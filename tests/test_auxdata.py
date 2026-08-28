import ftplib
import json
import os
import socket
import zipfile as zf
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from lxml import etree
from osgeo import gdal, osr

import pyroSAR.auxdata as auxdata
from pyroSAR.auxdata import (
    DEMHandler,
    ImplicitFTP_TLS,
    dem_autoload,
    dem_create,
    get_dem_options,
    get_egm_lookup,
    getasse30_hdr,
    vrt_check_sources,
)
from spatialist import Raster, bbox

DEM_TYPES = sorted([
    'AW3D30',
    'Copernicus 10m EEA DEM',
    'Copernicus 30m Global DEM',
    'Copernicus 30m Global DEM II',
    'Copernicus 90m Global DEM',
    'Copernicus 90m Global DEM II',
    'GETASSE30',
    'SRTM 1Sec HGT',
    'SRTM 3Sec',
])

DEM_TYPES_AUTH = sorted([
    'Copernicus 10m EEA DEM',
    'Copernicus 30m Global DEM II',
    'Copernicus 90m Global DEM II',
])

DEM_TYPES_NO_AUTH = sorted(set(DEM_TYPES) - set(DEM_TYPES_AUTH))


class _Response:
    def __init__(self, *, content=b'', text='', json_data=None, status_code=200):
        self.content = content
        self.text = text
        self._json_data = json_data
        self.status_code = status_code
        self.closed = False
    
    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f'HTTP {self.status_code}')
    
    def close(self):
        self.closed = True
    
    def json(self):
        return self._json_data


class _FakeVrtRaster:
    def __init__(self, xml):
        self.raster = SimpleNamespace(GetMetadata=lambda domain: [xml])
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        return None


@pytest.fixture
def raster_factory(tmp_path):
    """Create small georeferenced GeoTIFFs and return them as Raster objects."""
    rasters = []
    
    def create(
            name='dummy.tif',
            extent=(11, 51.98, 11.02, 52),
            shape=(2, 2),
            epsg=4326,
            dtype=gdal.GDT_Int16,
            nodata=None,
            value=1,
            array=None,
    ):
        filename = tmp_path / name
        filename.parent.mkdir(parents=True, exist_ok=True)
        
        if isinstance(extent, dict):
            xmin = extent['xmin']
            ymin = extent['ymin']
            xmax = extent['xmax']
            ymax = extent['ymax']
        else:
            xmin, ymin, xmax, ymax = extent
        
        if array is not None:
            array = np.asarray(array)
            rows, cols = array.shape
        else:
            rows, cols = shape
        
        dataset = gdal.GetDriverByName('GTiff').Create(
            str(filename),
            cols,
            rows,
            1,
            dtype,
        )
        dataset.SetGeoTransform((
            xmin,
            (xmax - xmin) / cols,
            0,
            ymax,
            0,
            -(ymax - ymin) / rows,
        ))
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(epsg)
        dataset.SetProjection(srs.ExportToWkt())
        
        band = dataset.GetRasterBand(1)
        if nodata is not None:
            band.SetNoDataValue(nodata)
        if array is None:
            band.Fill(value)
        else:
            band.WriteArray(array)
        band.FlushCache()
        dataset.FlushCache()
        band = None
        dataset = None
        
        raster = Raster(str(filename))
        rasters.append(raster)
        return raster
    
    yield create
    
    for raster in rasters:
        raster.close()


def _build_vrt(filename, sources):
    sources = [
        source.filename if isinstance(source, Raster) else str(source)
        for source in sources
    ]
    dataset = gdal.BuildVRT(str(filename), sources)
    assert dataset is not None
    dataset.FlushCache()
    dataset = None
    return str(filename)


def _dataset_bounds(dataset):
    gt = dataset.GetGeoTransform()
    xmin = gt[0]
    ymax = gt[3]
    xmax = xmin + dataset.RasterXSize * gt[1]
    ymin = ymax + dataset.RasterYSize * gt[5]
    return xmin, ymin, xmax, ymax


def _capture_gdalwarp(monkeypatch):
    """Patch gdalwarp and return a dictionary containing the call arguments."""
    captured = {}
    
    def fake_gdalwarp(*, src, dst, pbar=False, **kwargs):
        captured['src'] = src
        captured['dst'] = dst
        captured['pbar'] = pbar
        captured['kwargs'] = kwargs
    
    monkeypatch.setattr(auxdata, 'gdalwarp', fake_gdalwarp)
    return captured


def _dem_vrt(
        tmp_path,
        raster_factory,
        *,
        dem_type='SRTM 3Sec',
        product='dem',
        dtype=None,
        nodata=None,
        extent=(11, 51.98, 11.02, 52),
):
    """Create a test DEM tile with a recognizable filename and build a VRT."""
    basenames = {
        'SRTM 3Sec': {
            'dem': 'srtm_39_02.tif',
        },
        'Copernicus 30m Global DEM': {
            'dem': 'Copernicus_DSM_COG_10_N51_00_E011_00_DEM.tif',
            'wbm': 'Copernicus_DSM_COG_10_N51_00_E011_00_WBM.tif',
        },
        'GETASSE30': {
            'dem': '45N000E.GETASSE30',
        },
    }
    if dtype is None:
        dtype = gdal.GDT_Int16
    
    raster = raster_factory(
        name=f'{dem_type}/{basenames[dem_type][product]}',
        extent=extent,
        dtype=dtype,
        nodata=nodata,
    )
    vrt = tmp_path / f'{dem_type.replace(" ", "_")}_{product}.vrt'
    _build_vrt(vrt, [raster])
    return raster, str(vrt)


# -----------------------------------------------------------------------------
# Convenience wrappers


def test_dem_autoload_forwards_arguments(monkeypatch):
    geometry = object()
    captured = {}
    
    class FakeHandler:
        def __init__(self, vectorobject, buffer):
            captured['init'] = (vectorobject, buffer)
        
        def __enter__(self):
            return self
        
        def __exit__(self, exc_type, exc_val, exc_tb):
            return None
        
        def load(self, **kwargs):
            captured['load'] = kwargs
            return ['tile.tif']
    
    monkeypatch.setattr(auxdata, 'DEMHandler', FakeHandler)
    
    result = dem_autoload(
        vectorobject=geometry,
        demType='SRTM 3Sec',
        vrt='test.vrt',
        buffer=0.5,
        username='user',
        password='password',
        product='dem',
        crop=False,
        lock_timeout=12,
        offline=True,
        return_fname=False,
    )
    
    assert result == ['tile.tif']
    assert captured['init'] == (geometry, 0.5)
    assert captured['load'] == {
        'dem_type': 'SRTM 3Sec',
        'username': 'user',
        'password': 'password',
        'vrt': 'test.vrt',
        'product': 'dem',
        'crop': False,
        'lock_timeout': 12,
        'offline': True,
        'return_fname': False,
    }


def test_dem_create_forwards_arguments(monkeypatch):
    geometry = object()
    captured = {}
    
    class FakeHandler:
        def __init__(self, vectorobject, buffer):
            captured['init'] = (vectorobject, buffer)
        
        def __enter__(self):
            return self
        
        def __exit__(self, exc_type, exc_val, exc_tb):
            return None
        
        def create(self, **kwargs):
            captured['create'] = kwargs
    
    monkeypatch.setattr(auxdata, 'DEMHandler', FakeHandler)
    
    dem_create(
        vectorobject=geometry,
        src='test.vrt',
        dst='test.tif',
        buffer=0.25,
        t_srs=32632,
        tr=(30, 30),
        threads=2,
        geoid_convert=False,
        nodata=-32767,
        resampleAlg='cubic',
        dtype='Float32',
        pbar=True,
        creationOptions=['COMPRESS=LZW'],
    )
    
    assert captured['init'] == (geometry, 0.25)
    assert captured['create'] == {
        'src': 'test.vrt',
        'dst': 'test.tif',
        't_srs': 32632,
        'tr': (30, 30),
        'threads': 2,
        'geoid_convert': False,
        'nodata': -32767,
        'resampleAlg': 'cubic',
        'dtype': 'Float32',
        'pbar': True,
        'creationOptions': ['COMPRESS=LZW'],
    }


# -----------------------------------------------------------------------------
# DEMHandler initialization and configuration


def test_handler_init_geometry():
    extent = {'xmin': 11, 'xmax': 12, 'ymin': 51, 'ymax': 52}
    with bbox(coordinates=extent, crs=4326) as box:
        with DEMHandler(vectorobject=box) as handler:
            assert handler.extent == extent
            assert isinstance(handler.auxdatapath, str)


def test_handler_init_buffer():
    extent = {'xmin': 11, 'xmax': 12, 'ymin': 51, 'ymax': 52}
    expected = {'xmin': 10.5, 'xmax': 12.5, 'ymin': 50.5, 'ymax': 52.5}
    with bbox(coordinates=extent, crs=4326) as box:
        with DEMHandler(vectorobject=box, buffer=0.5) as handler:
            assert handler.extent == pytest.approx(expected)


def test_handler_init_global_extent():
    with DEMHandler() as handler:
        assert handler.extent == {
            'xmin': -180,
            'xmax': 180,
            'ymin': -90,
            'ymax': 90,
        }


def test_handler_init_invalid_geometry():
    with pytest.raises(
            RuntimeError,
            match='vectorobject.*must be of type Vector or None',
    ):
        DEMHandler('foobar')


def test_handler_init_invalid_crs():
    extent = {
        'xmin': -955867,
        'xmax': -915536,
        'ymin': -5915518,
        'ymax': -5863678,
    }
    with bbox(coordinates=extent, crs=32632) as box:
        with pytest.raises(RuntimeError, match="input vector object's CRS must be WGS84 LatLon"):
            DEMHandler(box)


def test_handler_init_auxdatapath_fallback(monkeypatch):
    class FakeExamineSnap:
        pass
    
    monkeypatch.setattr(auxdata, 'ExamineSnap', FakeExamineSnap)
    with DEMHandler() as handler:
        expected = os.path.join(os.path.expanduser('~'), '.snap', 'auxdata')
        assert handler.auxdatapath == expected


def test_handler_config():
    with DEMHandler() as handler:
        config = handler.config
    
    assert sorted(config) == DEM_TYPES
    required = {
        'url',
        'ocean_fill_value',
        'nodata',
        'resolution',
        'tilesize',
        'area_or_point',
        'vsi',
        'pattern',
        'datatype',
        'authentication',
        'vertical_datum',
    }
    for properties in config.values():
        assert required <= properties.keys()
        products = set(properties['pattern'])
        assert products == set(properties['datatype'])
        assert products == set(properties['nodata'])
        assert products == set(properties['ocean_fill_value'])
        assert isinstance(properties['authentication'], bool)


@pytest.mark.parametrize(
    'require_auth, expected',
    [
        (None, DEM_TYPES),
        (False, DEM_TYPES_NO_AUTH),
        (True, DEM_TYPES_AUTH),
    ],
)
def test_get_dem_options(require_auth, expected):
    assert get_dem_options(require_auth=require_auth) == expected


# -----------------------------------------------------------------------------
# Coordinate handling and resolution


@pytest.mark.parametrize(
    'extent, step, expected',
    [
        ({'xmin': 11, 'xmax': 12, 'ymin': 51, 'ymax': 51.5}, 1, ([51], [11])),
        ({'xmin': 11, 'xmax': 12, 'ymin': 51, 'ymax': 51.5}, 5, ([50], [10])),
        ({'xmin': 11, 'xmax': 12, 'ymin': 51, 'ymax': 51.5}, 15, ([45], [0])),
        ({'xmin': 179, 'xmax': -179, 'ymin': 51, 'ymax': 51.5}, 1, ([51], [179, -180])),
        ({'xmin': 175, 'xmax': -175, 'ymin': 51, 'ymax': 51.5}, 5, ([50], [175, -180])),
        ({'xmin': 165, 'xmax': -165, 'ymin': 51, 'ymax': 51.5}, 15, ([45], [165, -180])),
        ({'xmin': -59, 'xmax': -58, 'ymin': -52, 'ymax': -51.5}, 1, ([-52], [-59])),
    ],
)
def test_intrange(extent, step, expected):
    with bbox(coordinates=extent, crs=4326) as box:
        with DEMHandler(vectorobject=box) as handler:
            assert handler.intrange(step=step) == expected


@pytest.mark.parametrize(
    'latitude, expected',
    [
        (0, (1 / 3600, 1 / 3600)),
        (50, (1 / 3600, 1 / 3600)),
        (55, (1.5 / 3600, 1 / 3600)),
        (-65, (2 / 3600, 1 / 3600)),
        (75, (3 / 3600, 1 / 3600)),
        (82, (5 / 3600, 1 / 3600)),
        (87, (10 / 3600, 1 / 3600)),
        (90, (10 / 3600, 1 / 3600)),
    ],
)
def test_get_resolution(latitude, expected):
    with DEMHandler() as handler:
        result = handler._DEMHandler__get_resolution(
            dem_type='Copernicus 30m Global DEM',
            y=latitude,
        )
    assert result == pytest.approx(expected)


def test_get_resolution_outside_range():
    with DEMHandler() as handler:
        with pytest.raises(RuntimeError, match='could not get resolution'):
            handler._DEMHandler__get_resolution(
                dem_type='Copernicus 30m Global DEM',
                y=91,
            )


# -----------------------------------------------------------------------------
# Filename inference


@pytest.mark.parametrize(
    'filename, expected',
    [
        ('/tmp/AW3D30/N51E011_DSM.tif', ('AW3D30', 'dem')),
        ('/tmp/AW3D30/N51E011_MSK.tif', ('AW3D30', 'msk')),
        ('/tmp/AW3D30/N51E011_STK.tif', ('AW3D30', 'stk')),
        (
                '/tmp/Copernicus 30m Global DEM/Copernicus_DSM_N51_00_E011_00_DEM.tif',
                ('Copernicus 30m Global DEM', 'dem'),
        ),
        (
                '/tmp/Copernicus 30m Global DEM/Copernicus_DSM_N51_00_E011_00_EDM.tif',
                ('Copernicus 30m Global DEM', 'edm'),
        ),
        (
                '/tmp/Copernicus 30m Global DEM/Copernicus_DSM_N51_00_E011_00_FLM.tif',
                ('Copernicus 30m Global DEM', 'flm'),
        ),
        (
                '/tmp/Copernicus 30m Global DEM/Copernicus_DSM_N51_00_E011_00_HEM.tif',
                ('Copernicus 30m Global DEM', 'hem'),
        ),
        (
                '/tmp/Copernicus 30m Global DEM/Copernicus_DSM_N51_00_E011_00_WBM.tif',
                ('Copernicus 30m Global DEM', 'wbm'),
        ),
        ('/tmp/GETASSE30/45N000E.GETASSE30', ('GETASSE30', 'dem')),
        ('/tmp/SRTM 1Sec HGT/N51E011.hgt', ('SRTM 1Sec HGT', 'dem')),
        ('/tmp/SRTM 3Sec/srtm_39_03.tif', ('SRTM 3Sec', 'dem')),
        (r'C:\tmp\SRTM 1Sec HGT\N51E011.hgt', ('SRTM 1Sec HGT', 'dem')),
    ],
)
def test_info_from_filenames(filename, expected):
    with DEMHandler() as handler:
        assert handler.info_from_filenames([filename]) == expected


@pytest.mark.parametrize(
    'dem_type',
    ['Copernicus 30m Global DEM II', 'Copernicus 90m Global DEM II'],
)
def test_info_from_filenames_common_prefix(dem_type):
    filename = f'/tmp/{dem_type}/N51E011_DEM.tif'
    with DEMHandler() as handler:
        assert handler.info_from_filenames([filename]) == (dem_type, 'dem')


def test_info_from_filenames_multiple_files():
    filenames = [
        '/tmp/SRTM 1Sec HGT/N51E011.hgt',
        '/tmp/SRTM 1Sec HGT/N51E012.hgt',
    ]
    with DEMHandler() as handler:
        assert handler.info_from_filenames(filenames) == ('SRTM 1Sec HGT', 'dem')


@pytest.mark.parametrize(
    'filenames, match',
    [
        (['/tmp/N51E011.hgt'], 'could not infer DEM type'),
        (['/tmp/SRTM 1Sec HGT/N51E011.tif'], 'could not infer product type'),
        (
                ['/tmp/AW3D30/N51E011_DSM.tif', '/tmp/SRTM 1Sec HGT/N51E011.hgt'],
                'multiple DEM types found',
        ),
        (
                ['/tmp/AW3D30/N51E011_DSM.tif', '/tmp/AW3D30/N51E011_MSK.tif'],
                'multiple products found',
        ),
    ],
)
def test_info_from_filenames_errors(filenames, match):
    with DEMHandler() as handler:
        with pytest.raises(ValueError, match=match):
            handler.info_from_filenames(filenames)


# -----------------------------------------------------------------------------
# Remote product IDs


@pytest.mark.parametrize(
    'cases_fixture',
    [
        'auxdata_dem_cases_northern',
        'auxdata_dem_cases_northern_antimeridian',
        'auxdata_dem_cases_southern',
    ],
)
def test_remote_ids(cases_fixture, request):
    extent, cases = request.getfixturevalue(cases_fixture)
    with bbox(coordinates=extent, crs=4326) as box:
        with DEMHandler(box) as handler:
            for dem_type, reference in cases:
                assert handler.remote_ids(dem_type=dem_type) == reference


def test_remote_ids_invalid_dem_type():
    with DEMHandler() as handler:
        with pytest.raises(RuntimeError, match='not supported'):
            handler.remote_ids(dem_type='foobar')


def test_remote_ids_aw3d30():
    extent = {'xmin': 11.5, 'xmax': 11.9, 'ymin': 51, 'ymax': 51.5}
    with bbox(coordinates=extent, crs=4326) as box:
        with DEMHandler(box) as handler:
            assert handler.remote_ids('AW3D30') == ['N050E010/N051E011.tar.gz']


@pytest.mark.parametrize(
    'dem_type, index, expected',
    [
        ('SRTM 1Sec HGT', {'N51': {'E011': {'dem': 'srtm1'}}}, ['srtm1']),
        ('Copernicus 30m Global DEM', {'N51': {'E011': {'dem': 'cop30'}}}, ['cop30']),
        ('Copernicus 90m Global DEM', {'N51': {'E011': {'dem': 'cop90'}}}, ['cop90']),
        ('GETASSE30', {'45N': {'000E': {'dem': 'getasse'}}}, ['getasse']),
        ('SRTM 3Sec', {'02': {'39': {'dem': 'srtm3'}}}, ['srtm3']),
    ],
)
def test_remote_ids_from_local_index(monkeypatch, dem_type, index, expected):
    extent = {'xmin': 11.5, 'xmax': 11.9, 'ymin': 51, 'ymax': 51.5}
    with bbox(coordinates=extent, crs=4326) as box:
        with DEMHandler(box) as handler:
            monkeypatch.setattr(
                handler,
                '_DEMHandler__local_index',
                lambda dem_type: index,
            )
            result = handler.remote_ids(dem_type)
    assert result == expected


def test_remote_ids_missing_index_entry(monkeypatch):
    extent = {'xmin': 11.5, 'xmax': 11.9, 'ymin': 51, 'ymax': 51.5}
    with bbox(coordinates=extent, crs=4326) as box:
        with DEMHandler(box) as handler:
            monkeypatch.setattr(
                handler,
                '_DEMHandler__local_index',
                lambda dem_type: {},
            )
            assert handler.remote_ids('SRTM 1Sec HGT') == []


# -----------------------------------------------------------------------------
# Local product index


def test_local_index_existing(tmp_path, monkeypatch):
    with DEMHandler() as handler:
        handler.auxdatapath = str(tmp_path)
        path = tmp_path / 'dem' / 'SRTM 1Sec HGT' / 'index.json'
        path.parent.mkdir(parents=True)
        reference = {'N51': {'E011': {'dem': 'https://example.test/N51E011.zip'}}}
        path.write_text(json.dumps(reference))
        monkeypatch.setattr(
            auxdata.requests,
            'get',
            lambda *args, **kwargs: pytest.fail('network access not expected'),
        )
        result = handler._DEMHandler__local_index('SRTM 1Sec HGT')
    assert result == reference


@pytest.mark.parametrize(
    'dem_type, item, lat, lon',
    [
        ('GETASSE30', '45N000E.zip', '45N', '000E'),
        ('SRTM 1Sec HGT', 'N51E011.SRTMGL1.hgt.zip', 'N51', 'E011'),
        ('SRTM 3Sec', 'srtm_39_02.zip', '02', '39'),
    ],
)
def test_local_index_from_directory_listing(
        tmp_path,
        monkeypatch,
        dem_type,
        item,
        lat,
        lon,
):
    html = f'<html><body><a href="../">../</a><a href="{item}">{item}</a></body></html>'
    monkeypatch.setattr(auxdata.requests, 'get', lambda url: _Response(text=html))
    
    with DEMHandler() as handler:
        handler.auxdatapath = str(tmp_path)
        result = handler._DEMHandler__local_index(dem_type)
        url = handler.config[dem_type]['url'].rstrip('/') + '/' + item
    
    assert result[lat][lon]['dem'] == url
    assert (tmp_path / 'dem' / dem_type / 'index.json').is_file()


def test_local_index_copernicus_stac(tmp_path, monkeypatch):
    dem_type = 'Copernicus 30m Global DEM'
    catalog = 'dem_cop_30.json'
    item = 'Copernicus_DSM_COG_10_N51_00_E011_00.json'
    listing = f'''\
<ListBucketResult>
    <IsTruncated>false</IsTruncated>
    <Contents><Key>{catalog}</Key></Contents>
    <Contents><Key>{item}</Key></Contents>
</ListBucketResult>
'''.encode()
    asset_href = (
        'https://copernicus-dem-30m.s3.eu-central-1.amazonaws.com/'
        'Copernicus_DSM_COG_10_N51_00_E011_00_DEM/'
        'Copernicus_DSM_COG_10_N51_00_E011_00_DEM.tif'
    )
    
    def fake_get(url, params=None):
        if url.endswith(item):
            return _Response(json_data={'assets': {'elevation': {'href': asset_href}}})
        return _Response(content=listing)
    
    monkeypatch.setattr(auxdata.requests, 'get', fake_get)
    
    with DEMHandler() as handler:
        handler.auxdatapath = str(tmp_path)
        result = handler._DEMHandler__local_index(dem_type)
    
    tile = result['N51']['E011']
    assert tile['dem'].endswith(
        'Copernicus_DSM_COG_10_N51_00_E011_00_DEM/'
        'Copernicus_DSM_COG_10_N51_00_E011_00_DEM.tif'
    )
    assert tile['edm'].endswith('_DEM/AUXFILES/Copernicus_DSM_COG_10_N51_00_E011_00_EDM.tif')
    assert tile['flm'].endswith('_DEM/AUXFILES/Copernicus_DSM_COG_10_N51_00_E011_00_FLM.tif')
    assert tile['hem'].endswith('_DEM/AUXFILES/Copernicus_DSM_COG_10_N51_00_E011_00_HEM.tif')
    assert tile['wbm'].endswith('_DEM/AUXFILES/Copernicus_DSM_COG_10_N51_00_E011_00_WBM.tif')


def test_local_index_unsupported(tmp_path):
    with DEMHandler() as handler:
        handler.auxdatapath = str(tmp_path)
        with pytest.raises(RuntimeError, match='local indexing is not supported'):
            handler._DEMHandler__local_index('AW3D30')


# -----------------------------------------------------------------------------
# HTTP and FTP retrieval


def test_retrieve_empty(tmp_path, monkeypatch):
    monkeypatch.setattr(
        auxdata.requests,
        'get',
        lambda *args, **kwargs: pytest.fail('network access not expected'),
    )
    assert DEMHandler._DEMHandler__retrieve([], str(tmp_path)) == []


def test_retrieve_offline_existing(tmp_path):
    one = tmp_path / 'one.tif'
    two = tmp_path / 'two.tif'
    one.write_bytes(b'one')
    two.write_bytes(b'two')
    urls = [
        'https://example.test/two.tif',
        'https://example.test/one.tif',
        'https://example.test/one.tif',
    ]
    result = DEMHandler._DEMHandler__retrieve(urls, str(tmp_path), offline=True)
    assert result == sorted([str(one), str(two)])


def test_retrieve_offline_missing(tmp_path):
    with pytest.raises(RuntimeError, match='file not found locally'):
        DEMHandler._DEMHandler__retrieve(
            ['https://example.test/missing.tif'],
            str(tmp_path),
            offline=True,
        )


def test_retrieve_download(tmp_path, monkeypatch):
    base = _Response()
    remote = _Response(content=b'content')
    calls = []
    
    def fake_get(url):
        calls.append(url)
        return base if url == 'https://example.test' else remote
    
    monkeypatch.setattr(auxdata.requests, 'get', fake_get)
    result = DEMHandler._DEMHandler__retrieve(
        ['https://example.test/tile.tif'],
        str(tmp_path),
    )
    
    local = tmp_path / 'tile.tif'
    assert result == [str(local)]
    assert local.read_bytes() == b'content'
    assert calls == ['https://example.test', 'https://example.test/tile.tif']
    assert base.closed
    assert remote.closed


def test_retrieve_404_is_skipped(tmp_path, monkeypatch):
    def fake_get(url):
        if url == 'https://example.test':
            return _Response()
        return _Response(status_code=404)
    
    monkeypatch.setattr(auxdata.requests, 'get', fake_get)
    result = DEMHandler._DEMHandler__retrieve(
        ['https://example.test/missing.tif'],
        str(tmp_path),
    )
    assert result == []


def test_retrieve_ftp_offline_existing(tmp_path):
    local = tmp_path / 'tile.zip'
    local.write_bytes(b'content')
    result = DEMHandler._DEMHandler__retrieve_ftp(
        url='ftp://example.test/data',
        filenames=['tile.zip'],
        outdir=str(tmp_path),
        username=None,
        password=None,
        offline=True,
    )
    assert result == [str(local)]


@pytest.mark.parametrize('scheme', ['ftps', 'ftpes'])
def test_retrieve_ftp_authentication_required(tmp_path, scheme):
    with pytest.raises(ValueError, match='username or password'):
        DEMHandler._DEMHandler__retrieve_ftp(
            url=f'{scheme}://example.test/data',
            filenames=['tile.zip'],
            outdir=str(tmp_path),
            username=None,
            password=None,
        )


def test_retrieve_ftp_download(tmp_path, monkeypatch):
    instances = []
    
    class FakeFTP:
        def __init__(self, host, timeout):
            self.host = host
            self.timeout = timeout
            self.cwd_path = None
            self.closed = False
            instances.append(self)
        
        def login(self):
            return None
        
        def cwd(self, path):
            self.cwd_path = path
        
        def nlst(self, remote):
            return [remote]
        
        def retrbinary(self, command, callback):
            assert command == 'RETR tile.zip'
            callback(b'content')
        
        def close(self):
            self.closed = True
    
    monkeypatch.setattr(auxdata.ftplib, 'FTP', FakeFTP)
    result = DEMHandler._DEMHandler__retrieve_ftp(
        url='ftp://example.test/data',
        filenames=['tile.zip'],
        outdir=str(tmp_path),
        username=None,
        password=None,
    )
    
    assert result == [str(tmp_path / 'tile.zip')]
    assert (tmp_path / 'tile.zip').read_bytes() == b'content'
    assert instances[0].cwd_path == '/data'
    assert instances[0].closed


def test_retrieve_ftp_temporary_error(tmp_path, monkeypatch):
    class FakeFTP:
        def __init__(self, host, timeout):
            pass
        
        def login(self):
            return None
        
        def cwd(self, path):
            return None
        
        def nlst(self, remote):
            raise ftplib.error_temp('temporary failure')
        
        def close(self):
            return None
    
    monkeypatch.setattr(auxdata.ftplib, 'FTP', FakeFTP)
    result = DEMHandler._DEMHandler__retrieve_ftp(
        url='ftp://example.test/data',
        filenames=['tile.zip'],
        outdir=str(tmp_path),
        username=None,
        password=None,
    )
    assert result == []


# -----------------------------------------------------------------------------
# Dummy DEM creation and VRT construction


def test_create_dummy_dem_memory():
    with DEMHandler() as handler:
        handler.extent = {'xmin': 11, 'xmax': 12, 'ymin': 51, 'ymax': 52}
        dataset = handler._DEMHandler__create_dummy_dem(filename=None, fill_value=7)
    
    assert isinstance(dataset, gdal.Dataset)
    assert dataset.RasterXSize == 1
    assert dataset.RasterYSize == 1
    assert dataset.RasterCount == 1
    assert dataset.GetGeoTransform() == pytest.approx((11, 1, 0, 52, 0, -1))
    srs = osr.SpatialReference()
    srs.ImportFromWkt(dataset.GetProjection())
    assert srs.GetAuthorityCode(None) == '4326'
    band = dataset.GetRasterBand(1)
    assert band.DataType == gdal.GDT_Byte
    assert band.GetNoDataValue() == 255
    np.testing.assert_array_equal(band.ReadAsArray(), np.array([[7]], dtype=np.uint8))
    dataset = None


def test_create_dummy_dem_file(tmp_path):
    dst = tmp_path / 'dummy.tif'
    with DEMHandler() as handler:
        handler.extent = {'xmin': 11, 'xmax': 12, 'ymin': 51, 'ymax': 52}
        result = handler._DEMHandler__create_dummy_dem(filename=str(dst), fill_value=3)
    
    assert result is None
    assert dst.is_file()
    dataset = gdal.Open(str(dst))
    assert dataset.GetGeoTransform() == pytest.approx((11, 1, 0, 52, 0, -1))
    assert dataset.GetRasterBand(1).GetNoDataValue() == 255
    np.testing.assert_array_equal(dataset.ReadAsArray(), np.array([[3]], dtype=np.uint8))
    dataset = None


def test_create_dummy_dem_antimeridian():
    with DEMHandler() as handler:
        handler.extent = {'xmin': 179, 'xmax': -179, 'ymin': 51, 'ymax': 52}
        datasets = handler._DEMHandler__create_dummy_dem(filename=None, fill_value=5)
    
    assert isinstance(datasets, list)
    assert len(datasets) == 2
    assert datasets[0].GetGeoTransform() == pytest.approx((179, 1, 0, 52, 0, -1))
    assert datasets[1].GetGeoTransform() == pytest.approx((-180, 1, 0, 52, 0, -1))
    for dataset in datasets:
        np.testing.assert_array_equal(dataset.ReadAsArray(), np.array([[5]], dtype=np.uint8))
    datasets = None


def test_create_dummy_dem_antimeridian_file(tmp_path):
    with DEMHandler() as handler:
        handler.extent = {'xmin': 179, 'xmax': -179, 'ymin': 51, 'ymax': 52}
        with pytest.raises(ValueError, match='crossing the antimeridian'):
            handler._DEMHandler__create_dummy_dem(
                filename=str(tmp_path / 'dummy.tif'),
                fill_value=0,
            )


def test_create_dummy_dem_custom_extent():
    extent = {'xmin': 20, 'xmax': 21, 'ymin': 40, 'ymax': 41}
    with DEMHandler() as handler:
        dataset = handler._DEMHandler__create_dummy_dem(
            filename=None,
            fill_value=4,
            extent=extent,
        )
    
    assert dataset.GetGeoTransform() == pytest.approx((20, 1, 0, 41, 0, -1))
    np.testing.assert_array_equal(dataset.ReadAsArray(), np.array([[4]], dtype=np.uint8))
    dataset = None


def test_buildvrt_regular(tmp_path, raster_factory):
    res = 5 / 6000
    source = raster_factory(
        extent=(11, 52 - 2 * res, 11 + 2 * res, 52),
    )
    vrt = tmp_path / 'test.vrt'
    
    with DEMHandler() as handler:
        handler.extent = {'xmin': 11, 'xmax': 11.01, 'ymin': 51.99, 'ymax': 52}
        handler._DEMHandler__buildvrt(
            tiles=[source.filename],
            vrt=str(vrt),
            dem_type='SRTM 3Sec',
            product='dem',
            crop=True,
        )
    
    dataset = gdal.Open(str(vrt))
    assert _dataset_bounds(dataset) == pytest.approx((11, 51.99, 11.01, 52))
    dataset = None


def test_buildvrt_crop_false_expands_to_tile_grid(tmp_path, raster_factory):
    source = raster_factory()
    vrt = tmp_path / 'test.vrt'
    
    with DEMHandler() as handler:
        handler.extent = {'xmin': 11.1, 'xmax': 11.9, 'ymin': 51.1, 'ymax': 51.9}
        handler._DEMHandler__buildvrt(
            tiles=[source.filename],
            vrt=str(vrt),
            dem_type='SRTM 3Sec',
            product='dem',
            crop=False,
        )
    
    dataset = gdal.Open(str(vrt))
    assert _dataset_bounds(dataset) == pytest.approx((10, 50, 15, 55))
    dataset = None


def test_buildvrt_point_registration_shift(tmp_path, raster_factory):
    res = 1 / 3600
    source = raster_factory(
        extent=(11, 52 - 2 * res, 11 + 2 * res, 52),
    )
    vrt = tmp_path / 'test.vrt'
    
    with DEMHandler() as handler:
        handler.extent = {
            'xmin': 11,
            'xmax': 11 + 4 * res,
            'ymin': 52 - 4 * res,
            'ymax': 52,
        }
        handler._DEMHandler__buildvrt(
            tiles=[source.filename],
            vrt=str(vrt),
            dem_type='SRTM 1Sec HGT',
            product='dem',
            crop=True,
        )
    
    shift = res / 2
    dataset = gdal.Open(str(vrt))
    assert _dataset_bounds(dataset) == pytest.approx(
        (
            11 - shift,
            52 - 4 * res + shift,
            11 + 4 * res - shift,
            52 + shift,
        )
    )
    dataset = None


def test_buildvrt_empty_tiles_creates_dummy(tmp_path):
    vrt = tmp_path / 'test.vrt'
    with DEMHandler() as handler:
        handler.extent = {'xmin': 11, 'xmax': 11.01, 'ymin': 51.99, 'ymax': 52}
        handler._DEMHandler__buildvrt(
            tiles=[],
            vrt=str(vrt),
            dem_type='SRTM 3Sec',
            product='dem',
            crop=True,
        )
    
    assert vrt.is_file()
    assert Path(str(vrt).replace('.vrt', '_tmp.tif')).is_file()
    tree = etree.parse(str(vrt))
    assert tree.find('VRTRasterBand').attrib['dataType'] == 'Int16'


def test_buildvrt_antimeridian():
    with DEMHandler() as handler:
        handler.extent = {'xmin': 179, 'xmax': -179, 'ymin': 51, 'ymax': 52}
        with pytest.raises(RuntimeError, match='crossing the antimeridian'):
            handler._DEMHandler__buildvrt(
                tiles=[],
                vrt='/vsimem/test.vrt',
                dem_type='SRTM 3Sec',
                product='dem',
                crop=True,
            )


# -----------------------------------------------------------------------------
# DEM loading


def test_load_invalid_dem_type():
    with DEMHandler() as handler:
        with pytest.raises(RuntimeError, match='not supported'):
            handler.load('foobar')


def test_load_invalid_product():
    with DEMHandler() as handler:
        with pytest.raises(RuntimeError, match="Product 'foobar' is not available"):
            handler.load('SRTM 3Sec', product='foobar')


def test_load_https_returns_files(monkeypatch, tmp_path):
    local = str(tmp_path / 'tile.tif')
    with DEMHandler() as handler:
        handler.auxdatapath = str(tmp_path)
        monkeypatch.setattr(handler, 'remote_ids', lambda **kwargs: ['remote'])
        monkeypatch.setattr(handler, '_DEMHandler__retrieve', lambda **kwargs: [local])
        result = handler.load('Copernicus 30m Global DEM')
    assert result == [local]


def test_load_ftp_returns_vsi_paths(monkeypatch, tmp_path):
    local = str(tmp_path / 'N051E011.tar.gz')
    inside = local + '/N051E011/N051E011_DSM.tif'
    with DEMHandler() as handler:
        handler.auxdatapath = str(tmp_path)
        monkeypatch.setattr(handler, 'remote_ids', lambda **kwargs: ['remote.tar.gz'])
        monkeypatch.setattr(handler, '_DEMHandler__retrieve_ftp', lambda **kwargs: [local])
        monkeypatch.setattr(auxdata, 'finder', lambda *args, **kwargs: [inside])
        result = handler.load('AW3D30')
    assert result == ['/vsitar/' + inside]


def test_load_return_archive_names(monkeypatch, tmp_path):
    local = str(tmp_path / 'N051E011.tar.gz')
    inside = local + '/N051E011/N051E011_DSM.tif'
    with DEMHandler() as handler:
        handler.auxdatapath = str(tmp_path)
        monkeypatch.setattr(handler, 'remote_ids', lambda **kwargs: ['remote.tar.gz'])
        monkeypatch.setattr(handler, '_DEMHandler__retrieve_ftp', lambda **kwargs: [local])
        monkeypatch.setattr(auxdata, 'finder', lambda *args, **kwargs: [inside])
        result = handler.load('AW3D30', return_fname=False)
    assert result == [local]


def test_load_getasse_creates_header(monkeypatch, tmp_path):
    local = str(tmp_path / '45N000E.zip')
    calls = []
    with DEMHandler() as handler:
        handler.auxdatapath = str(tmp_path)
        monkeypatch.setattr(handler, 'remote_ids', lambda **kwargs: ['remote.zip'])
        monkeypatch.setattr(handler, '_DEMHandler__retrieve', lambda **kwargs: [local])
        monkeypatch.setattr(auxdata, 'getasse30_hdr', calls.append)
        monkeypatch.setattr(auxdata, 'finder', lambda *args, **kwargs: [])
        handler.load('GETASSE30')
    assert calls == [local]


def test_load_builds_vrt(monkeypatch, tmp_path):
    local = str(tmp_path / 'tile.tif')
    captured = {}
    with DEMHandler() as handler:
        handler.auxdatapath = str(tmp_path)
        monkeypatch.setattr(handler, 'remote_ids', lambda **kwargs: ['remote'])
        monkeypatch.setattr(handler, '_DEMHandler__retrieve', lambda **kwargs: [local])
        monkeypatch.setattr(
            handler,
            '_DEMHandler__buildvrt',
            lambda **kwargs: captured.update(kwargs),
        )
        result = handler.load(
            'Copernicus 30m Global DEM',
            vrt='test.vrt',
            crop=False,
            product='wbm',
        )
    
    assert result is None
    assert captured == {
        'tiles': [local],
        'vrt': 'test.vrt',
        'crop': False,
        'dem_type': 'Copernicus 30m Global DEM',
        'product': 'wbm',
    }


def test_load_forwards_ftp_options(monkeypatch, tmp_path):
    captured = {}
    with DEMHandler() as handler:
        handler.auxdatapath = str(tmp_path)
        monkeypatch.setattr(handler, 'remote_ids', lambda **kwargs: ['remote.tar.gz'])
        
        def fake_retrieve_ftp(**kwargs):
            captured.update(kwargs)
            return []
        
        monkeypatch.setattr(handler, '_DEMHandler__retrieve_ftp', fake_retrieve_ftp)
        handler.load(
            'AW3D30',
            username='user',
            password='password',
            lock_timeout=12,
            offline=True,
        )
    
    assert captured['url'].startswith('ftp://')
    assert captured['filenames'] == ['remote.tar.gz']
    assert captured['username'] == 'user'
    assert captured['password'] == 'password'
    assert captured['lock_timeout'] == 12
    assert captured['offline'] is True
    assert captured['port'] == 0


# -----------------------------------------------------------------------------
# DEM creation


def test_create_string_source_requires_vrt(tmp_path):
    with DEMHandler() as handler:
        with pytest.raises(RuntimeError, match='must be a VRT file'):
            handler.create(src='source.tif', dst=str(tmp_path / 'out.tif'))


def test_create_defaults(monkeypatch, tmp_path, raster_factory):
    _, vrt = _dem_vrt(tmp_path, raster_factory)
    captured = _capture_gdalwarp(monkeypatch)
    
    with DEMHandler() as handler:
        handler.create(
            src=vrt,
            dst=str(tmp_path / 'out.tif'),
            geoid_convert=False,
        )
    
    kwargs = captured['kwargs']
    assert kwargs['srcNodata'] is None
    assert kwargs['resampleAlg'] == 'bilinear'
    assert kwargs['targetAlignedPixels'] is True
    assert kwargs['outputType'] == gdal.GDT_Int16
    assert kwargs['dstNodata'] == -32768
    assert kwargs['xRes'] == pytest.approx(0.01)
    assert kwargs['yRes'] == pytest.approx(0.01)
    assert kwargs['srcSRS'] == 'EPSG:4326'
    assert kwargs['dstSRS'] == 'EPSG:4326'
    assert kwargs['multithread'] is True
    assert kwargs['format'] == 'GTiff'


def test_create_byte_defaults(monkeypatch, tmp_path, raster_factory):
    _, vrt = _dem_vrt(
        tmp_path,
        raster_factory,
        dem_type='Copernicus 30m Global DEM',
        product='wbm',
        dtype=gdal.GDT_Byte,
    )
    captured = _capture_gdalwarp(monkeypatch)
    
    with DEMHandler() as handler:
        handler.create(src=vrt, dst=str(tmp_path / 'out.tif'))
    
    kwargs = captured['kwargs']
    assert kwargs['resampleAlg'] == 'mode'
    assert kwargs['dstNodata'] == 255
    assert kwargs['outputType'] == gdal.GDT_Byte
    assert kwargs['srcSRS'] == 'EPSG:4326'


def test_create_dtype_override(monkeypatch, tmp_path, raster_factory):
    _, vrt = _dem_vrt(tmp_path, raster_factory)
    captured = _capture_gdalwarp(monkeypatch)
    
    with DEMHandler() as handler:
        handler.create(
            src=vrt,
            dst=str(tmp_path / 'out.tif'),
            dtype='Float32',
            geoid_convert=False,
        )
    
    assert captured['kwargs']['outputType'] == gdal.GDT_Float32


@pytest.mark.parametrize(
    'threads, multithread',
    [(1, False), (2, True), ('ALL_CPUS', True)],
)
def test_create_threads(monkeypatch, tmp_path, raster_factory, threads, multithread):
    _, vrt = _dem_vrt(tmp_path, raster_factory)
    captured = _capture_gdalwarp(monkeypatch)
    monkeypatch.setattr(auxdata.gdal, '__version__', '3.12.1')
    
    with DEMHandler() as handler:
        handler.create(
            src=vrt,
            dst=str(tmp_path / 'out.tif'),
            threads=threads,
            geoid_convert=False,
        )
    
    assert captured['kwargs']['multithread'] is multithread
    assert captured['kwargs']['warpOptions'] == {'NUM_THREADS': str(threads)}


@pytest.mark.parametrize(
    'threads, exception, match',
    [
        ('2', ValueError, 'unsupported value'),
        (0, ValueError, 'must be >= 1'),
        (-1, ValueError, 'must be >= 1'),
        (1.5, TypeError, 'must be of type int, str or None'),
    ],
)
def test_create_invalid_threads(
        tmp_path,
        raster_factory,
        threads,
        exception,
        match,
):
    _, vrt = _dem_vrt(tmp_path, raster_factory)
    with DEMHandler() as handler:
        with pytest.raises(exception, match=match):
            handler.create(
                src=vrt,
                dst=str(tmp_path / 'out.tif'),
                threads=threads,
                geoid_convert=False,
            )


def test_create_rejects_locked_kwargs(tmp_path, raster_factory):
    _, vrt = _dem_vrt(tmp_path, raster_factory)
    with DEMHandler() as handler:
        with pytest.raises(RuntimeError, match="argument 'xRes' cannot be set"):
            handler.create(
                src=vrt,
                dst=str(tmp_path / 'out.tif'),
                geoid_convert=False,
                xRes=30,
            )


def test_create_forwards_additional_kwargs(monkeypatch, tmp_path, raster_factory):
    _, vrt = _dem_vrt(tmp_path, raster_factory)
    captured = _capture_gdalwarp(monkeypatch)
    
    with DEMHandler() as handler:
        handler.create(
            src=vrt,
            dst=str(tmp_path / 'out.tif'),
            geoid_convert=False,
            creationOptions=['COMPRESS=LZW'],
        )
    
    assert captured['kwargs']['creationOptions'] == ['COMPRESS=LZW']


def test_create_existing_output_is_not_overwritten(monkeypatch, tmp_path, raster_factory):
    _, vrt = _dem_vrt(tmp_path, raster_factory)
    dst = tmp_path / 'out.tif'
    dst.write_bytes(b'existing')
    monkeypatch.setattr(
        auxdata,
        'gdalwarp',
        lambda **kwargs: pytest.fail('gdalwarp must not be called'),
    )
    
    with DEMHandler() as handler:
        handler.create(src=vrt, dst=str(dst), geoid_convert=False)
    
    assert dst.read_bytes() == b'existing'


def test_create_removes_partial_output_on_error(monkeypatch, tmp_path, raster_factory):
    _, vrt = _dem_vrt(tmp_path, raster_factory)
    dst = tmp_path / 'out.tif'
    
    def fake_gdalwarp(*, dst, **kwargs):
        Path(dst).write_bytes(b'partial')
        raise RuntimeError('warp failed')
    
    monkeypatch.setattr(auxdata, 'gdalwarp', fake_gdalwarp)
    
    with DEMHandler() as handler:
        with pytest.raises(RuntimeError, match='warp failed'):
            handler.create(src=vrt, dst=str(dst), geoid_convert=False)
    
    assert not dst.exists()


def test_create_restores_gdal_num_threads(monkeypatch, tmp_path, raster_factory):
    _, vrt = _dem_vrt(tmp_path, raster_factory)
    _capture_gdalwarp(monkeypatch)
    monkeypatch.setattr(auxdata.gdal, '__version__', '3.12.1')
    
    original = gdal.GetConfigOption('GDAL_NUM_THREADS')
    gdal.SetConfigOption('GDAL_NUM_THREADS', '7')
    try:
        with DEMHandler() as handler:
            handler.create(
                src=vrt,
                dst=str(tmp_path / 'out.tif'),
                threads=2,
                geoid_convert=False,
            )
        assert gdal.GetConfigOption('GDAL_NUM_THREADS') == '7'
    finally:
        gdal.SetConfigOption('GDAL_NUM_THREADS', original)


def test_create_list_source_adds_dummy_and_bounds(monkeypatch, tmp_path, raster_factory):
    source = raster_factory(
        name='SRTM 3Sec/srtm_39_02.tif',
        extent=(11, 51, 12, 52),
    )
    sources = [source.filename]
    captured = _capture_gdalwarp(monkeypatch)
    
    with DEMHandler() as handler:
        handler.extent = {'xmin': 11, 'xmax': 12, 'ymin': 51, 'ymax': 52}
        handler.create(
            src=sources,
            dst=str(tmp_path / 'out.tif'),
            tr=(0.5, 0.5),
            geoid_convert=False,
        )
    
    assert sources == [source.filename]
    assert isinstance(captured['src'][0], gdal.Dataset)
    assert captured['src'][1:] == [source.filename]
    assert _dataset_bounds(captured['src'][0]) == pytest.approx((11, 51, 12, 52))
    assert captured['kwargs']['outputBounds'] == pytest.approx([11, 51, 12, 52])


def test_create_list_source_output_bounds_override(monkeypatch, tmp_path, raster_factory):
    source = raster_factory(
        name='SRTM 3Sec/srtm_39_02.tif',
        extent=(11, 51, 12, 52),
    )
    bounds = [10, 50, 13, 53]
    captured = _capture_gdalwarp(monkeypatch)
    
    with DEMHandler() as handler:
        handler.extent = {'xmin': 11, 'xmax': 12, 'ymin': 51, 'ymax': 52}
        handler.create(
            src=[source.filename],
            dst=str(tmp_path / 'out.tif'),
            tr=(0.5, 0.5),
            outputBounds=bounds,
            geoid_convert=False,
        )
    
    assert captured['kwargs']['outputBounds'] == bounds


def test_create_list_source_no_intersection(tmp_path, raster_factory):
    source = raster_factory(
        name='SRTM 3Sec/srtm_39_02.tif',
        extent=(11, 51, 12, 52),
    )
    
    with DEMHandler() as handler:
        handler.extent = {'xmin': 20, 'xmax': 21, 'ymin': 60, 'ymax': 61}
        with pytest.raises(RuntimeError, match='does not intersect'):
            handler.create(
                src=[source.filename],
                dst=str(tmp_path / 'out.tif'),
                geoid_convert=False,
            )


def test_create_list_source_antimeridian(monkeypatch, tmp_path, raster_factory):
    source = raster_factory(
        name='SRTM 3Sec/srtm_72_02.tif',
        extent=(179, 51, 180, 52),
    )
    
    class Intersection:
        extent = {'xmin': 179, 'xmax': -179, 'ymin': 51, 'ymax': 52}
        
        def reproject(self, crs):
            return None
        
        def close(self):
            return None
    
    monkeypatch.setattr(auxdata, 'intersect', lambda *args, **kwargs: Intersection())
    
    with DEMHandler() as handler:
        with pytest.raises(RuntimeError, match='crossing the antimeridian'):
            handler.create(
                src=[source.filename],
                dst=str(tmp_path / 'out.tif'),
                geoid_convert=False,
            )


def test_create_geoid_conversion(monkeypatch, tmp_path, raster_factory):
    _, vrt = _dem_vrt(tmp_path, raster_factory)
    captured = _capture_gdalwarp(monkeypatch)
    calls = []
    monkeypatch.setattr(
        auxdata,
        'get_egm_lookup',
        lambda **kwargs: calls.append(kwargs),
    )
    
    with DEMHandler() as handler:
        handler.create(
            src=vrt,
            dst=str(tmp_path / 'out.tif'),
            geoid_convert=True,
        )
    
    assert calls == [{'geoid': 'EGM96', 'software': 'PROJ'}]
    assert captured['kwargs']['srcSRS'] == 'EPSG:4326+5773'


def test_create_geoid_lookup_failure(monkeypatch, tmp_path, raster_factory):
    _, vrt = _dem_vrt(tmp_path, raster_factory)
    
    def fail(**kwargs):
        raise OSError('lookup unavailable')
    
    monkeypatch.setattr(auxdata, 'get_egm_lookup', fail)
    
    with DEMHandler() as handler:
        with pytest.raises(RuntimeError, match='lookup unavailable'):
            handler.create(
                src=vrt,
                dst=str(tmp_path / 'out.tif'),
                geoid_convert=True,
            )


@pytest.mark.parametrize(
    'dem_type, product, dtype',
    [
        ('Copernicus 30m Global DEM', 'wbm', gdal.GDT_Byte),
        ('GETASSE30', 'dem', gdal.GDT_Int16),
    ],
)
def test_create_disables_unneeded_geoid_conversion(
        monkeypatch,
        tmp_path,
        raster_factory,
        dem_type,
        product,
        dtype,
):
    _, vrt = _dem_vrt(
        tmp_path,
        raster_factory,
        dem_type=dem_type,
        product=product,
        dtype=dtype,
    )
    captured = _capture_gdalwarp(monkeypatch)
    calls = []
    monkeypatch.setattr(
        auxdata,
        'get_egm_lookup',
        lambda **kwargs: calls.append(kwargs),
    )
    
    with DEMHandler() as handler:
        handler.create(
            src=vrt,
            dst=str(tmp_path / 'out.tif'),
            geoid_convert=True,
        )
    
    assert calls == []
    assert captured['kwargs']['srcSRS'] == 'EPSG:4326'


def test_create_vrt_multithreading_workaround(monkeypatch, tmp_path, raster_factory):
    _, vrt = _dem_vrt(tmp_path, raster_factory)
    translated = object()
    translate_calls = []
    captured = _capture_gdalwarp(monkeypatch)
    
    monkeypatch.setattr(auxdata.gdal, '__version__', '3.11.0')
    monkeypatch.setattr(
        auxdata.psutil,
        'virtual_memory',
        lambda: SimpleNamespace(available=10 ** 9),
    )
    
    def fake_translate(*, destName, srcDS, format):
        translate_calls.append((destName, srcDS, format))
        return translated
    
    monkeypatch.setattr(auxdata.gdal, 'Translate', fake_translate)
    
    with DEMHandler() as handler:
        handler.create(
            src=vrt,
            dst=str(tmp_path / 'out.tif'),
            threads=2,
            geoid_convert=False,
        )
    
    assert translate_calls == [('', vrt, 'MEM')]
    assert captured['src'] is translated


def test_dem_create_local_vrt(tmp_path, raster_factory):
    source = raster_factory(
        name='SRTM 3Sec/srtm_39_02.tif',
        nodata=-32768,
        array=np.array([[1, 2], [3, 4]], dtype=np.int16),
    )
    vrt = tmp_path / 'source.vrt'
    out = tmp_path / 'out.tif'
    _build_vrt(vrt, [source])
    
    extent = {'xmin': 11, 'xmax': 11.02, 'ymin': 51.98, 'ymax': 52}
    with bbox(coordinates=extent, crs=4326) as box:
        dem_create(
            vectorobject=box,
            src=str(vrt),
            dst=str(out),
            t_srs=4326,
            tr=(0.01, 0.01),
            geoid_convert=False,
        )
    
    assert out.is_file()
    dataset = gdal.Open(str(out))
    assert dataset is not None
    assert dataset.RasterCount == 1
    assert dataset.GetRasterBand(1).DataType == gdal.GDT_Int16
    dataset = None


# -----------------------------------------------------------------------------
# GETASSE30 header creation


@pytest.mark.parametrize(
    'filename, longitude, upper_latitude',
    [
        ('45N000E.zip', '0.0', '60.0'),
        ('45S015W.zip', '-15.0', '-30.0'),
    ],
)
def test_getasse30_hdr(tmp_path, filename, longitude, upper_latitude):
    path = tmp_path / filename
    with zf.ZipFile(path, 'w') as archive:
        archive.writestr(filename.replace('.zip', '.GETASSE30'), b'')
    
    getasse30_hdr(str(path))
    
    hdr_name = filename.replace('.zip', '.hdr')
    with zf.ZipFile(path) as archive:
        assert hdr_name in archive.namelist()
        text = archive.read(hdr_name).decode()
    
    assert 'samples = 1800' in text
    assert 'lines = 1800' in text
    assert 'byte order = 1' in text
    assert 'data type = 2' in text
    assert 'Geographic Lat/Lon' in text
    assert longitude in text
    assert upper_latitude in text
    assert 'WGS-84' in text


def test_getasse30_hdr_existing(tmp_path):
    path = tmp_path / '45N000E.zip'
    with zf.ZipFile(path, 'w') as archive:
        archive.writestr('45N000E.GETASSE30', b'')
    getasse30_hdr(str(path))
    getasse30_hdr(str(path))
    with zf.ZipFile(path) as archive:
        assert archive.namelist().count('45N000E.hdr') == 1


# -----------------------------------------------------------------------------
# EGM lookup files


def test_get_egm_lookup_proj_existing(tmp_path, monkeypatch):
    monkeypatch.setenv('PROJ_DATA', str(tmp_path))
    monkeypatch.delenv('PROJ_LIB', raising=False)
    (tmp_path / 'us_nga_egm96_15.tif').write_bytes(b'existing')
    monkeypatch.setattr(
        auxdata.requests,
        'get',
        lambda *args, **kwargs: pytest.fail('network access not expected'),
    )
    get_egm_lookup(geoid='EGM96', software='PROJ')


def test_get_egm_lookup_proj_download(tmp_path, monkeypatch):
    monkeypatch.setenv('PROJ_DATA', str(tmp_path))
    monkeypatch.delenv('PROJ_LIB', raising=False)
    response = _Response(content=b'grid')
    monkeypatch.setattr(auxdata.requests, 'get', lambda url: response)
    get_egm_lookup(geoid='EGM2008', software='PROJ')
    assert (tmp_path / 'us_nga_egm08_25.tif').read_bytes() == b'grid'
    assert response.closed


def test_get_egm_lookup_proj_lib_fallback(tmp_path, monkeypatch):
    monkeypatch.delenv('PROJ_DATA', raising=False)
    monkeypatch.setenv('PROJ_LIB', str(tmp_path))
    (tmp_path / 'us_nga_egm96_15.tif').write_bytes(b'existing')
    get_egm_lookup(geoid='EGM96', software='PROJ')


def test_get_egm_lookup_proj_environment_missing(monkeypatch):
    monkeypatch.delenv('PROJ_DATA', raising=False)
    monkeypatch.delenv('PROJ_LIB', raising=False)
    with pytest.raises(
            RuntimeError,
            match="Neither environment variable 'PROJ_DATA' nor 'PROJ_LIB'",
    ):
        get_egm_lookup(geoid='EGM96', software='PROJ')


def test_get_egm_lookup_proj_not_writable(tmp_path, monkeypatch):
    monkeypatch.setenv('PROJ_DATA', str(tmp_path))
    monkeypatch.delenv('PROJ_LIB', raising=False)
    monkeypatch.setattr(auxdata.os, 'access', lambda path, mode: False)
    with pytest.raises(OSError, match='cannot write'):
        get_egm_lookup(geoid='EGM96', software='PROJ')


def test_get_egm_lookup_snap_download(tmp_path, monkeypatch):
    class FakeExamineSnap:
        auxdatapath = str(tmp_path)
    
    response = _Response(content=b'grid')
    monkeypatch.setattr(auxdata, 'ExamineSnap', FakeExamineSnap)
    monkeypatch.setattr(auxdata.requests, 'get', lambda url: response)
    get_egm_lookup(geoid='EGM96', software='SNAP')
    expected = tmp_path / 'dem' / 'egm96' / 'ww15mgh_b.zip'
    assert expected.read_bytes() == b'grid'
    assert response.closed


def test_get_egm_lookup_invalid_software():
    with pytest.raises(TypeError, match="software must be either 'SNAP' or 'PROJ'"):
        get_egm_lookup(geoid='EGM96', software='foobar')


# -----------------------------------------------------------------------------
# Implicit FTPS socket handling


def test_implicit_ftp_tls_sock_none():
    ftp = ImplicitFTP_TLS()
    assert ftp.sock is None
    ftp.sock = None
    assert ftp.sock is None


def test_implicit_ftp_tls_wraps_socket():
    ftp = ImplicitFTP_TLS()
    wrapped = object()
    ftp.context = SimpleNamespace(wrap_socket=lambda value: wrapped)
    raw = socket.socket()
    try:
        ftp.sock = raw
        assert ftp.sock is wrapped
    finally:
        raw.close()


# -----------------------------------------------------------------------------
# VRT source validation


def test_vrt_check_sources(tmp_path, raster_factory):
    src = raster_factory(name='source.tif')
    vrt = tmp_path / 'test.vrt'
    _build_vrt(vrt, [src])
    
    tree = etree.parse(str(vrt))
    element = tree.find('.//SourceFilename')
    element.text = src.filename
    element.attrib['relativeToVRT'] = '0'
    tree.write(str(vrt))
    
    assert vrt_check_sources(str(vrt)) == [src.filename]


def test_vrt_check_sources_relative(tmp_path, raster_factory):
    src = raster_factory(name='source.tif')
    vrt = tmp_path / 'test.vrt'
    _build_vrt(vrt, [src])
    
    tree = etree.parse(str(vrt))
    element = tree.find('.//SourceFilename')
    element.text = Path(src.filename).name
    element.attrib['relativeToVRT'] = '1'
    tree.write(str(vrt))
    
    assert vrt_check_sources(str(vrt)) == ['source.tif']


def test_vrt_check_sources_multiple(tmp_path, raster_factory):
    src1 = raster_factory(name='source1.tif')
    src2 = raster_factory(
        name='source2.tif',
        extent=(11.02, 51.98, 11.04, 52),
    )
    vrt = tmp_path / 'test.vrt'
    _build_vrt(vrt, [src1, src2])
    
    tree = etree.parse(str(vrt))
    sources = [src1.filename, src2.filename]
    for element, source in zip(tree.findall('.//SourceFilename'), sources):
        element.text = source
        element.attrib['relativeToVRT'] = '0'
    tree.write(str(vrt))
    
    assert vrt_check_sources(str(vrt)) == sources


def test_vrt_check_sources_missing(tmp_path, monkeypatch):
    src = tmp_path / 'missing.tif'
    xml = f'''\
<VRTDataset>
    <VRTRasterBand>
        <SimpleSource><SourceFilename>{src}</SourceFilename></SimpleSource>
    </VRTRasterBand>
</VRTDataset>
'''
    monkeypatch.setattr(auxdata, 'Raster', lambda filename: _FakeVrtRaster(xml))
    with pytest.raises(RuntimeError, match='missing VRT source file'):
        vrt_check_sources('test.vrt')


def test_vrt_check_sources_invalid(tmp_path):
    with pytest.raises(RuntimeError, match='cannot read'):
        vrt_check_sources(str(tmp_path / 'missing.vrt'))


def test_vrt_check_sources_none(monkeypatch):
    xml = '''\
<VRTDataset>
    <VRTRasterBand>
        <SimpleSource><SourceFilename></SourceFilename></SimpleSource>
    </VRTRasterBand>
</VRTDataset>
'''
    monkeypatch.setattr(auxdata, 'Raster', lambda filename: _FakeVrtRaster(xml))
    with pytest.raises(ValueError, match='None value as source file name'):
        vrt_check_sources('test.vrt')


def test_vrt_check_sources_invalid_archive_path(monkeypatch):
    xml = '''\
<VRTDataset>
    <VRTRasterBand>
        <SimpleSource>
            <SourceFilename>/vsizip//tmp/archive/image.tif</SourceFilename>
        </SimpleSource>
    </VRTRasterBand>
</VRTDataset>
'''
    monkeypatch.setattr(auxdata, 'Raster', lambda filename: _FakeVrtRaster(xml))
    with pytest.raises(RuntimeError, match='could not match archive path'):
        vrt_check_sources('test.vrt')


def test_vrt_check_sources_vsizip(tmp_path, raster_factory):
    src = raster_factory(name='source.tif')
    archive = tmp_path / 'source.zip'
    vrt = tmp_path / 'test.vrt'
    source_path = Path(src.filename)
    with zf.ZipFile(archive, 'w') as obj:
        obj.write(source_path, arcname=source_path.name)
    
    vsi = '/vsizip/' + str(archive).replace('\\', '/') + '/' + source_path.name
    _build_vrt(vrt, [vsi])
    assert vrt_check_sources(str(vrt)) == [vsi]


def test_vrt_check_sources_vsizip_missing_archive(tmp_path, monkeypatch):
    archive = tmp_path / 'missing.zip'
    vsi = '/vsizip/' + str(archive).replace('\\', '/') + '/source.tif'
    xml = f'''\
<VRTDataset>
    <VRTRasterBand>
        <SimpleSource><SourceFilename>{vsi}</SourceFilename></SimpleSource>
    </VRTRasterBand>
</VRTDataset>
'''
    monkeypatch.setattr(auxdata, 'Raster', lambda filename: _FakeVrtRaster(xml))
    with pytest.raises(RuntimeError, match='missing VRT source file'):
        vrt_check_sources('test.vrt')


# -----------------------------------------------------------------------------
# Online integration tests for unauthenticated services


def test_autoload_aw3d30(travis):
    if travis:
        pytest.skip('Travis CI does not support FTP access')
    with bbox(
            {'xmin': 11.5, 'xmax': 11.9, 'ymin': 51, 'ymax': 51.5},
            crs=4326,
    ) as box:
        files = dem_autoload(box, 'AW3D30')
    assert len(files) == 1


def test_autoload_aw3d30_stk(travis):
    if travis:
        pytest.skip('Travis CI does not support FTP access')
    with bbox(
            {'xmin': 11.5, 'xmax': 11.9, 'ymin': 51, 'ymax': 51.5},
            crs=4326,
    ) as box:
        files = dem_autoload(box, 'AW3D30', product='stk')
    assert len(files) == 1


def test_autoload_srtm1_online_and_offline():
    with bbox(
            {'xmin': 11.5, 'xmax': 11.9, 'ymin': 51, 'ymax': 51.5},
            crs=4326,
    ) as box:
        online = dem_autoload(box, 'SRTM 1Sec HGT')
        offline = dem_autoload(box, 'SRTM 1Sec HGT', offline=True)
    assert len(online) == 1
    assert offline == online


def test_autoload_srtm3():
    with bbox(
            {'xmin': 11.5, 'xmax': 11.9, 'ymin': 51, 'ymax': 51.5},
            crs=4326,
    ) as box:
        files = dem_autoload(box, 'SRTM 3Sec')
    assert len(files) == 1


def test_autoload_ocean():
    with bbox(
            {'xmin': -30, 'xmax': -29, 'ymin': 40, 'ymax': 41},
            crs=4326,
    ) as box:
        files = dem_autoload(box, 'SRTM 1Sec HGT')
    assert files == []
