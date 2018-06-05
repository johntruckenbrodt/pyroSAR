import os
import shutil
import pytest
import numpy as np
from osgeo import ogr
from pyroSAR import identify
from pyroSAR.spatial import crsConvert, haversine, Raster, stack, ogr2ogr, gdal_translate, gdal_rasterize, dtypes
from pyroSAR.ancillary import finder


def test_crsConvert():
    assert crsConvert(crsConvert(4326, 'wkt'), 'proj4') == '+proj=longlat +datum=WGS84 +no_defs '
    assert crsConvert(crsConvert(4326, 'prettyWkt'), 'opengis') == 'http://www.opengis.net/def/crs/EPSG/0/4326'
    assert crsConvert('+proj=longlat +datum=WGS84 +no_defs ', 'epsg') == 4326
    assert crsConvert('http://www.opengis.net/def/crs/EPSG/0/4326', 'epsg') == 4326
    assert crsConvert(crsConvert('http://www.opengis.net/def/crs/EPSG/0/4326', 'osr'), 'epsg') == 4326
    with pytest.raises(TypeError):
        crsConvert('xyz', 'epsg')
    with pytest.raises(ValueError):
        crsConvert(4326, 'xyz')


def test_haversine():
    assert haversine(50, 10, 51, 10) == 111194.92664455889


def test_Vector():
    scene = identify('pyroSAR/tests/data/S1A_IW_GRDH_1SDV_20150222T170750_20150222T170815_004739_005DD8_3768.zip')
    bbox = scene.bbox()
    assert bbox.getArea() == 7.573045244595988
    assert bbox.extent == {'ymax': 52.183979, 'ymin': 50.295261, 'xmin': 8.017178, 'xmax': 12.0268}
    assert bbox.nlayers == 1
    assert bbox.getProjection('epsg') == 4326
    assert bbox.proj4 == '+proj=longlat +datum=WGS84 +no_defs'
    assert isinstance(bbox.getFeatureByIndex(0), ogr.Feature)
    with pytest.raises(IndexError):
        bbox.getFeatureByIndex(1)
    bbox.reproject(32633)
    assert bbox.proj4 == '+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs'


def test_Raster():
    with pytest.raises(OSError):
        ras = Raster('foobar')
    ras = Raster('pyroSAR/tests/data/S1A__IW___A_20150309T173017_VV_grd_mli_geo_norm_db.tif')
    assert ras.bands == 1
    assert ras.proj4 == '+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs '
    assert ras.cols == 268
    assert ras.rows == 217
    assert ras.dim == [217, 268, 1]
    assert ras.dtype == 'Float32'
    assert dtypes(ras.dtype) == 6
    assert ras.epsg == 32631
    assert ras.format == 'GTiff'
    assert ras.geo == {'ymax': 4830114.70107, 'rotation_y': 0.0, 'rotation_x': 0.0, 'xmax': 625408.241204, 'xres': 20.0,
                       'xmin': 620048.241204, 'ymin': 4825774.70107, 'yres': -20.0}
    assert ras.typemap() == {'int32': 5, 'int16': 3, 'float64': 7, 'complex128': 11, 'uint8': 1, 'uint16': 2,
                             'complex64': 10, 'uint32': 4, 'int8': 1, 'float32': 6}
    assert ras.geogcs == 'WGS 84'
    assert ras.is_valid() is True
    assert ras.proj4args == {'units': 'm', 'no_defs': None, 'datum': 'WGS84', 'proj': 'utm', 'zone': '31'}
    assert ras.allstats == [[-26.65471076965332, 1.4325850009918213, -12.124929534450377, 4.738273594738293]]
    assert ras.bbox().getArea() == 23262400.0
    assert len(ras.layers()) == 1
    ras.load()
    mat = ras.matrix()
    assert isinstance(mat, np.ndarray)
    ras.assign(mat)
    # ras.reduce()
    ras.rescale(lambda x: 10 * x)
    ras.write('pyroSAR/tests/data/S1A__IW___A_20150309T173017_VV_grd_mli_geo_norm_db_new.tif')
    for item in finder('pyroSAR/tests/data', ['S1A*_new*']):
        os.remove(item)


def test_stack():
    name = 'pyroSAR/tests/data/S1A__IW___A_20150309T173017_VV_grd_mli_geo_norm_db.tif'
    outname = 'pyroSAR/tests/data/S1A__IW___A_20150309T173017_VV_grd_mli_geo_norm_db_sub'
    tr = (30, 30)
    with pytest.raises(IOError):
        stack(srcfiles=[], resampling='near', targetres=tr,
              srcnodata=-99, dstnodata=-99, dstfile=outname)

    with pytest.raises(IOError):
        stack(srcfiles=[name, name], resampling='near', targetres=tr,
              srcnodata=-99, dstnodata=-99, dstfile=outname, layernames=['a'])

    with pytest.raises(RuntimeError):
        stack(srcfiles=[name, name], resampling='near', targetres=30,
              srcnodata=-99, dstnodata=-99, dstfile=outname)

    with pytest.raises(RuntimeError):
        stack(srcfiles=[name, name], resampling='near', targetres=(30, 30, 30),
              srcnodata=-99, dstnodata=-99, dstfile=outname)

    with pytest.raises(IOError):
        stack(srcfiles=[name, name], resampling='foobar', targetres=tr,
              srcnodata=-99, dstnodata=-99, dstfile=outname)

    with pytest.raises(OSError):
        stack(srcfiles=[name, 'foobar'], resampling='near', targetres=tr,
              srcnodata=-99, dstnodata=-99, dstfile=outname)

    stack(srcfiles=[name, name], resampling='near', targetres=tr, overwrite=True,
          srcnodata=-99, dstnodata=-99, dstfile=outname)
    for item in finder('pyroSAR/tests/data', ['S1A*_sub*']):
        os.remove(item)


def test_auxil():
    dir = 'pyroSAR/tests/data/subdir'
    if not os.path.exists(dir):
        os.makedirs(dir)
    ras = Raster('pyroSAR/tests/data/S1A__IW___A_20150309T173017_VV_grd_mli_geo_norm_db.tif')
    bbox = os.path.join(dir, 'bbox.shp')
    ras.bbox(bbox)
    ogr2ogr(bbox, os.path.join(dir, 'bbox.gml'), {'format': 'GML'})
    gdal_translate(ras.raster, os.path.join(dir, 'test'), {'format': 'ENVI'})
    gdal_rasterize(bbox, os.path.join(dir, 'test2'), {'format': 'GTiff', 'xRes': 20, 'yRes': 20})
    shutil.rmtree(dir)
