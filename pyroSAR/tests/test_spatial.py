import os
import pytest
import numpy as np
from osgeo import ogr
from pyroSAR import identify
from pyroSAR.spatial import crsConvert, haversine, Raster, stack, ogr2ogr, gdal_translate, gdal_rasterize, dtypes, bbox
from pyroSAR.spatial.vector import feature2vector, dissolve, Vector, intersect


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


def test_Vector(testdata):
    scene = identify(testdata['s1'])
    bbox1 = scene.bbox()
    assert bbox1.getArea() == 7.573045244595988
    assert bbox1.extent == {'ymax': 52.183979, 'ymin': 50.295261, 'xmin': 8.017178, 'xmax': 12.0268}
    assert bbox1.nlayers == 1
    assert bbox1.getProjection('epsg') == 4326
    assert bbox1.proj4 == '+proj=longlat +datum=WGS84 +no_defs'
    assert isinstance(bbox1.wkt, str)
    assert isinstance(bbox1.getFeatureByIndex(0), ogr.Feature)
    with pytest.raises(IndexError):
        bbox1.getFeatureByIndex(1)
    bbox1.reproject(32633)
    assert bbox1.proj4 == '+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs'
    assert isinstance(bbox1['id=1'], Vector)
    with pytest.raises(RuntimeError):
        test = bbox1[0.1]
    assert bbox1.fieldnames == ['id']
    assert bbox1.getUniqueAttributes('id') == [1]
    feat = bbox1.getFeatureByAttribute('id', 1)
    assert isinstance(feat, ogr.Feature)
    bbox2 = feature2vector(feat, ref=bbox1)
    bbox2.close()
    feat.Destroy()
    with pytest.raises(KeyError):
        select = bbox1.getFeatureByAttribute('foo', 'bar')
    with pytest.raises(RuntimeError):
        vec = Vector(driver='foobar')
    bbox1.close()


def test_dissolve(tmpdir, travis, testdata):
    scene = identify(testdata['s1'])
    bbox1 = scene.bbox()
    # retrieve extent and shift its coordinates by one degree
    ext = bbox1.extent
    for key in ext.keys():
        ext[key] += 1
    # create new bbox shapefile with modified extent
    bbox2_name = os.path.join(str(tmpdir), 'bbox2.shp')
    bbox(ext, bbox1.srs, bbox2_name)
    # assert intersection between the two bboxes and combine them into one
    with Vector(bbox2_name) as bbox2:
        assert intersect(bbox1, bbox2) is not None
        bbox1.addvector(bbox2)
        # write combined bbox into new shapefile
        bbox3_name = os.path.join(str(tmpdir), 'bbox3.shp')
        bbox1.write(bbox3_name)

    if not travis:
        # dissolve the geometries in bbox3 and write the result to new bbox4
        # this test is currently disabled for Travis as the current sqlite3 version on Travis seems to not support
        # loading gdal as extension; Travis CI setup: Ubuntu 14.04 (Trusty), sqlite3 version 3.8.2 (2018-06-04)
        bbox4_name = os.path.join(str(tmpdir), 'bbox4.shp')
        dissolve(bbox3_name, bbox4_name, field='id')
        assert os.path.isfile(bbox4_name)


def test_Raster(tmpdir, testdata):
    with pytest.raises(OSError):
        ras = Raster('foobar')
    with Raster(testdata['tif']) as ras:
        print(ras)
        assert ras.bands == 1
        assert ras.proj4 == '+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs '
        assert ras.cols == 268
        assert ras.rows == 217
        assert ras.dim == [217, 268, 1]
        assert ras.dtype == 'Float32'
        assert ras.epsg == 32631
        assert ras.format == 'GTiff'
        assert ras.geo == {'ymax': 4830114.70107, 'rotation_y': 0.0, 'rotation_x': 0.0, 'xmax': 625408.241204,
                           'xres': 20.0, 'xmin': 620048.241204, 'ymin': 4825774.70107, 'yres': -20.0}
        assert ras.geogcs == 'WGS 84'
        assert ras.is_valid() is True
        assert ras.proj4args == {'units': 'm', 'no_defs': None, 'datum': 'WGS84', 'proj': 'utm', 'zone': '31'}
        assert ras.allstats == [[-26.65471076965332, 1.4325850009918213, -12.124929534450377, 4.738273594738293]]
        assert ras.bbox().getArea() == 23262400.0
        assert len(ras.layers()) == 1
        assert ras.projcs == 'WGS 84 / UTM zone 31N'
        assert ras.res == [20.0, 20.0]
        ras.load()
        mat = ras.matrix()
        assert isinstance(mat, np.ndarray)
        ras.assign(mat, index=0)
        # ras.reduce()
        ras.rescale(lambda x: 10 * x)

        ras.write(os.path.join(str(tmpdir), 'test'), format='GTiff', compress_tif=True)
        with pytest.raises(RuntimeError):
            ras.write(os.path.join(str(tmpdir), 'test.tif'), format='GTiff')

    # test writing a raster subset file with no data in memory
    outname = os.path.join(str(tmpdir), 'test_sub.tif')
    with Raster(testdata['tif']) as ras:
        ras.write(outname, format='GTiff', dim=(0, 0, 100, 100))
    with Raster(outname) as ras:
        assert ras.cols == 100
        assert ras.rows == 100


def test_Raster_extract(testdata):
    with Raster(testdata['tif']) as ras:
        assert ras.extract(px=624000, py=4830000, radius=5) == -10.48837461270875
        with pytest.raises(RuntimeError):
            ras.extract(1, 4830000)
        with pytest.raises(RuntimeError):
            ras.extract(624000, 1)

        # ensure corner extraction capability
        assert ras.extract(px=ras.geo['xmin'], py=ras.geo['ymax']) == -10.147890090942383
        assert ras.extract(px=ras.geo['xmin'], py=ras.geo['ymin']) == -14.640368461608887
        assert ras.extract(px=ras.geo['xmax'], py=ras.geo['ymax']) == -9.599242210388182
        assert ras.extract(px=ras.geo['xmax'], py=ras.geo['ymin']) == -9.406558990478516

        # test nodata handling capability and correct indexing
        mat = ras.matrix()
        mat[0:10, 0:10] = ras.nodata
        mat[207:217, 258:268] = ras.nodata
        ras.assign(mat, index=0)
        assert ras.extract(px=ras.geo['xmin'], py=ras.geo['ymax'], radius=5) == ras.nodata
        assert ras.extract(px=ras.geo['xmax'], py=ras.geo['ymin'], radius=5) == ras.nodata


def test_dtypes():
    assert dtypes('Float32') == 6
    with pytest.raises(ValueError):
        dtypes('foobar')


def test_stack(tmpdir, testdata):
    name = testdata['tif']
    outname = os.path.join(str(tmpdir), 'test')
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


def test_auxil(tmpdir, testdata):
    dir = str(tmpdir)
    ras = Raster(testdata['tif'])
    bbox = os.path.join(dir, 'bbox.shp')
    ras.bbox(bbox)
    ogr2ogr(bbox, os.path.join(dir, 'bbox.gml'), {'format': 'GML'})
    gdal_translate(ras.raster, os.path.join(dir, 'test'), {'format': 'ENVI'})
    gdal_rasterize(bbox, os.path.join(dir, 'test2'), {'format': 'GTiff', 'xRes': 20, 'yRes': 20})
