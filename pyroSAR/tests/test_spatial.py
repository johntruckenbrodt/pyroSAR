
import pytest
from osgeo import ogr
from pyroSAR import identify
from pyroSAR.spatial import crsConvert, haversine


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
