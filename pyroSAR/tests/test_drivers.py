import pyroSAR
from pyroSAR.spatial import crsConvert, haversine
import pytest
import os

testdir = os.getenv('TESTDATA_DIR', 'pyroSAR/tests/data/')

testcases = [
    #SAFE
    {'path': os.path.join('pyroSAR/tests/data/', 'S1A_IW_GRDH_1SDV_20150222T170750_20150222T170815_004739_005DD8_3768.zip'),
     'acquisition_mode': 'IW',
     'bbox_area': 7.573045244595988,
     'compression': 'zip',
     'corners': {'ymax': 52.183979, 'ymin': 50.295261, 'xmin': 8.017178, 'xmax': 12.0268},
     'lines': 16685,
     'outname': 'S1A__IW___A_20150222T170750',
     'orbit': 'A',
     'polarizations': ['VV', 'VH'],
     'product': 'GRD',
     'samples': 25368,
     'sensor': 'S1A',
     'spacing': (10.0, 9.998647),
     'start': '20150222T170750',
     'stop': '20150222T170815'
     }
]


@pytest.fixture
def scene(case):
    case['pyro'] = pyroSAR.identify(case['path'])
    return case


@pytest.mark.parametrize('case', testcases)
class Test_Metadata():
    def test_attributes(self, scene):
        assert scene['pyro'].acquisition_mode == scene['acquisition_mode']
        assert scene['pyro'].compression == scene['compression']
        assert scene['pyro'].getCorners() == scene['corners']
        assert scene['pyro'].lines == scene['lines']
        assert scene['pyro'].outname_base() == scene['outname']
        assert scene['pyro'].orbit == scene['orbit']
        assert scene['pyro'].polarizations == scene['polarizations']
        assert scene['pyro'].product == scene['product']
        assert scene['pyro'].samples == scene['samples']
        assert scene['pyro'].start == scene['start']
        assert scene['pyro'].stop == scene['stop']
        assert scene['pyro'].sensor == scene['sensor']
        assert scene['pyro'].spacing == scene['spacing']
        assert scene['pyro'].is_processed('data/') is False
        assert scene['pyro'].bbox().getArea() == scene['bbox_area']


def test_identify_fail():
    with pytest.raises(IOError):
        pyroSAR.identify(os.path.join(testdir, 'foobar'))


def test_export2dict():
    pass


def test_archive():
    scene = 'pyroSAR/tests/data/S1A_IW_GRDH_1SDV_20150222T170750_20150222T170815_004739_005DD8_3768.zip'
    dbfile = os.path.join('pyroSAR/tests/data/', 'scenes.db')
    if os.path.isfile(dbfile):
        os.remove(dbfile)
    with pyroSAR.Archive(dbfile) as db:
        db.insert(scene, verbose=True)
        assert db.size == (1, 0)


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
