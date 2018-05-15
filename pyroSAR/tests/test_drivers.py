import pyroSAR
from pyroSAR.spatial import crsConvert, haversine
from pyroSAR.ancillary import finder
import pytest
import shutil
import sys
import os
from osgeo import ogr

testdir = os.getenv('TESTDATA_DIR', 'pyroSAR/tests/data/')

testcases = [
    #SAFE
    {'path': os.path.join('pyroSAR/tests/data/', 'S1A_IW_GRDH_1SDV_20150222T170750_20150222T170815_004739_005DD8_3768.zip'),
     'acquisition_mode': 'IW',
     'bbox_area': 7.573045244595988,
     'compression': 'zip',
     'corners': {'ymax': 52.183979, 'ymin': 50.295261, 'xmin': 8.017178, 'xmax': 12.0268},
     'hgt_len': 15,
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
        assert scene['pyro'].bbox().getArea() == scene['bbox_area']
        assert scene['pyro'].bbox().extent == scene['corners']
        assert scene['pyro'].bbox().nlayers == 1
        assert scene['pyro'].bbox().getProjection('epsg') == 4326
        assert scene['pyro'].bbox().proj4 == '+proj=longlat +datum=WGS84 +no_defs'
        assert isinstance(scene['pyro'].bbox().getFeatureByIndex(0), ogr.Feature)
        with pytest.raises(IndexError):
            scene['pyro'].bbox().getFeatureByIndex(1)
        bbox = scene['pyro'].bbox()
        bbox.reproject(32633)
        assert bbox.proj4 == '+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs'
        assert len(scene['pyro'].getHGT()) == scene['hgt_len']


def test_identify_fail():
    with pytest.raises(IOError):
        pyroSAR.identify(os.path.join(testdir, 'foobar'))


def test_export2dict():
    pass


def test_scene():
    scene = 'pyroSAR/tests/data/S1A_IW_GRDH_1SDV_20150222T170750_20150222T170815_004739_005DD8_3768.zip'
    dbfile = os.path.join('pyroSAR/tests/data/', 'scenes.db')
    with pyroSAR.Archive(dbfile) as db:
        db.insert(scene, verbose=True)
        assert db.size == (1, 0)
    id = pyroSAR.identify(scene)
    test_dir = 'pyroSAR/tests/data/test'
    os.makedirs(test_dir)
    id.bbox(outname='pyroSAR/tests/data/test/bbox_test.shp')
    assert id.is_processed(test_dir) is False
    id.unpack('pyroSAR/tests/data/test')
    assert id.compression is None
    os.remove(dbfile)
    id.export2sqlite(dbfile)
    with pytest.raises(IOError):
        id.getGammaImages()
    assert id.getGammaImages(id.scene) == []
    osvdir = os.path.join(id.scene, 'osv')
    if sys.version_info >= (2, 7, 9):
        id.getOSV(osvdir)

        with pyroSAR.OSV(osvdir) as osv:
            with pytest.raises(IOError):
                osv.catch(osvtype='XYZ')
            res = osv.catch(osvtype='RES', start=osv.mindate('POE'), stop=osv.maxdate('POE'))
            osv.retrieve(res)

            assert len(osv.getLocals('POE')) == 3
            assert len(osv.getLocals('RES')) == 21
            assert osv.match(id.start, 'POE') is not None
            assert osv.match(id.start, 'RES') is None
            osv.clean_res()
    else:
        with pytest.raises(RuntimeError):
            id.getOSV(osvdir, osvType='POE')

    shutil.rmtree(test_dir)
    os.remove(dbfile)


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
