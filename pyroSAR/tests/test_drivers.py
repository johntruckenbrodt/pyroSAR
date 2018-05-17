import pyroSAR
import pytest
import shutil
import tarfile as tf
import sys
import os
from osgeo import ogr
from datetime import datetime

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
        assert len(scene['pyro'].getHGT()) == scene['hgt_len']


def test_identify_fail():
    with pytest.raises(OSError):
        pyroSAR.identify(os.path.join(testdir, 'foobar'))
    with pytest.raises(RuntimeError):
        pyroSAR.identify('pyroSAR/tests/data/S1A__IW___A_20150309T173017_VV_grd_mli_geo_norm_db.tif')


def test_identify_many_fail():
    filename = 'pyroSAR/tests/data/S1A__IW___A_20150309T173017_VV_grd_mli_geo_norm_db.tif'
    assert pyroSAR.identify_many([filename]) == []


def test_filter_processed():
    filename = 'pyroSAR/tests/data/S1A_IW_GRDH_1SDV_20150222T170750_20150222T170815_004739_005DD8_3768.zip'
    scene = pyroSAR.identify(filename)
    outdir = 'pyroSAR/tests'
    assert len(pyroSAR.filter_processed([scene], outdir)) == 1


def test_parse_date():
    with pytest.raises(ValueError):
        print(pyroSAR.parse_date(1))
    with pytest.raises(ValueError):
        print(pyroSAR.parse_date('foobar'))
    assert pyroSAR.parse_date(datetime(2006, 11, 21)) == '20061121T000000'


def test_export2dict():
    pass


def test_getFileObj():
    filename = 'pyroSAR/tests/data/S1A_IW_GRDH_1SDV_20150222T170750_20150222T170815_004739_005DD8_3768.zip'
    scene = pyroSAR.identify(filename)
    tmpdir = 'pyroSAR/tests/data/tmp'
    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)
    os.makedirs(tmpdir)
    try:
        scene.unpack(tmpdir)
    except RuntimeError:
        pass
    item = scene.findfiles('manifest.safe')[0]
    assert os.path.basename(item) == 'manifest.safe'
    assert isinstance(scene.getFileObj(item).read(), (bytes, str))

    filename = os.path.join(tmpdir, os.path.basename(filename.replace('zip', 'tar.gz')))
    with tf.open(filename, 'w:gz') as tar:
        tar.add(scene.scene, arcname=os.path.basename(scene.scene))
    scene = pyroSAR.identify(filename)
    item = scene.findfiles('manifest.safe')[0]
    assert isinstance(scene.getFileObj(item).read(), (bytes, str))
    shutil.rmtree(tmpdir)
    with pytest.raises(IOError):
        pyroSAR.getFileObj('foo', 'bar')


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
