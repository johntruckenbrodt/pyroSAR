import pyroSAR
import pytest
import tarfile as tf
import sys
import os
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


testdata = 'pyroSAR/tests/data/'
testfile1 = os.path.join(testdata, 'S1A_IW_GRDH_1SDV_20150222T170750_20150222T170815_004739_005DD8_3768.zip')
testfile2 = os.path.join(testdata, 'S1A__IW___A_20150309T173017_VV_grd_mli_geo_norm_db.tif')


def test_identify_fail():
    with pytest.raises(OSError):
        pyroSAR.identify(os.path.join(testdir, 'foobar'))
    with pytest.raises(RuntimeError):
        pyroSAR.identify(testfile2)


def test_identify_many_fail():
    assert pyroSAR.identify_many([testfile2]) == []


def test_filter_processed():
    scene = pyroSAR.identify(testfile1)
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


def test_getFileObj(tmpdir):
    scene = pyroSAR.identify(testfile1)
    scene.unpack(str(tmpdir))
    scene = pyroSAR.identify(scene.scene)
    item = scene.findfiles('manifest.safe')[0]
    assert os.path.basename(item) == 'manifest.safe'
    assert isinstance(scene.getFileObj(item).read(), (bytes, str))

    filename = os.path.join(str(tmpdir), os.path.basename(testfile1.replace('zip', 'tar.gz')))
    with tf.open(filename, 'w:gz') as tar:
        tar.add(scene.scene, arcname=os.path.basename(scene.scene))
    scene = pyroSAR.identify(filename)
    assert scene.compression == 'tar'
    item = scene.findfiles('manifest.safe')[0]
    assert isinstance(scene.getFileObj(item).read(), (bytes, str))
    with pytest.raises(IOError):
        pyroSAR.getFileObj('foo', 'bar')


def test_scene(tmpdir):
    dbfile = os.path.join(str(tmpdir), 'scenes.db')
    id = pyroSAR.identify(testfile1)
    assert isinstance(id.export2dict(), dict)
    with pytest.raises(RuntimeError):
        assert isinstance(id.gdalinfo(), dict)
    id.summary()
    id.bbox(outname=os.path.join(str(tmpdir), 'bbox_test.shp'), overwrite=True)
    assert id.is_processed(str(tmpdir)) is False
    id.unpack(str(tmpdir), overwrite=True)
    assert id.compression is None
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


def test_archive(tmpdir):
    id = pyroSAR.identify(testfile1)
    dbfile = os.path.join(str(tmpdir), 'scenes.db')
    db = pyroSAR.Archive(dbfile)
    db.insert(testfile1, verbose=False)
    assert db.is_registered(testfile1) is True
    assert len(db.get_unique_directories()) == 1
    assert db.select_duplicates() == []
    assert db.select_duplicates(outname_base='S1A__IW___A_20150222T170750', scene='scene.zip') == []
    assert len(db.select(mindate='20141001T192312', maxdate='20201001T192312')) == 1
    assert len(db.select(polarizations=['VV'])) == 1
    assert len(db.select(vectorobject=id.bbox())) == 1
    with pytest.raises(IOError):
        db.filter_scenelist([1])
    db.close()
    # separately test the with statement
    with pyroSAR.Archive(dbfile) as db:
        assert db.size == (1, 0)
    os.remove(dbfile)
    dbfile_old = os.path.join(testdata, 'archive_outdated.csv')
    with pytest.raises(OSError):
        with pyroSAR.Archive(dbfile) as db:
            db.import_outdated(dbfile_old)
