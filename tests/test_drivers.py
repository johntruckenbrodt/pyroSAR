import pyroSAR

import pytest
import platform
import tarfile as tf
import os
from datetime import datetime
from spatialist import Vector
from sqlalchemy import Table, MetaData, Column, Integer, String
from geoalchemy2 import Geometry

metadata = MetaData()

mytable = Table('mytable', metadata,
                Column('mytable_id', Integer, primary_key=True),
                Column('value', String(50)),
                Column('shape', Geometry('POLYGON', management=True, srid=4326)))


@pytest.fixture()
def testcases():
    cases = {
        's1': {
            'acquisition_mode': 'IW',
            'bbox_area': 7.573045244595988,
            'compression': 'zip',
            'corners': {'ymax': 52.183979, 'ymin': 50.295261, 'xmin': 8.017178, 'xmax': 12.0268},
            'hgt_len': 15,
            'lines': 16685,
            'orbit': 'A',
            'outname': 'S1A__IW___A_20150222T170750',
            'polarizations': ['VV', 'VH'],
            'product': 'GRD',
            'samples': 25368,
            'sensor': 'S1A',
            'spacing': (10.0, 9.998647),
            'start': '20150222T170750',
            'stop': '20150222T170815'
        },
        'psr2': {
            'acquisition_mode': 'FBD',
            'compression': 'zip',
            'corners': {'xmin': -62.9005207, 'xmax': -62.1629744, 'ymin': -11.4233051, 'ymax': -10.6783401},
            'hgt_len': 2,
            'lines': 13160,
            'orbit': 'A',
            'outname': 'PSR2_FBD__A_20140909T043342',
            'polarizations': ['HH', 'HV'],
            'product': '1.5',
            'samples': 12870,
            'sensor': 'PSR2',
            'spacing': (6.25, 6.25),
            'start': '20140909T043342',
            'stop': '20140909T043352'
        }
    }
    return cases


@pytest.fixture
def scene(testcases, testdata, request):
    case = testcases[request.param]
    case['pyro'] = pyroSAR.identify(testdata[request.param])
    return case


class Test_Metadata():
    @pytest.mark.parametrize('scene', ['s1', 'psr2'], indirect=True)
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


def test_identify_fail(testdir, testdata):
    with pytest.raises(OSError):
        pyroSAR.identify(os.path.join(testdir, 'foobar'))
    with pytest.raises(RuntimeError):
        pyroSAR.identify(testdata['tif'])


def test_identify_many_fail(testdata):
    assert pyroSAR.identify_many([testdata['tif']]) == []


def test_filter_processed(tmpdir, testdata):
    scene = pyroSAR.identify(testdata['s1'])
    assert len(pyroSAR.filter_processed([scene], str(tmpdir))) == 1


def test_parse_date():
    with pytest.raises(ValueError):
        print(pyroSAR.parse_date(1))
    with pytest.raises(ValueError):
        print(pyroSAR.parse_date('foobar'))
    assert pyroSAR.parse_date(datetime(2006, 11, 21)) == '20061121T000000'


def test_export2dict():
    pass


def test_getFileObj(tmpdir, testdata):
    scene = pyroSAR.identify(testdata['s1'])
    if platform.system() == 'Windows':
        directory = u'\\\\?\\' + str(tmpdir)
    else:
        directory = str(tmpdir)
    scene.unpack(directory)
    scene = pyroSAR.identify(scene.scene)
    item = scene.findfiles('manifest.safe')[0]
    assert os.path.basename(item) == 'manifest.safe'
    assert isinstance(scene.getFileObj(item).read(), (bytes, str))
    
    filename = os.path.join(str(tmpdir), os.path.basename(testdata['s1'].replace('zip', 'tar.gz')))
    with tf.open(filename, 'w:gz') as tar:
        tar.add(scene.scene, arcname=os.path.basename(scene.scene))
    # test error if scene is not a directory, zip or tar
    with pytest.raises(RuntimeError):
        pyroSAR.getFileObj(scene=os.path.join(scene.scene, 'manifest.safe'), filename='bar')
    scene = pyroSAR.identify(filename)
    assert scene.compression == 'tar'
    item = scene.findfiles('manifest.safe')[0]
    assert isinstance(scene.getFileObj(item).read(), (bytes, str))
    with pytest.raises(RuntimeError):
        pyroSAR.getFileObj('foo', 'bar')


def test_scene(tmpdir, testdata):
    dbfile = os.path.join(str(tmpdir), 'scenes.db')
    id = pyroSAR.identify(testdata['s1'])
    assert isinstance(id.export2dict(), dict)
    with pytest.raises(RuntimeError):
        assert isinstance(id.gdalinfo(), dict)
    id.summary()
    id.bbox(outname=os.path.join(str(tmpdir), 'bbox_test.shp'), overwrite=True)
    assert id.is_processed(str(tmpdir)) is False
    id.unpack(str(tmpdir), overwrite=True)
    assert id.compression is None
    id.export2sqlite(dbfile)
    with pytest.raises(RuntimeError):
        id.getGammaImages()
    assert id.getGammaImages(id.scene) == []
    id = pyroSAR.identify(testdata['psr2'])
    assert id.getCorners() == {'xmax': -62.1629744, 'xmin': -62.9005207,
                               'ymax': -10.6783401, 'ymin': -11.4233051}


def test_archive(tmpdir, testdata):
    id = pyroSAR.identify(testdata['s1'])
    dbfile = os.path.join(str(tmpdir), 'scenes.db')
    db = pyroSAR.Archive(dbfile)
    db.insert(testdata['s1'])
    assert all(isinstance(x, str) for x in db.get_tablenames())
    assert all(isinstance(x, str) for x in db.get_colnames())
    assert db.is_registered(testdata['s1']) is True
    assert len(db.get_unique_directories()) == 1
    assert db.select_duplicates() == []
    assert db.select_duplicates(outname_base='S1A__IW___A_20150222T170750', scene='scene.zip') == []
    assert len(db.select(mindate='20141001T192312', maxdate='20201001T192312')) == 1
    assert len(db.select(polarizations=['VV'])) == 1
    assert len(db.select(vectorobject=id.bbox())) == 1
    assert len(db.select(sensor='S1A', vectorobject='foo', processdir=str(tmpdir))) == 1
    assert len(db.select(sensor='S1A', mindate='foo', maxdate='bar', foobar='foobar')) == 1
    out = db.select(vv=1, acquisition_mode=('IW', 'EW'))
    assert len(out) == 1
    assert isinstance(out[0], str)
    
    db.insert(testdata['s1_3'])
    db.insert(testdata['s1_4'])
    db.drop_element(testdata['s1_3'])
    assert db.size == (2, 0)
    db.drop_element(testdata['s1_4'])
    
    db.add_tables(mytable)
    assert 'mytable' in db.get_tablenames()
    with pytest.raises(TypeError):
        db.filter_scenelist([1])
    db.close()


def test_archive2(tmpdir, testdata):
    dbfile = os.path.join(str(tmpdir), 'scenes.db')
    with pyroSAR.Archive(dbfile) as db:
        db.insert(testdata['s1'])
        assert db.size == (1, 0)
        shp = os.path.join(str(tmpdir), 'db.shp')
        db.export2shp(shp)
    
    os.remove(dbfile)
    assert not os.path.isfile(dbfile)
    assert Vector(shp).nfeatures == 1
    
    with pytest.raises(OSError):
        with pyroSAR.Archive(dbfile) as db:
            db.import_outdated(testdata['archive_old_csv'])
    
    # the archive_old_bbox database contains a relative file name for the scene
    # so that it can be reimported into the new database. The working directory
    # is changed temporarily so that the scene can be found.
    cwd = os.getcwd()
    folder = os.path.dirname(os.path.realpath(__file__))
    os.chdir(os.path.join(folder, 'data'))
    with pyroSAR.Archive(dbfile) as db:
        with pyroSAR.Archive(testdata['archive_old_bbox'], legacy=True) as db_old:
            db.import_outdated(db_old)
    os.chdir(cwd)
    
    with pytest.raises(RuntimeError):
        db = pyroSAR.Archive(testdata['archive_old_csv'])
    with pytest.raises(RuntimeError):
        db = pyroSAR.Archive(testdata['archive_old_bbox'])


def test_archive_postgres(tmpdir, testdata):
    pguser = os.environ.get('PGUSER')
    pgpassword = os.environ.get('PGPASSWORD')
    pgport = os.environ.get('PGPORT')
    if pgport is not None:
        pgport = int(pgport)
    else:
        pgport = 5432
    
    id = pyroSAR.identify(testdata['s1'])
    db = pyroSAR.Archive('test', postgres=True, port=pgport, user=pguser, password=pgpassword)
    db.insert(testdata['s1'])
    assert all(isinstance(x, str) for x in db.get_tablenames())
    assert all(isinstance(x, str) for x in db.get_colnames())
    assert db.is_registered(testdata['s1']) is True
    assert len(db.get_unique_directories()) == 1
    assert db.select_duplicates() == []
    assert db.select_duplicates(outname_base='S1A__IW___A_20150222T170750', scene='scene.zip') == []
    assert len(db.select(mindate='20141001T192312', maxdate='20201001T192312')) == 1
    assert len(db.select(polarizations=['VV'])) == 1
    assert len(db.select(vectorobject=id.bbox())) == 1
    assert len(db.select(sensor='S1A', vectorobject='foo', processdir=str(tmpdir))) == 1
    assert len(db.select(sensor='S1A', mindate='foo', maxdate='bar', foobar='foobar')) == 1
    out = db.select(vv=1, acquisition_mode=('IW', 'EW'))
    assert len(out) == 1
    assert isinstance(out[0], str)
    db.add_tables(mytable)
    assert 'mytable' in db.get_tablenames()
    with pytest.raises(TypeError):
        db.filter_scenelist([1])
    db.close()
    with pyroSAR.Archive('test', postgres=True, port=pgport,
                         user=pguser, password=pgpassword) as db:
        assert db.size == (1, 0)
        shp = os.path.join(str(tmpdir), 'db.shp')
        db.export2shp(shp)
        pyroSAR.drop_archive(db)
    assert Vector(shp).nfeatures == 1
    
    with pyroSAR.Archive('test', postgres=True, port=pgport,
                         user=pguser, password=pgpassword) as db:
        with pytest.raises(OSError):
            db.import_outdated(testdata['archive_old_csv'])
        pyroSAR.drop_archive(db)
    
    # the archive_old_bbox database contains a relative file name for the scene
    # so that it can be reimported into the new database. The working directory
    # is changed temporarily so that the scene can be found.
    cwd = os.getcwd()
    folder = os.path.dirname(os.path.realpath(__file__))
    os.chdir(os.path.join(folder, 'data'))
    with pyroSAR.Archive('test', postgres=True, port=pgport,
                         user=pguser, password=pgpassword) as db:
        with pyroSAR.Archive(testdata['archive_old_bbox'], legacy=True) as db_old:
            db.import_outdated(db_old)
        pyroSAR.drop_archive(db)
    os.chdir(cwd)
    
    dbfile = os.path.join(str(tmpdir), 'scenes.db')
    with pyroSAR.Archive('test', postgres=True, port=pgport,
                         user=pguser, password=pgpassword) as db:
        with pyroSAR.Archive(dbfile, legacy=True) as db_sqlite:
            db.import_outdated(db_sqlite)
        pyroSAR.drop_archive(db)
    
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        pyroSAR.Archive('test', postgres=True, user='hello_world', port=7080)
    assert pytest_wrapped_e.type == SystemExit


datasets = ['asar', 'ers1_esa', 'ers1_ceos', 'psr2', 's1']


@pytest.mark.parametrize('dataset', datasets)
def test_geometry(testdata, dataset):
    scene = pyroSAR.identify(testdata[dataset])
    with scene.geometry() as geom:
        assert isinstance(geom, Vector)
