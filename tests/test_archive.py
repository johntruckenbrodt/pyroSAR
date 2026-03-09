import os

from sqlalchemy import Table, MetaData, Column, Integer, String
from geoalchemy2 import Geometry

from pyroSAR.drivers import identify
from pyroSAR.archive import Archive, drop_archive

from spatialist.vector import Vector

from shapely import wkt

import pytest

metadata = MetaData()

mytable = Table('mytable', metadata,
                Column('mytable_id', Integer, primary_key=True),
                Column('value', String(50)),
                Column('shape', Geometry(geometry_type='POLYGON',
                                         management=True, srid=4326)))


def test_archive(tmpdir, testdata):
    id = identify(testdata['s1'])
    dbfile = os.path.join(str(tmpdir), 'scenes.db')
    db = Archive(dbfile)
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
    out = db.select(vv=1, acquisition_mode=['IW', 'EW'])
    assert len(out) == 1
    assert isinstance(out[0], str)
    
    out = db.select(vv=1, return_value=['mindate', 'geometry_wkt', 'geometry_wkb'])
    assert len(out) == 1
    assert isinstance(out[0], tuple)
    assert out[0][0] == '20150222T170750'
    geom = wkt.loads('POLYGON(('
                     '8.505644 50.295261, 12.0268 50.688881, '
                     '11.653832 52.183979, 8.017178 51.788181, '
                     '8.505644 50.295261))')
    assert wkt.loads(out[0][1]) == geom
    assert out[0][2] == geom.wkb
    
    with pytest.raises(ValueError):
        out = db.select(vv=1, return_value=['foobar'])
    
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
    with Archive(dbfile) as db:
        db.insert(testdata['s1'])
        assert db.size == (1, 0)
        shp = os.path.join(str(tmpdir), 'db.shp')
        db.export2shp(shp)
    
    os.remove(dbfile)
    assert not os.path.isfile(dbfile)
    assert Vector(shp).nfeatures == 1
    
    with Archive(dbfile) as db:
        with pytest.raises(OSError):
            db.import_outdated(testdata['archive_old_csv'])
        with pytest.raises(RuntimeError):
            db.import_outdated('foobar')
    
    # the archive_old_bbox database contains a relative file name for the scene
    # so that it can be reimported into the new database. The working directory
    # is changed temporarily so that the scene can be found.
    cwd = os.getcwd()
    folder = os.path.dirname(os.path.realpath(__file__))
    os.chdir(os.path.join(folder, 'data'))
    with Archive(dbfile) as db:
        with Archive(testdata['archive_old_bbox'], legacy=True) as db_old:
            db.import_outdated(db_old)
    os.chdir(cwd)
    
    with pytest.raises(RuntimeError):
        db = Archive(testdata['archive_old_csv'])
    with pytest.raises(RuntimeError):
        db = Archive(testdata['archive_old_bbox'])


def test_archive_postgres(tmpdir, testdata):
    pguser = os.environ.get('PGUSER')
    pgpassword = os.environ.get('PGPASSWORD')
    pgport = os.environ.get('PGPORT')
    if pgport is not None:
        pgport = int(pgport)
    else:
        pgport = 5432
    
    id = identify(testdata['s1'])
    db = Archive('test', postgres=True, port=pgport, user=pguser, password=pgpassword)
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
    
    out = db.select(vv=1, return_value=['scene', 'start'])
    assert len(out) == 1
    assert isinstance(out[0], tuple)
    assert out[0][1] == '20150222T170750'
    
    with pytest.raises(ValueError):
        out = db.select(vv=1, return_value=['foobar'])
    
    db.add_tables(mytable)
    assert 'mytable' in db.get_tablenames()
    with pytest.raises(TypeError):
        db.filter_scenelist([1])
    db.close()
    with Archive('test', postgres=True, port=pgport,
                 user=pguser, password=pgpassword) as db:
        assert db.size == (1, 0)
        shp = os.path.join(str(tmpdir), 'db.shp')
        db.export2shp(shp)
        drop_archive(db)
    assert Vector(shp).nfeatures == 1
    
    with Archive('test', postgres=True, port=pgport,
                 user=pguser, password=pgpassword) as db:
        with pytest.raises(OSError):
            db.import_outdated(testdata['archive_old_csv'])
        drop_archive(db)
    
    # the archive_old_bbox database contains a relative file name for the scene
    # so that it can be reimported into the new database. The working directory
    # is changed temporarily so that the scene can be found.
    cwd = os.getcwd()
    folder = os.path.dirname(os.path.realpath(__file__))
    os.chdir(os.path.join(folder, 'data'))
    with Archive('test', postgres=True, port=pgport,
                 user=pguser, password=pgpassword) as db:
        with Archive(testdata['archive_old_bbox'], legacy=True) as db_old:
            db.import_outdated(db_old)
        drop_archive(db)
    os.chdir(cwd)
    
    dbfile = os.path.join(str(tmpdir), 'scenes.db')
    with Archive('test', postgres=True, port=pgport,
                 user=pguser, password=pgpassword) as db:
        with Archive(dbfile, legacy=True) as db_sqlite:
            db.import_outdated(db_sqlite)
        drop_archive(db)
    
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        Archive('test', postgres=True, user='hello_world', port=7080)
    assert pytest_wrapped_e.type == SystemExit
