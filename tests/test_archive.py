import pytest

import sys
import textwrap
import platform
from uuid import uuid4
import subprocess as sp
from pathlib import Path
from dataclasses import dataclass

from geoalchemy2 import Geometry
from shapely import wkt
from sqlalchemy import Table, MetaData, Column, Integer, String

from pyroSAR.archive import Archive, drop_archive
from pyroSAR.drivers import identify
from spatialist.vector import Vector


@dataclass(frozen=True)
class ArchiveConfig:
    driver: str
    dbfile: str
    kwargs: dict
    
    def open(self, **kwargs):
        options = {**self.kwargs, **kwargs}
        return Archive(self.dbfile, **options)


@pytest.fixture
def sqlite_archive_config(tmp_path):
    return ArchiveConfig(
        driver='sqlite',
        dbfile=str(tmp_path / 'scenes.db'),
        kwargs={},
    )


@pytest.fixture
def postgres_archive_config(pg_conn):
    config = ArchiveConfig(
        driver='postgresql',
        dbfile=f'pytest_{uuid4().hex}',
        kwargs={
            'postgres': True,
            'host': pg_conn.info.host,
            'port': pg_conn.info.port,
            'user': pg_conn.info.user,
            'password': pg_conn.info.password,
        },
    )
    
    yield config
    
    # Reopen the database so that it can also be cleaned up if the test
    # itself has already closed its Archive instance.
    db = config.open(cleanup=False)
    drop_archive(db)


@pytest.fixture(
    params=['sqlite_archive_config', 'postgres_archive_config'],
    ids=['sqlite', 'postgresql'],
)
def archive_config(request):
    """
    Run backend-independent tests against both SQLite/SpatiaLite
    and PostgreSQL/PostGIS.
    """
    return request.getfixturevalue(request.param)


@pytest.fixture
def archive(archive_config, testdata):
    """
    Create an Archive containing one Sentinel-1 scene.
    """
    with archive_config.open() as db:
        db.insert(testdata['s1'])
        yield db


@pytest.fixture
def extra_table():
    """
    Create a fresh SQLAlchemy table for each test.

    Archive.add_tables() modifies the table metadata, so sharing the same
    Table instance between tests should be avoided.
    """
    metadata = MetaData()
    
    return Table(
        'mytable',
        metadata,
        Column('mytable_id', Integer, primary_key=True),
        Column('value', String(50)),
        Column(
            'shape',
            Geometry(
                geometry_type='POLYGON',
                srid=4326,
            ),
        ),
    )


@pytest.fixture
def extra_tables(extra_table):
    metadata = extra_table.metadata
    
    table2 = Table(
        'mytable2',
        metadata,
        Column('mytable2_id', Integer, primary_key=True),
        Column('value', String(50)),
        Column(
            'shape',
            Geometry(
                geometry_type='POLYGON',
                srid=4326,
            ),
        ),
    )
    
    return [extra_table, table2]


@pytest.fixture(
    params=[
        ('extra_table', ['mytable']),
        ('extra_tables', ['mytable', 'mytable2']),
    ],
    ids=['single-table', 'two-tables'],
)
def extra_table_case(request):
    fixture_name, expected = request.param
    tables = request.getfixturevalue(fixture_name)
    return tables, expected


def test_schema(archive):
    assert all(isinstance(x, str) for x in archive.get_tablenames())
    assert all(isinstance(x, str) for x in archive.get_colnames())


def test_registration(archive, testdata):
    assert archive.is_registered(testdata['s1'])
    assert len(archive.get_unique_directories()) == 1


def test_select_duplicates(archive):
    assert archive.select_duplicates() == []
    
    assert archive.select_duplicates(
        outname_base='S1A__IW___A_20150222T170750',
        scene='scene.zip',
    ) == []


@pytest.mark.parametrize(
    'kwargs',
    [
        {
            'mindate': '20141001T192312',
            'maxdate': '20201001T192312',
        },
        {
            'polarizations': ['VV'],
        },
        {
            'vv': 1,
            'acquisition_mode': ['IW', 'EW'],
        },
    ],
    ids=[
        'date',
        'polarization',
        'acquisition-mode',
    ],
)
def test_select(archive, kwargs):
    out = archive.select(**kwargs)
    
    assert len(out) == 1
    assert isinstance(out[0], str)


def test_select_vectorobject(archive, testdata):
    scene = identify(testdata['s1'])
    
    out = archive.select(vectorobject=scene.bbox())
    
    assert len(out) == 1


def test_select_processdir(archive, tmp_path):
    out = archive.select(
        sensor='S1A',
        vectorobject='foo',
        processdir=str(tmp_path),
    )
    
    assert len(out) == 1


def test_select_invalid_filters(archive):
    out = archive.select(
        sensor='S1A',
        mindate='foo',
        maxdate='bar',
        foobar='foobar',
    )
    
    assert len(out) == 1


def test_select_return_values(archive):
    out = archive.select(
        vv=1,
        return_value=[
            'mindate',
            'geometry_wkt',
            'geometry_wkb',
        ],
    )
    
    assert len(out) == 1
    assert isinstance(out[0], tuple)
    assert out[0][0] == '20150222T170750'
    
    geom = wkt.loads(
        'POLYGON(('
        '8.505644 50.295261, '
        '12.0268 50.688881, '
        '11.653832 52.183979, '
        '8.017178 51.788181, '
        '8.505644 50.295261'
        '))'
    )
    
    assert wkt.loads(out[0][1]) == geom
    assert out[0][2] == geom.wkb


def test_select_invalid_return_value(archive):
    with pytest.raises(ValueError):
        archive.select(
            vv=1,
            return_value=['foobar'],
        )


def test_get_tablenames(archive):
    assert set(archive.get_tablenames()) == {
        'data',
        'duplicates',
    }

    assert {
        '_managed_tables',
        'data',
        'duplicates',
    } <= set(
        archive.get_tablenames(return_all=True)
    )


def test_managed_tables_persist(
        archive_config,
        extra_table,
):
    with archive_config.open() as db:
        db.add_tables(extra_table)

        assert 'mytable' in db.get_tablenames()

    with archive_config.open() as db:
        assert 'mytable' in db.get_tablenames()


def test_drop_element(archive, testdata):
    archive.insert(testdata['s1_3'])
    archive.insert(testdata['s1_4'])
    
    archive.drop_element(testdata['s1_3'])
    
    assert archive.size == (2, 0)
    
    archive.drop_element(testdata['s1_4'])


def test_add_tables(archive, extra_table_case):
    tables, expected = extra_table_case

    archive.add_tables(tables)

    tablenames = archive.get_tablenames()

    assert {'data', 'duplicates'} <= set(tablenames)
    assert set(expected) <= set(tablenames)


def test_drop_table(
        archive,
        extra_table,
):
    archive.add_tables(extra_table)

    assert 'mytable' in archive.get_tablenames()

    archive.drop_table('mytable')

    assert 'mytable' not in archive.get_tablenames()
    assert 'mytable' not in archive.get_tablenames(
        return_all=True
    )


def test_drop_managed_tables_registry(archive):
    with pytest.raises(ValueError):
        archive.drop_table('_managed_tables')


def test_filter_scenelist_invalid_input(archive):
    with pytest.raises(TypeError):
        archive.filter_scenelist([1])


def test_persistence_and_export(
        archive_config,
        testdata,
        tmp_path,
):
    with archive_config.open() as db:
        db.insert(testdata['s1'])
    
    with archive_config.open() as db:
        assert db.size == (1, 0)
        
        shp = tmp_path / 'db.shp'
        db.export2shp(str(shp))
    
    assert Vector(str(shp)).nfeatures == 1


def test_import_outdated_csv(
        archive_config,
        testdata,
):
    with archive_config.open() as db:
        with pytest.raises(OSError):
            db.import_outdated(
                testdata['archive_old_csv']
            )


def test_import_outdated_invalid_input(
        archive_config,
):
    with archive_config.open() as db:
        with pytest.raises(RuntimeError):
            db.import_outdated('foobar')


def test_import_outdated_archive(
        archive_config,
        testdata,
        monkeypatch,
):
    """
    Test importing an outdated Archive instance.

    The archive_old_bbox database contains a relative scene filename.
    Temporarily change the working directory so the scene can be found.
    """
    folder = Path(__file__).parent / 'data'
    monkeypatch.chdir(folder)
    
    with archive_config.open() as db:
        with Archive(
                testdata['archive_old_bbox'],
                legacy=True,
        ) as db_old:
            db.import_outdated(db_old)


@pytest.mark.parametrize(
    'key',
    [
        'archive_old_csv',
        'archive_old_bbox',
    ],
)
def test_open_outdated_without_legacy(
        testdata,
        key,
):
    with pytest.raises(RuntimeError):
        Archive(testdata[key])


def test_sqlite_deleted_after_close(
        sqlite_archive_config,
        testdata,
):
    """
    Check that closing a SQLite archive releases the database file.
    """
    with sqlite_archive_config.open() as db:
        db.insert(testdata['s1'])
    
    dbfile = Path(sqlite_archive_config.dbfile)
    
    if platform.system() == 'Windows':
        with pytest.raises(PermissionError):
            dbfile.unlink()
    else:
        dbfile.unlink()
        assert not dbfile.exists()


def test_postgres_import_from_sqlite(
        postgres_archive_config,
        tmp_path,
):
    """
    Check that a SQLite Archive can be passed to import_outdated()
    when the target archive uses PostgreSQL.
    """
    dbfile = tmp_path / 'scenes.db'
    
    with postgres_archive_config.open() as db:
        with Archive(
                str(dbfile),
                legacy=True,
        ) as db_sqlite:
            db.import_outdated(db_sqlite)


def test_postgres_invalid_connection():
    with pytest.raises(
            SystemExit,
            match='Server not found!',
    ):
        Archive(
            'test',
            postgres=True,
            user='hello_world',
            port=7080,
        )


def test_sqlite_close_does_not_break_lxml():
    """
    Confirm the current erroneous behavior of SpatiaLite on Windows.
    """
    code = textwrap.dedent("""
        import tempfile
        from pathlib import Path

        from lxml import html
        from pyroSAR import Archive

        with tempfile.TemporaryDirectory() as tmpdir:
            dbfile = Path(tmpdir) / "archive.db"

            archive = Archive(str(dbfile))

            test = html.fromstring("<p>before</p>")

            archive.close()

            test = html.fromstring("<p>after</p>")
            print(test.text)
    """)
    
    result = sp.run(
        args=[
            sys.executable,
            '-X',
            'faulthandler',
            '-c',
            code,
        ],
        capture_output=True,
        text=True,
    )
    
    if platform.system() == 'Windows':
        assert result.returncode != 0
    else:
        assert result.returncode == 0
    
    assert result.stdout.strip() == 'after'
