###############################################################################
# Scene database tools for pyroSAR
# Copyright (c) 2016-2026, the pyroSAR Developers.

# This file is part of the pyroSAR Project. It is subject to the
# license terms in the LICENSE.txt file found in the top-level
# directory of this distribution and at
# https://github.com/johntruckenbrodt/pyroSAR/blob/master/LICENSE.txt.
# No part of the pyroSAR project, including this file, may be
# copied, modified, propagated, or distributed except according
# to the terms contained in the LICENSE.txt file.
###############################################################################
from __future__ import annotations

import os
import re
import gc
import shutil
import sys
import socket
import time
import platform
import logging
import csv
import inspect
from datetime import datetime

import progressbar as pb

from types import TracebackType
from typing import Any, Protocol, runtime_checkable, Literal

from osgeo import gdal

from spatialist import sqlite3
from spatialist.vector import Vector
from spatialist.ancillary import finder

from pyroSAR.drivers import identify, identify_many, ID

from sqlalchemy import create_engine, Table, MetaData, Column, Integer, String, exc
from sqlalchemy import inspect as sql_inspect
from sqlalchemy.event import listen
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql import select, func
from sqlalchemy.engine.url import URL
from sqlalchemy.ext.automap import automap_base
from sqlalchemy_utils import database_exists, create_database, drop_database
from geoalchemy2 import Geometry

log = logging.getLogger(__name__)

gdal.UseExceptions()


@runtime_checkable
class SceneArchive(Protocol):
    """
    Common interface for scene archive backends.

    Implementations may represent local databases, STAC catalogs, remote APIs,
    or other scene repositories, but should expose a consistent `select`
    method and support context-manager usage.
    """
    
    def __enter__(self) -> SceneArchive:
        """
        Enter the archive context.
        """
        ...
    
    def __exit__(
            self,
            exc_type: type[BaseException] | None,
            exc_val: BaseException | None,
            exc_tb: TracebackType | None,
    ) -> None:
        """
        Exit the archive context and release resources if necessary.
        """
        ...
    
    def close(self) -> None:
        """
        Release open resources.

        Implementations that do not hold resources may implement this as a no-op.
        """
        ...
    
    def select(
            self,
            sensor: str | list[str] | None = None,
            product: str | list[str] | None = None,
            acquisition_mode: str | list[str] | None = None,
            mindate: str | datetime | None = None,
            maxdate: str | datetime | None = None,
            vectorobject: Vector | None = None,
            date_strict: bool = True,
            return_value: str | list[str] = "scene",
            **kwargs: Any,
    ) -> list[Any]:
        """
        Select scenes matching the query parameters.

        Parameters
        ----------
        sensor:
            One sensor or a list of sensors.
        product:
            One product type or a list of product types.
        acquisition_mode:
            One acquisition mode or a list of acquisition modes.
        mindate:
            Minimum acquisition date/time.
        maxdate:
            Maximum acquisition date/time.
        vectorobject:
            Spatial search geometry.
        date_strict:
            Whether date filtering should be strict.
        return_value:
            One return field or a list of return fields.
        **kwargs:
            Backend-specific optional query arguments.

        Returns
        -------
            The query result. Implementations may return a list of scalar values or
            tuples depending on `return_value`.
        """
        ...


class Archive(SceneArchive):
    """
    Utility for storing SAR image metadata in a database

    Parameters
    ----------
    dbfile:
        the filename for the SpatiaLite database. This might either point to an
        existing database or will be created otherwise. If postgres is set to True,
        this will be the name for the PostgreSQL database.
    custom_fields:
        a dictionary containing additional non-standard database column names and data types;
        the names must be attributes of the SAR scenes to be inserted (i.e. id.attr) or keys in their meta attribute
        (i.e. id.meta['attr'])
    postgres:
        enable postgres driver for the database. Default: False
    user:
        required for postgres driver: username to access the database. Default: 'postgres'
    password:
        required for postgres driver: password to access the database. Default: '1234'
    host:
        required for postgres driver: host where the database is hosted. Default: 'localhost'
    port:
        required for postgres driver: port number to the database. Default: 5432
    cleanup:
        check whether all registered scenes exist and remove missing entries?
    legacy:
        open an outdated database in legacy mode to import into a new database.
        Opening an outdated database without legacy mode will throw a RuntimeError.

    Examples
    ----------
    Ingest all Sentinel-1 scenes in a directory and its subdirectories into the database:

    >>> from pyroSAR import Archive, identify
    >>> from spatialist.ancillary import finder
    >>> dbfile = '/.../scenelist.db'
    >>> archive_s1 = '/.../sentinel1/GRD'
    >>> scenes_s1 = finder(archive_s1, [r'^S1.*.zip'], regex=True, recursive=True)
    >>> with Archive(dbfile) as archive:
    >>>     archive.insert(scenes_s1)

    select all Sentinel-1 A/B scenes stored in the database, which

     * overlap with a test site
     * were acquired in Ground-Range-Detected (GRD) Interferometric Wide Swath (IW) mode before 2018
     * contain a VV polarization image
     * have not been processed to directory `outdir` before

    >>> from pyroSAR import Archive
    >>> from spatialist import Vector
    >>> archive = Archive('/.../scenelist.db')
    >>> site = Vector('/path/to/site.shp')
    >>> outdir = '/path/to/processed/results'
    >>> maxdate = '20171231T235959'
    >>> selection_proc = archive.select(vectorobject=site, processdir=outdir,
    >>>                                 maxdate=maxdate, sensor=['S1A', 'S1B'],
    >>>                                 product='GRD', acquisition_mode='IW', vv=1)
    >>> archive.close()

    Alternatively, the `with` statement can be used.
    In this case to just check whether one particular scene is already registered in the database:

    >>> from pyroSAR import identify, Archive
    >>> scene = identify('S1A_IW_SLC__1SDV_20150330T170734_20150330T170801_005264_006A6C_DA69.zip')
    >>> with Archive('/.../scenelist.db') as archive:
    >>>     print(archive.is_registered(scene.scene))

    When providing 'postgres' as driver, a PostgreSQL database will be created at a given host.
    Additional arguments are required.

    >>> from pyroSAR import Archive, identify
    >>> from spatialist.ancillary import finder
    >>> dbfile = 'scenelist_db'
    >>> archive_s1 = '/.../sentinel1/GRD'
    >>> scenes_s1 = finder(archive_s1, [r'^S1.*.zip'], regex=True, recursive=True)
    >>> with Archive(dbfile, driver='postgres', user='user', password='password', host='host', port=5432) as archive:
    >>>     archive.insert(scenes_s1)

    Importing an old database:

    >>> from pyroSAR import Archive
    >>> db_new = 'scenes.db'
    >>> db_old = 'scenes_old.db'
    >>> with Archive(db_new) as db:
    >>>     with Archive(db_old, legacy=True) as db_old:
    >>>         db.import_outdated(db_old)
    """
    
    def __init__(
            self,
            dbfile: str,
            custom_fields: dict[str, Any] | None = None,
            postgres: bool = False,
            user: str = 'postgres',
            password: str = '1234',
            host: str = 'localhost',
            port: int = 5432,
            cleanup: bool = True,
            legacy: bool = False
    ):
        
        if dbfile.endswith('.csv'):
            raise RuntimeError("Please create a new Archive database and import the"
                               "CSV file using db.import_outdated('<file>.csv').")
        # check for driver, if postgres then check if server is reachable
        if not postgres:
            self.driver = 'sqlite'
            dirname = os.path.dirname(os.path.abspath(dbfile))
            w_ok = os.access(dirname, os.W_OK)
            if not w_ok:
                raise RuntimeError('cannot write to directory {}'.format(dirname))
            # catch if .db extension is missing
            root, ext = os.path.splitext(dbfile)
            if len(ext) == 0:
                dbfile = root + '.db'
        else:
            self.driver = 'postgresql'
            if not self.__check_host(host, port):
                sys.exit('Server not found!')
        
        connect_args = {}
        
        # create dict, with which a URL to the db is created
        if self.driver == 'sqlite':
            self.url_dict = {'drivername': self.driver,
                             'database': dbfile,
                             'query': {'charset': 'utf8'}}
        if self.driver == 'postgresql':
            self.url_dict = {'drivername': self.driver,
                             'username': user,
                             'password': password,
                             'host': host,
                             'port': port,
                             'database': dbfile}
            connect_args = {
                'keepalives': 1,
                'keepalives_idle': 30,
                'keepalives_interval': 10,
                'keepalives_count': 5}
        
        # create engine, containing URL and driver
        log.debug('starting DB engine for {}'.format(URL.create(**self.url_dict)))
        self.url = URL.create(**self.url_dict)
        # https://www.postgresql.org/docs/current/libpq-connect.html#LIBPQ-PARAMKEYWORDS
        self.engine = create_engine(url=self.url, echo=False,
                                    connect_args=connect_args)
        
        # call to __load_spatialite() for sqlite, to load mod_spatialite via event handler listen()
        if self.driver == 'sqlite':
            log.debug('loading spatialite extension')
            listen(target=self.engine, identifier='connect', fn=self.__load_spatialite)
            # check if loading was successful
            try:
                with self.engine.begin() as conn:
                    version = conn.execute('SELECT spatialite_version();')
            except exc.OperationalError:
                raise RuntimeError('could not load spatialite extension')
        
        # if database is new, (create postgres-db and) enable spatial extension
        if not database_exists(self.engine.url):
            if self.driver == 'postgresql':
                log.debug('creating new PostgreSQL database')
                create_database(self.engine.url)
            log.debug('enabling spatial extension for new database')
            with self.engine.begin() as conn:
                if self.driver == 'sqlite':
                    conn.execute(select([func.InitSpatialMetaData(1)]))
                else:
                    conn.exec_driver_sql('CREATE EXTENSION IF NOT EXISTS postgis;')
        # create Session (ORM) and get metadata
        self.Session = sessionmaker(bind=self.engine)
        self.meta = MetaData(self.engine)
        self.custom_fields = custom_fields
        
        # load or create tables
        self.__init_data_table()
        self.__init_duplicates_table()
        
        msg = ("the 'data' table is missing {}. Please create a new database "
               "and import the old one opened in legacy mode using "
               "Archive.import_outdated.")
        pk = sql_inspect(self.data_schema).primary_key
        if 'product' not in pk.columns.keys() and not legacy:
            raise RuntimeError(msg.format("a primary key 'product'"))
        
        if 'geometry' not in self.get_colnames() and not legacy:
            raise RuntimeError(msg.format("the 'geometry' column"))
        
        self.Base = automap_base(metadata=self.meta)
        self.Base.prepare(self.engine, reflect=True)
        self.Data = self.Base.classes.data
        self.Duplicates = self.Base.classes.duplicates
        self.dbfile = dbfile
        
        if cleanup:
            log.info('checking for missing scenes')
            self.cleanup()
            sys.stdout.flush()
    
    def add_tables(
            self,
            tables: Table | list[Table],
    ) -> None:
        """
        Add tables to the database per :class:`sqlalchemy.schema.Table`
        Tables provided here will be added to the database.

        .. note::

            Columns using Geometry must have setting management=True for SQLite,
            for example: ``geometry = Column(Geometry('POLYGON', management=True, srid=4326))``

        Parameters
        ----------
        tables:
            The table(s) to be added to the database.
        """
        created = []
        if isinstance(tables, list):
            for table in tables:
                table.metadata = self.meta
                if not sql_inspect(self.engine).has_table(str(table)):
                    table.create(self.engine)
                    created.append(str(table))
        else:
            table = tables
            table.metadata = self.meta
            if not sql_inspect(self.engine).has_table(str(table)):
                table.create(self.engine)
                created.append(str(table))
        log.info('created table(s) {}.'.format(', '.join(created)))
        self.Base = automap_base(metadata=self.meta)
        self.Base.prepare(self.engine, reflect=True)
    
    def __init_data_table(self) -> None:
        if sql_inspect(self.engine).has_table('data'):
            self.data_schema = Table('data', self.meta, autoload_with=self.engine)
            return
        
        log.debug("creating DB table 'data'")
        
        self.data_schema = Table('data', self.meta,
                                 Column('sensor', String),
                                 Column('orbit', String),
                                 Column('orbitNumber_abs', Integer),
                                 Column('orbitNumber_rel', Integer),
                                 Column('cycleNumber', Integer),
                                 Column('frameNumber', Integer),
                                 Column('acquisition_mode', String),
                                 Column('start', String),
                                 Column('stop', String),
                                 Column('product', String, primary_key=True),
                                 Column('samples', Integer),
                                 Column('lines', Integer),
                                 Column('outname_base', String, primary_key=True),
                                 Column('scene', String),
                                 Column('hh', Integer),
                                 Column('vv', Integer),
                                 Column('hv', Integer),
                                 Column('vh', Integer),
                                 Column('geometry', Geometry(geometry_type='POLYGON',
                                                             management=True, srid=4326)))
        # add custom fields
        if self.custom_fields is not None:
            for key, val in self.custom_fields.items():
                if val in ['Integer', 'integer', 'int']:
                    self.data_schema.append_column(Column(key, Integer))
                elif val in ['String', 'string', 'str']:
                    self.data_schema.append_column(Column(key, String))
                else:
                    log.info('Value in dict custom_fields must be "integer" or "string"!')
        
        self.data_schema.create(self.engine)
    
    def __init_duplicates_table(self) -> None:
        # create tables if not existing
        if sql_inspect(self.engine).has_table('duplicates'):
            self.duplicates_schema = Table('duplicates', self.meta, autoload_with=self.engine)
            return
        
        log.debug("creating DB table 'duplicates'")
        
        self.duplicates_schema = Table('duplicates', self.meta,
                                       Column('outname_base', String, primary_key=True),
                                       Column('scene', String, primary_key=True))
        self.duplicates_schema.create(self.engine)
    
    @staticmethod
    def __load_spatialite(dbapi_conn: sqlite3.Connection, connection_record: Any) -> None:
        """
        loads the spatialite extension for SQLite, not to be used outside the init()

        Parameters
        ----------
        dbapi_conn:
            db engine
        connection_record:
            not sure what it does, but it is needed by :func:`sqlalchemy.event.listen`
        """
        dbapi_conn.enable_load_extension(True)
        # check which platform and use according mod_spatialite
        if platform.system() == 'Linux':
            for option in ['mod_spatialite', 'mod_spatialite.so']:
                try:
                    dbapi_conn.load_extension(option)
                except sqlite3.OperationalError:
                    continue
        elif platform.system() == 'Darwin':
            for option in ['mod_spatialite.so', 'mod_spatialite.7.dylib',
                           'mod_spatialite.dylib']:
                try:
                    dbapi_conn.load_extension(option)
                except sqlite3.OperationalError:
                    continue
        else:
            dbapi_conn.load_extension('mod_spatialite')
    
    def __prepare_insertion(self, scene: str | ID) -> Any:
        """
        read scene metadata and parse a string for inserting it into the database

        Parameters
        ----------
        scene:
            a SAR scene

        Returns
        -------
            object of class Data
        """
        id = scene if isinstance(scene, ID) else identify(scene)
        pols = [x.lower() for x in id.polarizations]
        # insertion as an object of Class Data (reflected in the init())
        insertion = self.Data()
        colnames = self.get_colnames()
        for attribute in colnames:
            if attribute == 'geometry':
                geom = id.geometry()
                geom.reproject(4326)
                geom = geom.convert2wkt(set3D=False)[0]
                geom = 'SRID=4326;' + str(geom)
                # set attributes of the Data object according to input
                setattr(insertion, 'geometry', geom)
            elif attribute in ['hh', 'vv', 'hv', 'vh']:
                setattr(insertion, attribute, int(attribute in pols))
            else:
                if hasattr(id, attribute):
                    attr = getattr(id, attribute)
                elif attribute in id.meta.keys():
                    attr = id.meta[attribute]
                else:
                    raise AttributeError('could not find attribute {}'.format(attribute))
                value = attr() if inspect.ismethod(attr) else attr
                setattr(insertion, str(attribute), value)
        
        return insertion  # return the Data object
    
    def __select_missing(self, table: str) -> list[str]:
        """

        Parameters
        ----------
        table:
            the name of the table

        Returns
        -------
            the names of all scenes, which are no longer stored in their registered location
        """
        with self.Session() as session:
            if table == 'data':
                # using ORM query to get all scenes locations
                scenes = session.query(self.Data.scene)
            elif table == 'duplicates':
                scenes = session.query(self.Duplicates.scene)
            else:
                raise ValueError("parameter 'table' must either be 'data' or 'duplicates'")
        files = [self.to_str(x[0]) for x in scenes]
        return [x for x in files if not os.path.isfile(x)]
    
    def insert(
            self,
            scene_in: str | ID | list[str | ID],
            pbar: bool = False,
            test: bool = False
    ) -> None:
        """
        Insert one or many scenes into the database

        Parameters
        ----------
        scene_in:
            a SAR scene or a list of scenes to be inserted
        pbar:
            show a progress bar?
        test:
            should the insertion only be tested or directly be committed to the database?
        """
        
        if isinstance(scene_in, (ID, str)):
            scene_in = [scene_in]
        if not isinstance(scene_in, list):
            raise RuntimeError('scene_in must either be a string pointing to a file, a pyroSAR.ID object '
                               'or a list containing several of either')
        
        log.info('filtering scenes by name')
        scenes = self.filter_scenelist(scene_in)
        if len(scenes) == 0:
            log.info('...nothing to be done')
            return
        log.info('identifying scenes and extracting metadata')
        scenes = identify_many(scenes, pbar=pbar)
        
        if len(scenes) == 0:
            log.info('all scenes are already registered')
            return
        
        counter_regulars = 0
        counter_duplicates = 0
        list_duplicates = []
        
        message = 'inserting {0} scene{1} into database'
        log.info(message.format(len(scenes), '' if len(scenes) == 1 else 's'))
        log.debug('testing changes in temporary database')
        if pbar:
            progress = pb.ProgressBar(max_value=len(scenes))
        else:
            progress = None
        insertions = []
        with self.Session() as session:
            for i, id in enumerate(scenes):
                basename = id.outname_base()
                if not self.is_registered(id):
                    insertion = self.__prepare_insertion(id)
                    insertions.append(insertion)
                    counter_regulars += 1
                    log.debug('regular:   {}'.format(id.scene))
                elif not self.__is_registered_in_duplicates(id):
                    insertion = self.Duplicates(outname_base=basename,
                                                scene=id.scene)
                    insertions.append(insertion)
                    counter_duplicates += 1
                    log.debug('duplicate: {}'.format(id.scene))
                else:
                    list_duplicates.append(id.outname_base())
                
                if progress is not None:
                    progress.update(i + 1)
            
            if progress is not None:
                progress.finish()
            
            session.add_all(insertions)
            
            if not test:
                log.debug('committing transactions to permanent database')
                # commit changes of the session
                session.commit()
            else:
                log.info('rolling back temporary database changes')
                # roll back changes of the session
                session.rollback()
        
        message = '{0} scene{1} registered regularly'
        log.info(message.format(counter_regulars, '' if counter_regulars == 1 else 's'))
        message = '{0} duplicate{1} registered'
        log.info(message.format(counter_duplicates, '' if counter_duplicates == 1 else 's'))
    
    def is_registered(self, scene: str | ID) -> bool:
        """
        Simple check if a scene is already registered in the database.

        Parameters
        ----------
        scene:
            the SAR scene

        Returns
        -------
            is the scene already registered?
        """
        id = scene if isinstance(scene, ID) else identify(scene)
        with self.Session() as session:
            # ORM query, where scene equals id.scene, return first
            exists_data = session.query(self.Data.outname_base).filter_by(
                outname_base=id.outname_base(), product=id.product).first()
            exists_duplicates = session.query(self.Duplicates.outname_base).filter(
                self.Duplicates.outname_base == id.outname_base()).first()
        in_data = False
        in_dup = False
        if exists_data:
            in_data = len(exists_data) != 0
        if exists_duplicates:
            in_dup = len(exists_duplicates) != 0
        return in_data or in_dup
    
    def __is_registered_in_duplicates(self, scene: str | ID) -> bool:
        """
        Simple check if a scene is already registered in the database.

        Parameters
        ----------
        scene:
            the SAR scene

        Returns
        -------
            is the scene already registered?
        """
        id = scene if isinstance(scene, ID) else identify(scene)
        with self.Session() as session:
            # ORM query as in is registered
            exists_duplicates = session.query(self.Duplicates.outname_base).filter(
                self.Duplicates.outname_base == id.outname_base()).first()
        in_dup = False
        if exists_duplicates:
            in_dup = len(exists_duplicates) != 0
        return in_dup
    
    def cleanup(self) -> None:
        """
        Remove all scenes from the database, which are no longer stored in their registered location
        """
        missing = self.__select_missing('data')
        for scene in missing:
            log.info('Removing missing scene from database tables: {}'.format(scene))
            self.drop_element(scene, with_duplicates=True)
    
    @staticmethod
    def to_str(string: str | bytes, encoding: str = 'utf-8') -> str:
        if isinstance(string, bytes):
            return string.decode(encoding)
        else:
            return string
    
    def export2shp(self, path: str, table: str = 'data') -> None:
        """
        export the database to a shapefile

        Parameters
        ----------
        path:
            the path of the shapefile to be written.
            This will overwrite other files with the same name.
            If a folder is given in path it is created if not existing.
            If the file extension is missing '.shp' is added.
        table:
            the table to write to the shapefile; either 'data' (default) or 'duplicates'
        """
        if table not in self.get_tablenames():
            log.warning('Only data and duplicates can be exported!')
            return
        
        # add the .shp extension if missing
        if not path.endswith('.shp'):
            path += '.shp'
        
        # creates folder if not present, adds .shp if not within the path
        dirname = os.path.dirname(path)
        os.makedirs(dirname, exist_ok=True)
        
        launder_names = {'acquisition_mode': 'acq_mode',
                         'orbitNumber_abs': 'orbit_abs',
                         'orbitNumber_rel': 'orbit_rel',
                         'cycleNumber': 'cycleNr',
                         'frameNumber': 'frameNr',
                         'outname_base': 'outname'}
        
        sel_tables = ', '.join([f'"{s}" as {launder_names[s]}' if s in launder_names else s
                                for s in self.get_colnames(table)])
        
        if self.driver == 'sqlite':
            srcDS = self.dbfile
        elif self.driver == 'postgresql':
            srcDS = """PG:host={host} port={port} user={username}
            dbname={database} password={password} active_schema=public""".format(**self.url_dict)
        else:
            raise RuntimeError('unknown archive driver')
        
        gdal.VectorTranslate(destNameOrDestDS=path, srcDS=srcDS,
                             format='ESRI Shapefile',
                             SQLStatement=f'SELECT {sel_tables} FROM {table}',
                             SQLDialect=self.driver)
    
    def filter_scenelist(self, scenelist: list[str | ID]) -> list[str | ID]:
        """
        Filter a list of scenes by file names already registered in the database.

        Parameters
        ----------
        scenelist:
            the scenes to be filtered

        Returns
        -------
            The objects of `scenelist` for all scenes whose basename
            is not yet registered in the database.

        """
        for item in scenelist:
            if not isinstance(item, (ID, str)):
                raise TypeError("items in scenelist must be of type 'str' or 'pyroSAR.ID'")
        
        with self.Session() as session:
            # ORM query, get all scenes locations
            scenes_data = session.query(self.Data.scene)
            registered = [os.path.basename(self.to_str(x[0])) for x in scenes_data]
            scenes_duplicates = session.query(self.Duplicates.scene)
        duplicates = [os.path.basename(self.to_str(x[0])) for x in scenes_duplicates]
        names = [item.scene if isinstance(item, ID) else item for item in scenelist]
        filtered = [x for x, y in zip(scenelist, names)
                    if os.path.basename(y) not in registered + duplicates]
        return filtered
    
    def get_colnames(self, table: str = 'data') -> list[str]:
        """
        Return the names of all columns of a table.

        Returns
        -------
            the column names of the chosen table
        """
        # get all columns of `table`, but shows geometry columns not correctly
        table_info = Table(table, self.meta, autoload=True, autoload_with=self.engine)
        col_names = table_info.c.keys()
        
        return sorted([self.to_str(x) for x in col_names])
    
    def get_tablenames(self, return_all: bool = False) -> list[str]:
        """
        Return the names of all tables in the database

        Parameters
        ----------
        return_all:
            only gives tables data and duplicates on default.
            Set to True to get all other tables and views created automatically.

        Returns
        -------
            the table names
        """
        #  TODO: make this dynamic
        #  the method was intended to only return user generated tables by default, as well as data and duplicates
        all_tables = ['ElementaryGeometries', 'SpatialIndex', 'geometry_columns', 'geometry_columns_auth',
                      'geometry_columns_field_infos', 'geometry_columns_statistics', 'geometry_columns_time',
                      'spatial_ref_sys', 'spatial_ref_sys_aux', 'spatialite_history', 'sql_statements_log',
                      'sqlite_sequence', 'views_geometry_columns', 'views_geometry_columns_auth',
                      'views_geometry_columns_field_infos', 'views_geometry_columns_statistics',
                      'virts_geometry_columns', 'virts_geometry_columns_auth', 'virts_geometry_columns_field_infos',
                      'virts_geometry_columns_statistics', 'data_licenses', 'KNN']
        # get tablenames from metadata
        tables = sorted([self.to_str(x) for x in self.meta.tables.keys()])
        if return_all:
            return tables
        else:
            ret = []
            for i in tables:
                if i not in all_tables and 'idx_' not in i:
                    ret.append(i)
            return ret
    
    def get_unique_directories(self) -> list[str]:
        """
        Get a list of directories containing registered scenes

        Returns
        -------
            the directory names
        """
        with self.Session() as session:
            # ORM query, get all directories
            scenes = session.query(self.Data.scene)
        registered = [os.path.dirname(self.to_str(x[0])) for x in scenes]
        return list(set(registered))
    
    def import_outdated(self, dbfile: str | Archive) -> None:
        """
        import an older database

        Parameters
        ----------
        dbfile:
            the old database. If this is a string, the name of a CSV file is expected.
        """
        if isinstance(dbfile, str) and dbfile.endswith('csv'):
            with open(dbfile) as csvfile:
                text = csvfile.read()
                csvfile.seek(0)
                dialect = csv.Sniffer().sniff(text)
                reader = csv.DictReader(csvfile, dialect=dialect)
                scenes = []
                for row in reader:
                    scenes.append(row['scene'])
                self.insert(scenes)
        elif isinstance(dbfile, Archive):
            with self.engine.begin() as conn:
                scenes = conn.exec_driver_sql('SELECT scene from data')
                scenes = [s.scene for s in scenes]
            self.insert(scenes)
            reinsert = dbfile.select_duplicates(value='scene')
            if reinsert is not None:
                self.insert(reinsert)
        else:
            raise RuntimeError("'dbfile' must either be a CSV file name or an Archive object")
    
    def move(self, scenelist: list[str], directory: str, pbar: bool = False) -> None:
        """
        Move a list of files while keeping the database entries up to date.
        If a scene is registered in the database (in either the data or duplicates table),
        the scene entry is directly changed to the new location.

        Parameters
        ----------
        scenelist:
            the file locations
        directory:
            a folder to which the files are moved
        pbar:
            show a progress bar?
        """
        if not os.path.isdir(directory):
            os.mkdir(directory)
        if not os.access(directory, os.W_OK):
            raise RuntimeError('directory cannot be written to')
        failed = []
        double = []
        if pbar:
            progress = pb.ProgressBar(max_value=len(scenelist)).start()
        else:
            progress = None
        
        for i, scene in enumerate(scenelist):
            new = os.path.join(directory, os.path.basename(scene))
            if os.path.isfile(new):
                double.append(new)
                continue
            try:
                shutil.move(scene, directory)
            except shutil.Error:
                failed.append(scene)
                continue
            finally:
                if progress is not None:
                    progress.update(i + 1)
            if self.select(scene=scene) != 0:
                table = 'data'
            else:
                # using core connection to execute SQL syntax (as was before)
                query = '''SELECT scene FROM duplicates WHERE scene='{0}' '''.format(scene)
                with self.engine.begin() as conn:
                    query_duplicates = conn.exec_driver_sql(query)
                if len(query_duplicates) != 0:
                    table = 'duplicates'
                else:
                    table = None
            if table:
                # using core connection to execute SQL syntax (as was before)
                query = '''UPDATE {0} SET scene= '{1}' WHERE scene='{2}' '''.format(table, new, scene)
                with self.engine.begin() as conn:
                    conn.exec_driver_sql(query)
        if progress is not None:
            progress.finish()
        
        if len(failed) > 0:
            log.info('The following scenes could not be moved:\n{}'.format('\n'.join(failed)))
        if len(double) > 0:
            log.info('The following scenes already exist at the target location:\n{}'.format('\n'.join(double)))
    
    def select(
            self,
            sensor: str | list[str] | None = None,
            product: str | list[str] | None = None,
            acquisition_mode: str | list[str] | None = None,
            mindate: str | datetime | None = None,
            maxdate: str | datetime | None = None,
            vectorobject: Vector | None = None,
            date_strict: bool = True,
            processdir: str | None = None,
            recursive: bool = False,
            polarizations: list[str] | None = None,
            return_value: str | list[str] = "scene",
            **kwargs: Any
    ) -> list[str | bytes] | list[tuple[str | bytes]]:
        """
        select scenes from the database

        Parameters
        ----------
        sensor:
            the satellite sensor(s)
        product:
            the product type(s)
        acquisition_mode:
            the sensor's acquisition mode(s)
        mindate:
            the minimum acquisition date; strings must be in format YYYYmmddTHHMMSS; default: None
        maxdate:
            the maximum acquisition date; strings must be in format YYYYmmddTHHMMSS; default: None
        vectorobject:
            a geometry with which the scenes need to overlap. The object may only contain one feature.
        date_strict: bool
            treat dates as strict limits or also allow flexible limits to incorporate scenes
            whose acquisition period overlaps with the defined limit?

            - strict: start >= mindate & stop <= maxdate
            - not strict: stop >= mindate & start <= maxdate
        processdir:
            A directory to be scanned for already processed scenes;
            the selected scenes will be filtered to those that have not yet been processed. Default: None
        recursive:
            (only if `processdir` is not None) should also the subdirectories of the `processdir` be scanned?
        polarizations:
            a list of polarization strings, e.g. ['HH', 'VV']
        return_value:
            the query return value(s). Options:

            - `geometry_wkb`: the scene's footprint geometry formatted as WKB
            - `geometry_wkt`: the scene's footprint geometry formatted as WKT
            - `mindate`: the acquisition start datetime in UTC formatted as YYYYmmddTHHMMSS
            - `maxdate`: the acquisition end datetime in UTC formatted as YYYYmmddTHHMMSS
            - all further database column names (see :meth:`~Archive.get_colnames()`)

        **kwargs:
            any further arguments (columns), which are registered in the database. See :meth:`~Archive.get_colnames()`

        Returns
        -------
            If a single return_value is specified: list of values for that attribute.
            If multiple return_values are specified: list of tuples containing the requested attributes.
            The return value type is bytes for `geometry_wkb` and str for all others.
        """
        # Convert return_value to list if it's a string
        if isinstance(return_value, str):
            return_values = [return_value]
        else:
            return_values = return_value
        
        return_values_sql = []
        for val in return_values:
            if val == 'mindate':
                return_values_sql.append('start')
            elif val == 'maxdate':
                return_values_sql.append('stop')
            elif val == 'geometry_wkt':
                prefix = 'ST_' if self.driver == 'postgresql' else ''
                return_values_sql.append(f'{prefix}AsText(geometry) as geometry_wkt')
            elif val == 'geometry_wkb':
                prefix = 'ST_' if self.driver == 'postgresql' else ''
                return_values_sql.append(f'{prefix}AsBinary(geometry) as geometry_wkb')
            else:
                return_values_sql.append(val)
        
        # Validate that all requested return values exist in the database
        valid_columns = self.get_colnames()
        extra = ['mindate', 'maxdate', 'geometry_wkt', 'geometry_wkb']
        normal_returns = [x for x in return_values if x not in extra]
        invalid_returns = [x for x in normal_returns if x not in valid_columns]
        if invalid_returns:
            invalid_str = ', '.join(invalid_returns)
            msg = (f"The following options are not supported as "
                   f"return values: {invalid_str}")
            raise ValueError(msg)
        
        arg_valid = [x for x in kwargs.keys() if x in self.get_colnames()]
        arg_invalid = [x for x in kwargs.keys() if x not in self.get_colnames()]
        if len(arg_invalid) > 0:
            log.info(f"the following arguments will be ignored as they are not "
                     f"registered in the data base: {', '.join(arg_invalid)}")
        
        def convert_general(k: str, v: Any) -> str:
            if isinstance(v, (float, int, str)):
                return f"""{k}='{v}'"""
            elif isinstance(v, (tuple, list)):
                v_str = "', '".join(map(str, v))
                return f"""{k} IN ('{v_str}')"""
            else:
                raise TypeError(f"unsupported type for '{k}': {type(v)}")
        
        arg_format = []
        vals = []
        for key in arg_valid:
            if key == 'scene':
                arg_format.append('''scene LIKE '%%{0}%%' '''.format(os.path.basename(kwargs[key])))
            else:
                arg_format.append(convert_general(key, kwargs[key]))
        
        if sensor:
            arg_format.append(convert_general('sensor', sensor))
        
        if product:
            arg_format.append(convert_general('product', product))
        
        if acquisition_mode:
            arg_format.append(convert_general('acquisition_mode', acquisition_mode))
        
        if mindate:
            if isinstance(mindate, datetime):
                mindate = mindate.strftime('%Y%m%dT%H%M%S')
            if re.search('[0-9]{8}T[0-9]{6}', mindate):
                if date_strict:
                    arg_format.append('start>=?')
                else:
                    arg_format.append('stop>=?')
                vals.append(mindate)
            else:
                log.info('WARNING: argument mindate is ignored, must be in format YYYYmmddTHHMMSS')
        
        if maxdate:
            if isinstance(maxdate, datetime):
                maxdate = maxdate.strftime('%Y%m%dT%H%M%S')
            if re.search('[0-9]{8}T[0-9]{6}', maxdate):
                if date_strict:
                    arg_format.append('stop<=?')
                else:
                    arg_format.append('start<=?')
                vals.append(maxdate)
            else:
                log.info('WARNING: argument maxdate is ignored, must be in format YYYYmmddTHHMMSS')
        
        if polarizations:
            for pol in polarizations:
                if pol in ['HH', 'VV', 'HV', 'VH']:
                    arg_format.append('{}=1'.format(pol.lower()))
        
        if vectorobject:
            if isinstance(vectorobject, Vector):
                if vectorobject.nfeatures > 1:
                    raise RuntimeError("'vectorobject' contains more than one feature.")
                with vectorobject.clone() as vec:
                    vec.reproject(4326)
                    site_geom = vec.convert2wkt(set3D=False)[0]
                # postgres has a different way to store geometries
                if self.driver == 'postgresql':
                    statement = f"st_intersects(geometry, 'SRID=4326; {site_geom}')"
                    arg_format.append(statement)
                else:
                    arg_format.append('st_intersects(GeomFromText(?, 4326), geometry) = 1')
                    vals.append(site_geom)
            else:
                log.info('WARNING: argument vectorobject is ignored, must be of type spatialist.vector.Vector')
        
        if len(arg_format) > 0:
            subquery = ' WHERE {}'.format(' AND '.join(arg_format))
        else:
            subquery = ''
        
        # Modify the query to select the requested return values
        query = 'SELECT {}, outname_base FROM data{}'.format(', '.join(return_values_sql), subquery)
        
        # the query gets assembled stepwise here
        for val in vals:
            query = query.replace('?', """'{0}'""", 1).format(val)
        log.debug(query)
        
        # core SQL execution
        with self.engine.begin() as conn:
            query_rs = conn.exec_driver_sql(query)
            
            if processdir and os.path.isdir(processdir):
                scenes = [x for x in query_rs
                          if len(finder(processdir, [x[-1]],
                                        regex=True, recursive=recursive)) == 0]
            else:
                scenes = query_rs
            
            ret = []
            for x in scenes:
                # If only one return value was requested, append just that value
                if len(return_values) == 1:
                    ret.append(self.to_str(x[0]))
                else:
                    # If multiple return values were requested, append a tuple of all values
                    values = []
                    for k, v in zip(return_values, x[:-1]):  # Exclude outname_base
                        if k == 'geometry_wkb':
                            values.append(v)
                        else:
                            values.append(self.to_str(v))
                    ret.append(tuple(values))
        return ret
    
    def select_duplicates(
            self,
            outname_base: str | None = None,
            scene: str | None = None,
            value: Literal["id", "scene"] = "id"
    ) -> list[str]:
        """
        Select scenes from the duplicates table. In case both `outname_base` and `scene` are set to None all scenes in
        the table are returned, otherwise only those that match the attributes `outname_base` and `scene` if they are not None.

        Parameters
        ----------
        outname_base:
            the basename of the scene
        scene:
            the scene name
        value:
            the return value; either 'id' or 'scene'

        Returns
        -------
            the selected scene(s)
        """
        if value == 'id':
            key = 0
        elif value == 'scene':
            key = 1
        else:
            raise ValueError("argument 'value' must be either 0 or 1")
        
        with self.engine.begin() as conn:
            if not outname_base and not scene:
                # core SQL execution
                scenes = conn.exec_driver_sql('SELECT * from duplicates')
            else:
                cond = []
                arg = []
                if outname_base:
                    cond.append('outname_base=?')
                    arg.append(outname_base)
                if scene:
                    cond.append('scene=?')
                    arg.append(scene)
                query = 'SELECT * from duplicates WHERE {}'.format(' AND '.join(cond))
                for a in arg:
                    query = query.replace('?', ''' '{0}' ''', 1).format(a)
                # core SQL execution
                scenes = conn.exec_driver_sql(query)
            
            ret = []
            for x in scenes:
                ret.append(self.to_str(x[key]))
        
        return ret
    
    @property
    def size(self) -> tuple[int, int]:
        """
        get the number of scenes registered in the database

        Returns
        -------
            the number of scenes in (1) the main table and (2) the duplicates table
        """
        # ORM query
        with self.Session() as session:
            r1 = session.query(self.Data.outname_base).count()
            r2 = session.query(self.Duplicates.outname_base).count()
        return r1, r2
    
    def __enter__(self) -> Archive:
        return self
    
    def close(self) -> None:
        """
        close the database connection
        """
        self.engine.dispose()
        gc.collect(generation=2)  # this was added as a fix for win PermissionError when deleting sqlite.db files.
    
    def __exit__(
            self,
            exc_type: type[BaseException] | None,
            exc_val: BaseException | None,
            exc_tb: TracebackType | None
    ) -> None:
        self.close()
    
    def drop_element(self, scene: str, with_duplicates: bool = False) -> None:
        """
        Drop a scene from the data table.
        If the duplicates table contains a matching entry, it will be moved to the data table.

        Parameters
        ----------
        scene:
            a SAR scene
        with_duplicates:
            True: delete matching entry in duplicates table
            False: move matching entry from duplicates into data table
        """
        # save outname_base from to be deleted entry
        search = self.data_schema.select().where(self.data_schema.c.scene == scene)
        entry_data_outname_base = []
        with self.engine.begin() as conn:
            for rowproxy in conn.execute(search):
                entry_data_outname_base.append((rowproxy[12]))
        # log.info(entry_data_outname_base)
        
        # delete entry in data table
        delete_statement = self.data_schema.delete().where(self.data_schema.c.scene == scene)
        with self.engine.begin() as conn:
            conn.execute(delete_statement)
        
        return_sentence = 'Entry with scene-id: \n{} \nwas dropped from data'.format(scene)
        
        # with_duplicates == True, delete entry from duplicates
        if with_duplicates:
            delete_statement_dup = self.duplicates_schema.delete().where(
                self.duplicates_schema.c.outname_base == entry_data_outname_base[0])
            with self.engine.begin() as conn:
                conn.execute(delete_statement_dup)
            
            log.info(return_sentence + ' and duplicates!'.format(scene))
            return
        
        # else select scene info matching outname_base from duplicates
        select_in_duplicates_statement = self.duplicates_schema.select().where(
            self.duplicates_schema.c.outname_base == entry_data_outname_base[0])
        entry_duplicates_scene = []
        with self.engine.begin() as conn:
            for rowproxy in conn.execute(select_in_duplicates_statement):
                entry_duplicates_scene.append((rowproxy[1]))
        
        # check if there is a duplicate
        if len(entry_duplicates_scene) == 1:
            # remove entry from duplicates
            delete_statement_dup = self.duplicates_schema.delete().where(
                self.duplicates_schema.c.outname_base == entry_data_outname_base[0])
            with self.engine.begin() as conn:
                conn.execute(delete_statement_dup)
            
            # insert scene from duplicates into data
            self.insert(entry_duplicates_scene[0])
            
            return_sentence += ' and entry with outname_base \n{} \nand scene \n{} \n' \
                               'was moved from duplicates into data table'.format(
                entry_data_outname_base[0], entry_duplicates_scene[0])
        
        log.info(return_sentence + '!')
    
    def drop_table(self, table: str) -> None:
        """
        Drop a table from the database.

        Parameters
        ----------
        table:
            the table name
        """
        if table in self.get_tablenames(return_all=True):
            # this removes the idx tables and entries in geometry_columns for sqlite databases
            if self.driver == 'sqlite':
                with self.engine.begin() as conn:
                    query = "SELECT f_table_name FROM geometry_columns"
                    tab_with_geom = [rowproxy[0] for rowproxy
                                     in conn.exec_driver_sql(query)]
                    if table in tab_with_geom:
                        conn.exec_driver_sql("SELECT DropGeoTable('" + table + "')")
            else:
                table_info = Table(table, self.meta, autoload=True, autoload_with=self.engine)
                table_info.drop(self.engine)
            log.info('table {} dropped from database.'.format(table))
        else:
            raise ValueError("table {} is not registered in the database!".format(table))
        self.Base = automap_base(metadata=self.meta)
        self.Base.prepare(self.engine, reflect=True)
    
    @staticmethod
    def __is_open(ip: str, port: str | int) -> bool:
        """
        Checks server connection, from Ben Curtis (github: Fmstrat)

        Parameters
        ----------
        ip:
            ip of the server
        port:
            port of the server

        Returns
        -------
            is the server reachable?

        """
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.settimeout(3)
        try:
            s.connect((ip, int(port)))
            s.shutdown(socket.SHUT_RDWR)
            return True
        except:
            return False
        finally:
            s.close()
    
    def __check_host(self, ip: str, port: str | int) -> bool:
        """
        Calls __is_open() on ip and port, from Ben Curtis (github: Fmstrat)

        Parameters
        ----------
        ip:
            ip of the server
        port:
            port of the server

        Returns
        -------
            is the server reachable?
        """
        ipup = False
        for i in range(2):
            if self.__is_open(ip, port):
                ipup = True
                break
            else:
                time.sleep(5)
        return ipup


def drop_archive(archive: Archive) -> None:
    """
    drop (delete) a scene database

    Parameters
    ----------
    archive:
        the database to be deleted

    See Also
    --------
    :func:`sqlalchemy_utils.functions.drop_database()`

    Examples
    --------
    >>> pguser = os.environ.get('PGUSER')
    >>> pgpassword = os.environ.get('PGPASSWORD')

    >>> db = Archive('test', postgres=True, port=5432, user=pguser, password=pgpassword)
    >>> drop_archive(db)
    """
    if archive.driver == 'postgresql':
        url = archive.url
        archive.close()
        drop_database(url)
    else:
        raise RuntimeError('this function only works for PostgreSQL databases.'
                           'For SQLite databases it is recommended to just delete the DB file.')
