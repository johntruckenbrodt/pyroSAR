import os
import re

from ctypes.util import find_library


def check_loading():
    try:
        conn = sqlite3.connect(':memory:')
        conn.enable_load_extension(True)
    except (sqlite3.OperationalError, AttributeError):
        raise RuntimeError


errormessage = 'sqlite3 does not support loading extensions and {}; ' \
               'please refer to the pyroSAR installation instructions'
try:
    import sqlite3

    check_loading()
except RuntimeError:
    try:
        from pysqlite2 import dbapi2 as sqlite3

        check_loading()
    except ImportError:
        raise RuntimeError(errormessage.format('pysqlite2 does not exist as alternative'))
    except RuntimeError:
        raise RuntimeError(errormessage.format('neither does pysqlite2'))


def sqlite_setup(driver=':memory:', extensions=None):
    """
    setup a sqlite3 connection and load extensions to it

    Parameters
    ----------
    driver: str
        the database file or (by default) an in-memory database
    extensions: list
        a list of extensions to load

    Returns
    -------
    sqlite3.Connection
        the database connection
    """
    conn = __Handler(driver, extensions)
    return conn.conn


class __Handler(object):
    def __init__(self, driver=':memory:', extensions=None):
        self.conn = sqlite3.connect(driver)
        self.conn.enable_load_extension(True)
        self.extensions = []
        if isinstance(extensions, list):
            for ext in extensions:
                self.load_extension(ext)
        elif extensions is not None:
            raise RuntimeError('extensions must either be a list or None')
        print('using sqlite version {}'.format(self.version['sqlite']))
        if 'spatialite' in self.version.keys():
            print('using spatialite version {}'.format(self.version['spatialite']))

    @property
    def version(self):
        out = {'sqlite': sqlite3.sqlite_version}
        cursor = self.conn.cursor()
        try:
            cursor.execute('SELECT spatialite_version()')
            spatialite_version = cursor.fetchall()[0][0].encode('ascii')
            out['spatialite'] = spatialite_version
        except sqlite3.OperationalError:
            pass
        return out

    def get_tablenames(self):
        cursor = self.conn.cursor()
        cursor.execute('SELECT * FROM sqlite_master WHERE type="table"')
        return [x[1].encode('ascii') for x in cursor.fetchall()]

    def load_extension(self, extension):
        if re.search('spatialite', extension):
            select = None
            # first try to load the dedicated mod_spatialite adapter
            for option in ['mod_spatialite', 'mod_spatialite.so']:
                try:
                    self.conn.load_extension(option)
                    select = option
                    self.extensions.append(option)
                    print('loading extension {0} as {1}'.format(extension, option))
                    break
                except sqlite3.OperationalError:
                    continue

            # if loading mod_spatialite fails try to load libspatialite directly
            if select is None:
                self.__load_regular('spatialite')

            # initialize spatial support
            if 'spatial_ref_sys' not in self.get_tablenames():
                cursor = self.conn.cursor()
                if select is None:
                    # libspatialite extension
                    cursor.execute('SELECT InitSpatialMetaData();')
                else:
                    # mod_spatialite extension
                    cursor.execute('SELECT InitSpatialMetaData(1);')
                self.conn.commit()

        else:
            self.__load_regular(extension)

    def __load_regular(self, extension):
        options = []

        # create an extension library option starting with 'lib' without extension suffices;
        # e.g. 'libgdal' but not 'gdal.so'
        ext_base = self.__split_ext(extension)
        if not ext_base.startswith('lib'):
            ext_base = 'lib' + ext_base
        options.append(ext_base)

        # get the full extension library name; e.g. 'libgdal.so.20'
        ext_mod = find_library(extension.replace('lib', ''))
        if ext_mod is None:
            raise RuntimeError('no library found for extension {}'.format(extension))
        options.append(ext_mod)

        # loop through extension library name options and try to load them
        success = False
        for option in options:
            try:
                self.conn.load_extension(option)
                self.extensions.append(option)
                print('loading extension {0} as {1}'.format(extension, option))
                success = True
                break
            except sqlite3.OperationalError:
                continue

        if not success:
            raise RuntimeError('failed to load extension {}'.format(extension))

    def __split_ext(self, extension):
        base = extension
        while re.search('\.', base):
            base = os.path.splitext(base)[0]
        return base
