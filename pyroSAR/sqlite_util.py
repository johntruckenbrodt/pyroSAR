import os
import re

from ctypes.util import find_library


def check_loading():
    try:
        conn = sqlite3.connect(':memory:')
        conn.enable_load_extension(True)
    except sqlite3.OperationalError:
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
        out = {'sqlite': sqlite3.version}
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
        ext = os.path.splitext(os.path.basename(extension))[0]
        if re.search('spatialite', ext):
            select = None
            for option in ['mod_spatialite', 'mod_spatialite.so']:
                try:
                    self.conn.load_extension(option)
                    select = option
                    self.extensions.append(option)
                    print('loading extension {0} as {1}'.format(extension, option))
                    break
                except sqlite3.OperationalError:
                    continue

            if select is None:
                self.__load_regular('spatialite')

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
            self.__load_regular(ext)

    def __load_regular(self, extension):
        ext_mod = find_library(extension.replace('lib', ''))
        if ext_mod is None:
            raise RuntimeError('no library found for extension {}'.format(extension))
        print('loading extension {0} as {1}'.format(extension, ext_mod))
        try:
            self.conn.load_extension(ext_mod)
            self.extensions.append(ext_mod)
        except sqlite3.OperationalError:
            raise RuntimeError('failed to load extension {}'.format(ext_mod))
