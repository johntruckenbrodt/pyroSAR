import os
import re

try:
    from pysqlite2 import dbapi2 as sqlite3
except ImportError:
    import sqlite3

from ctypes.util import find_library


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
    def __init__(self, driver, extensions):
        self.conn = sqlite3.connect(driver)
        self.conn.enable_load_extension(True)
        for ext in extensions:
            self.load_extension(ext)

    def get_tablenames(self):
        cursor = self.conn.execute('''SELECT * FROM sqlite_master WHERE type="table"''')
        return [x[1].encode('ascii') for x in cursor.fetchall()]

    def load_extension(self, extension):
        ext = os.path.splitext(os.path.basename(extension))[0]
        if re.search('spatialite', ext):
            select = None
            for option in ['mod_spatialite', 'mod_spatialite.so', 'libspatialite', 'libspatialite.so']:
                try:
                    self.conn.load_extension(option)
                    select = option
                except sqlite3.OperationalError:
                    continue
            if select is None:
                raise RuntimeError('failed to load extension {}'.format(ext))
            else:
                if 'spatial_ref_sys' not in self.get_tablenames():
                    param = 1 if re.search('spatialite', select) else ''
                    self.conn.execute('SELECT InitSpatialMetaData({});'.format(param))
        else:
            ext_mod = find_library(ext.replace('lib', ''))
            if ext_mod is None:
                raise RuntimeError('no library found for extension {}'.format(extension))
            print('loading extension {0} as {1}'.format(extension, ext_mod))
            try:
                self.conn.load_extension(ext_mod)
            except sqlite3.OperationalError as e:
                print(e)
                raise RuntimeError('failed to load extension {}'.format(ext))
