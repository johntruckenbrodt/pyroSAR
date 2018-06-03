import os
try:
    from pysqlite2 import dbapi2 as sqlite3
except ImportError:
    import sqlite3


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
    conn = sqlite3.connect(driver)
    conn.enable_load_extension(True)
    for item in extensions:
        ext = os.path.splitext(os.path.basename(item))[0]
        # see https://www.gaia-gis.it/fossil/libspatialite/wiki?name=mod_spatialite
        if sqlite3.version < '3.7.17':
            ext += '.so'
        conn.load_extension(ext)
    return conn
