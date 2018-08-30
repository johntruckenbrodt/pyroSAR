import os
import math
from spatialist import sqlite_setup
from spatialist.ancillary import HiddenPrints
from datetime import datetime, timedelta


def passdb_create(ers1passes, ers2passes, dbname):
    """
    create a sqlite database from ERS pass tables
    downloaded from http://www.deos.tudelft.nl/ers/phases/starttimes.html.
    There you can also find additional information on the file structure and background.
    The fields `phase`, `cycle`, `pass`, `starttime` and `endtime` are read from the table.
    The latter two are converted to format YYYY-MM-DD HH:MM:SS.SSS.
    The fields `cycle` and `pass` are converted to integer.
    All five fields plus the name of the sensor (`ERS1` or `ERS2`) are then stored to the database.

    Parameters
    ----------
    ers1passes: str
        the name of the ERS-1 pass table
    ers2passes: str
        the name of the ERS-2 pass table
    dbname: str
        the name of the database to write the results to

    Returns
    -------

    """
    columns = {'satellite': 'TEXT',
               'phase': 'TEXT',
               'cycleNumber': 'INTEGER',
               'passNumber': 'INTEGER',
               'starttime': 'TEXT',
               'endtime': 'TEXT'}
    
    con = sqlite_setup(driver=dbname)
    
    create_string = '''CREATE TABLE if not exists data ({})'''.format(
        ', '.join([' '.join(x) for x in columns.items()]))
    cursor = con.cursor()
    cursor.execute(create_string)
    
    def time_convert(timestring):
        dt = datetime(1985, 1, 1) + timedelta(seconds=float(timestring))
        return dt.strftime('%Y-%m-%d %H:%M:%S.%f')
    
    insert_string = '''INSERT INTO data({0}) VALUES({1})''' \
        .format(', '.join(columns.keys()),
                ', '.join(['?'] * len(columns.keys())))
    
    for satellite, filename in [('ERS1', ers1passes), ('ERS2', ers2passes)]:
        with open(filename, 'r') as table:
            for line in table:
                phase, cycle, passNumber, starttime, endtime = line.split()[0:5]
                insertion = [satellite, phase,
                             int(cycle), int(passNumber),
                             time_convert(starttime), time_convert(endtime)]
                if satellite == 'ERS1':
                    print(tuple(insertion))
                cursor.execute(insert_string, tuple(insertion))
    con.commit()
    con.close()


def passdb_query(satellite, acqtime, dbname=None):
    """
    query the orbit information for an ERS acquisition

    Parameters
    ----------
    satellite: {'ERS1', 'ERS2'}
        the name of the satellite
    acqtime: datetime.datetime
        the acquisition of the satellite image
    dbname: str, None
        the name of the database as created by :func:`passdb_create`. If None, the default database delivered with
        pyroSAR is used

    Returns
    -------

    """
    if satellite == 'ERS1':
        # the last timestamp for which specific ERS-1 orbit information is present,
        # afterwards that of ERS-2 is used
        last = datetime.strptime('1996-06-02 21:59:26.618659', '%Y-%m-%d %H:%M:%S.%f')
        sat = 'ERS2' if acqtime > last else 'ERS1'
    elif satellite == 'ERS2':
        sat = 'ERS2'
    else:
        raise ValueError("satellite must be either 'ERS1' or 'ERS2'")
    
    if dbname is None:
        dbname = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', 'erspasses.db')
    with HiddenPrints():
        con = sqlite_setup(driver=dbname)
    
    cursor = con.cursor()
    acqtime_str = acqtime.strftime('%Y-%m-%d %H:%M:%S.%f')
    query = '''SELECT * FROM data WHERE satellite = ? AND starttime <= ? AND endtime >= ?'''
    cursor.execute(query, (sat, acqtime_str, acqtime_str))
    
    fetch = cursor.fetchall()
    if len(fetch) == 0:
        cursor.execute(query, ('ERS2', acqtime_str, acqtime_str))
        fetch = cursor.fetchall()
    
    result = dict(zip(['satellite', 'phase', 'cycleNumber', 'passNumber'], fetch[0][0:4]))
    result['satellite'] = satellite
    result['orbitNumber_rel'] = int(math.ceil(result['passNumber'] / 2.))
    return result
