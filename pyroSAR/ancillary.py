##############################################################
# ancillary routines for software pyroSAR
# John Truckenbrodt 2014-2018
##############################################################
"""
This script gathers central functions and classes for general pyroSAR applications.
"""
import re
from datetime import datetime
from._dev_config import product_pattern

try:
    import pathos.multiprocessing as mp
except ImportError:
    pass

from spatialist.ancillary import dissolve, finder, multicore, parse_literal, rescale, run, union, \
    urlQueryParser, which


def groupbyTime(images, function, time):
    """
    function to group images by their acquisition time difference

    Parameters
    ----------
    images: list of str
        a list of image names
    function: function
        a function to derive the time from the image names; see e.g. :func:`seconds`
    time: int or float
        a time difference in seconds by which to group the images

    Returns
    -------
    list
        a list of sub-lists containing the grouped images
    """
    # sort images by time stamp
    srcfiles = sorted(images, key=function)

    groups = [[srcfiles[0]]]
    group = groups[0]

    for i in range(1, len(srcfiles)):
        item = srcfiles[i]
        if 0 < abs(function(item) - function(group[-1])) <= time:
            group.append(item)
        else:
            groups.append([item])
            group = groups[-1]
    return [x[0] if len(x) == 1 else x for x in groups]


def seconds(filename):
    """
    function to extract time in seconds from a file name.
    the format must follow a fixed pattern: YYYYmmddTHHMMSS
    Images processed with pyroSAR functionalities via module snap or gamma will contain this information.

    Parameters
    ----------
    filename: str
        the name of a file from which to extract the time from

    Returns
    -------
    float
        the difference between the time stamp in filename and Jan 01 1900 in seconds
    """
    # return mktime(strptime(re.findall('[0-9T]{15}', filename)[0], '%Y%m%dT%H%M%S'))
    td = datetime.strptime(re.findall('[0-9T]{15}', filename)[0], '%Y%m%dT%H%M%S') - datetime(1900, 1, 1)
    return td.total_seconds()


def parse_productname(name, parse_date=False):
    """
    Parse the name of a pyroSAR processing product and extract its metadata components as dictionary
    
    Parameters
    ----------
    name: str
        the name of the file to be parsed
    parse_date: bool
        parse the start date to a :class:`~datetime.datetime` object or just return the string?

    Returns
    -------
    dict
        the metadata attributes
    
    Examples
    --------
    >>> meta = parse_productname('S1A__IW___A_20150309T173017_VV_grd_mli_geo_norm_db.tif')
    >>> print(list(meta.keys()))
    ['sensor', 'acquisition_mode', 'orbit', 'start', 'extensions', 'polarization', 'proc_steps']
    """
              
    match = re.match(re.compile(product_pattern), name)
    if not match:
        return
    out = match.groupdict()
    if out['extensions'] == '':
        out['extensions'] = None
    if len(out['proc_steps']) > 0:
        out['proc_steps'] = out['proc_steps'].split('_')
    else:
        out['proc_steps'] = None
    if parse_date:
        out['start'] = datetime.strptime(out['start'], '%Y%m%dT%H%M%S')
    return out
