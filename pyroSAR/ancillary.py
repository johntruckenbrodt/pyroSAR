##############################################################
# core routines for software pyroSAR
# John Truckenbrodt 2014-2018
##############################################################
"""
This script gathers central functions and object instances for general applications
Please refer to the descriptions of the individual functions/instances for details
"""
import re
from datetime import datetime

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
        a function to derive the time from the image names
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
