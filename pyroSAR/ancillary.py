###############################################################################
# ancillary routines for software pyroSAR

# Copyright (c) 2014-2024, the pyroSAR Developers.

# This file is part of the pyroSAR Project. It is subject to the
# license terms in the LICENSE.txt file found in the top-level
# directory of this distribution and at
# https://github.com/johntruckenbrodt/pyroSAR/blob/master/LICENSE.txt.
# No part of the pyroSAR project, including this file, may be
# copied, modified, propagated, or distributed except according
# to the terms contained in the LICENSE.txt file.
###############################################################################
"""
This module gathers central functions and classes for general pyroSAR applications.
"""
import os
import re
import time
import uuid
from pathlib import Path
from math import sin, radians
import inspect
from datetime import datetime
from . import patterns

from spatialist.ancillary import finder


def groupby(images, attribute):
    """
    group a list of images by a metadata attribute
    
    Parameters
    ----------
    images: list[str]
        the names of the images to be sorted
    attribute: str
        the name of the attribute used for sorting;
        see :func:`parse_datasetname` for options
    
    Returns
    -------
    list[list[str]]
        a list of sub-lists containing the grouped images
    """
    images_sort = sorted(images, key=lambda x: re.search(patterns.pyrosar, x).group(attribute))
    out_meta = [[parse_datasetname(images_sort.pop(0))]]
    while len(images_sort) > 0:
        filename = images_sort.pop(0)
        meta = parse_datasetname(filename)
        
        if out_meta[-1][0][attribute] == meta[attribute]:
            out_meta[-1].append(meta)
        else:
            out_meta.append([meta])
    out = [[x['filename'] for x in y] for y in out_meta]
    return out


def groupbyTime(images, function, time):
    """
    function to group images by their acquisition time difference

    Parameters
    ----------
    images: list[str]
        a list of image names
    function: function
        a function to derive the time from the image names; see e.g. :func:`seconds`
    time: int or float
        a time difference in seconds by which to group the images

    Returns
    -------
    list[list[str]]
        a list of sub-lists containing the grouped images
    """
    # sort images by time stamp
    srcfiles = sorted(images, key=function)
    
    groups = [[srcfiles[0]]]
    group = groups[0]
    
    for i in range(1, len(srcfiles)):
        item = srcfiles[i]
        timediff = abs(function(item) - function(group[-1]))
        if timediff <= time:
            group.append(item)
        else:
            groups.append([item])
            group = groups[-1]
    return [x[0] if len(x) == 1 else x for x in groups]


def multilook_factors(source_rg, source_az, target, geometry, incidence):
    """
    compute multi-looking factors to approximate a square pixel with defined target ground range pixel spacing.
    
    Parameters
    ----------
    source_rg: int or float
        the range pixel spacing
    source_az: int or float
        the azimuth pixel spacing
    target: int or float
        the target pixel spacing of an approximately square pixel
    geometry: str
        the imaging geometry; either 'SLANT_RANGE' or 'GROUND_RANGE'
    incidence: int or float
        the angle of incidence

    Returns
    -------
    tuple[int]
        the multi-looking factors as (range looks, azimuth looks)
    
    Examples
    --------
    >>> from pyroSAR.ancillary import multilook_factors
    >>> rlks, azlks = multilook_factors(source_rg=2, source_az=13, target=10,
    >>>                                 geometry='SLANT_RANGE', incidence=39)
    >>> print(rlks, azlks)
    4 1
    """
    azlks = int(round(float(target) / source_az))
    azlks = azlks if azlks > 0 else 1
    if geometry == 'SLANT_RANGE':
        rlks = float(azlks) * source_az * sin(radians(incidence)) / source_rg
    elif geometry == 'GROUND_RANGE':
        rlks = float(azlks) * source_az / source_rg
    else:
        raise ValueError("parameter 'geometry' must be either 'SLANT_RANGE' or 'GROUND_RANGE'")
    
    rlks = int(round(rlks))
    return rlks, azlks


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


def parse_datasetname(name, parse_date=False):
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
    >>> meta = parse_datasetname('S1A__IW___A_20150309T173017_VV_grd_mli_geo_norm_db.tif')
    >>> print(sorted(meta.keys()))
    ['acquisition_mode', 'extensions', 'filename', 'orbit',
    'outname_base', 'polarization', 'proc_steps', 'sensor', 'start']
    """
    
    filename = os.path.abspath(name) if os.path.isfile(name) else name
    
    match = re.match(re.compile(patterns.pyrosar), filename)
    if not match:
        return
    out = match.groupdict()
    if out['extensions'] == '':
        out['extensions'] = None
    if out['proc_steps'] is not None:
        out['proc_steps'] = out['proc_steps'].split('_')
    if parse_date:
        out['start'] = datetime.strptime(out['start'], '%Y%m%dT%H%M%S')
    out['filename'] = filename
    out['outname_base'] = out['outname_base'].strip('_')
    return out


def find_datasets(directory, recursive=False, **kwargs):
    """
    find pyroSAR datasets in a directory based on their metadata
    
    Parameters
    ----------
    directory: str
        the name of the directory to be searched
    recursive: bool
        search the directory recursively into subdirectories?
    kwargs:
        Metadata attributes for filtering the scene list supplied as `key=value`. e.g. `sensor='S1A'`.
        Multiple allowed options can be provided in tuples, e.g. `sensor=('S1A', 'S1B')`.
        Any types other than tuples require an exact match, e.g. `proc_steps=['grd', 'mli', 'geo', 'norm', 'db']`
        will be matched only if these processing steps are contained in the product name in this exact order.
        The special attributes `start` and `stop` can be used for time filtering where `start<=value<=stop`.
        See function :func:`parse_datasetname` for further options.
    
    Returns
    -------
    list of str
        the file names found in the directory and filtered by metadata attributes
    
    Examples
    --------
    >>> selection = find_datasets('path/to/files', sensor=('S1A', 'S1B'), polarization='VV')
    """
    files = finder(directory, [patterns.pyrosar], regex=True, recursive=recursive)
    selection = []
    for file in files:
        meta = parse_datasetname(file)
        matches = []
        for key, val in kwargs.items():
            if key == 'start':
                match = val <= meta['start']
            elif key == 'stop':
                match = val >= meta['start']  # only the start time stamp is contained in the filename
            elif isinstance(val, tuple):
                match = meta[key] in val
            else:
                match = meta[key] == val
            matches.append(match)
        if all(matches):
            selection.append(file)
    return selection


def getargs(func):
    """
    get the arguments of a function
    
    Parameters
    ----------
    func: function
        the function to be checked

    Returns
    -------
    list or str
        the argument names
    """
    return sorted(inspect.getfullargspec(func).args)


def hasarg(func, arg):
    """
    simple check whether a function takes a parameter as input
    
    Parameters
    ----------
    func: function
        the function to be checked
    arg: str
        the argument name to be found

    Returns
    -------
    bool
        does the function take this as argument?
    """
    return arg in getargs(func)


def windows_fileprefix(func, path, exc_info):
    """
    Helper function for :func:`shutil.rmtree` to exceed Windows' file name length limit of 256 characters.
    See `here <https://stackoverflow.com/questions/36219317/pathname-too-long-to-open>`_ for details.

    Parameters
    ----------
    func: function
        the function to be executed, i.e. :func:`shutil.rmtree`
    path: str
        the path to be deleted
    exc_info: tuple
        execution info as returned by :func:`sys.exc_info`

    Returns
    -------

    Examples
    --------
    >>> import shutil
    >>> from pyroSAR.ancillary import windows_fileprefix
    >>> shutil.rmtree('/path', onerror=windows_fileprefix)
    """
    func(u'\\\\?\\' + path)


class Lock(object):
    """
    File and folder locking mechanism.
    This mechanism creates lock files indicating whether a file/folder
    
     1. is being modified (`target`.lock),
     2. is being used/read (`target`.used_<uuid.uuid4>) or
     3. was damaged during modification (`target`.error).
    
    Although these files will not prevent locking by other mechanisms (UNIX
    locks are generally only advisory), this mechanism is respected across
    any running instances. I.e., if such a lock file exists, no process
    trying to acquire a lock using this class will succeed if a lock file
    intending to prevent it exists. This was implemented because other existing
    solutions like `filelock <https://github.com/tox-dev/filelock>`_ or
    `fcntl <https://docs.python.org/3/library/fcntl.html>`_ do not implement
    effective solutions for parallel jobs in HPC systems.
    
    Hard locks prevent any usage of the data. Damage/error locks work like hard
    locks except that `timeout` is ignored and a `RuntimeError` is raised immediately.
    Error locks are created if an error occurs whilst a hard lock is acquired and
    `target` exists (by renaming the hard lock file).
    Infinite usage locks may exist, each with a different random UUID. No hard
    lock may be acquired whilst usage locks exist. On error usage locks are simply
    deleted.
    
    It may happen that lock files remain when a process is killed by HPC schedulers
    like Slurm because in this case the process is not ended by Python. Optimally,
    hard locks should be renamed to error lock files and usage lock files should be
    deleted. This has to be done separately.
    
    Examples
    --------
    >>> from pyroSAR.ancillary import Lock
    >>> target = 'test.txt'
    >>> with Lock(target=target):
    >>>     with open(target, 'w') as f:
    >>>         f.write('Hello World!')

    Parameters
    ----------
    target: str
        the file/folder to lock
    soft: bool
        lock the file/folder only for reading (and not for modification)?
    timeout: int
        the time in seconds to retry acquiring a lock
    """
    def __init__(self, target, soft=False, timeout=7200):
        if os.path.isdir(target) and not os.path.exists(target):
            raise OSError('target does not exist: {}'.format(target))
        self.target = target
        used_id = str(uuid.uuid4())
        self.lock = self.target + '.lock'
        self.error = self.target + '.error'
        self.used = self.target + f'.used_{used_id}'
        self.soft = soft
        if os.path.isfile(self.error):
            msg = 'cannot acquire lock on damaged target: {}'
            raise RuntimeError(msg.format(self.target))
        end = time.time() + timeout
        while True:
            if time.time() > end:
                msg = 'could not acquire lock due to timeout: {}'
                raise RuntimeError(msg.format(self.target))
            try:
                if self.soft and not os.path.isfile(self.lock):
                    Path(self.used).touch(exist_ok=False)
                    break
                if not self.soft and not self.is_used():
                    Path(self.lock).touch(exist_ok=False)
                    break
            except FileExistsError:
                pass
            time.sleep(1)
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        if not self.soft and exc_type is not None:
            if os.path.exists(self.target):
                os.rename(self.lock, self.error)
        else:
            self.remove()
    
    def is_used(self):
        """
        Does any usage lock exist?
        
        Returns
        -------
        bool
        """
        base = os.path.basename(self.target)
        folder = os.path.dirname(self.target)
        files = list(Path(folder).glob(base + '.used*'))
        return len(files) > 0
    
    def remove(self):
        """
        Remove the acquired soft/hard lock
        
        Returns
        -------

        """
        if self.soft:
            os.remove(self.used)
        else:
            os.remove(self.lock)


class LockCollection(object):
    """
    Like :class:`Lock` but for multiple files/folders.

    Parameters
    ----------
    targets: list[str]
        the files/folders to lock
    soft: bool
        lock the files/folders only for reading (and not for modification)?
    timeout: int
        the time in seconds to retry acquiring a lock
    """
    def __init__(self, targets, soft=False, timeout=7200):
        self.locks = [Lock(x, soft=soft, timeout=timeout) for x in targets]
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        for lock in self.locks:
            lock.remove()
