##############################################################
# core routines for software pyroSAR
# John Truckenbrodt 2014-2018
##############################################################
"""
This script gathers central functions and object instances for general applications
Please refer to the descriptions of the individual functions/instances for details
"""
import sys

if sys.version_info >= (3, 0):
    from io import StringIO
    from urllib.parse import urlparse, urlunparse, urlencode
    from builtins import str
else:
    from urllib import urlencode
    from StringIO import StringIO
    from urlparse import urlparse, urlunparse

import re
import sys
import fnmatch
import inspect
import itertools
import os
import subprocess as sp
from datetime import datetime

try:
    import pathos.multiprocessing as mp
except ImportError:
    pass


class HiddenPrints:
    """
    suppress console stdout prints, i.e. redirect them to a temporary string object
    adapted from https://stackoverflow.com/questions/8391411/suppress-calls-to-print-python

    Examples
    --------
    >>> with HiddenPrints():
    >>>     print('foobar')
    >>> print('foobar')
    """

    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = StringIO()

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout = self._original_stdout


def dictmerge(x, y):
    """
    merge two dictionaries
    """
    z = x.copy()
    z.update(y)
    return z


# todo consider using itertools.chain like in function finder
def dissolve(inlist):
    """
    list and tuple flattening; e.g. [[1, 2], [3, 4]] -> [1, 2, 3, 4]
    """
    out = []
    for i in inlist:
        i = list(i) if isinstance(i, tuple) else i
        out.extend(dissolve(i)) if isinstance(i, list) else out.append(i)
    return out


def finder(folder, matchlist, foldermode=0, regex=False, recursive=True):
    """
    function for finding files/folders in folders and their subdirectories

    Parameters
    ----------
    folder: str or list of str
        the directory or al ist of directories to be searched
    matchlist: list
        a list of search patterns
    foldermode: int
        0: only files;
        1: files and folders;
        2: only folders
    regex: bool
        are the search patterns in matchlist regular expressions or unix shell standard (default)?
    recursive: bool
        search folder recursively into all subdirectories or only in the top level?

    Returns
    -------
    list of str
        the absolute names of files matching the patterns
    """
    # match patterns
    if isinstance(folder, str):
        pattern = r'|'.join(matchlist if regex else [fnmatch.translate(x) for x in matchlist])
        if recursive:
            out = dissolve([[os.path.join(root, x) for x in dirs + files if re.search(pattern, x)]
                            for root, dirs, files in os.walk(folder)])
        else:
            out = [os.path.join(folder, x) for x in os.listdir(folder) if re.search(pattern, x)]
        # exclude directories
        if foldermode == 0:
            out = [x for x in out if not os.path.isdir(x)]
        if foldermode == 2:
            out = [x for x in out if os.path.isdir(x)]
        return sorted(out)
    elif isinstance(folder, list):
        groups = [finder(x, matchlist, foldermode, regex, recursive) for x in folder]
        return list(itertools.chain(*groups))
    else:
        raise TypeError("parameter 'folder' must be of type str or list")


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


def multicore(function, cores, multiargs, **singleargs):
    """
    wrapper for multicore process execution

    Parameters
    ----------
    function
        individual function to be applied to each process item
    cores: int
        the number of subprocesses started/CPUs used;
        this value is reduced in case the number of subprocesses is smaller
    multiargs: dict
        a dictionary containing sub-function argument names as keys and lists of arguments to be
        distributed among the processes as values
    singleargs
        all remaining arguments which are invariant among the subprocesses

    Returns
    -------
    None or list
        the return of or function for all the subprocesses

    Notes
    -----
    - all multiargs value lists must be of same length, i.e. all argument keys must be explicitly defined for each
    subprocess
    - all function arguments passed via singleargs must be provided with the full argument name and its value
    (i.e. argname=argval); default function args are not accepted
    - if the processes return anything else than None, this function will return a list of results
    - if all processes return None, this function will be of type void

    Examples
    --------
    >>> def add(x, y, z):
    >>>     return x+y+z
    >>> multicore(add, cores=2, multiargs={'x': [1, 2]}, y=5, z=9)
    [15, 16]
    >>> multicore(add, cores=2, multiargs={'x': [1, 2], 'y': [5, 6]}, z=9)
    [15, 17]
    """

    # compare the function arguments with the multi and single arguments and raise errors if mismatches occur
    if sys.version_info >= (3, 0):
        check = inspect.getfullargspec(function)
        varkw = check.varkw
    else:
        check = inspect.getargspec(function)
        varkw = check.keywords

    if not check.varargs and not varkw:
        multiargs_check = [x for x in multiargs if x not in check.args]
        singleargs_check = [x for x in singleargs if x not in check.args]
        if len(multiargs_check) > 0:
            raise AttributeError('incompatible multi arguments: {0}'.format(', '.join(multiargs_check)))
        if len(singleargs_check) > 0:
            raise AttributeError('incompatible single arguments: {0}'.format(', '.join(singleargs_check)))

    # compare the list lengths of the multi arguments and raise errors if they are of different length
    arglengths = list(set([len(multiargs[x]) for x in multiargs]))
    if len(arglengths) > 1:
        raise AttributeError('multi argument lists of different length')

    # prevent starting more threads than necessary
    cores = cores if arglengths[0] >= cores else arglengths[0]

    # create a list of dictionaries each containing the arguments for individual function calls to be passed to the multicore processes
    processlist = [dictmerge(dict([(arg, multiargs[arg][i]) for arg in multiargs]), singleargs)
                   for i in range(len(multiargs[list(multiargs.keys())[0]]))]

    # block printing of the executed function
    result = None
    with HiddenPrints():
        # start pool of processes and do the work
        try:
            pool = mp.Pool(processes=cores)
        except NameError:
            raise ImportError("package 'pathos' could not be imported")
        result = pool.imap(lambda x: function(**x), processlist)
        pool.close()
        pool.join()

    # evaluate the return of the processing function;
    # if any value is not None then the whole list of results is returned
    result = list(result)
    eval = [x for x in result if x is not None]
    if len(eval) == 0:
        return None
    else:
        return result


def parse_literal(x):
    """
    return the smallest possible data type for a string

    :param x: a string to be parsed
    :return a value of type int, float or str
    """
    if isinstance(x, list):
        return [parse_literal(y) for y in x]
    elif isinstance(x, (bytes, str)):
        try:
            return int(x)
        except ValueError:
            try:
                return float(x)
            except ValueError:
                return x
    else:
        raise IOError('input must be a string or a list of strings')


class Queue(object):
    """
    classical queue implementation
    """

    def __init__(self, inlist=None):
        self.stack = [] if inlist is None else inlist

    def empty(self):
        return len(self.stack) == 0

    def length(self):
        return len(self.stack)

    def push(self, x):
        self.stack.append(x)

    def pop(self):
        if not self.empty():
            val = self.stack[0]
            del self.stack[0]
            return val


def rescale(inlist, newrange=(0, 1)):
    """
    rescale the values in a list between the values in newrange (a tuple with the new minimum and maximum)
    """
    OldMax = max(inlist)
    OldMin = min(inlist)

    if OldMin == OldMax:
        raise RuntimeError('list contains of only one unique value')

    OldRange = OldMax - OldMin
    NewRange = newrange[1] - newrange[0]
    result = [(((float(x) - OldMin) * NewRange) / OldRange) + newrange[0] for x in inlist]
    return result


def run(cmd, outdir=None, logfile=None, inlist=None, void=True, errorpass=False):
    """
    wrapper for subprocess execution including logfile writing and command prompt piping
    """
    cmd = [str(x) for x in dissolve(cmd)]
    if outdir is None:
        outdir = os.getcwd()
    log = sp.PIPE if logfile is None else open(logfile, 'a')
    proc = sp.Popen(cmd, stdin=sp.PIPE, stdout=log, stderr=sp.PIPE, cwd=outdir, universal_newlines=True)
    instream = None if inlist is None else ''.join([str(x) + '\n' for x in inlist])
    out, err = proc.communicate(instream)
    if not errorpass and proc.returncode != 0:
        raise sp.CalledProcessError(proc.returncode, cmd, err)
    # add line for separating log entries of repeated function calls
    if logfile:
        log.write('#####################################################################\n')
        log.close()
    if not void:
        return out, err


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


class Stack(object):
    """
    classical stack implementation
    input can be a list, a single value or None (i.e. Stack())
    """

    def __init__(self, inlist=None):
        if isinstance(inlist, list):
            self.stack = inlist
        elif inlist is None:
            self.stack = []
        else:
            self.stack = [inlist]

    # check whether stack is empty
    def empty(self):
        return len(self.stack) == 0

    # empty the stack
    def flush(self):
        self.stack = []

    # print the length of the stack
    def length(self):
        return len(self.stack)

    # append items to the stack; input can be a single value or a list
    def push(self, x):
        if isinstance(x, list):
            for item in x:
                self.stack.append(item)
        else:
            self.stack.append(x)

    # return the last stack element and delete it fro mthe list
    def pop(self):
        if not self.empty():
            val = self.stack[-1]
            del self.stack[-1]
            return val


def union(a, b):
    """
    union of two lists
    """
    return list(set(a) & set(b))


def urlQueryParser(url, querydict):
    """
    parse a url query
    """
    address_parse = urlparse(url)
    return urlunparse(address_parse._replace(query=urlencode(querydict)))


def which(program):
    """
    mimics UNIX's which
    taken from this post: http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    can be replaced by shutil.which() in Python 3.3
    """

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ['PATH'].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None
