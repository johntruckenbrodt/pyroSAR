##############################################################
# core routines for software pyroSAR
# John Truckenbrodt 2014-2017
##############################################################
"""
This script gathers central functions and object instances for general applications
Please refer to the descriptions of the individual functions/instances for details
"""

import re
import sys
import fnmatch
import inspect
import os
import subprocess as sp
from StringIO import StringIO
from urllib import urlencode
from urlparse import urlparse, urlunparse

try:
    import pathos.multiprocessing as mp
except ImportError:
    pass


def blockPrint():
    """
    suppress console stdout prints, i.e. redirect them to a temporary string object
    call enablePrint() to reverse the effect
    """
    if isinstance(sys.stdout, file):
        sys.stdout = StringIO()


def dictmerge(x, y):
    """
    merge two dictionaries
    """
    z = x.copy()
    z.update(y)
    return z


def dissolve(inlist):
    """
    list and tuple flattening; e.g. [[1, 2], [3, 4]] -> [1, 2, 3, 4]
    """
    out = []
    for i in inlist:
        i = list(i) if isinstance(i, tuple) else i
        out.extend(dissolve(i)) if isinstance(i, list) else out.append(i)
    return out


def enablePrint():
    """
    reverse console stdout print suppression from function blockPrint
    """
    if not isinstance(sys.stdout, file):
        sys.stdout = sys.__stdout__


def finder(folder, matchlist, foldermode=0, regex=False, recursive=True):
    """
    function for finding files in a folder and its subdirectories
    search patterns must be given in a list
    patterns can follow the unix shell standard (default) or regular expressions (via argument regex)
    the argument recursive decides whether search is done recursively in all subdirectories or in the defined directory only
    foldermodes:
    0: no folders
    1: folders included
    2: only folders
    """
    # match patterns
    pattern = r'|'.join(matchlist if regex else [fnmatch.translate(x) for x in matchlist])
    if recursive:
        out = dissolve([[os.path.join(root, x) for x in dirs+files if re.search(pattern, x)] for root, dirs, files in os.walk(folder)])
    else:
        out = [os.path.join(folder, x) for x in os.listdir(folder) if re.search(pattern, x)]
    # exclude directories
    if foldermode == 0:
        out = [x for x in out if not os.path.isdir(x)]
    if foldermode == 2:
        out = [x for x in out if os.path.isdir(x)]
    return sorted(out)


def multicore(function, cores, multiargs, **singleargs):
    """
    wrapper for multicore process execution
    function: individual function to be applied to each process item
    cores: the number of subprocesses started/CPUs used; this value is reduced in case the number of subprocesses is smaller
    multiargs: a dictionary containing subfunction argument names as keys and lists of arguments to be distributed among the processes as values
    singleargs: all remaining arguments which are invariant among the subprocesses
    important:
    -all multiargs value lists must be of same length, i.e. all argument keys must be explicitly defined for each subprocess
    -all function arguments passed via singleargs must be provided with the full argument name and its value (i.e. argname=argval); default function args are not accepted
    if the processes return anything else than None, this function will return a list of results
    if all processes return None, this function will be of type void
    example:
    def add(x, y, z):
        return x+y+z
    multicore(add, cores=2, multiargs={'x': [1, 2]}, y=5, z=9)
    -> returns [15, 16]
    multicore(add, cores=2, multiargs={'x': [1, 2], 'y': [5, 6]}, z=9)
    -> returns [15, 17]
    """

    # compare the function arguments with the multi and single arguments and raise errors if mismatches occur
    function_check = inspect.getargspec(function)
    multiargs_check = [x for x in multiargs if x not in function_check[0]]
    singleargs_check = [x for x in singleargs if x not in function_check[0]]
    if len(multiargs_check) > 0:
        raise AttributeError('incompatible multi arguments: {0}'.format(', '.join(multiargs_check)))
    if len(singleargs_check) > 0:
        raise AttributeError('incompatible single arguments: {0}'.format(', '.join(singleargs_check)))

    # compare the list lengths of the multi arguments and raise errors if they are of different length
    arglengths = list(set([len(x) for x in multiargs]))
    if len(arglengths) > 1:
        raise AttributeError('multi argument lists of different length')

    # prevent starting more threads than necessary
    cores = cores if arglengths[0] >= cores else arglengths[0]

    # create a list of dictionaries each containing the arguments for individual function calls to be passed to the multicore processes
    processlist = [dictmerge(dict([(arg, multiargs[arg][i]) for arg in multiargs]), singleargs) for i in range(len(multiargs[multiargs.keys()[0]]))]

    # block printing of the executed function
    blockPrint()

    # start pool of processes and do the work
    pool = mp.Pool(processes=cores)
    result = pool.imap_unordered(lambda x: function(**x), processlist)
    pool.close()
    pool.join()

    # unblock the printing
    enablePrint()

    # return result

    # # evaluate the return of the processing function; if any value is not None then the whole list of results is returned
    # eval = [x for x in result if x is None]
    # if len(eval) != len(result):
    #     return result


def parse_literal(x):
    """
    return the smallest possible data type
    """
    if isinstance(x, list):
        return [parse_literal(y) for y in x]
    try:
        return int(x)
    except ValueError:
        try:
            return float(x)
        except ValueError:
            return str(x)


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


class ReadPar(object):
    """
    read processing parameter text files
    """
    def __init__(self, filename, splits='[:|\t|=\s]', type=''):
        if type == 'exe':
            splits = '[\t\n]'
        with open(filename, 'r') as infile:
            self.index = []
            for line in infile:
                if not line.startswith('#'):
                    items = filter(None, re.split(splits, line))
                    if len(items) > 1:
                        if len(items) > 2:
                            entry = items[1] if items[2] in ['m', 'decimal', 'arc-sec', 'degrees'] else items[1:]
                            entry = [x for x in items if x not in ['m', 'decimal', 'arc-sec', 'degrees']]
                            if ''.join(entry) == 'WGS1984':
                                entry = 'WGS84'
                        else:
                            entry = items[1]
                        if type == 'exe':
                            items[0] = items[0].replace(' ', '_')
                        setattr(self, items[0], entry)
                        self.index.append(items[0])


def rescale(inlist, newrange=(0, 1)):
    """
    rescale the values in a list between the values in newrange (a tuple with the new minimum and maximum)
    """
    OldMax = max(inlist)
    OldMin = min(inlist)
    OldRange = OldMax - OldMin
    NewRange = newrange[1] - newrange[0]
    return [(((float(x) - OldMin) * NewRange) / OldRange) + newrange[0] for x in inlist]


def run(cmd, outdir=None, logfile=None, inlist=None, void=True, errorpass=False):
    """
    wrapper for subprocess execution including logfile writing and command prompt piping
    """
    cmd = [str(x) for x in dissolve(cmd)]
    if outdir is None:
        outdir = os.getcwd()
    log = sp.PIPE if logfile is None else open(logfile, 'a')
    proc = sp.Popen(cmd, stdin=sp.PIPE, stdout=log, stderr=sp.PIPE, cwd=outdir, universal_newlines=True)
    instream = None if inlist is None else ''.join([str(x)+'\n' for x in inlist])
    out, err = proc.communicate(instream)
    if not errorpass and proc.returncode != 0:
        raise sp.CalledProcessError(proc.returncode, cmd, err)
    # add line for separating log entries of repeated function calls
    if logfile:
        log.write('#####################################################################\n')
        log.close()
    if not void:
        return out, err


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


def writer(filename, arguments, strfill=True):
    """
    write parameter textfile
    """
    argstring = '\n'.join(['\t'.join(x) for x in arguments])+'\n'
    if strfill:
        argstring = argstring.replace(' ', '_')
    with open(filename, 'w') as out:
        out.write(argstring)


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
