##############################################################
# core routines for software pythonland
# John Truckenbrodt 2014-2015
# last update 2015-12-08
##############################################################
"""
This script gathers central functions and object instances for processing SAR images using the software GAMMA within the GUI
Please refer to the descriptions of the individual functions/instances for details
"""

import inspect
import math
import os
import subprocess as sp
from StringIO import StringIO
from Tkinter import *
from glob import glob
from gamma.error import gammaErrorHandler

from osgeo import osr

try:
    import pathos.multiprocessing as mp
except ImportError:
    pass

osr.UseExceptions()


# suppress console stdout prints, i.e. redirect them to a temporary string object
# call enablePrint() to reverse the effect
def blockPrint():
    if isinstance(sys.stdout, file):
        sys.stdout = StringIO()


# convert between epsg, wkt and proj4 spatial references
# crsText must be a osr.SpatialReference object or a string of type WKT, PROJ4 or EPSG
# crsOut must be either "wkt", "proj4", "epsg" or "osr"
# if type "osr" is selected the function will return a spatial reference object of type osr.SpatialReference()
def crsConvert(crsIn, crsOut):
    if isinstance(crsIn, osr.SpatialReference):
        srs = crsIn.Clone()
    else:
        srs = osr.SpatialReference()
        try:
            srs.ImportFromEPSG(crsIn)
        except (TypeError, RuntimeError):
            try:
                srs.ImportFromWkt(crsIn)
            except (TypeError, RuntimeError):
                try:
                    srs.ImportFromProj4(crsIn)
                except (TypeError, RuntimeError):
                    raise TypeError("crsText not recognized; must be of type WKT, PROJ4 or EPSG")
    if crsOut == "wkt":
        return srs.ExportToWkt()
    elif crsOut == "proj4":
        return srs.ExportToProj4()
    elif crsOut == "epsg":
        srs.AutoIdentifyEPSG()
        return int(srs.GetAuthorityCode(None))
    elif crsOut == "osr":
        return srs
    else:
        raise ValueError("crsOut not recognized; must be either wkt, proj4 or epsg")


# merge two dictionaries
def dictmerge(x, y):
    z = x.copy()
    z.update(y)
    return z


# list and tuple flattening
def dissolve(inlist):
    out = []
    for i in inlist:
        i = list(i) if type(i) is tuple else i
        out.extend(dissolve(i)) if type(i) is list else out.append(i)
    return out


# reverse console stdout print suppression from function blockPrint
def enablePrint():
    if not isinstance(sys.stdout, file):
        sys.stdout = sys.__stdout__


# function for finding files in a folder and its subdirectories
# search patterns must be given in a list
# patterns can follow the unix shell standard (default) or regular expressions (via argument regex)
# the argument recursive decides whether search is done recursively in all subdirectories or in the defined directory only
# foldermodes:
# 0: no folders
# 1: folders included
# 2: only folders
def finder(folder, matchlist, foldermode=0, regex=False, recursive=True):

    # match patterns
    if regex:
        if recursive:
            out = dissolve([[os.path.join(group[0], x) for x in dissolve(group) if re.search(pattern, x)] for group in os.walk(folder) for pattern in matchlist])
        else:
            out = dissolve([[os.path.join(folder, x) for x in os.listdir(folder) if re.search(pattern, x)] for pattern in matchlist])
    else:
        if recursive:
            out = list(set([f for files in [glob(os.path.join(item[0], pattern)) for item in os.walk(folder) for pattern in matchlist] for f in files]))
        else:
            out = dissolve([glob(os.path.join(folder, pattern)) for pattern in matchlist])
    # exclude directories
    if foldermode == 0:
        out = [x for x in out if not os.path.isdir(x)]
    if foldermode == 2:
        out = [x for x in out if os.path.isdir(x)]
    return sorted(out)


# compute distance in meters between two points in latlon
def haversine(lat1, lon1, lat2, lon2):
    radius = 6371000
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    a = math.sin((lat2-lat1)/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin((lon2-lon1)/2)**2
    c = 2 * math.asin(math.sqrt(a))
    return radius * c


# wrapper for multicore process execution
# function: individual function to be applied to each process item
# cores: the number of subprocesses started/CPUs used; this value is reduced in case the number of subprocesses is smaller
# multiargs: a dictionary containing subfunction argument names as keys and lists of arguments to be distributed among the processes as values
# singleargs: all remaining arguments which are invariant among the subprocesses
# important:
# -all multiargs value lists must be of same length, i.e. all argument keys must be explicitly defined for each subprocess
# -all function arguments passed via singleargs must be provided with the full argument name and its value (i.e. argname=argval); default function args are not accepted
# if the processes return anything else than None, this function will return a list of results
# if all processes return None, this function will be of type void
# example:
# def add(x, y, z):
#     return x+y+z
# multicore(add, cores=2, multiargs={"x": [1, 2]}, y=5, z=9)
# -> returns [15, 16]
# multicore(add, cores=2, multiargs={"x": [1, 2], "y": [5, 6]}, z=9)
# -> returns [15, 17]
def multicore(function, cores, multiargs, **singleargs):

    # compare the function arguments with the multi and single arguments and raise errors if mismatches occur
    function_check = inspect.getargspec(function)
    multiargs_check = [x for x in multiargs if x not in function_check[0]]
    singleargs_check = [x for x in singleargs if x not in function_check[0]]
    if len(multiargs_check) > 0:
        raise AttributeError("incompatible multi arguments: {0}".format(", ".join(multiargs_check)))
    if len(singleargs_check) > 0:
        raise AttributeError("incompatible single arguments: {0}".format(", ".join(singleargs_check)))

    # compare the list lengths of the multi arguments and raise errors if they are of different length
    arglengths = list(set([len(x) for x in multiargs]))
    if len(arglengths) > 1:
        raise AttributeError("multi argument lists of different length")

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

    return result

    # # evaluate the return of the processing function; if any value is not None then the whole list of results is returned
    # eval = [x for x in result if x is None]
    # if len(eval) != len(result):
    #     return result


# return the smallest possible data type
def parse_literal(x):
    if isinstance(x, list):
        return [parse_literal(y) for y in x]
    try:
        return int(x)
    except ValueError:
        try:
            return float(x)
        except ValueError:
            return str(x)


# classical queue implementation
class Queue(object):
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


# read processing parameter text files
class ReadPar(object):
    def __init__(self, filename, splits="[:|\t|=\s]", type=""):
        if type == "exe":
            splits = "[\t\n]"
        with open(filename, "r") as infile:
            self.index = []
            for line in infile:
                if not line.startswith("#"):
                    items = filter(None, re.split(splits, line))
                    if len(items) > 1:
                        if len(items) > 2:
                            entry = items[1] if items[2] in ["m", "decimal", "arc-sec", "degrees"] else items[1:]
                            entry = [x for x in items if x not in ["m", "decimal", "arc-sec", "degrees"]]
                            if "".join(entry) == "WGS1984":
                                entry = "WGS84"
                        else:
                            entry = items[1]
                        if type == "exe":
                            items[0] = items[0].replace(" ", "_")
                        setattr(self, items[0], entry)
                        self.index.append(items[0])


# wrapper for subprocess execution including logfile writing and command prompt piping
def run(cmd, outdir=None, logpath=None, inlist=None):
    cmd = [str(x) for x in dissolve(cmd)]
    if outdir is None:
        outdir = os.getcwd()
    if logpath is None:
        log = sp.PIPE
    else:
        index = 1 if cmd[0] in [sys.executable, "Rscript"] else 0
        logfile = os.path.join(logpath, os.path.splitext(cmd[index])[0]+".log")
        log = open(logfile, "a")

    if inlist is None:
        # proc = sp.Popen(cmd, stdin=sp.PIPE, stdout=log, stderr=err, cwd=outdir)
        proc = sp.Popen(cmd, stdin=sp.PIPE, stdout=log, stderr=sp.PIPE, cwd=outdir)
        gammaErrorHandler(proc)
    else:
        out, err = sp.Popen(cmd, stdin=sp.PIPE, stdout=log, stderr=sp.PIPE, cwd=outdir, universal_newlines=True, shell=False).communicate("".join([str(x)+"\n" for x in inlist]))
    # add line for separating log entries of repeated function calls
    if logpath is not None:
        log.write("#####################################################################\n")
        log.close()


# classical stack implementation
# input can be a list, a single value or None (i.e. Stack())
class Stack(object):
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


# union of two lists
def union(a, b):
    return list(set(a) & set(b))


# write parameter textfile
def writer(filename, arguments, strfill=True):
    with open(filename, "w") as out:
        for i in range(0, len(arguments)):
            if strfill:
                out.write(arguments[i][0].replace(" ", "_") + "\t" + arguments[i][1] + "\n")
            else:
                out.write(arguments[i][0] + "\t" + arguments[i][1] + "\n")
