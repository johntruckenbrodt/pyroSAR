##############################################################
# preparation of dem data for use in gamma
# module of software gammaGUI
# John Truckenbrodt 2014-15
##############################################################

"""
The following tasks are performed by executing this script:
-read the header file of the provided envi-formatted (!) dem
-if necessary:
--convert data type from 8Bit unsigned integer or 16Bit signed integer to 32Bit float
--swap the bytes of the file to big endian (in place!)
--adjust the envi header file according to prior conversion steps
--create a gamma parameter file from the header information
"""

import sys

import os
import subprocess as sp

from ancillary import run
from envi import hdr, HDRobject

# define (and create) directories for processing results and logfile
path_dem = os.path.join(os.getcwd(), "DEM/")
path_log = os.path.join(os.getcwd(), "LOG/GEO/")
for path in [path_log, path_dem]:
    if not os.path.exists(path_log):
        os.makedirs(path_log)

# get dem name from script call
dem = sys.argv[1]

# load header file
try:
    header = HDRobject(os.path.splitext(dem)[0]+".hdr")
except IOError:
    print "header file not found"

if "UTM" not in header.map_info:
    raise ValueError("currently only UTM projection is supported")

if header.data_type not in ["1", "2", "4"]:
    raise ValueError("data type not supported; provide int1u, int2s or float4")

process = 0

# swap byte order from LSF (little endian) to MSF (big endian)
if header.byte_order == "0":
    print "swapping bytes"
    sp.check_call(["swap_bytes", dem, dem+"_s", {"1": "4", "2": "2", "4": "4"}[header.data_type]], stdout=sp.PIPE)
    os.remove(dem)
    os.rename(dem+"_s", dem)
    # adjust header file
    header.byte_order = "1"
    process += 1

# convert file to 32Bit float if data type is 8Bit unsigned integer or 16Bit signed integer
if header.data_type in ["1", "2"]:
    print "converting data type"
    sp.check_call([{"1": "uchar2float", "2": "short2float"}[header.data_type], dem, dem+"_f"], stdout=sp.PIPE)
    os.remove(dem)
    os.rename(dem+"_f", dem)
    header.data_type = "4"
    process += 1

if process > 0:
    hdr(header)

# create gamma parameter file
if not os.path.isfile(dem + ".par"):
    print "creating parameter file"
    false_northing = "0" if "North" in header.map_info else "10000000"
    posting = "-"+header.map_info[6]+" "+header.map_info[5]
    topleft = header.map_info[4]+" "+header.map_info[3]
    dempar = ["UTM", "WGS84", "1", header.map_info[7], false_northing, os.path.basename(dem), "", "", "", header.samples, header.lines, posting, topleft]
    run(["create_dem_par", dem + ".par"], path_dem, path_log, inlist=dempar)
    process += 1

if process == 0:
    print "nothing to be done"
