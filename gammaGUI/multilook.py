##############################################################
# SAR image multilooking of slc and rslc files
# module of software gammaGUI
# Elias Wolf, John Truckenbrodt 2014-2015
##############################################################
import sys

import os

from gamma.util import ISPPar, Spacing, gamma
from envi import hdr

# receive input file
data = sys.argv[2]
meta = data + ".par"

print data

# read parameter file and compute multilooking parameters
par = ISPPar(meta)
mlk = Spacing(par, sys.argv[1])

# define (and create) directories for logfile
path_log = os.path.join(os.getcwd(), "LOG/ISP/")
if not os.path.exists(path_log):
    os.makedirs(path_log)

print "slc range pixel spacing (slant, ground):", par.range_pixel_spacing, mlk.groundRangePS
print "slc azimuth pixel spacing:", par.azimuth_pixel_spacing
print "number of looks looks (range, azimuth):", mlk.rlks, mlk.azlks
print "mli range pixel spacing (slant, ground):", int(mlk.rlks)*float(par.range_pixel_spacing), int(mlk.rlks)*mlk.groundRangePS
print "mli azimuth pixel spacing:", int(mlk.azlks)*float(par.azimuth_pixel_spacing)
print "-----------"

# concatenate output names
out_data = data+"_mli"
out_meta = out_data + ".par"

# set scaling factor
scale = 0.000001 if "ERS" in par.sensor else 1.0

# perform gamma command
gamma(["multi_look", data, meta, out_data, out_meta, mlk.rlks, mlk.azlks, "-", "-", scale], os.getcwd(), path_log)
hdr(out_meta)
