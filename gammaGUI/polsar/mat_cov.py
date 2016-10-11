##############################################################
# Calculate covariance matrix C elements from HH, HV, and VV SLC data
# module of software gammaGUI
# John Truckenbrodt 2015-16
##############################################################

import os
import re

from ancillary import finder, ReadPar
from gammaGUI.auxiliary import grouping
import gamma

path_log = os.path.join(os.getcwd(), "LOG/LAT/")
if not os.path.exists(path_log):
    os.makedirs(path_log)

par = ReadPar(os.path.join(os.getcwd(), "PAR/mat_cov.par"))

# create list of scene tuple objects
tuples = grouping()

print "#############################################"
print "creating covariance matrices..."

for scene in tuples:
    print scene.basename
    try:
        hh_slc = scene.getTop("HH_slc(?:cal_|)$")
        vv_slc = scene.getTop("VV_slc(?:cal_|)$")
        hv_slc = scene.getTop("HV_slc(?:cal_|)$")
        hh_mli = scene.getTop("HH_(?:slc_|)(?:cal_|)mli$")
    except ValueError:
        print "...skipped; no appropriate files found"
        continue

    rlks = ReadPar(hh_mli+".par").range_looks
    azlks = ReadPar(hh_mli+".par").azimuth_looks
    path_out = os.path.join(os.path.dirname(hh_slc), "POL/")
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    base = os.path.basename(hh_slc.replace("HH_", ""))
    gamma.process(["polcovar", hh_slc, hv_slc, vv_slc, hh_slc + ".par", hv_slc + ".par", vv_slc + ".par", base,
             base + "_mat_cov.par", rlks, azlks], path_out, path_log)

# rename files to consistent pattern
for filename in finder(os.getcwd(), ["*.c*"]):
    os.rename(filename, filename.replace(".c", "_c"))

print "...done"
