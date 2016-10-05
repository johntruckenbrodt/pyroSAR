##############################################################
# Calculate Pauli polarimetric decomposition from HH, VV, and HV SLC images
# module of software gammaGUI
# John Truckenbrodt 2015-16
##############################################################

import os

from ancillary import finder
from gammaGUI.auxiliary import grouping
from gamma.util import gamma

path_log = os.path.join(os.getcwd(), "LOG/LAT/")
if not os.path.exists(path_log):
    os.makedirs(path_log)

# create list of scene tuple objects
# all images following the defined patterns within the same folder (i.e. the same acquisition) will be grouped together
tuples = grouping()

print "#############################################"
print "creating pauli decomposition..."

for scene in tuples:
    print scene.basename
    try:
        hh_slc = scene.getTop("HH_slc(?:_cal|)$")
        vv_slc = scene.getTop("VV_slc(?:_cal|)$")
        hv_slc = scene.getTop("HV_slc(?:_cal|)$")
    except ValueError:
        print "...skipped; no appropriate files found"
        continue

    path_out = os.path.join(os.path.dirname(hh_slc), "POL/")
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    name_out = os.path.join(path_out, os.path.basename(hh_slc.replace("HH_", "")+"_pauli"))
    gamma(["pauli", hh_slc, vv_slc, hv_slc, hh_slc + ".par", vv_slc + ".par", hv_slc + ".par", name_out], os.getcwd(), path_log)

# rename files to consistent pattern
for filename in finder(os.getcwd(), ["*.slc*"]):
    os.rename(filename, filename.replace(".slc", ""))

print "...done"
