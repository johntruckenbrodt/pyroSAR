##############################################################
# Huynen decomposition to generate equivalent single target coherency matrix values
# module of software gammaGUI
# John Truckenbrodt 2015-16
##############################################################

import os

from ancillary import finder, ReadPar
from gammaGUI.auxiliary import grouping
import gamma

path_log = os.path.join(os.getcwd(), "LOG/LAT/")
if not os.path.exists(path_log):
    os.makedirs(path_log)

# create list of scene tuple objects
tuples = grouping()

print "#############################################"
print "creating huynen decomposition..."
counter = 0
for scene in tuples:
    print scene.basename
    try:
        hh_slc = scene.getTop("HH_slc(?:cal_|)$")
        vv_slc = scene.getTop("VV_slc(?:cal_|)$")
        hv_slc = scene.getTop("HV_slc(?:cal_|)$")
        t11 = scene.getTop("t11$")
        t12 = scene.getTop("t12$")
        t13 = scene.getTop("t13$")
    except ValueError:
        print "...skipped; no appropriate files found"
        continue

    samples = ReadPar(hh_slc+".par").range_samples
    base = os.path.basename(hh_slc).replace("HH_", "")

    for i in range(1, 4):
        gamma.process(["HUYNEN_DEC", hh_slc, hv_slc, vv_slc, t11, t12, t13, samples, base, i], os.path.dirname(t11), path_log)

else:
    # rename files to consistent pattern
    for pattern in ["*.t*", "*.im", "*.re"]:
        for filename in finder(os.getcwd(), [pattern]):
            os.rename(filename, filename.replace(pattern.strip("*"), pattern.strip("*").replace(".", "_")))

print "...done"
