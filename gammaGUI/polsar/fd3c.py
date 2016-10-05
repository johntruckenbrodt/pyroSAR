##############################################################
# Target decomposition based on Freeman-Durden 3-component algorithm
# module of software gammaGUI
# John Truckenbrodt 2015
##############################################################
import os

from ancillary import finder, ReadPar
from envi import hdr
from gammaGUI.auxiliary import grouping
from gamma.util import gamma

path_log = os.path.join(os.getcwd(), "LOG/LAT/")
if not os.path.exists(path_log):
    os.makedirs(path_log)

# create list of scene tuple objects
tuples = grouping()

print "#############################################"
print "creating fd3c decomposition..."

for scene in tuples:
    print scene.basename
    try:
        hh_mli = scene.getTop("HH_(?:slc_|)(?:cal_|)mli$")
        hh_slc = scene.getTop("HH_slc(?:_cal|)$")
        vv_slc = scene.getTop("VV_slc(?:_cal|)$")
        hv_slc = scene.getTop("HV_slc(?:_cal|)$")
        t13 = scene.getTop("t13")
    except ValueError:
        print "...skipped; no appropriate files found"
        continue

    rlks = ReadPar(hh_mli+".par").range_looks
    azlks = ReadPar(hh_mli+".par").azimuth_looks
    samples = ReadPar(scene.HH_slc+".par").range_samples
    base = os.path.basename(hh_slc.replace("HH_", ""))

    gamma(["FD3C_DEC", hh_slc, hv_slc, vv_slc, t13, samples, base, rlks, azlks], os.path.dirname(t13), path_log)
    for tag in ["_fdd_pd", "_fdd_ps", "_fdd_pv"]:
        hdr(hh_mli+".par", os.path.join(os.path.dirname(t13), base)+tag+".hdr")

# rename files to consistent pattern
for pattern in ["*.fdd*"]:
    for filename in finder(os.getcwd(), [pattern]):
        os.rename(filename, filename.replace(pattern.strip("*"), pattern.strip("*").replace(".", "_")))

print "...done"
