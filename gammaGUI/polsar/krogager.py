##############################################################
# Calculate Helix and Diplane composition from RR and LL circular components
# Multilook RL circular component (Sphere component)
# Sphere, Helix and Diplane represent the Krogager decomposition
# module of software gammaGUI
# John Truckenbrodt 2015-16
##############################################################

import os

from ancillary import ReadPar
from envi import hdr
from gammaGUI.auxiliary import grouping
import gamma

path_log = os.path.join(os.getcwd(), "LOG/LAT/")

# create list of scene tuple objects
tuples = grouping()

print "#############################################"
print "creating krogager decomposition..."

for scene in tuples:
    print scene.basename
    try:
        hh_mli = scene.getTop("HH_(?:slc_|)(?:cal_|)mli$")
        hh_slc = scene.getTop("HH_slc(?:_cal|)$")
        rl = scene.getTop("rl$")
        ll = scene.getTop("ll$")
        rr = scene.getTop("rr$")
    except (ValueError, AttributeError):
        print "...skipped; no appropriate files found"
        continue

    path_out = os.path.dirname(rl)
    mlipar = hh_mli+".par"
    rlks = ReadPar(mlipar).range_looks
    azlks = ReadPar(mlipar).azimuth_looks

    base = os.path.basename(hh_slc.replace("HH_", ""))

    # define placeholder parameter file names, which are required by the single commands but will not be needed in the future and are thus deleted once gamma commands have been executed
    dummypar1 = os.path.join(path_out, base+"_dummy1")
    dummypar2 = os.path.join(path_out, base+"_dummy2")
    gamma.process(["multi_look", rl, hh_slc + ".par", base + "_sphere", dummypar1, rlks, azlks], path_out, path_log)
    try:
        gamma.process(["diplane_helix", ll, rr, hh_slc + ".par", base + "_diplane", base + "_helix", dummypar2, rlks, azlks, "-", "-", "-"], path_out, path_log)
    # catch strange behaviour of diplane_helix command to throw an error; by execution in the shell, no error is shown, via python an error is forwarded
    except RuntimeError:
        pass
    os.remove(dummypar1)
    os.remove(dummypar2)

    for tag in ["_sphere", "_helix", "_diplane"]:
        hdr(mlipar, os.path.join(path_out, base)+tag+".hdr")

print "...done"
