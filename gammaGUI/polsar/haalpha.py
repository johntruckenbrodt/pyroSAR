##############################################################
# Calculate H/A/alpha (Entropy, Anisotropy, and alpha) decomposition from the 3D Pauli feature vector
# module of software gammaGUI
# John Truckenbrodt 2015-16
##############################################################

import os

from ancillary import ReadPar
from envi import hdr
from gammaGUI.auxiliary import grouping
from gamma.util import gamma

path_log = os.path.join(os.getcwd(), "LOG/LAT/")

# create list of scene tuple objects
tuples = grouping()

print "#############################################"
print "creating entropy-anisotropy-alpha decomposition..."

for scene in tuples:
    print scene.basename
    try:
        hh_mli = scene.getTop("HH_(?:slc_|)(?:cal_|)mli$")
        hh_slc = scene.getTop("HH_slc(?:_cal|)$")
        pauli_alpha = scene.getTop("pauli_alpha$")
        pauli_beta = scene.getTop("pauli_beta$")
        pauli_gamma = scene.getTop("pauli_gamma$")
    except (ValueError, AttributeError):
        print "...skipped; no appropriate files found"
        continue

    path_out = os.path.dirname(pauli_alpha)
    mlipar = hh_mli+".par"
    rlks = ReadPar(mlipar).range_looks
    azlks = ReadPar(mlipar).azimuth_looks

    base = os.path.basename(hh_slc.replace("HH_", ""))

    gamma(["haalpha", pauli_alpha, pauli_beta, pauli_gamma, hh_slc + ".par", base + "_cpd_A", base + "_cpd_alpha",
         base + "_cpd_H", base + "_cpd_l1", base + "_cpd_l2", base + "_cpd_l3", mlipar, rlks, azlks], path_out, path_log)

    for tag in ["_cpd_A", "_cpd_alpha", "_cpd_H"]:
        hdr(mlipar, os.path.join(path_out, base)+tag+".hdr")

print "...done"
