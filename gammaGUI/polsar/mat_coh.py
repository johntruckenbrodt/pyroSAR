##############################################################
# Calculate coherence matrix T elements from Pauli decomposition alpha, beta, gamma
# module of software gammaGUI
# John Truckenbrodt 2015-16
##############################################################
import os
import re

from ancillary import finder, ReadPar
from gammaGUI.auxiliary import grouping
from gamma.util import gamma

path_log = os.path.join(os.getcwd(), "LOG/LAT/")
if not os.path.exists(path_log):
    os.makedirs(path_log)

par = ReadPar(os.path.join(os.getcwd(), "PAR/mat_coh.par"))

tuples = grouping()

print "#############################################"
print "creating coherence matrices..."

for scene in tuples:
    print scene.basename
    try:
        hh_mli = scene.getTop("HH_(?:slc_|)(?:cal_|)mli$")
        pauli_alpha = scene.getTop("pauli_alpha$")
        pauli_beta = scene.getTop("pauli_beta$")
        pauli_gamma = scene.getTop("pauli_gamma$")
    except ValueError:
        print "...skipped; no appropriate files found"
        continue

    rlks = ReadPar(hh_mli+".par").range_looks
    azlks = ReadPar(hh_mli+".par").azimuth_looks
    base = pauli_alpha.replace("_pauli_alpha", "")
    gamma(["polcoh", pauli_alpha, pauli_beta, pauli_gamma, pauli_alpha + ".par", pauli_beta + ".par", pauli_gamma + ".par",
         base, base + "_mat_coh.par", rlks, azlks], os.path.dirname(pauli_alpha), path_log)

# rename files to consistent pattern
for filename in finder(os.getcwd(), ["*.t*"]):
    os.rename(filename, filename.replace(".t", "_t"))

print "...done"
