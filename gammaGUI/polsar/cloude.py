##############################################################
# Cloude target decomposition from elements of scattering and coherency matrix
# module of software gammaGUI
# John Truckenbrodt 2015-16
##############################################################

import os

from ancillary import finder, ReadPar
from envi import hdr, HDRobject
from gammaGUI.auxiliary import grouping
from gamma.util import gamma

path_log = os.path.join(os.getcwd(), "LOG/LAT/")
if not os.path.exists(path_log):
    os.makedirs(path_log)

par = ReadPar(os.path.join(os.getcwd(), "PAR/cloude.par"))

# create list of scene tuple objects
tuples = grouping()

print "#############################################"
print "creating cloude decomposition..."

for scene in tuples:
    print scene.basename
    try:
        hh_mli = scene.getTop("HH_(?:slc_|)(?:cal_|)mli$")
        hh_slc = scene.getTop("HH_slc(?:_cal|)$")
        vv_slc = scene.getTop("VV_slc(?:_cal|)$")
        hv_slc = scene.getTop("HV_slc(?:_cal|)$")
        t12 = scene.getTop("t12$")
        t13 = scene.getTop("t13$")
    except ValueError:
        print "...skipped; no appropriate files found"
        continue

    rlks = ReadPar(hh_mli+".par").range_looks
    azlks = ReadPar(hh_mli+".par").azimuth_looks
    samples = ReadPar(hh_slc+".par").range_samples

    base = os.path.basename(hh_slc.replace("HH_", ""))
    gamma(["CLOUDE_DEC", hh_slc, hv_slc, vv_slc, t12, t13, samples, base, rlks, azlks], os.path.dirname(t12), path_log)

    # create envi header files (note: number of lines must be reduced by 1 for import into envi)
    header = HDRobject(hh_mli+".par")
    header.lines -= 1
    for i in range(1, 4):
        hdr(header, os.path.join(os.path.dirname(t12), base)+"_ctd_"+str(i)+"_mag.hdr")

# rename files to consistent pattern
for pattern in ["*.ctd*", "*.mag", "*.pha"]:
    for filename in finder(os.getcwd(), [pattern]):
        os.rename(filename, filename.replace(pattern.strip("*"), pattern.strip("*").replace(".", "_")))

print "...done"
