##############################################################
# conversion to circular polarization basis
# module of software gammaGUI
# John Truckenbrodt 2015-16
##############################################################
import os

from ancillary import ReadPar
from gammaGUI.auxiliary import grouping
import gamma

path_log = os.path.join(os.getcwd(), "LOG/LAT/")
if not os.path.exists(path_log):
    os.makedirs(path_log)

par = ReadPar(os.path.join(os.getcwd(), "PAR/circular.par"))

# create list of scene tuple objects
tuples = grouping()

print "#############################################"
print "transforming scenes..."

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

    base = os.path.basename(hh_slc.replace("HH_", ""))
    samples = ReadPar(hh_slc+".par").range_samples

    gamma.process(["lin_comb_cpx", 3, hh_slc, vv_slc, hv_slc, par.constant_r, par.constant_i, par.factorHH_r, par.factorHH_i, par.factorVV_r, par.factorVV_i,
             par.factorHV_r, par.factorHV_i, base + "_rr", samples, "", "", par.pixav_x, par.pixav_y, 1], path_out, path_log)

    gamma.process(["lin_comb_cpx", 3, vv_slc, hh_slc, hv_slc, par.constant_r, par.constant_i, par.factorVV_r, par.factorVV_i, par.factorHH_r, par.factorHH_i,
             par.factorHV_r, par.factorHV_i, base + "_ll", samples, "", "", par.pixav_x, par.pixav_y, 1], path_out, path_log)

    gamma.process(["lin_comb_cpx", 2, hh_slc, vv_slc, par.constant_r, par.constant_i, par.factorHH_r, par.factorHH_i, par.factorVV_r, par.factorVV_i,
             base + "_rl", samples, "", "", par.pixav_x, par.pixav_y, 1], path_out, path_log)

print "...done"
