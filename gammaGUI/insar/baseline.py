##############################################################
# Interferometric Baseline Estimation
# module of software gammaGUI
# John Truckenbrodt 2015
##############################################################

"""
input: SLC and RSLC file (to be passed by executing the script, i.e. python baseline.py SLC RSLC
The following tasks are performed by executing this script:
-reading of a parameter file baseline.par
--see object par for necessary values; file is automatically created by starting the script via the GUI
-if necessary, creation of output and logfile directories
-check whether corresponding coregistration offset file and interferogram exist
-execution of initial baseline estimation and calculation of perpendicular and parallel components
"""

import re
import os

from ancillary import finder, ReadPar
from gamma.util import gamma

# read parameter file
par = ReadPar(os.path.join(os.getcwd(), "PAR/baseline.par"), type="exe")

# define (and create) directories for processing results and logfile
path_log = os.path.join(os.getcwd(), "LOG/ISP/")
path_out = os.path.join(os.getcwd(), "ISP/")
for path in [path_log, path_out]:
    if not os.path.exists(path):
        os.mkdir(path)

interferograms = finder(path_out, ["*int"])

if len(interferograms) > 0:
    print "#############################################"
    print "baseline estimation started..."
    for name_int in interferograms:
        name_off = name_int[:-3]+"off"
        name_base = name_int[:-3]+"base_init"
        print os.path.basename(name_base)

        slc = finder(os.getcwd(), [re.findall("[A-Z0-9_]{9,10}[0-9T]{15}_[HV]{2}_slc(?:_cal|)", name_int)[0]])[0]
        # rslc = finder(os.getcwd(), [re.findall("[A-Z0-9_]{10}[0-9T]{15}_[HV]{2}_slc(?:_cal)", name_int)[1]+"_reg"])[0]
        rslc = name_int[:-3]+"reg"

        gamma(["base_init", slc + ".par", rslc + ".par", name_off, name_int, name_base, par.method_flag, par.nrfft, par.nazfft, par.r_samp, par.az_line], path_out, path_log)

    print "...done"
else:
    print "#############################################"
    print "no interferograms found"
