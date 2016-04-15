##############################################################
# Interferogram Generation
# module of software gammaGUI
# John Truckenbrodt 2015
##############################################################

"""
input: SLC and RSLC file (to be passed by executing the script, i.e. python interferogram.py SLC RSLC
The following tasks are performed by executing this script:
-reading of a parameter file interferogram.par
--see object par for necessary values; file is automatically created by starting the script via the GUI
-if necessary, creation of output and logfile directories
-check whether corresponding coregistration offset file  exists
-interferogram generation
"""

import re
import os

from ancillary import finder, ReadPar, run

par = ReadPar(os.path.join(os.getcwd(), "PAR/interferogram.par"), type="exe")

# define (and create) directories for processing results and logfile
path_log = os.path.join(os.getcwd(), "LOG/ISP/")
path_out = os.path.join(os.getcwd(), "ISP/")
for path in [path_log, path_out]:
    if not os.path.exists(path):
        os.mkdir(path)

offsets = finder(path_out, ["*off"])

if len(offsets) > 0:
    print "#############################################"
    print "interferogram generation started..."
    for name_off in offsets:
        name_int = name_off[:-3]+"int"
        print os.path.basename(name_int)

        slc = finder(os.getcwd(), [re.findall("[A-Z0-9_]{9,10}[0-9T]{15}_[HV]{2}_slc(?:_cal|)", name_int)[0]])[0]
        # rslc = finder(os.getcwd(), [re.findall("[A-Z0-9_]{10}[0-9T]{15}_[HV]{2}_slc(?:_cal)", name_int)[1]+"_reg"])[0]
        rslc = name_off[:-3]+"reg"

        # read the master mli parameter file for multilooking factors
        try:
            par_mli = ReadPar(slc+"_mli.par")
        except IOError:
            print "MLI for the primary SLC missing"

        run(["SLC_intf", slc, rslc, slc + ".par", rslc + ".par", name_off, name_int, par_mli.range_looks, par_mli.azimuth_looks, "0", "-", par.sps_flg, par.azf_flg,
             par.rp1_flg, par.rp2_flg], path_out, path_log)

    print "...done"
else:
    print "#############################################"
    print "no coregistration offset files found"
