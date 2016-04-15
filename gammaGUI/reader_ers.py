##############################################################
# ERS-1/2-specific SLC data import and calibration
# module of software gammaGUI
# John Truckenbrodt 2015
##############################################################

"""
The following tasks are performed by executing this script:
-scan the defined import directory for the ERS-specific leader file LEA_01.001
-convert the data to gamma format
-create scene folders and rename the newly created files based on the sensor, frame id and acquisition time stamp
-scan the leader file for meta information relevant for calibration
-select a antenna gain correction lookup file from the extracted information; the files are stored in a subfolder CAL of this software package
-calibrate the SLC files based on the extracted meta information
-remove the uncalibrated SLCs
"""

import os
import sys
import time
import pyroSAR
from gamma.util import ISPPar
from ancillary import finder, run

orbit_correct = True if sys.argv[-1] == "True" else False

# path to delft orbit files
path_delft = "/pvdata2/john/ancillary/ERS/ORBIT/delft"

# path to antenna correction files
path_cal = "/pvdata2/john/ancillary/ERS/CAL/ERS_antenna"

# define (and create) directory for logfile
path_log = os.path.join(os.getcwd(), "LOG/IMP/")
if not os.path.exists(path_log):
    os.makedirs(path_log)

scenes = [os.path.dirname(x) for x in finder(sys.argv[1], ["*LEA_01.001"])]

if len(scenes) == 0:
    raise IOError("no appropriate file found")

for scene in scenes:
    print "----------"

    # read scene metadata
    id = pyroSAR.ERS(scene)

    tempname = os.path.join(os.getcwd(), "temp")
    print "importing..."
    run(["par_ESA_ERS", "LEA_01.001", tempname+".par", "DAT_01.001", tempname], scene, path_log, [""])
    par = ISPPar(tempname+".par")

    date = "".join([format(int(x), "02d") for x in par.date[0:3]])
    timestamp = date+"T"+time.strftime("%H%M%S", time.gmtime(round(par.center_time)))
    outname = par.sensor+"_"+id.frame+"_"+timestamp+"_VV_slc"
    path_out = os.path.join(os.getcwd(), outname[:-7])
    if not os.path.exists(path_out):
        print outname
        os.makedirs(path_out)
        os.rename(tempname, os.path.join(path_out, outname))
        os.rename(tempname+".par", os.path.join(path_out, outname+".par"))
    else:
        print "scene", outname, "already imported; removing temporary files"
        os.remove(tempname)
        os.remove(tempname+".par")
    outname = os.path.join(path_out, outname)

    antenna = "antenna_"+par.sensor

    print "sensor:", par.sensor
    print "date:", int(date)
    print "processing facility:", id.proc_fac
    print "processing system:", id.proc_sys
    print "processing version:", id.proc_vrs
    print "antenna correction flag:", id.antenna_flag
    print "antenna correction file:", antenna
    print "calibration constant K [dB]:", id.cal

    sc_db = 59.61 if par.sensor == "ERS1" else 60
    antenna_corr = 1 if id.antenna_flag == "0" else 0
    antenna = os.path.join(path_cal, antenna)

    if orbit_correct:
        print "...correcting orbits"
        try:
            if par.sensor == "ERS1":
                path_delft_target = os.path.join(path_delft, par.sensor, "dgm-e04" if int(date) <= 19960601 else "dgm-e04.fd")
            else:
                path_delft_target = os.path.join(path_delft, par.sensor, "dgm-e04")
            run(["DELFT_vec2", outname+".par", path_delft_target], path_out, path_log)
        except:
            pass

    # perform radiometric calibration
    print "...calibrating"
    run(["radcal_SLC", outname, outname+".par", outname+"_cal", outname+"_cal.par", 4, antenna, 1, antenna_corr, 1, sc_db, id.cal], path_out, path_log)

    os.remove(os.path.join(path_out, outname))
    os.remove(os.path.join(path_out, outname)+".par")
