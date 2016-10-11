#!/usr/bin/env python
##############################################################
# SLC calibration
# module of software gammaGUI
# John Truckenbrodt 2015
# last update 2015-11-19
##############################################################

"""
The script is intended for initial calibration of SLC images imported into the gammaGUI.
The central GAMMA command is radcal_SLC, which is parametrized specific to individual sensors.
The gammaGUI working directory is first searched for SLC files (suffix "_slc"). From the name of the file the sensor is matched and
calibration is performed. If the option "replace" was set in the GUI dialog, the original SLC files are deleted.
"""
import os
from ancillary import finder, ReadPar
import gamma

# read parameter file
par = ReadPar(os.path.join(os.getcwd(), "PAR/cal_slc.par"))

# define (and create) directories for processing results and logfile
path_log = os.path.join(os.getcwd(), "LOG/GEO/")
path_out = os.path.join(os.getcwd(), "ISP/")
for path in [path_log, path_out]:
    if not os.path.exists(path):
        os.makedirs(path)

list_K_dB = {"PSR1": -115.0, "PSR2": -115.0}

list_slc = finder(os.getcwd(), ["*_slc"])

rejected = []

if len(list_slc) > 0:
    print "#############################################"
    print "calibration started..."

    for name_slc in list_slc:
        sensor = os.path.basename(name_slc).split("_")[0]
        if sensor in list_K_dB:
            print name_slc
            K_dB = list_K_dB[sensor]
            name_cslc = name_slc+"_cal"
            gamma.process(["radcal_SLC", name_slc, name_slc + ".par", name_cslc, name_cslc + ".par", "-", "-", "-", "-", "-", "-", K_dB], path_out, path_log)
            if par.replace == "True":
                os.remove(name_slc)
                os.remove(name_slc+".par")
        else:
            if sensor not in rejected:
                print "calibration for sensor {} not implemented (yet)".format(sensor)
                rejected.append(sensor)
            continue
    print "...done"
else:
    print "#############################################"
    print "no SLCs found"
