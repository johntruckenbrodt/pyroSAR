##############################################################
# Interferogram simulation and differential interferogram computation
# module of software gammaGUI
# John Truckenbrodt 2015
##############################################################

import os
import re

from ancillary import finder, ReadPar
from gamma.util import gamma

# read parameter file
par = ReadPar(os.path.join(os.getcwd(), "PAR/baseline.par"))

# define (and create) directories for processing results and logfile
path_log = os.path.join(os.getcwd(), "LOG/ISP/")
path_out = os.path.join(os.getcwd(), "ISP/")
for path in [path_log, path_out]:
    if not os.path.exists(path):
        os.mkdir(path)

interferograms = finder(path_out, ["*int"])

if len(interferograms) > 0:
    print "#############################################"
    print "interferogram flattening started..."
    for name_int in interferograms:
        print os.path.basename(name_int)

        # retrieve full name of primary and secondary SLC files
        scenes = re.findall("[A-Z0-9_]{10}[0-9T]{15}_[HV]{2}_slc(?:_cal|)", name_int)
        primary = finder(os.getcwd(), [scenes[0]+"$"], regex=True)[0]
        # secondary = finder(os.getcwd(), [scenes[1]+"_reg"], regex=True)[0]
        secondary = name_int[:-3]+"reg"

        # collect geocoding lookup tables
        lut_list = [x for x in finder(os.getcwd(), [re.findall("[A-Z0-9_]{9,10}[0-9T]{15}", name_int)[0]], regex=True) if "map_to_rdc" in x]

        if len(lut_list) > 0:
            lut_fine = lut_list[0]
        else:
            print "geocoding lookup table missing"
            continue

        cc_sm = name_int[:-3]+"cc_sm"
        dem_map = lut_fine[:-10]+"dem"
        dem_rdc = dem_map+"_rdc"
        diff_int = name_int+"_diff"
        diff_int_sm = diff_int+"_sm"
        name_diffpar = diff_int[:-3]+"diff_par"
        name_off = name_int[:-3]+"off"
        name_base_init = name_int[:-3]+"base_init"
        name_base_refine = name_int[:-3]+"base_refine"
        name_base_res = name_int[:-3]+"base_res"
        par_dem = ReadPar(dem_map+".par")
        par_mli = ReadPar(primary+"_mli.par")
        ph_sim = name_int[:-3]+"ph_sim"

        if not os.path.isfile(dem_map):
            print "DEM", dem_map, "missing"
            continue

        if not os.path.isfile(dem_rdc):
            print "...transforming DEM to range-doppler coordinates"
            gamma(["geocode", lut_fine, dem_map, par_dem.width, dem_rdc, par_mli.range_samples, par_mli.azimuth_lines, "2"], os.path.dirname(dem_map), path_log)
        else:
            print "...found DEM in range-doppler coordinates"

        print "...initial DEM interferogram simulation"
        gamma(["phase_sim", primary + ".par", name_off, name_base_init, dem_rdc, ph_sim], os.path.dirname(name_int), path_log)

        print "...initial differential interferogram generation"
        gamma(["SLC_diff_intf", primary, secondary, primary + ".par", secondary + ".par", name_off, ph_sim, diff_int,
             par_mli.range_looks, par_mli.azimuth_looks], os.path.dirname(name_int), path_log)

        print "...refining baseline"
        gamma(["base_init", primary + ".par", secondary + ".par", name_off, diff_int, name_base_res, "4"], path_out, path_log)
        gamma(["base_add", name_base_init, name_base_res, name_base_refine, "1"], path_out, path_log)

        print "...refined DEM interferogram simulation"
        gamma(["phase_sim", primary + ".par", name_off, name_base_refine, dem_rdc, ph_sim], os.path.dirname(name_int), path_log)

        print "...refined differential interferogram generation"
        gamma(["SLC_diff_intf", primary, secondary, primary + ".par", secondary + ".par", name_off, ph_sim, diff_int,
             par_mli.range_looks, par_mli.azimuth_looks], os.path.dirname(name_int), path_log)
        print "----------"

    print "...done"
else:
    print "#############################################"
    print "no interferograms found"
