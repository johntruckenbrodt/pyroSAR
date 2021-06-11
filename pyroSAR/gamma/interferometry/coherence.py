# -*- coding: utf-8 -*-
#created by Marcel Felix & Friederike Metz
import pyroSAR as ps
from pyroSAR import gamma
import os
import math
#import spatialist
from time import perf_counter 
from pyroSAR.gamma.api import diff,isp,disp,lat
import auxmodules



# Initialisation #
auxmodules.Script_init()

# Path Variables #
project_path, data_path = auxmodules.Path_init()

## DEM Preparation ##
dem_type = auxmodules.DEM_Select()
if dem_type == "srtm":
    dem_name = auxmodules.SRTMDEM(data_path)
if dem_type == "palsar":
    dem_name = auxmodules.PALSARDEM(data_path)
dem_res = auxmodules.DEM_res()

# Processing parameters #
n_step, ccthresh, frac_thresh, phstdev_thresh, boxmin, boxmax = auxmodules.Processing()


## Last Hints ##
auxmodules.LastHints()
# timer start #
t_start = perf_counter()

### Scene Identification ###
os.chdir(data_path)
scene1, scene2 = auxmodules.Scene_ident(data_path)
scene1_id = ps.identify(scene1)
scene2_id = ps.identify(os.path.join(data_path,scene2))

## unzip ##
scene1_id.unpack(data_path,overwrite=True)
scene2_id.unpack(data_path,overwrite=True)

#auxmodules.Scene_unpack(data_path,scene1,scene2)

## new identification ##
scene1, scene2 = auxmodules.Scene_ident2(data_path)
scene1_id = ps.identify(os.path.join(data_path,scene1))
scene2_id = ps.identify(os.path.join(data_path,scene2))

## Set variables ##
date1 = scene1_id.start[0:8]
date2 = scene2_id.start[0:8]
dir1 = (os.path.join(data_path,scene1))
dir2 = (os.path.join(data_path,scene2))

## UNPACKING SWATHS AND CREATING PARAMETER AND SLC FILES ##
print("""### UNPACKING SWATHS AND CREATING PARAMETER AND SLC FILES ###\n
          ...""")
pol = "vv"
lst = ["iw1","iw2","iw3"]
def find_slc(pol,lst,directory):
    for file in os.listdir(os.path.join(directory,"measurement")):
        if (pol in file) and (lst in file):
            return file
def find_ann(pol,lst,directory):
    for file in os.listdir(os.path.join(directory,"annotation")):
        if (pol in file) and (lst in file):
            return file
def find_cal(pol,lst,directory):
    for file in os.listdir(os.path.join(directory,"annotation/calibration")):
        if (pol in file) and (lst in file) and ("calibration" in file):
            return file    
def find_noi(pol,lst,directory):
    for file in os.listdir(os.path.join(directory,"annotation/calibration")):
        if (pol in file) and (lst in file) and ("noise" in file):
            return file                              

for i in lst:
    slc1 = os.path.join(dir1,"measurement/")+find_slc(pol,i,dir1)
    ann1 = os.path.join(dir1,"annotation/")+find_ann(pol,i,dir1)
    cal1 = os.path.join(dir1,"annotation/calibration/")+find_cal(pol,i,dir1)
    noi1 = os.path.join(dir1,"annotation/calibration/")+find_noi(pol,i,dir1)
    isp.par_S1_SLC(slc1, ann1, cal1, noi1, 
                   date1+"_"+i+"_"+pol+"_slc.par", 
                   date1+"_"+i+"_"+pol+"_slc", 
                   date1+"_"+i+"_"+pol+"_slc.tops_par", 
                   dtype='-', sc_dB='-',noise_pwr='-', 
                   logpath=None, outdir=None, shellscript=None)
    slc2 = os.path.join(dir2,"measurement/")+find_slc(pol,i,dir2)
    ann2 = os.path.join(dir2,"annotation/")+find_ann(pol,i,dir2)
    cal2 = os.path.join(dir2,"annotation/calibration/")+find_cal(pol,i,dir2)
    noi2 = os.path.join(dir2,"annotation/calibration/")+find_noi(pol,i,dir2)
    isp.par_S1_SLC(slc2, ann2, cal2, noi2, 
                                 date2+"_"+i+"_"+pol+"_slc.par", 
                                 date2+"_"+i+"_"+pol+"_slc", 
                                 date2+"_"+i+"_"+pol+"_slc.tops_par", 
                                 dtype='-', sc_dB='-',noise_pwr='-', 
                                 logpath=None, outdir=None, shellscript=None)


## Burst Rectangles ##
print("""### BURST RECTANGLES ###\n
          ...""")
isp.ScanSAR_burst_corners(SLC_par=date1+"_iw1_"+pol+"_slc.par",
                          TOPS_par=date1+"_iw1_"+pol+"_slc.tops_par",
                          KML='-', logpath=None, outdir=None, shellscript=None)

## Creating the tab-files ##
tab1 = open("slc1_tab.txt", "w")
tab1.write(date1+"_iw1_"+pol+"_slc " +date1+"_iw1_"+pol+"_slc.par " +date1+"_iw1_"+pol+"_slc.tops_par\n")
tab1.write(date1+"_iw2_"+pol+"_slc " +date1+"_iw2_"+pol+"_slc.par " +date1+"_iw2_"+pol+"_slc.tops_par\n")
tab1.write(date1+"_iw3_"+pol+"_slc " +date1+"_iw3_"+pol+"_slc.par " +date1+"_iw3_"+pol+"_slc.tops_par\n")
tab1.close()

tab2 = open("slc2_tab.txt", "w")
tab2.write(date2+"_iw1_"+pol+"_slc " +date2+"_iw1_"+pol+"_slc.par " +date2+"_iw1_"+pol+"_slc.tops_par\n")
tab2.write(date2+"_iw2_"+pol+"_slc " +date2+"_iw2_"+pol+"_slc.par " +date2+"_iw2_"+pol+"_slc.tops_par\n")
tab2.write(date2+"_iw3_"+pol+"_slc " +date2+"_iw3_"+pol+"_slc.par " +date2+"_iw3_"+pol+"_slc.tops_par\n")
tab2.close()

## calculating the real range ##
# reading variables #
print("""### CALCULATING REAL RANGE AND VARIABLES FOR MULTILOOK ###\n
          ...""")
par = open(date1+"_iw1_"+pol+"_slc.par","r")
range_pixel_spacing = ""
azimuth_pixel_spacing = ""
incidence_angle = ""

for line in par:
    if "range_pixel_spacing" in line:
        linesplit = line.split()
        range_pixel_spacing = linesplit[len(linesplit)-2]
    if "azimuth_pixel_spacing" in line:
        linesplit = line.split()
        azimuth_pixel_spacing = linesplit[len(linesplit)-2]
    if "incidence_angle" in line:
        linesplit = line.split()
        incidence_angle = linesplit[len(linesplit)-2]
        break           
par.close()
range_pixel_spacing = float(range_pixel_spacing)
azimuth_pixel_spacing = float(azimuth_pixel_spacing)
incidence_angle = float(incidence_angle)

# Calculation #
incidence_angle = math.radians(incidence_angle) #in radians
range_pixel_spacing = range_pixel_spacing / math.sin(incidence_angle) #real pixel_spacing

ml_rg = dem_res/range_pixel_spacing
ml_az = dem_res/azimuth_pixel_spacing

## Mosaicing & Multilook ##
print("""### MOSAICING AND MULTILOOKING ###\n
          ...""")
isp.SLC_mosaic_S1_TOPS("slc1_tab.txt", date1+".slc", date1+".slc.par", ml_rg, ml_az,logpath=None, outdir=None,
                       shellscript=None)
isp.multi_look(date1+".slc", date1+".slc.par", date1+".mli", date1+".mli.par", ml_rg, ml_az)

## reading lines und columns ##
print("""### ACQUIRING MLI WIDTH AND LINES ###\n
          ...""")
par = open(date1+".mli.par","r") #Datei oeffnen
range_samples = ""
azimuth_lines = ""

for line in par:
    if "range_samples" in line:
        linesplit = line.split()
        range_samples = linesplit[len(linesplit)-1]
    if "azimuth_lines" in line:
        linesplit = line.split()
        azimuth_lines = linesplit[len(linesplit)-1]
        break
par.close()
mli1_width=int(range_samples)
mli1_lines=int(azimuth_lines)
print(""" MLI lines: %s  /  MLI width %s""" % (mli1_lines,mli1_width)) 

## LUT creation ## (creates EQA DEM)
print("""### CREATING LUT AND EQA.DEM ###\n
          ...""") 
diff.gc_map(MLI_par=date1+".mli.par", OFF_par="-", DEM_par=dem_name+"_par", DEM=dem_name, DEM_seg_par="EQA.dem_par",
       DEM_seg="EQA.dem", lookup_table=date1+".lt",
       lat_ovr='3', lon_ovr='3', sim_sar='-', u='-', v='-', inc=date1+'.inc', psi='-', pix='-', ls_map=date1+'.ls_map',
       frame='-', ls_mode='-', r_ovr='-', logpath=None, outdir=None, shellscript=None) #frame = 8; ls_mode = 2 in default

## reading lines and columns from the EQA DEM ##
print("""### ACQUIRING DEM WIDTH AND LINES ###\n
          ...""") 
par = open("EQA.dem_par","r")
width = ""
nlines = ""
for line in par:
    if "width" in line:
        linesplit = line.split()
        width = linesplit[len(linesplit)-1]
    if "nlines" in line:
        linesplit = line.split()
        nlines = linesplit[len(linesplit)-1]
        break
par.close()
dem_width=int(width)
dem_lines=int(nlines)
print(""" DEM lines: %s  /  DEM width %s""" % (dem_lines,dem_width)) 

## pixel area ##
print("""### CALCULATING PIXEL AREA ###\n
          ...""") 
diff.pixel_area(MLI_par=date1+".mli.par", DEM_par="EQA.dem_par", DEM="EQA.dem", lookup_table=date1+".lt",
                ls_map=date1+".ls_map", inc_map=date1+".inc", pix_sigma0="pix_sigma0",
                pix_gamma0='-', nstep=n_step, area_fact='-', sigma0_ratio='pix', gamma0_ratio='-' #nstep = default = 10 passt für unsere Anwendung
                , logpath=None, outdir=None, shellscript=None)

## converting DEM into radar geometry ##
print("""### CONVERT DEM INTO RADAR GEOMTERY ###\n
          ...""") 
diff.geocode(lookup_table=date1+".lt", data_in="EQA.dem", width_in=dem_width, data_out=date1+".hgt",
        width_out=mli1_width, nlines_out=mli1_lines, interp_mode='2', dtype='0',
        lr_in='-', lr_out='-', n_ovr='-', rad_max='-', nintr='-', logpath=None, outdir=None, shellscript=None)

## Creating SLC3_tab ## (for Coregistration)
tab3 = open("slc3_tab.txt", "w")
tab3.write(date2+"_iw1_rslc " +date2+"_iw1_rslc.par " +date2+"_iw1_rslc.tops_par\n")
tab3.write(date2+"_iw2_rslc " +date2+"_iw2_rslc.par " +date2+"_iw2_rslc.tops_par\n")
tab3.write(date2+"_iw3_rslc " +date2+"_iw3_rslc.par " +date2+"_iw3_rslc.tops_par\n")
tab3.close()

## Coregistration ##
print("""### COREGISTRATION ###\n
          ...""") 
diff.S1_coreg_TOPS(SLC1_tab="slc1_tab.txt", SLC1_ID=date1, SLC2_tab="slc2_tab.txt", SLC2_ID=date2, RSLC2_tab="slc3_tab.txt",
                   hgt=date1+'.hgt', rlks=round(ml_rg), azlks=round(ml_az), poly1='-', poly2='-',
                   cc_thresh=ccthresh, fraction_thresh=frac_thresh, ph_stdev_thresh=phstdev_thresh, 
                   cleaning=1, flag1='-', RSLC3_tab='-', RSLC3_ID='-',
                   logpath=None, outdir=None, shellscript=None)

# to check for plausibility #
#disp.disras(ras=date1+"_"+date2+".diff.bmp", mag='-', win_sz='-', logpath=None, outdir=None, shellscript=None)

## Multilook of second Scene ##
print("""### MULTILOOKING OF SECOND SCENE ###\n
          ...""") 
isp.multi_look(date2+".rslc", date2+".rslc.par", date2+".rmli", date2+".rmli.par", ml_rg, ml_az)
# look at it #
# disp.dis2pwr(date1+".mli", date2+".rmli", mli1_width, mli1_width)

## Baseline calculation ##
print("""### CALCULATING BASELINE ###\n
          ...""") 
isp.base_orbit(SLC1_par=date1+".slc.par", SLC2_par=date2+".rslc.par", baseline=date1+"_"+date2+".base",
               logpath=None, outdir=None, shellscript=None)

## Multiplication with backscatter intensities (topographic normalisation) ##
print("""### MULTIPLICATION WITH BACKSCATTER INTENSITY ###\n
          ...""") 
lat.product(data_1=date1+".mli", data_2="pix", product=date1+".cmli", width=mli1_width, bx=5, by=5, wgt_flag=1,
            logpath=None, outdir=None, shellscript=None)    
lat.product(data_1=date2+".rmli", data_2="pix", product=date2+".crmli", width=mli1_width, bx=5, by=5, wgt_flag=1,
            logpath=None, outdir=None, shellscript=None)    
# look at it # 
# disp.dis2pwr(date1+".mli", date1+".cmli", mli1_width, mli1_width)

## calculating coherence ##
print("""### CALCULATING COHERENCE ###\n
          ...""") 
lat.cc_ad(interf=date1+"_"+date2+".diff", pwr1=date1+".cmli", pwr2="-", slope="-", texture="-",
          cc_ad=date1+"_"+date2+".cc_ad", width=mli1_width, box_min=boxmin, box_max=boxmax, wgt_flag=1,
          loff='-', nl='-', logpath=None, outdir=None, shellscript=None)
            # box_min und box_max are default 3 and 9

## geocode back ##
print("""### GEOCODING BACK TO REAL GEOM ###\n
          ...""") 
diff.geocode_back(date1+"_"+date2+".cc_ad", mli1_width, date1+".lt", date1+"_"+date2+".geo.cc_ad",
                  dem_width, nlines_out='-',
                  interp_mode='-', dtype='-', lr_in='-', lr_out='-', order='-', e_flag='-',
                  logpath=None, outdir=None, shellscript=None)

## Conversion to geotiff ##
print("""### CREATING GEOTIFF FILE ###\n
          ...""") 
disp.data2geotiff("EQA.dem_par", date1+"_"+date2+".geo.cc_ad", 2, date1+"_"+date2+"_geo_cc_ad.tif",
                  nodata='-', logpath=None, outdir=None, shellscript=None)


# look at it #
# disp.dis_linear(date1+"_"+date2+".cc_ad", mli1_width, 5000)

# timer stop #
t_stop = perf_counter()
auxmodules.timer(t_start,t_stop)
