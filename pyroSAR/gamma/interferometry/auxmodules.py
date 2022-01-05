# -*- coding: utf-8 -*-
# Auxmodules for the coherence processing tool
# created by Marcel Felix & Friederike Metz

import os
import math
import zipfile
import spatialist
from pyroSAR.gamma.api import diff, isp, disp, lat


class HaltException(Exception):
    pass


def Script_init():
    print("""
          ################### WELCOME ######################### \n\n
          Welcome to the automated processing script for Sentinel 1 SAR Coherence pairs.
          Please read the following insructions and make necessary adjustments to your data.\n\n
          ################### CAPABILITIES #################### \n\n
          This script processes one coherence pair per execution.
          It needs to be called again for every coherence pair you want to process.\n\n
          ################### INSTRUCTIONS ####################\n
          ## folder structure ##
          Before you continue make sure the data structure in your project folder is as follows: \n
          -Main_project_folder
              --Daten
                  ---Your_area_of_interest_1 (optional)
                  ---..._2
              --other_folders \n
          all your data (DEM and SAR-scenes) should be stored in the "Daten" folder,
          OR in a specified sub-folder, in case you have multiple pairs to process.
          This is your working directory.\n
          
          ## DEM ##\n
          You need to provide a DEM. this script does not download any DEM automatically.
          This script can handle the >ALOS PALSAR DSM< and the >SRTM DEM<.
          
          # SRTM #
          If you have SRTM tiles, please put them directly in the working directory.\n
          # PALSAR #
          If you have PALSAR tiles, please unpack the packages in which they come in.
          Then leave the files in their respective folders within the working directory.
          The result should be like this: .../Daten/N060E125_N065E130/PALSAR_DSM.tif
          or this .../Daten/area_of_interest_1/N060E125_N065E130/PALSAR_DSM.tif\n
          
          ## SAR scenes ##
          Put the zip-files of your SAR scenes in their corresponding working directory.
          Please make sure only two scenes are in the folder!\n\n
          ################### START ########################### \n\n
          Checked the above instructions?
          Then you can begin to process. 
          Remember you can cancel the execution at any given time with CTRL+D
          """)
    if input("please enter 'start' (without quotation marks):") == "start":
        pass
    else:
        return Script_init()


def Path_init():
    print("""
          ################### PATHS ###########################\n\n
          Please provide the Path to your Main Project folder.
          Like this: /home/user/project_name/
          """)
    project_path = input("Project path:")
    if os.path.isdir(project_path):
        pass
    else:
        print("### PROVIDED PATH DOES NOT EXIST! ###")
        return Path_init()
    print("""
          \n
          Please provide the name of your working directory / data folder.
          Like this: Data/     or this: Data/NewYork/
          Please type slashes only as shown.
          """)
    data_path = input("working directory path:")
    if os.path.isdir(os.path.join(project_path, data_path)):
        data_path = os.path.join(project_path, data_path)
        return project_path, data_path
    else:
        print("### PROVIDED PATH DOES NOT EXIST! ###")
        return Path_init()


def Scene_ident(data_path):
    os.chdir(data_path)
    print("""
          ################### SCENE IDENTIFICATION ##########\n
          Scenes are being identified and unpacked...
          """)
    filelist = []
    for file in os.scandir(data_path):
        if file.name.startswith("S1") and file.name.endswith(".zip"):
            filelist.append(file.name)
    if len(filelist) > 2:
        print("""
              ### EXECUTION ERROR ###
              MORE THAN TWO SCENES (.zip) IN THE WORKING DIRECTORY!
              PLEASE CHECK THE PROBLEM AND RESTART THE SCRIPT.
              """)
        raise HaltException("Script execution stopped")
    if len(filelist) < 2:
        print("""
              ### EXECUTION ERROR ###
              LESS THAN TWO SCENES (.zip) IN THE WORKING DIRECTORY!
              PLEASE CHECK THE PROBLEM AND RESTART THE SCRIPT.
              """)
        raise HaltException("Script execution stopped")
    else:
        filelist = sorted(filelist)
        print("""
              Following scenes have been found: 
                  %s
                  %s
              """ % (filelist[0], filelist[1]))
        return filelist[0], filelist[1]


def Scene_unpack(data_path, scene1, scene2):
    with zipfile.ZipFile(os.path.join(data_path, scene1), 'r') as zip_ref:
        zip_ref.extractall(data_path)
    with zipfile.ZipFile(os.path.join(data_path, scene2), 'r') as zip_ref:
        zip_ref.extractall(data_path)


def Scene_ident2(data_path):
    os.chdir(data_path)
    filelist = []
    for file in os.scandir(data_path):
        if file.name.startswith("S1") and file.name.endswith(".SAFE"):
            filelist.append(file.name)
    if len(filelist) > 2:
        print("""
              ### EXECUTION ERROR ###
              MORE THAN TWO SCENES (.SAFE) IN THE WORKING DIRECTORY!
              PLEASE CHECK THE PROBLEM AND RESTART THE SCRIPT.
              """)
        raise HaltException("Script execution stopped")
    if len(filelist) < 2:
        print("""
              ### EXECUTION ERROR ###
              LESS THAN TWO SCENES (.SAFE) IN THE WORKING DIRECTORY!
              PLEASE CHECK THE PROBLEM AND RESTART THE SCRIPT.
              """)
        raise HaltException("Script execution stopped")
    else:
        filelist = sorted(filelist)
        print("""
              Following scenes have been found: 
                  %s
                  %s
              """ % (filelist[0], filelist[1]))
        return filelist[0], filelist[1]


def DEM_Select():
    print("""
          ################### DEM PREPARATION ###############\n
          Please provide the type of DEM you are using.
          """)
    dem_type = input("Enter DEM type (PALSAR or SRTM):")
    if dem_type.lower() == "srtm":
        print("The DEMs are being identified and merged...")
        return "srtm"
    if dem_type.lower() == "palsar":
        print("The DEMs are being identified and merged...")
        return "palsar"
    else:
        print("### PROVIDED TYPE DOES NOT EXIST! ###")
        return DEM_Select()


def PALSARDEM(datapath):
    # this function requires for the PALSAR-tif Rasters to be in their unpacke folders and
    # thos folders need to be in the working directory
    os.chdir(datapath)
    folderlist = []
    for folder in os.scandir():
        if folder.is_dir():
            if ("SAFE" not in folder.name) and ("zip" not in folder.name) and (
                    (folder.name[4] == "E") or (folder.name[4] == "W")):
                folderlist.append(folder.name)
    demlist = []
    for folder in folderlist:
        for file in os.scandir("./" + folder):
            if file.name.endswith("DSM.tif"):
                demlist.append(file.path)
    spatialist.auxil.gdalbuildvrt(demlist, "palsar.vrt", options=None, void=True)
    diff.dem_import("palsar.vrt", "palsar.dem", "palsar.dem_par")
    return "palsar.dem"


def SRTMDEM(datapath):
    # the DEM files should be directly in the working directory
    os.chdir(datapath)
    demlist = []
    for file in os.scandir():
        if file.name.endswith(".tif") and file.name.startswith("srtm"):
            demlist.append(file.name)
    spatialist.auxil.gdalbuildvrt(demlist, "srtm.vrt", options=None, void=True)
    diff.dem_import("srtm.vrt", "srtm.dem", "srtm.dem_par")
    return "srtm.dem"


def DEM_res():
    res = input(
        "\nPlease enter your desired resolution for Mulitlooking which should be in concordance with your DEM resolution:")
    if res.isdigit():
        return float(res)
    else:
        print("### PROVIDED INPUT IS NOT A PROPER NUMBER! ###")
        return DEM_res()


def Processing():
    print("""
         ################### PROCESSING PARAMETERS ########\n
         Now you have the choice to either go with default parameters for processing or specifiy many of them.
         For even more in depth options, please enter the scripts manually
         """)
    # gc_map_ls_mode = lsmode()
    # if gc_map_ls_mode == "2":
    #     pixel_area_allow = "yes"
    # else: pixel_area_allow = "no"
    defaultyesno = input("Do you want to continue and process with ALL DEFAULT parameters? yes/no:")
    if defaultyesno.lower() not in ["yes", "no"]:
        print("### PLEASE ONLY WRITE YES OR NO ###")
        return Processing()
    
    if defaultyesno.lower() == "no":
        pixel_area_nstep = nstep()
        coreg_ccthresh = cc_threshold()
        coreg_frac_thresh = fraction_thresh()
        coreg_ph_stdev_thresh = ph_stdev_thresh()
        cc_ad_boxmin, cc_ad_boxmax = cc_ad_box()
        return pixel_area_nstep, coreg_ccthresh, coreg_frac_thresh, coreg_ph_stdev_thresh, cc_ad_boxmin, cc_ad_boxmax
    else:
        return "-" * 6


def cc_ad_box():
    print("""
          ### Coherence Filter Boxes ###\n
          boxmin smallest correlation average box size default: 3.0
          boxmax largest correlation average box size  default: 9.0
          
          Please only consider values between 1 and 11. Else change the script manually.
          """)
    boxmin = input("\nPlease enter the desired boxmin for cc_ad:")
    if boxmin not in ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"]:
        print("### ! PROVIDED INPUT IS NOT BETWEEN 1 AND 10 ! ###")
        return cc_ad_box()
    
    boxmax = input("\nPlease enter the desired boxmin for cc_ad:")
    if boxmax not in ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"]:
        print("### ! PROVIDED INPUT IS NOT BETWEEN 1 AND 10 ! ###")
        return cc_ad_box()
    
    if int(boxmax) < int(boxmin):
        print("### ! PROVIDED INPUT IS NOT BIGGER THAN BOXMIN ! ###")
        return cc_ad_box()
    else:
        return boxmin, boxmax


def ph_stdev_thresh():
    print("""
          ### phase_stdev_threshold ###\n
          phase standard deviation threshold
          default is 0.8 radian
          """)
    stdev_thresh = input("\nPlease enter the desired phase-standard-deviation threshold for S1_coreg_TOPS:")
    if stdev_thresh.replace('.', '', 1).isdigit():
        if float(stdev_thresh) <= 2 * math.pi and float(stdev_thresh) >= 0:  # Phase from zero to 2xpi
            return stdev_thresh
        else:
            print("### ! PROVIDED INPUT IS NOT A PROPER VALUE ! ###")
            return ph_stdev_thresh()
    else:
        print("### ! PROVIDED INPUT IS NOT A PROPER VALUE ! ###")
        return ph_stdev_thresh()


def fraction_thresh():
    print("""
          ### fraction_threshold ###\n
          minimum valid fraction of unwrapped phase values
          default is 0.01
          """)
    frac_thresh = input("\nPlease enter the desired fraction threshold for S1_coreg_TOPS:")
    if frac_thresh.replace('.', '', 1).isdigit():
        if float(frac_thresh) <= 1 and float(frac_thresh) >= 0:
            return frac_thresh
        else:
            print("### ! PROVIDED INPUT IS NOT A PROPER VALUE ! ###")
            return fraction_thresh()
    else:
        print("### ! PROVIDED INPUT IS NOT A PROPER VALUE ! ###")
        return fraction_thresh()


def cc_threshold():
    print("""
          ### cc_threshold ###\n
          0-1 
          default is 0.8
          """)
    cc_thresh = input("\nPlease enter the desired coherence threshold for S1_coreg_TOPS:")
    if cc_thresh.replace('.', '', 1).isdigit():
        if float(cc_thresh) <= 1 and float(cc_thresh) >= 0:
            return cc_thresh
        else:
            print("### ! PROVIDED INPUT IS NOT A PROPER VALUE ! ###")
            return cc_threshold()
    else:
        print("### ! PROVIDED INPUT IS NOT A PROPER VALUE ! ###")
        return cc_threshold()


def nstep():
    print("""
          ### nstep ###\n
          number of steps to divide each dimension of the map pixels (default: 10)
          If the lookup table is on a fine grid relative to the slant-range resolution, 
          the value of nstep can be reduced to 8 or 6. In the case were the map pixels are larger, 
          then nstep should be increased.
          """)
    n_step = input("\nPlease enter the desired nstep for pixel_area:")
    check = []
    for i in range(21):  # Values larger 20 are not considered
        check.append(str(i))
    if n_step not in check:
        print("### ! PROVIDED INPUT IS NOT A PROPER VALUE ! ###")
        return nstep()
    else:
        return n_step


def lsmode():
    print("""
          ### ls_mode ###\n
          
          0: set to (0.,0.)
          1: linear interpolation across these regions
          2: actual value (default)
          3: nn-thinned

          """)
    ls_mode = input("\nPlease enter the desired ls_mode for gc_map:")
    if ls_mode not in ["0", "1", "2", "3"]:
        print("### ! PROVIDED INPUT IS NOT A PROPER VALUE ! ###")
        return lsmode()
    else:
        return ls_mode


def LastHints():
    print("""
          ################### PROCESSING INFORMATION ########\n
          From now on the script will proceed it's execution without user input.
          It uses default values. If you want you can adjust the Values 'in coherence.py'.
          A timer will run and inform you about the absolut execution time needed.
          Please be aware of errors that might occur during the process and check the log from time to time.
          """)


def timer(t_start, t_stop):
    t_sec = int(t_stop) - int(t_start)
    t_min = t_sec // 60
    t_rsec = t_sec - t_min * 60
    print("""
          ### THE SCRIPT FINISHED SUCCESSFULLY ###\n
          Execution time: %s minutes %s seconds
              """ % (t_min, t_rsec))
