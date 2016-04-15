#!/usr/bin/env python
##############################################################
# main interface for software gammaGUI
# John Truckenbrodt 2014-2015
##############################################################

import os
import tkMessageBox
from Tkinter import *
from PIL import ImageTk, Image

from dialog import Dialog
from fileQuery import FileQuery
from auxiliary import Environment


# CENTRAL CLASS FOR GUI MAIN WINDOW
class Main(Frame):
    def __init__(self, *args, **kwargs):
        Frame.__init__(self, *args, **kwargs)

        # set up main window
        self.master.geometry("600x400")
        root.title("GAMMA GUI")
        header = Label(root, Environment.header_ops, text="Welcome to the GEO410 GAMMA GUI")
        header.pack()

        # DEFINE EXECUTION ARGUMENTS

        # each list of arguments must always contain four sub-lists:
        # -(description) , (interpreter/compiler), name of script/command
        # -input file descriptors
        # -input parameter names
        # -input parameter default values
        # in case no files/parameters are required for the command, an empty sub-list must be provided

        # #######################################################################################################
        # IMP
        # raw data import
        args_import = [["general data import", sys.executable, "reader_old.py"], ["import directory"], ["sensor", "deramp (S1 SLC only)", "mosaic (S1 SLC only)"], ["All", True, True]]
        args_ers = [["ers data import", sys.executable, "reader_ers.py"], ["import directory"], ["orbit correction"], [True]]
        # #######################################################################################################
        # ISP
        # multilooking
        args_mli = [["multilooking", sys.executable, "multilook.py"], [], ["target resolution [m]"], ["automatic"]]

        # calibration
        args_cal_slc = [["SLC calibration", sys.executable, "cal_slc.py"], [], ["replace"], [True]]

        # SLC coregistration
        args_cor = [["SLC coregistration", sys.executable, "coreg.py"], ["primary SLC", "secondary SLC"], ["algorithm", "oversampling", "polynomial"], ["Intensity CC", 1, 4]]
        args_cor2 = [["SLC coregistration (batched)", sys.executable, "coreg_batch.py"], ["image list"], ["algorithm", "oversampling", "polynomial"], ["Intensity CC", 1, 4]]

        # interferogram generation
        args_int = [["interferogram generation", sys.executable, "interferogram.py"], [],
                    ["loff", "sps_flg", "azf_flg", "rp1_flg", "rp2_flg"],
                    [0, 1, 1, 1, 1]]

        # baseline estimation
        args_ibl = [["baseline estimation", sys.executable, "baseline.py"], [],
                    ["method flag", "nrfft", "nazfft", "r_samp", "az_line"],
                    [2, 1024, 1024, "-", "-"]]

        # interferogram flattening
        args_flt = [["interferogram flattening: bsl", sys.executable, "ph_slope_base.py"], [], [], []]
        args_flt2 = [["interferogram flattening: dem", sys.executable, "int_sim.py"], [], [], []]

        # coherence estimation
        args_cca = [["cc_ad", sys.executable, "coherence_ad.py"], [], ["box_min", "box_max", "wgt_ad", "differential"], [3, 7, 1, True]]
        args_ccw = [["cc_wave", sys.executable, "coherence_wave.py"], [], ["bx", "by", "wgt_wave", "differential"], [5, 5, 2, True]]
        # #######################################################################################################
        # GEO
        # download of SRTM tiles

        # DEM preparation and Geocoding
        args_srtm_prep = [["srtm preparation", sys.executable, "srtm.py"], ["SRTM archive"],
                    ["targetres", "utm", "arcsec"], [20, True, 1]]

        args_dem = [["dem preparation", sys.executable, "dem_preparation.py"], ["DEM"], [], []]

        args_geo = [["geocoding", sys.executable, "geocoding.py"], ["DEM"],
                    ["frame", "interp_mode", "oversampling", "polynomial", "interp_mode_intensity", "interp_mode_coherence", "topographic normalization"],
                    [8, "1/dist", 2, 4, 2, 2, True]]
        # #######################################################################################################
        # LAT

        args_pauli = [["pauli decomposition", sys.executable, "pauli.py"], [], [], []]
        args_coh = [["coherence matrix", sys.executable, "mat_coh.py"], [], [], []]
        args_cov = [["covariance matrix", sys.executable, "mat_cov.py"], [], [], []]
        args_huy = [["huynen decomposition", sys.executable, "huynen.py"], [], [], []]
        args_clo = [["cloude decomposition", sys.executable, "cloude.py"], [], [], []]
        args_fd3c = [["cloude decomposition", sys.executable, "fd3c.py"], [], [], []]
        args_circ = [["circular polarization basis transformation", sys.executable, "circular.py"], [],
                     ["constant_r", "constant_i", "factorHH_r", "factorHH_i", "factorVV_r", "factorVV_i", "factorHV_r", "factorHV_i", "pixav_x", "pixav_y"],
                     [0, 0, .5, 0, .5, 0, .5, 0, 1, 1]]
        args_kro = [["krogager decomposition", sys.executable, "krogager.py"], [], [], []]
        args_haa = [["haalpha decomposition", sys.executable, "haalpha.py"], [], [], []]
        args_ken = [["kennaugh decomposition", sys.executable, "kennaugh.py"], [], [], []]

        args_tfilt = [["temporal filtering", sys.executable, "tempfilter.py"], [], ["waz", "wr", "wgt_tfilt", "topographic normalization"], [3, 3, 0, True]]

        # #######################################################################################################
        # DISP
        args_dispwr = [["dispwr"], ["File"], [], []]
        args_dis2pwr = [["dis2pwr"], ["File 1", "File 2"], [], []]
        args_disSLC = [["disSLC"], ["File"], ["Width"], [2500]]
        args_dis2SLC = [["dis2SLC"], ["File 1", "File 2"], ["Width 1", "Width 2"], [2500, 2500]]
        args_discc = [["discc"], ["coherence", "intensity"], [], []]
        args_dismph = [["dismph"], ["interferogram"], [], []]
        args_dis2mph = [["dis2mph"], ["interferogram1", "interferogram2"], [], []]
        args_dishgt = [["disdem_par"], ["DEM File"], [], []]
        args_disrmg = [["disrmg"], ["Phase", "Intensity (optional)"], [], []]
        # #######################################################################################################
        # AUX
        args_cpdel = [["export/delete files", sys.executable, "cpdel.py"], [], ["pattern", "mode", "regexp"], ["", "export", False]]
        args_stack = [["layer stacking", sys.executable, "layerstack.py"], ["import directory", "output file", "shapefile"],
                      ["pattern", "regexp", "resampling method", "target resolution (x,y)", "NA value"],
                      ["", False, "bilinear", "", -999]]
        # #######################################################################################################

        # CREATE MENU BAR

        # #######################################################################################################

        # main menu
        menubar = Menu(self.master, Environment.menu_ops)
        self.master.config(menu=menubar)
        fileMenu = Menu(menubar, tearoff=0)
        # #######################################################################################################
        # MSP
        impMenu = Menu(fileMenu, Environment.menu_ops)
        menubar.add_cascade(label="IMP", menu=impMenu)
        impMenu.add_command(label="General Data Import", command=lambda: self.callChild(args_import))
        impMenu.add_command(label="ERS Data Import", command=lambda: self.callChild(args_ers))
        # #######################################################################################################
        # IMP
        ispMenu = Menu(fileMenu, Environment.menu_ops)
        corMenu = Menu(ispMenu, Environment.menu_ops)
        fltMenu = Menu(ispMenu, Environment.menu_ops)
        cohMenu = Menu(ispMenu, Environment.menu_ops)
        calMenu = Menu(ispMenu, Environment.menu_ops)
        menubar.add_cascade(label="ISP", menu=ispMenu)

        ispMenu.add_cascade(label="SLC coregistration", menu=corMenu)
        corMenu.add_command(label="single", command=lambda: self.callChild(args_cor))
        corMenu.add_command(label="batched", command=lambda: self.callChild(args_cor2))
        ispMenu.add_command(label="Multilooking", command=lambda: self.callChild(args_mli))
        ispMenu.add_command(label="Interferogram Generation", command=lambda: self.callChild(args_int))
        ispMenu.add_command(label="Baseline Estimation", command=lambda: self.callChild(args_ibl))
        ispMenu.add_cascade(label="Interferogram Flattening", menu=fltMenu)
        fltMenu.add_command(label="Estimation from Baseline", command=lambda: self.callChild(args_flt))
        fltMenu.add_command(label="DEM Interferogram Simulation", command=lambda: self.callChild(args_flt2))
        ispMenu.add_cascade(label="Coherence Estimation", menu=cohMenu)
        cohMenu.add_command(label="CC_AD", command=lambda: self.callChild(args_cca))
        cohMenu.add_command(label="CC_WAVE", command=lambda: self.callChild(args_ccw))
        ispMenu.add_cascade(label="Calibration", menu=calMenu)
        calMenu.add_command(label="SLC", command=lambda: self.callChild(args_cal_slc))
        # #######################################################################################################
        # GEO
        geoMenu = Menu(fileMenu, Environment.menu_ops)
        menubar.add_cascade(label="GEO", menu=geoMenu)
        geoMenu.add_command(label="SRTM preparation", command=lambda: self.callChild(args_srtm_prep))
        geoMenu.add_command(label="DEM preparation", command=lambda: self.callChild(args_dem))
        geoMenu.add_command(label="Geocoding", command=lambda: self.callChild(args_geo))
        # #######################################################################################################
        # LAT
        latMenu = Menu(fileMenu, Environment.menu_ops)
        polMenu = Menu(latMenu, Environment.menu_ops)
        filtMenu = Menu(latMenu, Environment.menu_ops)
        menubar.add_cascade(label="LAT", menu=latMenu)
        latMenu.add_cascade(label="Polarimetry", menu=polMenu)
        latMenu.add_cascade(label="Filtering", menu=filtMenu)
        polMenu.add_command(label="Pauli Decompostion", command=lambda: self.callChild(args_pauli))
        polMenu.add_command(label="Coherence Matrix", command=lambda: self.callChild(args_coh))
        polMenu.add_command(label="Covariance Matrix", command=lambda: self.callChild(args_cov))
        polMenu.add_command(label="Huynen Decompostion", command=lambda: self.callChild(args_huy))
        polMenu.add_command(label="Cloude Decompostion", command=lambda: self.callChild(args_clo))
        polMenu.add_command(label="Freeman-Durden 3 Component Decompostion", command=lambda: self.callChild(args_fd3c))
        polMenu.add_command(label="Circular Polarization Basis Transformation", command=lambda: self.callChild(args_circ))
        polMenu.add_command(label="Krogager Decompostion", command=lambda: self.callChild(args_kro))
        polMenu.add_command(label="H-A-Alpha Decompostion", command=lambda: self.callChild(args_haa))
        polMenu.add_command(label="Kennaugh Decompostion", command=lambda: self.callChild(args_ken))
        filtMenu.add_command(label="Multitemporal Filtering", command=lambda: self.callChild(args_tfilt))
        # #######################################################################################################
        # DISP
        dispMenu = Menu(fileMenu, Environment.menu_ops)
        pwrMenu = Menu(dispMenu, Environment.menu_ops)
        phaMenu = Menu(dispMenu, Environment.menu_ops)
        slcMenu = Menu(dispMenu, Environment.menu_ops)
        intMenu = Menu(dispMenu, Environment.menu_ops)
        ccMenu = Menu(dispMenu, Environment.menu_ops)
        hgtMenu = Menu(dispMenu, Environment.menu_ops)
        menubar.add_cascade(label="DISP", menu=dispMenu)

        dispMenu.add_cascade(label="Multilook", menu=pwrMenu)
        pwrMenu.add_command(label=args_dispwr[0], command=lambda: self.callChild(args_dispwr))
        pwrMenu.add_command(label=args_dis2pwr[0], command=lambda: self.callChild(args_dis2pwr))

        # SLC image display
        dispMenu.add_cascade(label="Single Look Complex", menu=slcMenu)
        slcMenu.add_command(label=args_disSLC[0], command=lambda: self.callChild(args_disSLC))
        slcMenu.add_command(label=args_dis2SLC[0], command=lambda: self.callChild(args_dis2SLC))

        # correlation image display
        dispMenu.add_cascade(label="Coherence", menu=ccMenu)
        ccMenu.add_command(label=args_discc[0], command=lambda: self.callChild(args_discc))

        # interferogram display
        dispMenu.add_cascade(label="Interferogram", menu=intMenu)
        intMenu.add_command(label=args_dismph[0], command=lambda: self.callChild(args_dismph))
        intMenu.add_command(label=args_dis2mph[0], command=lambda: self.callChild(args_dis2mph))

        # phase image display
        dispMenu.add_cascade(label="Phase", menu=phaMenu)
        phaMenu.add_command(label=args_disrmg[0], command=lambda: self.callChild(args_disrmg))

        # height image display
        dispMenu.add_cascade(label="Height Images", menu=hgtMenu)
        hgtMenu.add_command(label=args_dishgt[0], command=lambda: self.callChild(args_dishgt))

        # dialog for defining the working directory
        self.frame = Frame(self, bg="white", height=30, bd=2, width=300)
        self.frame.pack()
        Environment.workdir = FileQuery(self.frame, "Working Directory", 2).file
        # #######################################################################################################
        # AUX
        auxMenu = Menu(fileMenu, Environment.menu_ops)
        menubar.add_cascade(label="AUX", menu=auxMenu)
        auxMenu.add_command(label="Export/Delete", command=lambda: self.callChild(args_cpdel))
        auxMenu.add_command(label="Layer Stacking", command=lambda: self.callChild(args_stack))
        # #######################################################################################################

    # CREATE CHILD WINDOWS
    def callChild(self, arguments):
        if len(Environment.workdir.get()) == 0:
            tkMessageBox.showerror("Directory Missing", "please define working directory")
        else:
            Dialog(arguments)
# #######################################################################################################
# INITIALIZE GUI
if __name__ == "__main__":
    root = Tk()
    root.resizable(width=FALSE, height=FALSE)
    # load image as main window background
    basedir = os.path.dirname(os.path.abspath(__file__))

    image = Image.open(os.path.join(basedir,  "background.jpg"))
    tkimage = ImageTk.PhotoImage(image)
    Label(root, image=tkimage).place({"x": 0, "y": 0, "relwidth": 1, "relheight": 1})
    # call main window
    Main(root, Environment.bcolor).pack({"side": "top", "fill": "both"})
    root.mainloop()
