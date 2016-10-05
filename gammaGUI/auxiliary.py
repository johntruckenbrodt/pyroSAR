
import os
import shutil
from ancillary import finder, writer
import subprocess as sp
from Tkinter import *
from gamma.util import ISPPar, gamma

import tkMessageBox


# wrapper for batch-processing
def batch(args, tags_search, tag_exclude):
    # find processing candidates
    candidates = finder(Environment.workdir.get(), tags_search)
    for item in candidates:
        items = list(args)
        items.append(item)
        # check whether processing results already exist and execute command if they do not
        if not os.path.isfile(item+tag_exclude):
            sp.check_call(items, cwd=Environment.workdir.get())


# gui design and structure dictionaries
class Environment(object):

    # place holder for the later defined working directory
    workdir = ""

    # general GUI design specifications
    bcolor = {"bg": 'black'}
    label_ops = {"bg": "black", "fg": "white", "font": ("system", 10, "bold")}
    header_ops = {"bg": "black", "fg": "white", "font": ("system", 12, "bold"), "pady": 5}
    menu_ops = {"bg": "black", "fg": "white", "font": ("system", 10, "bold"), "activebackground": "white", "activeforeground": "black", "bd": 2}
    button_ops = {"bg": "white", "fg": "black", "font": ("system", 10, "bold"), "activebackground": "black",
                  "activeforeground": "white", "highlightbackground": "white"}
    checkbutton_ops = {"bg": "black", "fg": "black", "font": ("system", 10, "bold"), "activebackground": "black", "width": 21, "height": 1, "bd": .5,
                  "activeforeground": "black", "highlightbackground": "black"}

    # define extra text printed in the dialog window
    args_extra = {"multilooking": "Compute MLIs from all SLCs in all subdirectories",
                  "cc_ad": "estimate coherence from all generated interferograms",
                  "cc_wave": "estimate coherence from all generated interferograms",
                  "interferogram flattening": "flatten all existing interferograms",
                  "srtm download": "download SRTM tiles covering all existing SLCs",
                  "geocoding": "geocode and topo-normalize all MLI, coherence and decomposition images",
                  "SLC calibration": "calibrate all existing SLCs",
                  "SLC coregistration (batched)": "batched coregistration with a tab-separated two-column file list",
                  "export/delete files": "remove files matching a pattern or copy them to a subdirectory EXP",
                  "temporal filtering": "temporal filtering of all intensity images and polarimetric decompositions"}

    # path (relative to working directory) and name of parameter files for specified commands
    parnames = {"baseline estimation": "/PAR/baseline.par",
                "cc_ad": "/PAR/coherence_ad.par",
                "cc_wave": "/PAR/coherence_wave.par",
                "dem preparation": "/PAR/dem.par",
                "geocoding": "/PAR/scene.par",
                "interferogram flattening": "/PAR/ph_slope_base.par",
                "interferogram generation": "/PAR/interferogram.par",
                "kml generation": "/PAR/kml.par",
                "SLC coregistration": "/PAR/coreg.par"}

    # dropdown button lookup and options
    dropoptions = {"str": ["method flag", "mode", "resampling method", "sensor"],
                   "int": ["wgt_ad", "wgt_wave", "wgt_tfilt"],
                   "interp_mode": ["1/dist", "NN", "SQR(1/dist)", "constant", "Gauss"],
                   "method flag": ["0", "1", "2", "3", "4"],
                   "mode": ["export", "delete"],
                   "resampling method": ["near", "bilinear", "cubic", "cubicspline", "lanczos", "average", "mode", "max", "min", "med", "q1", "q3"],
                   "sensor": ["All", "PSR1", "PSR2", "CSK", "S1", "TDX", "TSX", "RS2"],
                   "wgt_ad": ["constant", "gaussian"],
                   "wgt_wave": ["constant", "triangular", "gaussian", "none (phase only)"],
                   "wgt_tfilt": ["uniform", "linear", "gaussian"],
                   "arcsec": ["1", "3"],
                   "algorithm": ["Intensity CC", "Fringe Vis."],
                   "polynomial": [1, 3, 4, 6],
                   "oversampling": [1, 2, 4, 8]}

    # check button lookup table
    checkbuttons = ["topographic normalization", "regexp", "deramp (S1 SLC only)", "mosaic (S1 SLC only)", "orbit correction", "differential", "utm", "replace"]


# execute child window command
def execute(action, fileList, arguments):

    text = []
    if len(action) > 1:
        if ".py" in action[1]:
            # add path to name of script
            action[1] = finder(os.environ["PYTHONLAND"], [action[1]])[0]
    # list all necessary arguments
    for arg in action:
        text.append(arg)

    for obj in fileList:
        text.append(obj.file.get())

    if len(action) > 1:
        if "srtm.py" in action[1]:
            arguments.append(["SRTM_archive", fileList[0].file])

    # evaluate processing parameters
    for i in range(0, len(arguments)):
        if arguments[i][0] in Environment.dropoptions["int"]:
            arguments[i][1] = str(Environment.dropoptions[arguments[i][0]].index(arguments[i][1].get()))
        else:
            try:
                arguments[i][1] = arguments[i][1].get()
            except:
                pass

        if arguments[i][1] != "SRTM_archive":
            text.append(arguments[i][1])

    # execute subprocess (GAMMA) commands
    try:
        if len(action) > 1:
            path_par = os.path.join(Environment.workdir.get(), "PAR/")
            if not os.path.exists(path_par):
                os.makedirs(path_par)
            basename_script = "coreg" if "coreg" in action[1] else os.path.splitext(os.path.basename(action[1]))[0]
            writer(path_par+basename_script+".par", arguments, strfill=False)
            if os.path.basename(action[1]) == "reader_old.py":
                print "#############################################"
                print "importing files..."
                del(text[0])
                importer(text)
                print "...done"
            elif os.path.basename(action[1]) == "multilook.py":
                print "#############################################"
                print "multilooking files..."
                batch(text, ["*_slc", "*_slc_cal", "*_slc_reg", "*_slc_cal_reg"], "_mli")
                print "...done"
            elif os.path.basename(action[1]) == "srtm.py":
                sp.check_call(text[:2], cwd=Environment.workdir.get())
            else:
                sp.check_call(text, cwd=Environment.workdir.get())

        elif "dis" in action[0]:
            if "dem_par" in action[0]:
                text.append(text[1]+".par")
            if "disrmg" in action[0]:
                if len(text[2]) == 0:
                    text[2] = "-"
                samples = ISPPar(finder(Environment.workdir.get(), [re.findall("[A-Z0-9_]{9,10}[0-9T]{15}_[HV]{2}_slc(?:_cal|)?", text[1])[0]])[0]+".par").range_samples
                text.append(samples)

            elif "mph" in action[0]:
                samples = ISPPar(finder(Environment.workdir.get(), [re.findall("[A-Z0-9_]{9,10}[0-9T]{15}_[HV]{2}_slc(?:_cal|)?", text[1])[0]])[0]+".par").range_samples
                text.append(samples)
                if "2" in action[0]:
                    samples = ISPPar(finder(Environment.workdir.get(), [re.findall("[A-Z0-9_]{9,10}[0-9T]{15}_[HV]{2}_slc(?:_cal|)?", text[2])[0]])[0]+".par").range_samples
                    text.append(samples)
            else:
                for i in range(1, len(text)):
                    try:
                        text.append(ISPPar(text[i] + ".par").range_samples)
                    except:
                        continue
            print "#############################################"
            sp.check_call([str(x) for x in text])
        if "dis" not in action[0]:
            tkMessageBox.showinfo("execution end", "processing successfully completed")
    except IOError:
        tkMessageBox.showerror("execution error", "incorrect or missing parameters")


# search working directory for scene folders and store image names into list of objects
def grouping(filepath=os.getcwd()):
    # list top level subfolders (i.e. no sub-subfolders) of the defined filepath if their names match the defined expression (i.e. sensor_timestamp)
    scenes = sorted([os.path.join(filepath, z) for z in os.listdir(filepath) if re.search("[0-9]{8}T[0-9]{6}", z)])
    scenes = finder(filepath, ["[0-9]{8}T[0-9]{6}"], regex=True, foldermode=2, recursive=False)
    # return list of scene objects
    return [Tuple(scene) for scene in scenes]


# wrapper for file import
# imported files are first moved to a temporary directory
# folders are created with names based on information retrieved from the newly generated gamma parameter files
# files in the temporary directory are then moved to the newly created folder and the temporary folder removed
def importer(text):
    # set all files and folders of the input directory as possible import candidates
    candidates = finder(text[1], ["*"], foldermode=1)
    # create directories
    path_tmp = os.path.join(Environment.workdir.get(), "Temp")
    path_log = os.path.join(Environment.workdir.get(), "LOG/IMP")
    for path in [path_tmp, path_log]:
        if not os.path.exists(path):
            os.makedirs(path)
    # lists for storing coordinates of imported files
    latitudes = []
    longitudes = []

    for item in candidates:
        # reduce list to scriptname
        args = [text[0]]

        # if set, add sensor
        if text[2] != "All":
            args.append("-m")
            args.append(text[2])

        # add filename to list of execution arguments
        args.append("-s")
        args.append(item)

        # try to import file, if not possible (i.e. not an accepted format) continue with the next file
        try:
            gamma(args, path_tmp, path_log)
            # find next best parameter file of the newly imported files
            name_par = finder(path_tmp, ["*.par"])[0]

            print re.sub("w[0-9]__", "wx__", os.path.basename(name_par)[:25])

            # read parameter file
            par = ISPPar(name_par)

            # extract time stamp and sensor from parameter file name
            time = re.findall("[0-9]{8}T[0-9]{6}", name_par)[0]
            sensor = os.path.basename(name_par).split("_")[0]

            # define path for the imported files to be stored
            path_out = os.path.join(Environment.workdir.get(), sensor+"_"+time)
            if not os.path.exists(path_out):
                os.makedirs(path_out)
            latitudes.append(par.center_latitude)
            longitudes.append(par.center_longitude)
            for x in finder(path_tmp, "*"):
                try:
                    shutil.move(x, path_out)
                except:
                    os.remove(x)
        except:
            continue
    os.rmdir(path_tmp)
    # give warning if lat/long differs by more than 1 degree
    # if max(longitudes) - min(longitudes) > 1 or max(latitudes) - min(latitudes) > 1:
    #     tkMessageBox.showwarning("working directory warning", "latitudes or longitudes differ by more than 1deg!")


# create an entry form for execution parameters
def makeform(self, fields, defaultentries):
    entries = []
    for i in range(0, len(fields)):
        row = Frame(self, Environment.bcolor)
        row.pack(side=TOP, padx=5, pady=5)
        lab = Label(row, Environment.label_ops, width=25, text=fields[i], anchor="w")
        lab.pack(side=LEFT)

        # add dropdown button
        if fields[i] in Environment.dropoptions:
            optionlist = Environment.dropoptions[fields[i]]
            var = StringVar(self)
            if fields[i] in Environment.dropoptions["int"]:
                default = Environment.dropoptions[fields[i]][int(defaultentries[i])]
            else:
                default = defaultentries[i]
            var.set(default)
            ent = OptionMenu(row, var, *optionlist)
            ent.pack(side=RIGHT, expand=YES)
            ent.config(Environment.button_ops, width=14, bd=.5)
            entries.append([fields[i], var])
        # add checkbutton
        elif fields[i] in Environment.checkbuttons:
            var = StringVar(self)
            var.set(str(defaultentries[i]))
            ent = Checkbutton(row, variable=var, onvalue="True", offvalue="False")
            ent.pack(side=RIGHT, expand=YES)
            ent.config(Environment.checkbutton_ops)
            entries.append([fields[i], var])
        # add text entry field
        else:
            ent = Entry(row)
            ent.pack(side=RIGHT, expand=YES)
            ent.insert(0, defaultentries[i])
            entries.append([fields[i], ent])
    return entries


# create index objects containing paths to all products of a particular scene
# a list of entries can be queried by calling Tuple.__dict__.keys(); Tuple.index contains only those entries that represent individual images or image groups
class Tuple(object):
    def __init__(self, scene):
        self.index = []
        self.main = scene
        self.basename = os.path.basename(scene)
        pattern = "[0-9]{8}T[0-9]{6}"
        for x in [y for y in finder(scene, "*") if ".par" not in y and ".hdr" not in y]:
            base = os.path.basename(x)
            if re.search(pattern, base):
                tag = base[re.search(pattern, base).end()+1:]
                setattr(self, tag, x)
                self.index.append(tag)
        # append interferometric products
        isp_full = [y for y in finder(os.path.join(os.path.split(scene)[0], "ISP"), "*") if re.search(pattern, y) and not re.search("(?:bmp|hdr|par|kml)$", y)]
        isp_sub = [z for z in isp_full if re.findall(pattern, os.path.basename(z))[0] in self.basename]
        isp_tags = list(set([x[re.search(re.findall(pattern+"_[HV]{2}_", os.path.basename(x))[1], x).end():] for x in isp_sub]))
        for tag in isp_tags:
            setattr(self, tag, [x for x in isp_sub if re.search(tag+"$", x)])
            self.index.append(tag)
        self.index.sort()

    def getTop(self, pattern):
        candidates = [x for x in self.index if re.search(pattern, x)]
        if len(candidates) > 0:
            return getattr(self, max(candidates, key=len))
        else:
            raise AttributeError("no entries matching the defined pattern")
