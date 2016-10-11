##############################################################
# multi-temporal filtering of SAR intensity images
# module of software gammaGUI
# John Truckenbrodt 2015
##############################################################

"""
The following tasks are performed by executing this script:
-for each polarization separately:
--search for processing candidates within the working directory
--crop these images to their common extent (using R package raster)
--swap bytes of files processed in R back to big endian for use in GAMMA
--write textfile table containing in one column the processing candidate and in the other the name of the filtered image
--perform multitemporal filtering (temp-filt) using the newly created textfile table
"""

import sys

import re
import os
import subprocess as sp

from ancillary import ReadPar, dissolve
from envi import HDRobject, hdr
from auxiliary import grouping
import gamma

# read processing parameter textfile
par = ReadPar(os.path.join(os.getcwd(), "PAR/tempfilter.par"), type="exe")

# use topographically normalized images?
tnorm = True if par.topographic_normalization == "True" else False

# define (and create) directory for logfile
path_log = os.path.join(os.getcwd(), "LOG/LAT/")
if not os.path.exists(path_log):
    os.makedirs(path_log)

# create list of scene tuple objects
tuples = grouping()

counter = 0
for tag in ["HH", "VV", "HV", "VH", "pauli_alpha", "pauli_beta", "pauli_gamma", "ctd_1_mag", "ctd_2_mag", "ctd_3_mag", "fdd_pd", "fdd_ps", "fdd_pv"]:

    # collect processing candidates
    processlist = dissolve([[getattr(scene, x) for x in scene.index if re.search(tag, x)] for scene in tuples])

    if tnorm:
        processlist = [x for x in processlist if re.search("_norm_geo$", x) and not os.path.isfile(x+"_tfilt")]
    else:
        processlist = [x for x in processlist if re.search("_geo$", x) and "norm" not in x and not os.path.isfile(x+"_tfilt")]

    if len(processlist) > 1:
        counter += 1
        if counter == 1:
            print "#############################################"
            print "multitemporal filtering started..."
        print tag
        # crop images to their common extent
        # for convenience this is outsourced to R (package raster required)
        sp.check_call(dissolve(["Rscript", os.path.join(os.path.dirname(sys.argv[0]), "multicrop.R"), processlist]))

        # update processing list with new names
        processlist = [x+"_sub1" for x in processlist]

        # swap byte order of R-processed files to big endian, create hdr file for the final product and remove redundant intermediate files
        for x in processlist:
            sp.check_call(["swap_bytes", x, x[:-1], "4"], stdout=sp.PIPE)
            hdrobj = HDRobject(x+".hdr")
            hdrobj.description = hdrobj.description[:-4]+"tfilt"
            hdrobj.byte_order = "1"
            hdr(hdrobj, x[:-4]+"tfilt.hdr")
            os.remove(x)
            os.remove(x+".hdr")
            os.remove(x+".aux.xml")

        # update processing list with new names
        processlist = [x[:-1] for x in processlist]

        # create table textfile for temporal filter
        name_ftab = os.path.join(os.getcwd(), "PAR/tab_tempfilt_"+tag)
        with open(name_ftab, "w") as out:
            for image in processlist:
                out.write(image+"\t"+image[:-3]+"tfilt\n")

        # read header file
        meta = HDRobject(processlist[0][:-3]+"tfilt.hdr")

        # perform temporal filtering
        gamma.process(["temp_filt", name_ftab, meta.samples, par.waz, par.wr, par.wgt_tfilt], os.getcwd(), path_log)

        for x in processlist:
            os.remove(x)
if counter == 0:
    print "#############################################"
    print "no candidates for filtering found"
    print "#############################################"
else:
    print "...done"
    print "#############################################"
