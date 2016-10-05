##############################################################
# SLC coregistration
# module of software gammaGUI
# John Truckenbrodt 2015
##############################################################

"""
The following tasks are performed by executing this script:
-reading of a parameter file coreg.par (in subfolder PAR of GUI working directory)
--see object par for necessary values; file is automatically created by starting the script via the GUI
-if necessary, creation of output and logfile directories
-creating the offset parameter file
-estimation of initial range and azimuth offsets
-precise estimation of offset polynomials
-generation of offsets polynomials
-resampling of slc file
"""

import os
import sys
from ancillary import ReadPar
from gamma.util import correlate, init_offset, gamma
from auxiliary import Environment

# retrieve additional arguments
slc1 = sys.argv[1]
slc2 = sys.argv[2]

# read processing parameter textfile
par = ReadPar(os.path.join(os.getcwd(), "PAR/coreg.par"), type="exe")
par.algorithm = Environment.dropoptions["algorithm"].index(par.algorithm)+1

# set SNR theshold (this should not be changed)
thres = 7.0

# define (and create) directories for processing results and logfile
path_log = os.path.join(os.getcwd(), "LOG/ISP/")
path_out = os.path.join(os.getcwd(), "ISP/")
for path in [path_log, path_out]:
    if not os.path.exists(path):
        os.makedirs(path)

# concatenate output names
name_base = os.path.join(path_out, os.path.basename(slc1)+"_"+os.path.basename(slc2)+"_")
name_coffs = name_base+"coffs"
name_coffsets = name_base+"coffsets"
name_off = name_base+"off"
name_offs = name_base+"offs"
name_offsets = name_base+"offsets"
name_snr = name_base+"snr"
name_reg = name_base+"reg"

print "#############################################"
print os.path.basename(slc1), "->", os.path.basename(slc2)
print "----------"
print "coregistration started..."

gamma(["create_offset", slc1 + ".par", slc2 + ".par", name_off, par.algorithm, 1, 1, 0], path_out, path_log)

print "...initial offset estimation"
init_offset(slc1, slc2, name_off, path_log, thres=thres)

print "...cross-correlation"
correlate(slc1, slc2, name_off, name_offs, name_offsets, name_snr, name_coffs, name_coffsets, path_log, maxwin=2048, overlap=.3, poly=par.polynomial, ovs=par.oversampling, thres=thres)

print "...resampling of SLC file"
gamma(["SLC_interp", slc2, slc1 + ".par", slc2 + ".par", name_off, name_reg, name_reg + ".par"], path_out, path_log)

print "...done"
