##############################################################
# batched SLC coregistration
# module of software gammaGUI
# John Truckenbrodt 2015
##############################################################

"""
The following tasks are performed by executing this script:
-read a provided two-column list (first and second column for primary and secondary SLC respectively)
--within this file '#' can be used for skipping coregistration couples
-print warning and continue with the next couple in case no file is found in the working directory or if the file descriptor fits multiple files
-call script coreg.py for each image couple to perform coregistration
"""

import sys

import os
import subprocess as sp

from ancillary import finder

with open(sys.argv[1], "r") as inlist:
    processlist = [line.split() for line in inlist if not line.startswith("#") and not line.strip() == ""]

for couple in processlist:
    slc1 = finder(os.getcwd(), [couple[0]+"$"], regex=True)
    slc2 = finder(os.getcwd(), [couple[1]+"$"], regex=True)

    if len(slc1) != 1:
        print "descriptor", slc1, "ambiguous or file not existing"
        continue
    elif len(slc2) != 1:
        print "descriptor", slc2, "ambiguous or file not existing"
        continue
    else:
        sp.check_call(["python", os.path.join(os.path.dirname(sys.argv[0]), "coreg.py"), slc1[0], slc2[0]], cwd=os.getcwd())
