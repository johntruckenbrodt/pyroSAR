##############################################################
# delete/remove files following defined patterns in the current directory and its subdirectories
# module of software gammaGUI
# John Truckenbrodt 2015
##############################################################

import sys

import os
import shutil

from ancillary import finder


# find all files matching the defined pattern(s)
items = finder(os.getcwd(), sys.argv[1].split(", "), regex=True if sys.argv[3] == "True" else False)

# exclude files in the export directory
items = [x for x in items if "/EXP/" not in x]

if len(items) > 0:
    path_exp = os.path.join(os.getcwd(), "EXP/")
    if sys.argv[2] == "export" and not os.path.exists(path_exp):
        os.makedirs(path_exp)
    print "the following files will be", {"export": "exported", "delete": "deleted"}[sys.argv[2]]+":"
    for item in items:
        print item
    decision = raw_input("proceed (y/n)?: ")
    if decision == "y":
        if sys.argv[2] == "export":
            for item in items:
                shutil.copy(item, path_exp)
        if sys.argv[2] == "delete":
            for item in items:
                os.remove(item)
    else:
        print "command aborted"
else:
    print "no files matching the defined pattern(s) found"



