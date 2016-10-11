##############################################################
# Kennaugh matrix derivation from scattering matrix elements
# module of software gammaGUI
# John Truckenbrodt 2015-16
##############################################################
import os
import re

from ancillary import dissolve, finder
from gammaGUI.auxiliary import grouping
import gamma
from envi import hdr

path_log = os.path.join(os.getcwd(), "LOG/LAT/")
if not os.path.exists(path_log):
    os.makedirs(path_log)

# create list of scene tuple objects
tuples = grouping()

print "#############################################"
print "creating kennaugh decomposition..."

for scene in tuples:
    print scene.basename

    try:
        hh_mli = scene.getTop("HH_(?:slc_|)(?:cal_|)mli$")
        nlines = gamma.ISPPar(hh_mli+".par").azimuth_lines
    except AttributeError:
        print "...skipped; no appropriate files found"
        continue

    elements = {1: ["t11", "t22", "t33"],
                2: ["t12"],
                3: ["t13"],
                4: ["t23"],
                5: ["t11", "t22", "t33"],
                6: ["t23"],
                7: ["t13"],
                8: ["t11", "t22", "t33"],
                9: ["t12"],
                10: ["t11", "t22", "t33"]}

    values = []
    keys = ["t11", "t22", "t33", "t12", "t13", "t23"]
    for key in keys:
        try:
            item = scene.getTop(key+"$")
            values.append(item)
        except AttributeError:
            values.append("-")
    valids = [x for x in values if x != "-"]
    if len(valids) == 0:
        print "...skipped; no appropriate files found"
        continue
    else:
        components = dict(zip(keys, values))
        base = re.sub("_t[123]{2}$", "", valids[0])

        for element in elements:
            requirements = [components[x] for x in elements[element]]
            if "-" not in requirements:
                gamma.process(dissolve(["KENNAUGH_MATRIX", values, nlines, base, element]), logpath=path_log)

    for element in finder(os.path.dirname(base), ["\.[tk][1-4]{2}$"], regex=True):
        new = re.sub("\.k", "_k", element) if re.search("\.k[1-4]{2}", element) else re.sub("\.t", "_k", element)
        os.rename(element, new)
        hdr(hh_mli+".par", new+".hdr")

print "...done"
