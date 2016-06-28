
import os
import re
from raster import stack
from ancillary import finder
from time import mktime, strptime

site = "Kilombero"
maindir = "/geonfs01_vol1/ve39vem/swos_process"

# store results in separate files
sep = False

srcfiles = finder(os.path.join(maindir, site, "proc_out"), ["*_VV_*db.tif"])

for i in range(6, 10):
    dstfile = os.path.join(maindir, site, "{0}_VV_dB_{1}".format(site, i))
    shp = "/geonfs01_vol1/ve39vem/swos_testsites/Test_Sites_project_v29_{0}_{1}.shp".format(site, i)
    if os.path.isfile(dstfile):
        raise IOError("dstfile already exists")

    def seconds(name): return mktime(strptime(re.findall("[0-9T]{15}", name)[0], "%Y%m%dT%H%M%S"))
    srcfiles = sorted(srcfiles, key=seconds)
    groups = []
    temp = []
    for item in srcfiles:
        if len(temp) == 0:
            temp.append(item)
        else:
            if 0 < abs(seconds(item)-seconds(temp[-1])) < 30:
                temp.append(item)
            else:
                groups.append(temp) if len(temp) > 1 else groups.append(temp[0])
                temp = [item]
    stack(groups, dstfile, "bilinear", [20, 20], srcnodata=-99, dstnodata=-99, shapefile=shp, sortfun=seconds, separate=sep)
