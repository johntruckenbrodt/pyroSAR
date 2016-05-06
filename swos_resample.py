
import re
from raster import stack
from ancillary import finder
from time import mktime, strptime

srcfiles = finder("/geonfs01_vol1/ve39vem/out", ["*_VV_*db.tif"])
dstfile = "/geonfs01_vol1/ve39vem/Camarque_VV_dB_update"
shp = "/geonfs01_vol1/ve39vem/S1/test_camarque/boundary/Test_Sites_project_v11_Camarque.shp"


def seconds(name): return mktime(strptime(re.findall("[0-9T]{15}", name)[0], "%Y%m%dT%H%M%S"))

srcfiles = sorted(srcfiles, key=seconds)

groups = []
temp = []

for item in srcfiles:
    if len(temp) == 0:
        temp.append(item)
    else:
        if abs(seconds(item)-seconds(temp[-1])) < 30:
            temp.append(item)
        else:
            groups.append(temp) if len(temp) > 1 else groups.append(temp[0])
            temp = [item]

stack(groups, dstfile, "bilinear", [20, 20], srcnodata=-99, dstnodata=-99, shapefile=shp, sortfun=seconds)
