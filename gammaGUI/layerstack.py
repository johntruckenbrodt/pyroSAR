##############################################################
# crop, stack and resample multiple raster layers
# module of software gammaGUI
# John Truckenbrodt 2015
##############################################################

"""
The following tasks are performed by executing this script:
-search raster layers by single or multiple search expression(s)
-stack raster layers into a single envi file including options to set:
--target resolution
--extent cropping by providing a rectangular shapefile
--resampling method
--NA value
-if options are left blank the defaults are taken
-see documentation of GDAL's gdalbuildvrt and gdalwarp for details
"""

import sys

import os
import re
import raster
import subprocess as sp

from ancillary import dissolve, finder
from envi import HDRobject, hdr

[path_in, file_out, shape, pattern, reg, resampling, targetres, nodata] = sys.argv[1:]

# find all files matching the defined pattern(s)
items = finder(path_in, pattern.split(", "), regex=True if reg == "True" else False)

# if a shapefile has been defined, retrieve its bounding box extent
if len(shape) > 0:
    if not re.search(".shp$", shape):
        raise IOError("no valid shapefile defined")
    shape_info = sp.Popen(["ogrinfo", "-ro", "-al", shape], stdout=sp.PIPE).stdout.read().split("\n")
    shape_extent = [re.findall("[0-9.]+", x) for x in shape_info if re.search("Extent", x)][0]
    extent = ["-te", shape_extent]

    # reduce list of rasters to those overlapping with the shapefile bounding box
    items = [x for x in items if raster.intersect(raster.Raster(x), raster.Extent(shape_extent))]
else:
    extent = []

# evaluate target resolution statement
if len(targetres) > 0:
    targetres = filter(None, re.split("[, ]", targetres))
    if len(targetres) == 2:
        try:
            int(targetres[0])
            int(targetres[1])
            targetres = ["-tr", targetres]
        except IOError:
            print "invalid resolution statement"
    else:
        raise IOError("invalid resolution statement")
else:
    targetres = []

nodata = ["-dstnodata", nodata]
format = ["-of", "ENVI"]
resampling = ["-r", resampling]


if len(items) > 0:
    path_out = os.path.dirname(file_out)
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    print "the following files will be stacked to file {0}:".format(file_out)
    for item in items:
        print item
    decision = raw_input("proceed (y/n)?: ")
    if decision == "y":
        vrt = file_out+".vrt"
        sp.check_call(dissolve(["gdalbuildvrt", "-q", "-overwrite", "-separate", extent, vrt, items]))
        sp.check_call(dissolve(["gdalwarp", "-q", "-overwrite", resampling, format, nodata, targetres, vrt, file_out]))
        os.remove(vrt)

        # add band names (basenames of the original files without extension) to the envi header file
        par = HDRobject(file_out+".hdr")
        par.band_names = [os.path.splitext(os.path.basename(x))[0] for x in items]
        hdr(par)
    else:
        print "command aborted"
else:
    print "no files matching the defined pattern(s) found"
