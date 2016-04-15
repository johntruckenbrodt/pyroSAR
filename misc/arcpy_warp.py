# general workflow for batched image warping using arcpy (-> requires ArcGIS-inherent Python interpreter)
# this can easily be extended to using more registration points and higher polynomial order warping
# produced files of type xml, ovr and tfw are considered redundant and will be automatically removed

import os
from glob import glob
from arcpy import Warp_management, env

path = "D:/HiWi/ForestBiomass/DEM_DSM/KAR/SAR/asc"

# enable overwriting
env.overwriteOutput = True

# define warp point for single image acquisition
pt_source = "'456128.283496 5433750.261664'"
pt_target = "'456233.270373 5433747.086658'"

# iterate over all files matching a certain pattern and warp them
for image in glob(os.path.join(path, "*140715.tif")):
    outname = os.path.join(path, os.path.basename(image)[:10]+"_r_temp.tif")
    Warp_management(image, pt_source, pt_target, outname, "POLYORDER0", "BILINEAR")
    # remove redundant files
    map(lambda x: os.remove(x), sum([glob(os.path.join(path, y)) for y in ["*.xml", "*.ovr", "*.tfw"]], []))

# define warp point for total dataset
pt_source = "'456241.346419 5433570.276603'"
pt_target = "'456430.788464 5433604.143338'"

# iterated warping
for image in glob(os.path.join(path, "*.tif")):
    outname = os.path.join(path, os.path.basename(image)[:10]+"_r.tif")
    Warp_management(image, pt_source, pt_target, outname, "POLYORDER0", "BILINEAR")
    map(lambda x: os.remove(x), sum([glob(os.path.join(path, y)) for y in ["*.xml", "*.ovr", "*.tfw"]], []))


# remove temporary files
map(lambda x: os.remove(x), glob(os.path.join(path, "*_temp*")))
