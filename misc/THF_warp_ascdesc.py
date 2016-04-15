__author__ = 'John'
import os
from glob import glob
from arcpy import Warp_management

gcp = "D:/HiWi/ForestBiomass/DEM_DSM/THF/SAR/gcp_asc_desc.txt"

source = []
target = []
for line in open(gcp, 'r'):
    items = line.split()
    source.append("'"+items[0]+" "+items[1]+"'")
    target.append("'"+items[2]+" "+items[3]+"'")
source = ";".join(source)
target = ";".join(target)

path = "D:/HiWi/ForestBiomass/DEM_DSM/THF/SAR/asc/"
outpath = "D:/HiWi/ForestBiomass/DEM_DSM/THF/SAR/asc/coreg_desc/"

for image in glob(os.path.join(path, "*.tif")):
    outname = outpath+"/"+os.path.basename(image)
    Warp_management(image, source, target, outname, "POLYORDER3", "BILINEAR")
map(lambda x: os.remove(x), sum([glob(os.path.join(outpath, y)) for y in ["*.xml", "*.ovr", "*.tfw"]], []))