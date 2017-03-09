"""
Demonstration script for the usage of the pyroSAR processing with gamma
"""
import os
from datetime import datetime
from ancillary import finder, multicore
from gamma.util import geocode
from pyroSAR import identify

#Path to the SRTM in the GAMMA format
srtm = "/path/to/srtmfile" 

#Path to the temporary folder, this will be emptied after the processing
tempdir = "/path/to/temporary" 	

#Path to the output folder
outdir = "/path/to/output" 

#path to the raw data folder
datadir = "/path/to/rawdata"

#Target resolution
targetres = 10

#Scaling of the output data could be 'linear' or 'db' or a list with both
scaling = ['linear','db']

# Should only one scene be processed? The temporal data would be kept. This is for testing purposes.
one_scene_test = False

################################################
# No further user input needed beyond this point
# ##############################################


func_geoback = 2
func_interp = 0


for dir in [tempdir, outdir]:
    if not os.path.isdir(dir):
        os.makedirs(dir)
        
scenes = finder(datadir, ["S1*GRDH*zip"])

selection = []
a = datetime.now()
for scene in scenes:
    try:
        id = identify(scene)
        if len(finder(outdir, [id.outname_base()], regex=True)) == 0:
           selection.append(id)
    except IOError:
        continue

b = datetime.now()
print b - a

with open(os.path.join(outdir, "scenelist"), "w") as outfile:
    outfile.write("\n".join([x.scene for x in selection]) + "\n")

scene = selection[0]
if one_scene_test:
    geocode(scene, srtm, tempdir, outdir, targetres, scaling, func_geoback, func_interp, sarsimulation=False, cleanup=False)
    exit()

print "start geocoding"
print datetime.now()
multicore(geocode, cores=4, multiargs={"scene": selection}, dem=srtm,
          tempdir=tempdir, outdir=outdir,
          targetres=targetres, scaling=scaling,
          func_geoback=func_geoback, func_interp=func_interp,
          sarsimulation=False,cleanup=True)
print datetime.now()
print "done"  
