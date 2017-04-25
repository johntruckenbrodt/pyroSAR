"""
Demonstration script for the usage of the pyroSAR processing with gamma
"""
import os
from spatial import vector
from datetime import datetime
from ancillary import finder, multicore
from gamma import geocode
from pyroSAR import identify, Archive

# Path to the SRTM in the GAMMA format
srtm = '/path/to/srtmfile'

# Path to the temporary folder, this will be emptied after the processing
tempdir = '/path/to/temporary'

# Path to the output folder
outdir = '/path/to/output'

# path to the raw data folder
datadir = '/path/to/rawdata'

# shapefile of the study area (same extent as the DEM file)
shp = '/path/to/shapefile'

# Target resolution in meters
targetres = 10

# scene archive csv textfile (here all filenames and metadata of the SAR scenes in datadir will be written to)
archive_file = '/path/to/archive.csv'

# search pattern for all the SAR scenes of interest in datadir
pattern = 'S1*GRDH*zip'

# Scaling of the output data could be 'linear' or 'db' or a list with both
scaling = ['linear', 'db']

# Should only one scene be processed? The temporal data would be kept. This is for testing purposes.
one_scene_test = True

# backward geocoding interpolation mode (see GAMMA command geocode_back or call 'help(geocode)')
func_geoback = 2

# output lookup table values in regions of layover, shadow, or DEM gaps (see GAMMA command gc_map or call 'help(geocode)')
func_interp = 0

# create the define directories of they do not yet exist
for dir in [tempdir, outdir]:
    if not os.path.isdir(dir):
        os.makedirs(dir)

# select the scenes in the datadir based on the define search pattern
scenes = finder(datadir, [pattern])

# load the shapefile
shpobj = vector.Vector(shp)

# create the archive and update it with scenes (the archive csv will either be newly created or updated with the file names which had not been registered before)
with Archive(archive_file, header=True) as archive:
    archive.update(scenes)
    # select all scenes which overlap with the shapefile of the study area (reprojection is performed automatically in-memory)
    selection_site = archive.select(shpobj)

# load all selected scenes to pyroSAR metadata objects (this might take a little time)
selection_proc = map(identify, selection_site)

# filter the selected scenes by additional acquisition parameters
# call 'selection_proc[0].summary()' to get an overview of all object parameter names which can be selected
selection_proc = [x for x in selection_proc if x.acquisition_mode == 'IW' and 'VV' in x.polarizations]

# process a single scene
scene = selection_proc[0]
if one_scene_test:
    geocode(scene, srtm, tempdir, outdir, targetres, scaling, func_geoback, func_interp, sarsimulation=False,
            cleanup=False)
    exit()

# process all scenes in parallel
# !!! WARNING: always try processing a single scene first to prevent system crashes!!!
# call 'help(geocode)' for a documentation and description of input parameters
print 'start geocoding'
print datetime.now()
multicore(geocode, cores=4, multiargs={'scene': selection_proc}, dem=srtm,
          tempdir=tempdir, outdir=outdir,
          targetres=targetres, scaling=scaling,
          func_geoback=func_geoback, func_interp=func_interp,
          sarsimulation=False, cleanup=True)
print datetime.now()
print 'done'
