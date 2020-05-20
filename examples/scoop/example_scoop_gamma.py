import os
import socket
from scoop import futures

from pyroSAR.S1 import OSV
from pyroSAR.gamma import geocode
from pyroSAR import Archive

from spatialist import vector
from spatialist.ancillary import finder, multicore

"""
This script is an example usage for processing Sentinel-1 scenes with GAMMA

Run this script by calling the 'start_gamma.sh' script.

The following tasks are performed:
- a directory is scanned for valid Sentinel-1 scenes
- the found scenes are ingested into a spatialite database
- orbit state vector (OSV) files are downloaded to a user-defined directory (these are needed for precise orbit information)
    - currently this is implemented to update a fixed directory in which all OSV files are stored
    - an empty directory will first be filled with all available OSV files on the server
- a cluster job is setup using package 'scoop', which assigns a list of testsites to different cluster nodes
- for each site:
    - query the SAR scenes, which overlap with your testsite and match certain criteria (e.g. sensor, acquisition mode etc.)
    - filter the selected scenes by those that have already been processed and saved to the defined output directory
    - do parallelized processing using package 'pathos'
"""

# the sites to be processed
# this is just an exemplary use case assuming a shapefile with different geometries for the test sites
sites = ['Egypt_Burullus', 'France_Camargue', 'Kenya_Lorian_Olbolossat', 'Sweden_Skogaryd', 'Sweden_Store-Mosse']

# the pyroSAR database file
dbfile = '/.../scenelist.db'

# the main directory for storing the processed results
maindir = '/.../swos_process'

# the directories for Sentinel-1 POE and RES orbit state vector files
# this is intended to be a fixed directory structure similar to that of ESA SNAP
# in the future all auxiliary data files will be stored in a structure defined by pyroSAR
# This is currently only used for Sentinel-1; within the processor two subdirectories will be created named
# POEORB and RESORB which will contain the respective orbit files
osvdir = '/.../.gamma/auxdata/Orbits/Sentinel-1'


def worker(sitename):
    #######################################################################################
    # setup general processing parameters

    resolution = 20

    # number of processes for Python pathos framework (multiple scenes in parallel)
    parallel1 = 6

    # number of parallel OpenMP threads; this is used by GAMMA internally
    parallel2 = 6
    os.environ['OMP_NUM_THREADS'] = str(parallel2)
    #######################################################################################
    # get the maximum date of the precise orbit files
    # as type also 'RES' can be selected. These files are not as precise as POE and thus geocoding might not be
    # quite as accurate
    with OSV(osvdir) as osv:
        maxdate = osv.maxdate(osvtype='POE', datetype='stop')
    #######################################################################################
    # define the directories for writing temporary and final results
    sitedir = os.path.join(maindir, sitename)
    tempdir = os.path.join(sitedir, 'proc_in')
    outdir = os.path.join(sitedir, 'proc_out')
    #######################################################################################
    # load the test site geometry into a vector object
    sites = vector.Vector('/.../testsites.shp')

    # query the test site by name; a column name 'Site_Name' must be saved in your shapefile
    site = sites["Site_Name='{}'".format(sitename)]
    #######################################################################################
    # query the database for scenes to be processed
    with Archive(dbfile) as archive:
        selection_proc = archive.select(vectorobject=site,
                                        processdir=outdir,
                                        maxdate=maxdate,
                                        sensor=('S1A', 'S1B'),
                                        product='GRD',
                                        acquisition_mode='IW',
                                        vv=1)

    print('{0}: {1} scenes found for site {2}'.format(socket.gethostname(), len(selection_proc), sitename))
    #######################################################################################
    # define the DEM file
    demfile = '{0}/{1}/DEM/{1}_srtm_utm'.format(maindir, sitename)
    if not os.path.isfile(demfile):
        print('DEM missing for site {}'.format(sitename))
        return
    #######################################################################################
    # call to processing utility
    if len(selection_proc) > 1:
        print('start processing')
    if len(selection_proc) > 1:
        if len(selection_proc) < parallel1:
            parallel1 = len(selection_proc)
        # run the function on multiple cores in parallel
        multicore(geocode, cores=parallel1, multiargs={'scene': selection_proc}, dem=demfile,
                  tempdir=tempdir, outdir=outdir,
                  targetres=resolution, scaling='db',
                  func_geoback=2, func_interp=0, sarSimCC=False, osvdir=osvdir, cleanup=True, allow_RES_OSV=False)
    elif len(selection_proc) == 1:
        scene = selection_proc[0]
        # run the function on a single core
        geocode(scene, dem=demfile,
                tempdir=tempdir, outdir=outdir,
                targetres=resolution, scaling='db',
                func_geoback=2, func_interp=0, sarSimCC=False, osvdir=osvdir, cleanup=True, allow_RES_OSV=False)
    return len(selection_proc)


if __name__ == '__main__':
    #######################################################################################
    # update Sentinel-1 GRD scene archive database

    # define a directory containing zipped scene archives and list all files starting with 'S1A' or 'S1B'
    archive_s1 = '/.../sentinel1/GRD'
    scenes_s1 = finder(archive_s1, ['^S1[AB]'], regex=True, recursive=False)

    with Archive(dbfile) as archive:
        archive.insert(scenes_s1)
    #######################################################################################
    # download the latest orbit state vector files
    with OSV(osvdir) as osv:
        osv.update()
    #######################################################################################
    # start the processing
    results = list(futures.map(worker, sites))
