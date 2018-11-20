import os
import socket
from scoop import futures

from pyroSAR.snap import geocode
from pyroSAR import Archive

from spatialist import vector
from spatialist.ancillary import finder

"""
This script is an example usage for processing Sentinel-1 scenes with SNAP

Run this script by calling the 'start_gamma.sh' scipt.

The following tasks are performed:
- a directory is scanned for valid Sentinel-1 scenes
- the found scenes are ingested into a spatialite database
- a cluster job is setup using package 'scoop', which assigns a list of testsites to different cluster nodes
- for each site:
    - query the SAR scenes, which overlap with your testsite and match certain criteria (e.g. sensor, acquisition mode etc.)
    - filter the selected scenes by those that have already been processed and saved to the defined output directory
    - run processing using ESA SNAP's Graph Processing Tool (GPT)
"""

# the sites to be processed
# this is just an exemplary use case assuming a shapefile with different geometries for the test sites
sites = ['Egypt_Burullus', 'France_Camargue', 'Kenya_Lorian_Olbolossat', 'Sweden_Skogaryd', 'Sweden_Store-Mosse']

# the pyroSAR database file
dbfile = '/.../scenelist.db'

# the main directory for storing the processed results
maindir = '/.../swos_process'


def worker(sitename):
    #######################################################################################
    # setup general processing parameters

    resolution = 20
    #######################################################################################
    # define the directories for writing temporary and final results
    sitedir = os.path.join(maindir, sitename)
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
                                        sensor=('S1A', 'S1B'),
                                        product='GRD',
                                        acquisition_mode='IW',
                                        vv=1)

    print('{0}: {1} scenes found for site {2}'.format(socket.gethostname(), len(selection_proc), sitename))
    #######################################################################################
    # call to processing utility
    if len(selection_proc) > 1:
        print('start processing')

        for scene in selection_proc:
            geocode(infile=scene, outdir=outdir, tr=resolution, scaling='db')
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
    # start the processing
    results = list(futures.map(worker, sites))
