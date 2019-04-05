##############################################################
# demonstration script for SRTM DEM preparation for use with GAMMA
# John Truckenbrodt 2017-2019
##############################################################

import os
import shutil

from pyroSAR import gamma
from pyroSAR.gamma import dem, par2hdr
from spatialist import vector

from spatialist.ancillary import run


def main():
    # define a shapefile for the study area
    shp = '/.../testsite.shp'
    
    # define the output name of the DEM (no file extension like .tif etc.!)
    outname = '/path/to/demfile'
    
    # define a buffer around the study area boundaries (in degrees)
    buffer = 0.01
    
    # load the defined shapefile
    with vector.Vector(shp) as site:
        # reproject the shapefile to latlon (in-memory, no file written or modified)
        site.reproject('+proj=longlat +datum=WGS84 +no_defs ')

        # extract the extent (bounding box) of the shapefile
        ext = site.extent
    
    # add the buffer to the bounding box
    ext['xmin'] -= buffer
    ext['ymin'] -= buffer
    ext['xmax'] += buffer
    ext['ymax'] += buffer
    
    # define a GDAL VRT file containing all SRTM tiles
    # this file has all hgt tiles in the same directory registered and is used for subsetting/mosaicing
    srtm_vrt = '/path/to/SRTM_1_HGT.vrt'
    
    # create a temporary directory for writing intermediate files (will be deleted at the end)
    tmpdir = outname + '__tmp'
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)
    
    # define a name for a temporary DEM file
    dem_tmp = os.path.join(tmpdir, 'srtm_tmp.tif')
    
    # create a DEM mosaic for the study site
    run(['gdalwarp', '-q', '-of', 'GTiff', '-te', ext['xmin'], ext['ymin'], ext['xmax'], ext['ymax'], srtm_vrt,
         dem_tmp])
    
    # transform the DEM to GAMMA format (including EGM96 geoid to WGS84 ellipsoid height reference correction)
    gamma.process(['srtm2dem', dem_tmp, outname, outname + '.par', 2, '-'], outdir=tmpdir)
    
    # create an ENVI header file
    par2hdr(outname + '.par', outname + '.hdr')
    
    # remove the temporary directory with all intermediate files
    shutil.rmtree(tmpdir)
    
    # optional: transform DEM to UTM projection
    # the UTM zone is automatically computed for the center of the DEM file
    dem.transform(outname, outname + '_utm', posting=20)


if __name__ == '__main__':
    main()
