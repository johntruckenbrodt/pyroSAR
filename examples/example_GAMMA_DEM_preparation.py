##############################################################
# demonstration script for SRTM DEM preparation for use with GAMMA
# John Truckenbrodt 2017-2019
##############################################################

from pyroSAR.gamma.dem import dem_autocreate
from spatialist import Vector


def main():
    # define a shapefile for the study area
    shp = '/.../testsite.shp'
    
    # define the output name of the DEM (no file extension like .tif etc.!)
    outname = '/path/to/demfile'
    
    # define a buffer around the study area boundaries (in degrees)
    buffer = 0.01
    
    # load the defined shapefile
    with Vector(shp) as site:
        
        # reproject the shapefile to latlon (in-memory, no file written or modified)
        site.reproject(4326)
        
        # create the DEM mosaic in Gamma format including tile download, warping to arbitrary CRS and resolution,
        # and conversion from geoid to ellipsoid heights
        # documentation: https://pyrosar.readthedocs.io/en/latest/pyroSAR.html#pyroSAR.gamma.dem.dem_autocreate
        dem_autocreate(geometry=site,
                       outfile=outname,
                       buffer=buffer,
                       demType='SRTM 1Sec HGT',
                       t_srs=32632,  # WGS84, UTM 32N
                       tr=(20, 20),  # 20 m target resolution in both x and y
                       geoid_mode='gamma',  # convert from EGM96 geoid to WGS84 ellipsoid heights using Gamma
                       resampling_method='bilinear')


if __name__ == '__main__':
    main()
