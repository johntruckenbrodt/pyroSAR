##############################################################
# demonstration script for resampling, stacking, mosaicking and subsetting SAR images
# John Truckenbrodt 2017
##############################################################
import os

from pyroSAR.ancillary import groupbyTime, seconds
from spatialist import stack
from spatialist.ancillary import finder


def main():
    # define input directory containing file sto be stacked
    dir_in = '/...'
    
    # define output file name
    dstfile = '/.../x'
    
    # shapefile (for stack boundaries)
    shp = '/../x.shp'
    
    # store results in separate files or one single stack file? If separate then dstfile is used as a directory.
    sep = True
    
    # list files to be resampled; those not overlapping with the shapefile geometry will excluded by function stack
    srcfiles = finder(dir_in, ['S1*_VV_*norm_db.tif'])
    
    # check whether dstfile is already a file
    if os.path.isfile(dstfile):
        raise IOError('dstfile already exists')
    
    # create groups of similar time stamps for mosaicking.
    # All images with a time stamp of less than 30s difference will be grouped
    groups = groupbyTime(srcfiles, seconds, 30)
    
    # final function call
    # groups will be mosaicked first
    # the resulting images will all have the same extent
    stack(srcfiles=groups, dstfile=dstfile, resampling='bilinear',
          targetres=[20, 20], srcnodata=-99, dstnodata=-99,
          shapefile=shp, sortfun=seconds, separate=sep, overwrite=False)


if __name__ == '__main__':
    main()
