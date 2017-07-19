##############################################################
# demonstration script for resampling, stacking, mosaicking and subsetting SAR images
# John Truckenbrodt 2017
##############################################################
import os

from pyroSAR.ancillary import finder, groupbyTime, seconds
from pyroSAR.spatial import stack, Vector


def main():
    # define input directory containing file sto be stacked
    dir_in = '/...'

    # define output file name
    dstfile = '/.../x'

    # shapefile (for stack boundaries)
    shp = '/../x.shp'

    # store results in separate files or one single stack file? If separate then dstfile is used as a directory.
    sep = True

    # define
    srcfiles = finder(dir_in, ['S1*_VV_*norm_db.tif'])

    if os.path.isfile(dstfile):
        raise IOError('dstfile already exists')

    site = Vector(shp)

    # create groups of similar time stamps for mosaicking. All images with a time stamp of less than 30s difference will be mosaicked
    groups = groupbyTime(srcfiles, seconds, 30)

    # final function call
    stack(srcfiles=groups, dstfile=dstfile, resampling='bilinear',
          targetres=[20, 20], srcnodata=-99, dstnodata=-99,
          shapefile=site, sortfun=seconds, separate=sep, overwrite=False)

if __name__ == '__main__':
    main()
