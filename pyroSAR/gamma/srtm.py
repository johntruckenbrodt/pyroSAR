#!/usr/bin/env python
##############################################################
# preparation of srtm data for use in gamma
# module of software pyroSAR
# John Truckenbrodt 2014-18
##############################################################

"""
The following tasks are performed by executing this script:
-reading of a parameter file dem.par
--see object par for necessary values; file is automatically created by starting the script via the GUI
-if necessary, creation of output and logfile directories
-generation of a DEM parameter file for each .hgt (SRTM) file in the working directory or its subdirectories
--the corresponding GAMMA command is create_dem_par, which is interactive. the list variables dempar and dempar2 are piped to the command line for automation
-if multiple files are found, mosaicing is performed
-replacement and interpolation of missing values
-transformation from equiangular (EQA) to UTM projection using a SLC parameter file
"""
import sys

if sys.version_info >= (3, 0):
    from urllib.request import urlopen
else:
    from urllib2 import urlopen

import os
import re
import shutil
import zipfile as zf

from ..spatial.envi import HDRobject, hdr
from ..spatial import raster
from . import ISPPar, process, UTM, slc_corners

import pyroSAR
from pyroSAR.ancillary import finder, run


def fill(dem, dem_out, logpath=None, replace=False):

    width = ISPPar(dem + '.par').width

    path_dem = os.path.dirname(dem_out)

    rpl_flg = 0
    dtype = 4

    # replace values
    value = 0
    new_value = 1
    process(['replace_values', dem, value, new_value, dem + '_temp', width, rpl_flg, dtype], path_dem, logpath)

    value = -32768
    new_value = 0
    process(['replace_values', dem + '_temp', value, new_value, dem + '_temp2', width, rpl_flg, dtype], path_dem, logpath)

    # interpolate missing values
    r_max = 9
    np_min = 40
    np_max = 81
    w_mode = 2
    process(['interp_ad', dem + '_temp2', dem_out, width, r_max, np_min, np_max, w_mode, dtype], path_dem, logpath)

    # remove temporary files
    os.remove(dem+'_temp')
    os.remove(dem+'_temp2')

    # duplicate parameter file for newly created dem
    shutil.copy(dem+'.par', dem_out+'.par')

    # create ENVI header file
    hdr(dem_out+'.par')

    if replace:
        for item in [dem+x for x in ['', '.par', '.hdr', '.aux.xml'] if os.path.isfile(dem+x)]:
            os.remove(item)


def transform(infile, outfile, posting=90):
    """
    transform SRTM DEM from EQA to UTM projection
    """
    # read DEM parameter file
    par = ISPPar(infile + '.par')

    # transform corner coordinate to UTM
    utm = UTM(infile + '.par')

    for item in [outfile, outfile+'.par']:
        if os.path.isfile(item):
            os.remove(item)

    # determine false northing from parameter file coordinates
    falsenorthing = 10000000. if par.corner_lat < 0 else 0

    # create new DEM parameter file with UTM projection details
    inlist = ['UTM', 'WGS84', 1, utm.zone, falsenorthing, os.path.basename(outfile), '', '', '', '', '', '-{0} {0}'.format(posting), '']
    process(['create_dem_par', outfile + '.par'], inlist=inlist)

    # transform dem
    process(['dem_trans', infile + '.par', infile, outfile + '.par', outfile, '-', '-', '-', 1])
    hdr(outfile+'.par')


def dempar(dem, logpath=None):
    """
    create GAMMA parameter text files for DEM files
    currently only EQA and UTM projections with WGS84 ellipsoid are supported
    """
    rast = raster.Raster(dem)

    # determine data type
    dtypes = {'Int16': 'INTEGER*2', 'UInt16': 'INTEGER*2', 'Float32': 'REAL*4'}
    if rast.dtype not in dtypes:
        raise IOError('data type not supported')
    else:
        dtype = dtypes[rast.dtype]

    # format pixel posting and top left coordinate
    posting = str(rast.geo['yres'])+' '+str(rast.geo['xres'])
    latlon = str(rast.geo['ymax'])+' '+str(rast.geo['xmin'])

    # evaluate projection
    projections = {'longlat': 'EQA', 'utm': 'UTM'}
    if rast.proj4args['proj'] not in projections:
        raise ValueError('projection not supported (yet)')
    else:
        projection = projections[rast.proj4args['proj']]

    # get ellipsoid
    ellipsoid = rast.proj4args['ellps'] if 'ellps' in rast.proj4args else rast.proj4args['datum']
    if ellipsoid != 'WGS84':
        raise ValueError('ellipsoid not supported (yet)')

    # create list for GAMMA command input
    if projection == 'UTM':
        zone = rast.proj4args['zone']
        falsenorthing = 10000000. if rast.geo['ymin'] < 0 else 0
        parlist = [projection, ellipsoid, 1, zone, falsenorthing, os.path.basename(dem), dtype, 0, 1, rast.cols, rast.rows, posting, latlon]
    else:
        parlist = [projection, ellipsoid, 1, os.path.basename(dem), dtype, 0, 1, rast.cols, rast.rows, posting, latlon]

    # execute GAMMA command
    process(['create_dem_par', os.path.splitext(dem)[0] + '.par'], os.path.dirname(dem), logpath, inlist=parlist)


def swap(data, outname):
    """
    byte swapping from small to big endian (as required by GAMMA)
    """
    rast = raster.Raster(data)
    dtype = rast.dtype
    if rast.format != 'ENVI':
        raise IOError('only ENVI format supported')
    dtype_lookup = {'Int16': 2, 'CInt16': 2, 'Int32': 4, 'Float32': 4, 'CFloat32': 4, 'Float64': 8}
    if dtype not in dtype_lookup:
        raise IOError('data type {} not supported'.format(dtype))
    process(['swap_bytes', data, outname, str(dtype_lookup[dtype])])
    header = HDRobject(data+'.hdr')
    header.byte_order = 1
    hdr(header, outname+'.hdr')


def mosaic(demlist, outname, byteorder=1, gammapar=True):
    """
    mosaicing of multiple DEMs
    """
    if len(demlist) < 2:
        raise IOError('length of demlist < 2')
    nodata = str(raster.Raster(demlist[0]).nodata)
    run(['gdalwarp', '-q', '-of', 'ENVI', '-srcnodata', nodata, '-dstnodata', nodata, demlist, outname])
    if byteorder == 1:
        swap(outname, outname+'_swap')
        for item in [outname, outname+'.hdr', outname+'.aux.xml']:
            os.remove(item)
        os.rename(outname+'_swap', outname)
        os.rename(outname+'_swap.hdr', outname+'.hdr')
    if gammapar:
        dempar(outname)


def hgt(parfiles):
    """
    concatenate hgt file names overlapping with multiple SAR scenes
    input is a list of GAMMA SAR scene parameter files
    this list is read for corner coordinates of which the next integer lower left latitude and longitude is computed
    hgt files are supplied in 1 degree equiangular format named e.g. N16W094.hgt (with pattern [NS][0-9]{2}[EW][0-9]{3}.hgt
    For north and east hemisphere the respective absolute latitude and longitude values are smaller than the lower left coordinate of the SAR image
    west and south coordinates are negative and hence the nearest lower left integer absolute value is going to be larger
    """

    lat = []
    lon = []
    for parfile in parfiles:
        if isinstance(parfile, pyroSAR.ID):
            corners = parfile.getCorners()
        elif parfile.endswith('.par'):
            corners = slc_corners(parfile)
        lat += [int(float(corners[x]) // 1) for x in ['ymin', 'ymax']]
        lon += [int(float(corners[x]) // 1) for x in ['xmin', 'xmax']]

    # add missing lat/lon values (and add an extra buffer of one degree)
    lat = range(min(lat), max(lat)+1)
    lon = range(min(lon), max(lon)+1)

    # convert coordinates to string with leading zeros and hemisphere identification letter
    lat = [str(x).zfill(2+len(str(x))-len(str(x).strip('-'))) for x in lat]
    lat = [x.replace('-', 'S') if '-' in x else 'N'+x for x in lat]

    lon = [str(x).zfill(3+len(str(x))-len(str(x).strip('-'))) for x in lon]
    lon = [x.replace('-', 'W') if '-' in x else 'E'+x for x in lon]

    # concatenate all formatted latitudes and longitudes with each other as final product
    return [x+y+'.hgt' for x in lat for y in lon]


def makeSRTM(scenes, srtmdir, outname):
    """
    Create a DEM from SRTM tiles
    Input is a list of pyroSAR.ID objects from which coordinates are read to determine the required DEM extent
    Mosaics SRTM DEM tiles, converts them to Gamma format and subtracts offset to WGS84 ellipsoid
    for DEMs downloaded from USGS http://gdex.cr.usgs.gov or CGIAR http://srtm.csi.cgiar.org
    """

    tempdir = outname+'___temp'
    os.makedirs(tempdir)

    hgt_options = hgt(scenes)

    hgt_files = finder(srtmdir, hgt_options)

    # todo: check if really needed
    nodatas = [str(int(raster.Raster(x).nodata)) for x in hgt_files]

    srtm_vrt = os.path.join(tempdir, 'srtm.vrt')
    srtm_temp = srtm_vrt.replace('.vrt', '_tmp')
    srtm_final = srtm_vrt.replace('.vrt', '')

    run(['gdalbuildvrt', '-overwrite', '-srcnodata', ' '.join(nodatas), srtm_vrt, hgt_files])

    run(['gdal_translate', '-of', 'ENVI', '-a_nodata', -32768, srtm_vrt, srtm_temp])

    process(['srtm2dem', srtm_temp, srtm_final, srtm_final + '.par', 2, '-'], outdir=tempdir)

    shutil.move(srtm_final, outname)
    shutil.move(srtm_final+'.par', outname+'.par')
    hdr(outname+'.par')

    shutil.rmtree(tempdir)


def hgt_collect(parfiles, outdir, demdir=None, arcsec=3):
    """
    automatic downloading and unpacking of srtm tiles
    base directory must contain SLC files in GAMMA format including their parameter files for reading coordinates
    additional dem directory may locally contain srtm files. This directory is searched for locally existing files, which are then copied to the current working directory
    """

    # concatenate required hgt tile names
    target_ids = hgt(parfiles)

    targets = []

    pattern = '[NS][0-9]{2}[EW][0-9]{3}'

    # if an additional dem directory has been defined, check this directory for required hgt tiles
    if demdir is not None:
        targets.extend(finder(demdir, target_ids))

    # check for additional potentially existing hgt tiles in the defined output directory
    extras = [os.path.join(outdir, x) for x in target_ids if os.path.isfile(os.path.join(outdir, x)) and not re.search(x, '\n'.join(targets))]
    targets.extend(extras)

    print('found {} relevant SRTM tiles...'.format(len(targets)))

    # search server for all required tiles, which were not found in the local directories
    if len(targets) < len(target_ids):
        print('searching for additional SRTM tiles on the server...')
        onlines = []

        if arcsec == 1:
            remotes = ['http://e4ftl01.cr.usgs.gov/SRTM/SRTMGL1.003/2000.02.11/']
            remotepattern = pattern+'.SRTMGL1.hgt.zip'
        elif arcsec == 3:
            server = 'http://dds.cr.usgs.gov/srtm/version2_1/SRTM3/'
            remotes = [os.path.join(server, x) for x in ['Africa', 'Australia', 'Eurasia', 'Islands', 'North_America', 'South_America']]
            remotepattern = pattern+'[.]hgt.zip'
        else:
            raise ValueError('argument arcsec must be of value 1 or 3')

        for remote in remotes:
            response = urlopen(remote).read()
            items = sorted(set(re.findall(remotepattern, response)))
            for item in items:
                outname = re.findall(pattern, item)[0]+'.hgt'
                if outname in target_ids and outname not in [os.path.basename(x) for x in targets]:
                    onlines.append(os.path.join(remote, item))

        # if additional tiles have been found online, download and unzip them to the local directory
        if len(onlines) > 0:
            print('downloading {} SRTM tiles...'.format(len(onlines)))
            for candidate in onlines:
                localname = os.path.join(outdir, re.findall(pattern, candidate)[0]+'.hgt')
                infile = urlopen(candidate)
                with open(localname+'.zip', 'wb') as outfile:
                    outfile.write(infile.read())
                infile.close()
                with zf.ZipFile(localname+'.zip', 'r') as z:
                    z.extractall(outdir)
                os.remove(localname+'.zip')
                targets.append(localname)
    return targets


if __name__ == '__main__':
    main()
