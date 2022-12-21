###############################################################################
# preparation of DEM data for use in GAMMA

# Copyright (c) 2014-2022, the pyroSAR Developers.

# This file is part of the pyroSAR Project. It is subject to the
# license terms in the LICENSE.txt file found in the top-level
# directory of this distribution and at
# https://github.com/johntruckenbrodt/pyroSAR/blob/master/LICENSE.txt.
# No part of the pyroSAR project, including this file, may be
# copied, modified, propagated, or distributed except according
# to the terms contained in the LICENSE.txt file.
################################################################################

"""
A collection of functions to handle digital elevation models in GAMMA
"""
from urllib.request import urlopen
import os
import re
import shutil
import zipfile as zf

from spatialist import raster, gdal_translate, gdalbuildvrt, gdalwarp, crsConvert
from spatialist.ancillary import finder
from spatialist.envi import HDRobject

from ..auxdata import dem_autoload, dem_create
from ..drivers import ID
from . import ISPPar, UTM, slc_corners, par2hdr
from pyroSAR.examine import ExamineGamma
from pyroSAR.ancillary import hasarg

import logging

log = logging.getLogger(__name__)

try:
    from .api import diff, disp, isp
except ImportError:
    pass


def fill(dem, dem_out, logpath=None, replace=False):
    """
    interpolate missing values in the SRTM DEM (value -32768)

    Parameters
    ----------
    dem: str
        the input DEM to be filled
    dem_out: str
        the name of the filled DEM
    logpath: str
        a directory to write logfiles to
    replace: bool
        delete `dem` once finished?

    Returns
    -------

    """
    width = ISPPar(dem + '.par').width
    
    path_dem = os.path.dirname(dem_out)
    
    rpl_flg = 0
    dtype = 4
    
    # replace values
    value = 0
    new_value = 1
    disp.replace_values(f_in=dem,
                        value=value,
                        new_value=new_value,
                        f_out=dem + '_temp',
                        width=width,
                        rpl_flg=rpl_flg,
                        dtype=dtype,
                        logpath=logpath)
    
    value = -32768
    new_value = 0
    disp.replace_values(f_in=dem + '_temp',
                        value=value,
                        new_value=new_value,
                        f_out=dem + '_temp2',
                        width=width,
                        rpl_flg=rpl_flg,
                        dtype=dtype,
                        outdir=path_dem,
                        logpath=logpath)
    
    # interpolate missing values
    isp.interp_ad(data_in=dem + '_temp2',
                  data_out=dem_out,
                  width=width,
                  r_max=9,
                  np_min=40,
                  np_max=81,
                  w_mode=2,
                  dtype=dtype,
                  outdir=path_dem,
                  logpath=logpath)
    
    # remove temporary files
    os.remove(dem + '_temp')
    os.remove(dem + '_temp2')
    
    # duplicate parameter file for newly created dem
    shutil.copy(dem + '.par', dem_out + '.par')
    
    # create ENVI header file
    par2hdr(dem_out + '.par', dem_out + '.hdr')
    
    if replace:
        for item in [dem + x for x in ['', '.par', '.hdr', '.aux.xml'] if os.path.isfile(dem + x)]:
            os.remove(item)


def transform(infile, outfile, posting=90):
    """
    transform SRTM DEM from EQA to UTM projection
    """
    # read DEM parameter file
    par = ISPPar(infile + '.par')
    
    # transform corner coordinate to UTM
    utm = UTM(infile + '.par')
    
    for item in [outfile, outfile + '.par']:
        if os.path.isfile(item):
            os.remove(item)
    
    # determine false northing from parameter file coordinates
    falsenorthing = 10000000. if par.corner_lat < 0 else 0
    
    # create new DEM parameter file with UTM projection details
    inlist = ['UTM', 'WGS84', 1, utm.zone, falsenorthing, os.path.basename(outfile), '', '', '', '', '',
              '-{0} {0}'.format(posting), '']
    
    diff.create_dem_par(DEM_par=outfile + '.par',
                        inlist=inlist)
    
    # transform dem
    diff.dem_trans(DEM1_par=infile + '.par',
                   DEM1=infile,
                   DEM2_par=outfile + '.par',
                   DEM2=outfile,
                   bflg=1)
    par2hdr(outfile + '.par', outfile + '.hdr')


def dem_autocreate(geometry, demType, outfile, buffer=None, t_srs=4326, tr=None, logpath=None,
                   username=None, password=None, geoid_mode='gamma', resampling_method='bilinear'):
    """
    | automatically create a DEM in GAMMA format for a defined spatial geometry.
    | The following steps will be performed:

    - collect all tiles overlapping with the geometry using :func:`pyroSAR.auxdata.dem_autoload`

      * if they don't yet exist locally they will automatically be downloaded
      * the tiles will be downloaded into the SNAP auxdata directory structure,
        e.g. ``$HOME/.snap/auxdata/dem/SRTM 3Sec``

    - create a mosaic GeoTIFF of the same spatial extent as the input geometry
      plus a defined buffer using :func:`pyroSAR.auxdata.dem_create`
    
    - if necessary, subtract the geoid-ellipsoid difference (see :func:`pyroSAR.auxdata.dem_autoload`
      for height references of different supported DEMs)
    
    - convert the result to GAMMA format
    
      * If ``t_srs`` is `4326` and the DEM's height reference is either `WGS84` ellipsoid or `EGM96` geoid,
        the command ``srtm2dem`` can be used. This is kept for backwards compatibility.
      * For all other cases the newer command ``dem_import`` can be used if it exists and if the command
        ``create_dem_par`` accepts a parameter `EPSG`.

    Parameters
    ----------
    geometry: spatialist.vector.Vector
        a vector geometry delimiting the output DEM size
    demType: str
        the type of DEM to be used; see :func:`~pyroSAR.auxdata.dem_autoload` for options
    outfile: str
        the name of the final DEM file
    buffer: float or None
        a buffer in degrees to create around the geometry
    t_srs: int, str or osgeo.osr.SpatialReference
        A target geographic reference system in WKT, EPSG, PROJ4 or OPENGIS format.
        See function :func:`spatialist.auxil.crsConvert()` for details.
        Default: `4326 <https://spatialreference.org/ref/epsg/4326/>`_.
    tr: tuple or None
        the target resolution as (xres, yres) in units of ``t_srs``; if ``t_srs`` is kept at its default value of 4326,
        ``tr`` does not need to be defined and the original resolution is preserved;
        in all other cases the default of None is rejected
    logpath: str
        a directory to write GAMMA logfiles to
    username: str or None
        (optional) the user name for services requiring registration;
        see :func:`~pyroSAR.auxdata.dem_autoload`
    password: str or None
        (optional) the password for the registration account
    geoid_mode: str
        the software to be used for converting geoid to ellipsoid heights (if necessary); options:
        
         - 'gamma'
         - 'gdal'
    resampling_method: str
        the gdalwarp resampling method; See `here <https://gdal.org/programs/gdalwarp.html#cmdoption-gdalwarp-r>`_
        for options.

    Returns
    -------

    """
    geometry = geometry.clone()
    
    epsg = crsConvert(t_srs, 'epsg') if t_srs != 4326 else t_srs
    
    if epsg != 4326:
        if not hasarg(diff.create_dem_par, 'EPSG'):
            raise RuntimeError('using a different CRS than 4326 is currently '
                               'not supported for this version of GAMMA')
        if 'dem_import' not in dir(diff):
            raise RuntimeError('using a different CRS than 4326 currently requires command '
                               'dem_import, which is not part of this version of GAMMA')
        if tr is None:
            raise RuntimeError('tr needs to be defined if t_srs is not 4326')
    
    if os.path.isfile(outfile):
        log.info('outfile already exists')
        return
    
    tmpdir = outfile + '__tmp'
    os.makedirs(tmpdir)
    
    try:
        if logpath is not None and not os.path.isdir(logpath):
            os.makedirs(logpath)
        
        vrt = os.path.join(tmpdir, 'dem.vrt')
        dem = os.path.join(tmpdir, 'dem.tif')
        
        if epsg == geometry.getProjection('epsg') and buffer is None:
            ext = geometry.extent
            bounds = [ext['xmin'], ext['ymin'], ext['xmax'], ext['ymax']]
        else:
            bounds = None
        geometry.reproject(4326)
        log.info('collecting DEM tiles')
        dem_autoload([geometry], demType, vrt=vrt, username=username,
                     password=password, buffer=buffer)
        
        # TanDEM-X DEM, GETASSE30 DEM: ellipsoidal heights,
        # Copernicus DEM: EGM2008 geoid, all others are EGM96 heights
        # GAMMA works only with ellipsoid heights and the offset needs to be corrected
        # starting from GDAL 2.2 the conversion can be done directly in GDAL; see docs of gdalwarp
        message = 'conversion to GAMMA format'
        geoid = None
        if demType not in ['TDX90m', 'GETASSE30']:
            message = 'geoid correction and conversion to GAMMA format'
            if re.search('Copernicus [139]0m', demType):
                geoid = 'EGM2008'
            elif demType in ['AW3D30', 'SRTM 1Sec HGT', 'SRTM 3Sec']:
                geoid = 'EGM96'
            else:
                raise RuntimeError("'demType' is not supported")
        
        if geoid_mode == 'gdal':
            gamma_geoid = None
            if geoid is not None:
                gdal_geoid = True
            else:
                gdal_geoid = False
        elif geoid_mode == 'gamma':
            gdal_geoid = False
            gamma_geoid = geoid
        else:
            raise RuntimeError("'geoid_mode' is not supported")
        
        dem_create(vrt, dem, t_srs=epsg, tr=tr, geoid_convert=gdal_geoid,
                   resampleAlg=resampling_method, outputBounds=bounds,
                   geoid=geoid)
        
        outfile_tmp = os.path.join(tmpdir, os.path.basename(outfile))
        
        log.info(message)
        
        dem_import(src=dem, dst=outfile_tmp, geoid=gamma_geoid,
                   logpath=logpath, outdir=tmpdir)
        
        for suffix in ['', '.par', '.hdr']:
            shutil.copyfile(outfile_tmp + suffix, outfile + suffix)
    
    except RuntimeError as e:
        raise e
    finally:
        shutil.rmtree(tmpdir)


def dem_import(src, dst, geoid=None, logpath=None, outdir=None):
    """
    convert an existing DEM in GDAL-readable format to GAMMA format including optional geoid-ellipsoid conversion.
    
    Parameters
    ----------
    src: str
        the input DEM
    dst: str
        the output DEM
    geoid: str or None
        the geoid height reference of `src`; supported options:
        
        - 'EGM96'
        - 'EGM2008'
        - None: assume WGS84 ellipsoid heights and do not convert heights
    logpath: str or None
        a directory to write logfiles to
    outdir: str or None
        the directory to execute the command in

    Returns
    -------

    """
    with raster.Raster(src) as ras:
        epsg = ras.epsg
    if epsg != 4326:
        if not hasarg(diff.create_dem_par, 'EPSG'):
            raise RuntimeError('using a different CRS than EPSG:4326 is currently '
                               'not supported for this version of GAMMA')
        if 'dem_import' not in dir(diff):
            raise RuntimeError('using a different CRS than 4326 currently requires command '
                               'dem_import, which is not part of this version of GAMMA')
    dst_base = os.path.splitext(dst)[0]
    if geoid is not None:
        # "Add interpolated geoid offset relative to the WGS84 datum;
        # NODATA are set to the interpolated geoid offset."
        gflg = 2
    else:
        # "No geoid offset correction, replace NODATA with a valid near-zero value."
        gflg = 0
    if epsg == 4326 and geoid == 'EGM96':
        # old approach for backwards compatibility
        diff.srtm2dem(SRTM_DEM=src,
                      DEM=dst,
                      DEM_par=dst + '.par',
                      gflg=gflg,
                      geoid='-',
                      logpath=logpath,
                      outdir=outdir)
    else:
        # new approach enabling an arbitrary target CRS EPSG code
        diff.create_dem_par(DEM_par=dst_base + '.par',
                            inlist=[''] * 9,
                            EPSG=epsg)
        dem_import_pars = {'input_DEM': src,
                           'DEM': dst,
                           'DEM_par': dst_base + '.par',
                           'logpath': logpath,
                           'outdir': outdir}
        if gflg == 2:
            home = ExamineGamma().home
            if geoid == 'EGM96':
                geoid_file = os.path.join(home, 'DIFF', 'scripts', 'egm96.dem')
            elif geoid == 'EGM2008':
                geoid_file = os.path.join(home, 'DIFF', 'scripts', 'egm2008-5.dem')
            else:
                raise RuntimeError('conversion of {} geoid is not supported by GAMMA'.format(geoid))
            dem_import_pars['geoid'] = geoid_file
            dem_import_pars['geoid_par'] = geoid_file + '_par'
        
        diff.dem_import(**dem_import_pars)
    
    par2hdr(dst_base + '.par', dst_base + '.hdr', nodata=0)


def dempar(dem, logpath=None):
    """
    create GAMMA parameter text files for DEM files

    currently only EQA and UTM projections with WGS84 ellipsoid are supported

    Parameters
    ----------
    dem: str
        the name of the DEM
    logpath: str
        a directory to write logfiles to

    Returns
    -------

    """
    rast = raster.Raster(dem)
    
    # determine data type
    dtypes = {'Int16': 'INTEGER*2', 'UInt16': 'INTEGER*2', 'Float32': 'REAL*4'}
    if rast.dtype not in dtypes:
        raise IOError('data type not supported')
    else:
        dtype = dtypes[rast.dtype]
    
    # format pixel posting and top left coordinate
    posting = str(rast.geo['yres']) + ' ' + str(rast.geo['xres'])
    latlon = str(rast.geo['ymax']) + ' ' + str(rast.geo['xmin'])
    
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
        parlist = [projection, ellipsoid, 1, zone, falsenorthing, os.path.basename(dem),
                   dtype, 0, 1, rast.cols, rast.rows, posting, latlon]
    else:
        parlist = [projection, ellipsoid, 1, os.path.basename(dem), dtype,
                   0, 1, rast.cols, rast.rows, posting, latlon]
    
    # execute GAMMA command
    diff.create_dem_par(DEM_par=os.path.splitext(dem)[0] + '.par',
                        inlist=parlist,
                        outdir=os.path.dirname(dem),
                        logpath=logpath)


def swap(data, outname):
    """
    byte swapping from small to big endian (as required by GAMMA)

    Parameters
    ----------
    data: str
        the DEM file to be swapped
    outname: str
        the name of the file to write

    Returns
    -------

    """
    with raster.Raster(data) as ras:
        dtype = ras.dtype
        ras_format = ras.format
    if ras_format != 'ENVI':
        raise IOError('only ENVI format supported')
    dtype_lookup = {'Int16': 2, 'CInt16': 2, 'Int32': 4, 'Float32': 4, 'CFloat32': 4, 'Float64': 8}
    if dtype not in dtype_lookup:
        raise IOError('data type {} not supported'.format(dtype))
    
    disp.swap_bytes(infile=data,
                    outfile=outname,
                    swap_type=dtype_lookup[dtype])
    
    with HDRobject(data + '.hdr') as header:
        header.byte_order = 1
        header.write(outname + '.hdr')


def mosaic(demlist, outname, byteorder=1, gammapar=True):
    """
    mosaicing of multiple DEMs

    Parameters
    ----------
    demlist: list[str]
        a list of DEM names to be mosaiced
    outname: str
        the name of the final mosaic file
    byteorder: {0, 1}
        the byte order of the mosaic

        - 0: small endian
        - 1: big endian

    gammapar: bool
        create a GAMMA parameter file for the mosaic?

    Returns
    -------

    """
    if len(demlist) < 2:
        raise IOError('length of demlist < 2')
    with raster.Raster(demlist[0]) as ras:
        nodata = ras.nodata
    
    par = {'format': 'ENVI',
           'srcNodata': nodata, ' dstNodata': nodata,
           'options': ['-q']}
    gdalwarp(src=demlist, dst=outname, **par)
    
    if byteorder == 1:
        swap(outname, outname + '_swap')
        for item in [outname, outname + '.hdr', outname + '.aux.xml']:
            os.remove(item)
        os.rename(outname + '_swap', outname)
        os.rename(outname + '_swap.hdr', outname + '.hdr')
    if gammapar:
        dempar(outname)


def hgt(parfiles):
    """
    concatenate hgt file names overlapping with multiple SAR scenes

    - this list is read for corner coordinates of which the next integer
      lower left latitude and longitude is computed
    - hgt files are supplied in 1 degree equiangular format named e.g.
      N16W094.hgt (with pattern [NS][0-9]{2}[EW][0-9]{3}.hgt
    - For north and east hemisphere the respective absolute latitude and longitude
      values are smaller than the lower left coordinate of the SAR image
    - west and south coordinates are negative and hence the nearest lower left
      integer absolute value is going to be larger

    Parameters
    ----------
    parfiles: list of str or pyroSAR.ID
        a list of GAMMA parameter files or pyroSAR ID objects

    Returns
    -------
    list
        the names of hgt files overlapping with the supplied parameter files/objects
    """
    
    lat = []
    lon = []
    for parfile in parfiles:
        if isinstance(parfile, ID):
            corners = parfile.getCorners()
        elif parfile.endswith('.par'):
            corners = slc_corners(parfile)
        else:
            raise RuntimeError('parfiles items must be of type pyroSAR.ID or GAMMA parfiles with suffix .par')
        lat += [int(float(corners[x]) // 1) for x in ['ymin', 'ymax']]
        lon += [int(float(corners[x]) // 1) for x in ['xmin', 'xmax']]
    
    # add missing lat/lon values (and add an extra buffer of one degree)
    lat = range(min(lat), max(lat) + 1)
    lon = range(min(lon), max(lon) + 1)
    
    # convert coordinates to string with leading zeros and hemisphere identification letter
    lat = [str(x).zfill(2 + len(str(x)) - len(str(x).strip('-'))) for x in lat]
    lat = [x.replace('-', 'S') if '-' in x else 'N' + x for x in lat]
    
    lon = [str(x).zfill(3 + len(str(x)) - len(str(x).strip('-'))) for x in lon]
    lon = [x.replace('-', 'W') if '-' in x else 'E' + x for x in lon]
    
    # concatenate all formatted latitudes and longitudes with each other as final product
    return [x + y + '.hgt' for x in lat for y in lon]


def makeSRTM(scenes, srtmdir, outname):
    """
    Create a DEM in GAMMA format from SRTM tiles

    - coordinates are read to determine the required DEM extent and select the necessary hgt tiles
    - mosaics SRTM DEM tiles, converts them to GAMMA format and subtracts offset to WGS84 ellipsoid

    intended for SRTM products downloaded from:

    - USGS: https://gdex.cr.usgs.gov/gdex/
    - CGIAR: https://srtm.csi.cgiar.org

    Parameters
    ----------
    scenes: list of str or pyroSAR.ID
        a list of Gamma parameter files or pyroSAR ID objects to read the DEM extent from
    srtmdir: str
        a directory containing the SRTM hgt tiles
    outname: str
        the name of the final DEM file

    Returns
    -------

    """
    
    tempdir = outname + '___temp'
    os.makedirs(tempdir)
    
    hgt_options = hgt(scenes)
    
    hgt_files = finder(srtmdir, hgt_options)
    
    nodatas = list(set([raster.Raster(x).nodata for x in hgt_files]))
    if len(nodatas) == 1:
        nodata = nodatas[0]
    else:
        raise RuntimeError('different nodata values are not permitted')
    
    srtm_vrt = os.path.join(tempdir, 'srtm.vrt')
    srtm_temp = srtm_vrt.replace('.vrt', '_tmp')
    srtm_final = srtm_vrt.replace('.vrt', '')
    
    gdalbuildvrt(src=hgt_files, dst=srtm_vrt, srcNodata=nodata, options=['-overwrite'])
    
    gdal_translate(src=srtm_vrt, dst=srtm_temp, format='ENVI', noData=nodata)
    
    diff.srtm2dem(SRTM_DEM=srtm_temp,
                  DEM=srtm_final,
                  DEM_par=srtm_final + '.par',
                  gflg=2,
                  geoid='-',
                  outdir=tempdir)
    
    shutil.move(srtm_final, outname)
    shutil.move(srtm_final + '.par', outname + '.par')
    par2hdr(outname + '.par', outname + '.hdr')
    
    shutil.rmtree(tempdir)


def hgt_collect(parfiles, outdir, demdir=None, arcsec=3):
    """
    automatic downloading and unpacking of srtm tiles

    Parameters
    ----------
    parfiles: list of str or pyroSAR.ID
        a list of Gamma parameter files or pyroSAR ID objects
    outdir: str
        a target directory to download the tiles to
    demdir: str or None
        an additional directory already containing hgt tiles
    arcsec: {1, 3}
        the spatial resolution to be used

    Returns
    -------
    list
        the names of all local hgt tiles overlapping with the parfiles
    """
    
    # concatenate required hgt tile names
    target_ids = hgt(parfiles)
    
    targets = []
    
    pattern = '[NS][0-9]{2}[EW][0-9]{3}'
    
    # if an additional dem directory has been defined, check this directory for required hgt tiles
    if demdir is not None:
        targets.extend(finder(demdir, target_ids))
    
    # check for additional potentially existing hgt tiles in the defined output directory
    extras = [os.path.join(outdir, x) for x in target_ids if
              os.path.isfile(os.path.join(outdir, x)) and not re.search(x, '\n'.join(targets))]
    targets.extend(extras)
    
    log.info('found {} relevant SRTM tiles...'.format(len(targets)))
    
    # search server for all required tiles, which were not found in the local directories
    if len(targets) < len(target_ids):
        log.info('searching for additional SRTM tiles on the server...')
        onlines = []
        
        if arcsec == 1:
            remotes = ['http://e4ftl01.cr.usgs.gov/SRTM/SRTMGL1.003/2000.02.11/']
            remotepattern = pattern + '.SRTMGL1.hgt.zip'
        elif arcsec == 3:
            server = 'https://dds.cr.usgs.gov/srtm/version2_1/SRTM3/'
            remotes = [os.path.join(server, x) for x in
                       ['Africa', 'Australia', 'Eurasia', 'Islands', 'North_America', 'South_America']]
            remotepattern = pattern + '[.]hgt.zip'
        else:
            raise ValueError('argument arcsec must be of value 1 or 3')
        
        for remote in remotes:
            response = urlopen(remote).read()
            items = sorted(set(re.findall(remotepattern, response)))
            for item in items:
                outname = re.findall(pattern, item)[0] + '.hgt'
                if outname in target_ids and outname not in [os.path.basename(x) for x in targets]:
                    onlines.append(os.path.join(remote, item))
        
        # if additional tiles have been found online, download and unzip them to the local directory
        if len(onlines) > 0:
            log.info('downloading {} SRTM tiles...'.format(len(onlines)))
            for candidate in onlines:
                localname = os.path.join(outdir, re.findall(pattern, candidate)[0] + '.hgt')
                infile = urlopen(candidate)
                with open(localname + '.zip', 'wb') as outfile:
                    outfile.write(infile.read())
                infile.close()
                with zf.ZipFile(localname + '.zip', 'r') as z:
                    z.extractall(outdir)
                os.remove(localname + '.zip')
                targets.append(localname)
    return targets
