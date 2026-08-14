###############################################################################
# preparation of DEM data for use in GAMMA

# Copyright (c) 2014-2026, the pyroSAR Developers.

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
from __future__ import annotations

import os
import re
import shutil

from spatialist import raster
from spatialist.auxil import gdalwarp, crsConvert
from spatialist.envi import HDRobject

from ..auxdata import dem_autoload, dem_create
from . import ISPPar, par2hdr
from pyroSAR.examine import ExamineGamma
from pyroSAR.ancillary import hasarg
from typing import Literal, TYPE_CHECKING

if TYPE_CHECKING:
    from spatialist.vector import Vector
    from osgeo.osr import SpatialReference

import logging

log = logging.getLogger(__name__)

try:
    from .api import diff, disp, isp
except ImportError:
    pass


def fill(dem: str, dem_out: str, logpath: str | None = None, replace: bool = False) -> None:
    """
    interpolate missing values in the SRTM DEM (value -32768)

    Parameters
    ----------
    dem
        the input DEM to be filled
    dem_out
        the name of the filled DEM
    logpath
        a directory to write logfiles to
    replace
        delete `dem` once finished?
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


def dem_autocreate(
        geometry: Vector,
        demType: str,
        outfile: str,
        buffer: int | float | None = None,
        t_srs: int | str | SpatialReference = 4326,
        tr: tuple[int | float, int | float] | None = None,
        logpath: str | None = None,
        username: str | None = None,
        password: str | None = None,
        geoid_mode: Literal['gamma', 'gdal'] = 'gamma',
        resampling_method: str = 'bilinear'
) -> None:
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
    geometry
        a vector geometry delimiting the output DEM size
    demType
        the type of DEM to be used; see :func:`~pyroSAR.auxdata.dem_autoload` for options
    outfile
        the name of the final DEM file
    buffer
        a buffer in degrees to create around the geometry
    t_srs
        A target geographic reference system in WKT, EPSG, PROJ4 or OPENGIS format.
        See function :func:`spatialist.auxil.crsConvert()` for details.
        Default: `4326 <https://spatialreference.org/ref/epsg/4326/>`_.
    tr
        the target resolution as (xres, yres) in units of ``t_srs``; if ``t_srs`` is kept at its default value of 4326,
        ``tr`` does not need to be defined and the original resolution is preserved;
        in all other cases the default of None is rejected
    logpath
        a directory to write GAMMA logfiles to
    username
        (optional) the user name for services requiring registration;
        see :func:`~pyroSAR.auxdata.dem_autoload`
    password
        (optional) the password for the registration account
    geoid_mode
        the software to be used for converting geoid to ellipsoid heights (if necessary)
    resampling_method
        the gdalwarp resampling method; See `here <https://gdal.org/programs/gdalwarp.html#cmdoption-gdalwarp-r>`_
        for options.
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
        dem_autoload(geometry=geometry, demType=demType,
                     vrt=vrt, username=username,
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
        
        dem_create(geometry=geometry,
                   src=vrt, dst=dem, t_srs=epsg, tr=tr, geoid_convert=gdal_geoid,
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


def dem_import(
        src: str,
        dst: str,
        geoid: Literal['EGM96', 'EGM2008'] | None = None,
        logpath: str | None = None,
        outdir: str | None = None,
        shellscript: str | None = None
) -> None:
    """
    convert an existing DEM in GDAL-readable format to GAMMA
    format including optional geoid-ellipsoid conversion.
    
    Parameters
    ----------
    src
        the input DEM
    dst
        the output DEM
    geoid
        the geoid height reference of `src`; supported options.
        Default `None`: assume WGS84 ellipsoid heights and do not convert heights.
    logpath
        a directory to write logfiles to
    outdir
        the directory to execute the command in
    shellscript
        a file to write the GAMMA commands to in shell format
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
                      outdir=outdir,
                      shellscript=shellscript)
    else:
        # new approach enabling an arbitrary target CRS EPSG code
        diff.create_dem_par(DEM_par=dst_base + '.par',
                            inlist=[''] * 9,
                            EPSG=epsg,
                            logpath=logpath,
                            outdir=outdir,
                            shellscript=shellscript)
        dem_import_pars = {'input_DEM': src,
                           'DEM': dst,
                           'DEM_par': dst_base + '.par',
                           'logpath': logpath,
                           'outdir': outdir,
                           'shellscript': shellscript}
        if gflg == 2:
            home = ExamineGamma().home
            if geoid == 'EGM96':
                geoid_file = os.path.join(home, 'DIFF', 'scripts', 'egm96.dem')
            elif geoid == 'EGM2008':
                geoid_file = os.path.join(home, 'DIFF', 'scripts', 'egm2008-5.dem')
            else:
                raise RuntimeError(f"conversion of '{geoid}' geoid is not supported by GAMMA")
            dem_import_pars['geoid'] = geoid_file
            dem_import_pars['geoid_par'] = geoid_file + '_par'
        
        diff.dem_import(**dem_import_pars)
    
    par2hdr(dst_base + '.par', dst_base + '.hdr', nodata=0)


def dempar(dem: str, logpath: str | None = None) -> None:
    """
    create GAMMA parameter text files for DEM files

    currently only EQA and UTM projections with WGS84 ellipsoid are supported

    Parameters
    ----------
    dem
        the name of the DEM
    logpath
        a directory to write logfiles to
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


def swap(data: str, outname: str) -> None:
    """
    byte swapping from small to big endian (as required by GAMMA)

    Parameters
    ----------
    data
        the DEM file to be swapped
    outname
        the name of the file to write
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


def mosaic(
        demlist: list[str],
        outname: str,
        byteorder: Literal[0, 1] = 1,
        gammapar: bool = True
) -> None:
    """
    mosaicing of multiple DEMs

    Parameters
    ----------
    demlist
        a list of DEM names to be mosaiced
    outname
        the name of the final mosaic file
    byteorder
        the byte order of the mosaic

        - 0: small endian
        - 1: big endian

    gammapar
        create a GAMMA parameter file for the mosaic?
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
