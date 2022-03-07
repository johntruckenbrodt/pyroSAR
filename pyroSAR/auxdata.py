###############################################################################
# tools for handling auxiliary data in software pyroSAR

# Copyright (c) 2019-2022, the pyroSAR Developers.

# This file is part of the pyroSAR Project. It is subject to the
# license terms in the LICENSE.txt file found in the top-level
# directory of this distribution and at
# https://github.com/johntruckenbrodt/pyroSAR/blob/master/LICENSE.txt.
# No part of the pyroSAR project, including this file, may be
# copied, modified, propagated, or distributed except according
# to the terms contained in the LICENSE.txt file.
###############################################################################
import io
import os
import re
import csv
import ssl
import ftplib
import requests
import zipfile as zf
from math import ceil, floor
from urllib.request import urlopen
from urllib.error import HTTPError
from urllib.parse import urlparse

from pyroSAR.examine import ExamineSnap
from spatialist.raster import Raster, Dtype
from spatialist.ancillary import dissolve, finder
from spatialist.auxil import gdalbuildvrt, crsConvert, gdalwarp
from spatialist.envi import HDRobject
from osgeo import gdal

import logging

log = logging.getLogger(__name__)


def dem_autoload(geometries, demType, vrt=None, buffer=None, username=None, password=None,
                 product='dem', nodata=None, hide_nodata=False):
    """
    obtain all relevant DEM tiles for selected geometries

    Parameters
    ----------
    geometries: list[spatialist.vector.Vector]
        a list of :class:`spatialist.vector.Vector` geometries to obtain DEM data for;
        CRS must be WGS84 LatLon (EPSG 4326)
    demType: str
        the type of DEM to be used; current options:

        - 'AW3D30' (ALOS Global Digital Surface Model "ALOS World 3D - 30m")

          * info: https://www.eorc.jaxa.jp/ALOS/en/aw3d30/index.htm
          * url: ftp://ftp.eorc.jaxa.jp/pub/ALOS/ext1/AW3D30/release_v1804
          * height reference: EGM96

        - 'Copernicus 10m EEA DEM' (Copernicus 10 m DEM available over EEA-39 countries)

          * registration: https://spacedata.copernicus.eu/web/cscda/data-access/registration
          * url: ftps://cdsdata.copernicus.eu/DEM-datasets/COP-DEM_EEA-10-DGED/2021_1
          * height reference: EGM2008

        - 'Copernicus 30m Global DEM'
          
          * info: https://copernicus-dem-30m.s3.amazonaws.com/readme.html
          * url: https://copernicus-dem-30m.s3.eu-central-1.amazonaws.com/
          * height reference: EGM2008

        - 'Copernicus 30m Global DEM II'
        
          * registration: https://spacedata.copernicus.eu/web/cscda/data-access/registration
          * url: ftps://cdsdata.copernicus.eu/DEM-datasets/COP-DEM_GLO-30-DGED/2021_1
          * height reference: EGM2008
        
        - 'Copernicus 90m Global DEM'
     
          * info: https://copernicus-dem-90m.s3.amazonaws.com/readme.html
          * url: https://copernicus-dem-90m.s3.eu-central-1.amazonaws.com/
          * height reference: EGM2008
        
        - 'Copernicus 90m Global DEM II'
        
          * registration: https://spacedata.copernicus.eu/web/cscda/data-access/registration
          * url: ftps://cdsdata.copernicus.eu/DEM-datasets/COP-DEM_GLO-90-DGED/2021_1
          * height reference: EGM2008
        
        - 'GETASSE30'
        
          * info: https://seadas.gsfc.nasa.gov/help-8.1.0/desktop/GETASSE30ElevationModel.html
          * url: https://step.esa.int/auxdata/dem/GETASSE30
          * height reference: WGS84
        
        - 'SRTM 1Sec HGT'

          * url: https://step.esa.int/auxdata/dem/SRTMGL1
          * height reference: EGM96

        - 'SRTM 3Sec'

          * url: https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF
          * height reference: EGM96
        
        - 'TDX90m'
        
          * registration:  https://geoservice.dlr.de/web/dataguide/tdm90
          * url: ftpes://tandemx-90m.dlr.de
          * height reference: WGS84

    vrt: str or None
        an optional GDAL VRT file created from the obtained DEM tiles
    buffer: int, float, None
        a buffer in degrees to add around the individual geometries
    username: str or None
        (optional) the user name for services requiring registration
    password: str or None
        (optional) the password for the registration account
    product: str
        the sub-product to extract from the DEM product.
        The following options are available for the respective DEM types:
        
        - 'AW3D30'
        
          * 'dem': the actual Digital Elevation Model
          * 'msk': mask information for each pixel (Cloud/Snow Mask, Land water and
            low correlation mask, Sea mask, Information of elevation dataset used
            for the void-filling processing)
          * 'stk': number of DSM-scene files which were used to produce the 5 m resolution DSM
        
        - 'Copernicus 10m EEA DEM'
        
          * 'dem': the actual Digital Elevation Model
          * 'edm': editing mask
          * 'flm': filling mask
          * 'hem': height error mask
          * 'wbm': water body mask
        
        - 'Copernicus 30m Global DEM'
        
          * 'dem': the actual Digital Elevation Model
        
        - 'Copernicus 30m Global DEM II'
        
          * 'dem': the actual Digital Elevation Model
          * 'edm': editing mask
          * 'flm': filling mask
          * 'hem': height error mask
          * 'wbm': water body mask
        
        - 'Copernicus 90m Global DEM'
        
          * 'dem': the actual Digital Elevation Model
        
        - 'Copernicus 90m Global DEM II'
        
          * 'dem': the actual Digital Elevation Model
          * 'edm': editing mask
          * 'flm': filling mask
          * 'hem': height error mask
          * 'wbm': water body mask
        
        - 'GETASSE30'
        
          * 'dem': the actual Digital Elevation Model
        
        - 'SRTM 1Sec HGT'
        
          * 'dem': the actual Digital Elevation Model
        
        - 'SRTM 3Sec'
        
          * 'dem': the actual Digital Elevation Model
        
        - 'TDX90m'
        
          * 'dem': the actual Digital Elevation Model
          * 'am2': Amplitude Mosaic representing the minimum value
          * 'amp': Amplitude Mosaic representing the mean value
          * 'com': Consistency Mask
          * 'cov': Coverage Map
          * 'hem': Height Error Map
          * 'lsm': Layover and Shadow Mask, based on SRTM C-band and Globe DEM data
          * 'wam': Water Indication Mask
    
    Returns
    -------
    list or None
        the names of the obtained files or None if a VRT file was defined
    
    Examples
    --------
    download all SRTM 1 arcsec DEMs overlapping with a Sentinel-1 scene and mosaic them to a single GeoTIFF file
    
    .. code-block:: python
        
        from pyroSAR import identify
        from pyroSAR.auxdata import dem_autoload
        from spatialist import gdalwarp
        
        # identify the SAR scene
        filename = 'S1A_IW_SLC__1SDV_20150330T170734_20150330T170801_005264_006A6C_DA69.zip'
        scene = identify(filename)
        
        # extract the bounding box as spatialist.Vector object
        bbox = scene.bbox()
        
        # download the tiles and virtually combine them in an in-memory
        # VRT file subsetted to the extent of the SAR scene plus a buffer of 0.01 degrees
        vrt = '/vsimem/srtm1.vrt'
        dem_autoload(geometries=[bbox], demType='SRTM 1Sec HGT',
                     vrt=vrt, buffer=0.01)
        
        # write the final GeoTIFF file
        outname = scene.outname_base() + 'srtm1.tif'
        gdalwarp(src=vrt, dst=outname, options={'format': 'GTiff'})
        
        # alternatively use function dem_create and warp the DEM to UTM
        # including conversion from geoid to ellipsoid heights
        from pyroSAR.auxdata import dem_create
        outname = scene.outname_base() + 'srtm1_ellp.tif'
        dem_create(src=vrt, dst=outname, t_srs=32632, tr=(30, 30),
                   geoid_convert=True, geoid='EGM96')
    """
    with DEMHandler(geometries) as handler:
        return handler.load(demType=demType,
                            username=username,
                            password=password,
                            vrt=vrt,
                            buffer=buffer,
                            product=product,
                            nodata=nodata,
                            hide_nodata=hide_nodata)


def dem_create(src, dst, t_srs=None, tr=None, resampling_method='bilinear', threads=None,
               geoid_convert=False, geoid='EGM96', outputBounds=None, dtype=None, pbar=False):
    """
    create a new DEM GeoTIFF file and optionally convert heights from geoid to ellipsoid
    
    Parameters
    ----------
    src: str
        the input dataset, e.g. a VRT from function :func:`dem_autoload`
    dst: str
        the output dataset
    t_srs: None, int, str or osr.SpatialReference
        A target geographic reference system in WKT, EPSG, PROJ4 or OPENGIS format.
        See function :func:`spatialist.auxil.crsConvert()` for details.
        Default (None): use the crs of ``src``.
    tr: None or tuple
        the target resolution as (xres, yres)
    resampling_method: str
        the gdalwarp resampling method; See `here <https://gdal.org/programs/gdalwarp.html#cmdoption-gdalwarp-r>`_
        for options.
    threads: int, str or None
        the number of threads to use. Possible values:
        
         - Default `None`: use the value of `GDAL_NUM_THREADS` without modification. If `GDAL_NUM_THREADS` is None,
           multi-threading is still turned on and two threads are used, one for I/O and one for computation.
         - integer value: temporarily modify `GDAL_NUM_THREADS` and reset it once done.
           If 1, multithreading is turned off.
         - `ALL_CPUS`: special string to use all cores/CPUs of the computer; will also temporarily
           modify `GDAL_NUM_THREADS`.
    geoid_convert: bool
        convert geoid heights?
    geoid: str
        the geoid model to be corrected, only used if ``geoid_convert == True``; current options:
        
         - 'EGM96'
         - 'EGM2008'
    outputBounds: list or None
        output bounds as [xmin, ymin, xmax, ymax] in target SRS
    dtype: str or None
        override the data type of the written file; Default None: use same type as source data.
        Data type notations of GDAL (e.g. `Float32`) and numpy (e.g. `int8`) are supported.
    pbar: bool
        add a progressbar?
    
    Returns
    -------

    """
    
    with Raster(src) as ras:
        nodata = ras.nodata
        epsg_in = ras.epsg
    
    if t_srs is None:
        epsg_out = epsg_in
    else:
        epsg_out = crsConvert(t_srs, 'epsg')
    
    threads_system = gdal.GetConfigOption('GDAL_NUM_THREADS')
    if threads is None:
        threads = threads_system
        try:
            threads = int(threads)
        except (ValueError, TypeError):
            pass
    if isinstance(threads, str):
        if threads != 'ALL_CPUS':
            raise ValueError("unsupported value for 'threads': '{}'".format(threads))
        else:
            multithread = True
            gdal.SetConfigOption('GDAL_NUM_THREADS', threads)
    elif isinstance(threads, int):
        if threads == 1:
            multithread = False
        elif threads > 1:
            multithread = True
            gdal.SetConfigOption('GDAL_NUM_THREADS', str(threads))
        else:
            raise ValueError("if 'threads' is of type int, it must be >= 1")
    elif threads is None:
        multithread = True
    else:
        raise TypeError("'threads' must be of type int, str or None. Is: {}".format(type(threads)))
    
    gdalwarp_args = {'format': 'GTiff', 'multithread': multithread,
                     'srcNodata': nodata, 'dstNodata': nodata,
                     'srcSRS': 'EPSG:{}'.format(epsg_in),
                     'dstSRS': 'EPSG:{}'.format(epsg_out),
                     'resampleAlg': resampling_method}
    
    if dtype is not None:
        gdalwarp_args['outputType'] = Dtype(dtype).gdalint
    
    if outputBounds is not None:
        gdalwarp_args['outputBounds'] = outputBounds
    
    if tr is not None:
        gdalwarp_args.update({'xRes': tr[0],
                              'yRes': tr[1],
                              'targetAlignedPixels': True})
    
    if geoid_convert:
        geoid_epsg = {'EGM96': 5773,
                      'EGM2008': 3855}
        if geoid in geoid_epsg.keys():
            epsg = geoid_epsg[geoid]
            gdalwarp_args['srcSRS'] += '+{}'.format(epsg)
            # the following line is a temporary workaround until compound EPSG codes can
            # directly be used for vertical CRS transformations
            # see https://github.com/OSGeo/gdal/pull/4639
            gdalwarp_args['srcSRS'] = crsConvert(gdalwarp_args['srcSRS'], 'proj4')
        else:
            raise RuntimeError('geoid model not yet supported')
        try:
            get_egm_lookup(geoid=geoid, software='PROJ')
        except OSError as e:
            errstr = str(e)
            addition = '\nplease refer to the following site for instructions ' \
                       'on how to use EGM GTX files (requires PROJ >= 5.0.0):\n' \
                       'https://gis.stackexchange.com/questions/258532/' \
                       'noaa-vdatum-gdal-variable-paths-for-linux-ubuntu'
            raise RuntimeError(errstr + addition)
    
    try:
        message = 'creating mosaic'
        crs = gdalwarp_args['dstSRS']
        if crs != 'EPSG:4326':
            message += ' and reprojecting to {}'.format(crs)
        message += ': {}'.format(dst)
        log.info(message)
        gdalwarp(src, dst, gdalwarp_args, pbar)
    except Exception:
        if os.path.isfile(dst):
            os.remove(dst)
        raise
    finally:
        gdal.SetConfigOption('GDAL_NUM_THREADS', threads_system)


class DEMHandler:
    """
    | An interface to obtain DEM data for selected geometries
    | The files are downloaded into the ESA SNAP auxdata directory structure
    
    Parameters
    ----------
    geometries: list of spatialist.vector.Vector
        a list of geometries
    """
    
    def __init__(self, geometries):
        if not isinstance(geometries, list):
            raise RuntimeError('geometries must be of type list')
        
        for geometry in geometries:
            if geometry.getProjection('epsg') != 4326:
                raise RuntimeError('input geometry CRS must be WGS84 LatLon (EPSG 4326)')
        
        self.geometries = geometries
        try:
            self.auxdatapath = ExamineSnap().auxdatapath
        except AttributeError:
            self.auxdatapath = os.path.join(os.path.expanduser('~'), '.snap', 'auxdata')
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        return
    
    @staticmethod
    def __applybuffer(extent, buffer):
        ext = dict(extent)
        if buffer is not None:
            ext['xmin'] -= buffer
            ext['xmax'] += buffer
            ext['ymin'] -= buffer
            ext['ymax'] += buffer
        return ext
    
    @staticmethod
    def __buildvrt(tiles, vrtfile, pattern, vsi, extent, nodata=None, hide_nodata=False):
        if vsi is not None:
            locals = [vsi + x for x in dissolve([finder(x, [pattern]) for x in tiles])]
        else:
            locals = tiles
        with Raster(locals[0]) as ras:
            if nodata is None:
                nodata = ras.nodata
            xres, yres = ras.res
        opts = {'outputBounds': (extent['xmin'], extent['ymin'],
                                 extent['xmax'], extent['ymax']),
                'srcNodata': nodata, 'targetAlignedPixels': True,
                'xRes': xres, 'yRes': yres, 'hideNodata': hide_nodata
                }
        gdalbuildvrt(src=locals, dst=vrtfile,
                     options=opts)
    
    def __commonextent(self, buffer=None):
        ext_new = {}
        for geo in self.geometries:
            if len(ext_new.keys()) == 0:
                ext_new = geo.extent
            else:
                for key in ['xmin', 'ymin']:
                    if geo.extent[key] > ext_new[key]:
                        ext_new[key] = geo.extent[key]
                for key in ['xmax', 'ymax']:
                    if geo.extent[key] < ext_new[key]:
                        ext_new[key] = geo.extent[key]
        ext_new = self.__applybuffer(ext_new, buffer)
        return ext_new
    
    @staticmethod
    def intrange(extent, step):
        """
        generate sequence of integer coordinates marking the tie points of the individual DEM tiles
        
        Parameters
        ----------
        extent: dict
            a dictionary with keys `xmin`, `xmax`, `ymin` and `ymax` with coordinates in EPSG:4326.
        step: int
            the sequence steps

        Returns
        -------
        tuple[range]
            the integer sequences as (latitude, longitude)
        """
        lat = range(floor(float(extent['ymin']) / step) * step,
                    ceil(float(extent['ymax']) / step) * step,
                    step)
        lon = range(floor(float(extent['xmin']) / step) * step,
                    ceil(float(extent['xmax']) / step) * step,
                    step)
        return lat, lon
    
    @staticmethod
    def __retrieve(url, filenames, outdir):
        files = list(set(filenames))
        os.makedirs(outdir, exist_ok=True)
        locals = []
        for file in files:
            infile = '{}/{}'.format(url, file)
            outfile = os.path.join(outdir, os.path.basename(file))
            if not os.path.isfile(outfile):
                try:
                    input = urlopen(infile)
                    log.info('{} <<-- {}'.format(outfile, infile))
                except HTTPError:
                    continue
                with open(outfile, 'wb') as output:
                    output.write(input.read())
                input.close()
            if os.path.isfile(outfile):
                locals.append(outfile)
        return sorted(locals)
    
    def __retrieve_ftp(self, url, filenames, outdir, username, password, port=0):
        files = list(set(filenames))
        os.makedirs(outdir, exist_ok=True)
        
        parsed = urlparse(url)
        
        if parsed.scheme == 'ftpes':
            ftp = ftplib.FTP_TLS(parsed.netloc)
            try:
                ftp.login(username, password)  # login anonymously before securing control channel
            except ftplib.error_perm as e:
                raise RuntimeError(str(e))
            ftp.prot_p()  # switch to secure data connection.. IMPORTANT! Otherwise, only the user and password is encrypted and not all the file data.
        elif parsed.scheme == 'ftps':
            ftp = ImplicitFTP_TLS()
            ftp.connect(host=parsed.netloc, port=port)
            ftp.login(username, password)
        else:
            ftp = ftplib.FTP(parsed.netloc, timeout=100)
            ftp.login()
        if parsed.path != '':
            ftp.cwd(parsed.path)
        locals = []
        for product_remote in files:
            product_local = os.path.join(outdir, os.path.basename(product_remote))
            if not os.path.isfile(product_local):
                try:
                    targetlist = ftp.nlst(product_remote)
                except ftplib.error_temp:
                    continue
                address = '{}://{}/{}{}'.format(parsed.scheme, parsed.netloc,
                                                parsed.path + '/' if parsed.path != '' else '', product_remote)
                log.info('{} <<-- {}'.format(product_local, address))
                with open(product_local, 'wb') as myfile:
                    ftp.retrbinary('RETR {}'.format(product_remote), myfile.write)
            if os.path.isfile(product_local):
                locals.append(product_local)
        ftp.close()
        return sorted(locals)
    
    @property
    def config(self):
        return {
            'AW3D30': {'url': 'ftp://ftp.eorc.jaxa.jp/pub/ALOS/ext1/AW3D30/release_v1804',
                       'nodata': -9999,
                       'vsi': '/vsitar/',
                       'pattern': {'dem': '*DSM.tif',
                                   'msk': '*MSK.tif',
                                   'stk': '*STK.tif'}
                       },
            'Copernicus 10m EEA DEM': {'url': 'ftps://cdsdata.copernicus.eu/DEM-datasets/COP-DEM_EEA-10-DGED/2021_1',
                                       'nodata': -32767.0,
                                       'vsi': '/vsitar/',
                                       'port': 990,
                                       'pattern': {'dem': '*DEM.tif',
                                                   'edm': '*EDM.tif',
                                                   'flm': '*FLM.tif',
                                                   'hem': '*HEM.tif',
                                                   'wbm': '*WBM.tif'}
                                       },
            'Copernicus 30m Global DEM': {'url': 'https://copernicus-dem-30m.s3.eu-central-1.amazonaws.com',
                                          'nodata': None,
                                          'vsi': None,
                                          'pattern': {'dem': '*DSM*'}
                                          },
            'Copernicus 30m Global DEM II': {
                'url': 'ftps://cdsdata.copernicus.eu/DEM-datasets/COP-DEM_GLO-30-DGED/2021_1',
                'nodata': -32767.0,
                'vsi': '/vsitar/',
                'port': 990,
                'pattern': {'dem': '*DEM.tif',
                            'edm': '*EDM.tif',
                            'flm': '*FLM.tif',
                            'hem': '*HEM.tif',
                            'wbm': '*WBM.tif'}
            },
            'Copernicus 90m Global DEM': {'url': 'https://copernicus-dem-90m.s3.eu-central-1.amazonaws.com',
                                          'nodata': None,
                                          'vsi': None,
                                          'pattern': {'dem': '*DSM*'}
                                          },
            'Copernicus 90m Global DEM II': {
                'url': 'ftps://cdsdata.copernicus.eu/DEM-datasets/COP-DEM_GLO-90-DGED/2021_1',
                'nodata': -32767.0,
                'vsi': '/vsitar/',
                'port': 990,
                'pattern': {'dem': '*DEM.tif',
                            'edm': '*EDM.tif',
                            'flm': '*FLM.tif',
                            'hem': '*HEM.tif',
                            'wbm': '*WBM.tif'}
            },
            'GETASSE30': {'url': 'https://step.esa.int/auxdata/dem/GETASSE30',
                          'nodata': None,
                          'vsi': '/vsizip/',
                          'pattern': {'dem': '*.GETASSE30'}
                          },
            'SRTM 1Sec HGT': {'url': 'https://step.esa.int/auxdata/dem/SRTMGL1',
                              'nodata': -32768.0,
                              'vsi': '/vsizip/',
                              'pattern': {'dem': '*.hgt'}
                              },
            'SRTM 3Sec': {'url': 'https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF',
                          'nodata': -32768.0,
                          'vsi': '/vsizip/',
                          'pattern': {'dem': 'srtm*.tif'}
                          },
            'TDX90m': {'url': 'ftpes://tandemx-90m.dlr.de',
                       'nodata': -32767.0,
                       'vsi': '/vsizip/',
                       'pattern': {'dem': '*_DEM.tif',
                                   'am2': '*_AM2.tif',
                                   'amp': '*_AMP.tif',
                                   'com': '*_COM.tif',
                                   'cov': '*_COV.tif',
                                   'hem': '*_HEM.tif',
                                   'lsm': '*_LSM.tif',
                                   'wam': '*_WAM.tif'}
                       }
        }
    
    def load(self, demType, vrt=None, buffer=None, username=None, password=None,
             product='dem', nodata=None, hide_nodata=False):
        """
        obtain DEM tiles for the given geometries
        
        Parameters
        ----------
        demType: str
            the type fo DEM to be used
        vrt: str or None
            an optional GDAL VRT file created from the obtained DEM tiles
        buffer: int, float, None
            a buffer in degrees to add around the individual geometries
        username: str or None
            the download account username
        password: str or None
            the download account password
        product: str
            the sub-product to extract from the DEM product
             - 'AW3D30'
             
              * 'dem': the actual Digital Elevation Model
              * 'msk': mask information for each pixel (Cloud/Snow Mask, Land water and
                low correlation mask, Sea mask, Information of elevation dataset used
                for the void-filling processing)
              * 'stk': number of DSM-scene files which were used to produce the 5m resolution DSM
             
             - 'Copernicus 10m EEA DEM'
             
              * 'dem': the actual Digital Elevation Model
              * 'edm': Editing Mask
              * 'flm': Filling Mask
              * 'hem': Height Error Mask
              * 'wbm': Water Body Mask
             
             - 'Copernicus 30m Global DEM'
             
              * 'dem': the actual Digital Elevation Model
             
             - 'Copernicus 30m Global DEM II'
             
              * 'dem': the actual Digital Elevation Model
              * 'edm': Editing Mask
              * 'flm': Filling Mask
              * 'hem': Height Error Mask
              * 'wbm': Water Body Mask
             
             - 'Copernicus 90m Global DEM'
             
              * 'dem': the actual Digital Elevation Model
             
             - 'Copernicus 90m Global DEM II'
             
              * 'dem': the actual Digital Elevation Model
              * 'edm': Editing Mask
              * 'flm': Filling Mask
              * 'hem': Height Error Mask
              * 'wbm': Water Body Mask
             
             - 'GETASSE30'
             
              * 'dem': the actual Digital Elevation Model
             
             - 'SRTM 1Sec HGT'
             
              * 'dem': the actual Digital Elevation Model
             
             - 'SRTM 3Sec'
             
              * 'dem': the actual Digital Elevation Model
             
             - 'TDX90m'
             
              * 'dem': the actual Digital Elevation Model
              * 'am2': Amplitude Mosaic representing the minimum value
              * 'amp': Amplitude Mosaic representing the mean value
              * 'com': Consistency Mask
              * 'cov': Coverage Map
              * 'hem': Height Error Map
              * 'lsm': Layover and Shadow Mask, based on SRTM C-band and Globe DEM data
              * 'wam': Water Indication Mask
        nodata: int or float or None
            the no data value to write in the VRT if it will be written.
            If `None`, the value of the source products is passed on.

        Returns
        -------
        list or None
            the names of the obtained files or None if a VRT file was defined
        """
        keys = self.config.keys()
        if demType not in keys:
            raise RuntimeError("demType '{}' is not supported\n  "
                               "possible options: '{}'"
                               .format(demType, "', '".join(keys)))
        
        products = self.config[demType]['pattern'].keys()
        if product not in products:
            raise RuntimeError("product '{0}' not available for demType '{1}'\n"
                               "  options: '{2}'".format(product, demType, "', '".join(products)))
        
        outdir = os.path.join(self.auxdatapath, 'dem', demType)
        
        remotes = []
        for geo in self.geometries:
            corners = self.__applybuffer(geo.extent, buffer)
            remotes.extend(self.remote_ids(corners, demType=demType,
                                           username=username, password=password))
        
        if demType in ['AW3D30', 'TDX90m', 'Copernicus 10m EEA DEM',
                       'Copernicus 30m Global DEM II', 'Copernicus 90m Global DEM II']:
            port = 0
            if 'port' in self.config[demType].keys():
                port = self.config[demType]['port']
            locals = self.__retrieve_ftp(self.config[demType]['url'], remotes, outdir,
                                         username=username, password=password, port=port)
        else:
            locals = self.__retrieve(self.config[demType]['url'], remotes, outdir)
        
        if demType == 'GETASSE30':
            for item in locals:
                getasse30_hdr(item)
        
        if nodata is None:
            if product == 'dem':
                nodata = self.config[demType]['nodata']
        
        if vrt is not None:
            self.__buildvrt(tiles=locals, vrtfile=vrt,
                            pattern=self.config[demType]['pattern'][product],
                            vsi=self.config[demType]['vsi'],
                            extent=self.__commonextent(buffer),
                            nodata=nodata, hide_nodata=hide_nodata)
            return None
        return locals
    
    def remote_ids(self, extent, demType, username=None, password=None):
        """
        parse the names of the remote files overlapping with an area of interest

        Parameters
        ----------
        extent: dict
            the extent of the area of interest with keys xmin, xmax, ymin, ymax
        demType: str
            the type fo DEM to be used
        username: str or None
            the download account user name
        password: str or None
            the download account password

        Returns
        -------
        str
            the sorted names of the remote files
        """
        
        def index(x=None, y=None, nx=3, ny=3, reverse=False):
            if reverse:
                pattern = '{c:0{n}d}{id}'
            else:
                pattern = '{id}{c:0{n}d}'
            if x is not None:
                xf = pattern.format(id='W' if x < 0 else 'E', c=abs(x), n=nx)
            else:
                xf = ''
            if y is not None:
                yf = pattern.format(id='S' if y < 0 else 'N', c=abs(y), n=ny)
            else:
                yf = ''
            return yf, xf
        
        def cop_dem_remotes(extent, arcsecs):
            lat, lon = self.intrange(extent, step=1)
            indices = [index(x, y, nx=3, ny=2)
                       for x in lon for y in lat]
            base = 'Copernicus_DSM_COG_{res}_{0}_00_{1}_00_DEM'
            skeleton = '{base}/{base}.tif'.format(base=base)
            remotes = [skeleton.format(res=arcsecs, *item) for item in indices]
            return remotes
        
        if demType == 'SRTM 1Sec HGT':
            lat, lon = self.intrange(extent, step=1)
            remotes = ['{0}{1}.SRTMGL1.hgt.zip'.format(*index(x, y, nx=3, ny=2))
                       for x in lon for y in lat]
        
        elif demType == 'GETASSE30':
            lat, lon = self.intrange(extent, step=15)
            remotes = ['{0}{1}.zip'.format(*index(x, y, nx=3, ny=2, reverse=True))
                       for x in lon for y in lat]
        
        elif demType == 'TDX90m':
            lat, lon = self.intrange(extent, step=1)
            remotes = []
            for x in lon:
                xr = abs(x) // 10 * 10
                for y in lat:
                    yf, xf = index(x=x, y=y, nx=3, ny=2)
                    remotes.append('90mdem/DEM/{y}/{hem}{xr:03d}/TDM1_DEM__30_{y}{x}.zip'
                                   .format(x=xf, xr=xr, y=yf, hem=xf[0]))
        
        elif demType == 'AW3D30':
            remotes = []
            lat, lon = self.intrange(extent, step=1)
            for x in lon:
                for y in lat:
                    remotes.append(
                        '{0}{1}/{2}{3}.tar.gz'.format(*index(x // 5 * 5, y // 5 * 5),
                                                      *index(x, y)))
        
        elif demType == 'SRTM 3Sec':
            lat = range(int((60 - float(extent['ymin'])) // 5) + 1,
                        int((60 - float(extent['ymax'])) // 5) + 2)
            lon = range(int((float(extent['xmin']) + 180) // 5) + 1,
                        int((float(extent['xmax']) + 180) // 5) + 2)
            remotes = ['srtm_{:02d}_{:02d}.zip'.format(x, y) for x in lon for y in lat]
        
        elif demType in ['Copernicus 10m EEA DEM',
                         'Copernicus 30m Global DEM II',
                         'Copernicus 90m Global DEM II']:
            lat, lon = self.intrange(extent, step=1)
            indices = [''.join(index(x, y, nx=3, ny=2))
                       for x in lon for y in lat]
            
            ftp = ImplicitFTP_TLS()
            parsed = urlparse(self.config[demType]['url'])
            host = parsed.netloc
            path = parsed.path
            ftp.connect(host=host, port=self.config[demType]['port'])
            ftp.login(username, password)
            ftp.cwd(path)
            
            obj = io.BytesIO()
            ftp.retrbinary('RETR mapping.csv', obj.write)
            obj = obj.getvalue().decode('utf-8').splitlines()
            
            ids = []
            stream = csv.reader(obj, delimiter=';')
            for row in stream:
                if row[1] + row[2] in indices:
                    ids.append(row[0])
            
            remotes = []
            remotes_base = []
            
            def ftp_search(target, files):
                pattern = '|'.join(files)
                if target.endswith('/'):
                    content = ftp.nlst(target)
                    for item in content:
                        ftp_search(target + '/' + item, files)
                else:
                    if target.endswith('.tar') and re.search(pattern, target):
                        base = os.path.basename(target)
                        if base not in remotes_base:
                            remotes.append(target)
                            remotes_base.append(base)
            
            ftp_search(path + '/', ids)
            ftp.quit()
        
        elif demType == 'Copernicus 30m Global DEM':
            remotes = cop_dem_remotes(extent, arcsecs=10)
        
        elif demType == 'Copernicus 90m Global DEM':
            remotes = cop_dem_remotes(extent, arcsecs=30)
        
        else:
            raise ValueError('unknown demType: {}'.format(demType))
        
        return sorted(remotes)


def getasse30_hdr(fname):
    """
    create an ENVI HDR file for zipped GETASSE30 DEM tiles
    
    Parameters
    ----------
    fname: str
        the name of the zipped tile

    Returns
    -------

    """
    basename = os.path.basename(fname)
    pattern = r'(?P<lat>[0-9]{2})' \
              '(?P<ns>[A-Z])' \
              '(?P<lon>[0-9]{3})' \
              '(?P<ew>[A-Z]).zip'
    match = re.search(pattern, basename).groupdict()
    
    lon = float(match['lon'])
    if match['ew'] == 'W':
        lon *= -1
    lat = float(match['lat'])
    if match['ns'] == 'S':
        lat *= -1
    posting = 30 / 3600  # 30 arc seconds
    pixels = 1800
    
    map_info = ['Geographic Lat/Lon', '1.0000', '1.0000',
                str(lon),
                str(lat + pixels * posting),
                str(posting),
                str(posting),
                'WGS-84', 'units=Degrees']
    
    with zf.ZipFile(fname, 'a') as zip:
        files = zip.namelist()
        hdr = basename.replace('.zip', '.hdr')
        if hdr not in files:
            with HDRobject() as obj:
                obj.samples = pixels
                obj.lines = pixels
                obj.byte_order = 1
                obj.data_type = 2
                obj.map_info = '{{{}}}'.format(','.join(map_info))
                obj.coordinate_system_string = crsConvert(4326, 'wkt')
                zip.writestr(hdr, str(obj))


def get_egm_lookup(geoid, software):
    """
    Download lookup tables for converting EGM geoid heights to WGS84 ellipsoid heights.
    
    Parameters
    ----------
    geoid: str
        the geoid model; current options:
        
        - SNAP: 'EGM96'
        - PROJ: 'EGM96', 'EGM2008'
    software: str
        the software for which to download the EGM lookup
        
        - SNAP: default directory: ``~/.snap/auxdata/dem/egm96``; URL:
        
          * https://step.esa.int/auxdata/dem/egm96/ww15mgh_b.zip
        - PROJ: requires ``PROJ_LIB`` environment variable to be set as download directory; URLs:
        
          * https://download.osgeo.org/proj/vdatum/egm96_15/egm96_15.gtx
          * https://download.osgeo.org/proj/vdatum/egm08_25/egm08_25.gtx

    Returns
    -------

    """
    if software == 'SNAP':
        try:
            auxdatapath = ExamineSnap().auxdatapath
        except AttributeError:
            auxdatapath = os.path.join(os.path.expanduser('~'), '.snap', 'auxdata')
        local = os.path.join(auxdatapath, 'dem', 'egm96', 'ww15mgh_b.zip')
        os.makedirs(os.path.dirname(local), exist_ok=True)
        if not os.path.isfile(local):
            remote = 'https://step.esa.int/auxdata/dem/egm96/ww15mgh_b.zip'
            log.info('{} <<-- {}'.format(local, remote))
            r = requests.get(remote)
            r.raise_for_status()
            with open(local, 'wb') as out:
                out.write(r.content)
    
    elif software == 'PROJ':
        gtx_lookup = {'EGM96': 'us_nga_egm96_15.tif',
                      'EGM2008': 'us_nga_egm08_25.tif'}
        gtx_remote = 'https://cdn.proj.org/' + gtx_lookup[geoid]
        
        proj_lib = os.environ.get('PROJ_LIB')
        if proj_lib is not None:
            gtx_local = os.path.join(proj_lib, os.path.basename(gtx_remote))
            if not os.path.isfile(gtx_local):
                if not os.access(proj_lib, os.W_OK):
                    raise OSError("cannot write to 'PROJ_LIB' path: {}".format(proj_lib))
                log.info('{} <<-- {}'.format(gtx_local, gtx_remote))
                r = requests.get(gtx_remote)
                r.raise_for_status()
                with open(gtx_local, 'wb') as out:
                    out.write(r.content)
        else:
            raise RuntimeError("environment variable 'PROJ_LIB' not set")
    else:
        raise TypeError("software must be either 'SNAP' or 'PROJ'")


class ImplicitFTP_TLS(ftplib.FTP_TLS):
    """
    FTP_TLS subclass that automatically wraps sockets in SSL to support implicit FTPS.
    taken from https://stackoverflow.com/a/36049814
    """
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._sock = None
    
    @property
    def sock(self):
        """Return the socket."""
        return self._sock
    
    @sock.setter
    def sock(self, value):
        """When modifying the socket, ensure that it is ssl wrapped."""
        if value is not None and not isinstance(value, ssl.SSLSocket):
            value = self.context.wrap_socket(value)
        self._sock = value
