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
import os
import re
import csv
import ssl
import numpy
import fnmatch
import ftplib
import requests
import zipfile as zf
from math import ceil, floor
from urllib.parse import urlparse
from lxml import etree as ET
from pyroSAR.examine import ExamineSnap
from spatialist.raster import Raster, Dtype
from spatialist.vector import bbox
from spatialist.ancillary import dissolve, finder
from spatialist.auxil import gdalbuildvrt, crsConvert, gdalwarp
from spatialist.envi import HDRobject
from osgeo import gdal

import logging

log = logging.getLogger(__name__)


def dem_autoload(geometries, demType, vrt=None, buffer=None, username=None, password=None,
                 product='dem', nodata=None, dst_nodata=None, hide_nodata=False, crop=True):
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
        (optional) the username for services requiring registration
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
    
    nodata: int or float or None
        the no data value of the source files.
    dst_nodata: int or float or None
        the nodata value of the VRT file.
    hide_nodata: bool
        hide the VRT no data value?
    crop: bool
        crop to the provided geometries (or return the full extent of the DEM tiles)?
        Argument `buffer` is ignored if set to `False`.
    
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
        return handler.load(dem_type=demType,
                            username=username,
                            password=password,
                            vrt=vrt,
                            buffer=buffer,
                            product=product,
                            nodata=nodata,
                            vrt_nodata=dst_nodata,
                            hide_nodata=hide_nodata,
                            crop=crop)


def dem_create(src, dst, t_srs=None, tr=None, resampling_method='bilinear', threads=None,
               geoid_convert=False, geoid='EGM96', outputBounds=None, nodata=None,
               dtype=None, pbar=False):
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
    nodata: int or float or str or None
        the no data value of the source and destination files.
        Can be used if no source nodata value can be read or to override it.
        A special string 'None' can be used to skip reading the value from the source file.
    dtype: str or None
        override the data type of the written file; Default None: use same type as source data.
        Data type notations of GDAL (e.g. `Float32`) and numpy (e.g. `int8`) are supported.
    pbar: bool
        add a progressbar?
    
    Returns
    -------

    """
    
    with Raster(src) as ras:
        if nodata is None:
            nodata = ras.nodata
        epsg_in = ras.epsg
    
    if nodata is None:
        raise RuntimeError('the nodata value could not be read from the source file. Please explicitly define it.')
    
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
    
    def __find_first(self, dem_type, product):
        outdir = os.path.join(self.auxdatapath, 'dem', dem_type)
        vsi = self.config[dem_type]['vsi']
        pattern = fnmatch.translate(self.config[dem_type]['pattern'][product])
        for root, dirs, files in os.walk(outdir):
            for file in files:
                if vsi is None:
                    if re.search(pattern, file):
                        return os.path.join(root, file)
                else:
                    if re.search(r'\.(?:zip|tar(\.gz)?)$', file):
                        fname = os.path.join(root, file)
                        content = finder(fname, [pattern], regex=True)
                        if len(content) > 0:
                            if dem_type == 'GETASSE30':
                                getasse30_hdr(fname)
                            return vsi + content[0]
    
    @staticmethod
    def __buildvrt(tiles, vrtfile, pattern, vsi, extent, src_nodata=None,
                   dst_nodata=None, hide_nodata=False, resolution=None, crop=True):
        """
        Build a VRT mosaic from DEM tiles. The VRT is cropped to the specified `extent` but the pixel grid
        of the source files is preserved and no resampling/shifting is applied.
        
        Parameters
        ----------
        tiles: list[str]
            a list of DEM files or compressed archives containing DEM files
        vrtfile: str
            the output VRT filename
        pattern: str
            the search pattern for finding DEM tiles in compressed archives
        vsi: str or None
            the GDAL VSI directive to prepend the DEM tile name, e.g. /vsizip/ or /vsitar/
        extent: dict:
            a dictionary with keys `xmin`, `ymin`, `xmax` and `ymax`
        src_nodata: int or float or None
            the nodata value of the source DEM tiles; default None: read the value from the first item in `tiles`
        dst_nodata: int or float or None
            the nodata value of the output VRT file; default None: do not define a nodata value
        hide_nodata: bool
            hide the nodata value of the ouptu VRT file?
        resolution: int or float or None
            the spatial resolution of the source DEM tiles; default None: read the value from the first item in `tiles`

        Returns
        -------

        """
        if vsi is not None and not tiles[0].endswith('.tif'):
            locals = [vsi + x for x in dissolve([finder(x, [pattern]) for x in tiles])]
        else:
            locals = tiles
        with Raster(locals[0]) as ras:
            if src_nodata is None:
                src_nodata = ras.nodata
            if resolution is None:
                xres, yres = ras.res
            else:
                xres, yres = resolution
        opts = {'srcNodata': src_nodata,
                'targetAlignedPixels': True,  # preserve source file pixel grid
                'xRes': xres, 'yRes': yres, 'hideNodata': hide_nodata
                }
        if dst_nodata is not None:
            opts['VRTNodata'] = dst_nodata
        if crop:
            opts['outputBounds'] = (extent['xmin'], extent['ymin'],
                                    extent['xmax'], extent['ymax'])
        gdalbuildvrt(src=locals, dst=vrtfile, options=opts)
    
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
    def __create_dummy_dem():
        path = os.path.join(os.path.expanduser('~'), '.pyrosar', 'auxdata')
        os.makedirs(name=path, exist_ok=True)
        filename = os.path.join(path, 'dummy_dem.tif')
        driver = gdal.GetDriverByName('GTiff')
        dataset = driver.Create(filename, 1, 1, 1, 1)
        geo = [-180, 360, 0, 90, 0, -180]
        dataset.SetGeoTransform(geo)
        dataset.SetProjection('EPSG:4326')
        band = dataset.GetRasterBand(1)
        band.SetNoDataValue(255)
        mat = numpy.zeros(shape=(1, 1))
        band.WriteArray(mat, 0, 0)
        band.FlushCache()
        del mat
        band = None
        dataset = None
        driver = None
        return filename
    
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
        # check that URL is reachable
        r = requests.get(url)
        r.raise_for_status()
        r.close()
        
        files = list(set(filenames))
        os.makedirs(outdir, exist_ok=True)
        locals = []
        n = len(files)
        for i, file in enumerate(files):
            remote = '{}/{}'.format(url, file)
            local = os.path.join(outdir, os.path.basename(file))
            if not os.path.isfile(local):
                msg = '[{i: >{w}}/{n}] {l} <<-- {r}'
                log.info(msg.format(i=i + 1, w=len(str(n)), n=n, l=local, r=remote))
                r = requests.get(remote)
                # a tile might not exist over ocean
                if r.status_code == 404:
                    r.close()
                    continue
                r.raise_for_status()
                with open(local, 'wb') as output:
                    output.write(r.content)
                r.close()
            else:
                msg = '[{i: >{w}}/{n}] found local file: {l}'
                log.info(msg.format(i=i + 1, w=len(str(n)), n=n, l=local))
            if os.path.isfile(local):
                locals.append(local)
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
        n = len(files)
        for i, product_remote in enumerate(files):
            product_local = os.path.join(outdir, os.path.basename(product_remote))
            if not os.path.isfile(product_local):
                try:
                    targetlist = ftp.nlst(product_remote)
                except ftplib.error_temp:
                    continue
                address = '{}://{}/{}{}'.format(parsed.scheme, parsed.netloc,
                                                parsed.path + '/' if parsed.path != '' else '', product_remote)
                msg = '[{i: >{w}}/{n}] {l} <<-- {r}'
                log.info(msg.format(i=i + 1, w=len(str(n)), n=n, l=product_local, r=address))
                with open(product_local, 'wb') as myfile:
                    ftp.retrbinary('RETR {}'.format(product_remote), myfile.write)
            else:
                msg = '[{i: >{w}}/{n}] found local file: {l}'
                log.info(msg.format(i=i + 1, w=len(str(n)), n=n, l=product_local))
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
                                   'stk': '*STK.tif'},
                       'authentication': False
                       },
            'Copernicus 10m EEA DEM': {'url': 'ftps://cdsdata.copernicus.eu/DEM-datasets/COP-DEM_EEA-10-DGED/2021_1',
                                       'nodata': -32767.0,
                                       'vsi': '/vsitar/',
                                       'port': 990,
                                       'pattern': {'dem': '*DEM.tif',
                                                   'edm': '*EDM.tif',
                                                   'flm': '*FLM.tif',
                                                   'hem': '*HEM.tif',
                                                   'wbm': '*WBM.tif'},
                                       'authentication': True
                                       },
            'Copernicus 30m Global DEM': {'url': 'https://copernicus-dem-30m.s3.eu-central-1.amazonaws.com',
                                          'nodata': None,
                                          'vsi': None,
                                          'pattern': {'dem': '*DSM*'},
                                          'authentication': False
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
                            'wbm': '*WBM.tif'},
                'authentication': True
            },
            'Copernicus 90m Global DEM': {'url': 'https://copernicus-dem-90m.s3.eu-central-1.amazonaws.com',
                                          'nodata': None,
                                          'vsi': None,
                                          'pattern': {'dem': '*DSM*'},
                                          'authentication': False
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
                            'wbm': '*WBM.tif'},
                'authentication': True
            },
            'GETASSE30': {'url': 'https://step.esa.int/auxdata/dem/GETASSE30',
                          'nodata': None,
                          'vsi': '/vsizip/',
                          'pattern': {'dem': '*.GETASSE30'},
                          'authentication': False
                          },
            'SRTM 1Sec HGT': {'url': 'https://step.esa.int/auxdata/dem/SRTMGL1',
                              'nodata': -32768.0,
                              'vsi': '/vsizip/',
                              'pattern': {'dem': '*.hgt'},
                              'authentication': False
                              },
            'SRTM 3Sec': {'url': 'https://download.esa.int/step/auxdata/dem/SRTM90/tiff',
                          'nodata': -32768.0,
                          'vsi': '/vsizip/',
                          'pattern': {'dem': 'srtm*.tif'},
                          'authentication': False
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
                                   'wam': '*_WAM.tif'},
                       'authentication': True
                       }
        }
    
    def load(self, dem_type, vrt=None, buffer=None, username=None, password=None,
             product='dem', nodata=None, vrt_nodata=None, hide_nodata=False, crop=True):
        """
        obtain DEM tiles for the given geometries and either return the file names in a list
        or combine them into a VRT mosaic. The VRT is cropped to the combined extent of the geometries
        but the pixel grid of the source files is preserved and no resampling/shifting is applied.
        
        Parameters
        ----------
        dem_type: str
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
            overrides the nodata values of the source files.
            If `None`, the value of the source products is read passed on.
        vrt_nodata: int or float or None
            the nodata value of the output VRT file; default None: do not define a nodata value.
        hide_nodata: bool
            hide the nodata value of the ouptut VRT file?
        crop: bool
            crop to the provided geometries (or return the full extent of the DEM tiles)?
            Argument `buffer` is ignored if set to `False`.
        
        Returns
        -------
        list[str] or None
            the names of the obtained files or None if a VRT file was defined
        """
        keys = self.config.keys()
        if dem_type not in keys:
            raise RuntimeError("demType '{}' is not supported\n  "
                               "possible options: '{}'"
                               .format(dem_type, "', '".join(keys)))
        
        products = self.config[dem_type]['pattern'].keys()
        if product not in products:
            raise RuntimeError("product '{0}' not available for demType '{1}'\n"
                               "  options: '{2}'".format(product, dem_type, "', '".join(products)))
        
        outdir = os.path.join(self.auxdatapath, 'dem', dem_type)
        
        candidates = []
        for geo in self.geometries:
            corners = self.__applybuffer(geo.extent, buffer)
            candidates.extend(self.remote_ids(corners, demType=dem_type,
                                              username=username, password=password))
        
        if dem_type in ['AW3D30', 'TDX90m', 'Copernicus 10m EEA DEM',
                        'Copernicus 30m Global DEM II', 'Copernicus 90m Global DEM II']:
            port = 0
            if 'port' in self.config[dem_type].keys():
                port = self.config[dem_type]['port']
            locals = self.__retrieve_ftp(url=self.config[dem_type]['url'],
                                         filenames=candidates,
                                         outdir=outdir, username=username,
                                         password=password, port=port)
        else:
            locals = self.__retrieve(url=self.config[dem_type]['url'],
                                     filenames=candidates, outdir=outdir)
        
        if len(locals) == 0:
            if product != 'dem':
                msg = 'could not find {0} {1} tiles for the area of interest'
                raise RuntimeError(msg.format(dem_type, product))
            locals = [self.__create_dummy_dem()]
            ref = self.__find_first(dem_type=dem_type, product=product)
            with Raster(ref) as ras:
                resolution = ras.res
        else:
            resolution = None
        
        if dem_type == 'GETASSE30':
            for item in locals:
                getasse30_hdr(item)
        
        if nodata is None:
            if product == 'dem':
                nodata = self.config[dem_type]['nodata']
        
        if vrt is not None:
            self.__buildvrt(tiles=locals, vrtfile=vrt,
                            pattern=self.config[dem_type]['pattern'][product],
                            vsi=self.config[dem_type]['vsi'],
                            extent=self.__commonextent(buffer),
                            src_nodata=nodata, dst_nodata=vrt_nodata,
                            hide_nodata=hide_nodata,
                            resolution=resolution, crop=crop)
        else:
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
            the download account username
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
            candidates = [skeleton.format(res=arcsecs, *item) for item in indices]
            remotes = []
            for candidate in candidates:
                response = requests.get(self.config[demType]['url'],
                                        params={'prefix': candidate})
                xml = ET.fromstring(response.content)
                content = xml.findall('.//Contents', namespaces=xml.nsmap)
                if len(content) > 0:
                    remotes.append(candidate)
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
            
            outdir = os.path.join(self.auxdatapath, 'dem', demType)
            mapping = os.path.join(outdir, 'mapping.csv')
            mapping2 = os.path.join(outdir, 'mapping_append.csv')
            
            def ftp_search(ftp, target):
                out = []
                if target.endswith('/'):
                    print(target)
                    content = ftp.nlst(target)
                    for item in content:
                        out.extend(ftp_search(ftp, target + item))
                else:
                    if target.endswith('DEM.tar'):
                        out.append(target.encode('latin-1').decode('utf-8'))
                return out
            
            def ftp_connect(host, path, username, password, port=990):
                ftp = ImplicitFTP_TLS()
                ftp.connect(host=host, port=port)
                ftp.login(username, password)
                ftp.cwd(path)
                return ftp
            
            if not os.path.isfile(mapping2):
                parsed = urlparse(self.config[demType]['url'])
                host = parsed.netloc
                path = parsed.path
                ftp = None
                os.makedirs(outdir, exist_ok=True)
                if not os.path.isfile(mapping):
                    print('downloading mapping.csv')
                    ftp = ftp_connect(host, path, username, password,
                                      port=self.config[demType]['port'])
                    with open(mapping, 'wb') as myfile:
                        ftp.retrbinary('RETR mapping.csv', myfile.write)
                print('searching FTP server')
                if ftp is None:
                    ftp = ftp_connect(host, path, username, password,
                                      port=self.config[demType]['port'])
                files = ftp_search(ftp, path + '/')
                files_base = [os.path.basename(x) for x in files]
                if ftp is not None:
                    ftp.quit()
                print('matching found files with mapping.csv')
                with open(mapping) as obj:
                    reader = csv.reader(obj, delimiter=';')
                    with open(mapping2, 'w', newline='') as out:
                        writer = csv.writer(out, delimiter=';')
                        writer.writerow(next(reader))  # write header
                        for row in reader:
                            index = files_base.index(row[0])
                            row.append(files[index])
                            del files_base[index]
                            del files[index]
                            writer.writerow(row)
            remotes = []
            with open(mapping2) as obj:
                stream = csv.reader(obj, delimiter=';')
                for row in stream:
                    if row[1] + row[2] in indices:
                        remotes.append(row[-1])
        
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


def get_dem_options(require_auth=None):
    """
    Get the names of all supported DEM type options.
    
    Parameters
    ----------
    require_auth: bool or None
        only return options that do/don't require authentication. Default None: return all options.

    Returns
    -------
    list[str]
        the names of the DEM options
    """
    out = []
    # create a dummy vector geometry for initializing the DEMHandler
    ext = {'xmin': -44, 'xmax': -43, 'ymin': 30, 'ymax': 31}
    with bbox(coordinates=ext, crs=4326) as vec:
        with DEMHandler(geometries=[vec]) as handler:
            for key, properties in handler.config.items():
                if require_auth is None:
                    out.append(key)
                else:
                    if require_auth == properties['authentication']:
                        out.append(key)
            return sorted(out)


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
        - PROJ: requires the ``PROJ_LIB`` environment variable to be set as download directory; URLs:
        
          * https://cdn.proj.org/us_nga_egm96_15.tif
          * https://cdn.proj.org/us_nga_egm08_25.tif

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
            r.close()
    
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
                r.close()
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
