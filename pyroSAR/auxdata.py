###############################################################################
# tools for handling auxiliary data in software pyroSAR

# Copyright (c) 2019-2025, the pyroSAR Developers.

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
from lxml import etree
from math import ceil, floor
from urllib.parse import urlparse
from lxml import etree as ET
from packaging import version
from pyroSAR.examine import ExamineSnap
from pyroSAR.ancillary import Lock
from spatialist.raster import Raster, Dtype
from spatialist.vector import bbox
from spatialist.ancillary import dissolve, finder
from spatialist.auxil import gdalbuildvrt, crsConvert, gdalwarp
from spatialist.envi import HDRobject
from osgeo import gdal

import logging

log = logging.getLogger(__name__)


def dem_autoload(geometries, demType, vrt=None, buffer=None, username=None,
                 password=None, product='dem', crop=True, lock_timeout=600):
    """
    obtain all relevant DEM tiles for selected geometries and optionally mosaic them in a VRT.

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

          * url: https://download.esa.int/step/auxdata/dem/SRTM90/tiff
          * height reference: EGM96

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
          * 'edm': Editing Mask
          * 'flm': Filling Mask
          * 'hem': Height Error Mask
          * 'wbm': Water Body Mask
        
        - 'Copernicus 30m Global DEM II'
        
          * 'dem': the actual Digital Elevation Model
          * 'edm': editing mask
          * 'flm': filling mask
          * 'hem': height error mask
          * 'wbm': water body mask
        
        - 'Copernicus 90m Global DEM'
        
          * 'dem': the actual Digital Elevation Model
          * 'edm': Editing Mask
          * 'flm': Filling Mask
          * 'hem': Height Error Mask
          * 'wbm': Water Body Mask
        
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
    
    crop: bool
        crop to the provided geometries (or return the full extent of the DEM tiles)?
    lock_timeout: int
        how long to wait to acquire a lock on the downloaded files?
    
    Returns
    -------
    list[str] or None
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
                            crop=crop,
                            lock_timeout=lock_timeout)


def dem_create(src, dst, t_srs=None, tr=None, threads=None,
               geoid_convert=False, geoid='EGM96', nodata=None,
               resampleAlg='bilinear', dtype=None, pbar=False,
               **kwargs):
    """
    Create a new DEM GeoTIFF file and optionally convert heights from geoid to ellipsoid.
    This is basically a convenience wrapper around :func:`osgeo.gdal.Warp` via :func:`spatialist.auxil.gdalwarp`.
    The following argument defaults deviate from those of :func:`osgeo.gdal.WarpOptions`:
    
    - `format` is set to 'GTiff'
    - `resampleAlg` is set to 'bilinear'
    - `targetAlignedPixels` is set to 'True'
    
    
    Parameters
    ----------
    src: str
        the input dataset, e.g. a VRT from function :func:`dem_autoload`
    dst: str
        the output dataset
    t_srs: None, int, str or osgeo.osr.SpatialReference
        A target geographic reference system in WKT, EPSG, PROJ4 or OPENGIS format.
        See function :func:`spatialist.auxil.crsConvert()` for details.
        Default (None): use the crs of ``src``.
    tr: None or tuple[int or float]
        the target resolution as (xres, yres)
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
    nodata: int or float or str or None
        the no data value of the source and destination files.
        Can be used if no source nodata value can be read or to override it.
        A special string 'None' can be used to skip reading the value from the source file.
    resampleAlg: str
        the resampling algorithm tu be used. See here for options:
        https://gdal.org/programs/gdalwarp.html#cmdoption-gdalwarp-r
    dtype: str or None
        override the data type of the written file; Default None: use same type as source data.
        Data type notations of GDAL (e.g. `Float32`) and numpy (e.g. `int8`) are supported.
        See :class:`spatialist.raster.Dtype`.
    pbar: bool
        add a progressbar?
    **kwargs
        additional keyword arguments to be passed to :func:`spatialist.auxil.gdalwarp`.
        See :func:`osgeo.gdal.WarpOptions` for options. The following arguments cannot
        be set as they are controlled internally:
        
        - `xRes`, `yRes`: controlled via argument `tr`
        - `srcSRS`, `dstSRS`: controlled via the CRS of `src` and arguments `t_srs`, `geoid`, `geoid_convert`
        - `srcNodata`, `dstNodata`: controlled via argument `nodata`
        - `outputType`: controlled via argument `dtype`
        - `multithread` controlled via argument `threads`
    
    Returns
    -------

    """
    
    vrt_check_sources(src)
    
    with Raster(src) as ras:
        if nodata is None:
            nodata = ras.nodata
        if tr is None:
            tr = ras.res
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
                     'resampleAlg': resampleAlg,
                     'xRes': tr[0], 'yRes': tr[1],
                     'targetAlignedPixels': True}
    
    if dtype is not None:
        gdalwarp_args['outputType'] = Dtype(dtype).gdalint
    
    if geoid_convert:
        geoid_epsg = {'EGM96': 5773,
                      'EGM2008': 3855}
        if geoid in geoid_epsg.keys():
            epsg = geoid_epsg[geoid]
            gdalwarp_args['srcSRS'] += '+{}'.format(epsg)
            # the following line is a workaround for older GDAL versions that did not
            # support compound EPSG codes. See https://github.com/OSGeo/gdal/pull/4639.
            if version.parse(gdal.__version__) < version.parse('3.4.0'):
                gdalwarp_args['srcSRS'] = crsConvert(gdalwarp_args['srcSRS'], 'proj4')
        else:
            raise RuntimeError('geoid model not yet supported')
        try:
            get_egm_lookup(geoid=geoid, software='PROJ')
        except OSError as e:
            errstr = str(e)
            raise RuntimeError(errstr)
    
    locked = ['xRes', 'yRes', 'srcSRS', 'dstSRS', 'srcNodata',
              'dstNodata', 'outputType', 'multithread']
    for key, val in kwargs.items():
        if key not in locked:
            gdalwarp_args[key] = val
        else:
            msg = "argument '{}' cannot be set via kwargs as it is set internally."
            raise RuntimeError(msg.format(key))
    try:
        if not os.path.isfile(dst):
            message = 'creating mosaic'
            crs = gdalwarp_args['dstSRS']
            if crs != 'EPSG:4326':
                message += ' and reprojecting to {}'.format(crs)
            log.info(f'{message}: {dst}')
            gdalwarp(src=src, dst=dst, pbar=pbar, **gdalwarp_args)
        else:
            log.info(f'mosaic already exists: {dst}')
    except Exception:
        if os.path.isfile(dst):
            os.remove(dst)
        raise
    finally:
        gdal.SetConfigOption('GDAL_NUM_THREADS', threads_system)


class DEMHandler:
    """
    An interface to obtain DEM data for selected geometries.
    The files are downloaded into the ESA SNAP auxiliary data directory structure.
    This class is the foundation for the convenience function :func:`~pyroSAR.auxdata.dem_autoload`.
    
    Parameters
    ----------
    geometries: list[spatialist.vector.Vector]
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
                   dst_nodata=None, hide_nodata=False, resolution=None,
                   tap=True, dst_datatype=None):
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
        extent: dict
            a dictionary with keys `xmin`, `ymin`, `xmax` and `ymax`
        src_nodata: int or float or None
            the nodata value of the source DEM tiles; default None: read the value from the first item in `tiles`
        dst_nodata: int or float or None
            the nodata value of the output VRT file.
            Default None: do not define a nodata value and use `src_nodata` instead.
        hide_nodata: bool
            hide the nodata value of the output VRT file?
        resolution: tuple[int or float] or None
            the spatial resolution (X, Y) of the source DEM tiles.
            Default None: read the value from the first item in `tiles`
        tap: bool
            align target pixels?
        dst_datatype: int or str or None
            the VRT data type as supported by :class:`spatialist.raster.Dtype`.
            Default None: use the same data type as the source files.
        
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
                'targetAlignedPixels': tap,
                'xRes': xres, 'yRes': yres, 'hideNodata': hide_nodata
                }
        if dst_nodata is not None:
            opts['VRTNodata'] = dst_nodata
        opts['outputBounds'] = (extent['xmin'], extent['ymin'],
                                extent['xmax'], extent['ymax'])
        
        gdalbuildvrt(src=locals, dst=vrtfile, **opts)
        if dst_datatype is not None:
            datatype = Dtype(dst_datatype).gdalstr
            tree = etree.parse(source=vrtfile)
            band = tree.find(path='VRTRasterBand')
            band.attrib['dataType'] = datatype
            tree.write(file=vrtfile, pretty_print=True,
                       xml_declaration=False, encoding='utf-8')
    
    def __commonextent(self, buffer=None):
        """
        
        Parameters
        ----------
        buffer: int or float or None

        Returns
        -------
        dict
        """
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
    def __create_dummy_dem(filename, extent):
        """
        Create a dummy file which spans the given extent and
        is 1x1 pixels large to be as small as possible.
        This file is used to create dummy DEMs over ocean.
        """
        driver = gdal.GetDriverByName('GTiff')
        dataset = driver.Create(filename, 1, 1, 1, 1)
        geo = [
            extent['xmin'],
            extent['xmax'] - extent['xmin'],
            0,
            extent['ymax'],
            0,
            extent['ymin'] - extent['ymax']  # negative
        ]
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
    
    def __get_resolution(self, dem_type, y):
        """
        
        Parameters
        ----------
        dem_type: str
            the DEM type
        y: int or float
            the latitude for which to get the resolution

        Returns
        -------
        tuple
            (xres, yres)
        """
        for key, val in self.config[dem_type]['resolution'].items():
            ymin, ymax = [int(y) for y in key.split('-')]
            if ymin <= abs(y) <= ymax:
                return val
    
    @staticmethod
    def __retrieve(url, filenames, outdir, lock_timeout=600):
        # check that base URL is reachable
        url_parse = urlparse(url)
        url_base = url_parse.scheme + '://' + url_parse.netloc
        r = requests.get(url_base)
        r.raise_for_status()
        r.close()
        
        files = list(set(filenames))
        os.makedirs(outdir, exist_ok=True)
        locals = []
        n = len(files)
        for i, file in enumerate(files):
            remote = '{}/{}'.format(url, file)
            local = os.path.join(outdir, os.path.basename(file))
            with Lock(local, timeout=lock_timeout):
                if not os.path.isfile(local):
                    r = requests.get(remote)
                    # a tile might not exist over ocean
                    if r.status_code == 404:
                        r.close()
                        continue
                    msg = '[{i: >{w}}/{n}] {l} <<-- {r}'
                    log.info(msg.format(i=i + 1, w=len(str(n)), n=n, l=local, r=remote))
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
    
    @staticmethod
    def __retrieve_ftp(url, filenames, outdir, username, password,
                       port=0, offline=False, lock_timeout=600):
        files = list(set(filenames))
        os.makedirs(outdir, exist_ok=True)
        
        parsed = urlparse(url)
        timeout = 100
        if not offline:
            if parsed.scheme == 'ftpes':
                ftp = ftplib.FTP_TLS(host=parsed.netloc, timeout=timeout)
                try:
                    ftp.login(username, password)  # login anonymously before securing control channel
                except ftplib.error_perm as e:
                    raise RuntimeError(str(e))
                ftp.prot_p()  # switch to secure data connection.. IMPORTANT! Otherwise, only the user and password is encrypted and not all the file data.
            elif parsed.scheme == 'ftps':
                ftp = ImplicitFTP_TLS()
                ftp.connect(host=parsed.netloc, timeout=timeout, port=port)
                ftp.login(username, password)
            else:
                ftp = ftplib.FTP(host=parsed.netloc, timeout=timeout)
                ftp.login()
            if parsed.path != '':
                ftp.cwd(parsed.path)
        else:
            ftp = None
        locals = []
        n = len(files)
        for i, remote in enumerate(files):
            local = os.path.join(outdir, os.path.basename(remote))
            with Lock(local, timeout=lock_timeout):
                if not os.path.isfile(local) and not offline:
                    try:
                        targetlist = ftp.nlst(remote)
                    except ftplib.error_temp:
                        continue
                    address = '{}://{}/{}{}'.format(parsed.scheme, parsed.netloc,
                                                    parsed.path + '/' if parsed.path != '' else '',
                                                    remote)
                    msg = '[{i: >{w}}/{n}] {l} <<-- {r}'
                    log.info(msg.format(i=i + 1, w=len(str(n)), n=n, l=local, r=address))
                    with open(local, 'wb') as myfile:
                        ftp.retrbinary('RETR {}'.format(remote), myfile.write)
                else:
                    msg = '[{i: >{w}}/{n}] found local file: {l}'
                    log.info(msg.format(i=i + 1, w=len(str(n)), n=n, l=local))
            if os.path.isfile(local):
                locals.append(local)
        if ftp is not None:
            ftp.close()
        return sorted(locals)
    
    @property
    def config(self):
        return {
            'AW3D30': {'url': 'ftp://ftp.eorc.jaxa.jp/pub/ALOS/ext1/AW3D30/release_v1804',
                       'nodata': {'dem': -9999,
                                  'msk': 3,
                                  'stk': 0},
                       'resolution': {'0-90': (1 / 3600, 1 / 3600)},
                       'tilesize': 1,
                       'area_or_point': 'area',
                       'vsi': '/vsitar/',
                       'pattern': {'dem': '*DSM.tif',
                                   'msk': '*MSK.tif',
                                   'stk': '*STK.tif'},
                       'datatype': {'dem': 'Int16',
                                    'msk': 'Byte',
                                    'stk': 'Byte'},
                       'authentication': False
                       },
            'Copernicus 10m EEA DEM': {'url': 'ftps://cdsdata.copernicus.eu/DEM-datasets/COP-DEM_EEA-10-DGED/2021_1',
                                       'nodata': {'dem': -32767.0,
                                                  'edm': 8,
                                                  'flm': 1,
                                                  'hem': -32767.0,
                                                  'wbm': 1},
                                       'resolution': {'0-50': (1 / 9000, 1 / 9000),
                                                      '50-60': (1.5 / 9000, 1 / 9000),
                                                      '60-70': (2 / 9000, 1 / 9000),
                                                      '70-80': (3 / 9000, 1 / 9000),
                                                      '80-85': (5 / 9000, 1 / 9000),
                                                      '85-90': (10 / 9000, 1 / 9000)},
                                       'tilesize': 1,
                                       'area_or_point': 'point',
                                       'vsi': '/vsitar/',
                                       'port': 990,
                                       'pattern': {'dem': '*DEM.tif',
                                                   'edm': '*EDM.tif',
                                                   'flm': '*FLM.tif',
                                                   'hem': '*HEM.tif',
                                                   'wbm': '*WBM.tif'},
                                       'datatype': {'dem': 'Float32',
                                                    'edm': 'Byte',
                                                    'flm': 'Byte',
                                                    'hem': 'Float32',
                                                    'wbm': 'Byte'},
                                       'authentication': True
                                       },
            'Copernicus 30m Global DEM': {'url': 'https://copernicus-dem-30m.s3.eu-central-1.amazonaws.com',
                                          'nodata': {'dem': -32767.0,
                                                     'edm': 8,
                                                     'flm': 1,
                                                     'hem': -32767.0,
                                                     'wbm': 1},
                                          'resolution': {'0-50': (1 / 3600, 1 / 3600),
                                                         '50-60': (1.5 / 3600, 1 / 3600),
                                                         '60-70': (2 / 3600, 1 / 3600),
                                                         '70-80': (3 / 3600, 1 / 3600),
                                                         '80-85': (5 / 3600, 1 / 3600),
                                                         '85-90': (10 / 3600, 1 / 3600)},
                                          'tilesize': 1,
                                          'area_or_point': 'point',
                                          'vsi': None,
                                          'pattern': {'dem': '*DEM.tif',
                                                      'edm': '*EDM.tif',
                                                      'flm': '*FLM.tif',
                                                      'hem': '*HEM.tif',
                                                      'wbm': '*WBM.tif'},
                                          'datatype': {'dem': 'Float32',
                                                       'edm': 'Byte',
                                                       'flm': 'Byte',
                                                       'hem': 'Float32',
                                                       'wbm': 'Byte'},
                                          'authentication': False
                                          },
            'Copernicus 30m Global DEM II': {
                'url': 'ftps://cdsdata.copernicus.eu/DEM-datasets/COP-DEM_GLO-30-DGED/2021_1',
                'nodata': {'dem': -32767.0,
                           'edm': 8,
                           'flm': 1,
                           'hem': -32767.0,
                           'wbm': 1},
                'resolution': {'0-50': (1 / 3600, 1 / 3600),
                               '50-60': (1.5 / 3600, 1 / 3600),
                               '60-70': (2 / 3600, 1 / 3600),
                               '70-80': (3 / 3600, 1 / 3600),
                               '80-85': (5 / 3600, 1 / 3600),
                               '85-90': (10 / 3600, 1 / 3600)},
                'tilesize': 1,
                'area_or_point': 'point',
                'vsi': '/vsitar/',
                'port': 990,
                'pattern': {'dem': '*DEM.tif',
                            'edm': '*EDM.tif',
                            'flm': '*FLM.tif',
                            'hem': '*HEM.tif',
                            'wbm': '*WBM.tif'},
                'datatype': {'dem': 'Float32',
                             'edm': 'Byte',
                             'flm': 'Byte',
                             'hem': 'Float32',
                             'wbm': 'Byte'},
                'authentication': True
            },
            'Copernicus 90m Global DEM': {'url': 'https://copernicus-dem-90m.s3.eu-central-1.amazonaws.com',
                                          'nodata': {'dem': -32767.0,
                                                     'edm': 8,
                                                     'flm': 1,
                                                     'hem': -32767.0,
                                                     'wbm': 1},
                                          'resolution': {'0-50': (1 / 1200, 1 / 1200),
                                                         '50-60': (1.5 / 1200, 1 / 1200),
                                                         '60-70': (2 / 1200, 1 / 1200),
                                                         '70-80': (3 / 1200, 1 / 1200),
                                                         '80-85': (5 / 1200, 1 / 1200),
                                                         '85-90': (10 / 1200, 1 / 1200)},
                                          'tilesize': 1,
                                          'area_or_point': 'point',
                                          'vsi': None,
                                          'pattern': {'dem': '*DEM.tif',
                                                      'edm': '*EDM.tif',
                                                      'flm': '*FLM.tif',
                                                      'hem': '*HEM.tif',
                                                      'wbm': '*WBM.tif'},
                                          'datatype': {'dem': 'Float32',
                                                       'edm': 'Byte',
                                                       'flm': 'Byte',
                                                       'hem': 'Float32',
                                                       'wbm': 'Byte'},
                                          'authentication': False
                                          },
            'Copernicus 90m Global DEM II': {
                'url': 'ftps://cdsdata.copernicus.eu/DEM-datasets/COP-DEM_GLO-90-DGED/2021_1',
                'nodata': {'dem': -32767.0,
                           'edm': 8,
                           'flm': 1,
                           'hem': -32767.0,
                           'wbm': 1},
                'resolution': {'0-50': (1 / 1200, 1 / 1200),
                               '50-60': (1.5 / 1200, 1 / 1200),
                               '60-70': (2 / 1200, 1 / 1200),
                               '70-80': (3 / 1200, 1 / 1200),
                               '80-85': (5 / 1200, 1 / 1200),
                               '85-90': (10 / 1200, 1 / 1200)},
                'tilesize': 1,
                'area_or_point': 'point',
                'vsi': '/vsitar/',
                'port': 990,
                'pattern': {'dem': '*DEM.tif',
                            'edm': '*EDM.tif',
                            'flm': '*FLM.tif',
                            'hem': '*HEM.tif',
                            'wbm': '*WBM.tif'},
                'datatype': {'dem': 'Float32',
                             'edm': 'Byte',
                             'flm': 'Byte',
                             'hem': 'Float32',
                             'wbm': 'Byte'},
                'authentication': True
            },
            'GETASSE30': {'url': 'https://step.esa.int/auxdata/dem/GETASSE30',
                          'nodata': {'dem': None},
                          'resolution': {'0-90': (15 / 1800, 15 / 1800)},
                          'tilesize': 15,
                          'area_or_point': 'area',
                          'vsi': '/vsizip/',
                          'pattern': {'dem': '*.GETASSE30'},
                          'datatype': {'dem': 'Int16'},
                          'authentication': False
                          },
            'SRTM 1Sec HGT': {'url': 'https://step.esa.int/auxdata/dem/SRTMGL1',
                              'nodata': {'dem': -32768.0},
                              'resolution': {'0-90': (1 / 3600, 1 / 3600)},
                              'tilesize': 1,
                              'area_or_point': 'point',
                              'vsi': '/vsizip/',
                              'pattern': {'dem': '*.hgt'},
                              'datatype': {'dem': 'Int16'},
                              'authentication': False
                              },
            'SRTM 3Sec': {'url': 'https://download.esa.int/step/auxdata/dem/SRTM90/tiff',
                          'nodata': {'dem': -32768.0},
                          'resolution': {'0-90': (5 / 6000, 5 / 6000)},
                          'tilesize': 5,
                          'area_or_point': 'area',
                          'vsi': '/vsizip/',
                          'pattern': {'dem': 'srtm*.tif'},
                          'datatype': {'dem': 'Int16'},
                          'authentication': False
                          },
            # 'TDX90m': {'url': 'ftpes://tandemx-90m.dlr.de',
            #            'nodata': {'dem': -32767.0,
            #                       'am2': 0,
            #                       'amp': 0,
            #                       'com': 0,
            #                       'cov': 0,
            #                       'hem': -32767.0,
            #                       'lsm': 0,
            #                       'wam': 0},
            #            'resolution': {'0-50': (1 / 1200, 1 / 1200),
            #                           '50-60': (1.5 / 1200, 1 / 1200),
            #                           '60-70': (2 / 1200, 1 / 1200),
            #                           '70-80': (3 / 1200, 1 / 1200),
            #                           '80-85': (5 / 1200, 1 / 1200),
            #                           '85-90': (10 / 1200, 1 / 1200)},
            #            'tilesize': 1,
            #            'area_or_point': 'point',
            #            'vsi': '/vsizip/',
            #            'pattern': {'dem': '*_DEM.tif',
            #                        'am2': '*_AM2.tif',
            #                        'amp': '*_AMP.tif',
            #                        'com': '*_COM.tif',
            #                        'cov': '*_COV.tif',
            #                        'hem': '*_HEM.tif',
            #                        'lsm': '*_LSM.tif',
            #                        'wam': '*_WAM.tif'},
            #            'datatype': {'dem': 'Float32',
            #                         'am2': 'UInt16',
            #                         'amp': 'UInt16',
            #                         'com': 'Byte',
            #                         'cov': 'Byte',
            #                         'hem': 'Float32',
            #                         'lsm': 'Byte',
            #                         'wam': 'Byte'},
            #            'authentication': True
            #            }
        }
    
    def load(self, dem_type, vrt=None, buffer=None, username=None,
             password=None, product='dem', crop=True, lock_timeout=600):
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
        buffer: int or float or None
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
              * 'edm': Editing Mask
              * 'flm': Filling Mask
              * 'hem': Height Error Mask
              * 'wbm': Water Body Mask
             
             - 'Copernicus 30m Global DEM II'
             
              * 'dem': the actual Digital Elevation Model
              * 'edm': Editing Mask
              * 'flm': Filling Mask
              * 'hem': Height Error Mask
              * 'wbm': Water Body Mask
             
             - 'Copernicus 90m Global DEM'
             
              * 'dem': the actual Digital Elevation Model
              * 'edm': Editing Mask
              * 'flm': Filling Mask
              * 'hem': Height Error Mask
              * 'wbm': Water Body Mask
             
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
        crop: bool
            If a VRT is created, crop it to  spatial extent of the provided geometries
            or return the full extent of the DEM tiles? In the latter case, the common
            bounding box of the geometries is expanded so that the coordinates are
            multiples of the tile size of the respective DEM option.
        lock_timeout: int
            how long to wait to acquire a lock on the downloaded files?
        
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
            corners = self.__applybuffer(extent=geo.extent, buffer=buffer)
            candidates.extend(self.remote_ids(extent=corners, dem_type=dem_type,
                                              username=username, password=password,
                                              product=product))
        
        if self.config[dem_type]['url'].startswith('ftp'):
            port = 0
            if 'port' in self.config[dem_type].keys():
                port = self.config[dem_type]['port']
            locals = self.__retrieve_ftp(url=self.config[dem_type]['url'],
                                         filenames=candidates,
                                         outdir=outdir, username=username,
                                         password=password, port=port,
                                         lock_timeout=lock_timeout)
        else:
            locals = self.__retrieve(url=self.config[dem_type]['url'],
                                     filenames=candidates, outdir=outdir,
                                     lock_timeout=lock_timeout)
        
        resolution = None
        datatype = None
        src_nodata = None
        dst_nodata = None
        tap = False
        extent = self.__commonextent(buffer=buffer)
        aop = self.config[dem_type]['area_or_point']
        res = self.__get_resolution(dem_type=dem_type, y=extent['ymin'])
        if not crop:
            f = self.config[dem_type]['tilesize']
            extent['xmin'] = floor(extent['xmin'] / f) * f
            extent['ymin'] = floor(extent['ymin'] / f) * f
            extent['xmax'] = ceil(extent['xmax'] / f) * f
            extent['ymax'] = ceil(extent['ymax'] / f) * f
        if aop == 'point':
            shift_x = res[0] / 2
            shift_y = res[1] / 2
            extent['xmin'] -= shift_x
            extent['ymin'] += shift_y
            extent['xmax'] -= shift_x
            extent['ymax'] += shift_y
        
        # special case where no DEM tiles were found because the AOI is completely over ocean
        if len(locals) == 0 and vrt is not None:
            # define a dummy file as source file
            # his file contains one pixel with a value of 0
            # nodata value is 255
            tif = vrt.replace('.vrt', '_tmp.tif')
            self.__create_dummy_dem(filename=tif, extent=extent)
            locals = [tif]
            datatype = self.config[dem_type]['datatype'][product]
            src_nodata = 0  # define the data value as nodata, so it can be overwritten in the VRT
            if product == 'dem':
                dst_nodata = 0
            else:
                dst_nodata = self.config[dem_type]['nodata'][product]
            # determine the target resolution based on minimum latitude
            resolution = self.__get_resolution(dem_type=dem_type, y=extent['ymin'])
        
        # make sure all GETASSE30 tiles get an ENVI HDR file so that they are GDAL-readable
        if dem_type == 'GETASSE30':
            for item in locals:
                getasse30_hdr(item)
        
        if vrt is not None:
            if src_nodata is None:
                src_nodata = self.config[dem_type]['nodata'][product]
            if dst_nodata is None:
                dst_nodata = 0 if product == 'dem' else None
            
            self.__buildvrt(tiles=locals, vrtfile=vrt,
                            pattern=self.config[dem_type]['pattern'][product],
                            vsi=self.config[dem_type]['vsi'],
                            extent=extent,
                            src_nodata=src_nodata, dst_nodata=dst_nodata,
                            hide_nodata=True,
                            resolution=resolution,
                            tap=tap, dst_datatype=datatype)
        else:
            return locals
    
    def remote_ids(self, extent, dem_type, product='dem', username=None, password=None):
        """
        parse the names of the remote files overlapping with an area of interest

        Parameters
        ----------
        extent: dict
            the extent of the area of interest with keys xmin, xmax, ymin, ymax
        dem_type: str
            the type fo DEM to be used
        product: str
            the sub-product to extract from the DEM product. Only needed for DEM options 'Copernicus 30m Global DEM'
            and 'Copernicus 90m Global DEM' and ignored otherwise.
        username: str or None
            the download account username
        password: str or None
            the download account password

        Returns
        -------
        str
            the sorted names of the remote files
        """
        keys = self.config.keys()
        if dem_type not in keys:
            raise RuntimeError("demType '{}' is not supported\n  "
                               "possible options: '{}'"
                               .format(dem_type, "', '".join(keys)))
        
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
        
        def cop_dem_remotes(extent, arcsecs, product='dem'):
            lat, lon = self.intrange(extent, step=1)
            indices = [index(x, y, nx=3, ny=2)
                       for x in lon for y in lat]
            base = 'Copernicus_DSM_COG_{res}_{0}_00_{1}_00'
            skeleton = '{base}_DEM/{sub}{base}_{product}.tif'
            sub = '' if product == 'dem' else 'AUXFILES/'
            base = skeleton.format(base=base, sub=sub, product=product.upper())
            candidates = [base.format(res=arcsecs, *item) for item in indices]
            remotes = []
            for candidate in candidates:
                response = requests.get(self.config[dem_type]['url'],
                                        params={'prefix': candidate})
                xml = ET.fromstring(response.content)
                content = xml.findall('.//Contents', namespaces=xml.nsmap)
                if len(content) > 0:
                    remotes.append(candidate)
            return remotes
        
        if dem_type == 'SRTM 1Sec HGT':
            lat, lon = self.intrange(extent, step=1)
            remotes = ['{0}{1}.SRTMGL1.hgt.zip'.format(*index(x, y, nx=3, ny=2))
                       for x in lon for y in lat]
        
        elif dem_type == 'GETASSE30':
            lat, lon = self.intrange(extent, step=15)
            remotes = ['{0}{1}.zip'.format(*index(x, y, nx=3, ny=2, reverse=True))
                       for x in lon for y in lat]
        
        elif dem_type == 'TDX90m':
            lat, lon = self.intrange(extent, step=1)
            remotes = []
            for x in lon:
                xr = abs(x) // 10 * 10
                for y in lat:
                    yf, xf = index(x=x, y=y, nx=3, ny=2)
                    remotes.append('DEM/{y}/{hem}{xr:03d}/TDM1_DEM__30_{y}{x}.zip'
                                   .format(x=xf, xr=xr, y=yf, hem=xf[0]))
        
        elif dem_type == 'AW3D30':
            remotes = []
            lat, lon = self.intrange(extent, step=1)
            for x in lon:
                for y in lat:
                    remotes.append(
                        '{0}{1}/{2}{3}.tar.gz'.format(*index(x // 5 * 5, y // 5 * 5),
                                                      *index(x, y)))
        
        elif dem_type == 'SRTM 3Sec':
            lat = range(int((60 - float(extent['ymin'])) // 5) + 1,
                        int((60 - float(extent['ymax'])) // 5) + 2)
            lon = range(int((float(extent['xmin']) + 180) // 5) + 1,
                        int((float(extent['xmax']) + 180) // 5) + 2)
            remotes = ['srtm_{:02d}_{:02d}.zip'.format(x, y) for x in lon for y in lat]
        
        elif dem_type in ['Copernicus 10m EEA DEM',
                          'Copernicus 30m Global DEM II',
                          'Copernicus 90m Global DEM II']:
            lat, lon = self.intrange(extent, step=1)
            indices = [''.join(index(x, y, nx=3, ny=2))
                       for x in lon for y in lat]
            
            outdir = os.path.join(self.auxdatapath, 'dem', dem_type)
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
                parsed = urlparse(self.config[dem_type]['url'])
                host = parsed.netloc
                path = parsed.path
                ftp = None
                os.makedirs(outdir, exist_ok=True)
                if not os.path.isfile(mapping):
                    print('downloading mapping.csv')
                    ftp = ftp_connect(host, path, username, password,
                                      port=self.config[dem_type]['port'])
                    with open(mapping, 'wb') as myfile:
                        ftp.retrbinary('RETR mapping.csv', myfile.write)
                print('searching FTP server')
                if ftp is None:
                    ftp = ftp_connect(host, path, username, password,
                                      port=self.config[dem_type]['port'])
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
        
        elif dem_type == 'Copernicus 30m Global DEM':
            remotes = cop_dem_remotes(extent=extent, arcsecs=10, product=product)
        
        elif dem_type == 'Copernicus 90m Global DEM':
            remotes = cop_dem_remotes(extent=extent, arcsecs=30, product=product)
        
        else:
            raise ValueError('unknown demType: {}'.format(dem_type))
        
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
        - PROJ: requires ``PROJ_DATA`` or ``PROJ_LIB`` environment variable to be set as download directory; URLs:
        
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
        lookup = {'EGM96': 'us_nga_egm96_15.tif',
                  'EGM2008': 'us_nga_egm08_25.tif'}
        remote = 'https://cdn.proj.org/' + lookup[geoid]
        
        # starting with PROJ 9.1, the PROJ_DATA variable is used.
        # Earlier versions make use of PROJ_LIB.
        var = 'PROJ_DATA'
        proj_dir = os.environ.get(var)
        if proj_dir is None:
            var = 'PROJ_LIB'
            proj_dir = os.environ.get(var)
        if proj_dir is not None:
            local = os.path.join(proj_dir, os.path.basename(remote))
            if not os.path.isfile(local):
                if not os.access(proj_dir, os.W_OK):
                    raise OSError("cannot write to '{0}' path: {1}".format(var, proj_dir))
                log.info('{} <<-- {}'.format(local, remote))
                r = requests.get(remote)
                r.raise_for_status()
                with open(local, 'wb') as out:
                    out.write(r.content)
                r.close()
        else:
            raise RuntimeError("Neither environment variable 'PROJ_DATA' nor 'PROJ_LIB' are set")
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


def vrt_check_sources(fname):
    """
    check the sanity of all source files of a given VRT.
    Currently does not check in-memory VRTs.
    
    Parameters
    ----------
    fname: str
        the VRT file name

    Returns
    -------
    
    Raises
    ------
    RuntimeError
    """
    if os.path.isfile(fname):
        tree = etree.parse(fname)
        sources = [x.text for x in tree.findall('.//SourceFilename')]
        for source in sources:
            if not os.path.isabs(source):
                base_dir = os.path.dirname(fname)
                source = os.path.normpath(os.path.join(base_dir, source))
            if not os.path.isfile(source):
                raise RuntimeError(f'missing VRT source file: {source}')
