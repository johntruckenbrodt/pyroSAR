###############################################################################
# tools for handling auxiliary data in software pyroSAR

# Copyright (c) 2019-2026, the pyroSAR Developers.

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
import socket
import json
import numpy as np
import fnmatch
import ftplib
import requests
import psutil
import zipfile as zf
from lxml import etree
from math import ceil, floor
from urllib.parse import urlparse
from collections import defaultdict
from packaging.version import Version
from pyroSAR.examine import ExamineSnap
from pyroSAR.ancillary import Lock
from spatialist.raster import Raster, Dtype
from spatialist.vector import bbox, Vector
from spatialist.ancillary import finder
from spatialist.auxil import gdalbuildvrt, crsConvert, gdalwarp
from spatialist.envi import HDRobject
from osgeo import gdal, osr

from typing import TypeAlias, Self, Any

import logging

log = logging.getLogger(__name__)

# typing
CRS: TypeAlias = int | str | osr.SpatialReference
EXT: TypeAlias = dict[str, int | float]


def dem_autoload(
        geometries: list[Vector] | None,
        demType: str,
        vrt: str | None = None,
        buffer: int | float | None = None,
        username: str | None = None,
        password: str | None = None,
        product: str = 'dem',
        crop: bool = True,
        lock_timeout: int = 600,
        offline: bool = False,
        return_fname: bool = True
) -> list[str] | None:
    """
    Obtain all relevant DEM tiles for selected geometries.
    The tiles are optionally mosaicked into a VRT file.
    Otherwise, the tiles are returned as file path list.
    With `return_fname` it can be controlled whether the zip file path
    or a path inside the zip including GDAL VSI directive is returned.

    Parameters
    ----------
    geometries
        a list of :class:`spatialist.vector.Vector` geometries to obtain DEM data for;
        CRS must be WGS84 LatLon (EPSG 4326). Can be set to `None` for global extent.
    demType
        the type of DEM to be used. Options:

        - 'AW3D30' (ALOS Global Digital Surface Model "ALOS World 3D - 30m")

          * info: https://www.eorc.jaxa.jp/ALOS/en/aw3d30/index.htm
          * url: ftp://ftp.eorc.jaxa.jp/pub/ALOS/ext1/AW3D30/release_v1804
          * height reference: EGM96

        - 'Copernicus 10m EEA DEM' (Copernicus 10 m DEM available over EEA-39 countries)

          * registration: https://spacedata.copernicus.eu/web/cscda/data-access/registration
          * url: ftps://cdsdata.copernicus.eu/DEM-datasets/COP-DEM_EEA-10-DGED/2021_1
          * height reference: EGM2008

        - 'Copernicus 30m Global DEM'
          
          * info: https://registry.opendata.aws/copernicus-dem
          * url: https://copernicus-dem-30m-stac.s3.amazonaws.com
          * height reference: EGM2008

        - 'Copernicus 30m Global DEM II'
        
          * registration: https://spacedata.copernicus.eu/web/cscda/data-access/registration
          * url: ftps://cdsdata.copernicus.eu/DEM-datasets/COP-DEM_GLO-30-DGED/2021_1
          * height reference: EGM2008
        
        - 'Copernicus 90m Global DEM'
     
          * info: https://registry.opendata.aws/copernicus-dem
          * url: https://copernicus-dem-90m-stac.s3.amazonaws.com
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

          * url: https://step.esa.int/auxdata/dem/SRTM90/tiff
          * height reference: EGM96

    vrt
        an optional GDAL VRT file created from the obtained DEM tiles.
        Setting this to `None` lets the function return the file paths.
        
        .. NOTE::
            VRTs do not support antimeridian splitting.
            For this, return the file paths and pass them to :func:`dem_create` directly.
        
        .. NOTE::
            If the geometries are entirely over ocean it might be possible
            that no DEM tiles are found. In this case a sidecar file is created,
            which covers the geometries' extent and contains one pixel with
            a value of 0 (Naming: ``vrt.replace('.vrt', '_tmp.tif')``).
            This file is then used as source in the VRT.
    
    buffer
        a buffer in degrees to add around the individual geometries
    username
        (optional) the username for services requiring registration
    password
        (optional) the password for the registration account
    product
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
    
    crop
        crop to the provided geometries (or return the full extent of the DEM tiles)?
        Only applies if `vrt!=None`.
    lock_timeout
        how long to wait to acquire a lock on the downloaded files?
    offline
        work offline? If `True`, only locally existing files are considered
        and no online check is performed. If a file is missing, an error is
        raised. For this to work, the function needs to be run in `online`
        mode once to create a local index.
    return_fname: bool
        return the file name including GDAL VSI directive
        (or just the path to the downloaded product)?
        E.g. `/vsizip/srtm_72_02.zip/srtm_72_02.tif` vs. `/srtm_72_02.zip`.
        Only applies if `vrt=None`.
    
    Returns
    -------
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
                   geoid_convert=True)
    """
    with DEMHandler(geometries) as handler:
        return handler.load(
            dem_type=demType,
            username=username,
            password=password,
            vrt=vrt,
            buffer=buffer,
            product=product,
            crop=crop,
            lock_timeout=lock_timeout,
            offline=offline,
            return_fname=return_fname
        )


def dem_create(
        geometries: list[Vector] | None,
        demType: str,
        product: str,
        src: str | list[str],
        dst: str,
        t_srs: CRS | None = None,
        tr: tuple[int | float, int | float] | None = None,
        threads: int | str | None = None,
        geoid_convert: bool = True,
        nodata: int | float | str | None = None,
        resampleAlg: str | None = None,
        dtype: str | None = None,
        pbar: bool = False,
        **kwargs
) -> None:
    """
    Create a new DEM GeoTIFF file and optionally convert heights from geoid to ellipsoid.
    This is basically a convenience wrapper around :func:`osgeo.gdal.Warp`
    via :func:`spatialist.auxil.gdalwarp`.

    Parameters
    ----------
    geometries
        a list of :class:`spatialist.vector.Vector` geometries to obtain DEM data for;
        CRS must be WGS84 LatLon (EPSG 4326). Can be set to `None` for global extent.
    demType
        The input DEM product type. See :func:`dem_autoload` for options.
    product
        The input DEM sub-product type. See :func:`dem_autoload` for options.
    src
        the input dataset(s) as returned by :func:`dem_autoload`.
    dst
        the output GeoTIFF file name.
    t_srs
        A target geographic reference system in WKT, EPSG, PROJ4 or OPENGIS format.
        See function :func:`spatialist.auxil.crsConvert()` for details.
        Default ``None``: use the crs of ``src``.
    tr
        the target resolution as (xres, yres).
        Default ``None``: use the resolution of ``src``.
    threads
        the number of threads to use. Possible values:

         - Default ``None``: use the value of ``GDAL_NUM_THREADS`` without modification.
           If ``GDAL_NUM_THREADS`` is None, multi-threading is still turned on and two
           threads are used, one for I/O and one for computation.
         - integer value: temporarily modify ``GDAL_NUM_THREADS`` and reset it once done.
           If 1, multithreading is turned off.
         - ``ALL_CPUS``: special string to use all cores/CPUs of the computer;
           will also temporarily modify ``GDAL_NUM_THREADS``.
    geoid_convert
        Convert geoid heights to WGS84 ellipsoid heights?
        Only applied if the vertical datum isn't already WGS84 and ``product='dem'``.
        The geoid model is inferred from ``demType``.
        See :func:`dem_autoload` for details.
    nodata
        the nodata value of `dst`. Default ``None``: use these values per data type:

        - ``Byte``: 255
        - ``Int16``: -32768
        - ``Float32``: -9999.0

    resampleAlg
        the resampling algorithm to be used. See here for options:
        https://gdal.org/programs/gdalwarp.html#cmdoption-gdalwarp-r.
        Default ``None``: use ``mode`` if the data type of ``src`` is ``Byte`` (for categorical
        value masks) and ``bilinear`` (for DEM and floating point error masks) otherwise.
    dtype
        the data type of ``dst``. Default ``None``: use the value of the source file(s).
        Data type notations of GDAL (e.g. ``Float32``) and numpy (e.g. ``int8``) are supported.
        See :class:`spatialist.raster.Dtype`.
    pbar
        add a progressbar?
    **kwargs
        additional keyword arguments to be passed to :func:`spatialist.auxil.gdalwarp`.
        See :func:`osgeo.gdal.WarpOptions` for options. The following arguments cannot
        be set as they are controlled internally:

        - ``xRes``, ``yRes``: controlled via argument ``tr``
        - ``srcSRS``, ``dstSRS``: controlled via arguments ``t_srs`` and ``geoid_convert``
        - ``srcNodata``: controlled via nodata value of ``src``
        - ``dstNodata``: controlled via argument ``nodata``
        - ``outputType``: controlled via argument ``dtype``
        - ``multithread``: controlled via argument ``threads``
        - ``format``: fixed to ``GTiff``
        - ``targetAlignedPixels``: fixed to ``True``
        - ``warpOptions``: currently used for setting the number of threads.
          Can be exposed if necessary.
    """
    with DEMHandler(geometries=geometries) as handler:
        handler.create(
            dem_type=demType,
            product=product,
            src=src,
            dst=dst,
            t_srs=t_srs,
            tr=tr,
            threads=threads,
            geoid_convert=geoid_convert,
            nodata=nodata,
            resampleAlg=resampleAlg,
            dtype=dtype,
            pbar=pbar,
            **kwargs
        )


class DEMHandler:
    """
    An interface to obtain DEM data for selected geometries.
    The files are downloaded into the ESA SNAP auxiliary data directory structure.
    This class is the foundation for the convenience functions
    :func:`~pyroSAR.auxdata.dem_autoload` and :func:`~pyroSAR.auxdata.dem_create`.
    
    Parameters
    ----------
    geometries
        a list of geometries. Default `None`: use the global extent.
    """
    
    def __init__(self, geometries: list[Vector] | None = None) -> None:
        if not (isinstance(geometries, list) or geometries is None):
            raise RuntimeError('geometries must be of type list')
        
        if geometries is not None:
            for geometry in geometries:
                if geometry.getProjection('epsg') != 4326:
                    raise RuntimeError('input geometry CRS must be WGS84 LatLon (EPSG 4326)')
        self.geometries = geometries
        try:
            self.auxdatapath = ExamineSnap().auxdatapath
        except AttributeError:
            self.auxdatapath = os.path.join(os.path.expanduser('~'), '.snap', 'auxdata')
    
    def __enter__(self) -> Self:
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        return
    
    @staticmethod
    def __applybuffer(extent: EXT, buffer: int | float | None) -> EXT:
        ext = dict(extent)
        if buffer is not None:
            ext['xmin'] = max(-180., ext['xmin'] - buffer)
            ext['xmax'] = min(180., ext['xmax'] + buffer)
            ext['ymin'] = max(-90., ext['ymin'] - buffer)
            ext['ymax'] = min(90., ext['ymax'] + buffer)
        return ext
    
    def __find_first(self, dem_type: str, product: str) -> str | None:
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
    
    def __buildvrt(
            self,
            tiles: list[str],
            vrt: str,
            dem_type: str,
            product: str,
            crop: bool
    ) -> None:
        """
        Build a VRT mosaic from DEM tiles.
        The VRT is cropped to the specified `extent` but the pixel grid
        of the source files is preserved and no resampling/shifting is applied.
        
        Parameters
        ----------
        tiles
            A list of DEM files or compressed archives containing DEM files.
        vrt
            The output VRT filename.
        dem_type
            The input DEM product type. See :func:`dem_autoload` for options.
        product
            The input DEM sub-product type. See :func:`dem_autoload` for options.
        crop
            Crop the VRT to the spatial extent of the provided geometries
            or return the full extent of the DEM tiles? In the latter case, the common
            bounding box of the geometries is expanded so that the coordinates are
            multiples of the tile size of the respective DEM option.
        """
        
        resolution = None
        dst_datatype = None
        dst_nodata = self.config[dem_type]['ocean_fill_value'][product]
        tap = False
        extent = self.__commonextent()
        if extent['xmin'] > extent['xmax']:
            raise RuntimeError('The extent is crossing the antimeridian. '
                               'It is not possible to create a VRT in this case.')
        
        aop = self.config[dem_type]['area_or_point']
        res = self.__get_resolution(dem_type=dem_type, y=extent['ymin'])
        
        # expand the extent to multiples of the DEM tile size
        if not crop:
            f = self.config[dem_type]['tilesize']
            extent['xmin'] = floor(extent['xmin'] / f) * f
            extent['ymin'] = floor(extent['ymin'] / f) * f
            extent['xmax'] = ceil(extent['xmax'] / f) * f
            extent['ymax'] = ceil(extent['ymax'] / f) * f
        
        # shift coordinates from upper left corner (area) to center (point)
        if aop == 'point':
            shift_x = res[0] / 2
            shift_y = res[1] / 2
            extent['xmin'] -= shift_x
            extent['ymin'] += shift_y
            extent['xmax'] -= shift_x
            extent['ymax'] += shift_y
        
        # special case where no DEM tiles were found because the AOI is completely over ocean
        if len(tiles) == 0:
            # define a dummy file as source file
            # this file contains one pixel with a value of 0
            # nodata value is 255
            tif = vrt.replace('.vrt', '_tmp.tif')
            self.__create_dummy_dem(filename=tif, extent=extent, fill_value=dst_nodata)
            tiles = [tif]
            dst_datatype = self.config[dem_type]['datatype'][product]
            # determine the target resolution based on minimum latitude
            resolution = self.__get_resolution(dem_type=dem_type, y=extent['ymin'])
        
        src_nodata = self.config[dem_type]['nodata'][product]
        
        with Raster(tiles[0]) as ras:
            if src_nodata is None:
                src_nodata = ras.nodata
            if resolution is None:
                xres, yres = ras.res
            else:
                xres, yres = resolution
        opts = {
            'srcNodata': src_nodata,
            'targetAlignedPixels': tap,
            'xRes': xres, 'yRes': yres,
            'hideNodata': True,
            'VRTNodata': dst_nodata,
            'outputBounds': (extent['xmin'], extent['ymin'],
                             extent['xmax'], extent['ymax'])
        }
        
        gdalbuildvrt(src=tiles, dst=vrt, **opts)
        
        if dst_datatype is not None:
            dst_datatype = Dtype(dst_datatype).gdalstr
            tree = etree.parse(source=vrt)
            band = tree.find(path='VRTRasterBand')
            band.attrib['dataType'] = dst_datatype
            tree.write(file=vrt, pretty_print=True,
                       xml_declaration=False, encoding='utf-8')
    
    def __commonextent(self, buffer: int | float | None = None) -> EXT:
        """
        
        Parameters
        ----------
        buffer
            a buffer to add to the common extent

        Returns
        -------
            the common extent of all geometries
        """
        if self.geometries is None:
            return self.__extent_global
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
    def __create_dummy_dem(
            filename: str | None,
            extent: EXT,
            fill_value: int | float
    ) -> gdal.Dataset | list[gdal.Dataset] | None:
        """
        Create dummy dataset(s) which span the given extent and
        are 1x1 pixels large to be as small as possible.
        These dataset(s) are used to create dummy DEMs over ocean.
        
        Parameters
        ----------
        filename
            The filename to be used if the dataset is to be written to a GeoTIFF file.
            If `None`, a :class:`gdal.Dataset` object is returned.
            If the extent is crossing the antimeridian, `filename` must be set
            to `None` because two datasets will be created.
        extent
            The extent for which to create the dataset(s).
        fill_value
            The value to use for extrapolating areas over ocean where no DEM tile exists.
        """
        
        def create_file(
                filename: str | None,
                extent: EXT,
                fill_value: int | float
        ) -> gdal.Dataset | None:
            if filename is None:
                filename = ''
                driver_name = 'MEM' if Version(gdal.__version__) >= Version('3.11') else 'Memory'
            else:
                driver_name = 'GTiff'
            driver = gdal.GetDriverByName(driver_name)
            dataset = driver.Create(utf8_path=filename, xsize=1, ysize=1,
                                    bands=1, eType=gdal.GDT_Byte)
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
            band.SetNoDataValue(value=255)
            arr = np.full(shape=(1, 1), fill_value=fill_value, dtype=np.uint8)
            band.WriteArray(array=arr, xoff=0, yoff=0)
            band.FlushCache()
            del arr
            band = None
            driver = None
            if len(filename) > 0:
                dataset = None
                return None
            else:
                return dataset
        
        # create two sub-datasets if the extent is crossing the antimeridian
        if extent['xmin'] > extent['xmax']:
            if filename is not None:
                raise ValueError("the extent is crossing the antimeridian."
                                 "In this case 'filename' must be None.")
            extents = [
                {
                    'xmin': extent['xmin'],
                    'xmax': 180.,
                    'ymin': extent['ymin'],
                    'ymax': extent['ymax']
                },
                {
                    'xmin': -180.,
                    'xmax': extent['xmax'],
                    'ymin': extent['ymin'],
                    'ymax': extent['ymax']
                }
            ]
            out = [
                create_file(
                    filename=filename,
                    extent=extent_sub,
                    fill_value=fill_value)
                for extent_sub in extents
            ]
        else:
            out = create_file(
                filename=filename,
                extent=extent,
                fill_value=fill_value)
        return out
    
    @property
    def __extent_global(self) -> EXT:
        return {'xmin': -180, 'xmax': 180, 'ymin': -90, 'ymax': 90}
    
    @staticmethod
    def intrange(extent: EXT, step: int) -> tuple[list[int], list[int]]:
        """
        generate a sequence of integer coordinates marking
        the tie points of the individual DEM tiles.
        
        Parameters
        ----------
        extent
            a dictionary with keys `xmin`, `xmax`, `ymin` and `ymax`
            with coordinates in EPSG:4326.
        step
            the sequence steps

        Returns
        -------
            the integer sequences as (latitude, longitude)
        """
        lat = list(range(floor(float(extent['ymin']) / step) * step,
                         ceil(float(extent['ymax']) / step) * step,
                         step))
        if extent['xmin'] > extent['xmax']:
            lon1 = range(floor(float(extent['xmin']) / step) * step,
                         ceil(180. / step) * step,
                         step)
            lon2 = range(floor(-180. / step) * step,
                         ceil(float(extent['xmax']) / step) * step,
                         step)
            lon = list(lon1) + list(lon2)
        else:
            lon = list(range(floor(float(extent['xmin']) / step) * step,
                             ceil(float(extent['xmax']) / step) * step,
                             step))
        return lat, lon
    
    def __get_resolution(
            self,
            dem_type: str, y: int | float
    ) -> tuple[float, float]:
        """
        
        Parameters
        ----------
        dem_type
            the DEM type
        y
            the latitude for which to get the resolution

        Returns
        -------
            (xres, yres)
        """
        for key, val in self.config[dem_type]['resolution'].items():
            ymin, ymax = [int(yr) for yr in key.split('-')]
            if ymin <= abs(y) <= ymax:
                return val
        raise RuntimeError(f"could not get resolution for DEM type "
                           f"'{dem_type}' and latitude '{y}'.")
    
    def __local_index(self, dem_type: str) -> dict[str, dict[str, dict[str, str]]]:
        path = os.path.join(self.auxdatapath, 'dem', dem_type, 'index.json')
        os.makedirs(os.path.dirname(path), exist_ok=True)
        if not os.path.isfile(path):
            with Lock(str(path)):
                if dem_type in ['Copernicus 30m Global DEM',
                                'Copernicus 90m Global DEM']:
                    log.debug(f"building local index for DEM type '{dem_type}'")
                    res = re.search('[39]0', dem_type).group()
                    catalog_json = f"dem_cop_{res}.json"
                    URL_STAC = self.config[dem_type]['url']
                    marker = None
                    out = defaultdict(lambda: defaultdict(dict[str, str]))
                    while True:
                        params = {}
                        if marker:
                            params["marker"] = marker
                        r = requests.get(URL_STAC, params=params)
                        root = etree.fromstring(r.content)
                        is_truncated = root.find(path="./IsTruncated",
                                                 namespaces=root.nsmap).text == "true"
                        items = [x.text for x in root.findall(path="./Contents/Key",
                                                              namespaces=root.nsmap)]
                        if marker is None:
                            del items[items.index(catalog_json)]
                        marker = items[-1]
                        items = sorted([URL_STAC + '/' + x for x in items if x is not None])
                        URL = None
                        for item in items:
                            if URL is None:
                                content = requests.get(item).json()
                                href = content['assets']['elevation']['href']
                                URL = 'https://' + urlparse(href).netloc
                            base = os.path.basename(item).replace('.json', '')
                            lat = re.search('[NS][0-9]{2}', base).group()
                            lon = re.search('[EW][0-9]{3}', base).group()
                            prefix = f"{URL}/{base}_DEM"
                            sub = {
                                "dem": f"{prefix}/{base}_DEM.tif",
                                "edm": f"{prefix}/AUXFILES/{base}_EDM.tif",
                                "flm": f"{prefix}/AUXFILES/{base}_FLM.tif",
                                "wbm": f"{prefix}/AUXFILES/{base}_WBM.tif",
                                "hem": f"{prefix}/AUXFILES/{base}_HEM.tif"
                            }
                            out[lat][lon] = sub
                        if not is_truncated:
                            break
                elif dem_type in ['GETASSE30', 'SRTM 1Sec HGT', 'SRTM 3Sec']:
                    url = self.config[dem_type]['url']
                    response = requests.get(url)
                    response.raise_for_status()
                    items = re.findall(r'href="([^"]+)"', response.text)
                    out = defaultdict(lambda: defaultdict(dict[str, str]))
                    patterns = {
                        'GETASSE30': '(?P<lat>[0-9]{2}[NS])(?P<lon>[0-9]{3}[EW])',
                        'SRTM 1Sec HGT': '(?P<lat>[NS][0-9]{2})(?P<lon>[EW][0-9]{3})',
                        'SRTM 3Sec': '(?P<lon>[0-9]{2})_(?P<lat>[0-9]{2})'
                    }
                    for item in items:
                        if item == '../':
                            continue
                        link = url.rstrip('/') + '/' + item
                        coord = re.search(patterns[dem_type], item).groupdict()
                        out[coord['lat']][coord['lon']] = {'dem': link}
                else:
                    raise RuntimeError(f"local indexing is not supported "
                                       f"for DEM type {dem_type}")
                with open(path, 'w') as f:
                    json.dump(out, f, indent=4)
        with Lock(str(path), soft=True):
            with open(path, 'r') as f:
                index = json.load(f)
        return index
    
    @staticmethod
    def __retrieve(
            urls: list[str],
            outdir: str,
            offline: bool = False,
            lock_timeout: int = 600
    ) -> list[str]:
        if len(urls) == 0:
            return []
        # check that base URL is reachable
        if not offline:
            url_parse = urlparse(urls[0])
            url_base = url_parse.scheme + '://' + url_parse.netloc
            r = requests.get(url_base)
            r.raise_for_status()
            r.close()
        
        urls = list(set(urls))
        os.makedirs(outdir, exist_ok=True)
        locals = []
        n = len(urls)
        for i, remote in enumerate(urls):
            local = os.path.join(outdir, os.path.basename(remote))
            if not os.path.isfile(local):
                if offline:
                    raise RuntimeError(f'file not found locally: {local}')
                else:
                    with Lock(local, timeout=lock_timeout):
                        r = requests.get(remote)
                        # a tile might not exist over the ocean
                        if r.status_code == 404:
                            r.close()
                            continue
                        msg = '[{i: >{w}}/{n}] {l} <<-- {r}'
                        log.info(msg.format(i=i + 1, w=len(str(n)),
                                            n=n, l=local, r=remote))
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
    def __retrieve_ftp(
            url: str,
            filenames: list[str],
            outdir: str,
            username: str | None,
            password: str | None,
            port: int = 0,
            offline: bool = False,
            lock_timeout: int = 600
    ) -> list[str]:
        files = list(set(filenames))
        os.makedirs(outdir, exist_ok=True)
        
        parsed = urlparse(url)
        timeout = 100
        if not offline:
            if parsed.scheme == 'ftpes':
                if username is None or password is None:
                    raise ValueError('Either username or password are set to None')
                ftp = ftplib.FTP_TLS(host=parsed.netloc, timeout=timeout)
                try:
                    ftp.login(username, password)  # login anonymously before securing control channel
                except ftplib.error_perm as e:
                    raise RuntimeError(str(e))
                ftp.prot_p()  # switch to secure data connection.. IMPORTANT! Otherwise, only the user and password is encrypted and not all the file data.
            elif parsed.scheme == 'ftps':
                if username is None or password is None:
                    raise ValueError('Either username or password are set to None')
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
    def config(self) -> dict[str, Any]:
        """
        Get DEM configuration options.
        """
        cop_dem_constants = {
            'area_or_point': 'point',
            'datatype': {
                'dem': 'Float32',
                'edm': 'Byte',
                'flm': 'Byte',
                'hem': 'Float32',
                'wbm': 'Byte'
            },
            'nodata': {
                'dem': None,
                'edm': 0,  # 0 = "Void (no data)"
                'flm': 0,  # 0 = "Void (no data)"
                'hem': -32767.0,
                'wbm': None
            },
            'ocean_fill_value': {
                'dem': 0,
                'edm': 8,  # editing mask; 8 = "Ocean Pixels"
                'flm': 1,  # filling mask; 1 = "Edited (except filled pixels)"
                'hem': -32767.0,  # height error mask; -32767 = Nodata
                'wbm': 1  # water body mask; 1 = "Ocean"
            },
            'pattern': {
                'dem': '*DEM.tif',
                'edm': '*EDM.tif',
                'flm': '*FLM.tif',
                'hem': '*HEM.tif',
                'wbm': '*WBM.tif'
            },
            'resolution_10m': {
                '0-50': (1 / 9000, 1 / 9000),
                '50-60': (1.5 / 9000, 1 / 9000),
                '60-70': (2 / 9000, 1 / 9000),
                '70-80': (3 / 9000, 1 / 9000),
                '80-85': (5 / 9000, 1 / 9000),
                '85-90': (10 / 9000, 1 / 9000)
            },
            'resolution_30m': {
                '0-50': (1 / 3600, 1 / 3600),
                '50-60': (1.5 / 3600, 1 / 3600),
                '60-70': (2 / 3600, 1 / 3600),
                '70-80': (3 / 3600, 1 / 3600),
                '80-85': (5 / 3600, 1 / 3600),
                '85-90': (10 / 3600, 1 / 3600)
            },
            'resolution_90m': {
                '0-50': (1 / 1200, 1 / 1200),
                '50-60': (1.5 / 1200, 1 / 1200),
                '60-70': (2 / 1200, 1 / 1200),
                '70-80': (3 / 1200, 1 / 1200),
                '80-85': (5 / 1200, 1 / 1200),
                '85-90': (10 / 1200, 1 / 1200)
            },
            'tilesize': 1,
            'vertical_datum': 'EGM2008'
        }
        
        return {
            'AW3D30': {
                'url': 'ftp://ftp.eorc.jaxa.jp/pub/ALOS/ext1/AW3D30/release_v1804',
                'ocean_fill_value': {
                    'dem': 0,
                    'msk': 3,  # mask file; 3 = 00000011 binary = "Sea mask"
                    'stk': 0  # stacking number file; 0 = no input files
                },
                'nodata': {
                    'dem': -9999,
                    'msk': None,
                    'stk': None
                },
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
                'authentication': False,
                'vertical_datum': 'EGM96',
            },
            'Copernicus 10m EEA DEM': {
                'url': 'ftps://cdsdata.copernicus.eu/DEM-datasets/COP-DEM_EEA-10-DGED/2021_1',
                'ocean_fill_value': cop_dem_constants['ocean_fill_value'],
                'nodata': cop_dem_constants['nodata'],
                'resolution': cop_dem_constants['resolution_10m'],
                'tilesize': cop_dem_constants['tilesize'],
                'area_or_point': cop_dem_constants['area_or_point'],
                'vsi': '/vsitar/',
                'port': 990,
                'pattern': cop_dem_constants['pattern'],
                'datatype': cop_dem_constants['datatype'],
                'authentication': True,
                'vertical_datum': cop_dem_constants['vertical_datum'],
            },
            'Copernicus 30m Global DEM': {
                'url': 'https://copernicus-dem-30m-stac.s3.amazonaws.com',
                'ocean_fill_value': cop_dem_constants['ocean_fill_value'],
                'nodata': cop_dem_constants['nodata'],
                'resolution': cop_dem_constants['resolution_30m'],
                'tilesize': cop_dem_constants['tilesize'],
                'area_or_point': cop_dem_constants['area_or_point'],
                'vsi': None,
                'pattern': cop_dem_constants['pattern'],
                'datatype': cop_dem_constants['datatype'],
                'authentication': False,
                'vertical_datum': cop_dem_constants['vertical_datum'],
            },
            'Copernicus 30m Global DEM II': {
                'url': 'ftps://cdsdata.copernicus.eu/DEM-datasets/COP-DEM_GLO-30-DGED/2021_1',
                'ocean_fill_value': cop_dem_constants['ocean_fill_value'],
                'nodata': cop_dem_constants['nodata'],
                'resolution': cop_dem_constants['resolution_30m'],
                'tilesize': cop_dem_constants['tilesize'],
                'area_or_point': cop_dem_constants['area_or_point'],
                'vsi': '/vsitar/',
                'port': 990,
                'pattern': cop_dem_constants['pattern'],
                'datatype': cop_dem_constants['datatype'],
                'authentication': True,
                'vertical_datum': cop_dem_constants['vertical_datum'],
            },
            'Copernicus 90m Global DEM': {
                'url': 'https://copernicus-dem-90m-stac.s3.amazonaws.com',
                'ocean_fill_value': cop_dem_constants['ocean_fill_value'],
                'nodata': cop_dem_constants['nodata'],
                'resolution': cop_dem_constants['resolution_90m'],
                'tilesize': cop_dem_constants['tilesize'],
                'area_or_point': cop_dem_constants['area_or_point'],
                'vsi': None,
                'pattern': cop_dem_constants['pattern'],
                'datatype': cop_dem_constants['datatype'],
                'authentication': False,
                'vertical_datum': cop_dem_constants['vertical_datum'],
            },
            'Copernicus 90m Global DEM II': {
                'url': 'ftps://cdsdata.copernicus.eu/DEM-datasets/COP-DEM_GLO-90-DGED/2021_1',
                'ocean_fill_value': cop_dem_constants['ocean_fill_value'],
                'nodata': cop_dem_constants['nodata'],
                'resolution': cop_dem_constants['resolution_90m'],
                'tilesize': cop_dem_constants['tilesize'],
                'area_or_point': cop_dem_constants['area_or_point'],
                'vsi': '/vsitar/',
                'port': 990,
                'pattern': cop_dem_constants['pattern'],
                'datatype': cop_dem_constants['datatype'],
                'authentication': True,
                'vertical_datum': cop_dem_constants['vertical_datum'],
            },
            'GETASSE30': {
                'url': 'https://step.esa.int/auxdata/dem/GETASSE30',
                'ocean_fill_value': {'dem': 0},
                'nodata': {'dem': None},
                'resolution': {'0-90': (15 / 1800, 15 / 1800)},
                'tilesize': 15,
                'area_or_point': 'area',
                'vsi': '/vsizip/',
                'pattern': {'dem': '*.GETASSE30'},
                'datatype': {'dem': 'Int16'},
                'authentication': False,
                'vertical_datum': 'WGS84',
            },
            'SRTM 1Sec HGT': {
                'url': 'https://step.esa.int/auxdata/dem/SRTMGL1',
                'ocean_fill_value': {'dem': 0},
                'nodata': {'dem': None},
                'resolution': {'0-90': (1 / 3600, 1 / 3600)},
                'tilesize': 1,
                'area_or_point': 'point',
                'vsi': '/vsizip/',
                'pattern': {'dem': '*.hgt'},
                'datatype': {'dem': 'Int16'},
                'authentication': False,
                'vertical_datum': 'EGM96',
            },
            'SRTM 3Sec': {
                'url': 'https://step.esa.int/auxdata/dem/SRTM90/tiff',
                'ocean_fill_value': {'dem': 0},
                'nodata': {'dem': None},
                'resolution': {'0-90': (5 / 6000, 5 / 6000)},
                'tilesize': 5,
                'area_or_point': 'area',
                'vsi': '/vsizip/',
                'pattern': {'dem': 'srtm*.tif'},
                'datatype': {'dem': 'Int16'},
                'authentication': False,
                'vertical_datum': 'EGM96',
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
    
    def create(
            self,
            dem_type: str,
            product: str,
            src: str | list[str],
            dst: str,
            t_srs: CRS | None = None,
            tr: tuple[int | float, int | float] | None = None,
            threads: int | str | None = None,
            geoid_convert: bool = True,
            nodata: int | float | str | None = None,
            resampleAlg: str | None = None,
            dtype: str | None = None,
            pbar: bool = False,
            **kwargs
    ) -> None:
        """
        Create a new DEM GeoTIFF file and optionally convert heights from geoid to ellipsoid.
        This is basically a convenience wrapper around :func:`osgeo.gdal.Warp`
        via :func:`spatialist.auxil.gdalwarp`.

        Parameters
        ----------
        dem_type
            The input DEM product type. See :func:`dem_autoload` for options.
        product
            The input DEM sub-product type. See :func:`dem_autoload` for options.
        src
            the input dataset(s) as returned by :func:`dem_autoload`.
            A string is expected to point to a VRT file.
            A list is interpreted as 0..n GDAL-readable datasets.
            :func:`dem_autoload` must be run with `return_fnames=True` to provide the correct format.
        dst
            the output GeoTIFF file name.
        t_srs
            A target geographic reference system in WKT, EPSG, PROJ4 or OPENGIS format.
            See function :func:`spatialist.auxil.crsConvert()` for details.
            Default ``None``: use the crs of ``src``.
        tr
            the target resolution as (xres, yres).
            Default ``None``: use the resolution of ``src``.
        threads
            the number of threads to use. Possible values:

             - Default ``None``: use the value of ``GDAL_NUM_THREADS`` without modification.
               If ``GDAL_NUM_THREADS`` is None, multi-threading is still turned on and two
               threads are used, one for I/O and one for computation.
             - integer value: temporarily modify ``GDAL_NUM_THREADS`` and reset it once done.
               If 1, multithreading is turned off.
             - ``ALL_CPUS``: special string to use all cores/CPUs of the computer;
               will also temporarily modify ``GDAL_NUM_THREADS``.
        geoid_convert
            Convert geoid heights to WGS84 ellipsoid heights?
            Only applied if the vertical datum isn't already WGS84 and ``product='dem'``.
            The geoid model is inferred from ``demType``.
            See :func:`dem_autoload` for details.
        nodata
            the nodata value of `dst`. Default ``None``: use these values per data type:

            - ``Byte``: 255
            - ``Int16``: -32768
            - ``Float32``: -9999.0

        resampleAlg
            the resampling algorithm to be used. See here for options:
            https://gdal.org/programs/gdalwarp.html#cmdoption-gdalwarp-r.
            Default ``None``: use ``mode`` if the data type of ``src`` is ``Byte`` (for categorical
            value masks) and ``bilinear`` (for DEM and floating point error masks) otherwise.
        dtype
            the data type of ``dst``. Default ``None``: use the value of the source file(s).
            Data type notations of GDAL (e.g. ``Float32``) and numpy (e.g. ``int8``) are supported.
            See :class:`spatialist.raster.Dtype`.
        pbar
            add a progressbar?
        **kwargs
            additional keyword arguments to be passed to :func:`spatialist.auxil.gdalwarp`.
            See :func:`osgeo.gdal.WarpOptions` for options. The following arguments cannot
            be set as they are controlled internally:

            - ``xRes``, ``yRes``: controlled via argument ``tr``
            - ``srcSRS``, ``dstSRS``: controlled via arguments ``t_srs`` and ``geoid_convert``
            - ``srcNodata``: controlled via nodata value of ``src``
            - ``dstNodata``: controlled via argument ``nodata``
            - ``outputType``: controlled via argument ``dtype``
            - ``multithread``: controlled via argument ``threads``
            - ``format``: fixed to ``GTiff``
            - ``targetAlignedPixels``: fixed to ``True``
            - ``warpOptions``: currently used for setting the number of threads.
              Can be exposed if necessary.
            
            The following arguments are set if they are not defined in `kwargs`:
            
            - ``outputBounds``: set to the common extent of ``geometries``
              projected to ``t_srs`` if ``isinstance(src, list)``.
        """
        
        if isinstance(src, list):
            src: list[str | gdal.Dataset] = src.copy()
        
        if isinstance(src, str):
            if not src.endswith('.vrt'):
                raise RuntimeError("If 'src' is a string, it must be a VRT file.")
            vrt_check_sources(src)
        
        src_nodata = self.config[dem_type]['nodata'][product]
        src_dtype = self.config[dem_type]['datatype'][product]
        vertical_datum = self.config[dem_type]['vertical_datum']
        fill_value = self.config[dem_type]['ocean_fill_value'][product]
        
        # data type handling
        if dtype is not None:
            dtype_obj = Dtype(dtype)
        else:
            dtype_obj = Dtype(src_dtype)
        
        # set dst nodata default if None
        if nodata is None:
            dst_nodata_lookup = {
                'Byte': 255,
                'Float32': -9999,
                'Int16': -32768
            }
            nodata = dst_nodata_lookup[src_dtype]
        
        if product != 'dem' or vertical_datum == 'WGS84':
            geoid_convert = False
        
        if resampleAlg is None:
            resampleAlg = 'mode' if src_dtype == 'Byte' else 'bilinear'
        ############################################################################
        # initial list of gdalwarp parameters
        
        gdalwarp_args = {
            'format': 'GTiff',
            'srcNodata': src_nodata,
            'resampleAlg': resampleAlg,
            'targetAlignedPixels': True,
            'outputType': dtype_obj.gdalint
        }
        ############################################################################
        # If the input is a list of DEM tiles, pass the user-defined extent directly
        # to gdalwarp (VRTs already contain this extent).
        # Also, add in-memory dummy dataset(s) to the file list so that the output layer
        # is extrapolated to areas where no DEM tile exists (over ocean).
        
        if isinstance(src, list):
            extent_4326 = self.__commonextent()
            if t_srs is not None:
                with bbox(coordinates=extent_4326, crs=4326) as vec:
                    vec.reproject(t_srs)
                    extent_out = vec.extent
            else:
                extent_out = extent_4326
            
            if 'outputBounds' not in kwargs.keys():
                gdalwarp_args['outputBounds'] = [extent_out['xmin'], extent_out['ymin'],
                                                 extent_out['xmax'], extent_out['ymax']]
            else:
                gdalwarp_args['outputBounds'] = kwargs['outputBounds']
                del kwargs['outputBounds']
            
            dummy = self.__create_dummy_dem(filename=None, extent=extent_4326,
                                            fill_value=fill_value)
            if isinstance(dummy, list):
                src = dummy + src
            else:
                src.insert(0, dummy)
        else:
            dummy = None
        ############################################################################
        # specify nodata, resolution and CRS
        
        with Raster(src[0] if isinstance(src, list) else src) as ras:
            src_format = ras.format
            if src_format == 'VRT':
                bytes = Dtype(ras.dtype).bytes
                expecteFileSize = ras.bands * ras.rows * ras.cols * bytes
            if nodata is None:
                nodata = ras.nodata
            if tr is None:
                tr = ras.res
            epsg_in = ras.epsg
        
        if t_srs is None:
            epsg_out = epsg_in
        else:
            epsg_out = crsConvert(t_srs, 'epsg')
        
        gdalwarp_args['dstNodata'] = nodata
        gdalwarp_args['xRes'] = tr[0]
        gdalwarp_args['yRes'] = tr[1]
        gdalwarp_args['srcSRS'] = f'EPSG:{epsg_in}'
        gdalwarp_args['dstSRS'] = f'EPSG:{epsg_out}'
        ############################################################################
        # determine number of threads to use
        
        threads_system = gdal.GetConfigOption('GDAL_NUM_THREADS')
        if isinstance(threads, str):
            if threads != 'ALL_CPUS':
                raise ValueError(f"unsupported value for 'threads': '{threads}'")
            else:
                multithread = True
        elif isinstance(threads, int):
            if threads == 1:
                multithread = False
            elif threads > 1:
                multithread = True
            else:
                raise ValueError("if 'threads' is of type int, it must be >= 1")
        elif threads is None:
            multithread = True
        else:
            raise TypeError(f"'threads' must be of type int, str or None. Is: {type(threads)}")
        
        gdalwarp_args['multithread'] = multithread
        gdalwarp_args['warpOptions'] = {"NUM_THREADS": f"{threads}"}
        ############################################################################
        # intermediate loading of VRT into memory to avoid GDAL bug
        
        c1 = (threads not in [1, None]) and (src_format == 'VRT')
        c2 = (Version(gdal.__version__) < Version('3.12.1'))
        if c1 and c2:
            log.info('Using multithreading for VRT warping is erroneous for smaller GDAL Versions. '
                     'See https://github.com/OSGeo/gdal/issues/13464.'
                     'VRT dataset is transformed to a memory dataset prior to warping.')
            # check free memory for TIFF file creation
            memory = psutil.virtual_memory()
            usedMemory = expecteFileSize * 100 / memory.available
            
            if usedMemory > 80:
                log.warning(f"Warning low memory for warping file "
                            f"{expecteFileSize} {memory.available} {usedMemory}")
            
            # Disable multithreaded gdal.Translate (GDAL_NUM_THREADS = None).
            # Prevents erroneous VRT treatment.
            gdal.SetConfigOption('GDAL_NUM_THREADS', None)
            
            driver_name = 'MEM' if Version(gdal.__version__) >= Version('3.11') else 'Memory'
            src = gdal.Translate(destName='', srcDS=src, format=driver_name)
        ############################################################################
        # update CRS with vertical datum
        
        if geoid_convert:
            geoid_epsg = {'EGM96': 5773,
                          'EGM2008': 3855}
            if vertical_datum in geoid_epsg.keys():
                epsg = geoid_epsg[vertical_datum]
                gdalwarp_args['srcSRS'] += f'+{epsg}'
                # the following line is a workaround for older GDAL versions that did not
                # support compound EPSG codes. See https://github.com/OSGeo/gdal/pull/4639.
                if Version(gdal.__version__) < Version('3.4.0'):
                    crs = crsConvert(gdalwarp_args['srcSRS'], 'proj4')
                    gdalwarp_args['srcSRS'] = crs
            else:
                raise RuntimeError('geoid model not (yet) supported')
            try:
                # download the geoid model
                get_egm_lookup(geoid=vertical_datum, software='PROJ')
            except OSError as e:
                errstr = str(e)
                raise RuntimeError(errstr)
        ############################################################################
        # update gdalwarp_args with user-defined kwargs
        
        locked = list(gdalwarp_args.keys())
        for key, val in kwargs.items():
            if key not in locked:
                gdalwarp_args[key] = val
            else:
                msg = f"argument '{key}' cannot be set via kwargs as it is set internally."
                raise RuntimeError(msg)
        ############################################################################
        # run warping
        
        try:
            if not os.path.isfile(dst):
                message = 'creating mosaic'
                crs = gdalwarp_args['dstSRS']
                if crs != 'EPSG:4326':
                    message += f' and reprojecting to {crs}'
                log.info(f'{message}: {dst}')
                gdalwarp(src=src, dst=dst, pbar=pbar, **gdalwarp_args)
            else:
                log.info(f'mosaic already exists: {dst}')
        except Exception:
            if os.path.isfile(dst):
                os.remove(dst)
            raise
        finally:
            src = dummy = None
            gdal.SetConfigOption('GDAL_NUM_THREADS', threads_system)
    
    def load(
            self,
            dem_type: str,
            vrt: str | None = None,
            buffer: int | float | None = None,
            username: str | None = None,
            password: str | None = None,
            product: str = 'dem',
            crop: bool = True,
            lock_timeout: int = 600,
            offline: bool = False,
            return_fname: bool = True
    ) -> list[str] | None:
        """
        Download DEM tiles. The result is either returned in a list of file
        names or combined into a VRT mosaic. The VRT is cropped to the combined
        extent of the geometries, but the pixel grid of the source files is
        preserved and no resampling/shifting is applied.
        
        Parameters
        ----------
        dem_type
            the type fo DEM to be used
        vrt
            an optional GDAL VRT file created from the obtained DEM tiles.
            NOTE: VRTs are not suited for geometries crossing the antimeridian.
        buffer
            a buffer in degrees to add around the individual geometries
        username
            the download account username
        password
            the download account password
        product
            the sub-product to extract from the DEM product. Options:
            
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
        crop
            If a VRT is created, crop it to the spatial extent of the provided geometries
            or return the full extent of the DEM tiles? In the latter case, the common
            bounding box of the geometries is expanded so that the coordinates are
            multiples of the tile size of the respective DEM option.
        lock_timeout
            how long to wait to acquire a lock on the downloaded files?
        offline
            work offline? If `True`, only locally existing files are considered
            and no online check is performed. If a file is missing, an error is
            raised. For this to work, the function needs to be run in `online`
            mode once to create a local index.
        return_fname: bool
            return the file name including GDAL VSI directive (or just the path to the downloaded product)?
            E.g. `/vsizip/srtm_72_02.zip/srtm_72_02.tif` vs. `/srtm_72_02.zip`.
            Only applies if `vrt=None`.
        
        Returns
        -------
            the names of the obtained files or None if a VRT file was defined with `vrt`.
        """
        keys = self.config.keys()
        if dem_type not in keys:
            options = ', '.join(keys)
            raise RuntimeError(f"DEM type '{dem_type}' is not supported.\n  "
                               f"possible options: '{options}'")
        
        products = self.config[dem_type]['pattern'].keys()
        if product not in products:
            options = ', '.join(products)
            raise RuntimeError(f"Product '{product}' is not available "
                               f"for DEM type '{dem_type}'.\n"
                               f"  options: '{options}'")
        
        outdir = os.path.join(self.auxdatapath, 'dem', dem_type)
        
        if self.geometries is not None:
            candidates = []
            for geo in self.geometries:
                corners = self.__applybuffer(extent=geo.get_extent(split_antimeridian=True),
                                             buffer=buffer)
                candidates.extend(self.remote_ids(extent=corners, dem_type=dem_type,
                                                  username=username, password=password,
                                                  product=product))
        else:
            candidates = self.remote_ids(extent=self.__extent_global, dem_type=dem_type,
                                         username=username, password=password,
                                         product=product)
        
        if self.config[dem_type]['url'].startswith('ftp'):
            port = 0
            if 'port' in self.config[dem_type].keys():
                port = self.config[dem_type]['port']
            locals = self.__retrieve_ftp(url=self.config[dem_type]['url'],
                                         filenames=candidates,
                                         outdir=outdir, username=username,
                                         password=password, port=port,
                                         lock_timeout=lock_timeout,
                                         offline=offline)
        else:
            locals = self.__retrieve(urls=candidates, outdir=outdir,
                                     lock_timeout=lock_timeout,
                                     offline=offline)
        
        # make sure all GETASSE30 tiles get an ENVI HDR file so that they are GDAL-readable
        if dem_type == 'GETASSE30':
            for item in locals:
                getasse30_hdr(item)
        
        # get the file paths inside zip/tar archives and prepend a GDAL VSI directive
        vsi = self.config[dem_type]['vsi']
        if vsi is not None:
            tiles = []
            pattern = self.config[dem_type]['pattern'][product]
            for local in locals:
                if not local.endswith('.tif'):
                    for tile in finder(local, [pattern]):
                        tiles.append(vsi + tile)
        else:
            tiles = locals
        
        if vrt is not None:
            self.__buildvrt(
                tiles=tiles,
                vrt=vrt,
                crop=crop,
                dem_type=dem_type,
                product=product
            )
        else:
            return tiles if return_fname else locals
    
    def remote_ids(
            self,
            extent: EXT,
            dem_type: str,
            product: str = 'dem',
            username: str | None = None,
            password: str | None = None
    ) -> list[str]:
        """
        parse the names/URLs of the remote files overlapping with an area of interest

        Parameters
        ----------
        extent
            the extent of the area of interest with keys xmin, xmax, ymin, ymax
            or `None` to not set any spatial filter.
        dem_type
            the type of DEM to be used
        product
            the sub-product to extract from the DEM product. Only needed for DEM options 'Copernicus 30m Global DEM'
            and 'Copernicus 90m Global DEM' and ignored otherwise.
        username
            the download account username
        password
            the download account password

        Returns
        -------
            the sorted names of the remote files
        """
        keys = self.config.keys()
        if dem_type not in keys:
            raise RuntimeError("demType '{}' is not supported\n  "
                               "possible options: '{}'"
                               .format(dem_type, "', '".join(keys)))
        
        def ids(
                x: int | None = None,
                y: int | None = None,
                nx: int = 3,
                ny: int = 3,
                reverse: bool = False
        ) -> tuple[str, str]:
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
        
        def remotes_from_index(
                indices: list[tuple[str, str]],
                product: str | None
        ) -> list[str]:
            lookup = self.__local_index(dem_type=dem_type)
            remotes = []
            for y, x in indices:
                try:
                    if product is None:
                        remotes.append(lookup[y][x])
                    else:
                        remotes.append(lookup[y][x][product])
                except KeyError:
                    pass
            return remotes
        
        if dem_type in ['Copernicus 30m Global DEM',
                        'Copernicus 90m Global DEM',
                        'SRTM 1Sec HGT']:
            lat, lon = self.intrange(extent, step=1)
            indices = [ids(x, y, nx=3, ny=2)
                       for x in lon for y in lat]
            remotes = remotes_from_index(indices, product=product)
        
        elif dem_type == 'GETASSE30':
            lat, lon = self.intrange(extent, step=15)
            indices = [ids(x, y, nx=3, ny=2, reverse=True)
                       for x in lon for y in lat]
            remotes = remotes_from_index(indices, product=product)
        
        elif dem_type == 'TDX90m':
            lat, lon = self.intrange(extent, step=1)
            remotes = []
            for x in lon:
                xr = abs(x) // 10 * 10
                for y in lat:
                    yf, xf = ids(x=x, y=y, nx=3, ny=2)
                    remotes.append('DEM/{y}/{hem}{xr:03d}/TDM1_DEM__30_{y}{x}.zip'
                                   .format(x=xf, xr=xr, y=yf, hem=xf[0]))
        
        elif dem_type == 'AW3D30':
            remotes = []
            lat, lon = self.intrange(extent, step=1)
            for x in lon:
                for y in lat:
                    remotes.append(
                        '{0}{1}/{2}{3}.tar.gz'.format(*ids(x // 5 * 5, y // 5 * 5),
                                                      *ids(x, y)))
        
        elif dem_type == 'SRTM 3Sec':
            lat = range(
                floor((60 - float(extent['ymax'])) / 5) + 1,
                ceil((60 - float(extent['ymin'])) / 5) + 1
            )
            
            if extent['xmin'] > extent['xmax']:
                lon1 = range(
                    floor((float(extent['xmin']) + 180) / 5) + 1,
                    ceil((180. + 180) / 5) + 1
                )
                lon2 = range(
                    floor((-180 + 180) / 5) + 1,
                    ceil((float(extent['xmax']) + 180) / 5) + 1
                )
                lon = list(lon1) + list(lon2)
            else:
                lon = range(
                    floor((float(extent['xmin']) + 180) / 5) + 1,
                    ceil((float(extent['xmax']) + 180) / 5) + 1
                )
            
            indices = [(f'{y:02d}', f'{x:02d}') for x in lon for y in lat]
            remotes = remotes_from_index(indices, product=product)
        
        elif dem_type in ['Copernicus 10m EEA DEM',
                          'Copernicus 30m Global DEM II',
                          'Copernicus 90m Global DEM II']:
            lat, lon = self.intrange(extent, step=1)
            indices = [''.join(ids(x, y, nx=3, ny=2))
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
                        out.append(target.to_str('latin-1').decode('utf-8'))
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
        else:
            raise ValueError('unknown demType: {}'.format(dem_type))
        
        return sorted(remotes)


def getasse30_hdr(fname: str) -> None:
    """
    create an ENVI HDR file for zipped GETASSE30 DEM tiles
    
    Parameters
    ----------
    fname
        the name of the zipped tile
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


def get_dem_options(require_auth: bool | None = None) -> list[str]:
    """
    Get the names of all supported DEM type options.
    
    Parameters
    ----------
    require_auth
        Only return options that do/don't require authentication.
        Default None: return all options.

    Returns
    -------
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


def get_egm_lookup(geoid: str, software: str) -> None:
    """
    Download lookup tables for converting EGM geoid heights to WGS84 ellipsoid heights.
    
    Parameters
    ----------
    geoid
        the geoid model; current options:
        
        - SNAP: 'EGM96'
        - PROJ: 'EGM96', 'EGM2008'
    software
        the software for which to download the EGM lookup
        
        - SNAP: default directory: ``~/.snap/auxdata/dem/egm96``; URL:
        
          * https://step.esa.int/auxdata/dem/egm96/ww15mgh_b.zip
        - PROJ: requires ``PROJ_DATA`` or ``PROJ_LIB`` environment variable to be set as download directory; URLs:
        
          * https://cdn.proj.org/us_nga_egm96_15.tif
          * https://cdn.proj.org/us_nga_egm08_25.tif
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
    
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self._sock: ssl.SSLSocket | None = None
    
    @property
    def sock(self) -> ssl.SSLSocket | None:
        """Return the socket."""
        return self._sock
    
    @sock.setter
    def sock(self, value: socket.socket | ssl.SSLSocket | None):
        """When modifying the socket, ensure that it is ssl wrapped."""
        if value is not None and not isinstance(value, ssl.SSLSocket):
            value = self.context.wrap_socket(value)
        self._sock = value


def vrt_check_sources(fname: str) -> None:
    """
    check the sanity of all source files of a given VRT.
    Currently, does not check in-memory VRTs.
    
    Parameters
    ----------
    fname
        the VRT file name
    
    Raises
    ------
    RuntimeError
    """
    
    def get_archive_path(path: str) -> str:
        """
        Extract the archive path from a GDAL VSI archive path.

        Examples
        --------
        /vsitar/C:\\tmp\\file.tar.gz\\folder\\image.tif
        -> C:\\tmp\\file.tar.gz

        /vsizip//tmp/file.zip/folder/image.tif
        -> /tmp/file.zip
        """
        pattern = (
            r'^/vsi[a-z]+/'  # VSI directive
            r'(?P<archive>.*?'
            r'\.(?:tar(?:\.gz)?|zip))'  # archive extension
        )
        
        match = re.match(pattern, path, flags=re.IGNORECASE)
        if not match:
            raise RuntimeError(f'could not match archive path: {path}')
        return match.group('archive')
    
    if os.path.isfile(fname):
        tree = etree.parse(fname)
        sources = [x.text for x in tree.findall('.//SourceFilename')]
        for source in sources:
            if source is None:
                raise ValueError('encountered None value as source file name')
            if not os.path.isabs(source):
                base_dir = os.path.dirname(fname)
                source = os.path.normpath(os.path.join(base_dir, source))
            if re.search('^/vsi', source):
                source = get_archive_path(source)
            if not os.path.isfile(source):
                raise RuntimeError(f'missing VRT source file: {source}')
