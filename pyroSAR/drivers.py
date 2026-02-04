###############################################################################
# Reading and Organizing system for SAR images
# Copyright (c) 2016-2025, the pyroSAR Developers.

# This file is part of the pyroSAR Project. It is subject to the
# license terms in the LICENSE.txt file found in the top-level
# directory of this distribution and at
# https://github.com/johntruckenbrodt/pyroSAR/blob/master/LICENSE.txt.
# No part of the pyroSAR project, including this file, may be
# copied, modified, propagated, or distributed except according
# to the terms contained in the LICENSE.txt file.
###############################################################################
"""
This is the core module of package pyroSAR.
It contains the drivers for the different SAR image formats and offers
functionality for retrieving metadata, unpacking images, downloading ancillary files like DEMs and
Orbit State Vector files as well as archiving scenes in a database.
The :class:`ID` class and its subclasses allow easy and standardized access to the metadata of
images from different SAR sensors.
"""

import sys
import gc

from builtins import str
from io import BytesIO

import abc
import ast
import csv
import inspect
import math
import os
import re
import shutil
import struct
import operator
import tarfile as tf
import xml.etree.ElementTree as ET
import zipfile as zf
from datetime import datetime, timezone, timedelta
from dateutil.parser import parse as dateparse
from time import strptime, strftime
from statistics import mean, median
from itertools import groupby
from PIL import Image

import progressbar as pb
from osgeo import gdal, osr, ogr
from osgeo.gdalconst import GA_ReadOnly
import numpy as np

from . import S1, patterns
from .config import __LOCAL__
from .ERS import passdb_query, get_resolution_nesz
from .xml_util import getNamespaces

from spatialist import crsConvert, sqlite3, Vector, bbox
from spatialist.ancillary import parse_literal, finder, multicore

from sqlalchemy import create_engine, Table, MetaData, Column, Integer, String, exc
from sqlalchemy import inspect as sql_inspect
from sqlalchemy.event import listen
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql import select, func
from sqlalchemy.engine.url import URL
from sqlalchemy.ext.automap import automap_base
from sqlalchemy_utils import database_exists, create_database, drop_database
from geoalchemy2 import Geometry
import socket
import time
import platform
import logging

log = logging.getLogger(__name__)


def identify(scene):
    """
    identify a SAR scene and return the appropriate metadata handler object

    Parameters
    ----------
    scene: str
        a file or directory name

    Returns
    -------
    pyroSAR.drivers.ID
        a pyroSAR metadata handler
    
    Examples
    --------

    >>> from pyroSAR import identify
    >>> filename = 'S1A_IW_GRDH_1SDV_20180829T170656_20180829T170721_023464_028DE0_F7BD.zip'
    >>> scene = identify(filename)
    >>> print(scene)
    pyroSAR ID object of type SAFE
    acquisition_mode: IW
    cycleNumber: 148
    frameNumber: 167392
    lines: 16703
    orbit: A
    orbitNumber_abs: 23464
    orbitNumber_rel: 117
    polarizations: ['VV', 'VH']
    product: GRD
    projection: +proj=longlat +datum=WGS84 +no_defs
    samples: 26056
    sensor: S1A
    spacing: (10.0, 10.0)
    start: 20180829T170656
    stop: 20180829T170721
    """
    if not os.path.exists(scene):
        raise OSError("No such file or directory: '{}'".format(scene))
    
    def get_subclasses(c):
        subclasses = c.__subclasses__()
        for subclass in subclasses.copy():
            subclasses.extend(get_subclasses(subclass))
        return list(set(subclasses))
    
    for handler in get_subclasses(ID):
        try:
            return handler(scene)
        except Exception:
            pass
    raise RuntimeError('data format not supported')


def identify_many(scenes, pbar=False, sortkey=None, cores=1):
    """
    wrapper function for returning metadata handlers of all valid scenes in a list,
    similar to function :func:`~pyroSAR.drivers.identify`.

    Parameters
    ----------
    scenes: list[str or ID]
        the file names of the scenes to be identified
    pbar: bool
        adds a progressbar if True
    sortkey: str or None
        sort the handler object list by an attribute
    cores: int
        the number of cores to parallelize identification
    
    Returns
    -------
    list[ID]
        a list of pyroSAR metadata handlers
    
    Examples
    --------
    >>> from pyroSAR import identify_many
    >>> files = finder('/path', ['S1*.zip'])
    >>> ids = identify_many(files, pbar=False, sortkey='start')
    """
    
    def handler(scene):
        if isinstance(scene, ID):
            return scene
        else:
            try:
                id = identify(scene)
                return id
            except RuntimeError:
                return None
            except PermissionError:
                log.warning("Permission denied: '{}'".format(scene))
    
    if cores == 1:
        idlist = []
        if pbar:
            progress = pb.ProgressBar(max_value=len(scenes)).start()
        else:
            progress = None
        for i, scene in enumerate(scenes):
            id = handler(scene)
            idlist.append(id)
            if progress is not None:
                progress.update(i + 1)
        if progress is not None:
            progress.finish()
    else:
        idlist = multicore(function=handler, multiargs={'scene': scenes},
                           pbar=pbar, cores=cores)
    if sortkey is not None:
        idlist.sort(key=operator.attrgetter(sortkey))
    idlist = list(filter(None, idlist))
    return idlist


def filter_processed(scenelist, outdir, recursive=False):
    """
    Filter a list of pyroSAR objects to those that have not yet been processed and stored in the defined directory.
    The search for processed scenes is either done in the directory only or recursively into subdirectories.
    The scenes must have been processed with pyroSAR in order to follow the right naming scheme.

    Parameters
    ----------
    scenelist: list[ID]
        a list of pyroSAR objects
    outdir: str
        the processing directory
    recursive: bool
        scan `outdir` recursively into subdirectories?

    Returns
    -------
    list[ID]
        a list of those scenes, which have not been processed yet
    """
    return [x for x in scenelist if not x.is_processed(outdir, recursive)]


class ID(object):
    """
    Abstract class for SAR meta data handlers
    """
    
    def __init__(self, metadict):
        """
        to be called by the __init__methods of the format drivers
        scans a metadata dictionary and registers entries with a standardized name as object attributes
        see __LOCAL__ for standard names. It must be ensured that each of these is actually read by the individual SAR format driver.

        :param metadict: a dictionary containing the metadata attributes of a SAR scene
        """
        self.locals = __LOCAL__
        for item in self.locals:
            setattr(self, item, metadict[item])
    
    def __getattr__(self, item):
        raise AttributeError("object has no attribute '{}'".format(item))
    
    def __str__(self):
        lines = ['pyroSAR ID object of type {}'.format(self.__class__.__name__)]
        for item in sorted(self.locals):
            value = getattr(self, item)
            if item == 'projection':
                value = crsConvert(value, 'proj4') if value is not None else None
            if value == -1:
                value = '<no global value per product>'
            line = '{0}: {1}'.format(item, value)
            lines.append(line)
        return '\n'.join(lines)
    
    def bbox(self, outname=None, driver=None, overwrite=True, buffer=None):
        """
        get the bounding box of a scene. The result is either returned as
        vector object or written to a file.

        Parameters
        ----------
        outname: str
            the name of the vector file to be written
        driver: str
            the output file format; needs to be defined if the format cannot
            be auto-detected from the filename extension
        overwrite: bool
            overwrite an existing vector file?
        buffer: None or int or float or tuple[int or float]
            a buffer to add around `coordinates`. Default None: do not add
            a buffer. A tuple is interpreted as (x buffer, y buffer).

        Returns
        -------
        ~spatialist.vector.Vector or None
            the vector object if `outname` is None and None otherwise
        
        See Also
        --------
        spatialist.vector.Vector.bbox
        """
        if outname is None:
            return bbox(coordinates=self.getCorners(), crs=self.projection,
                        buffer=buffer)
        else:
            bbox(coordinates=self.getCorners(), crs=self.projection,
                 outname=outname, driver=driver, overwrite=overwrite,
                 buffer=buffer)
    
    def geometry(self, outname=None, driver=None, overwrite=True):
        """
        get the footprint geometry of a scene either as a vector object or written to a file

        Parameters
        ----------
        outname: str
            the name of the vector file to be written
        driver: str
            the output file format; needs to be defined if the format cannot
            be auto-detected from the filename extension
        overwrite: bool
            overwrite an existing vector file?

        Returns
        -------
        ~spatialist.vector.Vector or None
            the vector object if `outname` is None, None otherwise
        
        See also
        --------
        spatialist.vector.Vector.write
        """
        if 'coordinates' not in self.meta.keys():
            raise NotImplementedError
        srs = crsConvert(self.projection, 'osr')
        points = ogr.Geometry(ogr.wkbMultiPoint)
        for lon, lat in self.meta['coordinates']:
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(lon, lat)
            points.AddGeometry(point)
        geom = points.ConvexHull()
        geom.FlattenTo2D()
        point = points = None
        exterior = geom.GetGeometryRef(0)
        if exterior.IsClockwise():
            points = list(exterior.GetPoints())
            exterior.Empty()
            for x, y in reversed(points):
                exterior.AddPoint(x, y)
            geom.CloseRings()
        exterior = points = None
        
        bbox = Vector(driver='Memory')
        bbox.addlayer('geometry', srs, geom.GetGeometryType())
        bbox.addfield('area', ogr.OFTReal)
        bbox.addfeature(geom, fields={'area': geom.Area()})
        geom = None
        if outname is None:
            return bbox
        else:
            bbox.write(outfile=outname, driver=driver, overwrite=overwrite)
    
    @property
    def compression(self):
        """
        check whether a scene is compressed into an tarfile or zipfile or not at all

        Returns
        -------
        str or None
            either 'zip', 'tar' or None
        """
        if os.path.isdir(self.scene):
            return None
        elif zf.is_zipfile(self.scene):
            return 'zip'
        elif tf.is_tarfile(self.scene):
            return 'tar'
        else:
            return None
    
    def export2dict(self):
        """
        Return the uuid and the metadata that is defined in `self.locals` as a dictionary
        """
        metadata = {item: self.meta[item] for item in self.locals}
        sq_file = os.path.basename(self.file)
        title = os.path.splitext(sq_file)[0]
        metadata['uuid'] = title
        return metadata
    
    def export2sqlite(self, dbfile):
        """
        Export relevant metadata to an SQLite database

        Parameters
        ----------
        dbfile: str
            the database file

        """
        with Archive(dbfile) as archive:
            archive.insert(self)
    
    def examine(self, include_folders=False):
        """
        check whether any items in the SAR scene structure (i.e. files/folders) match the regular expression pattern
        defined by the class. On success the item is registered in the object as attribute `file`.

        Parameters
        ----------
        include_folders: bool
            also match folder (or just files)?

        Returns
        -------

        Raises
        -------
        RuntimeError
        """
        files = self.findfiles(self.pattern, include_folders=include_folders)
        if len(files) == 1:
            self.file = files[0]
        elif len(files) == 0:
            raise RuntimeError('scene does not match {} naming convention'.format(type(self).__name__))
        else:
            raise RuntimeError('file ambiguity detected:\n{}'.format('\n'.join(files)))
    
    def findfiles(self, pattern, include_folders=False):
        """
        find files in the scene archive, which match a pattern.

        Parameters
        ----------
        pattern: str
            the regular expression to match
        include_folders: bool
             also match folders (or just files)?
        Returns
        -------
        list[str]
            the matched file names
        
        See Also
        --------
        :func:`spatialist.ancillary.finder`
        """
        foldermode = 1 if include_folders else 0
        
        try:
            files = finder(target=self.scene, matchlist=[pattern],
                           foldermode=foldermode, regex=True)
        except RuntimeError:
            # Return the scene if only a file and not zip
            return self.scene
        
        if os.path.isdir(self.scene) \
                and re.search(pattern, os.path.basename(self.scene)) \
                and include_folders:
            files.append(self.scene)
        
        return files
    
    def gdalinfo(self):
        """
        read metadata directly from the GDAL SAR image drivers

        Returns
        -------
        dict
            the metadata attributes
        """
        files = self.findfiles(r'(?:\.[NE][12]$|DAT_01\.001$|product\.xml|manifest\.safe$)')
        # If only one file return the file in array
        if isinstance(files, str):
            files = [files]
        
        if len(files) == 1:
            prefix = {'zip': '/vsizip/', 'tar': '/vsitar/', None: ''}[self.compression]
            header = files[0]
        elif len(files) > 1:
            raise RuntimeError('file ambiguity detected')
        else:
            raise RuntimeError('file type not supported')
        
        meta = {}
        
        ext_lookup = {'.N1': 'ASAR', '.E1': 'ERS1', '.E2': 'ERS2'}
        extension = os.path.splitext(header)[1]
        if extension in ext_lookup:
            meta['sensor'] = ext_lookup[extension]
            info = gdal.Info(prefix + header, options=gdal.InfoOptions(allMetadata=True, format='json'))
            meta['extra'] = info
        
        img = gdal.Open(prefix + header, GA_ReadOnly)
        gdalmeta = img.GetMetadata()
        meta['samples'], meta['lines'], meta['bands'] = img.RasterXSize, img.RasterYSize, img.RasterCount
        meta['projection'] = img.GetGCPProjection()
        meta['gcps'] = [((x.GCPPixel, x.GCPLine), (x.GCPX, x.GCPY, x.GCPZ)) for x in img.GetGCPs()]
        img = None
        
        for item in gdalmeta:
            entry = [item, parse_literal(gdalmeta[item].strip())]
            
            try:
                entry[1] = self.parse_date(str(entry[1]))
            except ValueError:
                pass
            
            if re.search('LAT|LONG', entry[0]):
                entry[1] /= 1000000.
            meta[entry[0]] = entry[1]
        return meta
    
    def getCorners(self):
        """
        Get the bounding box corner coordinates

        Returns
        -------
        dict
            the corner coordinates as a dictionary with keys `xmin`, `ymin`, `xmax`, `ymax`
        """
        if 'coordinates' not in self.meta.keys():
            raise NotImplementedError
        coordinates = self.meta['coordinates']
        lat = [x[1] for x in coordinates]
        lon = [x[0] for x in coordinates]
        return {'xmin': min(lon), 'xmax': max(lon), 'ymin': min(lat), 'ymax': max(lat)}
    
    def getFileObj(self, filename):
        """
        Load a file into a readable file object.

        Parameters
        ----------
        filename: str
            the name of a file in the scene archive, easiest to get with method :meth:`~ID.findfiles`

        Returns
        -------
        io.BytesIO
            a file pointer object
        """
        return getFileObj(self.scene, filename)
    
    def getGammaImages(self, directory=None):
        """
        list all files processed by GAMMA

        Parameters
        ----------
        directory: str or None
            the directory to be scanned; if left empty the object attribute `gammadir` is scanned

        Returns
        -------
        list[str]
            the file names of the images processed by GAMMA

        Raises
        -------
        RuntimeError
        """
        if directory is None:
            if hasattr(self, 'gammadir'):
                directory = self.gammadir
            else:
                raise RuntimeError(
                    'directory missing; please provide directory to function or define object attribute "gammadir"')
        return [x for x in finder(directory, [self.outname_base()], regex=True) if
                not re.search(r'\.(?:par|hdr|aux\.xml|swp|sh)$', x)]
    
    def getHGT(self):
        """
        get the names of all SRTM HGT tiles overlapping with the SAR scene

        Returns
        -------
        list[str]
            names of the SRTM HGT tiles
        """
        
        corners = self.getCorners()
        
        # generate sequence of integer coordinates marking the tie points of the overlapping hgt tiles
        lat = range(int(float(corners['ymin']) // 1), int(float(corners['ymax']) // 1) + 1)
        lon = range(int(float(corners['xmin']) // 1), int(float(corners['xmax']) // 1) + 1)
        
        # convert coordinates to string with leading zeros and hemisphere identification letter
        lat = [str(x).zfill(2 + len(str(x)) - len(str(x).strip('-'))) for x in lat]
        lat = [x.replace('-', 'S') if '-' in x else 'N' + x for x in lat]
        
        lon = [str(x).zfill(3 + len(str(x)) - len(str(x).strip('-'))) for x in lon]
        lon = [x.replace('-', 'W') if '-' in x else 'E' + x for x in lon]
        
        # concatenate all formatted latitudes and longitudes with each other as final product
        return [x + y + '.hgt' for x in lat for y in lon]
    
    def is_processed(self, outdir, recursive=False):
        """
        check whether a scene has already been processed and stored in the defined output directory
        (and subdirectories if scanned recursively)

        Parameters
        ----------
        outdir: str
            the directory to be checked

        Returns
        -------
        bool
            does an image matching the scene pattern exist?
        """
        if os.path.isdir(outdir):
            # '{}.*tif$'.format(self.outname_base())
            return len(finder(outdir, [self.outname_base()], regex=True, recursive=recursive)) != 0
        else:
            return False
    
    def outname_base(self, extensions=None):
        """
        parse a string containing basic information about the scene in standardized format.
        Currently, this id contains the sensor (4 digits), acquisition mode (4 digits), orbit (1 digit)
        and acquisition start time (15 digits)., e.g. `S1A__IW___A_20150523T122350`.
        
        Parameters
        ----------
        extensions: list[str] or None
            the names of additional parameters to append to the basename, e.g. ``['orbitNumber_rel']``
        Returns
        -------
        str
            a standardized name unique to the scene
            
        """
        
        fields = ('{:_<4}'.format(self.sensor),
                  '{:_<4}'.format(self.acquisition_mode),
                  self.orbit,
                  self.start)
        out = '_'.join(fields)
        if isinstance(extensions, list) and len(extensions) is not None:
            ext = '_'.join([str(getattr(self, key)) for key in extensions])
            out += '_' + ext
        return out
    
    @staticmethod
    def parse_date(x):
        """
        this function gathers known time formats provided in the different SAR products and converts them to a common
        standard of the form YYYYMMDDTHHMMSS.

        Parameters
        ----------
        x: str
            the time stamp

        Returns
        -------
        str
            the converted time stamp in format YYYYmmddTHHMMSS
        """
        return parse_date(x)
    
    @abc.abstractmethod
    def quicklook(self, outname, format='kmz'):
        """
        export a quick look image of the scene

        Parameters
        ----------
        outname: str
            the name of the output file
        format: str
            the format of the file to write;
            currently only kmz is supported

        Returns
        -------

        Examples
        --------

        >>> from pyroSAR import identify
        >>> scene = identify('S1A_IW_GRDH_1SDV_20180101T170648_20180101T170713_019964_021FFD_DA78.zip')
        >>> scene.quicklook('S1A__IW___A_20180101T170648.kmz')
        """
        raise NotImplementedError
    
    @property
    def start_dt(self) -> datetime:
        """
        
        Returns
        -------
            the acquisition start time as timezone-aware datetime object
        """
        out = datetime.strptime(self.start, '%Y%m%dT%H%M%S')
        return out.replace(tzinfo=timezone.utc)
    
    @property
    def stop_dt(self) -> datetime:
        """
        
        Returns
        -------
            the acquisition stop time as timezone-aware datetime object
        """
        out = datetime.strptime(self.stop, '%Y%m%dT%H%M%S')
        return out.replace(tzinfo=timezone.utc)
    
    def summary(self):
        """
        print the set of standardized scene metadata attributes

        Returns
        -------

        """
        print(self.__str__())
    
    @abc.abstractmethod
    def scanMetadata(self):
        """
        scan SAR scenes for metadata attributes.
        The returned dictionary is registered as attribute `meta` by the class upon object initialization.
        This dictionary furthermore needs to return a set of standardized attribute keys,
        which are directly registered as object attributes.

        Returns
        -------
        dict
            the derived attributes

        """
        raise NotImplementedError
    
    @abc.abstractmethod
    def unpack(self, directory, overwrite=False, exist_ok=False):
        """
        Unpack the SAR scene into a defined directory.

        Parameters
        ----------
        directory: str
            the base directory into which the scene is unpacked
        overwrite: bool
            overwrite an existing unpacked scene?
        exist_ok: bool
            allow existing output files and do not create new ones?

        Returns
        -------

        """
        raise NotImplementedError
    
    def _unpack(self, directory, offset=None, overwrite=False, exist_ok=False):
        """
        general function for unpacking scene archives; to be called by implementations of ID.unpack.
        Will reset object attributes `scene` and `file` to point to the locations of the unpacked scene
        
        Parameters
        ----------
        directory: str
            the name of the directory in which the files are written
        offset: str
            an archive directory offset; to be defined if only a subdirectory is to be unpacked (see e.g. TSX.unpack)
        overwrite: bool
            should an existing directory be overwritten?
        exist_ok: bool
            do not attempt unpacking if the target directory already exists? Ignored if ``overwrite==True``
        
        Returns
        -------
        
        """
        do_unpack = True
        if os.path.isdir(directory):
            if overwrite:
                shutil.rmtree(directory)
            else:
                if exist_ok:
                    do_unpack = False
                else:
                    raise RuntimeError('target scene directory already exists: {}'.format(directory))
        os.makedirs(directory, exist_ok=True)
        
        if do_unpack:
            if tf.is_tarfile(self.scene):
                archive = tf.open(self.scene, 'r')
                names = archive.getnames()
                if offset is not None:
                    names = [x for x in names if x.startswith(offset)]
                header = os.path.commonprefix(names)
                
                if header in names:
                    if archive.getmember(header).isdir():
                        for item in sorted(names):
                            if item != header:
                                member = archive.getmember(item)
                                if offset is not None:
                                    member.name = member.name.replace(offset + '/', '')
                                archive.extract(member, directory)
                        archive.close()
                    else:
                        archive.extractall(directory)
                        archive.close()
            
            elif zf.is_zipfile(self.scene):
                archive = zf.ZipFile(self.scene, 'r')
                names = archive.namelist()
                header = os.path.commonprefix(names)
                if header.endswith('/'):
                    for item in sorted(names):
                        if item != header:
                            repl = item.replace(header, '', 1)
                            outname = os.path.join(directory, repl)
                            outname = outname.replace('/', os.path.sep)
                            if item.endswith('/'):
                                os.makedirs(outname, exist_ok=True)
                            else:
                                os.makedirs(os.path.dirname(outname), exist_ok=True)
                                try:
                                    with open(outname, 'wb') as outfile:
                                        outfile.write(archive.read(item))
                                except zf.BadZipfile:
                                    log.info('corrupt archive, unpacking failed')
                                    continue
                    archive.close()
                else:
                    archive.extractall(directory)
                    archive.close()
            else:
                log.info('unpacking is only supported for TAR and ZIP archives')
                return
        
        self.scene = directory
        main = os.path.join(self.scene, os.path.basename(self.file))
        self.file = main if os.path.isfile(main) else self.scene


class BEAM_DIMAP(ID):
    """
    Handler class for BEAM-DIMAP data

    Sensors:
        * SNAP supported sensors
    """
    
    def __init__(self, scene):
        
        if not scene.lower().endswith('.dim'):
            raise RuntimeError('Scene format is not BEAM-DIMAP')
        
        self.root = None
        self.scene = scene
        self.meta = self.scanMetadata()
        
        super(BEAM_DIMAP, self).__init__(self.meta)
    
    def scanMetadata(self):
        meta = dict()
        
        self.root = ET.parse(self.scene).getroot()
        
        def get_by_name(attr: list[str] | str, section: str = 'Abstracted_Metadata') -> str:
            msg = 'cannot get attribute "{}" from section "{}"'
            if isinstance(attr, list):
                for i, item in enumerate(attr):
                    try:
                        return get_by_name(item, section=section)
                    except RuntimeError:
                        continue
                raise RuntimeError(msg.format('|'.join(attr), section))
            else:
                element = self.root.find(f'.//MDElem[@name="{section}"]')
                out = element.find(f'.//MDATTR[@name="{attr}"]')
                if out is None or out.text in ['99999', '99999.0']:
                    raise RuntimeError(msg.format(attr, section))
                return out.text
        
        missions = {'ENVISAT': 'ASAR',
                    'ERS1': 'ERS1',
                    'ERS2': 'ERS2',
                    'SENTINEL-1A': 'S1A',
                    'SENTINEL-1B': 'S1B',
                    'SENTINEL-1C': 'S1C',
                    'SENTINEL-1D': 'S1D'}
        
        section = 'Abstracted_Metadata'
        meta['sensor'] = missions[get_by_name('MISSION', section=section)]
        if re.search('S1[A-Z]', meta['sensor']):
            meta['acquisition_mode'] = get_by_name('ACQUISITION_MODE', section=section)
            meta['product'] = self.root.find('.//PRODUCT_TYPE').text
        elif meta['sensor'] in ['ASAR', 'ERS1', 'ERS2']:
            product_type = get_by_name('PRODUCT_TYPE', section=section)
            meta['acquisition_mode'] = product_type[4:7]
            # product overview table: https://doi.org/10.5167/UZH-96146
            if meta['acquisition_mode'] in ['APS', 'IMS', 'WSS']:
                meta['product'] = 'SLC'
            elif meta['acquisition_mode'] in ['APP', 'IMP']:
                meta['product'] = 'PRI'
            elif meta['acquisition_mode'] in ['APM', 'IMM', 'WSM']:
                meta['product'] = 'MR'
            else:
                raise RuntimeError(f"unsupported acquisition mode: '{meta['acquisition_mode']}'")
        else:
            raise RuntimeError('unknown sensor {}'.format(meta['sensor']))
        
        meta['IPF_version'] = get_by_name('Processing_system_identifier', section=section)
        
        meta['orbit'] = get_by_name('PASS', section=section)[0]
        pols = [x.text for x in self.root.findall('.//MDATTR[@desc="Polarization"]')]
        pols = list(filter(None, pols))
        meta['polarizations'] = list(set([x for x in pols if '-' not in x]))
        meta['spacing'] = (round(float(get_by_name('range_spacing', section=section)), 6),
                           round(float(get_by_name('azimuth_spacing', section=section)), 6))
        meta['looks'] = (float(get_by_name('range_looks', section=section)),
                         float(get_by_name('azimuth_looks', section=section)))
        meta['samples'] = int(self.root.find('.//BAND_RASTER_WIDTH').text)
        meta['lines'] = int(self.root.find('.//BAND_RASTER_HEIGHT').text)
        meta['bands'] = int(self.root.find('.//NBANDS').text)
        meta['orbitNumber_abs'] = int(get_by_name('ABS_ORBIT', section=section))
        meta['orbitNumber_rel'] = int(get_by_name('REL_ORBIT', section=section))
        meta['cycleNumber'] = int(get_by_name(['orbit_cycle', 'CYCLE'], section=section))
        meta['frameNumber'] = int(get_by_name(['data_take_id', 'ABS_ORBIT'], section=section))
        
        meta['swath'] = get_by_name('SWATH', section=section)
        
        srgr = bool(int(get_by_name('srgr_flag', section=section)))
        meta['image_geometry'] = 'GROUND_RANGE' if srgr else 'SLANT_RANGE'
        #################################################################################
        # start, stop
        start = datetime.strptime(self.root.find('.//PRODUCT_SCENE_RASTER_START_TIME').text,
                                  '%d-%b-%Y %H:%M:%S.%f')
        meta['start'] = start.strftime('%Y%m%dT%H%M%S')
        stop = datetime.strptime(self.root.find('.//PRODUCT_SCENE_RASTER_STOP_TIME').text,
                                 '%d-%b-%Y %H:%M:%S.%f')
        meta['stop'] = stop.strftime('%Y%m%dT%H%M%S')
        #################################################################################
        # incident angle
        # the incident angle is not stored consistently so several options are tried
        while True:
            # may be missing or set to '99999.0'
            try:
                inc_near = get_by_name('incidence_near', section=section)
                inc_far = get_by_name('incidence_far', section=section)
                incidence = (float(inc_near) + float(inc_far)) / 2
                break
            except RuntimeError:
                pass
            # this attribute might only apply to Sentinel-1
            inc_elements = self.root.findall('.//MDATTR[@name="incidenceAngleMidSwath"]')
            if len(inc_elements) > 0:
                incidence = [float(x.text) for x in inc_elements]
                incidence = mean(incidence)
                break
            # the tie point grids are no longer present in geocoded products
            inc_grid = os.path.join(self.scene.replace('.dim', '.data'),
                                    'tie_point_grids', 'incident_angle.img')
            if os.path.isfile(inc_grid):
                ras = gdal.Open(inc_grid)
                arr = ras.ReadAsArray()
                incidence = np.mean(arr[arr != 0])
                ras = arr = None
                break
            raise ValueError('cannot read the incident angle')
        meta['incidence'] = incidence
        #################################################################################
        # projection
        if self.root.find('.//WKT') is not None:
            meta['projection'] = self.root.find('.//WKT').text.lstrip()
        else:
            meta['projection'] = crsConvert(4326, 'wkt')
        #################################################################################
        # coordinates
        keys = ['{}_{}_{}'.format(a, b, c)
                for a in ['first', 'last']
                for b in ['far', 'near']
                for c in ['lat', 'long']]
        coords = {key: float(get_by_name(key, section=section))
                  for key in keys}
        
        meta['coordinates'] = [(coords['first_near_long'], coords['first_near_lat']),
                               (coords['last_near_long'], coords['last_near_lat']),
                               (coords['last_far_long'], coords['last_far_lat']),
                               (coords['first_far_long'], coords['first_far_lat'])]
        #################################################################################
        return meta
    
    def unpack(self, directory, overwrite=False, exist_ok=False):
        raise RuntimeError('unpacking of BEAM-DIMAP products is not supported')


class CEOS_ERS(ID):
    """
    Handler class for ERS data in CEOS format
    
    Sensors:
        * ERS1
        * ERS2
    
    Reference:
        ER-IS-EPO-GS-5902-3: Annex C. ERS SAR.SLC/SLC-I. CCT and EXABYTE
        (`ESA 1998 <https://earth.esa.int/documents/10174/1597298/SAR05E.pdf>`_)
    """
    
    def __init__(self, scene):
        self.pattern = patterns.ceos_ers
        
        self.pattern_pid = r'(?P<sat_id>(?:SAR|ASA))_' \
                           r'(?P<image_mode>(?:IM(?:S|P|G|M|_)|AP(?:S|P|G|M|_)|WV(?:I|S|W|_)|WS(?:M|S|_)))_' \
                           r'(?P<processing_level>[012B][CP])'
        
        self.scene = os.path.realpath(scene)
        
        self.examine()
        
        self.meta = self.scanMetadata()
        
        # register the standardized meta attributes as object attributes
        super(CEOS_ERS, self).__init__(self.meta)
    
    def unpack(self, directory, overwrite=False, exist_ok=False):
        if self.sensor in ['ERS1', 'ERS2']:
            base_file = re.sub(r'\.PS$', '', os.path.basename(self.file))
            base_dir = os.path.basename(directory.strip('/'))
            
            outdir = directory if base_file == base_dir else os.path.join(directory, base_file)
            
            self._unpack(outdir, overwrite=overwrite, exist_ok=exist_ok)
        else:
            raise NotImplementedError('sensor {} not implemented yet'.format(self.sensor))
    
    def scanMetadata(self):
        meta = dict()
        
        match = re.match(re.compile(self.pattern), os.path.basename(self.file))
        match2 = re.match(re.compile(self.pattern_pid), match.group('product_id'))
        
        if re.search('IM__0', match.group('product_id')):
            raise RuntimeError('product level 0 not supported (yet)')
        
        meta['acquisition_mode'] = match2.group('image_mode')
        meta['product'] = 'SLC' if meta['acquisition_mode'] in ['IMS', 'APS', 'WSS'] else 'PRI'
        
        lea_obj = self.getFileObj(self.findfiles('LEA_01.001')[0])
        lea = lea_obj.read()
        lea_obj.close()
        fdr = lea[0:720]  # file descriptor record
        dss = lea[720:(720 + 1886)]  # data set summary record
        mpd = lea[(720 + 1886):(720 + 1886 + 1620)]  # map projection data record
        ppd_start = 720 + 1886 + 1620
        ppd_length = struct.unpack('>i', lea[ppd_start + 8: ppd_start + 12])[0]
        ppd = lea[ppd_start:ppd_length]  # platform position data record
        frd_start = 720 + 1886 + 1620 + ppd_length
        frd = lea[frd_start:(frd_start + 12288)]  # facility related data record
        
        meta['sensor'] = dss[396:412].strip().decode()
        meta['start'] = self.parse_date(str(dss[1814:1838].decode('utf-8')))
        meta['stop'] = self.parse_date(str(dss[1862:1886].decode('utf-8')))
        meta['polarizations'] = ['VV']
        looks_range = float(dss[1174:1190])
        looks_azimuth = float(dss[1190:1206])
        meta['looks'] = (looks_range, looks_azimuth)
        meta['heading'] = float(dss[468:476])
        meta['orbit'] = 'D' if meta['heading'] > 180 else 'A'
        orbitNumber, frameNumber = map(int, re.findall('[0-9]+', dss[36:68].decode('utf-8')))
        meta['orbitNumber_abs'] = orbitNumber
        meta['frameNumber'] = frameNumber
        orbitInfo = passdb_query(meta['sensor'], datetime.strptime(meta['start'], '%Y%m%dT%H%M%S'))
        meta['cycleNumber'] = orbitInfo['cycleNumber']
        meta['orbitNumber_rel'] = orbitInfo['orbitNumber_rel']
        spacing_azimuth = float(dss[1686:1702])
        spacing_range = float(dss[1702:1718])
        meta['spacing'] = (spacing_range, spacing_azimuth)
        meta['incidence_angle'] = float(dss[484:492])
        meta['proc_facility'] = dss[1045:1061].strip().decode()
        meta['proc_system'] = dss[1061:1069].strip().decode()
        meta['proc_version'] = dss[1069:1077].strip().decode()
        
        meta['antenna_flag'] = int(frd[658:662])
        meta['k_db'] = -10 * math.log(float(frd[662:678]), 10)
        meta['sc_db'] = {'ERS1': 59.61, 'ERS2': 60}[meta['sensor']]
        
        meta['samples'] = int(mpd[60:76])
        meta['lines'] = int(mpd[76:92])
        ul = (float(mpd[1088:1104]), float(mpd[1072:1088]))
        ur = (float(mpd[1120:1136]), float(mpd[1104:1120]))
        lr = (float(mpd[1152:1168]), float(mpd[1136:1152]))
        ll = (float(mpd[1184:1200]), float(mpd[1168:1184]))
        meta['coordinates'] = [ul, ur, lr, ll]
        meta['projection'] = crsConvert(4326, 'wkt')
        return meta
        
        # def correctAntennaPattern(self):
        # the following section is only relevant for PRI products and can be considered future work
        # select antenna gain correction lookup file from extracted meta information
        # the lookup files are stored in a subfolder CAL which is included in the pythonland software package
        # if sensor == 'ERS1':
        #     if date < 19950717:
        #         antenna = 'antenna_ERS1_x_x_19950716'
        #     else:
        #         if proc_sys == 'VMP':
        #             antenna = 'antenna_ERS2_VMP_v68_x' if proc_vrs >= 6.8 else 'antenna_ERS2_VMP_x_v67'
        #         elif proc_fac == 'UKPAF' and date < 19970121:
        #             antenna = 'antenna_ERS1_UKPAF_19950717_19970120'
        #         else:
        #             antenna = 'antenna_ERS1'
        # else:
        #     if proc_sys == 'VMP':
        #         antenna = 'antenna_ERS2_VMP_v68_x' if proc_vrs >= 6.8 else 'antenna_ERS2_VMP_x_v67'
        #     elif proc_fac == 'UKPAF' and date < 19970121:
        #         antenna = 'antenna_ERS2_UKPAF_x_19970120'
        #     else:
        #         antenna = 'antenna_ERS2'


class CEOS_PSR(ID):
    """
    Handler class for ALOS-PALSAR data in CEOS format
    
    Sensors:
        * PSR1
        * PSR2

    PALSAR-1:
        References:
            * NEB-01006: ALOS/PALSAR Level 1 Product Format Description
              (`JAXA 2006 <https://www.eorc.jaxa.jp/ALOS/en/doc/fdata/PALSAR_L10_J_ENa.zip>`_)
            * NEB-070062B: ALOS/PALSAR Level 1.1/1.5 Product Format Description
              (`JAXA 2009 <https://www.eorc.jaxa.jp/ALOS/en/doc/fdata/PALSAR_x_Format_EL.pdf>`_)
        Products / processing levels:
            * 1.0
            * 1.1
            * 1.5
        Acquisition modes:
            * AB: [SP][HWDPC]
            * A: supplemental remarks of the sensor type:
                * S: Wide observation mode
                * P: all other modes
            * B: observation mode
                * H: Fine mode
                * W: ScanSAR mode
                * D: Direct downlink mode
                * P: Polarimetry mode
                * C: Calibration mode
    
    PALSAR-2:
        Reference:
            ALOS-2/PALSAR-2 Level 1.1/1.5/2.1/3.1 CEOS SAR Product Format Description
            (`JAXA 2014 <https://www.eorc.jaxa.jp/ALOS-2/en/doc/fdata/PALSAR-2_xx_Format_CEOS_E_r.pdf>`_).
        Products / processing levels:
            * 1.0
            * 1.1
            * 1.5
        Acquisition modes:
            * SBS: Spotlight mode
            * UBS: Ultra-fine mode Single polarization
            * UBD: Ultra-fine mode Dual polarization
            * HBS: High-sensitive mode Single polarization
            * HBD: High-sensitive mode Dual polarization
            * HBQ: High-sensitive mode Full (Quad.) polarimetry
            * FBS: Fine mode Single polarization
            * FBD: Fine mode Dual polarization
            * FBQ: Fine mode Full (Quad.) polarimetry
            * WBS: Scan SAR nominal [14MHz] mode Single polarization
            * WBD: Scan SAR nominal [14MHz] mode Dual polarization
            * WWS: Scan SAR nominal [28MHz] mode Single polarization
            * WWD: Scan SAR nominal [28MHz] mode Dual polarization
            * VBS: Scan SAR wide mode Single polarization
            * VBD: Scan SAR wide mode Dual polarization
    """
    
    def __init__(self, scene):
        
        self.scene = os.path.realpath(scene)
        
        candidates = [patterns.ceos_psr1, patterns.ceos_psr2]
        
        for i, pattern in enumerate(candidates):
            self.pattern = pattern
            try:
                self.examine()
                break
            except RuntimeError as e:
                if i + 1 == len(candidates):
                    raise e
        
        self.meta = self.scanMetadata()
        
        # register the standardized meta attributes as object attributes
        super(CEOS_PSR, self).__init__(self.meta)
    
    def _getLeaderfileContent(self):
        led_obj = self.getFileObj(self.led_filename)
        led = led_obj.read()
        led_obj.close()
        return led
    
    def _img_get_coordinates(self):
        img_filename = self.findfiles('IMG')[0]
        img_obj = self.getFileObj(img_filename)
        imageFileDescriptor = img_obj.read(720)
        
        lineRecordLength = int(imageFileDescriptor[186:192])  # bytes per line + 412
        numberOfRecords = int(imageFileDescriptor[180:186])
        
        signalDataDescriptor1 = img_obj.read(412)
        img_obj.seek(720 + lineRecordLength * (numberOfRecords - 1))
        signalDataDescriptor2 = img_obj.read()
        
        img_obj.close()
        
        lat = [signalDataDescriptor1[192:196], signalDataDescriptor1[200:204],
               signalDataDescriptor2[192:196], signalDataDescriptor2[200:204]]
        
        lon = [signalDataDescriptor1[204:208], signalDataDescriptor1[212:216],
               signalDataDescriptor2[204:208], signalDataDescriptor2[212:216]]
        
        lat = [struct.unpack('>i', x)[0] / 1000000. for x in lat]
        lon = [struct.unpack('>i', x)[0] / 1000000. for x in lon]
        
        return list(zip(lon, lat))
    
    def _parseSummary(self):
        try:
            summary_file = self.getFileObj(self.findfiles('summary|workreport')[0])
        except IndexError:
            return {}
        text = summary_file.getvalue().decode('utf-8').strip()
        summary_file.close()
        summary = ast.literal_eval('{"' + re.sub(r'\s*=', '":', text).replace('\n', ',"') + '}')
        for x, y in summary.items():
            summary[x] = parse_literal(y)
        return summary
    
    @property
    def led_filename(self):
        return self.findfiles(self.pattern)[0]
    
    def scanMetadata(self):
        ################################################################################################################
        # read leader (LED) file
        led = self._getLeaderfileContent()
        
        # read summary text file
        meta = self._parseSummary()
        
        # read polarizations from image file names
        meta['polarizations'] = [re.search('[HV]{2}', os.path.basename(x)).group(0) for x in self.findfiles('^IMG-')]
        ################################################################################################################
        # read start and stop time
        
        try:
            meta['start'] = self.parse_date(meta['Img_SceneStartDateTime'])
            meta['stop'] = self.parse_date(meta['Img_SceneEndDateTime'])
        except (AttributeError, KeyError):
            try:
                start_string = re.search('Img_SceneStartDateTime[ ="0-9:.]*', led).group()
                stop_string = re.search('Img_SceneEndDateTime[ ="0-9:.]*', led).group()
                meta['start'] = self.parse_date(re.search(r'\d+\s[\d:.]+', start_string).group())
                meta['stop'] = self.parse_date(re.search(r'\d+\s[\d:.]+', stop_string).group())
            except AttributeError:
                raise IndexError('start and stop time stamps cannot be extracted; see file {}'
                                 .format(self.led_filename))
        ################################################################################################################
        # read file descriptor record
        p0 = 0
        p1 = struct.unpack('>i', led[8:12])[0]
        fileDescriptor = led[p0:p1]
        # dataSetSummary
        dss_n = int(fileDescriptor[180:186])
        dss_l = int(fileDescriptor[186:192])
        # mapProjectionData
        mpd_n = int(fileDescriptor[192:198])
        mpd_l = int(fileDescriptor[198:204])
        # platformPositionData
        ppd_n = int(fileDescriptor[204:210])
        ppd_l = int(fileDescriptor[210:216])
        # attitudeData
        adr_n = int(fileDescriptor[216:222])
        adr_l = int(fileDescriptor[222:228])
        # radiometricData
        rdr_n = int(fileDescriptor[228:234])
        rdr_l = int(fileDescriptor[234:240])
        # dataQualitySummary
        dqs_n = int(fileDescriptor[252:258])
        dqs_l = int(fileDescriptor[258:264])
        meta['sensor'] = {'AL1': 'PSR1', 'AL2': 'PSR2'}[fileDescriptor[48:51].decode('utf-8')]
        ################################################################################################################
        # read leader file name information
        
        match = re.match(re.compile(self.pattern), os.path.basename(self.led_filename))
        
        if meta['sensor'] == 'PSR1':
            meta['acquisition_mode'] = match.group('sub') + match.group('mode')
        else:
            meta['acquisition_mode'] = match.group('mode')
        meta['product'] = match.group('level')
        ################################################################################################################
        # read led records
        p0 = p1
        p1 += dss_l * dss_n
        dataSetSummary = led[p0:p1]
        
        if mpd_n > 0:
            p0 = p1
            p1 += mpd_l * mpd_n
            mapProjectionData = led[p0:p1]
        else:
            mapProjectionData = None
        
        p0 = p1
        p1 += ppd_l * ppd_n
        platformPositionData = led[p0:p1]
        
        p0 = p1
        p1 += adr_l * adr_n
        attitudeData = led[p0:p1]
        
        p0 = p1
        p1 += rdr_l * rdr_n
        radiometricData = led[p0:p1]
        
        p0 = p1
        p1 += dqs_l * dqs_n
        dataQualitySummary = led[p0:p1]
        
        facilityRelatedData = []
        while p1 < len(led):
            p0 = p1
            length = struct.unpack('>i', led[(p0 + 8):(p0 + 12)])[0]
            p1 += length
            facilityRelatedData.append(led[p0:p1])
        ################################################################################################################
        # read map projection data record
        
        if mapProjectionData is not None:
            lat = list(map(float, [mapProjectionData[1072:1088],
                                   mapProjectionData[1104:1120],
                                   mapProjectionData[1136:1152],
                                   mapProjectionData[1168:1184]]))
            lon = list(map(float, [mapProjectionData[1088:1104],
                                   mapProjectionData[1120:1136],
                                   mapProjectionData[1152:1168],
                                   mapProjectionData[1184:1200]]))
            meta['coordinates'] = list(zip(lon, lat))
            
            # https://github.com/datalyze-solutions/LandsatProcessingPlugin/blob/master/src/metageta/formats/alos.py
            
            src_srs = osr.SpatialReference()
            # src_srs.SetGeogCS('GRS 1980','GRS 1980','GRS 1980',6378137.00000,298.2572220972)
            src_srs.SetWellKnownGeogCS('WGS84')
            # Proj CS
            projdesc = mapProjectionData[412:444].strip()
            epsg = 0  # default
            if projdesc == 'UTM-PROJECTION':
                nZone = int(mapProjectionData[476:480])
                dfFalseNorthing = float(mapProjectionData[496:512])
                if dfFalseNorthing > 0.0:
                    bNorth = False
                    epsg = 32700 + nZone
                else:
                    bNorth = True
                    epsg = 32600 + nZone
                src_srs.ImportFromEPSG(epsg)
                # src_srs.SetUTM(nZone,bNorth) #generates WKT that osr.SpatialReference.AutoIdentifyEPSG() doesn't return an EPSG for
            elif projdesc == 'UPS-PROJECTION':
                dfCenterLon = float(mapProjectionData[624, 640])
                dfCenterLat = float(mapProjectionData[640, 656])
                dfScale = float(mapProjectionData[656, 672])
                src_srs.SetPS(dfCenterLat, dfCenterLon, dfScale, 0.0, 0.0)
            elif projdesc == 'MER-PROJECTION':
                dfCenterLon = float(mapProjectionData[736, 752])
                dfCenterLat = float(mapProjectionData[752, 768])
                src_srs.SetMercator(dfCenterLat, dfCenterLon, 0, 0, 0)
            elif projdesc == 'LCC-PROJECTION':
                dfCenterLon = float(mapProjectionData[736, 752])
                dfCenterLat = float(mapProjectionData[752, 768])
                dfStdP1 = float(mapProjectionData[768, 784])
                dfStdP2 = float(mapProjectionData[784, 800])
                src_srs.SetLCC(dfStdP1, dfStdP2, dfCenterLat, dfCenterLon, 0, 0)
            meta['projection'] = src_srs.ExportToWkt()
        
        else:
            coordinates = self._img_get_coordinates()
            if all([x == (0, 0) for x in coordinates]):
                meta['projection'] = None
            else:
                meta['coordinates'] = coordinates
                meta['projection'] = crsConvert(4326, 'wkt')
        ################################################################################################################
        # read data set summary record
        
        scene_id = dataSetSummary[20:52].decode('ascii')
        
        if meta['sensor'] == 'PSR1':
            pattern = r'(?P<sat_id>[A-Z]{2})' \
                      r'(?P<sensor_id>[A-Z]{3})' \
                      r'(?P<sensor_id_sub>[A-Z]{1})' \
                      r'(?P<orbitNumber>[0-9]{5})' \
                      r'(?P<frameNumber>[0-9]{4})'
        elif meta['sensor'] == 'PSR2':
            pattern = r'(?P<sat_id>[A-Z0-9]{5})' \
                      r'(?P<orbitNumber>[0-9]{5})' \
                      r'(?P<frameNumber>[0-9]{4})-' \
                      r'(?P<obs_day>[0-9]{6})[ ]{11}'
        else:
            raise ValueError('sensor must be either PSR1 or PSR2; is: {}'.format(meta['sensor']))
        
        match = re.match(re.compile(pattern), scene_id)
        
        orbitsPerCycle = {'PSR1': 671, 'PSR2': 207}[meta['sensor']]
        
        meta['orbitNumber_abs'] = int(match.group('orbitNumber'))
        meta['orbitNumber_rel'] = meta['orbitNumber_abs'] % orbitsPerCycle
        meta['cycleNumber'] = meta['orbitNumber_abs'] // orbitsPerCycle + 1
        meta['frameNumber'] = int(match.group('frameNumber'))
        
        try:
            meta['lines'] = int(dataSetSummary[324:332]) * 2
        except ValueError:
            if 'Pdi_NoOfLines' in meta.keys():
                meta['lines'] = meta['Pdi_NoOfLines']
            else:
                meta['lines'] = None
        try:
            meta['samples'] = int(dataSetSummary[332:340]) * 2
        except ValueError:
            if 'Pdi_NoOfPixels' in meta.keys():
                meta['samples'] = meta['Pdi_NoOfPixels']
            else:
                meta['samples'] = None
        meta['incidence'] = float(dataSetSummary[484:492])
        meta['wavelength'] = float(dataSetSummary[500:516]) * 100  # in cm
        meta['proc_facility'] = dataSetSummary[1046:1062].strip()
        meta['proc_system'] = dataSetSummary[1062:1070].strip()
        meta['proc_version'] = dataSetSummary[1070:1078].strip()
        
        try:
            azlks = float(dataSetSummary[1174:1190])
            rlks = float(dataSetSummary[1190:1206])
            meta['looks'] = (rlks, azlks)
        except ValueError:
            meta['looks'] = (None, None)
        
        meta['orbit'] = dataSetSummary[1534:1542].decode('utf-8').strip()[0]
        
        try:
            spacing_azimuth = float(dataSetSummary[1686:1702])
            spacing_range = float(dataSetSummary[1702:1718])
            meta['spacing'] = (spacing_range, spacing_azimuth)
        except ValueError:
            meta['spacing'] = (None, None)
        ################################################################################################################
        # read radiometric data record
        if len(radiometricData) > 0:
            meta['k_dB'] = float(radiometricData[20:36])
        else:
            meta['k_dB'] = None
        ################################################################################################################
        # additional notes
        
        # the following can be used to read platform position time from the led file
        # this covers a larger time frame than the actual scene sensing time
        # y, m, d, nd, s = platformPositionData[144:182].split()
        # start = datetime(int(y), int(m), int(d)) + timedelta(seconds=float(s))
        # npoints = int(platformPositionData[140:144])
        # interval = float(platformPositionData[182:204])
        # stop = start + timedelta(seconds=(npoints - 1) * interval)
        # parse_date(start)
        # parse_date(stop)
        
        return meta
    
    def unpack(self, directory, overwrite=False, exist_ok=False):
        outdir = os.path.join(directory, os.path.basename(self.file).replace('LED-', ''))
        self._unpack(outdir, overwrite=overwrite, exist_ok=exist_ok)


class EORC_PSR(ID):
    """
    Handler class for ALOS-2/PALSAR-2 data in EORC (Earth Observation Research Center) Path format
    
    Sensors:
        * PALSAR-2

    PALSAR-2:
        Reference: 
            NDX-150019: ALOS-2/PALSAR-2 EORC Path Product Format Description (JAXA 2016)
        Products / processing levels:
            * 1.5
        Acquisition modes:
            * FBD: Fine mode Dual polarization
            * WBD: Scan SAR nominal [14MHz] mode Dual polarization
    """
    
    def __init__(self, scene):
        
        self.scene = os.path.realpath(scene)
        
        self.pattern = patterns.eorc_psr
        
        self.examine()
        
        self.meta = self.scanMetadata()
        
        # register the standardized meta attributes as object attributes
        super(EORC_PSR, self).__init__(self.meta)
    
    def _getHeaderfileContent(self):
        head_obj = self.getFileObj(self.header_filename)
        head = head_obj.read().decode('utf-8')
        head = list(head.split('\n'))
        head_obj.close()
        return head
    
    def _img_get_coordinates(self):
        img_filename = self.findfiles('IMG')[0]
        img_obj = self.getFileObj(img_filename)
        imageFileDescriptor = img_obj.read(720)
        
        lineRecordLength = int(imageFileDescriptor[186:192])  # bytes per line + 412
        numberOfRecords = int(imageFileDescriptor[180:186])
        
        signalDataDescriptor1 = img_obj.read(412)
        img_obj.seek(720 + lineRecordLength * (numberOfRecords - 1))
        signalDataDescriptor2 = img_obj.read()
        
        img_obj.close()
        
        lat = [signalDataDescriptor1[192:196], signalDataDescriptor1[200:204],
               signalDataDescriptor2[192:196], signalDataDescriptor2[200:204]]
        
        lon = [signalDataDescriptor1[204:208], signalDataDescriptor1[212:216],
               signalDataDescriptor2[204:208], signalDataDescriptor2[212:216]]
        
        lat = [struct.unpack('>i', x)[0] / 1000000. for x in lat]
        lon = [struct.unpack('>i', x)[0] / 1000000. for x in lon]
        
        return list(zip(lon, lat))
    
    def _parseFacter_m(self):
        try:
            facter_file = self.findfiles('facter_m.dat')[0]
        except IndexError:
            return {}
        facter_obj = self.getFileObj(facter_file)
        facter_m = facter_obj.read().decode('utf-8')
        facter_m = list(facter_m.split('\n'))
        facter_obj.close()
        return facter_m
    
    @property
    def header_filename(self):
        return self.findfiles(self.pattern)[0]
    
    def scanMetadata(self):
        ################################################################################################################
        # read header (HDR) file
        header = self._getHeaderfileContent()
        header = [head.replace(" ", "") for head in header]
        
        # read summary text file
        facter_m = self._parseFacter_m()
        facter_m = [fact.replace(" ", "") for fact in facter_m]
        
        meta = {}
        
        # read polarizations from image file names
        meta['polarizations'] = [re.search('[HV]{2}', os.path.basename(x)).group(0) for x in self.findfiles('^sar.')]
        meta['product'] = header[3]
        ################################################################################################################
        # read start and stop time --> TODO: in what format is the start and stop time?
        
        try:
            start_time = facter_m[168].split('.')[0].zfill(2) + facter_m[168].split('.')[1][:4]
            stop_time = facter_m[170].split('.')[0].zfill(2) + facter_m[170].split('.')[1][:4]
        except (AttributeError):
            raise IndexError('start and stop time stamps cannot be extracted; see file facter_m.dat')
        
        meta['start'] = str(header[6])  # +'T'+start_time
        meta['stop'] = str(header[6])  # +'T'+stop_time
        ################################################################################################################
        # read file metadata
        meta['sensor'] = header[2]
        ################################################################################################################
        # read leader file name information
        meta['acquisition_mode'] = header[12]
        # ##############################################################################################################
        # read map projection data 
        
        lat = list(map(float, [header[33], header[35], header[37], header[39]]))
        lon = list(map(float, [header[34], header[36], header[38], header[40]]))
        
        if len(lat) == 0 or len(lon) == 0:
            meta['coordinates'] = self._img_get_coordinates()
        else:
            meta['coordinates'] = list(zip(lon, lat))
        
        meta['projection'] = crsConvert(4918, 'wkt')  # EPSG: 4918: ITRF97, GRS80
        ################################################################################################################
        # read data set summary record
        
        orbitsPerCycle = int(207)
        
        meta['orbitNumber_rel'] = int(header[7])
        meta['cycleNumber'] = header[5]
        meta['frameNumber'] = ''
        meta['orbitNumber_abs'] = int(orbitsPerCycle * (meta['cycleNumber'] - 1) + meta['orbitNumber_rel'])
        
        meta['lines'] = int(float(facter_m[51]))
        meta['samples'] = int(float(facter_m[50]))
        meta['incidence'] = float(facter_m[119])
        meta['proc_facility'] = header[73]
        
        meta['spacing'] = (float(header[51]), float(header[52]))
        
        meta['orbit'] = header[9]
        ################################################################################################################
        # read radiometric data record
        
        meta['k_dB'] = float(header[64])
        
        return meta
    
    def unpack(self, directory, overwrite=False, exist_ok=False):
        outdir = os.path.join(directory, os.path.basename(self.file).replace('LED-', ''))
        self._unpack(outdir, overwrite=overwrite, exist_ok=exist_ok)


class ESA(ID):
    """
    Handler class for SAR data in ESA format (Envisat ASAR, ERS-1/2)
    
    Sensors:
        * ASAR
        * ERS1
        * ERS2
    """
    
    def __init__(self, scene):
        
        self.pattern = patterns.esa
        self.pattern_pid = r'(?P<sat_id>(?:SAR|ASA))_' \
                           r'(?P<image_mode>(?:IM(?:S|P|G|M|_)|AP(?:S|P|G|M|_)|WV(?:I|S|W|_)|WS(?:M|S|_)))_' \
                           r'(?P<processing_level>[012B][CP])'
        
        self.scene = os.path.realpath(scene)
        
        if re.search('.[EN][12]$', self.scene):
            self.file = self.scene
        else:
            self.examine()
        
        self.meta = self.scanMetadata()
        
        # register the standardized meta attributes as object attributes
        super(ESA, self).__init__(self.meta)
    
    def scanMetadata(self):
        match = re.match(re.compile(self.pattern), os.path.basename(self.file))
        match2 = re.match(re.compile(self.pattern_pid), match.group('product_id'))
        
        if re.search('IM__0', match.group('product_id')):
            raise RuntimeError('product level 0 not supported (yet)')
        
        meta = dict()
        sensor_lookup = {'N1': 'ASAR', 'E1': 'ERS1', 'E2': 'ERS2'}
        meta['sensor'] = sensor_lookup[match.group('satellite_ID')]
        meta['acquisition_mode'] = match2.group('image_mode')
        
        meta['image_geometry'] = 'GROUND_RANGE'
        # product overview table: https://doi.org/10.5167/UZH-96146
        if meta['acquisition_mode'] in ['APS', 'IMS', 'WSS']:
            meta['product'] = 'SLC'
            meta['image_geometry'] = 'SLANT_RANGE'
        elif meta['acquisition_mode'] in ['APP', 'IMP']:
            meta['product'] = 'PRI'
        elif meta['acquisition_mode'] in ['APM', 'IMM', 'WSM']:
            meta['product'] = 'MR'
        else:
            raise RuntimeError(f"unsupported acquisition mode: '{meta['acquisition_mode']}'")
        
        def val_convert(val):
            try:
                out = int(val)
            except ValueError:
                try:
                    out = float(val)
                except ValueError:
                    if re.search('[0-9]{2}-[A-Z]{3}-[0-9]{2}', val):
                        out = dateparse(val)
                        out = out.replace(tzinfo=timezone.utc)
                    else:
                        out = val
            return out
        
        def decode(raw):
            pattern = r'(?P<key>[A-Z0-9_]+)\=(")?(?P<value>.*?)("|<|$)'
            out = {}
            coord_keys = [f'{x}_{y}_{z}'
                          for x in ['FIRST', 'LAST']
                          for y in ['NEAR', 'MID', 'FAR']
                          for z in ['LAT', 'LONG']]
            lines = raw.split('\n')
            for line in lines:
                match = re.match(pattern, line)
                if match:
                    matchdict = match.groupdict()
                    val = val_convert(str(matchdict['value']).strip())
                    if matchdict['key'] in coord_keys:
                        val *= 10 ** -6
                    out[matchdict['key']] = val
            return out
        
        with self.getFileObj(self.file) as obj:
            origin = {}
            mph = obj.read(1247).decode('ascii')
            origin['MPH'] = decode(mph)
            
            sph_size = origin['MPH']['SPH_SIZE']
            dsd_size = origin['MPH']['DSD_SIZE']
            dsd_num = origin['MPH']['NUM_DSD']
            sph_descr_size = sph_size - dsd_size * dsd_num
            
            sph = obj.read(sph_descr_size).decode('ascii')
            origin['SPH'] = decode(sph)
            
            datasets = {}
            for i in range(dsd_num):
                dsd = obj.read(dsd_size).decode('ascii')
                dataset = decode(dsd)
                datasets[dataset.pop('DS_NAME')] = dataset
            origin['DSD'] = datasets
            
            meta['origin'] = origin
            
            key = 'GEOLOCATION GRID ADS'
            ds_offset = origin['DSD'][key]['DS_OFFSET']
            ds_size = origin['DSD'][key]['DS_SIZE']
            dsr_size = origin['DSD'][key]['DSR_SIZE']
            obj.seek(ds_offset)
            geo = obj.read(ds_size)
        
        geo = [geo[i:i + dsr_size] for i in range(0, len(geo), dsr_size)]
        
        keys = ['first_zero_doppler_time', 'attach_flag', 'line_num',
                'num_lines', 'sub_sat_track', 'first_line_tie_points',
                'spare', 'last_zero_doppler_time', 'last_line_tie_points',
                'swath_number']
        lengths = [12, 1, 4, 4, 4, 220, 22, 12, 220, 3, 19]
        
        meta['origin']['GEOLOCATION_GRID_ADS'] = []
        for granule in geo:
            start = 0
            values = {}
            for i, key in enumerate(keys):
                value = granule[start:sum(lengths[:i + 1])]
                if key in ['first_zero_doppler_time', 'last_zero_doppler_time']:
                    unpack = dict(zip(('days', 'seconds', 'microseconds'),
                                      struct.unpack('>lLL', value)))
                    value = datetime(year=2000, month=1, day=1, tzinfo=timezone.utc)
                    value += timedelta(**unpack)
                elif key in ['attach_flag']:
                    value = struct.unpack('B', value)[0]
                elif key in ['line_num', 'num_lines']:
                    value = struct.unpack('>L', value)[0]
                elif key in ['sub_sat_track']:
                    value = struct.unpack('>f', value)[0]
                elif key in ['first_line_tie_points', 'last_line_tie_points']:
                    sample_numbers = struct.unpack('>' + 'L' * 11, value[0:44])
                    slant_range_times = struct.unpack('>' + 'f' * 11, value[44:88])
                    incident_angles = struct.unpack('>' + 'f' * 11, value[88:132])
                    latitudes = struct.unpack('>' + 'l' * 11, value[132:176])
                    latitudes = [x / 1000000. for x in latitudes]
                    longitudes = struct.unpack('>' + 'l' * 11, value[176:220])
                    longitudes = [x / 1000000. for x in longitudes]
                    value = []
                    for j in range(11):
                        value.append({'sample_number': sample_numbers[j],
                                      'slant_range_time': slant_range_times[j],
                                      'incident_angle': incident_angles[j],
                                      'latitude': latitudes[j],
                                      'longitude': longitudes[j]})
                elif key == 'swath_number':
                    value = value.decode('ascii').strip()
                if key != 'spare':
                    values[key] = value
                start += lengths[i]
            meta['origin']['GEOLOCATION_GRID_ADS'].append(values)
        
        lat = []
        lon = []
        for granule in meta['origin']['GEOLOCATION_GRID_ADS']:
            for group in ['first', 'last']:
                for i in range(11):
                    lat.append(granule[f'{group}_line_tie_points'][i]['latitude'])
                    lon.append(granule[f'{group}_line_tie_points'][i]['longitude'])
        
        meta['coordinates'] = list(zip(lon, lat))
        
        if meta['sensor'] == 'ASAR':
            pols = [y for x, y in origin['SPH'].items() if 'TX_RX_POLAR' in x]
            pols = [x.replace('/', '') for x in pols if len(x) == 3]
            meta['polarizations'] = sorted(pols)
        elif meta['sensor'] in ['ERS1', 'ERS2']:
            meta['polarizations'] = ['VV']
        
        meta['orbit'] = origin['SPH']['PASS'][0]
        meta['start'] = origin['MPH']['SENSING_START'].strftime('%Y%m%dT%H%M%S')
        meta['stop'] = origin['MPH']['SENSING_STOP'].strftime('%Y%m%dT%H%M%S')
        meta['spacing'] = (origin['SPH']['RANGE_SPACING'], origin['SPH']['AZIMUTH_SPACING'])
        meta['looks'] = (origin['SPH']['RANGE_LOOKS'], origin['SPH']['AZIMUTH_LOOKS'])
        meta['samples'] = origin['SPH']['LINE_LENGTH']
        meta['lines'] = origin['DSD']['MDS1']['NUM_DSR']
        
        meta['orbitNumber_abs'] = origin['MPH']['ABS_ORBIT']
        meta['orbitNumber_rel'] = origin['MPH']['REL_ORBIT']
        meta['cycleNumber'] = origin['MPH']['CYCLE']
        meta['frameNumber'] = origin['MPH']['ABS_ORBIT']
        
        incident_angles = []
        for item in meta['origin']['GEOLOCATION_GRID_ADS']:
            for key in ['first', 'last']:
                pts = item[f'{key}_line_tie_points']
                for pt in pts:
                    incident_angles.append(pt['incident_angle'])
        
        meta['incidence_nr'] = min(incident_angles)
        meta['incidence_fr'] = max(incident_angles)
        meta['incidence'] = (meta['incidence_nr'] + meta['incidence_fr']) / 2
        
        resolution_rg, resolution_az, nesz_nr, nesz_fr = \
            get_resolution_nesz(sensor=meta['sensor'], mode=meta['acquisition_mode'],
                                swath_id=origin['SPH']['SWATH'], date=meta['start'])
        
        meta['resolution'] = (resolution_rg, resolution_az)
        meta['nesz'] = (nesz_nr, nesz_fr)
        
        meta['projection'] = crsConvert(4326, 'wkt')
        
        return meta
    
    def geo_grid(self, outname=None, driver=None, overwrite=True):
        """
        get the geo grid as vector geometry

        Parameters
        ----------
        outname: str
            the name of the vector file to be written
        driver: str
            the output file format; needs to be defined if the format cannot
            be auto-detected from the filename extension
        overwrite: bool
            overwrite an existing vector file?

        Returns
        -------
        spatialist.vector.Vector or None
            the vector object if `outname` is None, None otherwise

        See also
        --------
        spatialist.vector.Vector.write
        """
        vec = Vector(driver='Memory')
        vec.addlayer('geogrid', 4326, ogr.wkbPoint)
        field_defs = [
            ("swath", ogr.OFTString),
            ("azimuthTime", ogr.OFTDateTime),
            ("slantRangeTime", ogr.OFTReal),
            ("line", ogr.OFTInteger),
            ("pixel", ogr.OFTInteger),
            ("incidenceAngle", ogr.OFTReal)
        ]
        for name, ftype in field_defs:
            field = ogr.FieldDefn(name, ftype)
            vec.layer.CreateField(field)
        
        for granule in self.meta['origin']['GEOLOCATION_GRID_ADS']:
            line_first = granule['line_num']
            line_last = granule['line_num'] + granule['num_lines'] - 1
            for group in ['first', 'last']:
                meta = {'swath': granule['swath_number'],
                        'azimuthTime': granule[f'{group}_zero_doppler_time'],
                        'line': line_first if group == 'first' else line_last}
                tp = granule[f'{group}_line_tie_points']
                for i in range(11):
                    x = tp[i]['longitude']
                    y = tp[i]['latitude']
                    geom = ogr.Geometry(ogr.wkbPoint)
                    geom.AddPoint(x, y)
                    geom.FlattenTo2D()
                    meta['slantRangeTime'] = tp[i]['slant_range_time']
                    meta['pixel'] = tp[i]['sample_number']
                    meta['incidenceAngle'] = tp[i]['incident_angle']
                    vec.addfeature(geom, fields=meta)
        geom = None
        if outname is None:
            return vec
        else:
            vec.write(outfile=outname, driver=driver, overwrite=overwrite)
    
    def unpack(self, directory, overwrite=False, exist_ok=False):
        base_file = os.path.basename(self.file).strip(r'\.zip|\.tar(?:\.gz|)')
        base_dir = os.path.basename(directory.strip('/'))
        
        outdir = directory if base_file == base_dir else os.path.join(directory, base_file)
        
        self._unpack(outdir, overwrite=overwrite, exist_ok=exist_ok)


class SAFE(ID):
    """
    Handler class for Sentinel-1 data
    
    Sensors:
        * S1A
        * S1B
        * S1C
        * S1D

    References:
        * S1-RS-MDA-52-7443 Sentinel-1 IPF Auxiliary Product Specification
        * MPC-0243 Masking "No-value" Pixels on GRD Products generated by the Sentinel-1 ESA IPF
    """
    
    def __init__(self, scene):
        
        self.scene = os.path.realpath(scene)
        
        self.pattern = patterns.safe
        
        self.pattern_ds = r'^s1[abcd]-' \
                          r'(?P<swath>s[1-6]|iw[1-3]?|ew[1-5]?|wv[1-2]|n[1-6])-' \
                          r'(?P<product>slc|grd|ocn)-' \
                          r'(?P<pol>hh|hv|vv|vh)-' \
                          r'(?P<start>[0-9]{8}t[0-9]{6})-' \
                          r'(?P<stop>[0-9]{8}t[0-9]{6})-' \
                          r'(?:[0-9]{6})-(?:[0-9a-f]{6})-' \
                          r'(?P<id>[0-9]{3})' \
                          r'\.xml$'
        
        self.examine(include_folders=True)
        
        if not re.match(re.compile(self.pattern), os.path.basename(self.file)):
            raise RuntimeError('folder does not match S1 scene naming convention')
        
        # scan the metadata XML file and add selected attributes to a meta dictionary
        self.meta = self.scanMetadata()
        self.meta['projection'] = crsConvert(4326, 'wkt')
        
        # register the standardized meta attributes as object attributes
        super(SAFE, self).__init__(self.meta)
        
        self.gammafiles = {'slc': [], 'pri': [], 'grd': []}
    
    def removeGRDBorderNoise(self, method='pyroSAR'):
        """
        mask out Sentinel-1 image border noise.
        
        Parameters
        ----------
        method: str
            the border noise removal method to be applied; one of the following:
            
             - 'ESA': the pure implementation as described by ESA
             - 'pyroSAR': the ESA method plus the custom pyroSAR refinement

        Returns
        -------
        
        See Also
        --------
        :func:`~pyroSAR.S1.removeGRDBorderNoise`
        """
        S1.removeGRDBorderNoise(self, method=method)
    
    def geo_grid(self, outname=None, driver=None, overwrite=True):
        """
        get the geo grid as vector geometry

        Parameters
        ----------
        outname: str
            the name of the vector file to be written
        driver: str
            the output file format; needs to be defined if the format cannot
            be auto-detected from the filename extension
        overwrite: bool
            overwrite an existing vector file?

        Returns
        -------
        ~spatialist.vector.Vector or None
            the vector object if `outname` is None, None otherwise
        
        See also
        --------
        spatialist.vector.Vector.write
        """
        annotations = self.findfiles(self.pattern_ds)
        key = lambda x: re.search('-[vh]{2}-', x).group()
        groups = groupby(sorted(annotations, key=key), key=key)
        annotations = [list(value) for key, value in groups][0]
        
        vec = Vector(driver='Memory')
        vec.addlayer('geogrid', 4326, ogr.wkbPoint25D)
        field_defs = [
            ("swath", ogr.OFTString),
            ("azimuthTime", ogr.OFTDateTime),
            ("slantRangeTime", ogr.OFTReal),
            ("line", ogr.OFTInteger),
            ("pixel", ogr.OFTInteger),
            ("incidenceAngle", ogr.OFTReal),
            ("elevationAngle", ogr.OFTReal),
        ]
        for name, ftype in field_defs:
            field = ogr.FieldDefn(name, ftype)
            vec.layer.CreateField(field)
        
        for ann in annotations:
            with self.getFileObj(ann) as ann_xml:
                tree = ET.fromstring(ann_xml.read())
            swath = tree.find(".//adsHeader/swath").text
            points = tree.findall(".//geolocationGridPoint")
            for point in points:
                meta = {child.tag: child.text for child in point}
                meta["swath"] = swath
                x = float(meta.pop("longitude"))
                y = float(meta.pop("latitude"))
                z = float(meta.pop("height"))
                geom = ogr.Geometry(ogr.wkbPoint25D)
                geom.AddPoint(x, y, z)
                az_time = dateparse(meta["azimuthTime"])
                meta["azimuthTime"] = az_time.replace(tzinfo=timezone.utc)
                for key in ["slantRangeTime", "incidenceAngle", "elevationAngle"]:
                    meta[key] = float(meta[key])
                for key in ["line", "pixel"]:
                    meta[key] = int(meta[key])
                vec.addfeature(geom, fields=meta)
        geom = None
        if outname is None:
            return vec
        else:
            vec.write(outfile=outname, driver=driver, overwrite=overwrite)
    
    def getOSV(self, osvdir=None, osvType='POE', returnMatch=False, useLocal=True, timeout=300, url_option=1):
        """
        download Orbit State Vector files for the scene

        Parameters
        ----------
        osvdir: str
            the directory of OSV files; subdirectories POEORB and RESORB are created automatically;
            if no directory is defined, the standard SNAP auxdata location is used
        osvType: str or list[str]
            the type of orbit file either 'POE', 'RES' or a list of both;
            if both are selected, the best matching file will be retrieved. I.e., POE if available and RES otherwise
        returnMatch: bool
            return the best matching orbit file?
        useLocal: bool
            use locally existing files and do not search for files online if the right file has been found?
        timeout: int or tuple or None
            the timeout in seconds for downloading OSV files as provided to :func:`requests.get`
        url_option: int
            the OSV download URL option; see :meth:`pyroSAR.S1.OSV.catch` for options

        Returns
        -------
        str or None
            the best matching OSV file if `returnMatch` is True or None otherwise
        
        See Also
        --------
        :class:`pyroSAR.S1.OSV`
        """
        with S1.OSV(osvdir, timeout=timeout) as osv:
            if useLocal:
                match = osv.match(sensor=self.sensor, timestamp=self.start,
                                  osvtype=osvType)
                if match is not None:
                    return match if returnMatch else None
            
            if osvType in ['POE', 'RES']:
                files = osv.catch(sensor=self.sensor, osvtype=osvType,
                                  start=self.start, stop=self.stop,
                                  url_option=url_option)
            elif sorted(osvType) == ['POE', 'RES']:
                files = osv.catch(sensor=self.sensor, osvtype='POE',
                                  start=self.start, stop=self.stop,
                                  url_option=url_option)
                if len(files) == 0:
                    files = osv.catch(sensor=self.sensor, osvtype='RES',
                                      start=self.start, stop=self.stop,
                                      url_option=url_option)
            else:
                msg = "osvType must either be 'POE', 'RES' or a list of both"
                raise TypeError(msg)
            
            osv.retrieve(files)
            
            if returnMatch:
                match = osv.match(sensor=self.sensor, timestamp=self.start,
                                  osvtype=osvType)
                return match
    
    def quicklook(self, outname, format='kmz', na_transparent=True):
        """
        Write a quicklook file for the scene.
        
        Parameters
        ----------
        outname: str
            the file to write
        format: str
            the quicklook format. Currently supported options:
            
             - kmz
        na_transparent: bool
            make NA values transparent?

        Returns
        -------

        """
        if self.product not in ['GRD', 'SLC']:
            msg = 'this method has only been implemented for GRD and SLC, not {}'
            raise RuntimeError(msg.format(self.product))
        
        if format != 'kmz':
            raise RuntimeError('currently only kmz is supported as format')
        kml_name = self.findfiles('map-overlay.kml')[0]
        png_name = self.findfiles('quick-look.png')[0]
        with zf.ZipFile(outname, 'w') as out:
            with self.getFileObj(kml_name) as kml_in:
                kml = kml_in.getvalue().decode('utf-8')
                kml = kml.replace('Sentinel-1 Map Overlay', self.outname_base())
                out.writestr('doc.kml', data=kml)
            with self.getFileObj(png_name) as png_in:
                if na_transparent:
                    img = Image.open(png_in)
                    img = img.convert('RGBA')
                    datas = img.getdata()
                    newData = []
                    for item in datas:
                        if item[0] == 0 and item[1] == 0 and item[2] == 0:
                            newData.append((0, 0, 0, 0))
                        else:
                            newData.append(item)
                    img.putdata(newData)
                    buf = BytesIO()
                    img.save(buf, format='png')
                    out.writestr('quick-look.png', buf.getvalue())
                else:
                    out.writestr('quick-look.png', data=png_in.getvalue())
    
    def resolution(self):
        """
        Compute the mid-swath resolution of the Sentinel-1 product. For GRD products the resolution is expressed in
        ground range and in slant range otherwise.
        
        References:
            * https://sentinel.esa.int/web/sentinel/user-guides/sentinel-1-sar/resolutions/level-1-single-look-complex
            * https://sentinel.esa.int/web/sentinel/user-guides/sentinel-1-sar/resolutions/level-1-ground-range-detected
            * https://sentinel.esa.int/web/sentinel/user-guides/sentinel-1-sar/document-library/-/asset_publisher/1dO7RF5fJMbd/content/sentinel-1-product-definition
        
        Returns
        -------
        tuple[float]
            the resolution as (range, azimuth)
        """
        if 'resolution' in self.meta.keys():
            return self.meta['resolution']
        if self.product not in ['GRD', 'SLC']:
            msg = 'this method has only been implemented for GRD and SLC, not {}'
            raise RuntimeError(msg.format(self.product))
        
        annotations = self.findfiles(self.pattern_ds)
        key = lambda x: re.search('-[vh]{2}-', x).group()
        groups = groupby(sorted(annotations, key=key), key=key)
        annotations = [list(value) for key, value in groups][0]
        proc_pars = []  # processing parameters per sub-swath
        sp_az = []  # azimuth pixel spacings per sub-swath
        ti_az = []  # azimuth time intervals per sub-swath
        for ann in annotations:
            with self.getFileObj(ann) as ann_xml:
                tree = ET.fromstring(ann_xml.read())
                par = tree.findall('.//swathProcParams')
                proc_pars.extend(par)
                for i in range(len(par)):
                    sp_az.append(float(tree.find('.//azimuthPixelSpacing').text))
                    ti_az.append(float(tree.find('.//azimuthTimeInterval').text))
        c = 299792458.0  # speed of light
        # see Sentinel-1 product definition for Hamming window coefficients
        # and Impulse Response Width (IRW) broadening factors:
        coefficients = [0.52, 0.6, 0.61, 0.62, 0.63, 0.65, 0.70, 0.72, 0.73, 0.75]
        b_factors = [1.54, 1.32, 1.3, 1.28, 1.27, 1.24, 1.18, 1.16, 1.15, 1.13]
        resolutions_rg = []
        resolutions_az = []
        for i, par in enumerate(proc_pars):
            # computation of slant range resolution
            rg_proc = par.find('rangeProcessing')
            wrg = float(rg_proc.find('windowCoefficient').text)
            brg = float(rg_proc.find('processingBandwidth').text)
            lbrg = float(rg_proc.find('lookBandwidth').text)
            lrg = brg / lbrg
            kbrg = b_factors[coefficients.index(wrg)]
            resolutions_rg.append(0.886 * c / (2 * brg) * kbrg * lrg)
            
            # computation of azimuth resolution; yet to be checked for correctness
            az_proc = par.find('azimuthProcessing')
            waz = float(az_proc.find('windowCoefficient').text)
            baz = float(az_proc.find('processingBandwidth').text)
            lbaz = float(az_proc.find('lookBandwidth').text)
            laz = baz / lbaz
            kbaz = b_factors[coefficients.index(waz)]
            vsat = sp_az[i] / ti_az[i]
            resolutions_az.append(0.886 * vsat / baz * kbaz * laz)
        
        resolution_rg = median(resolutions_rg)
        resolution_az = median(resolutions_az)
        
        if self.meta['image_geometry'] == 'GROUND_RANGE':
            resolution_rg /= math.sin(math.radians(self.meta['incidence']))
        
        self.meta['resolution'] = resolution_rg, resolution_az
        return self.meta['resolution']
    
    def scanMetadata(self):
        with self.getFileObj(self.findfiles('manifest.safe')[0]) as input:
            manifest = input.getvalue()
        namespaces = getNamespaces(manifest)
        tree = ET.fromstring(manifest)
        
        meta = dict()
        key = 's1sarl1'
        obj_prod = tree.find('.//{}:productType'.format(key), namespaces)
        if obj_prod == None:
            key = 's1sarl2'
            obj_prod = tree.find('.//{}:productType'.format(key), namespaces)
        
        meta['product'] = obj_prod.text
        
        acqmode = tree.find('.//{}:mode'.format(key), namespaces).text
        if acqmode == 'SM':
            meta['acquisition_mode'] = tree.find('.//{}:swath'.format(key), namespaces).text
        else:
            meta['acquisition_mode'] = acqmode
        meta['acquisition_time'] = dict(
            [(x, tree.find('.//safe:{}Time'.format(x), namespaces).text) for x in ['start', 'stop']])
        meta['start'], meta['stop'] = (self.parse_date(meta['acquisition_time'][x]) for x in ['start', 'stop'])
        meta['coordinates'] = [tuple([float(y) for y in x.split(',')][::-1]) for x in
                               tree.find('.//gml:coordinates', namespaces).text.split()]
        meta['orbit'] = tree.find('.//s1:pass', namespaces).text[0]
        
        meta['orbitNumber_abs'] = int(tree.find('.//safe:orbitNumber[@type="start"]', namespaces).text)
        meta['orbitNumber_rel'] = int(tree.find('.//safe:relativeOrbitNumber[@type="start"]', namespaces).text)
        meta['cycleNumber'] = int(tree.find('.//safe:cycleNumber', namespaces).text)
        meta['frameNumber'] = int(tree.find('.//{}:missionDataTakeID'.format(key), namespaces).text)
        
        meta['orbitNumbers_abs'] = dict(
            [(x, int(tree.find('.//safe:orbitNumber[@type="{0}"]'.format(x), namespaces).text)) for x in
             ['start', 'stop']])
        meta['orbitNumbers_rel'] = dict(
            [(x, int(tree.find('.//safe:relativeOrbitNumber[@type="{0}"]'.format(x), namespaces).text)) for x in
             ['start', 'stop']])
        key_pol = './/{}:transmitterReceiverPolarisation'.format(key)
        meta['polarizations'] = [x.text for x in tree.findall(key_pol, namespaces)]
        meta['category'] = tree.find('.//{}:productClass'.format(key), namespaces).text
        family = tree.find('.//safe:familyName', namespaces).text.replace('ENTINEL-', '')
        number = tree.find('.//safe:number', namespaces).text
        meta['sensor'] = family + number
        meta['IPF_version'] = float(tree.find('.//safe:software', namespaces).attrib['version'])
        sliced = tree.find('.//{}:sliceProductFlag'.format(key), namespaces).text == 'true'
        if sliced:
            meta['sliceNumber'] = int(tree.find('.//{}:sliceNumber'.format(key), namespaces).text)
            meta['totalSlices'] = int(tree.find('.//{}:totalSlices'.format(key), namespaces).text)
        else:
            meta['sliceNumber'] = None
            meta['totalSlices'] = None
        
        if meta['product'] == 'OCN':
            meta['spacing'] = -1
            meta['samples'] = -1
            meta['lines'] = -1
        else:
            annotations = self.findfiles(self.pattern_ds)
            key = lambda x: re.search('-[vh]{2}-', x).group()
            groups = groupby(sorted(annotations, key=key), key=key)
            annotations = [list(value) for key, value in groups][0]
            ann_trees = []
            for ann in annotations:
                with self.getFileObj(ann) as ann_xml:
                    ann_trees.append(ET.fromstring(ann_xml.read()))
            
            sp_rg = [float(x.find('.//rangePixelSpacing').text) for x in ann_trees]
            sp_az = [float(x.find('.//azimuthPixelSpacing').text) for x in ann_trees]
            meta['spacing'] = (median(sp_rg), median(sp_az))
            
            looks_rg = [float(x.find('.//rangeProcessing/numberOfLooks').text) for x in ann_trees]
            looks_az = [float(x.find('.//azimuthProcessing/numberOfLooks').text) for x in ann_trees]
            meta['looks'] = (median(looks_rg), median(looks_az))
            
            samples = [x.find('.//imageAnnotation/imageInformation/numberOfSamples').text for x in ann_trees]
            meta['samples'] = sum([int(x) for x in samples])
            
            lines = [x.find('.//imageAnnotation/imageInformation/numberOfLines').text for x in ann_trees]
            meta['lines'] = sum([int(x) for x in lines])
            
            heading = median(float(x.find('.//platformHeading').text) for x in ann_trees)
            meta['heading'] = heading if heading > 0 else heading + 360
            
            incidence = [float(x.find('.//incidenceAngleMidSwath').text) for x in ann_trees]
            meta['incidence'] = median(incidence)
            
            meta['image_geometry'] = ann_trees[0].find('.//projection').text.replace(' ', '_').upper()
        
        return meta
    
    def unpack(self, directory, overwrite=False, exist_ok=False):
        outdir = os.path.join(directory, os.path.basename(self.file))
        self._unpack(outdir, overwrite=overwrite, exist_ok=exist_ok)


class TSX(ID):
    """
    Handler class for TerraSAR-X and TanDEM-X data
    
    Sensors:
        * TSX1
        * TDX1

    References:
        * TX-GS-DD-3302  TerraSAR-X Basic Product Specification Document
        * TX-GS-DD-3303  TerraSAR-X Experimental Product Description
        * TD-GS-PS-3028  TanDEM-X Experimental Product Description
        * TerraSAR-X Image Product Guide (Airbus Defence and Space)
    
    Acquisition modes:
        * ST:    Staring Spotlight
        * HS:    High Resolution SpotLight
        * HS300: High Resolution SpotLight 300 MHz
        * SL:    SpotLight
        * SM:    StripMap
        * SC:    ScanSAR
        * WS:    Wide ScanSAR
    
    Polarisation modes:
        * Single (S): all acquisition modes
        * Dual   (D): High Resolution SpotLight (HS), SpotLight (SL) and StripMap (SM)
        * Twin   (T): StripMap (SM) (experimental)
        * Quad   (Q): StripMap (SM) (experimental)
    
    Products:
        * SSC: Single Look Slant Range Complex
        * MGD: Multi Look Ground Range Detected
        * GEC: Geocoded Ellipsoid Corrected
        * EEC: Enhanced Ellipsoid Corrected
    """
    
    def __init__(self, scene):
        if isinstance(scene, str):
            self.scene = os.path.realpath(scene)
            
            self.pattern = patterns.tsx
            
            self.pattern_ds = r'^IMAGE_(?P<pol>HH|HV|VH|VV)_(?:SRA|FWD|AFT)_(?P<beam>[^\.]+)\.(cos|tif)$'
            self.examine(include_folders=False)
            
            if not re.match(re.compile(self.pattern), os.path.basename(self.file)):
                raise RuntimeError('folder does not match TSX scene naming convention')
            
            self.meta = self.scanMetadata()
            self.meta['projection'] = crsConvert(4326, 'wkt')
        
        super(TSX, self).__init__(self.meta)
    
    def scanMetadata(self):
        annotation = self.getFileObj(self.file).getvalue()
        namespaces = getNamespaces(annotation)
        tree = ET.fromstring(annotation)
        meta = dict()
        meta['sensor'] = tree.find('.//generalHeader/mission', namespaces).text.replace('-', '')
        meta['product'] = tree.find('.//orderInfo/productVariant', namespaces).text
        meta['orbit'] = tree.find('.//missionInfo/orbitDirection', namespaces).text[0]
        meta['polarizations'] = [x.text for x in
                                 tree.findall('.//acquisitionInfo/polarisationList/polLayer', namespaces)]
        
        meta['orbitNumber_abs'] = int(tree.find('.//missionInfo/absOrbit', namespaces).text)
        meta['orbitNumber_rel'] = int(tree.find('.//missionInfo/relOrbit', namespaces).text)
        meta['cycleNumber'] = int(tree.find('.//missionInfo/orbitCycle', namespaces).text)
        meta['frameNumber'] = int(tree.find('.//inputData/uniqueDataTakeID', namespaces).text)
        
        meta['acquisition_mode'] = tree.find('.//acquisitionInfo/imagingMode', namespaces).text
        meta['start'] = self.parse_date(tree.find('.//sceneInfo/start/timeUTC', namespaces).text)
        meta['stop'] = self.parse_date(tree.find('.//sceneInfo/stop/timeUTC', namespaces).text)
        spacing_row = float(tree.find('.//imageDataInfo/imageRaster/rowSpacing', namespaces).text)
        spacing_col = float(tree.find('.//imageDataInfo/imageRaster/columnSpacing', namespaces).text)
        meta['spacing'] = (spacing_col, spacing_row)
        meta['samples'] = int(tree.find('.//imageDataInfo/imageRaster/numberOfColumns', namespaces).text)
        meta['lines'] = int(tree.find('.//imageDataInfo/imageRaster/numberOfRows', namespaces).text)
        rlks = float(tree.find('.//imageDataInfo/imageRaster/rangeLooks', namespaces).text)
        azlks = float(tree.find('.//imageDataInfo/imageRaster/azimuthLooks', namespaces).text)
        meta['looks'] = (rlks, azlks)
        meta['incidence'] = float(tree.find('.//sceneInfo/sceneCenterCoord/incidenceAngle', namespaces).text)
        
        geocs = self.getFileObj(self.findfiles('GEOREF.xml')[0]).getvalue()
        tree = ET.fromstring(geocs)
        pts = tree.findall('.//gridPoint')
        lat = [float(x.find('lat').text) for x in pts]
        lon = [float(x.find('lon').text) for x in pts]
        # shift lon in case of west direction.
        lon = [x - 360 if x > 180 else x for x in lon]
        meta['coordinates'] = list(zip(lon, lat))
        
        return meta
    
    def unpack(self, directory, overwrite=False, exist_ok=False):
        match = self.findfiles(self.pattern, True)
        header = [x for x in match if not x.endswith('xml') and 'iif' not in x][0].replace(self.scene, '').strip('/')
        outdir = os.path.join(directory, os.path.basename(header))
        self._unpack(outdir, offset=header, overwrite=overwrite, exist_ok=exist_ok)


class TDM(TSX):
    """
    Handler class for TerraSAR-X and TanDEM-X experimental data
    
    Sensors:
        * TDM1

    References:
        * TD-GS-PS-3028  TanDEM-X Experimental Product Description
    
    Acquisition modes:
        * HS:    High Resolution SpotLight
        * SL:    SpotLight
        * SM:    StripMap
    
    Polarisation modes:
        * Single (S): all acquisition modes
        * Dual   (D): High Resolution SpotLight (HS), SpotLight (SL) and StripMap (SM)
        * Twin   (T): StripMap (SM) (experimental)
        * Quad   (Q): StripMap (SM) (experimental)
    
    Products:
        * CoSSCs: (bi-static) SAR co-registered single look slant range complex products (CoSSCs)


    Examples
    ----------
    Ingest all Tandem-X Bistatic scenes in a directory and its sub-directories into the database:

    >>> from pyroSAR import Archive, identify
    >>> from spatialist.ancillary import finder
    >>> dbfile = '/.../scenelist.db'
    >>> archive_tdm = '/.../TDM/'
    >>> scenes_tdm = finder(archive_tdm, [r'^TDM1.*'], foldermode=2, regex=True, recursive=True)
    >>> with Archive(dbfile) as archive:
    >>>     archive.insert(scenes_tdm)
    """
    
    def __init__(self, scene):
        self.scene = os.path.realpath(scene)
        
        self.pattern = patterns.tdm
        
        self.pattern_ds = r'^IMAGE_(?P<pol>HH|HV|VH|VV)_(?:SRA|FWD|AFT)_(?P<beam>[^\.]+)\.(cos|tif)$'
        self.examine(include_folders=False)
        
        if not re.match(re.compile(self.pattern), os.path.basename(self.file)):
            raise RuntimeError('folder does not match TDM scene naming convention')
        
        self.meta = self.scanMetadata()
        self.meta['projection'] = crsConvert(4326, 'wkt')
        
        super(TDM, self).__init__(self.meta)
    
    def scanMetadata(self):
        annotation = self.getFileObj(self.file).getvalue()
        namespaces = getNamespaces(annotation)
        tree = ET.fromstring(annotation)
        meta = dict()
        meta['sensor'] = tree.find('.//commonAcquisitionInfo/missionID', namespaces).text.replace('-', '')
        meta['product'] = tree.find('.//productInfo/productType', namespaces).text
        meta['SAT1'] = tree.find('.//commonAcquisitionInfo/satelliteIDsat1', namespaces).text
        meta['SAT2'] = tree.find('.//commonAcquisitionInfo/satelliteIDsat2', namespaces).text
        meta['inSARmasterID'] = tree.find('.//commonAcquisitionInfo/inSARmasterID', namespaces).text
        pattern = './/commonAcquisitionInfo/satelliteID{}'.format(meta['inSARmasterID'].lower())
        meta['inSARmaster'] = tree.find(pattern, namespaces).text.replace('-', '')
        
        pattern = './/commonAcquisitionInfo/operationsInfo/acquisitionItemID'
        meta['acquisitionItemID'] = int(tree.find(pattern, namespaces).text)
        
        meta['effectiveBaseline'] = float(tree.find('.//acquisitionGeometry/effectiveBaseline', namespaces).text)
        meta['heightOfAmbiguity'] = float(tree.find('.//acquisitionGeometry/heightOfAmbiguity', namespaces).text)
        meta['distanceActivePos'] = float(tree.find('.//acquisitionGeometry/distanceActivePos', namespaces).text)
        meta['distanceTracks'] = float(tree.find('.//acquisitionGeometry/distanceTracks', namespaces).text)
        
        meta['cooperativeMode'] = tree.find('.//commonAcquisitionInfo/cooperativeMode', namespaces).text
        
        if meta['cooperativeMode'].lower() == "bistatic":
            meta['bistatic'] = True
        else:
            meta['bistatic'] = False
        
        meta['orbit'] = tree.find('.//acquisitionGeometry/orbitDirection', namespaces).text[0]
        
        pattern = ".//productComponents/component[@componentClass='imageData']/file/location/name"
        elements = tree.findall(pattern, )
        self.primary_scene = os.path.join(self.scene, elements[0].text)
        self.secondary_scene = os.path.join(self.scene, elements[1].text)
        meta["SAT1"] = TSX(self.primary_scene).scanMetadata()
        meta["SAT2"] = TSX(self.secondary_scene).scanMetadata()
        
        meta['start'] = self.parse_date(tree.find('.//orbitHeader/firstStateTime/firstStateTimeUTC', namespaces).text)
        meta['stop'] = self.parse_date(tree.find('.//orbitHeader/lastStateTime/lastStateTimeUTC', namespaces).text)
        meta['samples'] = int(tree.find('.//coregistration/coregRaster/samples', namespaces).text)
        meta['lines'] = int(tree.find('.//coregistration/coregRaster/lines', namespaces).text)
        rlks = float(tree.find('.//processingInfo/inSARProcessing/looks/range', namespaces).text)
        azlks = float(tree.find('.//processingInfo/inSARProcessing/looks/azimuth', namespaces).text)
        meta['looks'] = (rlks, azlks)
        meta['incidence'] = float(tree.find('.//commonSceneInfo/sceneCenterCoord/incidenceAngle', namespaces).text)
        
        meta['orbit'] = meta[meta['inSARmasterID']]['orbit']
        meta['polarizations'] = meta[meta['inSARmasterID']]['polarizations']
        
        meta['orbitNumber_abs'] = meta[meta['inSARmasterID']]['orbitNumber_abs']
        meta['orbitNumber_rel'] = meta[meta['inSARmasterID']]['orbitNumber_rel']
        meta['cycleNumber'] = meta[meta['inSARmasterID']]['cycleNumber']
        meta['frameNumber'] = meta[meta['inSARmasterID']]['frameNumber']
        
        meta['acquisition_mode'] = meta[meta['inSARmasterID']]['acquisition_mode']
        meta['start'] = meta[meta['inSARmasterID']]['start']
        meta['stop'] = meta[meta['inSARmasterID']]['stop']
        meta['spacing'] = meta[meta['inSARmasterID']]['spacing']
        meta['samples'] = meta[meta['inSARmasterID']]['samples']
        meta['lines'] = meta[meta['inSARmasterID']]['lines']
        meta['looks'] = meta[meta['inSARmasterID']]['looks']
        meta['incidence'] = meta[meta['inSARmasterID']]['incidence']
        
        pts = tree.findall('.//sceneCornerCoord')
        lat = [float(x.find('lat').text) for x in pts]
        lon = [float(x.find('lon').text) for x in pts]
        # shift lon in case of west direction.
        lon = [x - 360 if x > 180 else x for x in lon]
        meta['coordinates'] = list(zip(lon, lat))
        
        return meta


class Archive(object):
    """
    Utility for storing SAR image metadata in a database

    Parameters
    ----------
    dbfile: str
        the filename for the SpatiaLite database. This might either point to an existing database or will be created otherwise.
        If postgres is set to True, this will be the name for the PostgreSQL database.
    custom_fields: dict or None
        a dictionary containing additional non-standard database column names and data types;
        the names must be attributes of the SAR scenes to be inserted (i.e. id.attr) or keys in their meta attribute
        (i.e. id.meta['attr'])
    postgres: bool
        enable postgres driver for the database. Default: False
    user: str
        required for postgres driver: username to access the database. Default: 'postgres'
    password: str
        required for postgres driver: password to access the database. Default: '1234'
    host: str
        required for postgres driver: host where the database is hosted. Default: 'localhost'
    port: int
        required for postgres driver: port number to the database. Default: 5432
    cleanup: bool
        check whether all registered scenes exist and remove missing entries?
    legacy: bool
        open an outdated database in legacy mode to import into a new database.
        Opening an outdated database without legacy mode will throw a RuntimeError.

    Examples
    ----------
    Ingest all Sentinel-1 scenes in a directory and its subdirectories into the database:

    >>> from pyroSAR import Archive, identify
    >>> from spatialist.ancillary import finder
    >>> dbfile = '/.../scenelist.db'
    >>> archive_s1 = '/.../sentinel1/GRD'
    >>> scenes_s1 = finder(archive_s1, [r'^S1[AB].*.zip'], regex=True, recursive=True)
    >>> with Archive(dbfile) as archive:
    >>>     archive.insert(scenes_s1)

    select all Sentinel-1 A/B scenes stored in the database, which
    
     * overlap with a test site
     * were acquired in Ground-Range-Detected (GRD) Interferometric Wide Swath (IW) mode before 2018
     * contain a VV polarization image
     * have not been processed to directory `outdir` before

    >>> from pyroSAR import Archive
    >>> from spatialist import Vector
    >>> archive = Archive('/.../scenelist.db')
    >>> site = Vector('/path/to/site.shp')
    >>> outdir = '/path/to/processed/results'
    >>> maxdate = '20171231T235959'
    >>> selection_proc = archive.select(vectorobject=site, processdir=outdir,
    >>>                                 maxdate=maxdate, sensor=('S1A', 'S1B'),
    >>>                                 product='GRD', acquisition_mode='IW', vv=1)
    >>> archive.close()

    Alternatively, the `with` statement can be used.
    In this case to just check whether one particular scene is already registered in the database:

    >>> from pyroSAR import identify, Archive
    >>> scene = identify('S1A_IW_SLC__1SDV_20150330T170734_20150330T170801_005264_006A6C_DA69.zip')
    >>> with Archive('/.../scenelist.db') as archive:
    >>>     print(archive.is_registered(scene.scene))

    When providing 'postgres' as driver, a PostgreSQL database will be created at a given host.
    Additional arguments are required.

    >>> from pyroSAR import Archive, identify
    >>> from spatialist.ancillary import finder
    >>> dbfile = 'scenelist_db'
    >>> archive_s1 = '/.../sentinel1/GRD'
    >>> scenes_s1 = finder(archive_s1, [r'^S1[AB].*.zip'], regex=True, recursive=True)
    >>> with Archive(dbfile, driver='postgres', user='user', password='password', host='host', port=5432) as archive:
    >>>     archive.insert(scenes_s1)
    
    Importing an old database:
    
    >>> from pyroSAR import Archive
    >>> db_new = 'scenes.db'
    >>> db_old = 'scenes_old.db'
    >>> with Archive(db_new) as db:
    >>>     with Archive(db_old, legacy=True) as db_old:
    >>>         db.import_outdated(db_old)
    """
    
    def __init__(self, dbfile, custom_fields=None, postgres=False, user='postgres',
                 password='1234', host='localhost', port=5432, cleanup=True,
                 legacy=False):
        
        if dbfile.endswith('.csv'):
            raise RuntimeError("Please create a new Archive database and import the"
                               "CSV file using db.import_outdated('<file>.csv').")
        # check for driver, if postgres then check if server is reachable
        if not postgres:
            self.driver = 'sqlite'
            dirname = os.path.dirname(os.path.abspath(dbfile))
            w_ok = os.access(dirname, os.W_OK)
            if not w_ok:
                raise RuntimeError('cannot write to directory {}'.format(dirname))
            # catch if .db extension is missing
            root, ext = os.path.splitext(dbfile)
            if len(ext) == 0:
                dbfile = root + '.db'
        else:
            self.driver = 'postgresql'
            if not self.__check_host(host, port):
                sys.exit('Server not found!')
        
        connect_args = {}
        
        # create dict, with which a URL to the db is created
        if self.driver == 'sqlite':
            self.url_dict = {'drivername': self.driver,
                             'database': dbfile,
                             'query': {'charset': 'utf8'}}
        if self.driver == 'postgresql':
            self.url_dict = {'drivername': self.driver,
                             'username': user,
                             'password': password,
                             'host': host,
                             'port': port,
                             'database': dbfile}
            connect_args = {
                'keepalives': 1,
                'keepalives_idle': 30,
                'keepalives_interval': 10,
                'keepalives_count': 5}
        
        # create engine, containing URL and driver
        log.debug('starting DB engine for {}'.format(URL.create(**self.url_dict)))
        self.url = URL.create(**self.url_dict)
        # https://www.postgresql.org/docs/current/libpq-connect.html#LIBPQ-PARAMKEYWORDS
        self.engine = create_engine(url=self.url, echo=False,
                                    connect_args=connect_args)
        
        # call to __load_spatialite() for sqlite, to load mod_spatialite via event handler listen()
        if self.driver == 'sqlite':
            log.debug('loading spatialite extension')
            listen(target=self.engine, identifier='connect', fn=self.__load_spatialite)
            # check if loading was successful
            try:
                with self.engine.begin() as conn:
                    version = conn.execute('SELECT spatialite_version();')
            except exc.OperationalError:
                raise RuntimeError('could not load spatialite extension')
        
        # if database is new, (create postgres-db and) enable spatial extension
        if not database_exists(self.engine.url):
            if self.driver == 'postgresql':
                log.debug('creating new PostgreSQL database')
                create_database(self.engine.url)
            log.debug('enabling spatial extension for new database')
            with self.engine.begin() as conn:
                if self.driver == 'sqlite':
                    conn.execute(select([func.InitSpatialMetaData(1)]))
                else:
                    conn.exec_driver_sql('CREATE EXTENSION IF NOT EXISTS postgis;')
        # create Session (ORM) and get metadata
        self.Session = sessionmaker(bind=self.engine)
        self.meta = MetaData(self.engine)
        self.custom_fields = custom_fields
        
        # load or create tables
        self.__init_data_table()
        self.__init_duplicates_table()
        
        msg = ("the 'data' table is missing {}. Please create a new database "
               "and import the old one opened in legacy mode using "
               "Archive.import_outdated.")
        pk = sql_inspect(self.data_schema).primary_key
        if 'product' not in pk.columns.keys() and not legacy:
            raise RuntimeError(msg.format("a primary key 'product'"))
        
        if 'geometry' not in self.get_colnames() and not legacy:
            raise RuntimeError(msg.format("the 'geometry' column"))
        
        self.Base = automap_base(metadata=self.meta)
        self.Base.prepare(self.engine, reflect=True)
        self.Data = self.Base.classes.data
        self.Duplicates = self.Base.classes.duplicates
        self.dbfile = dbfile
        
        if cleanup:
            log.info('checking for missing scenes')
            self.cleanup()
            sys.stdout.flush()
    
    def add_tables(self, tables):
        """
        Add tables to the database per :class:`sqlalchemy.schema.Table`
        Tables provided here will be added to the database.
        
        .. note::
        
            Columns using Geometry must have setting management=True for SQLite,
            for example: ``geometry = Column(Geometry('POLYGON', management=True, srid=4326))``
        
        Parameters
        ----------
        tables: :class:`sqlalchemy.schema.Table` or list[:class:`sqlalchemy.schema.Table`]
            The table(s) to be added to the database.
        """
        created = []
        if isinstance(tables, list):
            for table in tables:
                table.metadata = self.meta
                if not sql_inspect(self.engine).has_table(str(table)):
                    table.create(self.engine)
                    created.append(str(table))
        else:
            table = tables
            table.metadata = self.meta
            if not sql_inspect(self.engine).has_table(str(table)):
                table.create(self.engine)
                created.append(str(table))
        log.info('created table(s) {}.'.format(', '.join(created)))
        self.Base = automap_base(metadata=self.meta)
        self.Base.prepare(self.engine, reflect=True)
    
    def __init_data_table(self):
        if sql_inspect(self.engine).has_table('data'):
            self.data_schema = Table('data', self.meta, autoload_with=self.engine)
            return
        
        log.debug("creating DB table 'data'")
        
        self.data_schema = Table('data', self.meta,
                                 Column('sensor', String),
                                 Column('orbit', String),
                                 Column('orbitNumber_abs', Integer),
                                 Column('orbitNumber_rel', Integer),
                                 Column('cycleNumber', Integer),
                                 Column('frameNumber', Integer),
                                 Column('acquisition_mode', String),
                                 Column('start', String),
                                 Column('stop', String),
                                 Column('product', String, primary_key=True),
                                 Column('samples', Integer),
                                 Column('lines', Integer),
                                 Column('outname_base', String, primary_key=True),
                                 Column('scene', String),
                                 Column('hh', Integer),
                                 Column('vv', Integer),
                                 Column('hv', Integer),
                                 Column('vh', Integer),
                                 Column('geometry', Geometry(geometry_type='POLYGON',
                                                             management=True, srid=4326)))
        # add custom fields
        if self.custom_fields is not None:
            for key, val in self.custom_fields.items():
                if val in ['Integer', 'integer', 'int']:
                    self.data_schema.append_column(Column(key, Integer))
                elif val in ['String', 'string', 'str']:
                    self.data_schema.append_column(Column(key, String))
                else:
                    log.info('Value in dict custom_fields must be "integer" or "string"!')
        
        self.data_schema.create(self.engine)
    
    def __init_duplicates_table(self):
        # create tables if not existing
        if sql_inspect(self.engine).has_table('duplicates'):
            self.duplicates_schema = Table('duplicates', self.meta, autoload_with=self.engine)
            return
        
        log.debug("creating DB table 'duplicates'")
        
        self.duplicates_schema = Table('duplicates', self.meta,
                                       Column('outname_base', String, primary_key=True),
                                       Column('scene', String, primary_key=True))
        self.duplicates_schema.create(self.engine)
    
    @staticmethod
    def __load_spatialite(dbapi_conn, connection_record):
        """
        loads the spatialite extension for SQLite, not to be used outside the init()
        
        Parameters
        ----------
        dbapi_conn:
            db engine
        connection_record:
            not sure what it does, but it is needed by :func:`sqlalchemy.event.listen`
        """
        dbapi_conn.enable_load_extension(True)
        # check which platform and use according mod_spatialite
        if platform.system() == 'Linux':
            for option in ['mod_spatialite', 'mod_spatialite.so']:
                try:
                    dbapi_conn.load_extension(option)
                except sqlite3.OperationalError:
                    continue
        elif platform.system() == 'Darwin':
            for option in ['mod_spatialite.so', 'mod_spatialite.7.dylib',
                           'mod_spatialite.dylib']:
                try:
                    dbapi_conn.load_extension(option)
                except sqlite3.OperationalError:
                    continue
        else:
            dbapi_conn.load_extension('mod_spatialite')
    
    def __prepare_insertion(self, scene):
        """
        read scene metadata and parse a string for inserting it into the database

        Parameters
        ----------
        scene: str or ID
            a SAR scene

        Returns
        -------
        object of class Data, insert string
        """
        id = scene if isinstance(scene, ID) else identify(scene)
        pols = [x.lower() for x in id.polarizations]
        # insertion as an object of Class Data (reflected in the init())
        insertion = self.Data()
        colnames = self.get_colnames()
        for attribute in colnames:
            if attribute == 'geometry':
                geom = id.geometry()
                geom.reproject(4326)
                geom = geom.convert2wkt(set3D=False)[0]
                geom = 'SRID=4326;' + str(geom)
                # set attributes of the Data object according to input
                setattr(insertion, 'geometry', geom)
            elif attribute in ['hh', 'vv', 'hv', 'vh']:
                setattr(insertion, attribute, int(attribute in pols))
            else:
                if hasattr(id, attribute):
                    attr = getattr(id, attribute)
                elif attribute in id.meta.keys():
                    attr = id.meta[attribute]
                else:
                    raise AttributeError('could not find attribute {}'.format(attribute))
                value = attr() if inspect.ismethod(attr) else attr
                setattr(insertion, str(attribute), value)
        
        return insertion  # return the Data object
    
    def __select_missing(self, table):
        """

        Returns
        -------
        list[str]
            the names of all scenes, which are no longer stored in their registered location
        """
        with self.Session() as session:
            if table == 'data':
                # using ORM query to get all scenes locations
                scenes = session.query(self.Data.scene)
            elif table == 'duplicates':
                scenes = session.query(self.Duplicates.scene)
            else:
                raise ValueError("parameter 'table' must either be 'data' or 'duplicates'")
        files = [self.encode(x[0]) for x in scenes]
        return [x for x in files if not os.path.isfile(x)]
    
    def insert(self, scene_in, pbar=False, test=False):
        """
        Insert one or many scenes into the database

        Parameters
        ----------
        scene_in: str or ID or list[str or ID]
            a SAR scene or a list of scenes to be inserted
        pbar: bool
            show a progress bar?
        test: bool
            should the insertion only be tested or directly be committed to the database?
        """
        
        if isinstance(scene_in, (ID, str)):
            scene_in = [scene_in]
        if not isinstance(scene_in, list):
            raise RuntimeError('scene_in must either be a string pointing to a file, a pyroSAR.ID object '
                               'or a list containing several of either')
        
        log.info('filtering scenes by name')
        scenes = self.filter_scenelist(scene_in)
        if len(scenes) == 0:
            log.info('...nothing to be done')
            return
        log.info('identifying scenes and extracting metadata')
        scenes = identify_many(scenes, pbar=pbar)
        
        if len(scenes) == 0:
            log.info('all scenes are already registered')
            return
        
        counter_regulars = 0
        counter_duplicates = 0
        list_duplicates = []
        
        message = 'inserting {0} scene{1} into database'
        log.info(message.format(len(scenes), '' if len(scenes) == 1 else 's'))
        log.debug('testing changes in temporary database')
        if pbar:
            progress = pb.ProgressBar(max_value=len(scenes))
        else:
            progress = None
        insertions = []
        with self.Session() as session:
            for i, id in enumerate(scenes):
                basename = id.outname_base()
                if not self.is_registered(id):
                    insertion = self.__prepare_insertion(id)
                    insertions.append(insertion)
                    counter_regulars += 1
                    log.debug('regular:   {}'.format(id.scene))
                elif not self.__is_registered_in_duplicates(id):
                    insertion = self.Duplicates(outname_base=basename,
                                                scene=id.scene)
                    insertions.append(insertion)
                    counter_duplicates += 1
                    log.debug('duplicate: {}'.format(id.scene))
                else:
                    list_duplicates.append(id.outname_base())
                
                if progress is not None:
                    progress.update(i + 1)
            
            if progress is not None:
                progress.finish()
            
            session.add_all(insertions)
            
            if not test:
                log.debug('committing transactions to permanent database')
                # commit changes of the session
                session.commit()
            else:
                log.info('rolling back temporary database changes')
                # roll back changes of the session
                session.rollback()
        
        message = '{0} scene{1} registered regularly'
        log.info(message.format(counter_regulars, '' if counter_regulars == 1 else 's'))
        message = '{0} duplicate{1} registered'
        log.info(message.format(counter_duplicates, '' if counter_duplicates == 1 else 's'))
    
    def is_registered(self, scene):
        """
        Simple check if a scene is already registered in the database.

        Parameters
        ----------
        scene: str or ID
            the SAR scene

        Returns
        -------
        bool
            is the scene already registered?
        """
        id = scene if isinstance(scene, ID) else identify(scene)
        with self.Session() as session:
            # ORM query, where scene equals id.scene, return first
            exists_data = session.query(self.Data.outname_base).filter_by(
                outname_base=id.outname_base(), product=id.product).first()
            exists_duplicates = session.query(self.Duplicates.outname_base).filter(
                self.Duplicates.outname_base == id.outname_base()).first()
        in_data = False
        in_dup = False
        if exists_data:
            in_data = len(exists_data) != 0
        if exists_duplicates:
            in_dup = len(exists_duplicates) != 0
        return in_data or in_dup
    
    def __is_registered_in_duplicates(self, scene):
        """
        Simple check if a scene is already registered in the database.

        Parameters
        ----------
        scene: str or ID
            the SAR scene

        Returns
        -------
        bool
            is the scene already registered?
        """
        id = scene if isinstance(scene, ID) else identify(scene)
        with self.Session() as session:
            # ORM query as in is registered
            exists_duplicates = session.query(self.Duplicates.outname_base).filter(
                self.Duplicates.outname_base == id.outname_base()).first()
        in_dup = False
        if exists_duplicates:
            in_dup = len(exists_duplicates) != 0
        return in_dup
    
    def cleanup(self):
        """
        Remove all scenes from the database, which are no longer stored in their registered location

        Returns
        -------

        """
        missing = self.__select_missing('data')
        for scene in missing:
            log.info('Removing missing scene from database tables: {}'.format(scene))
            self.drop_element(scene, with_duplicates=True)
    
    @staticmethod
    def encode(string, encoding='utf-8'):
        if not isinstance(string, str) and hasattr(string, 'encode'):
            return string.encode(encoding)
        else:
            return string
    
    def export2shp(self, path, table='data'):
        """
        export the database to a shapefile

        Parameters
        ----------
        path: str
            the path of the shapefile to be written.
            This will overwrite other files with the same name.
            If a folder is given in path it is created if not existing.
            If the file extension is missing '.shp' is added.
        table: str
            the table to write to the shapefile; either 'data' (default) or 'duplicates'
        
        Returns
        -------
        """
        if table not in self.get_tablenames():
            log.warning('Only data and duplicates can be exported!')
            return
        
        # add the .shp extension if missing
        if not path.endswith('.shp'):
            path += '.shp'
        
        # creates folder if not present, adds .shp if not within the path
        dirname = os.path.dirname(path)
        os.makedirs(dirname, exist_ok=True)
        
        launder_names = {'acquisition_mode': 'acq_mode',
                         'orbitNumber_abs': 'orbit_abs',
                         'orbitNumber_rel': 'orbit_rel',
                         'cycleNumber': 'cycleNr',
                         'frameNumber': 'frameNr',
                         'outname_base': 'outname'}
        
        sel_tables = ', '.join([f'"{s}" as {launder_names[s]}' if s in launder_names else s
                                for s in self.get_colnames(table)])
        
        if self.driver == 'sqlite':
            srcDS = self.dbfile
        elif self.driver == 'postgresql':
            srcDS = """PG:host={host} port={port} user={username}
            dbname={database} password={password} active_schema=public""".format(**self.url_dict)
        else:
            raise RuntimeError('unknown archive driver')
        
        gdal.VectorTranslate(destNameOrDestDS=path, srcDS=srcDS,
                             format='ESRI Shapefile',
                             SQLStatement=f'SELECT {sel_tables} FROM {table}',
                             SQLDialect=self.driver)
    
    def filter_scenelist(self, scenelist):
        """
        Filter a list of scenes by file names already registered in the database.

        Parameters
        ----------
        scenelist: list[str or ID]
            the scenes to be filtered

        Returns
        -------
        list[ID]
            the file names of the scenes whose basename is not yet registered in the database

        """
        for item in scenelist:
            if not isinstance(item, (ID, str)):
                raise TypeError("items in scenelist must be of type 'str' or 'pyroSAR.ID'")
        
        with self.Session() as session:
            # ORM query, get all scenes locations
            scenes_data = session.query(self.Data.scene)
            registered = [os.path.basename(self.encode(x[0])) for x in scenes_data]
            scenes_duplicates = session.query(self.Duplicates.scene)
        duplicates = [os.path.basename(self.encode(x[0])) for x in scenes_duplicates]
        names = [item.scene if isinstance(item, ID) else item for item in scenelist]
        filtered = [x for x, y in zip(scenelist, names) if os.path.basename(y) not in registered + duplicates]
        return filtered
    
    def get_colnames(self, table='data'):
        """
        Return the names of all columns of a table.

        Returns
        -------
        list[str]
            the column names of the chosen table
        """
        # get all columns of one table, but shows geometry columns not correctly
        table_info = Table(table, self.meta, autoload=True, autoload_with=self.engine)
        col_names = table_info.c.keys()
        
        return sorted([self.encode(x) for x in col_names])
    
    def get_tablenames(self, return_all=False):
        """
        Return the names of all tables in the database
        
        Parameters
        ----------
        return_all: bool
            only gives tables data and duplicates on default.
            Set to True to get all other tables and views created automatically.

        Returns
        -------
        list[str]
            the table names
        """
        #  TODO: make this dynamic
        #  the method was intended to only return user generated tables by default, as well as data and duplicates
        all_tables = ['ElementaryGeometries', 'SpatialIndex', 'geometry_columns', 'geometry_columns_auth',
                      'geometry_columns_field_infos', 'geometry_columns_statistics', 'geometry_columns_time',
                      'spatial_ref_sys', 'spatial_ref_sys_aux', 'spatialite_history', 'sql_statements_log',
                      'sqlite_sequence', 'views_geometry_columns', 'views_geometry_columns_auth',
                      'views_geometry_columns_field_infos', 'views_geometry_columns_statistics',
                      'virts_geometry_columns', 'virts_geometry_columns_auth', 'virts_geometry_columns_field_infos',
                      'virts_geometry_columns_statistics', 'data_licenses', 'KNN']
        # get tablenames from metadata
        tables = sorted([self.encode(x) for x in self.meta.tables.keys()])
        if return_all:
            return tables
        else:
            ret = []
            for i in tables:
                if i not in all_tables and 'idx_' not in i:
                    ret.append(i)
            return ret
    
    def get_unique_directories(self):
        """
        Get a list of directories containing registered scenes

        Returns
        -------
        list[str]
            the directory names
        """
        with self.Session() as session:
            # ORM query, get all directories
            scenes = session.query(self.Data.scene)
        registered = [os.path.dirname(self.encode(x[0])) for x in scenes]
        return list(set(registered))
    
    def import_outdated(self, dbfile):
        """
        import an older database

        Parameters
        ----------
        dbfile: str or Archive
            the old database. If this is a string, the name of a CSV file is expected.

        Returns
        -------

        """
        if isinstance(dbfile, str) and dbfile.endswith('csv'):
            with open(dbfile) as csvfile:
                text = csvfile.read()
                csvfile.seek(0)
                dialect = csv.Sniffer().sniff(text)
                reader = csv.DictReader(csvfile, dialect=dialect)
                scenes = []
                for row in reader:
                    scenes.append(row['scene'])
                self.insert(scenes)
        elif isinstance(dbfile, Archive):
            with self.engine.begin() as conn:
                scenes = conn.exec_driver_sql('SELECT scene from data')
                scenes = [s.scene for s in scenes]
            self.insert(scenes)
            reinsert = dbfile.select_duplicates(value='scene')
            if reinsert is not None:
                self.insert(reinsert)
        else:
            raise RuntimeError("'dbfile' must either be a CSV file name or an Archive object")
    
    def move(self, scenelist, directory, pbar=False):
        """
        Move a list of files while keeping the database entries up to date.
        If a scene is registered in the database (in either the data or duplicates table),
        the scene entry is directly changed to the new location.

        Parameters
        ----------
        scenelist: list[str]
            the file locations
        directory: str
            a folder to which the files are moved
        pbar: bool
            show a progress bar?

        Returns
        -------
        """
        if not os.path.isdir(directory):
            os.mkdir(directory)
        if not os.access(directory, os.W_OK):
            raise RuntimeError('directory cannot be written to')
        failed = []
        double = []
        if pbar:
            progress = pb.ProgressBar(max_value=len(scenelist)).start()
        else:
            progress = None
        
        for i, scene in enumerate(scenelist):
            new = os.path.join(directory, os.path.basename(scene))
            if os.path.isfile(new):
                double.append(new)
                continue
            try:
                shutil.move(scene, directory)
            except shutil.Error:
                failed.append(scene)
                continue
            finally:
                if progress is not None:
                    progress.update(i + 1)
            if self.select(scene=scene) != 0:
                table = 'data'
            else:
                # using core connection to execute SQL syntax (as was before)
                query = '''SELECT scene FROM duplicates WHERE scene='{0}' '''.format(scene)
                with self.engine.begin() as conn:
                    query_duplicates = conn.exec_driver_sql(query)
                if len(query_duplicates) != 0:
                    table = 'duplicates'
                else:
                    table = None
            if table:
                # using core connection to execute SQL syntax (as was before)
                query = '''UPDATE {0} SET scene= '{1}' WHERE scene='{2}' '''.format(table, new, scene)
                with self.engine.begin() as conn:
                    conn.exec_driver_sql(query)
        if progress is not None:
            progress.finish()
        
        if len(failed) > 0:
            log.info('The following scenes could not be moved:\n{}'.format('\n'.join(failed)))
        if len(double) > 0:
            log.info('The following scenes already exist at the target location:\n{}'.format('\n'.join(double)))
    
    def select(self, vectorobject=None, mindate=None, maxdate=None, date_strict=True,
               processdir=None, recursive=False, polarizations=None, return_value="scene", **args):
        """
        select scenes from the database

        Parameters
        ----------
        vectorobject: :class:`~spatialist.vector.Vector` or None
            a geometry with which the scenes need to overlap. The object may only contain one feature.
        mindate: str or datetime.datetime or None
            the minimum acquisition date; strings must be in format YYYYmmddTHHMMSS; default: None
        maxdate: str or datetime.datetime or None
            the maximum acquisition date; strings must be in format YYYYmmddTHHMMSS; default: None
        date_strict: bool
            treat dates as strict limits or also allow flexible limits to incorporate scenes
            whose acquisition period overlaps with the defined limit?

            - strict: start >= mindate & stop <= maxdate
            - not strict: stop >= mindate & start <= maxdate
        processdir: str or None
            A directory to be scanned for already processed scenes;
            the selected scenes will be filtered to those that have not yet been processed. Default: None
        recursive: bool
            (only if `processdir` is not None) should also the subdirectories of the `processdir` be scanned?
        polarizations: list[str] or None
            a list of polarization strings, e.g. ['HH', 'VV']
        return_value: str or List[str]
            the query return value(s). Options:
            
            - `geometry_wkb`: the scene's footprint geometry formatted as WKB
            - `geometry_wkt`: the scene's footprint geometry formatted as WKT
            - `mindate`: the acquisition start datetime in UTC formatted as YYYYmmddTHHMMSS
            - `maxdate`: the acquisition end datetime in UTC formatted as YYYYmmddTHHMMSS
            - all further database column names (see :meth:`~Archive.get_colnames()`)
        **args:
            any further arguments (columns), which are registered in the database. See :meth:`~Archive.get_colnames()`

        Returns
        -------
        List[str] or List[tuple[str]]
            If a single return_value is specified: list of values for that attribute
            If multiple return_values are specified: list of tuples containing the requested attributes
        """
        # Convert return_value to list if it's a string
        if isinstance(return_value, str):
            return_values = [return_value]
        else:
            return_values = return_value
        
        return_values_sql = []
        for val in return_values:
            if val == 'mindate':
                return_values_sql.append('start')
            elif val == 'maxdate':
                return_values_sql.append('stop')
            elif val == 'geometry_wkt':
                prefix = 'ST_' if self.driver == 'postgresql' else ''
                return_values_sql.append(f'{prefix}AsText(geometry) as geometry_wkt')
            elif val == 'geometry_wkb':
                prefix = 'ST_' if self.driver == 'postgresql' else ''
                return_values_sql.append(f'{prefix}AsBinary(geometry) as geometry_wkb')
            else:
                return_values_sql.append(val)
        
        # Validate that all requested return values exist in the database
        valid_columns = self.get_colnames()
        extra = ['mindate', 'maxdate', 'geometry_wkt', 'geometry_wkb']
        normal_returns = [x for x in return_values if x not in extra]
        invalid_returns = [x for x in normal_returns if x not in valid_columns]
        if invalid_returns:
            invalid_str = ', '.join(invalid_returns)
            msg = (f"The following options are not supported as "
                   f"return values: {invalid_str}")
            raise ValueError(msg)
        
        arg_valid = [x for x in args.keys() if x in self.get_colnames()]
        arg_invalid = [x for x in args.keys() if x not in self.get_colnames()]
        if len(arg_invalid) > 0:
            log.info('the following arguments will be ignored as they are not registered in the data base: {}'.format(
                ', '.join(arg_invalid)))
        arg_format = []
        vals = []
        for key in arg_valid:
            if key == 'scene':
                arg_format.append('''scene LIKE '%%{0}%%' '''.format(os.path.basename(args[key])))
            else:
                if isinstance(args[key], (float, int, str)):
                    arg_format.append("""{0}='{1}'""".format(key, args[key]))
                elif isinstance(args[key], (tuple, list)):
                    arg_format.append("""{0} IN ('{1}')""".format(key, "', '".join(map(str, args[key]))))
        
        if mindate:
            if isinstance(mindate, datetime):
                mindate = mindate.strftime('%Y%m%dT%H%M%S')
            if re.search('[0-9]{8}T[0-9]{6}', mindate):
                if date_strict:
                    arg_format.append('start>=?')
                else:
                    arg_format.append('stop>=?')
                vals.append(mindate)
            else:
                log.info('WARNING: argument mindate is ignored, must be in format YYYYmmddTHHMMSS')
        
        if maxdate:
            if isinstance(maxdate, datetime):
                maxdate = maxdate.strftime('%Y%m%dT%H%M%S')
            if re.search('[0-9]{8}T[0-9]{6}', maxdate):
                if date_strict:
                    arg_format.append('stop<=?')
                else:
                    arg_format.append('start<=?')
                vals.append(maxdate)
            else:
                log.info('WARNING: argument maxdate is ignored, must be in format YYYYmmddTHHMMSS')
        
        if polarizations:
            for pol in polarizations:
                if pol in ['HH', 'VV', 'HV', 'VH']:
                    arg_format.append('{}=1'.format(pol.lower()))
        
        if vectorobject:
            if isinstance(vectorobject, Vector):
                if vectorobject.nfeatures > 1:
                    raise RuntimeError("'vectorobject' contains more than one feature.")
                with vectorobject.clone() as vec:
                    vec.reproject(4326)
                    site_geom = vec.convert2wkt(set3D=False)[0]
                # postgres has a different way to store geometries
                if self.driver == 'postgresql':
                    statement = f"st_intersects(geometry, 'SRID=4326; {site_geom}')"
                    arg_format.append(statement)
                else:
                    arg_format.append('st_intersects(GeomFromText(?, 4326), geometry) = 1')
                    vals.append(site_geom)
            else:
                log.info('WARNING: argument vectorobject is ignored, must be of type spatialist.vector.Vector')
        
        if len(arg_format) > 0:
            subquery = ' WHERE {}'.format(' AND '.join(arg_format))
        else:
            subquery = ''
        
        # Modify the query to select the requested return values
        query = 'SELECT {}, outname_base FROM data{}'.format(', '.join(return_values_sql), subquery)
        
        # the query gets assembled stepwise here
        for val in vals:
            query = query.replace('?', """'{0}'""", 1).format(val)
        log.debug(query)
        
        # core SQL execution
        with self.engine.begin() as conn:
            query_rs = conn.exec_driver_sql(query)
            
            if processdir and os.path.isdir(processdir):
                scenes = [x for x in query_rs
                          if len(finder(processdir, [x[-1]],
                                        regex=True, recursive=recursive)) == 0]
            else:
                scenes = query_rs
            
            ret = []
            for x in scenes:
                # If only one return value was requested, append just that value
                if len(return_values) == 1:
                    ret.append(self.encode(x[0]))
                else:
                    # If multiple return values were requested, append a tuple of all values
                    ret.append(tuple(self.encode(val) for val in x[:-1]))  # Exclude outname_base
        
        return ret
    
    def select_duplicates(self, outname_base=None, scene=None, value='id'):
        """
        Select scenes from the duplicates table. In case both `outname_base` and `scene` are set to None all scenes in
        the table are returned, otherwise only those that match the attributes `outname_base` and `scene` if they are not None.

        Parameters
        ----------
        outname_base: str
            the basename of the scene
        scene: str
            the scene name
        value: str
            the return value; either 'id' or 'scene'

        Returns
        -------
        list[str]
            the selected scene(s)
        """
        if value == 'id':
            key = 0
        elif value == 'scene':
            key = 1
        else:
            raise ValueError("argument 'value' must be either 0 or 1")
        
        with self.engine.begin() as conn:
            if not outname_base and not scene:
                # core SQL execution
                scenes = conn.exec_driver_sql('SELECT * from duplicates')
            else:
                cond = []
                arg = []
                if outname_base:
                    cond.append('outname_base=?')
                    arg.append(outname_base)
                if scene:
                    cond.append('scene=?')
                    arg.append(scene)
                query = 'SELECT * from duplicates WHERE {}'.format(' AND '.join(cond))
                for a in arg:
                    query = query.replace('?', ''' '{0}' ''', 1).format(a)
                # core SQL execution
                scenes = conn.exec_driver_sql(query)
            
            ret = []
            for x in scenes:
                ret.append(self.encode(x[key]))
        
        return ret
    
    @property
    def size(self):
        """
        get the number of scenes registered in the database

        Returns
        -------
        tuple[int]
            the number of scenes in (1) the main table and (2) the duplicates table
        """
        # ORM query
        with self.Session() as session:
            r1 = session.query(self.Data.outname_base).count()
            r2 = session.query(self.Duplicates.outname_base).count()
        return r1, r2
    
    def __enter__(self):
        return self
    
    def close(self):
        """
        close the database connection
        """
        self.engine.dispose()
        gc.collect(generation=2)  # this was added as a fix for win PermissionError when deleting sqlite.db files.
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
    
    def drop_element(self, scene, with_duplicates=False):
        """
        Drop a scene from the data table.
        If duplicates table contains matching entry, it will be moved to the data table.

        Parameters
        ----------
        scene: str
            a SAR scene
        with_duplicates: bool
            True: delete matching entry in duplicates table
            False: move matching entry from duplicates into data table

        Returns
        -------
        """
        # save outname_base from to be deleted entry
        search = self.data_schema.select().where(self.data_schema.c.scene == scene)
        entry_data_outname_base = []
        with self.engine.begin() as conn:
            for rowproxy in conn.execute(search):
                entry_data_outname_base.append((rowproxy[12]))
        # log.info(entry_data_outname_base)
        
        # delete entry in data table
        delete_statement = self.data_schema.delete().where(self.data_schema.c.scene == scene)
        with self.engine.begin() as conn:
            conn.execute(delete_statement)
        
        return_sentence = 'Entry with scene-id: \n{} \nwas dropped from data'.format(scene)
        
        # with_duplicates == True, delete entry from duplicates
        if with_duplicates:
            delete_statement_dup = self.duplicates_schema.delete().where(
                self.duplicates_schema.c.outname_base == entry_data_outname_base[0])
            with self.engine.begin() as conn:
                conn.execute(delete_statement_dup)
            
            log.info(return_sentence + ' and duplicates!'.format(scene))
            return
        
        # else select scene info matching outname_base from duplicates
        select_in_duplicates_statement = self.duplicates_schema.select().where(
            self.duplicates_schema.c.outname_base == entry_data_outname_base[0])
        entry_duplicates_scene = []
        with self.engine.begin() as conn:
            for rowproxy in conn.execute(select_in_duplicates_statement):
                entry_duplicates_scene.append((rowproxy[1]))
        
        # check if there is a duplicate
        if len(entry_duplicates_scene) == 1:
            # remove entry from duplicates
            delete_statement_dup = self.duplicates_schema.delete().where(
                self.duplicates_schema.c.outname_base == entry_data_outname_base[0])
            with self.engine.begin() as conn:
                conn.execute(delete_statement_dup)
            
            # insert scene from duplicates into data
            self.insert(entry_duplicates_scene[0])
            
            return_sentence += ' and entry with outname_base \n{} \nand scene \n{} \n' \
                               'was moved from duplicates into data table'.format(
                entry_data_outname_base[0], entry_duplicates_scene[0])
        
        log.info(return_sentence + '!')
    
    def drop_table(self, table):
        """
        Drop a table from the database.

        Parameters
        ----------
        table: str
            the table name

        Returns
        -------
        """
        if table in self.get_tablenames(return_all=True):
            # this removes the idx tables and entries in geometry_columns for sqlite databases
            if self.driver == 'sqlite':
                with self.engine.begin() as conn:
                    query = "SELECT f_table_name FROM geometry_columns"
                    tab_with_geom = [rowproxy[0] for rowproxy
                                     in conn.exec_driver_sql(query)]
                    if table in tab_with_geom:
                        conn.exec_driver_sql("SELECT DropGeoTable('" + table + "')")
            else:
                table_info = Table(table, self.meta, autoload=True, autoload_with=self.engine)
                table_info.drop(self.engine)
            log.info('table {} dropped from database.'.format(table))
        else:
            raise ValueError("table {} is not registered in the database!".format(table))
        self.Base = automap_base(metadata=self.meta)
        self.Base.prepare(self.engine, reflect=True)
    
    @staticmethod
    def __is_open(ip, port):
        """
        Checks server connection, from Ben Curtis (github: Fmstrat)

        Parameters
        ----------
        ip: str
            ip of the server
        port: str or int
            port of the server

        Returns
        -------
        bool:
            is the server reachable?
            
        """
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.settimeout(3)
        try:
            s.connect((ip, int(port)))
            s.shutdown(socket.SHUT_RDWR)
            return True
        except:
            return False
        finally:
            s.close()
    
    def __check_host(self, ip, port):
        """
        Calls __is_open() on ip and port, from Ben Curtis (github: Fmstrat)

        Parameters
        ----------
        ip: str
            ip of the server
        port: str or int
            port of the server

        Returns
        -------
        bool:
            is the server reachable?
        """
        ipup = False
        for i in range(2):
            if self.__is_open(ip, port):
                ipup = True
                break
            else:
                time.sleep(5)
        return ipup


def drop_archive(archive):
    """
    drop (delete) a scene database
    
    Parameters
    ----------
    archive: pyroSAR.drivers.Archive
        the database to be deleted

    Returns
    -------
    
    See Also
    --------
    :func:`sqlalchemy_utils.functions.drop_database()`
    
    Examples
    --------
    >>> pguser = os.environ.get('PGUSER')
    >>> pgpassword = os.environ.get('PGPASSWORD')
    
    >>> db = Archive('test', postgres=True, port=5432, user=pguser, password=pgpassword)
    >>> drop_archive(db)
    """
    if archive.driver == 'postgresql':
        url = archive.url
        archive.close()
        drop_database(url)
    else:
        raise RuntimeError('this function only works for PostgreSQL databases.'
                           'For SQLite databases it is recommended to just delete the DB file.')


def getFileObj(scene, filename):
    """
    Load a file in a SAR scene archive into a readable file object.

    Parameters
    ----------
    scene: str
        the scene archive. Can be either a directory or a compressed archive of type `zip` or `tar.gz`.
    filename: str
        the name of a file in the scene archive, easiest to get with method :meth:`~ID.findfiles`

    Returns
    -------
    ~io.BytesIO
        a file object
    """
    membername = filename.replace(scene, '').strip(r'\/')
    
    if not os.path.exists(scene):
        raise RuntimeError('scene does not exist')
    
    if os.path.isdir(scene):
        obj = BytesIO()
        with open(filename, 'rb') as infile:
            obj.write(infile.read())
        obj.seek(0)
    
    # the scene consists of a single file
    elif os.path.isfile(scene) and scene == filename:
        obj = BytesIO()
        with open(filename, 'rb') as infile:
            obj.write(infile.read())
        obj.seek(0)
    
    elif zf.is_zipfile(scene):
        obj = BytesIO()
        with zf.ZipFile(scene, 'r') as zip:
            obj.write(zip.open(membername).read())
        obj.seek(0)
    
    elif tf.is_tarfile(scene):
        obj = BytesIO()
        tar = tf.open(scene, 'r:gz')
        obj.write(tar.extractfile(membername).read())
        tar.close()
        obj.seek(0)
    else:
        raise RuntimeError('input must be either a file name or a location in an zip or tar archive')
    return obj


def parse_date(x):
    """
    this function gathers known time formats provided in the different SAR products and converts them to a common
    standard of the form YYYYMMDDTHHMMSS

    Parameters
    ----------
    x: str or ~datetime.datetime
        the time stamp to be converted

    Returns
    -------
    str
        the converted time stamp in format YYYYmmddTHHMMSS
    """
    if isinstance(x, datetime):
        return x.strftime('%Y%m%dT%H%M%S')
    elif isinstance(x, str):
        for timeformat in ['%d-%b-%Y %H:%M:%S.%f',
                           '%Y%m%d%H%M%S%f',
                           '%Y-%m-%dT%H:%M:%S.%f',
                           '%Y-%m-%dT%H:%M:%S.%fZ',
                           '%Y%m%d %H:%M:%S.%f']:
            try:
                return strftime('%Y%m%dT%H%M%S', strptime(x, timeformat))
            except (TypeError, ValueError):
                continue
        raise ValueError('unknown time format; check function parse_date')
    else:
        raise ValueError('input must be either a string or a datetime object')
