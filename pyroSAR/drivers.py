#!/usr/bin/env python2.7

##############################################################
# Reading and Organizing system for SAR images
# John Truckenbrodt, Felix Cremer 2016-2017
##############################################################
"""
this script is intended to contain several SAR scene identifier classes to read basic metadata from the scene
folders/files, convert to GAMMA format and do simple pre-processing
"""
from __future__ import print_function, division

import sys


if sys.version_info >= (3, 0):
    from io import StringIO
    from urllib.request import urlopen
else:
    from StringIO import StringIO
    from urllib2 import urlopen

import abc
import ast
import inspect
import math
import os
import re
import ssl
import shutil
import struct
import tarfile as tf
import xml.etree.ElementTree as ET
import zipfile as zf
from datetime import datetime, timedelta
from time import strptime, strftime

import numpy as np
import progressbar as pb
from osgeo import gdal, osr
from osgeo.gdalconst import GA_ReadOnly, GA_Update

####################################################
# Vielleicht in den __init__ rein?
try:
    from pysqlite2 import dbapi2 as sqlite3
except ImportError:
    import sqlite3
####################################################

from pyroSAR import linesimplify as ls
from pyroSAR import spatial
from pyroSAR.ancillary import finder, parse_literal, urlQueryParser, run
from pyroSAR.xml_util import getNamespaces

__LOCAL__ = ['sensor', 'projection', 'orbit', 'polarizations', 'acquisition_mode', 'start', 'stop', 'product',
             'spacing', 'samples', 'lines']


def identify(scene):
    """Return a metadata handler of the given scene."""
    for handler in ID.__subclasses__():
        try:
            return handler(scene)
        except (IOError, KeyError):
            pass
    raise IOError('data format not supported')


def identify_many(scenes):
    """
    return metadata handlers of all valid scenes in a list
    """
    idlist = []
    pbar = pb.ProgressBar(maxval=len(scenes)).start()
    for i, scene in enumerate(scenes):
        if isinstance(scene, ID):
            idlist.append(scene)
        else:
            try:
                id = identify(scene)
                idlist.append(id)
            except IOError:
                continue
        pbar.update(i + 1)
    pbar.finish()
    return idlist


def filter_processed(scenelist, outdir, recursive=False):
    """
    filter a list of pyroSAR objects to those that have not yet been processed and stored in the defined directory
    the search for processed scenes is either done in the directory only or recursively into subdirectories
    the scenes must have been processed with pyroSAR in order to follow the right naming scheme
    """
    return [x for x in scenelist if not x.is_processed(outdir, recursive)]


# todo: add bounding box info to init and summary methods
class ID(object):
    """Abstract class for SAR meta data handlers."""

    def __init__(self, metadict):
        # additional variables? looks, coordinates, ...
        self.locals = __LOCAL__
        for item in self.locals:
            setattr(self, item, metadict[item])

    def bbox(self, outname=None, overwrite=True):
        """
        get the bounding box of a scene either as an vector object or written to a shapefile

        :param outname: the name of the shapefile to be written, default: None
        :param overwrite: ...an existing shapefile
        :return: None if outname is None, otherwise an object of type pyroSAR.spatial.vector
        """
        if outname is None:
            return spatial.bbox(self.getCorners(), self.projection)
        else:
            spatial.bbox(self.getCorners(), self.projection, outname=outname, format='ESRI Shapefile',
                         overwrite=overwrite)

    @property
    def compression(self):
        """
        check whether a scene is compressed into an tarfile or zipfile or not at all

        :return: either 'zip', 'tar' or None
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
        Return the uuid and the metadata that is defined in self.locals as a dictionary
        """
        metadata = {item: self.meta[item] for item in self.locals}
        sq_file = os.path.basename(self.file)
        title = os.path.splitext(sq_file)[0]
        metadata['uuid'] = title
        return metadata

    def export2sqlite(self, dbfile):
        """
        Export relevant metadata to a sqlite database
        :param dbfile: the database file
        :return: None
        """
        with Archive(dbfile) as archive:
            archive.insert(self)

    def examine(self, include_folders=False):
        files = self.findfiles(self.pattern, include_folders=include_folders)
        if len(files) == 1:
            self.file = files[0]
        elif len(files) == 0:
            raise IOError('scene does not match {} naming convention'.format(type(self).__name__))
        else:
            raise IOError('file ambiguity detected:\n{}'.format('\n'.join(files)))

    def findfiles(self, pattern, include_folders=False):
        if os.path.isdir(self.scene):
            files = finder(self.scene, [pattern], regex=True, foldermode=1 if include_folders else 0)
            if re.search(pattern, os.path.basename(self.scene)) and include_folders:
                files.append(self.scene)
        elif zf.is_zipfile(self.scene):
            with zf.ZipFile(self.scene, 'r') as zip:
                files = [os.path.join(self.scene, x) for x in zip.namelist() if
                         re.search(pattern, os.path.basename(x.strip('/')))]
                if include_folders:
                    files = [x.strip('/') for x in files]
                else:
                    files = [x for x in files if not x.endswith('/')]
        elif tf.is_tarfile(self.scene):
            tar = tf.open(self.scene)
            files = [x for x in tar.getnames() if re.search(pattern, os.path.basename(x.strip('/')))]
            if not include_folders:
                files = [x for x in files if not tar.getmember(x).isdir()]
            tar.close()
            files = [os.path.join(self.scene, x) for x in files]
        else:
            files = [self.scene] if re.search(pattern, self.scene) else []
        return files

    def gdalinfo(self, scene):
        """
        Args:
            scene: an archive containing a SAR scene

        returns a dictionary of metadata attributes
        """
        self.scene = os.path.realpath(scene)
        files = self.findfiles('(?:\.[NE][12]$|DAT_01\.001$|product\.xml|manifest\.safe$)')

        if len(files) == 1:
            prefix = {'zip': '/vsizip/', 'tar': '/vsitar/', None: ''}[self.compression]
            header = files[0]
        elif len(files) > 1:
            raise IOError('file ambiguity detected')
        else:
            raise IOError('file type not supported')

        meta = {}

        ext_lookup = {'.N1': 'ASAR', '.E1': 'ERS1', '.E2': 'ERS2'}
        extension = os.path.splitext(header)[1]
        if extension in ext_lookup:
            meta['sensor'] = ext_lookup[extension]

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

            if re.search('(?:LAT|LONG)', entry[0]):
                entry[1] /= 1000000.
            meta[entry[0]] = entry[1]
        return meta

    @abc.abstractmethod
    def getCorners(self):
        raise NotImplementedError

    def getFileObj(self, filename):
        """
        load a file into a readable file object
        if the scene is unpacked this will be a regular 'file' object
        for a tarfile this is an object of type 'tarfile.ExtFile'
        for a zipfile this is an StringIO object (the zipfile.ExtFile object does not support setting file pointers via function 'seek', which is needed later on)
        """
        membername = filename.replace(self.scene, '').strip('/')

        if os.path.isdir(self.scene):
            obj = open(filename)
        elif zf.is_zipfile(self.scene):
            obj = StringIO()
            with zf.ZipFile(self.scene, 'r') as zip:
                obj.write(zip.open(membername).read())
            obj.seek(0)

        elif tf.is_tarfile(self.scene):
            obj = StringIO()
            tar = tf.open(self.scene, 'r:gz')
            obj.write(tar.extractfile(membername).read())
            tar.close()
        else:
            raise IOError('input must be either a file name or a location in an zip or tar archive')
        return obj

    def getGammaImages(self, directory=None):
        if directory is None:
            if hasattr(self, 'gammadir'):
                directory = self.gammadir
            else:
                raise IOError(
                    'directory missing; please provide directory to function or define object attribute "gammadir"')
        return [x for x in finder(directory, [self.outname_base()], regex=True) if
                not re.search('\.(?:par|hdr|aux\.xml)$', x)]

    def getHGT(self):
        """
        Returns: names of all SRTM hgt tiles overlapping with the SAR scene
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
        check whether a scene has already been processed and stored in the defined output directory (and subdirectories if recursive)
        """
        if os.path.isdir(outdir):
            # '{}.*tif$'.format(self.outname_base())
            return len(finder(outdir, [self.outname_base()], regex=True, recursive=recursive)) != 0
        else:
            return False

    def outname_base(self):
        fields = ('{:_<4}'.format(self.sensor),
                  '{:_<4}'.format(self.acquisition_mode),
                  self.orbit,
                  self.start)
        return '_'.join(fields)

    @staticmethod
    def parse_date(x):
        """
        this function gathers known time formats provided in the different SAR products and converts them to a common standard of the form YYYYMMDDTHHMMSS
        """
        # todo: check module time for more general approaches
        for timeformat in ['%d-%b-%Y %H:%M:%S.%f',
                           '%Y%m%d%H%M%S%f',
                           '%Y-%m-%dT%H:%M:%S.%f',
                           '%Y-%m-%dT%H:%M:%S.%fZ',
                           '%Y%m%d %H:%M:%S.%f']:
            try:
                return strftime('%Y%m%dT%H%M%S', strptime(x, timeformat))
            except (TypeError, ValueError):
                continue
        raise ValueError('unknown time format; check function ID.parse_date')

    def summary(self):
        for item in sorted(self.locals):
            print('{0}: {1}'.format(item, getattr(self, item)))

    @abc.abstractmethod
    def scanMetadata(self):
        raise NotImplementedError

    @abc.abstractmethod
    def unpack(self, directory):
        raise NotImplementedError

    # todo: replace with functionality from module archivist
    def _unpack(self, directory, offset=None, overwrite=False):
        """
        general function for unpacking scene archives
        :param directory: the name of the directory in which the files are written
        :param offset: an archive directory offset; to be define if only a subdirectory is to be unpacked (see e.g. TSX:unpack)
        :param overwrite: should an existing directory be overwritten?
        :return: None
        """
        if os.path.isdir(directory):
            if overwrite:
                shutil.rmtree(directory)
            else:
                raise RuntimeError('target scene directory already exists: {}'.format(directory))
        os.makedirs(directory)

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
                        outname = os.path.join(directory, item.replace(header, '', 1))
                        if item.endswith('/'):
                            os.makedirs(outname)
                        else:
                            try:
                                with open(outname, 'w') as outfile:
                                    outfile.write(archive.read(item))
                            except zf.BadZipfile:
                                print('corrupt archive, unpacking failed')
                                continue
                archive.close()
            else:
                archive.extractall(directory)
                archive.close()

        else:
            print('unpacking is only supported for TAR and ZIP archives')
            return

        self.scene = directory
        main = os.path.join(self.scene, os.path.basename(self.file))
        self.file = main if os.path.isfile(main) else self.scene


class CEOS_ERS(ID):
    """
    Handler class for ERS data in CEOS format
    
    References:
        ER-IS-EPO-GS-5902-3: Annex C. ERS SAR.SLC/SLC-I. CCT and EXABYTE (ESA 1998)
    """

    def __init__(self, scene):
        self.pattern = r'(?P<product_id>(?:SAR|ASA)_(?:IM(?:S|P|G|M|_)|AP(?:S|P|G|M|_)|WV(?:I|S|W|_)|WS(?:M|S|_))_[012B][CP])' \
                       r'(?P<processing_stage_flag>[A-Z])' \
                       r'(?P<originator_ID>[A-Z\-]{3})' \
                       r'(?P<start_day>[0-9]{8})_' \
                       r'(?P<start_time>[0-9]{6})_' \
                       r'(?P<duration>[0-9]{8})' \
                       r'(?P<phase>[0-9A-Z]{1})' \
                       r'(?P<cycle>[0-9]{3})_' \
                       r'(?P<relative_orbit>[0-9]{5})_' \
                       r'(?P<absolute_orbit>[0-9]{5})_' \
                       r'(?P<counter>[0-9]{4,})\.' \
                       r'(?P<satellite_ID>[EN][12])' \
                       r'(?P<extension>(?:\.zip|\.tar\.gz|\.PS|))$'

        self.pattern_pid = r'(?P<sat_id>(?:SAR|ASA))_' \
                           r'(?P<image_mode>(?:IM(?:S|P|G|M|_)|AP(?:S|P|G|M|_)|WV(?:I|S|W|_)|WS(?:M|S|_)))_' \
                           r'(?P<processing_level>[012B][CP])'

        self.scene = os.path.realpath(scene)

        self.examine()

        match = re.match(re.compile(self.pattern), os.path.basename(self.file))
        match2 = re.match(re.compile(self.pattern_pid), match.group('product_id'))

        if re.search('IM__0', match.group('product_id')):
            raise IOError('product level 0 not supported (yet)')

        self.meta = self.gdalinfo(self.scene)

        self.meta['acquisition_mode'] = match2.group('image_mode')
        self.meta['polarizations'] = ['VV']
        self.meta['product'] = 'SLC' if self.meta['acquisition_mode'] in ['IMS', 'APS', 'WSS'] else 'PRI'
        self.meta['spacing'] = (self.meta['CEOS_PIXEL_SPACING_METERS'], self.meta['CEOS_LINE_SPACING_METERS'])
        self.meta['sensor'] = self.meta['CEOS_MISSION_ID']
        self.meta['incidence_angle'] = self.meta['CEOS_INC_ANGLE']
        self.meta['k_db'] = -10 * math.log(float(self.meta['CEOS_CALIBRATION_CONSTANT_K']), 10)
        self.meta['sc_db'] = {'ERS1': 59.61, 'ERS2': 60}[self.meta['sensor']]

        # acquire additional metadata from the file LEA_01.001
        self.meta.update(self.scanMetadata())

        # register the standardized meta attributes as object attributes
        ID.__init__(self, self.meta)

    # todo: change coordinate extraction to the exact boundaries of the image (not outer pixel center points)
    def getCorners(self):
        lat = [x[1][1] for x in self.meta['gcps']]
        lon = [x[1][0] for x in self.meta['gcps']]
        return {'xmin': min(lon), 'xmax': max(lon), 'ymin': min(lat), 'ymax': max(lat)}

    def unpack(self, directory, overwrite=False):
        if self.sensor in ['ERS1', 'ERS2']:
            base_file = re.sub('\.PS$', '', os.path.basename(self.file))
            base_dir = os.path.basename(directory.strip('/'))

            outdir = directory if base_file == base_dir else os.path.join(directory, base_file)

            self._unpack(outdir, overwrite=overwrite)
        else:
            raise NotImplementedError('sensor {} not implemented yet'.format(self.sensor))

    def scanMetadata(self):
        """
        read the leader file and extract relevant metadata
        """
        lea_obj = self.getFileObj(self.findfiles('LEA_01.001')[0])
        lea = lea_obj.read()
        lea_obj.close()
        meta = dict()
        offset = 720
        looks_range = float(lea[(offset + 1174):(offset + 1190)])
        looks_azimuth = float(lea[(offset + 1190):(offset + 1206)])
        meta['looks'] = (looks_range, looks_azimuth)
        meta['heading'] = float(lea[(offset + 468):(offset + 476)])
        meta['orbit'] = 'D' if meta['heading'] > 180 else 'A'
        orbitNumber, frameNumber = map(int, re.findall('[0-9]+', lea[(offset + 36):(offset + 68)]))
        meta['orbitNumber'] = orbitNumber
        meta['frameNumber'] = frameNumber
        meta['start'] = self.parse_date(lea[(offset + 1814):(offset + 1838)])
        meta['stop'] = self.parse_date(lea[(offset + 1862):(offset + 1886)])
        # the following parameters are already read by gdalinfo
        # meta['sensor'] = lea[(offset+396):(offset+412)].strip()
        # spacing_azimuth = float(lea[(offset+1686):(offset+1702)])
        # spacing_range = float(lea[(offset+1702):(offset+1718)])
        # meta['spacing'] = (spacing_range, spacing_azimuth)
        # meta['incidence_angle'] = float(lea[(offset+484):(offset+492)])
        meta['proc_facility'] = lea[(offset + 1045):(offset + 1061)].strip()
        meta['proc_system'] = lea[(offset + 1061):(offset + 1069)].strip()
        meta['proc_version'] = lea[(offset + 1069):(offset + 1077)].strip()
        # text_subset = lea[re.search('FACILITY RELATED DATA RECORD \[ESA GENERAL TYPE\]', lea).start() - 13:]
        # meta['k_db'] = -10*math.log(float(text_subset[663:679].strip()), 10)
        # meta['antenna_flag'] = int(text_subset[659:663].strip())
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

    PALSAR-1:

        References:
            NEB-070062B: ALOS/PALSAR Level 1.1/1.5 product Format description (JAXA 2009)
    
        Products / processing levels:
            1.0
            1.1
            1.5
    
        Acquisition modes:
            AB: [SP][HWDPC]
            A: supplemental remarks of the sensor type:
                S: Wide observation mode
                P: all other modes
            B: observation mode
                H: Fine mode
                W: ScanSAR mode
                D: Direct downlink mode
                P: Polarimetry mode
                C: Calibration mode
    
    PALSAR-2:
    
        References:
            ALOS-2/PALSAR-2 Level 1.1/1.5/2.1/3.1 CEOS SAR Product Format Description
        
        Products / processing levels:
            1.0
            1.1
            1.5
        
        Acquisition modes:
            SBS: Spotlight mode 
            UBS: Ultra-fine mode Single polarization 
            UBD: Ultra-fine mode Dual polarization 
            HBS: High-sensitive mode Single polarization
            HBD: High-sensitive mode Dual polarization 
            HBQ: High-sensitive mode Full (Quad.) polarimetry 
            FBS: Fine mode Single polarization 
            FBD: Fine mode Dual polarization 
            FBQ: Fine mode Full (Quad.) polarimetry 
            WBS: Scan SAR nominal [14MHz] mode Single polarization 
            WBD: Scan SAR nominal [14MHz] mode Dual polarization 
            WWS: Scan SAR nominal [28MHz] mode Single polarization 
            WWD: Scan SAR nominal [28MHz] mode Dual polarization 
            VBS: Scan SAR wide mode Single polarization 
            VBD: Scan SAR wide mode Dual polarization
    """

    def __init__(self, scene):

        self.scene = os.path.realpath(scene)

        patterns = [r'^LED-ALPSR'
                    r'(?P<sub>P|S)'
                    r'(?P<orbit>[0-9]{5})'
                    r'(?P<frame>[0-9]{4})-'
                    r'(?P<mode>[HWDPC])'
                    r'(?P<level>1\.[015])'
                    r'(?P<proc>G|_)'
                    r'(?P<proj>[UPML_])'
                    r'(?P<orbit_dir>A|D)$',
                    r'^LED-ALOS2'
                    r'(?P<orbit>[0-9]{5})'
                    r'(?P<frame>[0-9]{4})-'
                    r'(?P<date>[0-9]{6})-'
                    r'(?P<mode>SBS|UBS|UBD|HBS|HBD|HBQ|FBS|FBD|FBQ|WBS|WBD|WWS|WWD|VBS|VBD)'
                    r'(?P<look_dir>L|R)'
                    r'(?P<level>1\.0|1\.1|1\.5|2\.1|3\.1)'
                    r'(?P<proc>[GR_])'
                    r'(?P<proj>[UPML_])'
                    r'(?P<orbit_dir>A|D)$']

        for i, pattern in enumerate(patterns):
            self.pattern = pattern
            try:
                self.examine()
                break
            except IOError as e:
                if i + 1 == len(patterns):
                    raise e
                else:
                    continue

        self.meta = self.scanMetadata()

        # register the standardized meta attributes as object attributes
        ID.__init__(self, self.meta)

    def getLeaderfile(self):
        led_filename = self.findfiles(self.pattern)[0]
        led_obj = self.getFileObj(led_filename)
        led = led_obj.read()
        led_obj.close()
        return led

    def parseSummary(self):
        try:
            summary_file = self.getFileObj(self.findfiles('summary|workreport')[0])
        except IndexError:
            return {}
        text = summary_file.read().strip()
        summary_file.close()
        summary = ast.literal_eval('{"' + re.sub('\s*=', '":', text).replace('\n', ',"') + '}')
        for x, y in summary.iteritems():
            summary[x] = parse_literal(y)
        return summary

    def scanMetadata(self):
        led_filename = self.findfiles(self.pattern)[0]
        led_obj = self.getFileObj(led_filename)
        led = led_obj.read()
        led_obj.close()

        meta = self.parseSummary()

        p0 = 0
        p1 = struct.unpack('>i', led[8:12])[0]
        fileDescriptor = led[p0:p1]
        dss_n = int(fileDescriptor[180:186])
        dss_l = int(fileDescriptor[186:192])
        mpd_n = int(fileDescriptor[192:198])
        mpd_l = int(fileDescriptor[198:204])
        ppd_n = int(fileDescriptor[204:210])
        ppd_l = int(fileDescriptor[210:216])
        adr_n = int(fileDescriptor[216:222])
        adr_l = int(fileDescriptor[222:228])
        rdr_n = int(fileDescriptor[228:234])
        rdr_l = int(fileDescriptor[234:240])
        dqs_n = int(fileDescriptor[252:258])
        dqs_l = int(fileDescriptor[258:264])
        meta['sensor'] = {'AL1': 'PSR1', 'AL2': 'PSR2'}[fileDescriptor[48:51]]

        p0 = p1
        p1 += dss_l * dss_n
        dataSetSummary = led[p0:p1]

        if mpd_n > 0:
            p0 = p1
            p1 += mpd_l * mpd_n
            mapProjectionData = led[p0:p1]

            lat = map(float, [mapProjectionData[1072:1088],
                              mapProjectionData[1104:1120],
                              mapProjectionData[1136:1152],
                              mapProjectionData[1168:1184]])
            lon = map(float, [mapProjectionData[1088:1104],
                              mapProjectionData[1120:1136],
                              mapProjectionData[1152:1168],
                              mapProjectionData[1184:1200]])
            meta['corners'] = {'xmin': min(lon), 'xmax': max(lon), 'ymin': min(lat), 'ymax': max(lat)}

            # https://github.com/datalyze-solutions/LandsatProcessingPlugin/blob/master/src/metageta/formats/alos.py

            src_srs = osr.SpatialReference()
            # src_srs.SetGeogCS('GRS 1980','GRS 1980','GRS 1980',6378137.00000,298.2572220972)
            src_srs.SetWellKnownGeogCS("WGS84")
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
            meta['projection'] = 'GEOGCS["WGS 84",' \
                                 'DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],' \
                                 'PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],' \
                                 'UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],' \
                                 'AUTHORITY["EPSG","4326"]]'

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

        # for i in range(0, 10):
        #     p0 = p1
        #     length = struct.unpack('>i', led[(p0 + 8):(p0 + 12)])[0]
        #     print length
        #     p1 += length
        #     facilityRelatedData[i] = led[p0:p1]
        #
        # facilityRelatedData[10] = led[p1:]

        meta['lines'] = int(dataSetSummary[324:332]) * 2
        meta['samples'] = int(dataSetSummary[332:340]) * 2
        meta['incidence'] = float(dataSetSummary[484:492])
        meta['wavelength'] = float(dataSetSummary[500:516]) * 100  # in cm
        meta['proc_facility'] = dataSetSummary[1046:1062].strip()
        meta['proc_system'] = dataSetSummary[1062:1070].strip()
        meta['proc_version'] = dataSetSummary[1070:1078].strip()

        azlks = float(dataSetSummary[1174:1190])
        rlks = float(dataSetSummary[1190:1206])
        meta['looks'] = (rlks, azlks)

        meta['orbit'] = dataSetSummary[1534:1542].strip()[0]

        spacing_azimuth = float(dataSetSummary[1686:1702])
        spacing_range = float(dataSetSummary[1702:1718])
        meta['spacing'] = (spacing_range, spacing_azimuth)

        match = re.match(re.compile(self.pattern), os.path.basename(led_filename))

        if meta['sensor'] == 'PSR1':
            meta['acquisition_mode'] = match.group('sub') + match.group('mode')
        else:
            meta['acquisition_mode'] = match.group('mode')
        meta['product'] = match.group('level')

        try:
            meta['start'] = self.parse_date(self.meta['Img_SceneStartDateTime'])
            meta['stop'] = self.parse_date(self.meta['Img_SceneEndDateTime'])
        except (AttributeError, KeyError):
            try:
                start_string = re.search('Img_SceneStartDateTime[ ="0-9:.]*', led).group()
                stop_string = re.search('Img_SceneEndDateTime[ ="0-9:.]*', led).group()
                meta['start'] = self.parse_date(re.search('\d+\s[\d:.]+', start_string).group())
                meta['stop'] = self.parse_date(re.search('\d+\s[\d:.]+', stop_string).group())
            except AttributeError:
                raise IndexError('start and stop time stamps cannot be extracted; see file {}'.format(led_filename))

        meta['polarizations'] = [re.search('[HV]{2}', os.path.basename(x)).group(0) for x in self.findfiles('^IMG-')]
        meta['k_dB'] = float(radiometricData[20:36])
        return meta

    def unpack(self, directory, overwrite=False):
        outdir = os.path.join(directory, os.path.basename(self.file).replace('LED-', ''))
        self._unpack(outdir, overwrite=overwrite)

    # todo: create summary/workreport file entries for coordinates if they were read from an IMG file
    def getCorners(self):
        if 'corners' not in self.meta.keys():
            lat = [y for x, y in self.meta.iteritems() if 'Latitude' in x]
            lon = [y for x, y in self.meta.iteritems() if 'Longitude' in x]
            if len(lat) == 0 or len(lon) == 0:
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

            self.meta['corners'] = {'xmin': min(lon), 'xmax': max(lon), 'ymin': min(lat), 'ymax': max(lat)}

        return self.meta['corners']


class ESA(ID):
    """
    Handler class for SAR data in ESA format (Envisat ASAR, ERS-1/2)
    """

    def __init__(self, scene):

        self.pattern = r'(?P<product_id>(?:SAR|ASA)_(?:IM(?:S|P|G|M|_)|AP(?:S|P|G|M|_)|WV(?:I|S|W|_)|WS(?:M|S|_))_[012B][CP])' \
                       r'(?P<processing_stage_flag>[A-Z])' \
                       r'(?P<originator_ID>[A-Z\-]{3})' \
                       r'(?P<start_day>[0-9]{8})_' \
                       r'(?P<start_time>[0-9]{6})_' \
                       r'(?P<duration>[0-9]{8})' \
                       r'(?P<phase>[0-9A-Z]{1})' \
                       r'(?P<cycle>[0-9]{3})_' \
                       r'(?P<relative_orbit>[0-9]{5})_' \
                       r'(?P<absolute_orbit>[0-9]{5})_' \
                       r'(?P<counter>[0-9]{4,})\.' \
                       r'(?P<satellite_ID>[EN][12])' \
                       r'(?P<extension>(?:\.zip|\.tar\.gz|))$'

        self.pattern_pid = r'(?P<sat_id>(?:SAR|ASA))_' \
                           r'(?P<image_mode>(?:IM(?:S|P|G|M|_)|AP(?:S|P|G|M|_)|WV(?:I|S|W|_)|WS(?:M|S|_)))_' \
                           r'(?P<processing_level>[012B][CP])'

        self.scene = os.path.realpath(scene)

        self.examine()

        match = re.match(re.compile(self.pattern), os.path.basename(self.file))
        match2 = re.match(re.compile(self.pattern_pid), match.group('product_id'))

        if re.search('IM__0', match.group('product_id')):
            raise IOError('product level 0 not supported (yet)')

        self.meta = self.gdalinfo(self.scene)

        self.meta['acquisition_mode'] = match2.group('image_mode')

        self.meta['product'] = 'SLC' if self.meta['acquisition_mode'] in ['IMS', 'APS', 'WSS'] else 'PRI'

        if self.meta['sensor'] == 'ASAR':
            self.meta['polarizations'] = [y.replace('/', '') for x, y in self.meta.iteritems() if
                                          'TX_RX_POLAR' in x and len(y) == 3]
        elif self.meta['sensor'] in ['ERS1', 'ERS2']:
            self.meta['polarizations'] = ['VV']

        self.meta['orbit'] = self.meta['SPH_PASS'][0]
        self.meta['start'] = self.meta['MPH_SENSING_START']
        self.meta['stop'] = self.meta['MPH_SENSING_STOP']
        self.meta['spacing'] = (self.meta['SPH_RANGE_SPACING'], self.meta['SPH_AZIMUTH_SPACING'])
        self.meta['looks'] = (self.meta['SPH_RANGE_LOOKS'], self.meta['SPH_AZIMUTH_LOOKS'])

        # register the standardized meta attributes as object attributes
        ID.__init__(self, self.meta)

    def getCorners(self):
        lon = [self.meta[x] for x in self.meta if re.search('LONG', x)]
        lat = [self.meta[x] for x in self.meta if re.search('LAT', x)]
        return {'xmin': min(lon), 'xmax': max(lon), 'ymin': min(lat), 'ymax': max(lat)}

    def unpack(self, directory, overwrite=False):
        base_file = os.path.basename(self.file).strip('\.zip|\.tar(?:\.gz|)')
        base_dir = os.path.basename(directory.strip('/'))

        outdir = directory if base_file == base_dir else os.path.join(directory, base_file)

        self._unpack(outdir, overwrite=overwrite)


class SAFE(ID):
    """
    Handler class for Sentinel-1 data

    References:
        S1-RS-MDA-52-7443 Sentinel-1 IPF Auxiliary Product Specification
        MPC-0243 Masking "No-value" Pixels on GRD Products generated by the Sentinel-1 ESA IPF
    """

    def __init__(self, scene):

        self.scene = os.path.realpath(scene)

        self.pattern = r'^(?P<sensor>S1[AB])_' \
                       r'(?P<beam>S1|S2|S3|S4|S5|S6|IW|EW|WV|EN|N1|N2|N3|N4|N5|N6|IM)_' \
                       r'(?P<product>SLC|GRD|OCN)(?:F|H|M|_)_' \
                       r'(?:1|2)' \
                       r'(?P<category>S|A)' \
                       r'(?P<pols>SH|SV|DH|DV|VV|HH|HV|VH)_' \
                       r'(?P<start>[0-9]{8}T[0-9]{6})_' \
                       r'(?P<stop>[0-9]{8}T[0-9]{6})_' \
                       r'(?P<orbitNumber>[0-9]{6})_' \
                       r'(?P<dataTakeID>[0-9A-F]{6})_' \
                       r'(?P<productIdentifier>[0-9A-F]{4})' \
                       r'\.SAFE$'

        self.pattern_ds = r'^s1[ab]-' \
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
            raise IOError('folder does not match S1 scene naming convention')

        # scan the manifest.safe file and add selected attributes to a meta dictionary
        self.meta = self.scanMetadata()
        self.meta['projection'] = 'GEOGCS["WGS 84",' \
                                  'DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],' \
                                  'PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],' \
                                  'UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],' \
                                  'AUTHORITY["EPSG","4326"]]'

        annotations = self.findfiles(self.pattern_ds)
        ann_xml = self.getFileObj(annotations[0])
        ann_tree = ET.fromstring(ann_xml.read())
        ann_xml.close()
        self.meta['spacing'] = tuple(
            [float(ann_tree.find('.//{}PixelSpacing'.format(dim)).text) for dim in ['range', 'azimuth']])
        self.meta['samples'] = int(ann_tree.find('.//imageAnnotation/imageInformation/numberOfSamples').text)
        self.meta['lines'] = int(ann_tree.find('.//imageAnnotation/imageInformation/numberOfLines').text)

        # register the standardized meta attributes as object attributes
        ID.__init__(self, self.meta)

        self.gammafiles = {'slc': [], 'pri': [], 'grd': []}

    def removeGRDBorderNoise(self):
        """
        mask out Sentinel-1 image border noise
        reference:
            'Masking "No-value" Pixels on GRD Products generated by the Sentinel-1 ESA IPF' (issue 1, June 2015)
            available online under 'https://sentinel.esa.int/web/sentinel/user-guides/sentinel-1-sar/document-library'
        """
        if self.compression is not None:
            raise RuntimeError('scene is not yet unpacked')

        blocksize = 2000

        # compute noise scaling factor
        if self.meta['IPF_version'] < 2.5:
            knoise = {'IW': 75088.7, 'EW': 56065.87}[self.acquisition_mode]
            cads = self.getFileObj(self.findfiles('calibration-s1[ab]-[ie]w-grd-(?:hh|vv)')[0])
            caltree = ET.fromstring(cads.read())
            cads.close()
            adn = float(caltree.find('.//calibrationVector/dn').text.split()[0])
            if self.meta['IPF_version'] < 2.34:
                scalingFactor = knoise * adn
            else:
                scalingFactor = knoise * adn * adn
        else:
            scalingFactor = 1

        # read noise vectors from corresponding annotation xml
        noisefile = self.getFileObj(self.findfiles('noise-s1[ab]-[ie]w-grd-(?:hh|vv)')[0])
        noisetree = ET.fromstring(noisefile.read())
        noisefile.close()
        noiseVectors = noisetree.findall('.//noiseVector')

        # define boundaries of image subsets to be masked (4x the first lines/samples of the image boundaries)
        subsets = [(0, 0, blocksize, self.lines),
                   (0, 0, self.samples, blocksize),
                   (self.samples - blocksize, 0, self.samples, self.lines),
                   (0, self.lines - blocksize, self.samples, self.lines)]

        # extract column indices of noise vectors
        yi = np.array([int(x.find('line').text) for x in noiseVectors])

        # create links to the tif files for a master co-polarization and all other polarizations as slaves
        master = self.findfiles('s1.*(?:vv|hh).*tiff')[0]
        ras_master = gdal.Open(master, GA_Update)
        ras_slaves = [gdal.Open(x, GA_Update) for x in self.findfiles('s1.*tiff') if x != master]

        outband_master = ras_master.GetRasterBand(1)
        outband_slaves = [x.GetRasterBand(1) for x in ras_slaves]

        # iterate over the four image subsets
        for subset in subsets:
            print(subset)
            xmin, ymin, xmax, ymax = subset
            xdiff = xmax - xmin
            ydiff = ymax - ymin
            # linear interpolation of noise vectors to array
            noise_interp = np.empty((ydiff, xdiff), dtype=float)
            for i in range(0, len(noiseVectors)):
                if ymin <= yi[i] <= ymax:
                    # extract row indices of noise vector
                    xi = map(int, noiseVectors[i].find('pixel').text.split())
                    # extract noise values
                    noise = map(float, noiseVectors[i].find('noiseLut').text.split())
                    # interpolate values along rows
                    noise_interp[yi[i] - ymin, :] = np.interp(range(0, xdiff), xi, noise)
            for i in range(0, xdiff):
                yi_t = yi[(ymin <= yi) & (yi <= ymax)] - ymin
                # interpolate values along columns
                noise_interp[:, i] = np.interp(range(0, ydiff), yi_t, noise_interp[:, i][yi_t])

            # read subset of image to array and subtract interpolated noise (denoising)
            mat_master = outband_master.ReadAsArray(*[xmin, ymin, xdiff, ydiff])
            denoisedBlock = mat_master.astype(float) ** 2 - noise_interp * scalingFactor
            # mask out all pixels with a value below 0.5 in the denoised block or 30 in the original block
            denoisedBlock[(denoisedBlock < 0.5) | (mat_master < 30)] = 0
            denoisedBlock = np.sqrt(denoisedBlock)

            # mask out negative values
            def helper1(x):
                return len(x) - np.argmax(x > 0)

            def helper2(x):
                return len(x) - np.argmax(x[::-1] > 0)

            if subset == (0, 0, blocksize, self.lines):
                border = np.apply_along_axis(helper1, 1, denoisedBlock)
                border = blocksize - np.array(ls.reduce(border))
                for j in range(0, ydiff):
                    denoisedBlock[j, :border[j]] = 0
                    denoisedBlock[j, border[j]:] = 1
            elif subset == (0, self.lines - blocksize, self.samples, self.lines):
                border = np.apply_along_axis(helper2, 0, denoisedBlock)
                border = ls.reduce(border)
                for j in range(0, xdiff):
                    denoisedBlock[border[j]:, j] = 0
                    denoisedBlock[:border[j], j] = 1
            elif subset == (self.samples - blocksize, 0, self.samples, self.lines):
                border = np.apply_along_axis(helper2, 1, denoisedBlock)
                border = ls.reduce(border)
                for j in range(0, ydiff):
                    denoisedBlock[j, border[j]:] = 0
                    denoisedBlock[j, :border[j]] = 1
            elif subset == (0, 0, self.samples, blocksize):
                border = np.apply_along_axis(helper1, 0, denoisedBlock)
                border = blocksize - np.array(ls.reduce(border))
                for j in range(0, xdiff):
                    denoisedBlock[:border[j], j] = 0
                    denoisedBlock[border[j]:, j] = 1

            mat_master[denoisedBlock == 0] = 0
            # write modified array back to original file
            outband_master.WriteArray(mat_master, xmin, ymin)
            outband_master.FlushCache()
            # perform reading, masking and writing for all other polarizations
            for outband in outband_slaves:
                mat = outband.ReadAsArray(*[xmin, ymin, xdiff, ydiff])
                mat[denoisedBlock == 0] = 0
                outband.WriteArray(mat, xmin, ymin)
                outband.FlushCache()
        # detach file links
        outband_master = None
        ras_master = None
        for outband in outband_slaves:
            outband = None
        for ras in ras_slaves:
            ras = None

    def getCorners(self):
        coordinates = self.meta['coordinates']
        lat = [x[0] for x in coordinates]
        lon = [x[1] for x in coordinates]
        return {'xmin': min(lon), 'xmax': max(lon), 'ymin': min(lat), 'ymax': max(lat)}

    def getOSV(self, outdir):
        date = datetime.strptime(self.start, '%Y%m%dT%H%M%S')

        before = (date - timedelta(days=1)).strftime('%Y-%m-%d')
        after = (date + timedelta(days=1)).strftime('%Y-%m-%d')

        query = dict()
        query['mission'] = self.sensor
        query['validity_start_time'] = '{0}..{1}'.format(before, after)

        remote_poe = 'https://qc.sentinel1.eo.esa.int/aux_poeorb/'

        pattern = 'S1[AB]_OPER_AUX_(?:POE|RES)ORB_OPOD_[0-9TV_]{48}\.EOF'

        sslcontext = ssl._create_unverified_context()

        subaddress = urlQueryParser(remote_poe, query)
        response = urlopen(subaddress, context=sslcontext).read()
        remotes = [os.path.join(remote_poe, x) for x in sorted(set(re.findall(pattern, response)))]

        if not os.access(outdir, os.W_OK):
            raise RuntimeError('insufficient directory permissions, unable to write')
        downloads = [x for x in remotes if not os.path.isfile(os.path.join(outdir, os.path.basename(x)))]
        for item in downloads:
            infile = urlopen(item, context=sslcontext)
            with open(os.path.join(outdir, os.path.basename(item)), 'wb') as outfile:
                outfile.write(infile.read())
            infile.close()

    def scanMetadata(self):
        """
        read the manifest.safe file and extract relevant metadata
        """
        manifest = self.getFileObj(self.findfiles('manifest.safe')[0])
        namespaces = getNamespaces(manifest)
        tree = ET.fromstring(manifest.read())
        manifest.close()

        meta = dict()
        meta['acquisition_mode'] = tree.find('.//s1sarl1:mode', namespaces).text
        meta['acquisition_time'] = dict(
            [(x, tree.find('.//safe:{}Time'.format(x), namespaces).text) for x in ['start', 'stop']])
        meta['start'], meta['stop'] = (self.parse_date(meta['acquisition_time'][x]) for x in ['start', 'stop'])
        meta['coordinates'] = [tuple([float(y) for y in x.split(',')]) for x in
                               tree.find('.//gml:coordinates', namespaces).text.split()]
        meta['orbit'] = tree.find('.//s1:pass', namespaces).text[0]
        meta['orbitNumbers_abs'] = dict(
            [(x, int(tree.find('.//safe:orbitNumber[@type="{0}"]'.format(x), namespaces).text)) for x in
             ['start', 'stop']])
        meta['orbitNumbers_rel'] = dict(
            [(x, int(tree.find('.//safe:relativeOrbitNumber[@type="{0}"]'.format(x), namespaces).text)) for x in
             ['start', 'stop']])
        meta['polarizations'] = [x.text for x in tree.findall('.//s1sarl1:transmitterReceiverPolarisation', namespaces)]
        meta['product'] = tree.find('.//s1sarl1:productType', namespaces).text
        meta['category'] = tree.find('.//s1sarl1:productClass', namespaces).text
        meta['sensor'] = tree.find('.//safe:familyName', namespaces).text.replace('ENTINEL-', '') + tree.find(
            './/safe:number', namespaces).text
        meta['IPF_version'] = float(tree.find('.//safe:software', namespaces).attrib['version'])

        return meta

    def unpack(self, directory, overwrite=False):
        outdir = os.path.join(directory, os.path.basename(self.file))
        self._unpack(outdir, overwrite=overwrite)


class TSX(ID):
    """
    Handler class for TerraSAR-X and TanDEM-X data

    References:
        TX-GS-DD-3302  TerraSAR-X Basic Product Specification Document
        TX-GS-DD-3303  TerraSAR-X Experimental Product Description
        TD-GS-PS-3028  TanDEM-X Experimental Product Description
        TerraSAR-X Image Product Guide (Airbus Defence and Space)
    
    Acquisition modes:
        ST:    Staring Spotlight
        HS:    High Resolution SpotLight
        HS300: High Resolution SpotLight 300 MHz
        SL:    SpotLight
        SM:    StripMap
        SC:    ScanSAR
        WS:    Wide ScanSAR
    
    Polarisation modes:
        Single (S): all acquisition modes
        Dual   (D): High Resolution SpotLight (HS), SpotLight (SL) and StripMap (SM)
        Twin   (T): StripMap (SM) (experimental)
        Quad   (Q): StripMap (SM) (experimental)
    
    Products:
        SSC: Single Look Slant Range Complex
        MGD: Multi Look Ground Range Detected
        GEC: Geocoded Ellipsoid Corrected
        EEC: Enhanced Ellipsoid Corrected
    """

    def __init__(self, scene):
        self.scene = os.path.realpath(scene)

        self.pattern = r'^(?P<sat>T[DS]X1)_SAR__' \
                       r'(?P<prod>SSC|MGD|GEC|EEC)_' \
                       r'(?P<var>____|SE__|RE__|MON1|MON2|BTX1|BRX2)_' \
                       r'(?P<mode>SM|SL|HS|HS300|ST|SC)_' \
                       r'(?P<pols>[SDTQ])_' \
                       r'(?:SRA|DRA)_' \
                       r'(?P<start>[0-9]{8}T[0-9]{6})_' \
                       r'(?P<stop>[0-9]{8}T[0-9]{6})(?:\.xml|)$'

        self.pattern_ds = r'^IMAGE_(?P<pol>HH|HV|VH|VV)_(?:SRA|FWD|AFT)_(?P<beam>[^\.]+)\.(cos|tif)$'
        self.examine(include_folders=False)

        if not re.match(re.compile(self.pattern), os.path.basename(self.file)):
            raise IOError('folder does not match TSX scene naming convention')

        self.meta = self.scanMetadata()
        self.meta['projection'] = 'GEOGCS["WGS 84",' \
                                  'DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],' \
                                  'PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],' \
                                  'UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],' \
                                  'AUTHORITY["EPSG","4326"]]'

        ID.__init__(self, self.meta)

    def scanMetadata(self):
        annotation = self.getFileObj(self.file)
        namespaces = getNamespaces(annotation)
        tree = ET.fromstring(annotation.read())
        annotation.close()
        meta = dict()
        meta['sensor'] = tree.find('.//generalHeader/mission', namespaces).text.replace('-', '')
        meta['product'] = tree.find('.//orderInfo/productVariant', namespaces).text
        meta['orbit'] = tree.find('.//missionInfo/orbitDirection', namespaces).text[0]
        meta['polarizations'] = [x.text for x in
                                 tree.findall('.//acquisitionInfo/polarisationList/polLayer', namespaces)]
        meta['orbit_abs'] = int(tree.find('.//missionInfo/absOrbit', namespaces).text)
        meta['orbit_rel'] = int(tree.find('.//missionInfo/relOrbit', namespaces).text)
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
        return meta

    def unpack(self, directory, overwrite=False):
        match = self.findfiles(self.pattern, True)
        header = [x for x in match if not x.endswith('xml') and 'iif' not in x][0].replace(self.scene, '').strip('/')
        outdir = os.path.join(directory, os.path.basename(header))
        self._unpack(outdir, offset=header, overwrite=overwrite)


class Archive(object):
    """
    Utility for storing SAR image metadata in a spatialite database
    """

    def __init__(self, dbfile):
        self.dbfile = dbfile
        self.conn = sqlite3.connect(self.dbfile)
        self.conn.enable_load_extension(True)
        try:
            self.conn.load_extension('mod_spatialite')
            if 'spatial_ref_sys' not in self.get_tablenames():
                self.conn.execute('SELECT InitSpatialMetaData(1);')
        except sqlite3.OperationalError:
            self.conn.load_extension('libspatialite')
            if 'spatial_ref_sys' not in self.get_tablenames():
                self.conn.execute('SELECT InitSpatialMetaData();')

        self.lookup = {'sensor': 'TEXT',
                       'orbit': 'TEXT',
                       'acquisition_mode': 'TEXT',
                       'start': 'TEXT',
                       'stop': 'TEXT',
                       'product': 'TEXT',
                       'samples': 'INTEGER',
                       'lines': 'INTEGER',
                       'outname_base': 'TEXT PRIMARY KEY',
                       'scene': 'TEXT',
                       'hh': 'INTEGER',
                       'vv': 'INTEGER',
                       'hv': 'INTEGER',
                       'vh': 'INTEGER'}
        create_string = '''CREATE TABLE if not exists data ({})'''.format(
            ', '.join([' '.join(x) for x in self.lookup.items()]))
        cursor = self.conn.cursor()
        cursor.execute(create_string)
        if 'bbox' not in self.get_colnames():
            cursor.execute('SELECT AddGeometryColumn("data","bbox" , 4326, "POLYGON", "XY", 0)')
        self.conn.commit()

        create_string = '''CREATE TABLE if not exists duplicates (outname_base TEXT, scene TEXT)'''
        cursor = self.conn.cursor()
        cursor.execute(create_string)
        self.conn.commit()

    def __prepare_insertion(self, scene):
        """
        read scene metadata and parse a string for inserting it into the database

        :param scene: a SAR scene
        :return: the actual insert string and a tuple containing parameters for the command, e.g.
        execute('''INSERT INTO data(a, b) VALUES(?, ?)''', (1, 2))
        where '?' is a placeholder for a value in the following tuple
        """
        id = scene if isinstance(scene, ID) else identify(scene)
        pols = [x.lower() for x in id.polarizations]
        insertion = []
        colnames = self.get_colnames()
        for attribute in colnames:
            if attribute == 'bbox':
                geom = id.bbox().convert2wkt(set3D=False)[0]
                insertion.append(geom)
            elif attribute in ['hh', 'vv', 'hv', 'vh']:
                insertion.append(int(attribute in pols))
            else:
                attr = getattr(id, attribute)
                value = attr() if inspect.ismethod(attr) else attr
                insertion.append(value)
        insert_string = '''INSERT INTO data({0}) VALUES({1})''' \
            .format(', '.join(colnames),
                    ', '.join(['GeomFromText(?, 4326)' if x == 'bbox' else '?' for x in colnames]))
        return insert_string, tuple(insertion)

    def insert(self, scene_in, verbose=False):
        """
        Insert one or many scenes into the database

        :param scene_in: a SAR scene or a list of scenes to be inserted
        :param verbose: should status information and a progress bar be printed into the console?
        :return: None
        """

        if isinstance(scene_in, (ID, str)):
            scenes = [scene_in if isinstance(scene_in, ID) else identify(scene_in)]
        elif isinstance(scene_in, list):
            if verbose:
                print('filtering scenes by name...')
            scenes = self.filter_scenelist(scene_in)
            if len(scenes) > 0:
                if verbose:
                    print('extracting scene metadata...')
                scenes = identify_many(scenes)
            else:
                print('all scenes are already registered')
                return
        else:
            raise RuntimeError('scene_in must either be a string pointing to a file, a pyroSAR.ID object '
                               'or a list containing several of either')
        duplicates_counter = 0
        if verbose:
            print('inserting scenes into the database...')
            pbar = pb.ProgressBar(maxval=len(scenes)).start()
        for i, id in enumerate(scenes):
            insert_string, insertion = self.__prepare_insertion(id)
            try:
                self.conn.execute(insert_string, insertion)
                self.conn.commit()
            except sqlite3.IntegrityError as e:
                if str(e) == 'UNIQUE constraint failed: data.outname_base':
                    self.conn.execute('INSERT INTO duplicates(outname_base, scene) VALUES(?, ?)',
                                      (id.outname_base(), id.scene))
                    self.conn.commit()
                    duplicates_counter += 1
                else:
                    raise e
            if verbose:
                pbar.update(i + 1)
        if verbose:
            pbar.finish()
        print('{} duplicates detected and registered'.format(duplicates_counter))

    def is_registered(self, scene):
        return len(self.select(scene=scene)) != 0 or len(self.select_duplicates(scene=scene)) != 0

    def export2shp(self, shp):
        """
        export the database to a shapefile

        :param shp: the name of the shapefile to be written
        :return: None
        """
        run(['ogr2ogr', '-f', '"ESRI Shapefile"', shp, self.dbfile])

    def filter_scenelist(self, scenelist):
        """
        filter a list of scenes by file names already registered in the database.

        Args:
            scenelist: a list of scenes (absolute path strings or pyroSAR.ID objects)

        Returns: a list which only contains files whose basename is not yet registered in the database

        """
        cursor = self.conn.execute('SELECT scene FROM data')
        registered = [os.path.basename(x[0].encode('ascii')) for x in cursor.fetchall()]
        cursor = self.conn.execute('SELECT scene FROM duplicates')
        duplicates = [os.path.basename(x[0].encode('ascii')) for x in cursor.fetchall()]
        return [x for x in scenelist if os.path.basename(x) not in registered+duplicates]

    def get_colnames(self):
        """
        return the names of the database table

        :return: a list containing the column names of the data table
        """
        cursor = self.conn.execute('''PRAGMA table_info(data)''')
        return [x[1].encode('ascii') for x in cursor.fetchall()]

    def get_tablenames(self):
        """
        return the names of all tables in the database

        :return: a list of table names
        """
        cursor = self.conn.execute('''SELECT * FROM sqlite_master WHERE type="table"''')
        return [x[1].encode('ascii') for x in cursor.fetchall()]

    def move(self, scenelist, directory):
        """
        move a list of files while keeping the database entries up to date
        if a scene is registered in the database (in either the data or duplicates table),
        the scene entry is directly changed to the new location

        :param scenelist: a list of file locations
        :param directory: a folder to which the files are moved
        :return: None
        """
        if not os.access(directory, os.W_OK):
            raise RuntimeError('directory cannot be written to')
        failed = []
        pbar = pb.ProgressBar(maxval=len(scenelist)).start()
        for i, scene in enumerate(scenelist):
            new = os.path.join(directory, os.path.basename(scene))
            try:
                shutil.move(scene, directory)
            except shutil.Error:
                failed.append(scene)
                continue
            finally:
                pbar.update(i + 1)
            if self.select(scene=scene) != 0:
                table = 'data'
            else:
                cursor = self.conn.execute('''SELECT scene FROM duplicates WHERE scene=?''', (scene,))
                if len(cursor.fetchall()) != 0:
                    table = 'duplicates'
                else:
                    table = None
            if table:
                self.conn.execute('''UPDATE {} SET scene=? WHERE scene=?'''.format(table), (new, scene))
                self.conn.commit()
        pbar.finish()
        if len(failed) > 0:
            print('the following scenes could not be moved:\n{}'.format('\n'.join(failed)))

    def select(self, vectorobject=None, mindate=None, maxdate=None, processdir=None, recursive=False, polarizations=None, **args):
        """
        select scenes from the database

        Args:
            vectorobject: an object of type spatial.vector.Vector
            mindate: a date string of format YYYYmmddTHHMMSS
            maxdate: a date string of format YYYYmmddTHHMMSS
            processdir: a directory to be scanned for already processed scenes; the selected scenes will be filtered to those that have not yet been processed
            recursive: should also the subdirectories of the processdir be scanned?
            **args: any further arguments (columns), which are registered in the database. See Archive.get_colnames()

        Returns: a list of strings pointing to the file locations of the selected scenes

        """
        arg_valid = [x for x in args.keys() if x in self.get_colnames()]
        arg_invalid = [x for x in args.keys() if x not in self.get_colnames()]
        if len(arg_invalid) > 0:
            print('the following arguments will be ignored as they are not registered in the data base: {}'.format(
                ', '.join(arg_invalid)))
        arg_format = []
        vals = []
        for key in arg_valid:
            if isinstance(args[key], (float, int, str)):
                arg_format.append('{0}="{1}"'.format(key, args[key]))
            elif isinstance(args[key], (tuple, list)):
                arg_format.append('{0} IN ("{1}")'.format(key, '", "'.join(map(str, args[key]))))
        if mindate:
            if re.search('[0-9]{8}T[0-9]{6}', mindate):
                arg_format.append('start>=?')
                vals.append(mindate)
            else:
                print('argument mindate is ignored, must be in format YYYYmmddTHHMMSS')
        if maxdate:
            if re.search('[0-9]{8}T[0-9]{6}', maxdate):
                arg_format.append('stop<=?')
                vals.append(maxdate)
            else:
                print('argument maxdate is ignored, must be in format YYYYmmddTHHMMSS')

        if polarizations:
            for pol in polarizations:
                if pol in ['HH', 'VV', 'HV', 'VH']:
                    arg_format.append('{}=1'.format(pol.lower()))

        if vectorobject:
            if isinstance(vectorobject, spatial.vector.Vector):
                vectorobject.reproject('+proj=longlat +datum=WGS84 +no_defs ')
                site_geom = vectorobject.convert2wkt(set3D=False)[0]
                arg_format.append('st_intersects(GeomFromText(?, 4326), bbox) = 1')
                vals.append(site_geom)
            else:
                print('argument vectorobject is ignored, must be of type spatial.vector.Vector')

        query = '''SELECT scene, outname_base FROM data WHERE {}'''.format(' AND '.join(arg_format))
        print(query)
        cursor = self.conn.execute(query, tuple(vals))
        if processdir:
            scenes = [x for x in cursor.fetchall()
                      if len(finder(processdir, [x[1]], regex=True, recursive=recursive)) == 0]
        else:
            scenes = cursor.fetchall()
        return [x[0].encode('ascii') for x in scenes]

    def select_duplicates(self, outname_base=None, scene=None):
        if not outname_base and not scene:
            cursor = self.conn.execute('SELECT * from duplicates')
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
            cursor = self.conn.execute(query, tuple(arg))
        return cursor.fetchall()

    @property
    def size(self):
        """
        get the number of scenes registered in the database

        :return: the number of scenes (integer)
        """
        cursor1 = self.conn.execute('''SELECT Count(*) FROM data''')
        cursor2 = self.conn.execute('''SELECT Count(*) FROM duplicates''')
        return (cursor1.fetchone()[0], cursor2.fetchone()[0])

    def __enter__(self):
        return self

    def close(self):
        self.conn.close()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.conn.close()
