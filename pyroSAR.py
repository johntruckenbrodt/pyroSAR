##############################################################
# Reading and Organizing system for SAR images
# John Truckenbrodt, Felix Cremer 2016-2017
##############################################################
"""
this script is intended to contain several SAR scene identifier classes to read basic metadata from the scene folders/files, convert to GAMMA format and do simple pre-processing
"""
import os
import re
import abc
import ssl
import struct
import math
import sqlite3
import StringIO
import numpy as np
import zipfile as zf
import tarfile as tf
import subprocess as sp
from urllib2 import urlopen, URLError
import xml.etree.ElementTree as ET
from time import strptime, strftime
from datetime import datetime, timedelta

from osgeo import gdal, ogr, osr
from osgeo.gdalconst import GA_ReadOnly, GA_Update

import envi
import gamma
import spatial
import linesimplify as ls
from xml_util import getNamespaces
from ancillary import finder, parse_literal, urlQueryParser, run


def identify(scene):
    """Return a metadata handler of the given scene."""
    for handler in ID.__subclasses__():
        try:
            return handler(scene)
        except (IOError, KeyError):
            pass
    raise IOError('data format not supported')


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
        locals = ['sensor', 'projection', 'orbit', 'polarizations', 'acquisition_mode', 'start', 'stop', 'product',
                  'spacing', 'samples', 'lines']
        for item in locals:
            setattr(self, item, metadict[item])

    def bbox(self, outname=None, overwrite=True):
        """Return the bounding box."""
        if outname is None:
            return spatial.bbox(self.getCorners(), self.projection)
        else:
            spatial.bbox(self.getCorners(), self.projection, outname=outname, format='ESRI Shapefile',
                         overwrite=overwrite)

    @abc.abstractmethod
    def calibrate(self, replace=False):
        raise NotImplementedError

    @property
    def compression(self):
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

    def export2sqlite(self, target=None):
        """
        Export the most important metadata in a sqlite database which is located in the same folder as the source file.
        """
        print 'Begin export'
        if self.compression is None:
            raise RuntimeError('Uncompressed data is not suitable for the metadata base')

        database = os.path.join(os.path.dirname(self.scene), 'data.db') if target is None else target
        conn = sqlite3.connect(database)
        # todo: following line doesn't work:
        conn.enable_load_extension(True)
        conn.execute('SELECT load_extension("libspatialite")')
        conn.execute('SELECT InitSpatialMetaData();')
        cursor = conn.cursor()
        create_string = '''CREATE TABLE if not exists data (
                            title TEXT NOT NULL,
                            file TEXT NOT NULL,
                            scene TEXT,
                            sensor TEXT,
                            projection TEXT,
                            orbit TEXT,
                            polarisation TEXT,
                            acquisition_mode TEXT,
                            start TEXT,
                            stop TEXT,
                            CONSTRAINT pk_data
                            PRIMARY KEY(file, scene)
                            )'''

        cursor.execute(create_string)
        cursor.execute('SELECT AddGeometryColumn("data","bbox" , 4326, "POLYGON", "XY", 0)')
        conn.commit()
        insert_string = '''
                        INSERT INTO data
                        (title, file, scene, sensor, projection, bbox, orbit, polarisation, acquisition_mode, start, stop)
                        VALUES( ?,?,?,?,?,GeomFromText(?, 4326),?,?,?,?,?)
                        '''
        geom = self.bbox().convert2wkt()[0]
        projection = spatial.crsConvert(self.projection, 'wkt')

        sq_file = os.path.basename(self.file)
        title = os.path.splitext(sq_file)[0]
        input = (title, sq_file, self.scene, self.sensor, projection, geom, self.orbit, 'polarisation', 'acquisition',
                 self.start, self.stop)
        try:
            cursor.execute(insert_string, input)
        except sqlite3.IntegrityError as e:
            print 'SQL error:', e

        conn.commit()
        conn.close()

    @abc.abstractmethod
    def convert2gamma(self, directory):
        raise NotImplementedError

    def examine(self, include_folders=False):
        files = self.findfiles(self.pattern, include_folders=include_folders)
        if len(files) == 1:
            self.file = files[0]
        elif len(files) == 0:
            raise IOError('folder does not match {} scene naming convention'.format(type(self).__name__))
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
        for a zipfile this is an StringIO.StringIO object (the zipfile.ExtFile object does not support setting file pointers via function 'seek', which is needed later on)
        """
        membername = filename.replace(self.scene, '').strip('/')

        if os.path.isdir(self.scene):
            obj = open(filename)
        elif zf.is_zipfile(self.scene):
            obj = StringIO.StringIO()
            with zf.ZipFile(self.scene, 'r') as zip:
                obj.write(zip.open(membername).read())
            obj.seek(0)

        elif tf.is_tarfile(self.scene):
            obj = StringIO.StringIO()
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
                           '%Y%m%d%H:%M:%S.%f']:
            try:
                return strftime('%Y%m%dT%H%M%S', strptime(x, timeformat))
            except (TypeError, ValueError):
                continue
        raise ValueError('unknown time format; check function ID.parse_date')

    def summary(self):
        for item in sorted(self.__dict__.keys()):
            if item != 'meta':
                print '{0}: {1}'.format(item, getattr(self, item))

    @abc.abstractmethod
    def scanMetadata(self):
        raise NotImplementedError

    @abc.abstractmethod
    def unpack(self, directory):
        raise NotImplementedError

    # todo: prevent unpacking if target files already exist
    # todo: replace with functionality from module archivist
    def _unpack(self, directory, offset=None):
        if not os.path.isdir(directory):
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
                            # note: the following is a pretty ugly workaround. Sentinel-1 Tiffs are occasionally provided with the wrong CRC-32 checksum although the file itself is intact.
                            # the command unzip unpacks the file first, then throws the error, while the Python zipfile module operates the other way around. If the Python module fails,
                            # the unzip command can still be operated and the file be used.
                            # todo: investigate how this can be done any better
                            except zf.BadZipfile:
                                cmd = ['unzip', '-j', '-qq', '-o', archive.filename, item, '-d',
                                       os.path.dirname(outname)]
                                try:
                                    sp.check_call(cmd)
                                except sp.CalledProcessError:
                                    continue
                archive.close()
            else:
                archive.extractall(directory)
                archive.close()
        self.scene = directory
        self.file = os.path.join(self.scene, os.path.basename(self.file))


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

    def unpack(self, directory):
        if self.sensor in ['ERS1', 'ERS2']:
            base_file = re.sub('\.PS$', '', os.path.basename(self.file))
            base_dir = os.path.basename(directory.strip('/'))

            outdir = directory if base_file == base_dir else os.path.join(directory, base_file)

            self._unpack(outdir)
        else:
            raise NotImplementedError('sensor {} not implemented yet'.format(self.sensor))

    def convert2gamma(self, directory):
        self.gammadir = directory
        if self.sensor in ['ERS1', 'ERS2']:
            if self.product == 'SLC' and self.meta['proc_system'] in ['PGS-ERS', 'VMP-ERS', 'SPF-ERS']:
                basename = '{}_{}_{}'.format(self.outname_base(), self.polarizations[0], self.product.lower())
                outname = os.path.join(directory, basename)
                if not os.path.isfile(outname):
                    lea = self.findfiles('LEA_01.001')[0]
                    dat = self.findfiles('DAT_01.001')[0]
                    title = re.sub('\.PS$', '', os.path.basename(self.file))
                    gamma.process(['par_ESA_ERS', lea, outname + '.par', dat, outname], inlist=[title])
                else:
                    print 'scene already converted'
            else:
                raise NotImplementedError(
                    'ERS {} product of {} processor in CEOS format not implemented yet'.format(self.product, self.meta[
                        'proc_system']))
        else:
            raise NotImplementedError('sensor {} in CEOS format not implemented yet'.format(self.sensor))

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


# id = identify('/geonfs01_vol1/ve39vem/archive/SAR/ERS/DRAGON/ERS1_0132_2529_20dec95.zip')


class CEOS_PSR(ID):
    """
    Handler class for ALOS PALSAR data in CEOS format

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
    """

    def __init__(self, scene):
        self.pattern = r'^LED-ALPSR' \
                       r'(?P<sub>P|S)' \
                       r'(?P<orbit>[0-9]{5})' \
                       r'(?P<frame>[0-9]{4})-' \
                       r'(?P<mode>[HWDPC])' \
                       r'(?P<level>1\.[015])' \
                       r'(?P<proc>G|_)' \
                       r'(?P<proj>[UPML_])' \
                       r'(?P<orbit_dir>A|D)$'

        self.scene = os.path.realpath(scene)

        self.examine()

        self.meta = self.scanMetadata()

        # register the standardized meta attributes as object attributes
        ID.__init__(self, self.meta)

    def calibrate(self, replace=False):
        for image in self.getGammaImages(self.scene):
            if image.endswith('_slc'):
                gamma.process(
                    ['radcal_SLC', image, image + '.par', image + '_cal', image + '_cal.par',
                     '-', '-', '-', '-', '-', '-', self.meta['k_dB']])
                envi.hdr(image + '_cal.par')

    def getLeaderfile(self):
        led_filename = self.findfiles(self.pattern)[0]
        led_obj = self.getFileObj(led_filename)
        led = led_obj.read()
        led_obj.close()
        return led

    def scanMetadata(self):
        led_filename = self.findfiles(self.pattern)[0]
        led_obj = self.getFileObj(led_filename)
        led = led_obj.read()
        led_obj.close()

        p0 = 0
        p1 = 720
        fileDescriptor = led[p0:p1]
        # dss_n = int(fileDescriptor[180:186])
        mpd_n = int(fileDescriptor[192:198])
        # ppd_n = int(fileDescriptor[204:210])
        # adr_n = int(fileDescriptor[216:222])
        # rdr_n = int(fileDescriptor[228:234])

        p0 = p1
        p1 += 4096
        dataSetSummary = led[p0:p1]

        if mpd_n == 1:
            p0 = p1
            p1 += 1620
            mapProjectionData = led[p0:p1]
        else:
            mapProjectionData = None

        p0 = p1
        p1 += 4680
        platformPositionData = led[p0:p1]
        p0 = p1
        p1 += 8192
        attitudeData = led[p0:p1]
        p0 = p1
        p1 += 9860
        radiometricData = led[p0:p1]
        p0 = p1
        p1 += 1620
        dataQualitySummary = led[p0:p1]

        facilityRelatedData = [0] * 11
        for i in range(0, 10):
            p0 = p1
            length = struct.unpack('>i', led[(p0 + 8):(p0 + 12)])[0]
            p1 += length
            facilityRelatedData[i] = led[p0:p1]

        facilityRelatedData[10] = led[p1:]

        meta = {}
        meta['sensor'] = 'PSR1'
        rlks = float(dataSetSummary[1190:1206])
        azlks = float(dataSetSummary[1174:1190])
        meta['looks'] = (rlks, azlks)
        meta['incidence'] = float(dataSetSummary[484:492])
        spacing_range = float(dataSetSummary[1702:1718])
        spacing_azimuth = float(dataSetSummary[1686:1702])
        meta['spacing'] = (spacing_range, spacing_azimuth)

        match = re.match(re.compile(self.pattern), os.path.basename(led_filename))
        meta['orbit'] = match.group('orbit_dir')
        meta['acquisition_mode'] = match.group('sub') + match.group('mode')
        meta['product'] = match.group('level')
        meta['proc_facility'] = dataSetSummary[1046:1062].strip()
        meta['proc_system'] = dataSetSummary[1062:1070].strip()
        meta['proc_version'] = dataSetSummary[1070:1078].strip()
        meta['samples'] = int(dataSetSummary[332:340]) * 2
        meta['lines'] = int(dataSetSummary[324:332]) * 2
        meta['wavelength'] = float(dataSetSummary[500:516]) * 100  # in cm
        try:
            start1, start2 = re.search('Img_SceneStartDateTime[ ="0-9:.]*', led).span()
            stop1, stop2 = re.search('Img_SceneEndDateTime[ ="0-9:.]*', led).span()
        except AttributeError:
            raise IndexError('start and stop time stamps cannot be extracted; see file {}'.format(led_filename))
        start_string = led[start1:start2]
        stop_string = led[stop1:stop2]
        meta['start'] = self.parse_date(re.sub('[a-zA-Z =_"]+', '', start_string))
        meta['stop'] = self.parse_date(re.sub('[a-zA-Z =_"]+', '', stop_string))
        meta['polarizations'] = [re.search('[HV]{2}', os.path.basename(x)).group(0) for x in self.findfiles('^IMG-')]
        meta['projection'] = 'GEOGCS["WGS 84",' \
                             'DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],' \
                             'PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],' \
                             'UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],' \
                             'AUTHORITY["EPSG","4326"]]'
        meta['k_dB'] = float(radiometricData[20:36])
        return meta

        # id = CEOS_PSR('/geonfs01_vol1/ve39vem/archive/SAR/PALSAR-1/ALPSRP224031000-H1.1__A.zip')

    def convert2gamma(self, directory):
        images = self.findfiles('^IMG-')
        if self.product == '1.0':
            raise RuntimeError('PALSAR-1 level 1.0 products are not supported')
        for image in images:
            polarization = re.search('[HV]{2}', os.path.basename(image)).group(0)
            if self.product == '1.1':
                outname_base = '{}_{}_slc'.format(self.outname_base(), polarization)
                outname = os.path.join(directory, outname_base)
                gamma.process(['par_EORC_PALSAR', self.file, outname + '.par', image, outname])
            else:
                outname_base = '{}_{}_mli_geo'.format(self.outname_base(), polarization)
                outname = os.path.join(directory, outname_base)
                gamma.process(
                    ['par_EORC_PALSAR_geo', self.file, outname + '.par', outname + '_dem.par', image, outname])
            envi.hdr(outname + '.par')

    def unpack(self, directory):
        outdir = os.path.join(directory, os.path.basename(self.file).replace('LED-', ''))
        self._unpack(outdir)

    def getCorners(self):

        img_filename = self.findfiles('IMG')[0]
        img_obj = self.getFileObj(img_filename)
        imageFileDescriptor = img_obj.read(720)

        lineRecordLength = int(imageFileDescriptor[186:192])  # byte per line + 412
        numberOfRecords = int(imageFileDescriptor[180:186])

        signalDataDescriptor1 = img_obj.read(412)
        img_obj.seek(720 + lineRecordLength * (numberOfRecords-1))
        signalDataDescriptor2 = img_obj.read()

        img_obj.close()

        y11 = struct.unpack('>i', signalDataDescriptor1[192:196])[0] / 1000000.
        y12 = struct.unpack('>i', signalDataDescriptor1[200:204])[0] / 1000000.
        x11 = struct.unpack('>i', signalDataDescriptor1[204:208])[0] / 1000000.
        x12 = struct.unpack('>i', signalDataDescriptor1[212:216])[0] / 1000000.

        y21 = struct.unpack('>i', signalDataDescriptor2[192:196])[0] / 1000000.
        y22 = struct.unpack('>i', signalDataDescriptor2[200:204])[0] / 1000000.
        x21 = struct.unpack('>i', signalDataDescriptor2[204:208])[0] / 1000000.
        x22 = struct.unpack('>i', signalDataDescriptor2[212:216])[0] / 1000000.

        lat = [y11, y12, y21, y22]
        lon = [x11, x12, x21, x22]

        return {'xmin': min(lon), 'xmax': max(lon), 'ymin': min(lat), 'ymax': max(lat)}


# id = CEOS_PSR('/geonfs01_vol1/ve39vem/archive/SAR/PALSAR-1/ALPSRP224031000-H1.1__A.zip')


class ESA(ID):
    """
    Handler class for SAR data in ESA format (Envisat, ERS)
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

    # todo: prevent conversion if target files already exists
    def convert2gamma(self, directory):
        """
        the command par_ASAR also accepts a K_dB argument in which case the resulting image names will carry the suffix GRD;
        this is not implemented here but instead in method calibrate
        """
        self.gammadir = directory
        outname = os.path.join(directory, self.outname_base())
        if len(self.getGammaImages(directory)) == 0:
            gamma.process(['par_ASAR', os.path.basename(self.file), outname], os.path.dirname(self.file))
            os.remove(outname + '.hdr')
            for item in finder(directory, [os.path.basename(outname)], regex=True):
                ext = '.par' if item.endswith('.par') else ''
                base = os.path.basename(item).strip(ext)
                base = base.replace('.', '_')
                base = base.replace('PRI', 'pri')
                base = base.replace('SLC', 'slc')
                newname = os.path.join(directory, base + ext)
                os.rename(item, newname)
                if newname.endswith('.par'):
                    envi.hdr(newname)
        else:
            raise IOError('scene already processed')

    def calibrate(self, replace=False):
        k_db = {'ASAR': 55., 'ERS1': 58.24, 'ERS2': 59.75}[self.sensor]
        inc_ref = 90. if self.sensor == 'ASAR' else 23.
        # candidates = [x for x in self.getGammaImages(self.gammadir) if not re.search('_(?:cal|grd)$', x)]
        candidates = [x for x in self.getGammaImages(self.gammadir) if re.search('_pri$', x)]
        for image in candidates:
            out = image.replace('pri', 'grd')
            gamma.process(['radcal_PRI', image, image + '.par', out, out + '.par', k_db, inc_ref])
            envi.hdr(out + '.par')
            if replace:
                for item in [image, image + '.par', image + '.hdr']:
                    if os.path.isfile(item):
                        os.remove(item)
                        # candidates = [x for x in self.getGammaImages(self.gammadir) if re.search('_slc$', x)]
                        # for image in candidates:
                        #     par = gamma.ISPPar(image+'.par')
                        #     out = image+'_cal'
                        #     fcase = 1 if par.image_format == 'FCOMPLEX' else 3
                        #     gamma.process(['radcal_SLC', image, image + '.par', out, out + '.par', fcase, '-', '-', '-', '-', '-', k_db, inc_ref])
                        #     envi.hdr(out + '.par')
                        #     if replace:
                        #         for item in [image, image+'.par', image+'.hdr']:
                        #             if os.path.isfile(item):
                        #                 os.remove(item)

    def unpack(self, directory):
        base_file = os.path.basename(self.file).strip('\.zip|\.tar(?:\.gz|)')
        base_dir = os.path.basename(directory.strip('/'))

        outdir = directory if base_file == base_dir else os.path.join(directory, base_file)

        self._unpack(outdir)


# class RS2(ID):
#     def __init__(self, scene):
#
#         raise IOError
#
#         self.pattern = r'^(?:RS2|RSAT2)_(?:OK[0-9]+)_(?:PK[0-9]+)_(?:DK[0-9]+)_' \
#                        r'(?P<beam>[0-9A-Z]+)_' \
#                        r'(?P<date>[0-9]{8})_' \
#                        r'(?P<time>[0-9]{6})_' \
#                        r'(?P<pols>[HV]{2}_' \
#                        r'(?P<level>SLC|SGX|SGF|SCN|SCW|SSG|SPG)$'
#
#         self.sensor = 'RS2'
#         self.scene = os.path.realpath(scene)
#         self.gdalinfo(self.scene)
#         self.start = self.ACQUISITION_START_TIME
#         self.incidence = (self.FAR_RANGE_INCIDENCE_ANGLE + self.NEAR_RANGE_INCIDENCE_ANGLE)/2
#         self.spacing = (self.PIXEL_SPACING, self.LINE_SPACING)
#         self.orbit = self.ORBIT_DIRECTION[0]
#
#     def getCorners(self):
#         lat = [x[1][1] for x in self.gcps]
#         lon = [x[1][0] for x in self.gcps]
#         return {'xmin': min(lon), 'xmax': max(lon), 'ymin': min(lat), 'ymax': max(lat)}

# id = identify('/geonfs01_vol1/ve39vem/RS2/RS2_OK53107_PK504800_DK448361_FQ1_20140606_055403_HH_VV_HV_VH_SLC.zip')


# todo: check self.file and self.scene assignment after unpacking
class SAFE(ID):
    """
    Handler class for Sentinel-1 data
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
            print subset
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

    def calibrate(self, replace=False):
        print 'calibration already performed during import'

    def convert2gamma(self, directory, noiseremoval=True):
        if self.compression is not None:
            raise RuntimeError('scene is not yet unpacked')
        if self.product == 'OCN':
            raise IOError('Sentinel-1 OCN products are not supported')
        if self.meta['category'] == 'A':
            raise IOError('Sentinel-1 annotation-only products are not supported')

        if not os.path.isdir(directory):
            os.makedirs(directory)

        for xml_ann in finder(os.path.join(self.scene, 'annotation'), [self.pattern_ds], regex=True):
            base = os.path.basename(xml_ann)
            match = re.compile(self.pattern_ds).match(base)

            tiff = os.path.join(self.scene, 'measurement', base.replace('.xml', '.tiff'))
            xml_cal = os.path.join(self.scene, 'annotation', 'calibration', 'calibration-' + base)

            product = match.group('product')

            # specify noise calibration file
            # L1 GRD product: thermal noise already subtracted, specify xml_noise to add back thermal noise
            # SLC products: specify noise file to remove noise
            # xml_noise = '-': noise file not specified
            if (noiseremoval and product == 'slc') or (not noiseremoval and product == 'grd'):
                xml_noise = os.path.join(self.scene, 'annotation', 'calibration', 'noise-' + base)
            else:
                xml_noise = '-'

            fields = (self.outname_base(),
                      match.group('pol').upper(),
                      product)
            name = os.path.join(directory, '_'.join(fields))

            if product == 'slc':
                swath = match.group('swath').upper()
                name = name.replace('{:_<{l}}'.format(self.acquisition_mode, l=len(swath)), swath)
                cmd = ['par_S1_SLC', tiff, xml_ann, xml_cal, xml_noise, name + '.par', name, name + '.tops_par']
            else:
                cmd = ['par_S1_GRD', tiff, xml_ann, xml_cal, xml_noise, name + '.par', name]

            gamma.process(cmd)
            envi.hdr(name + '.par')
            self.gammafiles[product].append(name)

    def correctOSV(self, osvdir=None):
        logdir = os.path.join(self.scene, 'logfiles')
        if not os.path.isdir(logdir):
            os.makedirs(logdir)
        if osvdir is None:
            osvdir = os.path.join(self.scene, 'osv')
        if not os.path.isdir(osvdir):
            os.makedirs(osvdir)
        try:
            self.getOSV(osvdir)
        except URLError:
            print '..no internet access'
        for image in self.getGammaImages(self.scene):
            gamma.process(['OPOD_vec', image + '.par', osvdir], outdir=logdir)

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

    def unpack(self, directory):
        outdir = os.path.join(directory, os.path.basename(self.file))
        self._unpack(outdir)


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
        Dual (D):   High Resolution SpotLight (HS), SpotLight (SL) and StripMap (SM)
        Twin (T):   StripMap (SM) (experimental)
        Quad (Q):   StripMap (SM) (experimental)
    
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

    def convert2gamma(self, directory):
        images = self.findfiles(self.pattern_ds)
        pattern = re.compile(self.pattern_ds)
        for image in images:
            pol = pattern.match(os.path.basename(image)).group('pol')
            outname = os.path.join(directory, self.outname_base() + '_' + pol)
            if self.product == 'SSC':
                outname += '_slc'
                gamma.process(['par_TX_SLC', self.file, image, outname + '.par', outname, pol])
            elif self.product == 'MGD':
                outname += '_mli'
                gamma.process(['par_TX_GRD', self.file, image, outname + '.par', outname, pol])
            else:
                outname += '_mli_geo'
                gamma.process(['par_TX_geo', self.file, image, outname + '.par', outname + '_dem.par', outname, pol])
            envi.hdr(outname + '.par')

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

    def unpack(self, directory):
        match = self.findfiles(self.pattern, True)
        header = [x for x in match if not x.endswith('xml') and 'iif' not in x][0].replace(self.scene, '').strip('/')
        outdir = os.path.join(directory, os.path.basename(header))
        self._unpack(outdir, header)


# id = identify('/geonfs01_vol1/ve39vem/archive/SAR/TerraSAR-X/TSX1_SAR__MGD_SE___SL_S_SRA_20110902T015248_20110902T015249.zip')

# todo: add method export2shp
class Archive(object):
    """
    Utility for storing relevant SAR image metadata in a CSV file based database
    """

    def __init__(self, scenelist, header=False, keys=None):
        self.scenelist = scenelist
        self.reg = {}
        if keys is None:
            self.keys = ['sensor', 'acquisition_mode', 'polarizations', 'scene', 'bbox']
        if os.path.isfile(self.scenelist):
            self.file = open(scenelist, 'a+', 0)
            if header:
                self.keys = self.file.readline().strip().split(';')
            for line in self.file:
                items = dict(zip(self.keys, line.strip().split(';')))
                base = os.path.basename(items['scene'])
                self.reg[base] = items
            self.file.seek(0)
        else:
            self.file = open(scenelist, 'w', 0)
            if header:
                self.file.write(';'.join(self.keys) + '\n')

    def __enter__(self):
        return self

    def update(self, scenes):
        for scene in scenes:
            base = os.path.basename(scene)
            if base not in self.reg:
                try:
                    id = identify(scene)
                except IOError:
                    print 'failed:', base
                    continue
                print base
                items = [id.sensor, id.acquisition_mode, ','.join(id.polarizations), id.scene,
                         id.bbox().convert2wkt()[0]]
                self.file.write(';'.join(items) + '\n')
                self.reg[base] = dict(zip(self.keys, items))

    def select(self, vectorobject):
        vectorobject.reproject('+proj=longlat +datum=WGS84 +no_defs ')
        site_geom = ogr.CreateGeometryFromWkt(vectorobject.convert2wkt()[0])
        selection_site = []
        for entry in self.reg:
            geom = ogr.CreateGeometryFromWkt(self.reg[entry]['bbox'])
            intersection = geom.Intersection(site_geom)
            if intersection.GetArea() > 0:
                selection_site.append(self.reg[entry]['scene'])
        return selection_site

    def delete(self, filename):
        for entry in self.reg:
            if self.reg[entry]['scene'] == filename:
                del self.reg[entry]
                break
        with open(self.scenelist, 'r+') as f:
            lines = f.readlines()
            f.seek(0)
            for line in lines:
                if not re.search(filename, line):
                    f.write(line)
            f.truncate()
        if os.path.isfile(filename):
            print filename
            os.remove(filename)

    @property
    def size(self):
        return len(self.reg)

    def close(self):
        self.file.close()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.file.close()
