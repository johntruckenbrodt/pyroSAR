##############################################################
# general utilities for Sentinel-1
# John Truckenbrodt 2016-2017
##############################################################
import os
import re
import ssl
import sys
import time
from datetime import datetime
from urllib import urlopen, urlencode
from urlparse import urlparse, urlunparse

from .. import gamma
from ..ancillary import finder, dissolve

try:
    import argparse
except ImportError:
    try:
        os.remove(os.path.join(os.path.dirname(sys.argv[0]), 'locale.pyc'))
    finally:
        import argparse


def init_parser():
    """
    initialize argument parser for S1 processing utilities
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--transform', action='store_true', help='transform the final DEM to UTM coordinates')
    parser.add_argument('-l', '--logfiles', action='store_true', help='create logfiles of the executed GAMMA commands')
    parser.add_argument('-i', '--intermediates', action='store_true', help='keep intermediate files')
    parser.add_argument('-q', '--quiet', action='store_true', help='suppress standard console prints')
    parser.add_argument('-tr', '--targetresolution', default=20, help='the target resolution in meters for x and y',
                        type=int)
    parser.add_argument('-fg', '--func_geoback', default=2, help='backward geocoding interpolation function; '
                                                                 '0 - Nearest Neighbor, 1 - Bicubic Spline, 2 - Bicubic Spline-Log; '
                                                                 'method 1: negative values possible (e.g. in urban areas) - use method 2 to avoid this',
                        type=int)
    parser.add_argument('-fi', '--func_interp', default=0,
                        help='function for interpolation of layover/shadow/foreshortening/DEM gaps; '
                             '0 - set to 0, 1 - linear interpolation, 2 - actual value, 3 - nn-thinned', type=int)
    parser.add_argument('-poe', '--poedir', default=None,
                        help='directory containing aux_poeorb (precise orbit ephemerides) orbit state vector files')
    parser.add_argument('-res', '--resdir', default=None,
                        help='directory containing aux_resorb (restituted orbit) orbit state vector files')
    parser.add_argument('zipfile', help='S1 zipped scene archive to be used')
    parser.add_argument('tempdir', help='temporary directory for intermediate files')
    parser.add_argument('outdir', help='output directory')
    parser.add_argument('srtmdir', help='directory containing SRTM hgt tiles (subdirectories possible)')
    return parser


def deburst(burst1, burst2, burst3, name_out, rlks=5, azlks=1, replace=False, path_log=None):
    """
    debursting of S1 SLC imagery
    the procedure consists of two steps. First antenna pattern deramping and then mosaicing of the single deramped bursts
    for mosaicing, the burst boundaries are calculated from the number of looks in range (rlks) and azimuth (azlks), in this case 5 range looks and 1 azimuth looks.
    Alternately 10 range looks and 2 azimuth looks could be used.
    if replace is set to True, the original files will be deleted
    """
    for burst in [burst1, burst2, burst3]:
        if not os.path.isfile(burst) or not os.path.isfile(burst + '.par') or not os.path.isfile(burst + '.tops_par'):
            raise IOError('input files missing; parameter files must be named e.g. {burst1}.par and {burst1}.tops_par')
    outpath = os.path.dirname(name_out)
    if not os.path.isdir(outpath):
        os.makedirs(outpath)
    tab_in = os.path.join(outpath, 'tab_deramp1')
    tab_out = os.path.join(outpath, 'tab_deramp2')
    with open(tab_in, 'w') as out1:
        with open(tab_out, 'w') as out2:
            for item in [burst1, burst2, burst3]:
                out1.write(item + '\t' + item + '.par\t' + item + '.tops_par\n')
                out2.write(item + '_drp\t' + item + '_drp.par\t' + item + '_drp.tops_par\n')
    gamma.process(['SLC_deramp_S1_TOPS', tab_in, tab_out, 0, 0], logpath=path_log)
    gamma.process(['SLC_mosaic_S1_TOPS', tab_out, name_out, name_out + '.par', rlks, azlks], logpath=path_log)

    if replace:
        for item in [burst1, burst2, burst3]:
            for subitem in [item + x for x in ['', '.par', '.tops_par']]:
                os.remove(subitem)
    for item in [burst1, burst2, burst3]:
        for subitem in [item + x for x in ['_drp', '_drp.par', '_drp.tops_par']]:
            os.remove(subitem)
    os.remove(tab_in)
    os.remove(tab_out)


class OSV(object):
    """
    interface for management of S1 Orbit State Vector (OSV) files
    input are two directories, one for Precise Orbit Ephemerides (POE) and one for Restituted Orbit (RES) files; these directories are created if they do not exist
    actions performed upon calling the main function 'update':
    -the ESA Quality Control (QC) server is checked for any POE files not in the local directory
    -POE  files on the server and not in the local directory are downloaded
    -RES files newer than the latest POE file are downloaded; POE files are approximately 18 days behind the actual date, thus RES files can be used instead
    -delete all RES files for whose date a POE file has become available
    using function 'match' the corresponding POE (priority) or RES file is returned for a timestamp
    timestamps are always handled in the format YYYYMMDDThhmmss
    """
    def __init__(self, outdir_poe, outdir_res):
        self.remote_poe = 'https://qc.sentinel1.eo.esa.int/aux_poeorb/'
        self.remote_res = 'https://qc.sentinel1.eo.esa.int/aux_resorb/'
        if outdir_poe == outdir_res:
            raise IOError('POE and RES directories must be different')
        self.outdir_poe = outdir_poe
        self.outdir_res = outdir_res
        self.pattern = 'S1[AB]_OPER_AUX_(?:POE|RES)ORB_OPOD_[0-9TV_]{48}\.EOF'
        self.pattern_fine = 'S1[AB]_OPER_AUX_(?P<type>(?:POE|RES)ORB)_OPOD_(?P<publish>[0-9]{8}T[0-9]{6})_V(?P<start>[0-9]{8}T[0-9]{6})_(?P<stop>[0-9]{8}T[0-9]{6})\.EOF'
        self.sslcontext = ssl._create_unverified_context()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return

    def _init_dir(self):
        """
        create directories if they don't exist yet
        """
        for dir in [self.outdir_poe, self.outdir_res]:
            if not os.path.isdir(dir):
                os.makedirs(dir)

    def _typeEvaluate(self, type):
        """
        evaluate the 'type' function argument and return the corresponding local directory
        """
        if type not in ['POE', 'RES']:
            raise IOError('type must be either "POE" or "RES"')
        if type == 'POE':
            return self.remote_poe, self.outdir_poe
        else:
            return self.remote_res, self.outdir_res

    def catch(self, type='POE', start=None, stop=None):
        """
        check a server for files
        """
        address, outdir = self._typeEvaluate(type)
        address_parse = urlparse(address)
        query = {'page': 1}
        files = []
        if start is not None:
            date_start = datetime.strptime(start, '%Y%m%dT%H%M%S').strftime('%Y-%m-%d')
        else:
            date_start = '2014-08-22'
        if stop is not None:
            date_stop = datetime.strptime(stop, '%Y%m%dT%H%M%S').strftime('%Y-%m-%d')
        else:
            date_stop = time.strftime('%Y-%m-%d')
        query['validity_start_time'] = '{0}..{1}'.format(date_start, date_stop)
        print 'searching for new {} files'.format(type)
        while True:
            subaddress = urlunparse(address_parse._replace(query=urlencode(query)))
            try:
                response = urlopen(subaddress, context=self.sslcontext).read()
                print subaddress
            except IOError as e:
                raise RuntimeError(e)
            remotes = [os.path.join(address, x) for x in sorted(set(re.findall(self.pattern, response)))]
            if start is not None:
                remotes = [x for x in remotes if self.date(x, 'stop') > start]
            selection = [x for x in remotes if x not in files]
            if len(selection) == 0:
                break
            else:
                files += selection
                query['page'] += 1
        if type == 'RES':
            files = [x for x in files if self.date(x, 'stop') > self.maxdate('POE', 'stop')]
        return files

    def date(self, file, type):
        """
        extract a date from an OSV file name; types: 'publish', 'start', 'stop'
        """
        return re.match(self.pattern_fine, os.path.basename(file)).group(type)

    def clean_res(self):
        """
        delete all RES files for whose date a POE file exists
        """
        maxdate_poe = self.maxdate('POE', 'stop')
        depreceated = [x for x in self.getLocals('RES') if self.date(x, 'stop') < maxdate_poe]
        print 'deleting {0} RES files'.format(len(depreceated))
        for item in depreceated:
            os.remove(item)

    def getLocals(self, type='POE'):
        """
        get a list of local files
        """
        address, directory = self._typeEvaluate(type)
        return finder(directory, [self.pattern], regex=True)

    def maxdate(self, type='POE', datetype='stop'):
        """
        return the latest date of POE/RES files; datetypes: 'publish', 'start', 'stop'
        """
        address, directory = self._typeEvaluate(type)
        files = finder(directory, [self.pattern], regex=True)
        return max([self.date(x, datetype) for x in files]) if len(files) > 0 else None

    def mindate(self, type='POE', datetype='start'):
        """
        return the latest date of POE/RES files; datetypes: 'publish', 'start', 'stop'
        """
        address, directory = self._typeEvaluate(type)
        files = finder(directory, [self.pattern], regex=True)
        return min([self.date(x, datetype) for x in files]) if len(files) > 0 else None

    def match(self, timestamp):
        """
        return the corresponding OSV file for the provided time stamp
        """
        for item in dissolve([self.getLocals('POE'), self.getLocals('RES')]):
            if self.date(item, 'start') <= timestamp <= self.date(item, 'stop'):
                return item
        return None

    def retrieve(self, files, type='POE'):
        """
        download the newest files
        """
        address, outdir = self._typeEvaluate(type)
        if not os.access(outdir, os.W_OK):
            raise RuntimeError('insufficient directory permissions')
        downloads = [x for x in files if not os.path.isfile(os.path.join(outdir, os.path.basename(x)))]
        print 'downloading {0} {1} files to {2}'.format(len(downloads), type, outdir)
        for item in downloads:
            infile = urlopen(item, context=self.sslcontext)
            with open(os.path.join(outdir, os.path.basename(item)), 'wb') as outfile:
                outfile.write(infile.read())
            infile.close()

    def update(self):
        """
        perform creating/updating operations for POE and RES files: 
        download newest POE and RES files, delete RES files which can be replaced by newly downloaded POE files
        """
        self._init_dir()
        try:
            files_poe = self.catch('POE', start=self.maxdate('POE', 'start'))
        except RuntimeError:
            print 'no internet connection'
            return
        self.retrieve(files_poe, 'POE')
        print '---------------------------------------------------------'
        files_res = self.catch('RES', start=self.maxdate('RES', 'start'))
        self.retrieve(files_res, 'RES')
        self.clean_res()
