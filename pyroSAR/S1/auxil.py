##############################################################
# general utilities for Sentinel-1
# John Truckenbrodt 2016-2018
##############################################################
import sys

if sys.version_info >= (3, 0):
    from urllib.request import urlopen
else:
    from urllib import urlopen

import os
import re
import ssl
import time
from datetime import datetime

from ..ancillary import finder, urlQueryParser

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


# todo check existence not by file name but by start and stop time; files are sometimes re-published
class OSV(object):
    """
    interface for management of S1 Orbit State Vector (OSV) files

    input is a directory which is supposed to contain, or already contains, OSV files.
    Two subdirectories are expected and created otherwise:
    one for Precise Orbit Ephemerides (POE) named POEORB and one for Restituted Orbit (RES) files named RESORB

    Using method :meth:`match` the corresponding POE (priority) or RES file is returned for a timestamp.
    Timestamps are always handled in the format YYYYmmddTHHMMSS.

    Parameters
    ----------
    osvdir: str
        the directory to write the orbit files to
    """
    def __init__(self, osvdir):
        self.remote_poe = 'https://qc.sentinel1.eo.esa.int/aux_poeorb/'
        self.remote_res = 'https://qc.sentinel1.eo.esa.int/aux_resorb/'
        self.outdir_poe = os.path.join(osvdir, 'POEORB')
        self.outdir_res = os.path.join(osvdir, 'RESORB')
        self.pattern = 'S1[AB]_OPER_AUX_(?:POE|RES)ORB_OPOD_[0-9TV_]{48}\.EOF'
        self.pattern_fine = 'S1[AB]_OPER_AUX_' \
                            '(?P<type>(?:POE|RES)ORB)_OPOD_' \
                            '(?P<publish>[0-9]{8}T[0-9]{6})_V' \
                            '(?P<start>[0-9]{8}T[0-9]{6})_' \
                            '(?P<stop>[0-9]{8}T[0-9]{6})\.EOF'
        if sys.version_info >= (2, 7, 9):
            self.sslcontext = ssl._create_unverified_context()
        else:
            raise RuntimeError('this functionality requires Python Version >=2.7.9')

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

    def _typeEvaluate(self, osvtype):
        """
        evaluate the 'osvtype' method argument and return the corresponding remote repository and local directory
        :param osvtype: the type of orbit files required; either 'POE' or 'RES'
        :return: the remote repository and local directory of the osv type
        """
        if osvtype not in ['POE', 'RES']:
            raise IOError('type must be either "POE" or "RES"')
        if osvtype == 'POE':
            return self.remote_poe, self.outdir_poe
        else:
            return self.remote_res, self.outdir_res

    def catch(self, osvtype='POE', start=None, stop=None):
        """
        check a server for files

        Parameters
        ----------
        osvtype: {'POE', 'RES'}
            the type of orbit files required
        start: str
            the date to start searching for files
        stop: str
            the date to stop searching for files

        Returns
        -------
        list
            the URLs of the remote OSV files
        """
        address, outdir = self._typeEvaluate(osvtype)
        # a dictionary for storing the url arguments
        query = {'page': 1}
        # the collection of files to be returned
        files = []
        # set the defined date or the date of the first existing OSV file otherwise
        if start is not None:
            date_start = datetime.strptime(start, '%Y%m%dT%H%M%S').strftime('%Y-%m-%d')
        else:
            date_start = '2014-08-22'
        # set the defined date or the current date otherwise
        if stop is not None:
            date_stop = datetime.strptime(stop, '%Y%m%dT%H%M%S').strftime('%Y-%m-%d')
        else:
            date_stop = time.strftime('%Y-%m-%d')

        # pattern for scanning urlopen response for links to OSV files
        pattern_url = 'http.*{}'.format(self.pattern)

        # append the time frame to the query dictionary
        query['validity_start'] = '{0}..{1}'.format(date_start, date_stop)
        print('searching for new {} files'.format(osvtype))
        # iterate through the url pages and look for files
        while True:
            # parse the url
            subaddress = urlQueryParser(address, query)
            # read the remote content
            try:
                response = urlopen(subaddress, context=self.sslcontext).read().decode('utf-8')
                print(subaddress)
            except IOError as e:
                raise RuntimeError(e)
            # list all osv files found on the page
            remotes = sorted(set(re.findall(pattern_url, response)))
            # do a more accurate filtering of the time stamps
            if start is not None:
                remotes = [x for x in remotes if self.date(x, 'stop') > start]
            # filter files already existing in the files collection
            selection = [x for x in remotes if x not in files]
            # stop the loop if no more files are found on the current url page
            if len(selection) == 0:
                break
            else:
                # append the found files to the collection and increment the url page
                files += selection
                query['page'] += 1
        # in case the type 'RES' is selected then only return those files covering
        # a time period not covered by any POE file
        if osvtype == 'RES':
            files = [x for x in files if self.date(x, 'stop') > self.maxdate('POE', 'stop')]
        return files

    def date(self, file, datetype):
        """
        extract a date from an OSV file name

        Parameters
        ----------
        datetype: {'publish', 'start', 'stop'}
            one of three possible date types contained in the OSV filename

        Returns
        -------
        str
            a time stamp in the format YYYYmmddTHHMMSS
        """
        return re.match(self.pattern_fine, os.path.basename(file)).group(datetype)

    def clean_res(self):
        """
        delete all RES files for whose date a POE file exists
        """
        maxdate_poe = self.maxdate('POE', 'stop')
        deprecated = [x for x in self.getLocals('RES') if self.date(x, 'stop') < maxdate_poe]
        print('deleting {0} RES files'.format(len(deprecated)))
        for item in deprecated:
            os.remove(item)

    def getLocals(self, osvtype='POE'):
        """
        get a list of local files of specific type
        Parameters
        ----------
        osvtype: {'POE', 'RES'}
            the type of orbit files required
        Returns
        -------
        list
            a selection of local OSV files
        """
        address, directory = self._typeEvaluate(osvtype)
        return finder(directory, [self.pattern], regex=True)

    def maxdate(self, osvtype='POE', datetype='stop'):
        """
        return the latest date of locally existing POE/RES files

        Parameters
        ----------
        osvtype: {'POE', 'RES'}
            the type of orbit files required
        datetype: {'publish', 'start', 'stop'}
            one of three possible date types contained in the OSV filename

        Returns
        -------
        str
            a timestamp in format YYYYmmddTHHMMSS
        """
        address, directory = self._typeEvaluate(osvtype)
        files = finder(directory, [self.pattern], regex=True)
        return max([self.date(x, datetype) for x in files]) if len(files) > 0 else None

    def mindate(self, osvtype='POE', datetype='start'):
        """
        return the earliest date of locally existing POE/RES files

        Parameters
        ----------
        osvtype: {'POE', 'RES'}
            the type of orbit files required
        datetype: {'publish', 'start', 'stop'}
            one of three possible date types contained in the OSV filename

        Returns
        -------
        str
            a timestamp in format YYYYmmddTHHMMSS
        """
        address, directory = self._typeEvaluate(osvtype)
        files = finder(directory, [self.pattern], regex=True)
        return min([self.date(x, datetype) for x in files]) if len(files) > 0 else None

    def match(self, timestamp, osvtype='POE'):
        """
        return the corresponding OSV file for the provided time stamp.
        The file returned is one which covers the acquisition time and, if multiple exist,
        the one which was published last.
        In case a list of options is provided as osvtype, the file of higher accuracy (i.e. POE over RES) is returned.

        Parameters
        ----------
        timestamp: str
            the time stamp in the format 'YYYmmddTHHMMSS'
        osvtype: {'POE', 'RES'} or list
            the type of orbit files required; either 'POE', 'RES' or a list of both

        Returns
        -------
        str
            the best matching orbit file (overlapping time plus latest publication date)
        """
        # list all locally existing files of the defined type
        if osvtype in ['POE', 'RES']:
            locals = self.getLocals(osvtype)
            # filter the files to those which contain data for the defined time stamp
            files = [x for x in locals if self.date(x, 'start') <= timestamp <= self.date(x, 'stop')]
            if len(files) > 0:
                # select the file which was published last
                best = self.sortByDate(files, 'publish')[-1]
                return best
            elif len(files) == 1:
                return files[0]
            return None
        elif sorted(osvtype) == ['POE', 'RES']:
            best = self.match(timestamp, 'POE')
            if not best:
                best = self.match(timestamp, 'RES')
            return best

    def retrieve(self, files):
        """
        download a list of remote files into the respective subdirectories, i.e. POEORB or RESORB

        Parameters
        ----------
        files: list
            a list of remotely existing OSV files as returned by method :meth:`catch`

        Returns
        -------
        """
        self._init_dir()
        for type in ['POE', 'RES']:
            address, outdir = self._typeEvaluate(type)
            downloads = [x for x in files
                         if re.search('{}ORB'.format(type), x) and
                         not os.path.isfile(os.path.join(outdir, os.path.basename(x)))]
            for item in downloads:
                infile = urlopen(item, context=self.sslcontext)
                with open(os.path.join(outdir, os.path.basename(item)), 'wb') as outfile:
                    outfile.write(infile.read())
                infile.close()

    def sortByDate(self, files, datetype='start'):
        """
        sort a list of OSV files by a specific date type

        Parameters
        ----------
        files: list
            some OSV files
        datetype: {'publish', 'start', 'stop'}
            one of three possible date types contained in the OSV filename

        Returns
        -------
        list
            the input OSV files sorted by the defined date
        """
        return sorted(files, key=lambda x: self.date(x, datetype))

    def update(self, update_res=True):
        """
        Caution! This method is intended for downloading all available POE files and all RES files for whose
        time span no POE file yet exists.
        This will be a data volume of several GB and is particularly suited for multi-node SAR processing where not
        all nodes might have internet access and thus all files have to be downloaded before starting the processing.

        If you want to download the OSV file for a single scene either use the respective methods
        of the SAR drivers (e.g. :meth:`pyroSAR.drivers.SAFE.getOSV`) or methods :meth:`catch` and :meth:`retrieve` in combination.

        Perform creating/updating operations for POE and RES files:
        download newest POE and RES files, delete RES files which can be replaced by newly downloaded POE files.

        actions performed:
         * the ESA Quality Control (QC) server is checked for any POE files not in the local directory
         * POE  files on the server and not in the local directory are downloaded
         * RES files newer than the latest POE file are downloaded; POE files are approximately 18 days behind the actual date, thus RES files can be used instead
         * delete all RES files for whose date a POE file now exists locally


        Parameters
        ----------
        update_res: bool
            should the RES files also be updated (or just the POE files)

        Returns
        -------

        """
        self._init_dir()
        try:
            files_poe = self.catch('POE', start=self.maxdate('POE', 'start'))
        except RuntimeError as e:
            raise e
        self.retrieve(files_poe)
        if update_res:
            print('---------------------------------------------------------')
            files_res = self.catch('RES', start=self.maxdate('RES', 'start'))
            self.retrieve(files_res)
            self.clean_res()
