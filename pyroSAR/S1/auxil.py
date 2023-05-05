###############################################################################
# general utilities for Sentinel-1

# Copyright (c) 2016-2021, the pyroSAR Developers.

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
import sys
import requests
import tempfile
import zipfile as zf
from io import BytesIO
from datetime import datetime, timedelta
from dateutil import parser as dateutil_parser
import xml.etree.ElementTree as ET
import numpy as np
from osgeo import gdal
from osgeo.gdalconst import GA_Update
from . import linesimplify as ls
from pyroSAR.examine import ExamineSnap
import progressbar as pb

from spatialist.ancillary import finder

import logging
log = logging.getLogger(__name__)

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
    timeout: int or tuple or None
        the timeout in seconds for downloading OSV files as provided to :func:`requests.get`
    
    See Also
    --------
    `requests timeouts <https://requests.readthedocs.io/en/master/user/advanced/#timeouts>`_
    """
    
    def __init__(self, osvdir=None, timeout=300):
        self.timeout = timeout
        if osvdir is None:
            try:
                auxdatapath = ExamineSnap().auxdatapath
            except AttributeError:
                auxdatapath = os.path.join(os.path.expanduser('~'), '.snap', 'auxdata')
            osvdir = os.path.join(auxdatapath, 'Orbits', 'Sentinel-1')
        self.outdir_poe = os.path.join(osvdir, 'POEORB')
        self.outdir_res = os.path.join(osvdir, 'RESORB')
        self.pattern = r'S1[AB]_OPER_AUX_(?:POE|RES)ORB_OPOD_[0-9TV_]{48}\.EOF'
        self.pattern_fine = r'(?P<sensor>S1[AB])_OPER_AUX_' \
                            r'(?P<type>(?:POE|RES)ORB)_OPOD_' \
                            r'(?P<publish>[0-9]{8}T[0-9]{6})_V' \
                            r'(?P<start>[0-9]{8}T[0-9]{6})_' \
                            r'(?P<stop>[0-9]{8}T[0-9]{6})\.EOF'
        self._init_dir()
        self._reorganize()
    
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
    
    def _parse(self, file):
        basename = os.path.basename(file)
        groups = re.match(self.pattern_fine, basename).groupdict()
        return groups
    
    def _reorganize(self):
        """
        compress and move EOF files into subdirectories

        Returns
        -------

        """
        message = True
        for subdir in [self.outdir_poe, self.outdir_res]:
            if not os.path.isdir(subdir):
                continue
            files = finder(subdir, [self.pattern], recursive=False, regex=True)
            for eof in files:
                base = os.path.basename(eof)
                target = os.path.join(self._subdir(eof), base + '.zip')
                os.makedirs(os.path.dirname(target), exist_ok=True)
                if not os.path.isfile(target):
                    if message:
                        log.info('compressing and reorganizing EOF files')
                        message = False
                    with zf.ZipFile(file=target,
                                    mode='w',
                                    compression=zf.ZIP_DEFLATED) as zip:
                        zip.write(filename=eof,
                                  arcname=base)
                os.remove(eof)
    
    def _typeEvaluate(self, osvtype):
        """
        evaluate the 'osvtype' method argument and return the corresponding remote repository and local directory

        Parameters
        ----------
        osvtype: str
            the type of orbit files required; either 'POE' or 'RES'

        Returns
        -------
        tuple of str
            the remote repository and local directory of the osv type
        """
        if osvtype not in ['POE', 'RES']:
            raise IOError('type must be either "POE" or "RES"')
        if osvtype == 'POE':
            return self.outdir_poe
        else:
            return self.outdir_res
    
    def __catch_aux_sentinel(self, sensor, start, stop, osvtype='POE'):
        url = 'http://aux.sentinel1.eo.esa.int'
        skeleton = '{url}/{osvtype}ORB/{year}/{month:02d}/{day:02d}/'
        
        files = []
        date_search = start
        busy = True
        while busy:
            url_sub = skeleton.format(url=url,
                                      osvtype=osvtype,
                                      year=date_search.year,
                                      month=date_search.month,
                                      day=date_search.day)
            response = requests.get(url_sub, timeout=self.timeout)
            response.raise_for_status()
            result = response.text
            files_sub = list(set(re.findall(self.pattern, result)))
            if len(files_sub) == 0:
                break
            for file in files_sub:
                match = re.match(self.pattern_fine, file)
                start2 = datetime.strptime(match.group('start'), '%Y%m%dT%H%M%S')
                stop2 = datetime.strptime(match.group('stop'), '%Y%m%dT%H%M%S')
                if sensor == match.group('sensor'):
                    if start2 < stop and stop2 > start:
                        log.info(url_sub)
                        files.append({'filename': file,
                                      'href': url_sub + '/' + file,
                                      'auth': None})
                if start2 >= stop:
                    busy = False
            date_search += timedelta(days=1)
        
        return files
    
    def __catch_step_auxdata(self, sensor, start, stop, osvtype='POE'):
        url = 'https://step.esa.int/auxdata/orbits/Sentinel-1'
        skeleton = '{url}/{osvtype}ORB/{sensor}/{year}/{month:02d}/'
        
        if isinstance(sensor, str):
            sensor = [sensor]
        
        files = []
        for sens in sensor:
            date_search = datetime(year=start.year,
                                   month=start.month,
                                   day=1)
            busy = True
            while busy:
                url_sub = skeleton.format(url=url,
                                          osvtype=osvtype,
                                          sensor=sens,
                                          year=date_search.year,
                                          month=date_search.month)
                log.info(url_sub)
                response = requests.get(url_sub, timeout=self.timeout)
                response.raise_for_status()
                result = response.text
                files_sub = list(set(re.findall(self.pattern, result)))
                if len(files_sub) == 0:
                    break
                for file in files_sub:
                    match = re.match(self.pattern_fine, file)
                    start2 = datetime.strptime(match.group('start'), '%Y%m%dT%H%M%S')
                    stop2 = datetime.strptime(match.group('stop'), '%Y%m%dT%H%M%S')
                    if start2 < stop and stop2 > start:
                        files.append({'filename': file,
                                      'href': url_sub + '/' + file + '.zip',
                                      'auth': None})
                    if start2 >= stop:
                        busy = False
                if date_search.month < 12:
                    date_search = date_search.replace(month=date_search.month + 1)
                else:
                    date_search = date_search.replace(year=date_search.year + 1, month=1)
        
        return files
    
    def __catch_gnss(self, sensor, start, stop, osvtype='POE'):
        url = 'https://scihub.copernicus.eu/gnss'
        redirect = 'https://dhusfeed.dhus.onda-dias.net/gnss'
        auth = ('gnssguest', 'gnssguest')
        # a dictionary for storing the url arguments
        query = {}
        
        if osvtype == 'POE':
            query['producttype'] = 'AUX_POEORB'
        elif osvtype == 'RES':
            query['producttype'] = 'AUX_RESORB'
        else:
            raise RuntimeError("osvtype must be either 'POE' or 'RES'")
        
        if sensor in ['S1A', 'S1B']:
            query['platformname'] = 'Sentinel-1'
            # filename starts w/ sensor
            query['filename'] = '{}*'.format(sensor)
        elif sorted(sensor) == ['S1A', 'S1B']:
            query['platformname'] = 'Sentinel-1'
        else:
            raise RuntimeError('unsupported input for parameter sensor')
        
        # the collection of files to be returned
        collection = []
        
        date_start = start.strftime('%Y-%m-%dT%H:%M:%SZ')
        date_stop = stop.strftime('%Y-%m-%dT%H:%M:%SZ')
        
        # append the time frame to the query dictionary
        query['beginPosition'] = '[{} TO {}]'.format(date_start, date_stop)
        query['endPosition'] = '[{} TO {}]'.format(date_start, date_stop)
        query_list = []
        for keyword, value in query.items():
            query_elem = '{}:{}'.format(keyword, value)
            query_list.append(query_elem)
        query_str = ' '.join(query_list)
        target = '{}/search?q={}&format=json'.format(url, query_str)
        log.info(target)
        
        def _parse_gnsssearch_json(search_dict):
            parsed_dict = {}
            # Will return ['entry'] as dict if only one item
            # If so just make a list
            if isinstance(search_dict, dict):
                search_dict = [search_dict]
            for entry in search_dict:
                id = entry['id']
                entry_dict = {}
                
                for key, value in entry.items():
                    if key == 'title':
                        entry_dict[key] = value
                    elif key == 'id':
                        entry_dict[key] = value
                    elif key == 'ondemand':
                        if value.lower() == 'true':
                            entry_dict[key] = True
                        else:
                            entry_dict[key] = False
                    elif key == 'str':
                        for elem in value:
                            entry_dict[elem['name']] = elem['content']
                    elif key == 'link':
                        for elem in value:
                            if 'rel' in elem.keys():
                                href_key = 'href_' + elem['rel']
                                entry_dict[href_key] = elem['href']
                            else:
                                entry_dict['href'] = elem['href']
                    elif key == 'date':
                        for elem in value:
                            entry_dict[elem['name']] = dateutil_parser.parse(elem['content'])
                
                parsed_dict[id] = entry_dict
            return parsed_dict
        
        def _parse_gnsssearch_response(response_json):
            if 'entry' in response_json.keys():
                search_dict = response_json['entry']
                parsed_dict = _parse_gnsssearch_json(search_dict)
            else:
                parsed_dict = {}
            return parsed_dict
        
        response = requests.get(target, auth=auth, timeout=self.timeout)
        response.raise_for_status()
        response_json = response.json()['feed']
        total_results = response_json['opensearch:totalResults']
        subquery = [link['href'] for link in response_json['link'] if link['rel'] == 'self'][0]
        subquery = subquery.replace(redirect, url.strip())
        if int(total_results) > 10:
            subquery = subquery.replace('rows=10', 'rows=100')
        while subquery:
            subquery_response = requests.get(subquery, auth=auth, timeout=self.timeout)
            subquery_response.raise_for_status()
            subquery_json = subquery_response.json()['feed']
            subquery_products = _parse_gnsssearch_response(subquery_json)
            items = list(subquery_products.values())
            for item in items:
                item['auth'] = auth
            collection += list(subquery_products.values())
            if 'next' in [link['rel'] for link in subquery_json['link']]:
                subquery = [link['href'] for link in subquery_json['link'] if link['rel'] == 'next'][0]
                subquery = subquery.replace(redirect, url.strip())
            else:
                subquery = None
        if osvtype == 'RES' and self.maxdate('POE', 'stop') is not None:
            collection = [x for x in collection
                          if self.date(x['filename'], 'start') > self.maxdate('POE', 'stop')]
        for item in collection:
            item['href'] = item['href'].replace(redirect, url)
        return collection
    
    def catch(self, sensor, osvtype='POE', start=None, stop=None, url_option=1):
        """
        check a server for files

        Parameters
        ----------
        sensor: str or list[str]
            The S1 mission(s):
            
             - 'S1A'
             - 'S1B'
             - ['S1A', 'S1B']
        osvtype: str or list[str]
            the type of orbit files required
        start: str or None
            the date to start searching for files in format YYYYmmddTHHMMSS
        stop: str or None
            the date to stop searching for files in format YYYYmmddTHHMMSS
        url_option: int
            the URL to query for OSV files
            
             - 1: https://scihub.copernicus.eu/gnss
             - 2: https://step.esa.int/auxdata/orbits/Sentinel-1

        Returns
        -------
        list[dict]
            the product dictionary of the remote OSV files, with href
        """
        
        log.info('searching for new {} files'.format(osvtype))
        
        if start is not None:
            start = datetime.strptime(start, '%Y%m%dT%H%M%S')
        else:
            start = datetime.strptime('2014-07-31', '%Y-%m-%d')
        # set the defined date or the current date otherwise
        if stop is not None:
            stop = datetime.strptime(stop, '%Y%m%dT%H%M%S')
        else:
            stop = datetime.now()
        
        if url_option == 1:
            items = self.__catch_gnss(sensor, start, stop, osvtype)
        elif url_option == 2:
            items = self.__catch_step_auxdata(sensor, start, stop, osvtype)
        else:
            raise ValueError("'url_option' must be either 1 or 2")
        
        if osvtype == 'RES' and self.maxdate('POE', 'stop') is not None:
            items = [x for x in items
                     if self.date(x['filename'], 'start') > self.maxdate('POE', 'stop')]
        log.info('found {} results'.format(len(items)))
        
        return items
    
    def date(self, file, datetype):
        """
        extract a date from an OSV file name

        Parameters
        ----------
        file: str
            the OSV file
        datetype: {'publish', 'start', 'stop'}
            one of three possible date types contained in the OSV filename

        Returns
        -------
        str
            a time stamp in the format YYYYmmddTHHMMSS
        """
        return self._parse(file)[datetype]
    
    def clean_res(self):
        """
        delete all RES files for whose date a POE file exists
        """
        maxdate_poe = self.maxdate('POE', 'stop')
        if maxdate_poe is not None:
            deprecated = [x for x in self.getLocals('RES') if self.date(x, 'stop') < maxdate_poe]
            log.info('deleting {} RES file{}'.format(len(deprecated), '' if len(deprecated) == 1 else 's'))
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
        list[str]
            a selection of local OSV files
        """
        directory = self._typeEvaluate(osvtype)
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
        directory = self._typeEvaluate(osvtype)
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
        directory = self._typeEvaluate(osvtype)
        files = finder(directory, [self.pattern], regex=True)
        return min([self.date(x, datetype) for x in files]) if len(files) > 0 else None
    
    def match(self, sensor, timestamp, osvtype='POE'):
        """
        return the corresponding OSV file for the provided sensor and time stamp.
        The file returned is one which covers the acquisition time and, if multiple exist,
        the one which was published last.
        In case a list of options is provided as osvtype, the file of higher accuracy (i.e. POE over RES) is returned.

        Parameters
        ----------
        sensor: str
            The S1 mission:
            
             - 'S1A'
             - 'S1B'
        timestamp: str
            the time stamp in the format 'YYYmmddTHHMMSS'
        osvtype: str or list[str]
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
            files = [x for x in files if os.path.basename(x).startswith(sensor)]
            if len(files) > 0:
                # select the file which was published last
                best = self.sortByDate(files, 'publish')[-1]
                return best
            elif len(files) == 1:
                return files[0]
            return None
        elif sorted(osvtype) == ['POE', 'RES']:
            best = self.match(sensor=sensor, timestamp=timestamp, osvtype='POE')
            if not best:
                best = self.match(sensor=sensor, timestamp=timestamp, osvtype='RES')
            return best
    
    def retrieve(self, products, pbar=False):
        """
        download a list of product dictionaries into the respective subdirectories, i.e. POEORB or RESORB

        Parameters
        ----------
        products: list[dict]
            a list of remotely existing OSV product dictionaries as returned by method :meth:`catch`
        pbar: bool
            add a progressbar?

        Returns
        -------
        """
        downloads = []
        for product in products:
            if all(key not in ['filename', 'href'] for key in product.keys()):
                raise RuntimeError("product dictionaries must contain 'filename' and 'href' keys")
            basename = product['filename']
            remote = product['href']
            auth = product['auth']
            
            outdir = self._subdir(basename)
            os.makedirs(outdir, exist_ok=True)
            local = os.path.join(outdir, basename) + '.zip'
            if not os.path.isfile(local):
                downloads.append((remote, local, basename, auth))
        if len(downloads) == 0:
            return
        log.info('downloading {} file{}'.format(len(downloads), '' if len(downloads) == 1 else 's'))
        if pbar:
            progress = pb.ProgressBar(max_value=len(downloads))
        else:
            progress = None
        i = 0
        for remote, local, basename, auth in downloads:
            response = requests.get(remote, auth=auth, timeout=self.timeout)
            response.raise_for_status()
            infile = response.content

            # use a tempfile to allow atomic writes in the case of
            # parallel executions dependent on the same orbit files
            fd, tmp_path = tempfile.mkstemp(prefix=os.path.basename(local), dir=os.path.dirname(local))
            os.close(fd)
            try:
                if remote.endswith('.zip'):
                    with zf.ZipFile(file=BytesIO(infile)) as tmp:
                        members = tmp.namelist()
                        target = [x for x in members if re.search(basename, x)][0]
                        with zf.ZipFile(tmp_path, 'w') as outfile:
                            outfile.write(filename=tmp.extract(target),
                                          arcname=basename)
                else:
                    with zf.ZipFile(file=tmp_path,
                                    mode='w',
                                    compression=zf.ZIP_DEFLATED) \
                            as outfile:
                        outfile.writestr(zinfo_or_arcname=basename,
                                         data=infile)
                os.rename(tmp_path, local)
            except Exception as e:
                os.unlink(tmp_path)
                raise

            if pbar:
                i += 1
                progress.update(i)
        if pbar:
            progress.finish()
        self.clean_res()
    
    def sortByDate(self, files, datetype='start'):
        """
        sort a list of OSV files by a specific date type

        Parameters
        ----------
        files: list[str]
            some OSV files
        datetype: {'publish', 'start', 'stop'}
            one of three possible date types contained in the OSV filename

        Returns
        -------
        list[str]
            the input OSV files sorted by the defined date
        """
        return sorted(files, key=lambda x: self.date(x, datetype))
    
    def _subdir(self, file):
        """
        | return the subdirectory in which to store the EOF file,
        | i.e. basedir/{type}ORB/{sensor}/{year}/{month}
        | e.g. basedir/POEORB/S1A/2018/12

        Parameters
        ----------
        file: str
            the EOF filename

        Returns
        -------
        str
            the target directory
        """
        attr = self._parse(file)
        outdir = self._typeEvaluate(attr['type'][:3])
        start = self.date(file, datetype='start')
        start = datetime.strptime(start, '%Y%m%dT%H%M%S')
        month = '{:02d}'.format(start.month)
        outdir = os.path.join(outdir, attr['sensor'],
                              str(start.year), month)
        return outdir


def removeGRDBorderNoise(scene, method='pyroSAR'):
    """
    Mask out Sentinel-1 image border noise. This function implements the method for removing GRD border noise as
    published by ESA :cite:`Miranda2018` and implemented in SNAP and additionally adds further refinement of the result using an image
    border line simplification approach. In this approach the border between valid and invalid pixels is first
    simplified using the poly-line vertex reduction method by Visvalingam and Whyatt :cite:`Visvalingam1993`.
    The line segments of the new border are then shifted until all pixels considered invalid before the simplification
    are again on one side of the line. See image below for further clarification.

    Parameters
    ----------
    scene: pyroSAR.drivers.SAFE
        the Sentinel-1 scene object
    method: str
        the border noise removal method to be applied; one of the following:
        
         - 'ESA': the pure implementation as described by ESA
         - 'pyroSAR': the ESA method plus the custom pyroSAR refinement


    .. figure:: figures/S1_bnr.png
        :scale: 30%

        Demonstration of the border noise removal for a vertical left image border. The area under the respective lines
        covers pixels considered valid, everything above will be masked out. The blue line is the result of the noise
        removal as recommended by ESA, in which a lot of noise is still present. The red line is the over-simplified
        result using the Visvalingam-Whyatt method. The green line is the final result after further correcting the
        VW-simplified result.

    """
    if scene.compression is not None:
        raise RuntimeError('scene is not yet unpacked')
    
    if method not in ['pyroSAR', 'ESA']:
        raise AttributeError("parameter 'method' must be either 'pyroSAR' or 'ESA'")
    
    blocksize = 2000
    
    # compute noise scaling factor
    if scene.meta['IPF_version'] >= 2.9:
        log.info('border noise removal not necessary for IPF version {}'.format(scene.meta['IPF_version']))
        return
    elif scene.meta['IPF_version'] <= 2.5:
        knoise = {'IW': 75088.7, 'EW': 56065.87}[scene.acquisition_mode]
        cads = scene.getFileObj(scene.findfiles('calibration-s1[ab]-[ie]w-grd-(?:hh|vv)')[0])
        caltree = ET.fromstring(cads.read())
        cads.close()
        adn = float(caltree.find('.//calibrationVector/dn').text.split()[0])
        if scene.meta['IPF_version'] < 2.34:
            scalingFactor = knoise * adn
        else:
            scalingFactor = knoise * adn * adn
    else:
        scalingFactor = 1
    
    # read noise vectors from corresponding annotation xml
    noisefile = scene.getFileObj(scene.findfiles('noise-s1[ab]-[ie]w-grd-(?:hh|vv)')[0])
    noisetree = ET.fromstring(noisefile.read())
    noisefile.close()
    noiseVectors = noisetree.findall('.//noiseVector')
    
    # define boundaries of image subsets to be masked (4x the first lines/samples of the image boundaries)
    subsets = [(0, 0, blocksize, scene.lines),
               (0, 0, scene.samples, blocksize),
               (scene.samples - blocksize, 0, scene.samples, scene.lines),
               (0, scene.lines - blocksize, scene.samples, scene.lines)]
    
    # extract column indices of noise vectors
    yi = np.array([int(x.find('line').text) for x in noiseVectors])
    
    # create links to the tif files for a master co-polarization and all other polarizations as slaves
    master = scene.findfiles('s1.*(?:vv|hh).*tiff')[0]
    ras_master = gdal.Open(master, GA_Update)
    ras_slaves = [gdal.Open(x, GA_Update) for x in scene.findfiles('s1.*tiff') if x != master]
    
    outband_master = ras_master.GetRasterBand(1)
    outband_slaves = [x.GetRasterBand(1) for x in ras_slaves]
    
    # iterate over the four image subsets
    for subset in subsets:
        log.info(subset)
        xmin, ymin, xmax, ymax = subset
        xdiff = xmax - xmin
        ydiff = ymax - ymin
        # linear interpolation of noise vectors to array
        noise_interp = np.empty((ydiff, xdiff), dtype=float)
        for i in range(0, len(noiseVectors)):
            if ymin <= yi[i] <= ymax:
                # extract row indices of noise vector
                xi = [int(x) for x in noiseVectors[i].find('pixel').text.split()]
                # extract noise values
                noise = [float(x) for x in noiseVectors[i].find('noiseLut').text.split()]
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
        
        if method == 'pyroSAR':
            # helper functions for masking out negative values
            def helper1(x):
                return len(x) - np.argmax(x > 0)
            
            def helper2(x):
                return len(x) - np.argmax(x[::-1] > 0)
            
            # mask out negative values and simplify borders (custom implementation)
            if subset == (0, 0, blocksize, scene.lines):
                border = np.apply_along_axis(helper1, 1, denoisedBlock)
                border = blocksize - ls.reduce(border)
                for j in range(0, ydiff):
                    denoisedBlock[j, :border[j]] = 0
                    denoisedBlock[j, border[j]:] = 1
            elif subset == (0, scene.lines - blocksize, scene.samples, scene.lines):
                border = np.apply_along_axis(helper2, 0, denoisedBlock)
                border = ls.reduce(border)
                for j in range(0, xdiff):
                    denoisedBlock[border[j]:, j] = 0
                    denoisedBlock[:border[j], j] = 1
            elif subset == (scene.samples - blocksize, 0, scene.samples, scene.lines):
                border = np.apply_along_axis(helper2, 1, denoisedBlock)
                border = ls.reduce(border)
                for j in range(0, ydiff):
                    denoisedBlock[j, border[j]:] = 0
                    denoisedBlock[j, :border[j]] = 1
            elif subset == (0, 0, scene.samples, blocksize):
                border = np.apply_along_axis(helper1, 0, denoisedBlock)
                border = blocksize - ls.reduce(border)
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
