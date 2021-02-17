###############################################################################
# general utilities for Sentinel-1

# Copyright (c) 2016-2020, the pyroSAR Developers.

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
import zipfile as zf
from datetime import datetime
import xml.etree.ElementTree as ET
import numpy as np
from osgeo import gdal
from osgeo.gdalconst import GA_Update
from . import linesimplify as ls
from pyroSAR.examine import ExamineSnap
import progressbar as pb

from spatialist.ancillary import finder, urlQueryParser

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
    
    def __init__(self, osvdir=None):
        if osvdir is None:
            try:
                auxdatapath = ExamineSnap().auxdatapath
            except AttributeError:
                auxdatapath = os.path.join(os.path.expanduser('~'), '.snap', 'auxdata')
            osvdir = os.path.join(auxdatapath, 'Orbits', 'Sentinel-1')
        self.url = 'https://qc.sentinel1.eo.esa.int/api/v1/'
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
                        print('compressing and reorganizing EOF files')
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
    
    def catch(self, sensor, osvtype='POE', start=None, stop=None):
        """
        check a server for files

        Parameters
        ----------
        sensor: str or list
            The S1 mission(s):
             - 'S1A'
             - 'S1B'
             - ['S1A', 'S1B']
        osvtype: {'POE', 'RES'}
            the type of orbit files required
        start: str
            the date to start searching for files in format YYYYmmddTHHMMSS
        stop: str
            the date to stop searching for files in format YYYYmmddTHHMMSS

        Returns
        -------
        list
            the URLs of the remote OSV files
        """
        # a dictionary for storing the url arguments
        query = {}
        
        if osvtype == 'POE':
            query['product_type'] = 'AUX_POEORB'
        elif osvtype == 'RES':
            query['product_type'] = 'AUX_RESORB'
        else:
            raise RuntimeError("osvtype must be either 'POE' or 'RES'")
        
        if sensor in ['S1A', 'S1B']:
            query['sentinel1__mission'] = sensor
        elif sorted(sensor) == ['S1A', 'S1B']:
            pass
        else:
            raise RuntimeError('unsupported input for parameter sensor')
        
        # the collection of files to be returned
        collection = []
        # set the defined date or the date of the first existing OSV file otherwise
        # two days are added/subtracted from the defined start and stop dates since the
        # online query does only allow for searching the start time; hence, if e.g.
        # the start date is 2018-01-01T000000, the query would not return the corresponding
        # file, whose start date is 2017-12-31 (V20171231T225942_20180102T005942)
        if start is not None:
            date_start = datetime.strptime(start, '%Y%m%dT%H%M%S').strftime('%Y-%m-%dT%H:%M:%S')
        else:
            date_start = '2014-07-31'
        # set the defined date or the current date otherwise
        if stop is not None:
            date_stop = datetime.strptime(stop, '%Y%m%dT%H%M%S').strftime('%Y-%m-%dT%H:%M:%S')
        else:
            date_stop = datetime.now().strftime('%Y-%m-%dT%H:%M:%S')
        
        # append the time frame to the query dictionary
        query['validity_start__gte'] = date_start
        query['validity_stop__lte'] = date_stop
        print('searching for new {} files'.format(osvtype))
        target = urlQueryParser(self.url, query).replace('%3A', ':')
        print(target)
        while target is not None:
            response = requests.get(target).json()
            remotes = [item['remote_url'] for item in response['results']]
            collection += remotes
            target = response['next']
        if osvtype == 'RES' and self.maxdate('POE', 'stop') is not None:
            collection = [x for x in collection
                          if self.date(x, 'start') > self.maxdate('POE', 'stop')]
        return collection
    
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
            print('deleting {} RES file{}'.format(len(deprecated), '' if len(deprecated) == 1 else 's'))
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
        osvtype: str or list
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
    
    def retrieve(self, files, pbar=False):
        """
        download a list of remote files into the respective subdirectories, i.e. POEORB or RESORB

        Parameters
        ----------
        files: list
            a list of remotely existing OSV files as returned by method :meth:`catch`
        pbar: bool
            add a progressbar?

        Returns
        -------
        """
        downloads = []
        for remote in files:
            outdir = self._subdir(remote)
            os.makedirs(outdir, exist_ok=True)
            basename = os.path.basename(remote)
            local = os.path.join(outdir, basename) + '.zip'
            if not os.path.isfile(local):
                downloads.append((remote, local, basename))
        if len(downloads) == 0:
            return
        print('downloading {} file{}'.format(len(downloads), '' if len(downloads) == 1 else 's'))
        if pbar:
            progress = pb.ProgressBar(max_value=len(downloads))
        i = 0
        for remote, local, basename in downloads:
            infile = requests.get(remote)
            with zf.ZipFile(file=local,
                            mode='w',
                            compression=zf.ZIP_DEFLATED) \
                    as outfile:
                outfile.writestr(zinfo_or_arcname=basename,
                                 data=infile.content)
            infile.close()
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
    recommended by ESA and implemented in SNAP and additionally adds further refinement of the result using an image
    border line simplification approach. In this approach the border between valid and invalid pixels is first
    simplified using the method by Visvalingam and Whyatt references below. The line segments of the new border are then
    shifted until all pixels considered invalid before the simplification are again on one side of the line.
    See image below for further clarification.

    References:
        - 'Masking "No-value" Pixels on GRD Products generated by the Sentinel-1 ESA IPF'
          (issue 2.1 Jan 29 2018); available online under
          https://sentinel.esa.int/web/sentinel/user-guides/sentinel-1-sar/document-library
        - Visvalingam, M. and Whyatt J.D. (1993):
          "Line Generalisation by Repeated Elimination of Points",
          Cartographic J., 30 (1), 46 - 51

    Parameters
    ----------
    scene: ~pyroSAR.drivers.SAFE
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
        result using the Visvalingam-Whyatt method of poly-line vertex reduction. The green line is the final result
        after further correcting the VW-simplified result.

    """
    if scene.compression is not None:
        raise RuntimeError('scene is not yet unpacked')
    
    if method not in ['pyroSAR', 'ESA']:
        raise AttributeError("parameter 'method' must be either 'pyroSAR' or 'ESA'")
    
    blocksize = 2000
    
    # compute noise scaling factor
    if scene.meta['IPF_version'] >= 2.9:
        print('border noise removal not necessary for IPF version {}'.format(scene.meta['IPF_version']))
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
        print(subset)
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
