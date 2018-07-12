import sys
import warnings

if sys.version_info >= (3, 0):
    from io import StringIO
    from urllib.request import urlopen
    from urllib.error import HTTPError
else:
    from cStringIO import StringIO
    from urllib2 import urlopen, HTTPError

import os
import re
import shutil
import subprocess as sp
import xml.etree.ElementTree as ET
import zipfile as zf
from ftplib import FTP
from time import strftime, gmtime
from xml.dom import minidom
from os.path import expanduser

import pyroSAR
from pyroSAR import identify
from pyroSAR.ancillary import dissolve, finder, which
from pyroSAR.spatial import gdal_translate
from .._dev_config import LOOKUP


def parse_recipe(name):
    name = name if name.endswith('.xml') else name + '.xml'
    absname = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'recipes', name)
    with open(absname, 'r') as workflow:
        tree = ET.fromstring(workflow.read())
    return tree


def parse_node(name):
    name = name if name.endswith('.xml') else name + '.xml'
    absname = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'recipes', 'nodes', name)
    with open(absname, 'r') as workflow:
        tree = ET.fromstring(workflow.read())
    return tree


def parse_suffix(workflow):
    nodes = workflow.findall('node')
    suffix = '_'.join(filter(None, [LOOKUP.snap.suffix[x] for x in [y.attrib['id'] for y in nodes]]))
    return suffix


def insert_node(workflow, predecessor_id, node):
    predecessor = workflow.find('.//node[@id="{}"]'.format(predecessor_id))
    position = list(workflow).index(predecessor) + 1
    workflow.insert(position, node)
    newnode = workflow[position]
    newnode.find('.//sources/sourceProduct').attrib['refid'] = predecessor.attrib['id']
    workflow[position + 1].find('.//sources/sourceProduct').attrib['refid'] = newnode.attrib['id']


def write_recipe(recipe, outfile):
    outfile = outfile if outfile.endswith('.xml') else outfile + '.xml'
    rough_string = ET.tostring(recipe, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    with open(outfile, 'w') as out:
        out.write(reparsed.toprettyxml(indent='\t', newl=''))


def getOrbitContentVersions(contentVersion):
    return dict(
        [re.split('\s*=\s*', x.strip('\r')) for x in contentVersion.read().split('\n') if re.search('^[0-9]{4}', x)])


class GetAuxdata:
    def __init__(self, datasets, scenes):
        self.datasets = datasets
        self.scenes = [identify(scene) if isinstance(scene, str) else scene for scene in scenes]
        self.sensors = list(set([scene.sensor for scene in scenes]))

        try:
            self.auxDataPath = os.path.join(os.environ['HOME'], '.snap/auxdata')
        except KeyError:
            self.auxDataPath = os.path.join(os.environ['USERPROFILE'], '.snap/auxdata')

    def srtm_1sec_hgt(self):
        pass
        # Wird nicht benutzt?


#        for dataset in self.datasets:
#            files = [x.replace('hgt', 'SRTMGL1.hgt.zip') for x in
#                     list(set(dissolve([scene.getHGT() for scene in self.scenes])))]

def getAuxdata(datasets, scenes):
    auxDataPath = os.path.join(expanduser("~"), '.snap/auxdata')

    scenes = [identify(scene) if isinstance(scene, str) else scene for scene in scenes]
    sensors = list(set([scene.sensor for scene in scenes]))
    for dataset in datasets:
        if dataset == 'SRTM 1Sec HGT':
            files = [x.replace('hgt', 'SRTMGL1.hgt.zip') for x in
                     list(set(dissolve([scene.getHGT() for scene in scenes])))]
            for file in files:
                infile = os.path.join('http://step.esa.int/auxdata/dem/SRTMGL1', file)
                outfile = os.path.join(auxDataPath, 'dem/SRTM 1Sec HGT', file)
                if not os.path.isfile(outfile):
                    print(infile)
                    try:
                        input = urlopen(infile)
                    except HTTPError:
                        print('-> not available')
                        continue
                    with open(outfile, 'wb') as output:
                        output.write(input.read())
                    input.close()
        elif dataset == 'POEORB':
            for sensor in sensors:
                if re.search('S1[AB]', sensor):

                    dates = [(scene.start[:4], scene.start[4:6]) for scene in scenes]
                    years = list(set([x[0] for x in dates]))

                    remote_contentVersion = urlopen(
                        'http://step.esa.int/auxdata/orbits/Sentinel-1/POEORB/remote_contentVersion.txt')
                    versions_remote = getOrbitContentVersions(remote_contentVersion)

                    for year in years:
                        dir_orb = os.path.join(auxDataPath, 'Orbits/Sentinel-1/POEORB', year)

                        if not os.path.isdir(dir_orb):
                            os.makedirs(dir_orb)
                        contentVersionFile = os.path.join(dir_orb, 'contentVersion.txt')

                        if os.path.isfile(contentVersionFile):
                            contentVersion = open(contentVersionFile, 'r+')
                            versions_local = getOrbitContentVersions(contentVersion)
                        else:
                            contentVersion = open(contentVersionFile, 'w')
                            versions_local = {}

                        combine = dict(set(versions_local.items()) & set(versions_remote.items()))

                        dates_select = [x for x in dates if x[0] == year]
                        months = list(set([x[1] for x in dates_select]))

                        orb_ids = sorted(
                            [x for x in ['{}-{}.zip'.format(year, month) for month in months] if not x in combine])

                        if len(orb_ids) > 0:
                            contentVersion.write('#\n#{}\n'.format(strftime('%a %b %d %H:%M:%S %Z %Y', gmtime())))

                            for orb_id in orb_ids:
                                orb_remote = urlopen(
                                    'http://step.esa.int/auxdata/orbits/Sentinel-1/POEORB/{}'.format(orb_id))
                                orb_remote_stream = zf.ZipFile(StringIO(orb_remote.read()), 'r')
                                orb_remote.close()

                                targets = [x for x in orb_remote_stream.namelist() if
                                           not os.path.isfile(os.path.join(dir_orb, x))]
                                orb_remote_stream.extractall(dir_orb, targets)
                                orb_remote_stream.close()

                                versions_local[orb_id] = versions_remote[orb_id]

                                for key, val in versions_local.iteritems():
                                    contentVersion.write('{}={}\n'.format(key, val))

                        contentVersion.close()
                    remote_contentVersion.close()
                else:
                    print('not implemented yet')
        elif dataset == 'Delft Precise Orbits':
            path_server = 'dutlru2.lr.tudelft.nl'
            subdirs = {'ASAR:': 'ODR.ENVISAT1/eigen-cg03c', 'ERS1': 'ODR.ERS-1/dgm-e04', 'ERS2': 'ODR.ERS-2/dgm-e04'}
            ftp = FTP(path_server)
            ftp.login()
            for sensor in sensors:
                if sensor in subdirs.keys():
                    path_target = os.path.join('pub/orbits', subdirs[sensor])
                    path_local = os.path.join(auxDataPath, 'Orbits/Delft Precise Orbits', subdirs[sensor])
                    ftp.cwd(path_target)
                    for item in ftp.nlst():
                        ftp.retrbinary('RETR ' + item, open(os.path.join(path_local, item), 'wb').write)
            ftp.quit()
        else:
            print('not implemented yet')


def gpt(xmlfile):
    """
    wrapper for ESA SNAP Graph Processing Tool GPT
    input is a readily formatted workflow xml file as created by function geocode in module snap.util
    """
    try:
        snap_exec = ExamineSnap().path
    except AttributeError:
        raise RuntimeError('could not find SNAP executable')
    gpt_exec = os.path.join(os.path.dirname(snap_exec), 'gpt')

    with open(xmlfile, 'r') as infile:
        workflow = ET.fromstring(infile.read())
    write = workflow.find('.//node[@id="Write"]')
    outname = write.find('.//parameters/file').text
    outdir = os.path.dirname(outname)
    format = write.find('.//parameters/formatName').text
    infile = workflow.find('.//node[@id="Read"]/parameters/file').text

    if format == 'GeoTiff-BigTIFF':
        cmd = [gpt_exec,
               # '-Dsnap.dataio.reader.tileWidth=*',
               # '-Dsnap.dataio.reader.tileHeight=1',
               '-Dsnap.dataio.bigtiff.tiling.width=256',
               '-Dsnap.dataio.bigtiff.tiling.height=256',
               # '-Dsnap.dataio.bigtiff.compression.type=LZW',
               # '-Dsnap.dataio.bigtiff.compression.quality=0.75',
               xmlfile]
    else:
        cmd = [gpt_exec, xmlfile]

    proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = proc.communicate()
    if proc.returncode != 0:
        if os.path.isfile(outname + '.tif'):
            os.remove(outname + '.tif')
        elif os.path.isdir(outname):
            shutil.rmtree(outname)
        print(out + err)
        print('failed: {}'.format(os.path.basename(infile)))
        err_match = re.search('Error: (.*)\n', out + err)
        errmessage = err_match.group(1) if err_match else err
        raise RuntimeError(errmessage)

    if format == 'ENVI':
        id = pyroSAR.identify(infile)
        suffix = parse_suffix(workflow)
        for item in finder(outname, ['*.img']):
            pol = re.search('[HV]{2}', item).group()
            name_new = os.path.join(outdir, '{}_{}_{}.tif'.format(id.outname_base(), pol, suffix))
            translateoptions = {'options': ['-q', '-co', 'INTERLEAVE=BAND', '-co', 'TILED=YES'], 'format': 'GTiff'}
            gdal_translate(item, name_new, translateoptions)
        shutil.rmtree(outname)


class ExamineSnap(object):
    """
    Class to check if snap is installed. This will be called with snap.__init__
    as snap_config. If you are running multiple snap versions or the package can
    not find the snap executable, you can set an path via: snap_config.set_path("path")
    """

    def __init__(self):

        candidates = ['snap64.exe', 'snap32.exe', 'snap.exe', 'snap']
        executables = list(filter(None, [which(x) for x in candidates]))

        if len(executables) == 0:

            self.path = None

            warnings.warn(
                'Cannot find SNAP installation. You can download it from '
                'http://step.esa.int/main/download/ or you can specify a path with '
                'snap_config.set_path(path_to_snap)', UserWarning)
        else:
            self.path = executables[0]

            if os.path.islink(self.path):
                self.path = os.path.realpath(self.path)

            if len(executables) > 1:
                raise warnings.warn(
                    'There are more than one instances of SNAP installed:\n{options}\n'
                    'Using the first option:\n{select}\n'
                    'You can define a different one with ExamineSnap.set_path()'
                    ''.format(options='\n'.join(executables), select=self.path), UserWarning)

            self.__get_etc()
            self.__read_config()

        self.auxdatapath = os.path.join(expanduser('~'), '.snap/auxdata')

    def __get_etc(self):
        try:
            self.etc = os.path.join(os.path.dirname(os.path.dirname(self.path)), 'etc')
            self.auxdata = os.listdir(self.etc)
            self.config_path = os.path.join(self.etc, [s for s in self.auxdata if 'snap.auxdata.properties' in s][0])

        except OSError:
            raise AssertionError('ETC directory is not existent.')

    def __read_config(self):
        with open(self.config_path) as config:
            self.config = []
            for line in config:
                self.config.append(line)
