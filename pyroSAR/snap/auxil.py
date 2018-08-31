import sys
import warnings
import pkg_resources

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
from pyroSAR import identify, ConfigHandler
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


def gpt(xmlfile, basename_extensions=None):
    """
    wrapper for ESA SNAP Graph Processing Tool GPT
    input is a readily formatted workflow xml file as created by function geocode in module snap.util
    """
    try:
        gpt_exec = ExamineSnap().gpt
    except AttributeError:
        raise RuntimeError('could not find SNAP GPT executable')

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
            name_new = os.path.join(outdir, '{}_{}_{}.tif'.format(id.outname_base(basename_extensions), pol, suffix))
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
        
        self.identifiers = ['path', 'gpt', 'etc', 'auxdata']
        
        self.__read_config()
        
        if not self.__is_identified():
            self.__identify_snap()

        self.__read_config_attr('auxdatapath', 'SNAP')
        if not hasattr(self, 'auxdatapath'):
            self.auxdatapath = os.path.join(os.path.expanduser('~'), '.snap', 'auxdata')

        self.__read_config_attr('properties', 'SNAP')
        if not hasattr(self, 'properties'):
            path = 'data/snap.auxdata.properties'
            self.properties = pkg_resources.resource_filename(__name__, path)

        if not hasattr(self, 'snap_properties'):
            self.__read_snap_properties()
        
        self.__update_config()
    
    def __is_identified(self):
        return len(list(filter(None, self.identifiers))) == 0
    
    def __identify_snap(self):
        defaults = ['snap64.exe', 'snap32.exe', 'snap.exe', 'snap']
        paths = os.environ['PATH'].split(os.path.pathsep)
        options = [os.path.join(path, option) for path in paths for option in defaults]
        options = [x for x in options if os.path.isfile(x)]
        
        if not hasattr(self, 'path') or not os.path.isfile(self.path):
            executables = options
        else:
            executables = [self.path] + options
        for path in executables:
            if os.path.islink(path):
                path = os.path.realpath(path)
            etc = os.path.join(os.path.dirname(os.path.dirname(path)), 'etc')
            if not os.path.isdir(etc):
                continue
            auxdata = os.listdir(etc)
            if 'snap.auxdata.properties' not in auxdata:
                continue
            else:
                auxdata_properties = os.path.join(etc, 'snap.auxdata.properties')
            gpt_candidates = finder(os.path.dirname(path), ['gpt', 'gpt.exe'])
            if len(gpt_candidates) == 0:
                continue
            else:
                gpt = gpt_candidates[0]

            self.path = path
            self.etc = etc
            self.gpt = gpt
            self.auxdata = auxdata
            self.properties = auxdata_properties
            return

    def __read_config(self):
        """
        This method reads the config.ini to examine the snap paths.
        If the snap paths are not in the config.ini or the paths are wrong they will be automatically created.

        Returns
        -------
        None

        """
        for attr in self.identifiers:
            self.__read_config_attr(attr, 'SNAP')

        snap_properties = {}
        if 'OUTPUT' in ConfigHandler.sections:
            for key in ConfigHandler.keys('OUTPUT'):
                snap_properties[key] = ConfigHandler.get('OUTPUT', key)
        if len(snap_properties.keys()) > 0:
            setattr(self, 'snap_properties', snap_properties)
    
    def __read_config_attr(self, attr, section):
        if section in ConfigHandler.sections:
            if attr in ConfigHandler.keys(section):
                val = ConfigHandler.get(section, attr)
                if attr in ['path', 'gpt']:
                    exist = os.path.isfile(val)
                elif attr == 'auxdata':
                    exist = isinstance(val, list)
                else:
                    exist = os.path.isdir(val)
                if exist:
                    setattr(self, attr, val)
    
    def __update_config(self):
        for section in ['SNAP', 'OUTPUT', 'URL']:
            if section not in ConfigHandler.sections:
                ConfigHandler.add_section(section)
        
        for key in self.identifiers + ['auxdatapath', 'properties']:
            ConfigHandler.set('SNAP', key, getattr(self, key), True)
        
        demPath = os.path.join(self.auxdatapath, 'dem')
        landCoverPath = os.path.join(self.auxdatapath, 'LandCover')
        
        for key, val in self.snap_properties.items():
            if val is None:
                value = None
            elif '${AuxDataPath}' in val:
                value = val.replace('${AuxDataPath}', self.auxdatapath)
            elif '${demPath}' in val:
                value = val.replace('${demPath}', demPath)
            elif '${landCoverPath}' in val:
                value = val.replace('${landCoverPath}', landCoverPath)
            else:
                continue

            if value is not None:
                value = value.replace('\\', '/')
            ConfigHandler.set('OUTPUT', key=key, value=value)

    def __read_snap_properties(self):
        """
        Read the properties file.

        Returns
        -------
        None

        """

        with open(self.properties) as config:
            lines = config.read().split('\n')
            self.snap_properties = {}
            for line in lines:
                if not line.startswith('#') and len(line) > 0:
                    try:
                        key, val = line.split(' = ')
                    except ValueError:
                        if '=' in line:
                            line = line.strip(' =')
                            key = line
                            val = None
                        else:
                            continue
                    self.snap_properties[key] = val
                else:
                    pass
