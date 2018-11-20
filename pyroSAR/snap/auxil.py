import sys
import ast
import json
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

from osgeo import gdal
from osgeo.gdalconst import GA_Update

from pyroSAR import identify, ConfigHandler
from .._dev_config import LOOKUP

from spatialist.auxil import gdal_translate
from spatialist.ancillary import dissolve, finder


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


def insert_node(workflow, node, before=None, after=None):
    if before and not after:
        predecessor = workflow.find('.//node[@id="{}"]'.format(before))
        position = list(workflow).index(predecessor) + 1
        workflow.insert(position, node)
        newnode = workflow[position]
        # set the source product for the new node
        newnode.find('.//sources/sourceProduct').attrib['refid'] = predecessor.attrib['id']
        # set the source product for the node after the new node
        successor = workflow[position + 1]
        if successor.tag == 'node':
            successor.find('.//sources/sourceProduct').attrib['refid'] = newnode.attrib['id']
    elif after and not before:
        successor = workflow.find('.//node[@id="{}"]'.format(after))
        position = list(workflow).index(successor)
        workflow.insert(position, node)
        newnode = workflow[position]
        # set the source product for the new node
        predecessor = workflow[position - 1]
        source_id = predecessor.attrib['id'] if predecessor.tag == 'node' else None
        newnode.find('.//sources/sourceProduct').attrib['refid'] = source_id
        # set the source product for the node after the new node
        successor.find('.//sources/sourceProduct').attrib['refid'] = newnode.attrib['id']
    else:
        raise RuntimeError('cannot insert node if both before and after are set')


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
                infile = os.path.join('https://step.esa.int/auxdata/dem/SRTMGL1', file)
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
                        'https://step.esa.int/auxdata/orbits/Sentinel-1/POEORB/remote_contentVersion.txt')
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
                                    'https://step.esa.int/auxdata/orbits/Sentinel-1/POEORB/{}'.format(orb_id))
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
    # print('- processing workflow {}'.format(os.path.basename(xmlfile)))
    proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = proc.communicate()
    out = out.decode('utf-8') if isinstance(out, bytes) else out
    err = err.decode('utf-8') if isinstance(err, bytes) else err
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
        # print('- converting to GTiff')
        suffix = parse_suffix(workflow)
        translateoptions = {'options': ['-q', '-co', 'INTERLEAVE=BAND', '-co', 'TILED=YES'],
                            'format': 'GTiff',
                            'noData': 0}
        for item in finder(outname, ['*.img']):
            pol = re.search('[HV]{2}', item).group()
            name_new = outname.replace(suffix, '{0}_{1}.tif'.format(pol, suffix))
            gdal_translate(item, name_new, translateoptions)
        shutil.rmtree(outname)
    elif format == 'GeoTiff-BigTIFF':
        ras = gdal.Open(outname + '.tif', GA_Update)
        for i in range(1, ras.RasterCount + 1):
            ras.GetRasterBand(i).SetNoDataValue(0)
        ras = None


class ExamineSnap(object):
    """
    Class to check if ESA SNAP is installed.
    Upon initialization, this class searches for relevant binaries and the accompanying
    relative directory structure, which uniquely identify an ESA SNAP installation on a system.
    First, all relevant file and folder names are read from the pyroSAR config file if it exists
    and their existence is verified.
    If this fails, a system check is performed to find relevant binaries in the system PATH variable and
    additional files and folders relative to them.
    Furthermore, a snap.auxdata.properties file is scanned for auxiliary data URLs and local storage location.
    This is used by SNAP to manage data from e.g. the SRTM mission. In case SNAP is not installed, the respective
    information is read from a default file delivered with pyroSAR. This has the advantage of using the SNAP download
    URLs and local directory structure without having SNAP installed such that it can be adapted by other SAR software.
    """
    
    def __init__(self):
        
        # define some attributes which identify SNAP
        self.identifiers = ['path', 'gpt', 'etc', 'auxdata']
        
        # try reading all necessary attributes from the config file
        # print('reading config..')
        self.__read_config()
        
        # if SNAP could not be identified from the config attributes, do a system search for it
        if not self.__is_identified():
            # print('identifying SNAP..')
            self.__identify_snap()
        
        # if the auxdatapath attribute was not yet set, create a default directory
        if not hasattr(self, 'auxdatapath'):
            self.auxdatapath = os.path.join(os.path.expanduser('~'), '.snap', 'auxdata')
            if not os.path.isdir(self.auxdatapath):
                os.makedirs(self.auxdatapath)
        
        # if the SNAP auxdata properties attribute was not yet identified,
        # point it to the default file delivered with pyroSAR
        if not hasattr(self, 'properties'):
            # print('using default properties file..')
            template = 'data/snap.auxdata.properties'
            self.properties = pkg_resources.resource_filename(__name__, template)
        
        # update the snap properties; this reads the 'properties' file and looks for any changes,
        # which are then updated for the object
        self.__update_snap_properties()
        
        # update the config file: this scans for config changes and re-writes the config file if any are found
        self.__update_config()
    
    def __is_identified(self):
        """
        Check if SNAP has been properly identified, i.e. all paths in self.identifiers
        have been detected and confirmed.
        
        Returns
        -------
        bool
        """
        return sum([hasattr(self, x) for x in self.identifiers]) == len(self.identifiers)
    
    def __identify_snap(self):
        """
        do a comprehensive search for an ESA SNAP installation
        
        Returns
        -------
        bool
            has the SNAP properties file been changed?
        """
        # create a list of possible SNAP executables
        defaults = ['snap64.exe', 'snap32.exe', 'snap.exe', 'snap']
        paths = os.environ['PATH'].split(os.path.pathsep)
        options = [os.path.join(path, option) for path in paths for option in defaults]
        options = [x for x in options if os.path.isfile(x)]
        
        if not hasattr(self, 'path') or not os.path.isfile(self.path):
            executables = options
        else:
            executables = [self.path] + options
        
        # for each possible SNAP executable, check whether additional files and directories exist relative to it
        # to confirm whether it actually is a ESA SNAP installation or something else like e.g. the Ubuntu App Manager
        for path in executables:
            if os.path.islink(path):
                path = os.path.realpath(path)
            
            # check whether a directory etc exists relative to the SNAP executable
            etc = os.path.join(os.path.dirname(os.path.dirname(path)), 'etc')
            if not os.path.isdir(etc):
                continue
            
            # check the content of the etc directory
            auxdata = os.listdir(etc)
            if 'snap.auxdata.properties' not in auxdata:
                continue
            else:
                auxdata_properties = os.path.join(etc, 'snap.auxdata.properties')
            
            # identify the gpt executable
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
        
        warnings.warn('SNAP could not be identified')
    
    def __read_config(self):
        """
        This method reads the config.ini to examine the snap paths.
        If the snap paths are not in the config.ini or the paths are wrong they will be automatically created.

        Returns
        -------

        """
        for attr in self.identifiers + ['auxdatapath', 'properties']:
            self.__read_config_attr(attr, 'SNAP')
        
        snap_properties = {}
        if 'OUTPUT' in ConfigHandler.sections:
            snap_properties = ConfigHandler['OUTPUT']
        if len(snap_properties.keys()) > 0:
            setattr(self, 'snap_properties', snap_properties)
    
    def __read_config_attr(self, attr, section):
        """
        read an attribute from the config file and set it as an object attribute
        
        Parameters
        ----------
        attr: str
            the attribute name
        section: str
            the config section to read the attribute from
        
        Returns
        -------
        
        """
        if section in ConfigHandler.sections:
            if attr in ConfigHandler[section].keys():
                val = ConfigHandler[section][attr]
                if attr in ['path', 'gpt', 'properties']:
                    exist = os.path.isfile(val)
                elif attr == 'auxdata':
                    val = ast.literal_eval(val)
                    exist = isinstance(val, list)
                else:
                    exist = os.path.isdir(val)
                if exist:
                    # print('setting attribute {}'.format(attr))
                    setattr(self, attr, val)
    
    def __update_config(self):
        for section in ['SNAP', 'OUTPUT', 'URL']:
            if section not in ConfigHandler.sections:
                # print('creating section {}..'.format(section))
                ConfigHandler.add_section(section)
        
        for key in self.identifiers + ['auxdatapath', 'properties']:
            if hasattr(self, key):
                self.__update_config_attr(key, getattr(self, key), 'SNAP')
        
        for key, value in self.snap_properties.items():
            self.__update_config_attr(key, value, 'OUTPUT')
    
    @staticmethod
    def __update_config_attr(attr, value, section):
        if isinstance(value, list):
            value = json.dumps(value)
        
        if attr not in ConfigHandler[section].keys() or ConfigHandler[section][attr] != value:
            # print('updating attribute {0}:{1}..'.format(section, attr))
            # print('  {0} -> {1}'.format(repr(ConfigHandler[section][attr]), repr(value)))
            ConfigHandler.set(section, key=attr, value=value, overwrite=True)
    
    def __update_snap_properties(self):
        """
        Read the snap.auxdata.properties file entries to object attributes

        Returns
        -------

        """
        pattern = r'^(?P<key>[\w\.]*)\s*=\s*(?P<value>.*)\n'
        
        if not hasattr(self, 'snap_properties'):
            self.snap_properties = {}
        
        demPath = os.path.join(self.auxdatapath, 'dem')
        landCoverPath = os.path.join(self.auxdatapath, 'LandCover')
        
        with open(self.properties, 'r') as prop:
            for line in prop:
                if re.search(pattern, line):
                    key, value = re.match(re.compile(pattern), line).groups()
                    value = value \
                        .replace('${AuxDataPath}', self.auxdatapath) \
                        .replace('${demPath}', demPath) \
                        .replace('${landCoverPath}', landCoverPath) \
                        .replace('\\', '/')
                    if not key in self.snap_properties.keys() or self.snap_properties[key] != value:
                        self.snap_properties[key] = value
