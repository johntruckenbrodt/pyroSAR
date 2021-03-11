###############################################################################
# Examination of SAR processing software

# Copyright (c) 2019-2020, the pyroSAR Developers.

# This file is part of the pyroSAR Project. It is subject to the
# license terms in the LICENSE.txt file found in the top-level
# directory of this distribution and at
# https://github.com/johntruckenbrodt/pyroSAR/blob/master/LICENSE.txt.
# No part of the pyroSAR project, including this file, may be
# copied, modified, propagated, or distributed except according
# to the terms contained in the LICENSE.txt file.
###############################################################################
import ast
import json
import os
import shutil

import re
import warnings

import pkg_resources

from pyroSAR._dev_config import ConfigHandler
from spatialist.ancillary import finder

config = ConfigHandler()


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
        
        # a list of relevant sections
        self.sections = ['SNAP', 'OUTPUT', 'SNAP_SUFFIX']
        
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
            os.makedirs(self.auxdatapath, exist_ok=True)
        
        # if the SNAP auxdata properties attribute was not yet identified,
        # point it to the default file delivered with pyroSAR
        if not hasattr(self, 'properties'):
            # print('using default properties file..')
            template = os.path.join('snap', 'data', 'snap.auxdata.properties')
            self.properties = pkg_resources.resource_filename(__name__, template)
        
        if not hasattr(self, 'suffices'):
            template = os.path.join('snap', 'data', 'snap.suffices.properties')
            fname_suffices = pkg_resources.resource_filename(__name__, template)
            with open(fname_suffices, 'r') as infile:
                content = infile.read().split('\n')
            self.__suffices = {k: v for k, v in [x.split('=') for x in content]}
        
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
        
        warnings.warn('SNAP could not be identified. If you have installed it please add the path to the SNAP '
                      'executables (bin subdirectory) to the PATH environment. '
                      'E.g. in the Linux .bashrc file add the following line:\nexport PATH=$PATH:path/to/snap/bin"')
    
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
        if 'OUTPUT' in config.sections:
            snap_properties = config['OUTPUT']
        if len(snap_properties.keys()) > 0:
            self.snap_properties = snap_properties
        
        suffices = {}
        if 'SNAP_SUFFIX' in config.sections:
            suffices = config['SNAP_SUFFIX']
        if len(suffices.keys()) > 0:
            self.__suffices = suffices
    
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
        if section in config.sections:
            if attr in config[section].keys():
                val = config[section][attr]
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
        for section in self.sections:
            if section not in config.sections:
                # print('creating section {}..'.format(section))
                config.add_section(section)
        
        for key in self.identifiers + ['auxdatapath', 'properties']:
            if hasattr(self, key):
                self.__update_config_attr(key, getattr(self, key), 'SNAP')
        
        for key, value in self.snap_properties.items():
            self.__update_config_attr(key, value, 'OUTPUT')
        
        for key in sorted(self.__suffices.keys()):
            self.__update_config_attr(key, self.__suffices[key], 'SNAP_SUFFIX')
    
    @staticmethod
    def __update_config_attr(attr, value, section):
        if isinstance(value, list):
            value = json.dumps(value)
        
        if attr not in config[section].keys() or config[section][attr] != value:
            # print('updating attribute {0}:{1}..'.format(section, attr))
            # print('  {0} -> {1}'.format(repr(config[section][attr]), repr(value)))
            config.set(section, key=attr, value=value, overwrite=True)
    
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
                    if key not in self.snap_properties.keys() or self.snap_properties[key] != value:
                        self.snap_properties[key] = value
    
    def get_suffix(self, operator):
        """
        get the file name suffix for an operator
        
        Parameters
        ----------
        operator: str
            the name of the operator

        Returns
        -------
        str
            the file suffix
        
        Examples
        --------
        >>> from pyroSAR.examine import ExamineSnap
        >>> config = ExamineSnap()
        >>> print(config.get_suffix('Terrain-Flattening'))
        'TF'
        """
        if operator in self.__suffices.keys():
            return self.__suffices[operator]
        else:
            return None


class ExamineGamma(object):
    def __init__(self):
        home_sys = os.environ.get('GAMMA_HOME')
        if home_sys is not None and not os.path.isdir(home_sys):
            warnings.warn('found GAMMA_HOME environment variable, but directory does not exist')
            home_sys = None
        
        self.__read_config()
        
        if hasattr(self, 'home'):
            if home_sys is not None and self.home != home_sys:
                print('the value of GAMMA_HOME is different to that in the pyroSAR configuration;\n'
                      '  was: {}\n'
                      '  is : {}\n'
                      'resetting the configuration and deleting parsed modules'
                      .format(self.home, home_sys))
                parsed = os.path.join(os.path.dirname(self.fname), 'gammaparse')
                shutil.rmtree(parsed)
                self.home = home_sys
        if not hasattr(self, 'home'):
            if home_sys is not None:
                setattr(self, 'home', home_sys)
            else:
                raise RuntimeError('could not read Gamma installation directory')
        self.version = re.search('GAMMA_SOFTWARE-(?P<version>[0-9]{8})',
                                 getattr(self, 'home')).group('version')
        
        if not hasattr(self, 'gdal_config'):
            gdal_config = '/usr/bin/gdal-configue'
            if not os.path.isfile(gdal_config):
                warnings.warn('could not find GDAL system installation under {}\n'
                              'please modify the pyroSAR configuration: {}'.format(gdal_config, self.fname))
            self.gdal_config = gdal_config
        self.__update_config()
    
    def __read_config(self):
        self.fname = config.file
        if 'GAMMA' in config.sections:
            attr = config['GAMMA']
            for key, value in attr.items():
                setattr(self, key, value)
    
    def __update_config(self):
        if 'GAMMA' not in config.sections:
            config.add_section('GAMMA')
        
        for attr in ['home', 'version', 'gdal_config']:
            self.__update_config_attr(attr, getattr(self, attr), 'GAMMA')
    
    @staticmethod
    def __update_config_attr(attr, value, section):
        if isinstance(value, list):
            value = json.dumps(value)
        
        if attr not in config[section].keys() or config[section][attr] != value:
            config.set(section, key=attr, value=value, overwrite=True)
