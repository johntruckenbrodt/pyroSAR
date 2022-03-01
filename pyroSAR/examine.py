###############################################################################
# Examination of SAR processing software
# Copyright (c) 2019-2021, the pyroSAR Developers.

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
import platform
import re
import warnings
import subprocess as sp
import pkg_resources

from pyroSAR._dev_config import ConfigHandler
from spatialist.ancillary import finder

import logging

log = logging.getLogger(__name__)

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
        self.__read_config()
        
        # if SNAP could not be identified from the config attributes, do a system search for it
        if not self.__is_identified():
            log.debug('identifying SNAP')
            self.__identify_snap()
        
        # if the auxdatapath attribute was not yet set, create a default directory
        if not hasattr(self, 'auxdatapath'):
            self.auxdatapath = os.path.join(os.path.expanduser('~'), '.snap', 'auxdata')
            os.makedirs(self.auxdatapath, exist_ok=True)
        
        # if the SNAP auxdata properties attribute was not yet identified,
        # point it to the default file delivered with pyroSAR
        if not hasattr(self, 'properties'):
            # log.info('using default properties file..')
            template = os.path.join('snap', 'data', 'snap.auxdata.properties')
            self.properties = pkg_resources.resource_filename(__name__, template)
        
        # if the SNAP suffices attribute was not yet identified,
        # point it to the default file delivered with pyroSAR
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
    
    def __getattr__(self, item):
        raise AttributeError("'ExamineSnap' object has no attribute '{}'".format(item))
    
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
        
        if len(executables) == 0:
            log.debug("could not detect any potential 'snap' executables")
        
        # for each possible SNAP executable, check whether additional files and directories exist relative to it
        # to confirm whether it actually is a ESA SNAP installation or something else like e.g. the Ubuntu App Manager
        for path in executables:
            log.debug('checking candidate {}'.format(path))
            if os.path.islink(path):
                path = os.path.realpath(path)
            
            # check whether a directory etc exists relative to the SNAP executable
            etc = os.path.join(os.path.dirname(os.path.dirname(path)), 'etc')
            if not os.path.isdir(etc):
                log.debug("could not find the 'etc' directory")
                continue
            
            # check the content of the etc directory
            auxdata = os.listdir(etc)
            if 'snap.auxdata.properties' not in auxdata:
                log.debug("could not find the 'snap.auxdata.properties' file")
                continue
            else:
                auxdata_properties = os.path.join(etc, 'snap.auxdata.properties')
            
            # identify the gpt executable
            gpt_candidates = finder(os.path.dirname(path), ['gpt', 'gpt.exe'])
            if len(gpt_candidates) == 0:
                log.debug("could not find the 'gpt' executable")
                continue
            else:
                gpt = gpt_candidates[0]
            
            self.path = path
            self.etc = etc
            self.gpt = gpt
            self.auxdata = auxdata
            self.properties = auxdata_properties
            return
        
        log.warning('SNAP could not be identified. If you have installed it please add the path to the SNAP '
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
                    # log.info('setting attribute {}'.format(attr))
                    setattr(self, attr, val)
    
    def __update_config(self):
        for section in self.sections:
            if section not in config.sections:
                # log.info('creating section {}..'.format(section))
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
            # log.info('updating attribute {0}:{1}..'.format(section, attr))
            # log.info('  {0} -> {1}'.format(repr(config[section][attr]), repr(value)))
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
    
    def get_version(self, module):
        """
        Read the version and date of different SNAP modules.
        This scans a file 'messages.log', which is re-written every time SNAP is started.
        
        Parameters
        ----------
        module: str
            one of the following
            
            - core
            - desktop
            - rstbx
            - s1tbx
            - s2tbx
            - s3tbx

        Returns
        -------
        dict
            a dictionary with keys 'version' and 'date'
        """
        # base search patterns for finding the right lines
        patterns = {'core': r'org\.esa\.snap\.snap\.core',
                    'desktop': r'org\.esa\.snap\.snap\.ui',
                    'rstb': r'org\.csa\.rstb\.rstb\.kit',
                    's1tbx': r'org\.esa\.s1tbx\.s1tbx\.kit',
                    's2tbx': r'org\.esa\.s2tbx\.s2tbx\.kit',
                    's3tbx': r'org\.esa\.s3tbx\.s3tbx\.kit'}
        
        if module in patterns.keys():
            pattern = patterns[module]
            pattern += r' \[(?P<version>[0-9.]+) [0-9.]+ (?P<date>[0-9]{12})'
        else:
            raise RuntimeError('module not supported')
        
        system = platform.system()
        if system in ['Linux', 'Darwin']:
            path = os.path.join(os.path.expanduser('~'), '.snap', 'system')
        elif system == 'Windows':
            path = os.path.join(os.environ['APPDATA'], 'SNAP')
        else:
            raise RuntimeError('operating system not supported')
        
        conda_env_path = os.environ.get('CONDA_PREFIX')
        if conda_env_path is not None and conda_env_path in self.gpt:
            fname = os.path.join(conda_env_path, 'snap', '.snap', 'system', 'var', 'log', 'messages.log')
        else:
            fname = os.path.join(path, 'var', 'log', 'messages.log')
        
        if not os.path.isfile(fname):
            try:
                # This will start SNAP and immediately stop it because of the invalid argument.
                # Currently this seems to be the only way to create the messages.log file if it does not exist.
                sp.check_call([self.path, '--nosplash', '--dummytest', '--console', 'suppress'])
            except sp.CalledProcessError:
                pass
        
        if not os.path.isfile(fname):
            raise RuntimeError("cannot find 'messages.log' to read SNAP module versions from.")
        
        with open(fname, 'r') as m:
            content = m.read()
        match = re.search(pattern, content)
        if match is None:
            raise RuntimeError('cannot read version information from {}.\nPlease restart SNAP.'.format(fname))
        return match.groupdict()


class ExamineGamma(object):
    def __init__(self):
        home_sys = os.environ.get('GAMMA_HOME')
        if home_sys is not None and not os.path.isdir(home_sys):
            warnings.warn('found GAMMA_HOME environment variable, but directory does not exist')
            home_sys = None
        
        self.__read_config()
        
        if hasattr(self, 'home'):
            if home_sys is not None and self.home != home_sys:
                log.info('the value of GAMMA_HOME is different to that in the pyroSAR configuration;\n'
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
                raise RuntimeError('could not read GAMMA installation directory')
        self.version = re.search('GAMMA_SOFTWARE-(?P<version>[0-9]{8})',
                                 getattr(self, 'home')).group('version')
        
        if not hasattr(self, 'gdal_config'):
            gdal_config = '/usr/bin/gdal-config'
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
