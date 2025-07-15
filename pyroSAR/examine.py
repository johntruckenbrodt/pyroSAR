###############################################################################
# Examination of SAR processing software
# Copyright (c) 2019-2024, the pyroSAR Developers.

# This file is part of the pyroSAR Project. It is subject to the
# license terms in the LICENSE.txt file found in the top-level
# directory of this distribution and at
# https://github.com/johntruckenbrodt/pyroSAR/blob/master/LICENSE.txt.
# No part of the pyroSAR project, including this file, may be
# copied, modified, propagated, or distributed except according
# to the terms contained in the LICENSE.txt file.
###############################################################################
import json
import os
import shutil
import platform
import re
import warnings
import subprocess as sp
import importlib.resources

from pyroSAR.config import ConfigHandler
from spatialist.ancillary import finder, run

import logging

log = logging.getLogger(__name__)

__config__ = ConfigHandler()


class ExamineSnap(object):
    """
    Class to check if ESA SNAP is installed.
    Upon initialization, this class searches for relevant binaries and the accompanying
    relative directory structure, which uniquely identify an ESA SNAP installation on a system.
    First, all relevant file and folder names are read from the pyroSAR config file if it exists
    and their existence is verified.
    If this fails, a system check is performed to find relevant binaries in the system PATH variable and
    additional files and folders relative to them.
    In case SNAP is not installed, a default `snap.auxdata.properties` file delivered with pyroSAR will be copied to
    `$HOME/.snap/etc` so that SNAP download URLS and local directory structure can be adapted by other software.
    
    SNAP configuration can be read and modified via the attribute `snap_properties` of type
    :class:`~pyroSAR.examine.SnapProperties` or the properties :attr:`~pyroSAR.examine.ExamineSnap.userpath` and
    :attr:`~pyroSAR.examine.ExamineSnap.auxdatapath`.
    """
    
    def __init__(self):
        # update legacy config files
        if 'OUTPUT' in __config__.sections:
            __config__.remove_section('OUTPUT')
        if 'SNAP' in __config__.sections:
            snap_keys = __config__.keys('SNAP')
            for key in ['auxdata', 'auxdatapath', 'properties']:
                if key in snap_keys:
                    __config__.remove_option(section='SNAP', key=key)
        
        # define some attributes which identify SNAP
        self.identifiers = ['path', 'gpt', 'etc']
        
        # a list of relevant sections
        self.sections = ['SNAP', 'SNAP_SUFFIX']
        
        # set attributes path, gpt, etc, __suffices
        self.__read_config()
        
        # if SNAP could not be identified from the config attributes, do a system search for it
        # sets attributes path, gpt, etc
        if not self.__is_identified():
            log.debug('identifying SNAP')
            self.__identify_snap()
        
        # if SNAP cannot be identified, copy the snap.auxdata.properties file to $HOME/.snap/etc
        if not self.__is_identified():
            self.etc = os.path.join(os.path.expanduser('~'), '.snap', 'etc')
            os.makedirs(self.etc, exist_ok=True)
            dst = os.path.join(self.etc, 'snap.auxdata.properties')
            if not os.path.isfile(dst):
                dir_data = importlib.resources.files('pyroSAR') / 'snap' / 'data'
                src = str(dir_data / 'snap.auxdata.properties')
                log.debug(f'creating {dst}')
                shutil.copyfile(src, dst)
        
        # if the SNAP suffices attribute was not yet identified,
        # point it to the default file delivered with pyroSAR
        if not hasattr(self, '__suffices'):
            dir_data = importlib.resources.files('pyroSAR') / 'snap' / 'data'
            fname_suffices = str(dir_data / 'snap.suffices.properties')
            with open(fname_suffices, 'r') as infile:
                content = infile.read().split('\n')
            self.__suffices = {k: v for k, v in [x.split('=') for x in content]}
        
        # SNAP property read/modification interface
        self.snap_properties = SnapProperties(path=os.path.dirname(self.etc))
        
        # update the config file: this scans for config changes and re-writes the config file if any are found
        self.__update_config()
    
    def __getattr__(self, item):
        if item in ['path', 'gpt']:
            msg = ('SNAP could not be identified. If you have installed it '
                   'please add the path to the SNAP executables (bin subdirectory) '
                   'to the PATH environment. E.g. in the Linux .bashrc file add '
                   'the following line:\nexport PATH=$PATH:path/to/snap/bin"')
        else:
            msg = "'ExamineSnap' object has no attribute '{}'".format(item)
        raise AttributeError(msg)
    
    def __is_identified(self):
        """
        Check if SNAP has been properly identified, i.e. all paths in `self.identifiers`
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
        # to confirm whether it actually is an ESA SNAP installation or something else like e.g. the Ubuntu App Manager
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
            config_files = os.listdir(etc)
            expected = ['snap.auxdata.properties', 'snap.clusters',
                        'snap.conf', 'snap.properties']
            for name in expected:
                if name not in config_files:
                    log.debug(f"could not find the '{name}' file")
                    continue
            
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
            return
    
    def __read_config(self):
        """
        This method reads the config.ini to examine the snap paths.
        If the snap paths are not in the config.ini or the paths are
        wrong they will be automatically created.

        Returns
        -------

        """
        for attr in self.identifiers:
            self.__read_config_attr(attr, section='SNAP')
        
        suffices = {}
        if 'SNAP_SUFFIX' in __config__.sections:
            suffices = __config__['SNAP_SUFFIX']
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
        if section in __config__.sections:
            if attr in __config__[section].keys():
                val = __config__[section][attr]
                if os.path.exists(val):
                    # log.info('setting attribute {}'.format(attr))
                    setattr(self, attr, val)
    
    def __update_config(self):
        for section in self.sections:
            if section not in __config__.sections:
                # log.info('creating section {}..'.format(section))
                __config__.add_section(section)
        
        for key in self.identifiers:
            if hasattr(self, key):
                self.__update_config_attr(key, getattr(self, key), 'SNAP')
        
        for key in sorted(self.__suffices.keys()):
            self.__update_config_attr(key, self.__suffices[key], 'SNAP_SUFFIX')
    
    @staticmethod
    def __update_config_attr(attr, value, section):
        if isinstance(value, list):
            value = json.dumps(value)
        
        if attr not in __config__[section].keys() or __config__[section][attr] != value:
            # log.info('updating attribute {0}:{1}..'.format(section, attr))
            # log.info('  {0} -> {1}'.format(repr(config[section][attr]), repr(value)))
            __config__.set(section, key=attr, value=value, overwrite=True)
    
    def get_suffix(self, operator):
        """
        get the file name suffix for an operator
        
        Parameters
        ----------
        operator: str
            the name of the operator

        Returns
        -------
        str or None
            the file suffix or None if unknown
        
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
            - rstb
            - opttbx
            - microwavetbx

        Returns
        -------
        dict
            a dictionary with keys 'version' and 'date'
        """
        # base search patterns for finding the right lines
        patterns = {'core': r'org\.esa\.snap\.snap\.core',
                    'desktop': r'org\.esa\.snap\.snap\.ui',
                    'rstb': r'org\.csa\.rstb\.rstb\.kit',
                    'opttbx': r'eu\.esa\.opt\.opttbx\.kit',
                    'microwavetbx': r'eu\.esa\.microwavetbx\.microwavetbx\.kit'}
        
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
                # Currently, this seems to be the only way to create the messages.log file if it does not exist.
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
    
    @property
    def auxdatapath(self):
        """
        Get/set the SNAP configuration for `AuxDataPath` in `snap.auxdata.properties`.
        
        Example
        -------
        >>> from pyroSAR.examine import ExamineSnap
        >>> config = ExamineSnap()
        >>> config.auxdatapath = '/path/to/snap/auxdata'
        # This is equivalent to
        >>> config.snap_properties['AuxDataPath'] = '/path/to/snap/auxdata'
        """
        out = self.snap_properties['AuxDataPath']
        if out is None:
            out = os.path.join(self.userpath, 'auxdata')
        return out
    
    @auxdatapath.setter
    def auxdatapath(self, value):
        self.snap_properties['AuxDataPath'] = value
    
    @property
    def userpath(self):
        """
        Get/set the SNAP configuration for `snap.userdir` in `snap.properties`.

        Example
        -------
        >>> from pyroSAR.examine import ExamineSnap
        >>> config = ExamineSnap()
        >>> config.userpath = '/path/to/snap/data'
        # This is equivalent to
        >>> config.snap_properties['snap.userdir'] = '/path/to/snap/data'
        """
        return self.snap_properties.userpath
    
    @userpath.setter
    def userpath(self, value):
        self.snap_properties.userpath = value


class ExamineGamma(object):
    """
    Class to check if GAMMA is installed.
    
    Examples
    --------
    >>> from pyroSAR.examine import ExamineGamma
    >>> config = ExamineGamma()
    >>> print(config.home)
    >>> print(config.version)
    
    """
    
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
        self.version = re.search('GAMMA_SOFTWARE[-/](?P<version>[0-9]{8})',
                                 getattr(self, 'home')).group('version')
        
        try:
            out, err = run(['which', 'gdal-config'], void=False)
            gdal_config = out.strip('\n')
            self.gdal_config = gdal_config
        except sp.CalledProcessError:
            raise RuntimeError('could not find command gdal-config.')
        self.__update_config()
    
    def __read_config(self):
        self.fname = __config__.file
        if 'GAMMA' in __config__.sections:
            attr = __config__['GAMMA']
            for key, value in attr.items():
                setattr(self, key, value)
    
    def __update_config(self):
        if 'GAMMA' not in __config__.sections:
            __config__.add_section('GAMMA')
        
        for attr in ['home', 'version']:
            self.__update_config_attr(attr, getattr(self, attr), 'GAMMA')
    
    @staticmethod
    def __update_config_attr(attr, value, section):
        if isinstance(value, list):
            value = json.dumps(value)
        
        if attr not in __config__[section].keys() or __config__[section][attr] != value:
            __config__.set(section, key=attr, value=value, overwrite=True)


class SnapProperties(object):
    """
    SNAP configuration interface. This class enables reading and modifying
    SNAP configuration in properties files. Modified properties are directly
    written to the files.
    Currently, the files `snap.properties` and `snap.auxdata.properties` are
    supported. These files can be found in two locations:
    
    - `<SNAP installation directory>/etc`
    - `<user directory>/.snap/etc`
    
    Configuration in the latter has higher priority and modified properties will
    always be written there so that the installation directory is not modified.

    Parameters
    ----------
    path: str
        SNAP installation directory path
    
    Examples
    --------
    >>> from pyroSAR.examine import ExamineSnap, SnapProperties
    >>> path = ExamineSnap().path
    >>> config = SnapProperties(path=path)
    >>> config['snap.userdir'] = '/path/to/snap/auxdata'
    """
    
    def __init__(self, path):
        self.pattern = r'^(?P<comment>#?)(?P<key>[\w\.]*)[ ]*=[ ]*(?P<value>.*)\n*'
        self.properties_path = os.path.join(path, 'etc', 'snap.properties')
        self.auxdata_properties_path = os.path.join(path, 'etc', 'snap.auxdata.properties')
        
        log.debug(f"reading {self.properties_path}")
        self.properties = self._to_dict(self.properties_path)
        self.auxdata_properties = self._to_dict(self.auxdata_properties_path)
        
        self._dicts = [self.properties, self.auxdata_properties]
        
        # some properties need to be read from the default user path to
        # be visible to SNAP
        pairs = [(self.userpath_properties, self.properties_path),
                 (self.userpath_auxdata_properties, self.auxdata_properties_path)]
        for default, defined in pairs:
            if default != defined:
                conf = self._to_dict(default)
                if len(conf.keys()) > 0:
                    log.debug(f"updating keys {list(conf.keys())} from {default}")
                    self.properties.update(conf)
    
    def __getitem__(self, key):
        """
        
        Parameters
        ----------
        key: str
        

        Returns
        -------

        """
        for section in self._dicts:
            if key in section:
                return section[key]
        raise KeyError(f'could not find key {key}')
    
    def __setitem__(self, key, value):
        """
        
        Parameters
        ----------
        key: str
        value: Any

        Returns
        -------

        """
        if value == self[key] and isinstance(value, type(self[key])):
            return
        if key in self.properties:
            self.properties[key] = value
        else:
            self.auxdata_properties[key] = value
        if value is not None:
            value = str(value).encode('unicode-escape').decode()
            value = value.replace(':', '\\:')
        if key in self.properties:
            path = self.userpath_properties
        elif key in self.auxdata_properties:
            path = self.userpath_auxdata_properties
        else:
            raise KeyError(f'unknown key {key}')
        if os.path.isfile(path):
            with open(path, 'r') as f:
                content = f.read()
        else:
            content = ''
        pattern = r'#?{}[ ]*=[ ]*(?P<value>.*)'.format(key)
        match = re.search(pattern, content)
        if match:
            repl = f'#{key} =' if value is None else f'{key} = {value}'
            content = content.replace(match.group(), repl)
        else:
            content += f'\n{key} = {value}'
        
        os.makedirs(os.path.dirname(path), exist_ok=True)
        log.debug(f"writing key '{key}' with value '{value}' to '{path}'")
        with open(path, 'w') as f:
            f.write(content)
    
    def _to_dict(self, path):
        """
        
        Parameters
        ----------
        path: str

        Returns
        -------
        dict
        """
        out = {}
        if os.path.isfile(path):
            with open(path, 'r') as f:
                for line in f:
                    if re.search(self.pattern, line):
                        match = re.match(re.compile(self.pattern), line)
                        comment, key, value = match.groups()
                        value = self._string_convert(value)
                        out[key] = value if comment == '' else None
        return out
    
    @staticmethod
    def _string_convert(string):
        if string.lower() == 'none':
            return None
        elif string.lower() == 'true':
            return True
        elif string.lower() == 'false':
            return False
        else:
            try:
                return int(string)
            except ValueError:
                try:
                    return float(string)
                except ValueError:
                    return string.replace('\\:', ':').replace('\\\\', '\\')
    
    def keys(self):
        """
        
        Returns
        -------
        list[str]
            all known SNAP property keys
        """
        keys = []
        for item in self._dicts:
            keys.extend(list(item.keys()))
        return sorted(keys)
    
    @property
    def userpath(self):
        key = 'snap.userdir'
        if key not in self.keys() or self[key] is None:
            return os.path.join(os.path.expanduser('~'), '.snap')
        else:
            return self[key]
    
    @userpath.setter
    def userpath(self, value):
        self['snap.userdir'] = value
    
    @property
    def userpath_auxdata_properties(self):
        return os.path.join(os.path.expanduser('~'), '.snap',
                            'etc', 'snap.auxdata.properties')
    
    @property
    def userpath_properties(self):
        return os.path.join(os.path.expanduser('~'), '.snap',
                            'etc', 'snap.properties')
