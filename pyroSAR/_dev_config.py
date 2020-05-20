# -*- coding: utf-8 -*-
###############################################################################
# pyroSAR configuration handling

# Copyright (c) 2018-2020, the pyroSAR Developers.

# This file is part of the pyroSAR Project. It is subject to the
# license terms in the LICENSE.txt file found in the top-level
# directory of this distribution and at
# https://github.com/johntruckenbrodt/pyroSAR/blob/master/LICENSE.txt.
# No part of the pyroSAR project, including this file, may be
# copied, modified, propagated, or distributed except according
# to the terms contained in the LICENSE.txt file.
###############################################################################
import os
import sys
import json

# Python 3 comparability

if sys.version_info >= (3, 0):
    import configparser as ConfigParser
else:
    import ConfigParser

__LOCAL__ = ['sensor', 'projection', 'orbit', 'polarizations', 'acquisition_mode',
             'start', 'stop', 'product', 'spacing', 'samples', 'lines']

# a pattern to search for pyroSAR processing products and extract metadata attributes from the file name
product_pattern = r'(?:.*[/\\]|)' \
                  r'(?P<outname_base>' \
                  r'(?P<sensor>[A-Z0-9]{1,4})_+' \
                  r'(?P<acquisition_mode>[A-Z0-9]{1,4})_+' \
                  r'(?P<orbit>[AD])_' \
                  r'(?P<start>[0-9T]{15})' \
                  r'(?:_(?P<extensions>\w*?)|)' \
                  r')_*' \
                  r'(?:(?P<polarization>[HV]{2})_' \
                  r'(?P<proc_steps>\w*)|)' \
                  r'(?P<filetype>(?:.tif|.nc|))$'


class Storage(dict):
    """
    Dict class with point access to store the lookups, pattern and URLs

    Attributes
    ----------
    STORAGE.LOOKUP : Storage
        All lookup table merged in a Storage class instance:
            * snap : SNAP process.
            * attributes : Attributes for different sensors.
    STORAGE.URL : dict (with point access)
        All URLs for DEMs, orbit files etc.:
            * dem : URL to download specific DEMs:
                * strm3
                * ace2
                * strm3_FTP
                * strm1HGT
                * ace
            * orbit : URL to download the orbit files:
                * ers1
                * ers2
                * s1_poe
                * s1_pres
                * doris
            * auxcal : URL to download the auxcal data:
                * s1
                * envisat
                * ers

    Note
    ----
    There may be additional attributes not listed above depending of the
    specific solver. Since this class is essentially a subclass of dict
    with attribute accessors, one can see which attributes are available
    using the `keys()` method.
    """
    
    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name)
    
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
    
    def __repr__(self):
        if self.keys():
            m = max(map(len, list(self.keys()))) + 1
            return '\n'.join([k.rjust(m) + ': ' + repr(v)
                              for k, v in sorted(self.items())])
        else:
            return self.__class__.__name__ + '()'
    
    def __dir__(self):
        return list(self.keys())


# ==============================================================================
# LOOKUP
# ==============================================================================
LOOKUP = Storage(attributes={'sensor': 'TEXT',
                             'orbit': 'TEXT',
                             'acquisition_mode': 'TEXT',
                             'start': 'TEXT',
                             'stop': 'TEXT',
                             'product': 'TEXT',
                             'samples': 'INTEGER',
                             'lines': 'INTEGER',
                             'outname_base': 'TEXT PRIMARY KEY',
                             'scene': 'TEXT',
                             'hh': 'INTEGER',
                             'vv': 'INTEGER',
                             'hv': 'INTEGER',
                             'vh': 'INTEGER'})

# ==============================================================================
# URL
# ==============================================================================
dem = Storage(ace2='https://step.esa.int/auxdata/dem/ACE2/5M/',
              ace='https://step.esa.int/auxdata/dem/ACE30/',
              srtm3_FTP='xftp.jrc.it',
              srtm3='http://srtm.csi.cgiar.org/SRT-ZIP/SRTM_V41/SRTM_Data_GeoTiff/',
              srtm1Hgt='https://step.esa.int/auxdata/dem/SRTMGL1/', )

orbit = Storage(doris='https://step.esa.int/auxdata/orbits/Doris/vor',
                ers1='https://step.esa.int/auxdata/orbits/ers_precise_orb/ERS1',
                ers2='https://step.esa.int/auxdata/orbits/ers_precise_orb/ERS2',
                s1_poe='https://step.esa.int/auxdata/orbits/Sentinel-1/POEORB/',
                s1_res='https://step.esa.int/auxdata/orbits/Sentinel-1/RESORB/')

auxcal = Storage(s1='https://step.esa.int/auxdata/auxcal/S1/',
                 envisat='https://step.esa.int/auxdata/auxcal/ENVISAT/',
                 ers='https://step.esa.int/auxdata/auxcal/ERS/')

URL = Storage(dem=dem,
              orbit=orbit,
              auxcal=auxcal)

# ==============================================================================
# Merge
# ==============================================================================
STORAGE = Storage(URL=URL,
                  LOOKUP=LOOKUP)


class Singleton(type):
    """
    Define an Instance operation that lets clients access its unique instance.
    https://sourcemaking.com/design_patterns/singleton/python/1
    """
    
    def __init__(cls, name, bases, attrs, **kwargs):
        super().__init__(name, bases, attrs)
        cls._instance = None
    
    def __call__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = super().__call__(*args, **kwargs)
        return cls._instance


class ConfigHandler(metaclass=Singleton):
    """
    ConfigHandler is a configuration handler for pyroSAR. It is intended to be called by a class's '__init__' and
    set or get the configuration parameters throughout an entire package.
    The primary goal with ConfigHandler is to load a single, consistent configuration environment to be passed 
    amongst ALL objects within a package.
        
    ConfigHandler is a SINGLETON, meaning once instantiated, THE SAME OBJECT
    will be returned to every class object calling it.

    Parameters
    ----------
    path : str or None
        A path where the .pyrosar directory will be created. If None (default) it will be created in the user home
        directory.
    config_fname : str
        Name of the config file. Default is 'config.ini'.
    
    Methods
    -------
    make_dir : Create a .pyrosar directory in home directory.
    create_config : Create a config.ini file in .pyrosar directory.
    open : Open the config.ini file.
    add_section : Create a new section in the configuration.
    set : Set an option in the configuration.
    remove_option : Remove an option in the configuration.

    Notes
    -----
    The syntax is the same as in ConfigParser. Here, keys are called options.

    """
    
    # Define __setter to control changeable keys (optional)
    # __setter = ["etc", "auxdata"]
    
    def __init__(self):
        path = os.path.join(os.path.expanduser('~'), '.pyrosar')
        
        self.__GLOBAL = {
            'path': path,
            'config_fname': 'config.ini',
            'config': os.path.join(path, 'config.ini'),
        }
        
        if not os.path.isfile(self.__GLOBAL['config']):
            self.__create_config()
        
        self.parser = ConfigParser.RawConfigParser(allow_no_value=True)
        self.parser.optionxform = str
        self.parser.read(self.__GLOBAL['config'])
    
    def __create_config(self):
        """
        Create a config.ini file in .pyrosar directory.

        Returns
        -------
        None
        """
        
        if not os.path.exists(self.__GLOBAL['path']):
            os.makedirs(self.__GLOBAL['path'])
        
        with open(self.__GLOBAL['config'], 'w'):
            pass
    
    def __str__(self):
        items = []
        for section in self.parser.sections():
            items.append('  Section: {0}\n'.format(section))
            
            for options in self.parser.options(section):
                items.append('    x {0} :: {1} :: {2}\n'
                             .format(options,
                                     self.parser.get(section, options),
                                     str(type(options))))
        out = 'Class    : Config\n' \
              'Path     : {0}\n' \
              'Sections : {1}\n' \
              'Contents : \n{2}' \
            .format(self.__GLOBAL['config'],
                    len(self.parser.sections()),
                    ''.join(items))
        
        return out
    
    def __getitem__(self, section):
        if not self.parser.has_section(section):
            raise AttributeError('Section {0} does not exist.'.format(str(section)))
        return dict(self.parser.items(section))
    
    @property
    def sections(self):
        return self.parser.sections()
    
    def keys(self, section):
        """
        Get all keys (options) of a section.

        Parameters
        ----------
        section : str
            Section name.

        Returns
        -------
        list : options (keys) of a section.

        """
        return self.parser.options(section)
    
    def open(self):
        """
        Open the config.ini file. This method will open the config.ini file in a external standard app (text editor).

        Returns
        -------
        os.startfile

        """
        
        os.startfile(self.__GLOBAL['config'])
    
    def add_section(self, section):
        """
        Create a new section in the configuration.

        Parameters
        ----------
        section : str
            Section name

        Returns
        -------
        None

        """
        if not self.parser.has_section(section):
            self.parser.add_section(section)
            self.write()
        else:
            raise RuntimeError('section already exists')
    
    @property
    def file(self):
        return self.__GLOBAL['config']
    
    def set(self, section, key, value, overwrite=False):
        """
        Set an option.

        Parameters
        ----------
        section : str
            Section name.
        key : str
            the attribute name
        value :
            the attribute value
        overwrite : bool
            If True and the defined key exists the value will be overwritten.

        Returns
        -------

        """
        if not self.parser.has_section(section):
            raise AttributeError('Section {0} does not exist.'.format(str(section)))
        
        if isinstance(value, list):
            value = json.dumps(value)
        
        if key in self.parser.options(section) and not overwrite:
            raise RuntimeError('Value already exists.')
        
        self.parser.set(section, key, value)
        self.write()
    
    def remove_option(self, section, key):
        """
        Remove an option and key.

        Parameters
        ----------
        section : str
            Section name.
        key : str
            Key value.
        
        Returns
        -------
        
        """
        if not self.parser.has_section(section):
            raise AttributeError('Section {0} does not exist.'.format(str(section)))
        
        if key not in self.parser.options(section):
            raise AttributeError('Key {0} does not exist.'.format(str(key)))
        
        self.parser.remove_option(section, key)
        self.write()
    
    def remove_section(self, section):
        """
        remove a section
        
        Parameters
        ----------
        section: str
            Section name.

        Returns
        -------

        """
        self.parser.remove_section(section)
        self.write()
    
    def write(self):
        if sys.version_info >= (3, 0):
            with open(self.__GLOBAL['config'], 'w', encoding='utf8') as out:
                self.parser.write(out)
        else:
            with open(self.__GLOBAL['config'], 'w') as out:
                self.parser.write(out)
