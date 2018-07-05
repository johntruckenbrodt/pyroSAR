# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 10:10:41 2017

@author: ibari
"""
import os
import warnings
from distutils.spawn import find_executable
from os.path import expanduser
import ConfigParser
import io

__LOCAL__ = ['sensor', 'projection', 'orbit', 'polarizations', 'acquisition_mode',
             'start', 'stop', 'product', 'spacing', 'samples', 'lines']

__CONFIG_NANE__ = "config.ini"
__PATH__ = os.path.join(expanduser("~"), '.pyrosar')

__CONFIG__ = os.path.join(__PATH__, __CONFIG_NANE__)


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
            return self.__class__.__name__ + "()"

    def __dir__(self):
        return list(self.keys())


# LOOKUP
snap_suffix = {'Apply-Orbit-File': 'Orb',
               'Calibration': 'Cal',
               'Cross-Correlation': '',
               'LinearToFromdB': 'dB',
               'Remove-GRD-Border-Noise': 'bnr',
               'SAR-Simulation': 'Sim',
               'SARSim-Terrain-Correction': 'TC',
               'Subset': '',
               'Terrain-Correction': 'TC',
               'Terrain-Flattening': 'TF',
               'Read': '',
               'Write': ''}

snap = Storage(suffix=snap_suffix)

LOOKUP = Storage(snap=snap,

                 attributes={'sensor': 'TEXT',
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

# URL
dem_url = Storage(ace2='http://step.esa.int/auxdata/dem/ACE2/5M/',
                  ace='http://step.esa.int/auxdata/dem/ACE30/',
                  srtm3_FTP='xftp.jrc.it',
                  srtm3='http://srtm.csi.cgiar.org/SRT-ZIP/SRTM_V41/SRTM_Data_GeoTiff/',
                  srtm1Hgt='http://step.esa.int/auxdata/dem/SRTMGL1/', )

orbit_url = Storage(doris='http://step.esa.int/auxdata/orbits/Doris/vor',
                    ers1='http://step.esa.int/auxdata/orbits/ers_precise_orb/ERS1',
                    ers2='http://step.esa.int/auxdata/orbits/ers_precise_orb/ERS2',
                    s1_poe='http://step.esa.int/auxdata/orbits/Sentinel-1/POEORB/',
                    s1_res='http://step.esa.int/auxdata/orbits/Sentinel-1/RESORB/')

auxcal_url = Storage(s1='http://step.esa.int/auxdata/auxcal/S1/',
                     envisat='http://step.esa.int/auxdata/auxcal/ENVISAT/',
                     ers='http://step.esa.int/auxdata/auxcal/ERS/')

URL = Storage(dem=dem_url,
              orbit=orbit_url,
              auxcal=auxcal_url)

# Set Variable Names
SECTION = Storage(snap='SNAP')

KEY_SNAP = Storage(etc_path='etc_path',
                   auxdata_path='auxdata_path',
                   config_path='config_path')

KEY = Storage(snap=KEY_SNAP)

CONFIG = Storage(section=SECTION,
                 key=KEY)

# Merge
STORAGE = Storage(URL=URL,
                  LOOKUP=LOOKUP,
                  CONFIG=CONFIG)


# Class Definitions
class ExamineExe(object):
    def __init__(self):
        # todo: Update Docstrings
        pass

    @staticmethod
    def examine(name):

        if isinstance(name, str):
            name = [name]

        executable_list = []

        for item in name:
            executable_temp = find_executable(item)
            executable_list.append(executable_temp)

            # Check True values
            True_values = [item for item in executable_list if item is not None]

        if len(True_values) > 1:
            raise ValueError(
                "There are more than one instances installed. Define which one you want to use with self.set_path(...)")

        else:
            status = any(item is not None for item in executable_list)

            try:
                return status, os.path.abspath(True_values[0])

            except IndexError:
                warnings.warn(
                    "One of the executables {0} should be installed. You can download it from http://step.esa.int/main/toolboxes/snap/ or you can specify a path with snap_config.set_path(path_to_snap)".format(
                        name), UserWarning)


class Config(object):
    def __init__(self, config_file=__CONFIG__, path=__PATH__):
        self.cfg = config_file
        self.cfgpath = path

        Config.__make_dir(self.cfgpath)
        self.__create_config()

    @staticmethod
    def __make_dir(path):
        if not os.path.exists(path):
            os.makedirs(path)

        else:
            pass

    def __create_config(self):
        if not os.path.isfile(self.cfg):
            with open(self.cfg, 'w'):
                pass
        else:
            pass

    def print_config(self):
        cfgp = ConfigParser.RawConfigParser()
        cfgp.read(self.cfg)

        print("List all contents:")

        for section in cfgp.sections():
            print("\n Section: {0}".format(section))

            for options in cfgp.options(section):
                print("\t x {0} :: {1} :: {2}".format(options, cfgp.get(section, options), str(type(options))))

    def add_section(self, section):
        cfgp = ConfigParser.RawConfigParser()
        cfgp.read(self.cfg)

        if section not in cfgp.sections():
            cfgp.add_section(section)
            with open(self.cfg, 'wb') as cfg:
                cfgp.write(cfg)

        else:
            pass
            # warnings.warn("Section already exists.")

    def set(self, section, variable, value, overwrite=True):
        cfgp = ConfigParser.RawConfigParser()
        cfgp.read(self.cfg)

        if section not in cfgp.sections():
            raise AttributeError(
                "The section name is not in the config file. You must define it with self.create_section().")
        else:
            if variable in cfgp.options(section):
                if overwrite:
                    cfgp.set(section, variable, value)
                    with open(self.cfg, 'wb') as cfg:
                        cfgp.write(cfg)
                else:
                    warnings.warn("Value already exists.")

            else:
                cfgp.set(section, variable, value)
                with open(self.cfg, 'wb') as cfg:
                    cfgp.write(cfg)

    def get(self, section, key=None):
        cfgp = ConfigParser.RawConfigParser()
        cfgp.read(self.cfg)

        if section not in cfgp.sections():
            raise AttributeError(
                "The section name is not in the config file. You must define it with self.create_section().")
        else:
            if key is None:
                return dict(cfgp.items(section))

            else:
                items = dict(cfgp.items(section))
                return items[key]

# config = Config()
#
# config.print_config()
# config.get('SNAP', CONFIG_VAR.snap.etc_path)
#
# config.add_section('SNAP')
# config.set('SNAP', 'etc_path', __PATH__)
# config.set('SNAP', 'auxdata_path', 'some/path/to/auxdata')
# config.set('SNAP', 'snap_path', __PATH__, overwrite=True)
# config.set('SNAP', 'etc_path', 'change', overwrite=False)
#
# config.add_section('ENVI')
# config.set('ENVI', 'path', 'some/path/to/envi')
