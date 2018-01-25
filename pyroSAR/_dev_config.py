# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 10:10:41 2017

@author: ibari
"""
import os
import platform
from distutils.spawn import find_executable

import numpy as np

OS_SYSTEM = platform.system()

__LOCAL__ = ['sensor', 'projection', 'orbit', 'polarizations', 'acquisition_mode',
             'start', 'stop', 'product', 'spacing', 'samples', 'lines']


class Storage(dict):
    # TODO: Update docstrings
    """
    Dict class with point access to store the lookups, pattern and URLs

    Attributes
    ----------
    STORAGE.LOOKUP : dict (with point access)
        All lookuptable merged in a dict:
            * suffix : SNAP process suffix names.
            * archive : ...

    STORAGE.PATTERN : dict (with point access)
        All pattern merged in a dict:
            * ceos_ers
            * ceos_ers_pid
            * ceos_psr
            * esa
            * esa_pid
            * tsx
            * projection

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
    To import the instance use: from radoptics.config import STORAGE
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


# ==============================================================================
# LOOKUP
# ==============================================================================
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
               'Write': ''},

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

# ==============================================================================
# URL
# ==============================================================================
dem = Storage(ace2='http://step.esa.int/auxdata/dem/ACE2/5M/',
              ace='http://step.esa.int/auxdata/dem/ACE30/',
              srtm3_FTP='xftp.jrc.it',
              srtm3='http://srtm.csi.cgiar.org/SRT-ZIP/SRTM_V41/SRTM_Data_GeoTiff/',
              srtm1Hgt='http://step.esa.int/auxdata/dem/SRTMGL1/', )

orbit = Storage(doris='http://step.esa.int/auxdata/orbits/Doris/vor',
                ers1='http://step.esa.int/auxdata/orbits/ers_precise_orb/ERS1',
                ers2='http://step.esa.int/auxdata/orbits/ers_precise_orb/ERS2',
                s1_poe='http://step.esa.int/auxdata/orbits/Sentinel-1/POEORB/',
                s1_res='http://step.esa.int/auxdata/orbits/Sentinel-1/RESORB/')

auxcal = Storage(s1='http://step.esa.int/auxdata/auxcal/S1/',
                 envisat='http://step.esa.int/auxdata/auxcal/ENVISAT/',
                 ers='http://step.esa.int/auxdata/auxcal/ERS/')

URL = Storage(dem=dem,
              orbit=orbit,
              auxcal=auxcal)

# ==============================================================================
# Merge
# ==============================================================================
STORAGE = Storage(URL=URL,
                  LOOKUP=LOOKUP)


# ==============================================================================
# Class Definitions
# ==============================================================================
class ExamineExe(object):
    def __init__(self):
        self.SNAP_EXECUTABLE = ['snap64.exe', 'snap32.exe', 'snap.exe', 'snap']

    def examine(name):

        executable_list = []
        if isinstance(name, (tuple, list)):
            for item in name:
                executable_temp = find_executable(item) is not None
                executable_list.append(executable_temp)

                # Check True values
                True_values = [item for item in executable_list if item]
                # True_values = executable[np.where(executable is True)]

            if len(True_values) > 1:
                raise ValueError(
                    "There are more than one instances installed. Define with one you want to \
                    use with self.set_path(...)")

            else:
                status = any(item == True for item in executable_list)

                try:
                    temp_loc = [item for item, executable_list in enumerate(executable_list) if executable_list][0]

                except IndexError:
                    raise ValueError("One of the executables {0} must be installed.".format(name))

                return status, os.path.abspath(find_executable(name[temp_loc]))

        else:
            status = find_executable(name) is not None
            if status is False:
                raise ValueError("The executables {0} must be installed.".format(name))

            return status, os.path.abspath(find_executable(name))