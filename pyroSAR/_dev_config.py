# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 10:10:41 2017

@author: ibari
"""
import os
import warnings
from distutils.spawn import find_executable

__LOCAL__ = ['sensor', 'projection', 'orbit', 'polarizations', 'acquisition_mode',
             'start', 'stop', 'product', 'spacing', 'samples', 'lines']


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
