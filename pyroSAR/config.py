# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 10:10:41 2017

@author: ibari
"""
import os
import platform
import numpy as np
from distutils.spawn import find_executable

OS_SYSTEM = platform.system()

__LOCAL__ = ['sensor', 'projection', 'orbit', 'polarizations', 'acquisition_mode', 
             'start', 'stop', 'product', 'spacing', 'samples', 'lines']

SNAP_EXECUTABLE = ['snap64.exe', 'snap32.exe', 'snap.exe', 'snap']
#==============================================================================
# LOOKUP
#==============================================================================
SUFFIX_LOOKUP = {'Apply-Orbit-File': 'Orb',
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

ARCHIVE_LOOKUP = {'sensor': 'TEXT',
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
                  'vh': 'INTEGER'}

#==============================================================================
# PATTERN
#==============================================================================

CEOS_ERS_PATTERN = r'(?P<product_id>(?:SAR|ASA)_(?:IM(?:S|P|G|M|_)|AP(?:S|P|G|M|_)|WV(?:I|S|W|_)|WS(?:M|S|_))_[012B][CP])' \
                   r'(?P<processing_stage_flag>[A-Z])' \
                   r'(?P<originator_ID>[A-Z\-]{3})' \
                   r'(?P<start_day>[0-9]{8})_' \
                   r'(?P<start_time>[0-9]{6})_' \
                   r'(?P<duration>[0-9]{8})' \
                   r'(?P<phase>[0-9A-Z]{1})' \
                   r'(?P<cycle>[0-9]{3})_' \
                   r'(?P<relative_orbit>[0-9]{5})_' \
                   r'(?P<absolute_orbit>[0-9]{5})_' \
                   r'(?P<counter>[0-9]{4,})\.' \
                   r'(?P<satellite_ID>[EN][12])' \
                   r'(?P<extension>(?:\.zip|\.tar\.gz|\.PS|))$'
                   
CEOS_ERS_PATTERN_PID = r'(?P<sat_id>(?:SAR|ASA))_' \
                       r'(?P<image_mode>(?:IM(?:S|P|G|M|_)|AP(?:S|P|G|M|_)|WV(?:I|S|W|_)|WS(?:M|S|_)))_' \
                       r'(?P<processing_level>[012B][CP])'
                       
CEOS_PSR_PATTERN = [r'^LED-ALPSR'
                    r'(?P<sub>P|S)'
                    r'(?P<orbit>[0-9]{5})'
                    r'(?P<frame>[0-9]{4})-'
                    r'(?P<mode>[HWDPC])'
                    r'(?P<level>1\.[015])'
                    r'(?P<proc>G|_)'
                    r'(?P<proj>[UPML_])'
                    r'(?P<orbit_dir>A|D)$',
                    r'^LED-ALOS2'
                    r'(?P<orbit>[0-9]{5})'
                    r'(?P<frame>[0-9]{4})-'
                    r'(?P<date>[0-9]{6})-'
                    r'(?P<mode>SBS|UBS|UBD|HBS|HBD|HBQ|FBS|FBD|FBQ|WBS|WBD|WWS|WWD|VBS|VBD)'
                    r'(?P<look_dir>L|R)'
                    r'(?P<level>1\.0|1\.1|1\.5|2\.1|3\.1)'
                    r'(?P<proc>[GR_])'
                    r'(?P<proj>[UPML_])'
                    r'(?P<orbit_dir>A|D)$']

ESA_PATTERN = r'(?P<product_id>(?:SAR|ASA)_(?:IM(?:S|P|G|M|_)|AP(?:S|P|G|M|_)|WV(?:I|S|W|_)|WS(?:M|S|_))_[012B][CP])' \
              r'(?P<processing_stage_flag>[A-Z])' \
              r'(?P<originator_ID>[A-Z\-]{3})' \
              r'(?P<start_day>[0-9]{8})_' \
              r'(?P<start_time>[0-9]{6})_' \
              r'(?P<duration>[0-9]{8})' \
              r'(?P<phase>[0-9A-Z]{1})' \
              r'(?P<cycle>[0-9]{3})_' \
              r'(?P<relative_orbit>[0-9]{5})_' \
              r'(?P<absolute_orbit>[0-9]{5})_' \
              r'(?P<counter>[0-9]{4,})\.' \
              r'(?P<satellite_ID>[EN][12])' \
              r'(?P<extension>(?:\.zip|\.tar\.gz|))$'
              
ESA_PATTERN_PID = r'(?P<sat_id>(?:SAR|ASA))_' \
                  r'(?P<image_mode>(?:IM(?:S|P|G|M|_)|AP(?:S|P|G|M|_)|WV(?:I|S|W|_)|WS(?:M|S|_)))_' \
                  r'(?P<processing_level>[012B][CP])'              

TSX_PATTERN = r'^(?P<sat>T[DS]X1)_SAR__' \
              r'(?P<prod>SSC|MGD|GEC|EEC)_' \
              r'(?P<var>____|SE__|RE__|MON1|MON2|BTX1|BRX2)_' \
              r'(?P<mode>SM|SL|HS|HS300|ST|SC)_' \
              r'(?P<pols>[SDTQ])_' \
              r'(?:SRA|DRA)_' \
              r'(?P<start>[0-9]{8}T[0-9]{6})_' \
              r'(?P<stop>[0-9]{8}T[0-9]{6})(?:\.xml|)$'
                       
PROJECTION = 'GEOGCS["WGS 84",' \
             'DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],' \
             'PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],' \
             'UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],' \
             'AUTHORITY["EPSG","4326"]]'              

class URL:
    def __init__(self):
        pass
    
    def dem(self):
        self.ace2_HTTP = 'http://step.esa.int/auxdata/dem/ACE2/5M/'
        self.ace_HTTP = 'http://step.esa.int/auxdata/dem/ACE30/'
        self.srtm3_FTP = 'xftp.jrc.it'
        self.srtm3_HTTP = 'http://srtm.csi.cgiar.org/SRT-ZIP/SRTM_V41/SRTM_Data_GeoTiff/'
        self.srtm1Hgt = 'http://step.esa.int/auxdata/dem/SRTMGL1/'

    def orbit(self):
        self.doris_HTTP = 'http://step.esa.int/auxdata/orbits/Doris/vor'
        self.ERS1_HTTP = 'http://step.esa.int/auxdata/orbits/ers_precise_orb/ERS1'
        self.ERS2_HTTP = 'http://step.esa.int/auxdata/orbits/ers_precise_orb/ERS2'
        self.sentinet1POE = 'http://step.esa.int/auxdata/orbits/Sentinel-1/POEORB/'
        self.sentinel1RES = 'http://step.esa.int/auxdata/orbits/Sentinel-1/RESORB/'

    def auxcal(self):
        self.sentinel1 = 'http://step.esa.int/auxdata/auxcal/S1/'
        self.ENVISAT = 'http://step.esa.int/auxdata/auxcal/ENVISAT/'
        self.ERS = 'http://step.esa.int/auxdata/auxcal/ERS/'

#==============================================================================
# Class Definitions
#==============================================================================           
class ExamineExe(object):
    def __init__(self):
        pass
    
    def check_status(self, name):
        """Check whether executable is on PATH."""
       
        if isinstance(name, tuple) or isinstance(name, list):
            executable_list = []

            for item in name:
                executable_temp = find_executable(item) is not None
                executable_list.append(executable_temp)
            
            # Check True values
            executable = np.asarray(executable_list)
            True_values = executable[np.where(executable==True)]
            
            if len(True_values) > 1:
                raise ValueError("There are more than one instances installed. Define with one you want to use with self.set_path(...)")
            
            else:  
                return np.any(np.asarray(executable_list) == True)
        
        else:
            return find_executable(name) is not None
        
    def get_path(self, name):
        status = self.check_status(name)
        
        if status:
            executable_list = []
            if isinstance(name, tuple) or isinstance(name, list):
                for item in name:
                    executable_temp = find_executable(item) is not None
                    executable_list.append(executable_temp)
                
                return os.path.abspath(find_executable(name[np.where(np.asarray(executable_list) == True)[0][0]]))
            
            else:
                return os.path.abspath(find_executable(name))             