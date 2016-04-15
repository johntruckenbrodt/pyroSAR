#!/usr/bin/env python2.7

# Copyright (c) 2015, Stefan Engelhardt
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL STEFAN ENGELHARDT BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


"""Dynamic SAR data import system for Gamma.

This module provides a dynamic import system for the Gamma Remote Sensing software. It automatically recognizes the type
of SAR data and calls the correct import tools (par_*). The recognition is based on naming conventions of the data
providers.

The implementation consists of a registry class and a data handler class for each input data format. The handlers are
providing a regular expression for their format and are registered at the registry. Based on the regular expression, the
registry class determines which data handler to use for the current input data. For the actual import each handler has
to have a method 'process' which is called by the registry.

Inheriting from the Handler base class is not mandatory, but recommended to enable automatic searches for handlers based
on Python's introspection capabilities. You can use any custom class which implements to following interface
    
    * string attribute 'mission' containing a abbreviation of the satellites name or mission
    * string attribute 'pattern' representing the regular expression to identify the data format
    * method 'process' taking two parameters, first is the dataset filename, second is the match object returned by the
      re.match function

Data handlers for several data formats are predefined within this module. Note that every handler requires a specific
file or directory of the SAR product for an unambiguous identification
    
    * TerraSAR-X: XML product annotation file
    * TanDEM-X: XML product annotation file
    * COSMO-SkyMed: HDF5 data file
    * RADARSAT-2: product directory which contains the 'product.xml' annotation
    * Sentinel-1: product directory with the extension .SAFE
    * PALSAR-1: CEOS leader file (LED...)
    * PALSAR-2 CEOS: CEOS leader file (LED...)
    * PALSAR-2 GeoTIFF: product KML file (ALOS2....kml)

The generated data and parameter files should follow the naming conventions hereafter
    
    AAAA_BBBB_CCCCCCCCTDDDDDD_EE
    
    Constituent ID   Constituent Name      Remarks
    ---------------------------------------------------------------------------
    AAAA             satellite name        should always contain the satellite
                                           number in multi-satellite
                                           constellations
    BBBB             variant               e.g. sub-swaths in ScanSAR modes or
                                           channels from TanDEM-X cooperative
                                           modes; ____ if not required
    CCCCCCCC         acquisition date      YYYYMMDD format
    DDDDDD           acquisition time      HHMMSS format
    EE               polarization channel  TxRx polarization (HH, HV, VH or VV)
    
    Processing levels shall be indicated as lower case three-character extensions separated by underscores, e.g.
    TSX1______20090405T054534_HH_slc_cal for calibrated single-look complex data. Used identifiers in this script are
    slc, mli and geo.
    
    Filename extension for Gamma ISP parameter files shall be .par and .dem_par for DIFF&GEO parameter files.

If used as script this module will try to import all data sets given as command-line parameters. For more details see
the main function in this module. Additionally there is a functionality to simplify the script usage within larger
software frameworks. This means on the one hand the possibility to reduce the enabled mission per run-normally all
available handlers are enabled. Mit the -m option you can restrict the set of enabled missions. A list of all defined
missions can be obtained by using the -l flag. The flag -s can be used to keep the output free of error messages due to
file not recognized by any handler, which can be used to recursively search directories for data sets.

This module was tested with Python version 2.7 only.
"""


# The table hereafter represents tests, applied to assess the software's
# functionality with real data sets.
#
# Product                                    Result
# -----------------------------------------------------------------------------
# TSX1 SM Dual-pol SSC                       passed
# TSX1 SM Single-pol MGD                     passed
# TSX1 SM Dual-pol EEC                       failed, error by par_TX_geo
# TSX1 ST Single-pol SSC                     passed
# TSX1 ST Single-pol MGD                     passed
# TSX1 ST Single-pol EEC                     failed, error by par_TX_geo
# TSX1 SC Single-pol GEC                     failed, error by par_TX_geo
# TSX1 SC Single-pol MGD                     passed
# TDM1 SM bistatic Single-pol CoSSC          passed
# CSK SM PingPong Dual-pol GEC               passed
# CSK SM Himage Single-pol SCS               passed
# CSK SM Himage Single-pol GEC               passed
# CSK WR Single-pol GEC                      failed, error by libhdf5
# CSK HR Single-pol GEC                      failed, error by libhdf5
# S1A IW SLC                                 passed
# S1A IW GRDH                                passed
# S1A SM GRDH                                passed
# RS2 Standard Dual-pol SGF                  passed
# RS2 Fine Quad-pol SLC                      passed
# RS2 Fine Dual-pol SGF                      passed
# RS2 ScanSAR Dual-pol SGF                   passed
# RS2 Ultra-Fine Single-pol SLC              passed
# PSR1 FBS Dual-pol L1.1                     passed
# PSR1 FBS Dual-pol L1.5                     passed
# PSR1 PLR Quad-pol L1.1                     passed
# PSR1 PLR Quad-pol L1.5                     passed
# PSR1 WB1 Single-pol L1.5                   passed
# PSR2 UBD Dual-pol L1.5 CEOS                passed
# PSR2 FBD Dual-pol L1.5 GeoTIFF             passed, desired error triggered
# PSR2 HBQ Quad-pol L1.5 GeoTIFF             passed, desired error triggered


from __future__ import print_function

import argparse
import glob
import os
import os.path
import random
import re
import string
import struct
import subprocess as sp
import sys
import warnings
from ancillary import finder


__all__ = ['HanderRegistry', 'Handler', 'TerraSARXL1BHandler', 'TanDEMXCoSSCHandler', 'CosmoSkyMedHandler',
           'Radarsat2Handler', 'Sentinel1Handler', 'EORCPalsar1Handler', 'EORCPalsar2CEOSHandler',
           'EORCPalsar2TIFFHandler']


def call_command(*args):
    """Spawns a new process and blocks until it returns."""
    args = [str(i) for i in args]
    sp.check_call(args)


class HandlerError(Exception):
    """Exception class indicating a generic error in a handler.
    
    Attributes:
        handler: Handler instance which triggered the exception.
    """
    
    def __init__(self, handler, message):
        Exception.__init__(self, '{0}: {1}'.format(type(handler).__name__, message))
        self.handler = handler


class MissingHandlerError(Exception):
    """Exception class indicating a missing data handler.
    
    Attributes:
        filename: Name of the data set which triggered the exception.
    """
    
    def __init__(self, filename):
        Exception.__init__(self, '{0}: not recognised as a supported file format'.format(filename))
        self.filename = filename


class Handler(object):
    """Base class for data handlers.
    
    The purpose of this class is to simplify the automatic detection of handler classes by inheriting from a common
    base class.
    """
    
    def process(self, filename, match):
        raise NotImplementedError('process method not implemented')


class CosmoSkyMedHandler(Handler):
    """Handler class for COSMO-SkyMed data.

    References
        ASI-CSM-ENG-RS-092-A  COSMO SkyMed SAR Products Handbook
    """

    mission = 'CSK'
    pattern = r'^CSKS(?P<sat>[1234])_(?P<type>RAW_B|SCS_B|SCS_U|DGM_B|GEC_B|GTC_B)_(?P<mode>HI|PP|WR|HR|S2)_(?P<swath>[A-Z0-9]{2})_(?P<pol>HH|VV|HV|VH|CO|CH|CV)_(?:L|R)(?P<orbit_dir>A|D)_(?:F|S)(?:N|F)_(?P<start>[0-9]{14})_(?P<stop>[0-9]{14})\.h5$'

    mapping_geo = [('{0}_{1}_{2}.geo', 'CSK{3}______{4}_{1}_slc_mli_geo'),
                   ('{0}_{1}_{2}.mli.par', 'CSK{3}______{4}_{1}_slc_mli_geo.par'),
                   ('{0}_{1}_{2}.dem_par', 'CSK{3}______{4}_{1}_slc_mli_geo.dem_par')]

    mapping_slc = [('{0}_{1}_{2}.slc', 'CSK{3}______{4}_{1}_slc'),
                   ('{0}_{1}_{2}.slc.par', 'CSK{3}______{4}_{1}_slc.par')]

    def process(self, filename, match):

        def random_string(length):
            return ''.join(random.choice(string.lowercase) for _ in range(length))

        # Determine required polarizations
        if match.group('pol') == 'CO':
            polarizations = ['HH', 'VV']
        elif match.group('pol') == 'CH':
            polarizations = ['HH', 'HV']
        elif match.group('pol') == 'CV':
            polarizations = ['VV', 'VH']
        else:
            polarizations = [match.group('pol')]  # Single-pol modes
        temp = random_string(6)
        move_command = 'move' if os.name == 'nt' else 'mv'
        if match.group('type').startswith('SCS'):
            call_command('par_CS_SLC', filename, temp)
            mapping = CosmoSkyMedHandler.mapping_slc
        elif match.group('type') == 'GEC_B' or match.group('type') == 'GTC_B':
            call_command('par_CS_geo', filename, temp)
            mapping = CosmoSkyMedHandler.mapping_geo
        else:
            raise HandlerError(self, 'CSK {0} products are not supported'.format(match.group('type')))
        # Rename all files to the generic naming conventions
        for pol in polarizations:
            fields = (temp, pol, match.group('swath'), match.group('sat'), match.group('start')[:8] + 'T' + match.group('start')[8:])
            for src, dst in mapping:
                call_command(move_command, src.format(*fields), dst.format(*fields))


class EORCPalsar1Handler(Handler):
    """Handler class for ALOS PALSAR-1 data from EORC/JAXA.

    References
        NEB-01006    ALOS PALSAR Level 1 Product Format Description Vol. 1
        NEB-070062B  ALOS PALSAR Level 1 Product Format Description Vol. 2
    """

    mission = 'PSR1'
    pattern = r'^LED-ALPSR(?P<sub>P|S)(?P<orbit>[0-9]{5})(?P<frame>[0-9]{4})-(?P<mode>[HWDPC])(?P<level>1\.[015])(?P<proc>G|_)(?P<proj>[UPML_])(?P<orbit_dir>A|D)$'

    @staticmethod
    def read_ceos_datetime(filename):
        """Reads the acquisition date from a CEOS leader file.

        The acquisition center time is read from the Data set Summary Record. For more information on the file
        format see NEB-070062B 'ALOS PALSAR Level 1 Product Format Description Vol. 2'

        Args:
            file: CEOS leader file (LED...).

        Returns:
            Acquisition center time in YYYYMMDDTHHMMSS format.
        """
        with open(filename, 'rb') as f:
            data = f.read(802)
        if len(data) != 802:
            raise EOFError('incomplete CEOS leader file')
        constituents = list(struct.unpack('720x68x8c6c', data))
        constituents.insert(8, 'T')
        return ''.join(constituents)

    def process(self, filename, match):

        if match.group('level') == '1.0':
            raise HandlerError(self, 'PALSAR-1 level 1.0 products are not supported')
        pattern_img = re.compile(
            r'^IMG-(?P<pol>HH|HV|VH|VV)-ALPSR(?P<sub>P|S)(?P<orbit>[0-9]{5})(?P<frame>[0-9]{4})-(?P<mode>[HWDPC])(?P<level>1\.[015])(?P<proc>G|_)(?P<map>[UPML_])(?P<orbit_dir>A|D)$')
        for i in os.listdir(os.path.dirname(filename)):
            match_img = pattern_img.match(i)
            if match_img:
                fields = 'PSR1', EORCPalsar1Handler.read_ceos_datetime(filename), match_img.group('pol')
                if match.group('level') == '1.1':
                    name = '{0}______{1}_{2}_slc'.format(*fields)
                    call_command('par_EORC_PALSAR', filename, name + '.par', os.path.join(os.path.dirname(filename), i), name)
                else:
                    name = '{0}______{1}_{2}_slc_mli_geo'.format(*fields)
                    call_command('par_EORC_PALSAR_geo', filename, name + '.par', name + '.dem_par', os.path.join(os.path.dirname(filename), i), name)


class EORCPalsar2CEOSHandler(Handler):
    """Handler class for ALOS PALSAR-2 CEOS data from EORC/JAXA.

    References:
        ALOS-2/PALSAR-2 Level 1.1/1.5/2.1/3.1 CEOS SAR Product Format Description
    """

    mission = 'PSR2'
    pattern = r'^LED-ALOS2(?P<orbit>[0-9]{5})(?P<frame>[0-9]{4})-(?P<date>[0-9]{6})-(?P<mode>SBS|UBS|UBD|HBS|HBD|HBQ|FBS|FBD|FBQ|WBS|WBD|WWS|WWD|VBS|VBD)(?P<look_dir>L|R)(?P<level>1\.0|1\.1|1\.5|2\.1|3\.1)(?P<proc>[GR_])(?P<map>[UPML_])(?P<orbit_dir>A|D)$'

    def process(self, filename, match):
        if not match.group('level') in ['1.1', '1.5', '2.1']:
            raise HandlerError(self, 'PALSAR-2 level {0} products are not supported'.format(match.group('level')))
        if match.group('level') == '2.1':
            warnings.warn('PALSAR-2 level 2.1 product support is experimental', RuntimeWarning)
        pattern_img = re.compile(
            r'^IMG-(?P<pol>HH|HV|VH|VV)-ALOS2(?P<orbit>[0-9]{5})(?P<frame>[0-9]{4})-(?P<date>[0-9]{6})-(?P<mode>SBS|UBS|UBD|HBS|HBD|HBQ|FBS|FBD|FBQ|WBS|WBD|WWS|WWD|VBS|VBD)(?P<look_dir>L|R)(?P<level>1\.0|1\.1|1\.5|2\.1|3\.1)(?P<proc>[GR_])(?P<map>[UPML_])(?P<orbit_dir>A|D)$')
        for i in os.listdir(os.path.dirname(filename)):
            match_img = pattern_img.match(i)
            if match_img:
                fields = 'PSR2', EORCPalsar1Handler.read_ceos_datetime(filename), match_img.group('pol')
                if match_img.group('level') == '1.1':
                    name = '{0}______{1}_{2}_slc'.format(*fields)
                    call_command('par_EORC_PALSAR', filename, name + '.par', os.path.join(os.path.dirname(filename), i), name)
                else:
                    name = '{0}______{1}_{2}_slc_mli_geo'.format(*fields)
                    call_command('par_EORC_PALSAR_geo', filename, name + '.par', name + '.dem_par', os.path.join(os.path.dirname(filename), i), name)


class EORCPalsar2TIFFHandler(Handler):
    """Handler class for ALOS PALSAR-2 GeoTIFF data from EORC/JAXA.

    References:
        ALOS-2/PALSAR-2 Level 1.1/1.5/2.1/3.1 GeoTiff Product Format Description
    """

    mission = 'PSR2'
    pattern = r'^ALOS2(?P<orbit>[0-9]{5})(?P<frame>[0-9]{4})-(?P<date>[0-9]{6})_(?P<mode>SBS|UBS|UBD|HBS|HBD|HBQ|FBS|FBD|FBQ|WBS|WBD|WWS|WWD|VBS|VBD)(?P<look_dir>L|R)(?P<level>1\.0|1\.1|1\.5|2\.1|3\.1)(?P<proc>[GR_])(?P<map>[UPML_])(?P<orbit_dir>A|D)\.kml$'

    def process(self, filename, match):
        raise HandlerError(self, 'PALSAR-2 GeoTIFF products are not supported')


class Radarsat2Handler(Handler):
    """Handler class for RADARSAT-2 data.

    References
        RN-RP-51-2713  RADARSAT-2 Product Format Definition
    """

    mission = 'RS2'
    pattern = r'^(?:RS2|RSAT2)_(?:OK[0-9]+)_(?:PK[0-9]+)_(?:DK[0-9]+)_(?P<beam>[0-9A-Z]+)_(?P<date>[0-9]{8})_(?P<time>[0-9]{6})_(?P<pols>(?:HH_|VV_|HV_|VH_){1,4})(?P<level>SLC|SGX|SGF|SCN|SCW|SSG|SPG)$'

    def process(self, filename, match):
        polarizations = match.group('pols')[:-1].split('_')
        assert len(polarizations) > 0
        if match.group('level') == 'SGF' or match.group('level') == 'SGX':
            for pol in polarizations:
                name = 'RS2______{0}T{1}_{2}_slc_mli'.format(match.group('date'), match.group('time'), pol)
                call_command('par_RSAT2_SG', os.path.join(filename, 'product.xml'), os.path.join(filename, 'lutSigma.xml'),
                             os.path.join(filename, 'imagery_{0}.tif'.format(pol)), pol, name + '.par', name)
        elif match.group('level') == 'SLC':
            for pol in polarizations:
                name = 'RS2______{0}T{1}_{2}_slc'.format(match.group('date'), match.group('time'), pol)

                # scan working system for all available GAMMA versions and iteratively execute the respective commands until data has successfully been imported
                cmd = sp.Popen(["which", "par_RSAT2_SLC"], stdout=sp.PIPE).stdout.read().strip("\n")
                cmd_base = re.sub("GAMMA_SOFTWARE(.*)", "", cmd)
                versions = [int(re.findall("[0-9]{8}", x)[0]) for x in finder(cmd_base, ["GAMMA_SOFTWARE*"], foldermode=2, recursive=False)]
                processed = False
                for version in sorted(versions, reverse=True):
                    try:
                        command = re.sub("[0-9]{8}", str(version), cmd)
                        call_command(command, os.path.join(filename, 'product.xml'), os.path.join(filename, 'lutSigma.xml'),
                                     os.path.join(filename, 'imagery_{0}.tif'.format(pol)), pol, name+'.par', name)
                        processed = True
                        break
                    except sp.CalledProcessError:
                        os.remove(name)
                        os.remove(name+'.par')
                        continue
                if not processed:
                    raise HandlerError(self, 'segmentation fault while executing par_RSAT2_SLC')
        else:
            raise HandlerError(self, 'Radarsat2 {0} products are not supported'.format(match.group('level')))


# class Sentinel1Handler(Handler):
#     """Handler class for Sentinel-1 data.
#
#     References
#         S1-RS-MDA-52-7441  Sentinel-1 Product Specification
#     """
#
#     mission = 'S1'
#     pattern = r'^(?P<sat>S1[AB])_(?P<beam>S1|S2|S3|S4|S5|S6|IW|EW|WV|EN|N1|N2|N3|N4|N5|N6|IM)_(?P<prod>SLC|GRD|OCN)(?:F|H|M|_)_(?:1|2)(?P<class>S|A)(?P<pols>SH|SV|DH|DV|HH|HV|VV|VH)_(?P<start>[0-9]{8}T[0-9]{6})_(?P<stop>[0-9]{8}T[0-9]{6})_(?:[0-9]{6})_(?:[0-9A-F]{6})_(?:[0-9A-F]{4})\.SAFE$'
#
#     def process(self, filename, match):
#         if match.group('prod') == 'OCN':
#             raise HandlerError(self, 'Sentinel-1 OCN products are not supported')
#         if match.group('class') == 'A':
#             raise HandlerError(self, 'Sentinel-1 annotation-only products are not supported')
#         # Find and import all datasets
#         pattern_ds = re.compile(r'^s1[ab]-(?P<swath>s[1-6]|iw[1-3]?|ew[1-5]?|wv[1-2]|n[1-6])-(?P<prod>slc|grd|ocn)-(?P<pol>hh|hv|vv|vh)-(?P<start>[0-9]{8}t[0-9]{6})-(?P<stop>[0-9]{8}t[0-9]{6})-(?:[0-9]{6})-(?:[0-9a-f]{6})-(?P<id>[0-9]{3})\.xml$')
#         for i in os.listdir(os.path.join(filename, 'annotation')):
#             match_ds = pattern_ds.match(i)
#             if match_ds is not None:
#                 annotation_xml = os.path.join(filename, 'annotation', i)
#                 measurement_tiff = os.path.join(filename, 'measurement', i.replace('.xml', '.tiff'))
#                 calibration_xml = os.path.join(filename, 'annotation', 'calibration', 'calibration-' + i)
#                 noise_xml = os.path.join(filename, 'annotation', 'calibration', 'noise-' + i)
#                 fields = (match.group('sat').upper(), match_ds.group('id'), match_ds.group('start').upper(), match_ds.group('pol').upper())
#                 if match.group('prod') == 'SLC':
#                     name = '{0}__{1}__{2}_{3}_slc'.format(*fields)
#                     call_command('par_S1_SLC', measurement_tiff, annotation_xml, calibration_xml, noise_xml, name + '.par', name, name + '.tops_par')
#                 else:
#                     name = '{0}__{1}__{2}_{3}_slc_mli'.format(*fields)
#                     call_command('par_S1_GRD', measurement_tiff, annotation_xml, calibration_xml, noise_xml, name + '.par', name)


class Sentinel1Handler(Handler):
    """Handler class for Sentinel-1 data.

    References
        S1-RS-MDA-52-7441  Sentinel-1 Product Specification
    """

    mission = 'S1'
    pattern = r'^(?P<sat>S1[AB])_(?P<beam>S1|S2|S3|S4|S5|S6|IW|EW|WV|EN|N1|N2|N3|N4|N5|N6|IM)_(?P<prod>SLC|GRD|OCN)(?:F|H|M|_)_(?:1|2)(?P<class>S|A)(?P<pols>SH|SV|DH|DV|HH|HV|VV|VH)_(?P<start>[0-9]{8}T[0-9]{6})_(?P<stop>[0-9]{8}T[0-9]{6})_(?:[0-9]{6})_(?:[0-9A-F]{6})_(?:[0-9A-F]{4})\.SAFE$'

    def process(self, filename, match):
        if match.group('prod') == 'OCN':
            raise HandlerError(self, 'Sentinel-1 OCN products are not supported')
        if match.group('class') == 'A':
            raise HandlerError(self, 'Sentinel-1 annotation-only products are not supported')
        # Find and import all datasets
        pattern_ds = re.compile(r'^s1[ab]-(?P<swath>s[1-6]|iw[1-3]?|ew[1-5]?|wv[1-2]|n[1-6])-(?P<prod>slc|grd|ocn)-(?P<pol>hh|hv|vv|vh)-(?P<start>[0-9]{8}t[0-9]{6})-(?P<stop>[0-9]{8}t[0-9]{6})-(?:[0-9]{6})-(?:[0-9a-f]{6})-(?P<id>[0-9]{3})\.xml$')
        for i in os.listdir(os.path.join(filename, 'annotation')):
            match_ds = pattern_ds.match(i)
            if match_ds is not None:
                annotation_xml = os.path.join(filename, 'annotation', i)
                measurement_tiff = os.path.join(filename, 'measurement', i.replace('.xml', '.tiff'))
                calibration_xml = os.path.join(filename, 'annotation', 'calibration', 'calibration-' + i)
                # the use of the noise xml file has been found to occasionally cause severe image artifacts of manifold nature and is thus excluded
                # the reason (GAMMA command error vs. bad ESA xml file entry) is yet to be discovered
                # noise_xml = os.path.join(filename, 'annotation', 'calibration', 'noise-' + i)
                noise_xml = "-"
                fields = (match.group('sat').upper(), match_ds.group('swath'), match_ds.group('start').upper(), match_ds.group('pol').upper())
                if match.group('prod') == 'SLC':
                    name = '{0}_{1}__{2}_{3}_slc'.format(*fields)
                    call_command('par_S1_SLC', measurement_tiff, annotation_xml, calibration_xml, noise_xml, name + '.par', name, name + '.tops_par')
                else:
                    name = '{0}______{2}_{3}_mli'.format(*fields)
                    call_command('par_S1_GRD', measurement_tiff, annotation_xml, calibration_xml, noise_xml, name + '.par', name)


class TanDEMXCoSSCHandler(Handler):
    """Handler class for experimental TanDEM-X CoSSC data.

    References:
        TD-GS-PS-3028  TanDEM-X Experimental Product Description
    """

    mission = 'TDX'
    pattern = r'^(?P<mission>TDM[0-9])_SAR__COS_(?P<coop>MONO|BIST|ALT1|ALT2|NONE)_(?P<mode>SM|SL|SL)_(?P<pols>[SDQ])_(?:SRA|DRA)_(?P<start>[0-9]{8}T[0-9]{6})_(?P<stop>[0-9]{8}T[0-9]{6})\.xml$'

    def process(self, filename, match):
        # Since the CoSSC products are essentially two L1B products combined, the existing handler can be used.
        pattern_l1b = re.compile(TerraSARXL1BHandler.pattern)
        handler_l1b = TerraSARXL1BHandler()
        for i in glob.glob(os.path.join(os.path.dirname(filename), '*', '*.xml')):
            match_l1b = pattern_l1b.match(os.path.basename(i))
            if match_l1b:
                handler_l1b.process(i, match_l1b)


class TerraSARXL1BHandler(Handler):
    """Handler class for TerraSAR-X data.
    
    References:
        TX-GS-DD-3302  TerraSAR-X Basic Product Specification Document
        TX-GS-DD-3303  TerraSAR-X Experimental Product Description
        TD-GS-PS-3028  TanDEM-X Experimental Product Description
    """
    
    mission = 'TSX'
    pattern = r'^(?P<sat>T[DS]X1)_SAR__(?P<prod>SSC|MGD|GEC|EEC)_(?P<var>____|SE__|RE__|MON1|MON2|BTX1|BRX2)_(?P<mode>SM|SL|HS|ST|SC)_(?P<pols>[SDTQ])_(?:SRA|DRA)_(?P<start>[0-9]{8}T[0-9]{6})_(?P<stop>[0-9]{8}T[0-9]{6})\.xml$'
    
    def process(self, filename, match):
        # Issue warning for experimental products
        if match.group('pols') == 'T' or match.group('pols') == 'Q':
            warnings.warn('TSX twin- and quad-pol mode are designated experimental', RuntimeWarning)
        # Find all image files and import them with par_TX commands
        image_folder = os.path.join(os.path.dirname(filename), 'IMAGEDATA')
        pattern_img = re.compile(r'^IMAGE_(?P<pol>HH|HV|VH|VV)_(?:SRA|FWD|AFT)_(?P<beam>[^\.]+)\.(cos|tif)$')
        for i in os.listdir(image_folder):
            match_img = pattern_img.match(i)
            if match_img:
                fields = match.group('sat'), match.group('var'), match.group('start'), match_img.group('pol')
                if match.group('prod') == 'SSC':
                    name = '{0}_{1}_{2}_{3}_slc'.format(*fields)
                    call_command('par_TX_SLC', filename, os.path.join(image_folder, i), name + '.par', name, match_img.group('pol'))
                elif match.group('prod') == 'MGD':
                    name = '{0}_{1}_{2}_{3}_slc_mli'.format(*fields)
                    call_command('par_TX_GRD', filename, os.path.join(image_folder, i), name + '.par', name, match_img.group('pol'))
                else:
                    name = '{0}_{1}_{2}_{3}_slc_mli_geo'.format(*fields)
                    call_command('par_TX_geo', filename, os.path.join(image_folder, i), name + '.par', name + '.dem_par', name, match_img.group('pol'))


class HandlerRegistry(object):
    """Generic interface for using handlers.
    
    Attributes:
        handlers: List of registered handlers.
    """
    
    handlers = property(lambda self: self._handlers.values(), doc="""List of registered handlers.""")
    
    def __init__(self, predefined=False):
        """Creates a new registry.
        
        Args:
            predefined: If True, add all data handler defined in the current
                namespace. Note that this will find only those handlers
                inheriting from the base class Handler. Defaults to False.
        """
        self._handlers = {}
        # Automatically search for predefined handlers based on their common super class Handler. This requires
        # new-style classes to work.
        if predefined:
            for handler_class in Handler.__subclasses__():
                try:
                    handler_instance = handler_class()
                    self.register_handler(handler_instance.pattern, handler_instance)
                except AttributeError:
                    warnings.warn('skipped handler {0} because of missing pattern definition'.format(handler_class.__name__), RuntimeWarning)

    def register_handler(self, pattern, handler):
        """Adds a new data handler to the registry.
        
        Args:
            pattern: Regular expression.
            handler: Instance of the handler class.
        
        Raises:
            TypeError: Handler instance is None.
        """
        if handler is None:
            raise TypeError('unexpected NoneType')
        self._handlers[re.compile(pattern)] = handler
    
    def handle(self, file):
        """Determines the handler and imports a data set.
        
        Args:
            file: Filename of the data set. For more information on expected
                names see the documentation of the module.
        
        Raises:
            MissingHandlerError: No matching handler for file found.
            HandlerError: A handler produced an error.
        """
        # Remove trailing path separators to force base name to assume it
        # represents a file rather than a directory. Otherwise basename would
        # return a zero-length string for directories, since it does not use
        # the actual file system.
        basename = os.path.basename(re.sub(r'[/\\]+$', '', file))
        for pattern, handler in self._handlers.items():
            match = pattern.match(basename)
            if match:
                handler.process(file, match)
                return  # Stop after the first pattern matches
        raise MissingHandlerError(file)


def main():
    parser = argparse.ArgumentParser(
        description=('This program provides a dynamic import system for the Gamma Remote Sensing '
                     'software. It automatically recognizes the type of SAR data and calls the '
                     'corresponding import tools (par_*). The recognition is based on naming conventions '
                     'of the data providers.'))
    parser.add_argument('-m', '--mission', nargs=1, metavar='ABBR', action='append', help='abbreviations of mission to be enabled; by default all missions are enabled')
    parser.add_argument('-l', '--list', action='store_true', help='list all available missions and their abbreviations usable with this program')
    parser.add_argument('-s', '--skip', action='store_true', help='if set the program skips files which do not match any handler without any error message')
    parser.add_argument('dataset', nargs='+', help='SAR dataset to be read')
    args = parser.parse_args()
    # List all handlers
    if args.list:
        registry = HandlerRegistry(True)
        # Print list of known handlers, using only one handler per mission
        for mission in set([handler.mission for handler in registry.handlers]):
            print(mission)
        sys.exit()
    # Register requested handlers
    if args.mission is None:
        registry = HandlerRegistry(True)
    else:
        registry = HandlerRegistry(False)
        for mission in args.mission:
            for handler_class in Handler.__subclasses__():
                if mission[0] == handler_class.mission:
                    registry.register_handler(handler_class.pattern, handler_class())
    # print('Enabled mission handler: ' + ', '.join(set([handler.mission for handler in registry.handlers])))
    # Read each passed dataset
    for dataset in args.dataset:
        try:
            registry.handle(dataset)
        except MissingHandlerError:
            if not args.skip:
                print('{0}: {1}: not recognised as a supported file format'.format(os.path.basename(sys.argv[0]), dataset), file=sys.stderr)
        except HandlerError, e:
            print('{0}: {1}: {2}'.format(os.path.basename(sys.argv[0]), dataset, str(e)), file=sys.stderr)

if __name__ == '__main__':
    main()
