##############################################################
# general GAMMA utilities
# Stefan Engelhardt, John Truckenbrodt 2014-2018
##############################################################
import math
import os
import re
import json
import subprocess as sp
from datetime import datetime

from spatialist.ancillary import parse_literal, run, union, dissolve
from spatialist.envi import hdr

from pyroSAR import ConfigHandler
from .error import gammaErrorHandler


class ISPPar(object):
    """
    Reader for ISP parameter files of the GAMMA software package

    This class allows to read all information from files in GAMMA's parameter file format.
    Each key-value pair is parsed and added as attribute. For instance if the parameter file
    contains the pair 'sensor:    TSX-1' an attribute named 'sensor' with the value 'TSX-1' will be available.

    The values are converted to native Python types, while unit identifiers like 'dB' or 'Hz' are removed.
    Please see the GAMMA reference manual for further information on the actual file format.
    
    Parameters
    ----------
    filename: str
        the Gamma parameter file
    
    Examples
    --------
    >>> from pyroSAR.gamma import ISPPar
    >>> with ISPPar('S1A__IW___A_20141115T181801_VH_grd.par') as par:
    ...     print(par) # print an overview of all available metadata
    ...     print(par.keys) # print all parameter names
    ...     for key, value in par.envidict().items():
    ...         print('{0}: {1}'.format(key, value)) # print the ENVI HDR compliant metadata
    
    Attributes
    ----------
    keys : list
        the names of all parameters
    """
    
    _re_kv_pair = re.compile(r'^(\w+):\s*(.+)\s*')
    _re_float_literal = re.compile(r'^[+-]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?')
    
    def __init__(self, filename):
        """Parses an ISP parameter file from disk.

        Args:
            filename: The filename or file object representing the ISP parameter file.
        """
        if isinstance(filename, str):
            par_file = open(filename, 'r')
        else:
            par_file = filename
        
        self.keys = []
        
        try:
            par_file.readline()  # Skip header line
            for line in par_file:
                match = ISPPar._re_kv_pair.match(line)
                if not match:
                    continue  # Skip malformed lines with no key-value pairs
                key = match.group(1)
                items = match.group(2).split()
                if len(items) == 0:
                    value = None
                elif len(items) == 1:
                    value = parse_literal(items[0])
                else:
                    if not ISPPar._re_float_literal.match(items[0]):
                        # Value is a string literal containing whitespace characters
                        value = match.group(2)
                    else:
                        # Evaluate each item and stop at the first non-float literal
                        value = []
                        for i in items:
                            match = ISPPar._re_float_literal.match(i)
                            if match:
                                value.append(parse_literal(match.group()))
                            else:
                                # If the first float literal is immediately followed by a non-float literal handle the
                                # first one as singular value, e.g. in '20.0970 dB'
                                if len(value) == 1:
                                    value = value[0]
                                break
                self.keys.append(key)
                setattr(self, key, value)
        finally:
            par_file.close()
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        return
    
    def __str__(self):
        maxlen = len(max(self.keys, key=len)) + 1
        return '\n'.join(['{key}:{sep}{value}'.format(key=key,
                                                      sep=(maxlen - len(key)) * ' ',
                                                      value=getattr(self, key)) for key in self.keys])
    
    def envidict(self):
        """
        export relevant metadata to a ENVI HDR file compliant format
        
        Returns
        -------
        dict
        """
        out = dict(bands=1,
                   header_offset=0,
                   file_type='ENVI Standard',
                   interleave='bsq',
                   sensor_type='Unknown',
                   byte_order=1,
                   wavelength_units='Unknown')
        
        out['samples'] = getattr(self, union(['width', 'range_samples', 'samples'], self.keys)[0])
        out['lines'] = getattr(self, union(['nlines', 'azimuth_lines', 'lines'], self.keys)[0])
        
        dtypes_lookup = {'FCOMPLEX': 6, 'FLOAT': 4, 'REAL*4': 4, 'INTEGER*2': 2, 'SHORT': 12}
        dtype = getattr(self, union(['data_format', 'image_format'], self.keys)[0])
        
        if dtype == 'SCOMPLEX':
            raise TypeError('unsupported data type: SCOMPLEX (2x16 bit complex)')
        
        out['data_type'] = dtypes_lookup[dtype]
        
        if out['data_type'] == 6:
            out['complex_function'] = 'Power'
        # projections = ['AEAC', 'EQA', 'LCC', 'LCC2', 'OMCH', 'PC', 'PS', 'SCH', 'TM', 'UTM']
        if hasattr(self, 'DEM_projection'):
            if self.DEM_projection == 'UTM':
                hem = 'North' if float(self.false_northing) == 0 else 'South'
                out['map_info'] = ['UTM', '1.0000', '1.0000',
                                   self.corner_east, self.corner_north,
                                   str(abs(float(self.post_east))), str(abs(float(self.post_north))),
                                   self.projection_zone, hem, 'WGS-84', 'units=Meters']
            elif self.DEM_projection == 'EQA':
                out['map_info'] = ['Geographic Lat/Lon', '1.0000', '1.0000',
                                   self.corner_lon, self.corner_lat,
                                   str(abs(float(self.post_lon))), str(abs(float(self.post_lat))),
                                   'WGS-84', 'units=Degrees']
            else:
                raise RuntimeError('unsupported projection')
        return out


def par2hdr(parfile, hdrfile):
    """
    Create an ENVI HDR file from a Gamma PAR file
    
    Parameters
    ----------
    parfile: str
        the Gamma parfile
    hdrfile: str
        the ENVI HDR file

    Returns
    -------

    """
    with ISPPar(parfile) as par:
        hdr(par.envidict(), hdrfile)


class UTM(object):
    """
    convert a gamma parameter file corner coordinate from EQA to UTM
    """
    
    def __init__(self, parfile):
        par = ISPPar(parfile)
        inlist = [str(x) for x in ['WGS84', 1, 'EQA', par.corner_lon, par.corner_lat, '', 'WGS84', 1, 'UTM', '']]
        proc = sp.Popen(['coord_trans'], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True,
                        shell=False).communicate(''.join([x + '\n' for x in inlist]))
        proc = [x for x in filter(None, proc[0].split('\n')) if ':' in x]
        self.index = []
        for item in proc:
            entry = item.split(': ')
            entry = [entry[0].replace(' ', '_'), entry[1].split()]
            if len(entry[1]) > 1:
                setattr(self, entry[0], entry[1])
            else:
                setattr(self, entry[0], entry[1][0])
            self.index.append(entry[0])
            if 'UTM' in entry[0]:
                self.zone, self.northing, self.easting = entry[1]
                self.index = list(set(self.index + ['zone', 'northing', 'easting']))


def process(cmd, outdir=None, logfile=None, logpath=None, inlist=None, void=True, shellscript=None):
    """
    wrapper function to execute GAMMA commands via module :mod:`subprocess`
    
    Parameters
    ----------
    cmd: list
        the command line arguments
    outdir: str
        the directory to execute the command in
    logfile: str
        a file to write the command log to; overrides parameter logpath
    logpath: str
        a directory to write logfiles to; the file will be named {GAMMA command}.log, e.g. gc_map.log;
        is overridden by parameter logfile
    inlist: list
        a list of values, which is passed as interactive inputs via stdin
    void: bool
        return the stdout and stderr messages?
    shellscript: str
        a file to write the Gamma commands to in shell format
    
    Returns
    -------
    tuple of str or None
        the stdout and stderr messages if void is False, otherwise None
    """
    if logfile is not None:
        log = logfile
    else:
        log = os.path.join(logpath, os.path.basename(cmd[0]) + '.log') if logpath else None
    if shellscript is not None:
        line = ' '.join([str(x) for x in dissolve(cmd)])
        if inlist is not None:
            line += ' <<< $"{}"'.format('\n'.join([str(x) for x in inlist]) + '\n')
        with open(shellscript, 'a+') as sh:
            if outdir is not None:
                first = sh.read(1)
                if not first:
                    ts = datetime.now().strftime('%a %b %d %H:%M:%S %Y')
                    sh.write('# this script was created automatically by pyroSAR on {}\n\n'.format(ts))
                    sh.write('export base={}\n\n'.format(outdir))
                line = line.replace(outdir, '$base')
            sh.write(line + '\n\n')
    out, err = run(cmd, outdir=outdir, logfile=log, inlist=inlist, void=False, errorpass=True)
    gammaErrorHandler(out, err)
    if not void:
        return out, err


class Spacing(object):
    def __init__(self, par, targetres='automatic'):
        """
        compute ground multilooking factors and pixel spacings from an ISPPar object for a defined target resolution
        """
        # compute ground range pixel spacing
        par = par if isinstance(par, ISPPar) else ISPPar(par)
        self.groundRangePS = par.range_pixel_spacing / (math.sin(math.radians(par.incidence_angle)))
        # compute initial multilooking factors
        if targetres == 'automatic':
            if self.groundRangePS > par.azimuth_pixel_spacing:
                ratio = self.groundRangePS / par.azimuth_pixel_spacing
                self.rlks = 1
                self.azlks = int(round(ratio))
            else:
                ratio = par.azimuth_pixel_spacing / self.groundRangePS
                self.rlks = int(round(ratio))
                self.azlks = 1
        else:
            self.rlks = int(round(float(targetres) / self.groundRangePS))
            self.azlks = int(round(float(targetres) / par.azimuth_pixel_spacing))


class Namespace(object):
    def __init__(self, directory, basename):
        self.__base = basename
        self.__outdir = directory
        self.__reg = []
    
    def appreciate(self, keys):
        for key in keys:
            setattr(self, key.replace('.', '_'), os.path.join(self.__outdir, self.__base + '_' + key))
            if key not in self.__reg:
                self.__reg.append(key.replace('.', '_'))
    
    def depreciate(self, keys):
        for key in keys:
            setattr(self, key.replace('.', '_'), '-')
            if key not in self.__reg:
                self.__reg.append(key.replace('.', '_'))
    
    def getall(self):
        out = {}
        for key in self.__reg:
            out[key] = getattr(self, key)
        return out
    
    def select(self, selection):
        return [getattr(self, key) for key in selection]
    
    def isregistered(self, key):
        return key in self.__reg
    
    def isappreciated(self, key):
        if self.isregistered(key):
            if self.get(key) != '-':
                return True
        return False
    
    def get(self, key):
        return getattr(self, key)


def slc_corners(parfile):
    """
    extract the corner coordinates of a SAR scene
    
    Parameters
    ----------
    parfile: str
        the GAMMA parameter file to read coordinates from

    Returns
    -------
    dict of float
        a dictionary with keys xmin, xmax, ymin, ymax
    """
    out, err = process(['SLC_corners', parfile], void=False)
    pts = {}
    for line in out.split('\n'):
        if line.startswith('min. latitude'):
            pts['ymin'], pts['ymax'] = [float(x) for x in
                                        re.findall('[0-9]+\.[0-9]+', line)]
        elif line.startswith('min. longitude'):
            pts['xmin'], pts['xmax'] = [float(x) for x in
                                        re.findall('[0-9]+\.[0-9]+', line)]
    return pts


class ExamineGamma(object):
    def __init__(self):
        if 'GAMMA' in ConfigHandler.sections:
            attr = ConfigHandler['GAMMA']
            for key, value in attr.items():
                setattr(self, key, value)
        if not hasattr(self, 'home'):
            try:
                setattr(self, 'home', os.environ['GAMMA_HOME'])
            except KeyError:
                raise RuntimeError('could not read Gamma installation directory')
        self.version = re.search('GAMMA_SOFTWARE-(?P<version>[0-9]{8})',
                                 getattr(self, 'home')).group('version')
        self.__update_config()
    
    def __update_config(self):
        if 'GAMMA' not in ConfigHandler.sections:
            ConfigHandler.add_section('GAMMA')
    
        for attr in ['home', 'version']:
            self.__update_config_attr(attr, getattr(self, attr), 'GAMMA')

    @staticmethod
    def __update_config_attr(attr, value, section):
        if isinstance(value, list):
            value = json.dumps(value)
    
        if attr not in ConfigHandler[section].keys() or ConfigHandler[section][attr] != value:
            ConfigHandler.set(section, key=attr, value=value, overwrite=True)
