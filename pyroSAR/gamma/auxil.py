##############################################################
# general GAMMA utilities
# Stefan Engelhardt, John Truckenbrodt 2014-2019
##############################################################
import math
import os
import re
import json
import string
import shutil
import codecs
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
            content = par_file.read().split('\n')
        except UnicodeDecodeError:
            par_file = codecs.open(filename, 'r', encoding='utf-8', errors='ignore')
            content = par_file.read()
            printable = set(string.printable)
            content = filter(lambda x: x in printable, content)
            content = ''.join(list(content)).split('\n')
        try:
            for line in content:
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

        if hasattr(self, 'date'):
            self.date = '{}-{:02d}-{:02d}T{:02d}:{:02d}:{:02f}'.format(*self.date)
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        return
    
    def __str__(self):
        maxlen = len(max(self.keys, key=len)) + 1
        return '\n'.join(['{key}:{sep}{value}'.format(key=key,
                                                      sep=(maxlen - len(key)) * ' ',
                                                      value=getattr(self, key)) for key in self.keys])
    
    def envidict(self, nodata=None):
        """
        export relevant metadata to a ENVI HDR file compliant format
        
        Parameters
        ----------
        nodata: int, float or None
            a no data value to write to the HDR file via attribute 'data ignore value'
        
        Returns
        -------
        dict
            a dictionary containing attributes translated to ENVI HDR naming
        """
        out = dict(bands=1,
                   header_offset=0,
                   file_type='ENVI Standard',
                   interleave='bsq',
                   sensor_type='Unknown',
                   byte_order=1,
                   wavelength_units='Unknown')
        
        if hasattr(self, 'date'):
            out['acquisition_time'] = self.date + 'Z'
        
        out['samples'] = getattr(self, union(['width', 'range_samples', 'samples'], self.keys)[0])
        out['lines'] = getattr(self, union(['nlines', 'azimuth_lines', 'lines'], self.keys)[0])
        
        dtypes_lookup = {'FCOMPLEX': 6, 'FLOAT': 4, 'REAL*4': 4, 'INTEGER*2': 2, 'SHORT': 12}
        dtype = getattr(self, union(['data_format', 'image_format'], self.keys)[0])
        
        if dtype not in dtypes_lookup.keys():
            raise TypeError('unsupported data type: {}'.format(dtype))
        
        out['data_type'] = dtypes_lookup[dtype]
        
        if nodata is not None:
            out['data_ignore_value'] = nodata
        
        if out['data_type'] == 6:
            out['complex_function'] = 'Power'
        # projections = ['AEAC', 'EQA', 'LCC', 'LCC2', 'OMCH', 'PC', 'PS', 'SCH', 'TM', 'UTM']
        # the corner coordinates are shifted by 1/2 pixel to the Northwest since GAMMA pixel
        # coordinates are defined for the pixel center while in ENVI it is the upper left
        if hasattr(self, 'DEM_projection'):
            if self.DEM_projection == 'UTM':
                hem = 'North' if float(self.false_northing) == 0 else 'South'
                out['map_info'] = ['UTM', '1.0000', '1.0000',
                                   self.corner_east - (abs(self.post_east) / 2),
                                   self.corner_north + (abs(self.post_north) / 2),
                                   str(abs(float(self.post_east))),
                                   str(abs(float(self.post_north))),
                                   self.projection_zone, hem, 'WGS-84', 'units=Meters']
            elif self.DEM_projection == 'EQA':
                out['map_info'] = ['Geographic Lat/Lon', '1.0000', '1.0000',
                                   self.corner_lon - (abs(self.post_lon) / 2),
                                   self.corner_lat + (abs(self.post_lat) / 2),
                                   str(abs(float(self.post_lon))),
                                   str(abs(float(self.post_lat))),
                                   'WGS-84', 'units=Degrees']
            else:
                raise RuntimeError('unsupported projection: {}'.format(self.DEM_projection))
        return out


def par2hdr(parfile, hdrfile, modifications=None, nodata=None):
    """
    Create an ENVI HDR file from a Gamma PAR file
    
    Parameters
    ----------
    parfile: str
        the Gamma parfile
    hdrfile: str
        the ENVI HDR file
    modifications: dict or None
        a dictionary containing value deviations to write to the HDR file
    nodata: int, float or None
        a no data value to write to the HDR file via attribute 'data ignore value'

    Returns
    -------
    
    Examples
    --------
    >>> from pyroSAR.gamma.auxil import par2hdr
    >>> par2hdr('dem_seg.par', 'inc.hdr')
    # write a HDR file for byte data based on a parfile of float data
    >>> par2hdr('dem_seg.par', 'ls_map.hdr', modifications={'data_type': 1})
    
    See Also
    --------
    :class:`spatialist.envi.HDRobject`
    :func:`spatialist.envi.hdr`
    """
    
    with ISPPar(parfile) as par:
        items = par.envidict(nodata)
        if modifications is not None:
            items.update(modifications)
        hdr(items, hdrfile)


class UTM(object):
    """
    convert a gamma parameter file corner coordinate from EQA to UTM
    
    Parameters
    ----------
    parfile: str
        the Gamma parameter file to read the coordinate from
    
    Example
    -------
    
    >>> from pyroSAR.gamma import UTM
    >>> print(UTM('gamma.par').zone)
    """
    
    def __init__(self, parfile):
        par = ISPPar(parfile)
        inlist = ['WGS84', 1, 'EQA', par.corner_lon, par.corner_lat, '', 'WGS84', 1, 'UTM', '']
        inlist = map(str, inlist)
        proc = sp.Popen(['coord_trans'], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE,
                        universal_newlines=True, shell=False)
        out, err = proc.communicate(''.join([x + '\n' for x in inlist]))
        out = [x for x in filter(None, out.split('\n')) if ':' in x]
        
        self.meta = dict()
        for line in out:
            key, value = re.split(r'\s*:\s*', line)
            value = value.split()
            value = map(parse_literal, value) if len(value) > 1 else value[0]
            self.meta[key] = value
        try:
            self.zone, self.northing, self.easting, self.altitude = \
                self.meta['UTM zone/northing/easting/altitude (m)']
        except KeyError:
            self.zone, self.northing, self.easting = \
                self.meta['UTM zone/northing/easting (m)']


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
        if not os.path.isfile(shellscript):
            with open(shellscript, 'w') as init:
                pass
        line = ' '.join([str(x) for x in dissolve(cmd)])
        if inlist is not None:
            line += ' <<< $"{}"'.format('\n'.join([str(x) for x in inlist]) + '\n')
        with open(shellscript, 'r+') as sh:
            if outdir is not None:
                content = sh.read()
                sh.seek(0)
                is_new = re.search('this script was created automatically by pyroSAR', content) is None
                if is_new:
                    ts = datetime.now().strftime('%a %b %d %H:%M:%S %Y')
                    sh.write('# this script was created automatically by pyroSAR on {}\n\n'.format(ts))
                    sh.write('export base={}\n\n'.format(outdir))
                    sh.write(content)
                line = line.replace(outdir, '$base')
            sh.seek(0, 2)  # set pointer to the end of the file
            sh.write(line + '\n\n')
    
    # create an environment containing the locations of all GAMMA submodules to be passed ot the subprocess calls
    gammaenv = {'GAMMA_HOME': ExamineGamma().home}
    for module in ['DIFF', 'DISP', 'IPTA', 'ISP', 'LAT']:
        loc = os.path.join(gammaenv['GAMMA_HOME'], module)
        if os.path.isdir(loc):
            gammaenv[module + '_HOME'] = loc
    
    # execute the command
    out, err = run(cmd, outdir=outdir, logfile=log, inlist=inlist, void=False, errorpass=True, env=gammaenv)
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
    
    def __getitem__(self, item):
        return getattr(self, item.replace('.', '_'))
    
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
    
    def isfile(self, key):
        return hasattr(self, key) and os.path.isfile(getattr(self, key))
    
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
                                        re.findall(r'[0-9]+\.[0-9]+', line)]
        elif line.startswith('min. longitude'):
            pts['xmin'], pts['xmax'] = [float(x) for x in
                                        re.findall(r'[0-9]+\.[0-9]+', line)]
    return pts


class ExamineGamma(object):
    def __init__(self):
        home_sys = os.environ.get('GAMMA_HOME')
        if 'GAMMA' in ConfigHandler.sections:
            attr = ConfigHandler['GAMMA']
            for key, value in attr.items():
                setattr(self, key, value)
        if hasattr(self, 'home'):
            if home_sys is not None and self.home != home_sys:
                print('the value of GAMMA_HOME is different to that in the pyroSAR configuration;\n'
                      '  was: {}\n'
                      '  is : {}\n'
                      'resetting the configuration and deleting parsed modules'
                      .format(self.home, home_sys))
                parsed = os.path.join(os.path.dirname(ConfigHandler.file), 'gammaparse')
                shutil.rmtree(parsed)
                self.home = home_sys
        if not hasattr(self, 'home'):
            if home_sys is not None:
                setattr(self, 'home', home_sys)
            else:
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
