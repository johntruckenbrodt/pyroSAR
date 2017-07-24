##############################################################
# ENVI header management
# John Truckenbrodt 2015-2017
##############################################################
"""
This script offers functionality for editing ENVI header files
The object HDRobject can be initialized from an ENVI header file or a GAMMA parameter file
writing of files is possible by
- calling the function hdr with either an HDR object, a GAMMA parameter file or an ENVI header file or
- calling the class function write on an existing HDRobject
example (adding band names to an existing hdr file):
obj = HDRobject('E:/test.hdr')
obj.band_names = ['one', 'two']
obj.write()
"""
import os
import re
from .ancillary import union, parse_literal


def hdr(data, filename='same'):
    """
    write ENVI header files
    """
    hdrobj = data if isinstance(data, HDRobject) else HDRobject(data)
    hdrobj.write(filename)


# todo: check whether this function is of benefit
# http://gis.stackexchange.com/questions/48618/how-to-read-write-envi-metadata-using-gdal
def get_envi_header_dict(hdrfile):
    with open(hdrfile, 'r') as infile:
        hdr = infile.read()
    # Get all 'key = {val}' type matches
    regex = re.compile(r'^(.+?)\s*=\s*({\s*.*?\n*.*?})$', re.M | re.I)
    matches = regex.findall(hdr)

    # Remove them from the header
    subhdr = regex.sub('', hdr)

    # Get all 'key = val' type matches
    regex = re.compile(r'^(.+?)\s*=\s*(.*?)$', re.M | re.I)
    matches.extend(regex.findall(subhdr))

    return dict(matches)


class HDRobject(object):
    """
    create ENVI hdr file object from existing .par or .hdr file
    for creating new headers from .par files currently only EQA and UTM projections with WGS-84 ellipsoid are supported
    """
    def __init__(self, parfile='None'):
        self.filename = 'None' if parfile == 'None' else parfile
        if re.search('.hdr$', parfile):
            with open(parfile, 'r') as infile:
                lines = infile.readlines()
                i = 0
                while i < len(lines):
                    line = lines[i].strip('\r\n')
                    if '=' in line:
                        if '{' in line and '}' not in line:
                            while '}' not in line:
                                i += 1
                                line += lines[i].strip('\n').lstrip()
                        line = filter(None, re.split('\s+=\s+', line))
                        line[1] = re.split(',[ ]*', line[1].strip('{}'))
                        setattr(self, line[0].replace(' ', '_'), line[1] if len(line[1]) > 1 else line[1][0])
                    i += 1
            if type(self.band_names) == str:
                self.band_names = [self.band_names]
        else:
            args = {'bands': 1,
                    'header_offset': 0,
                    'file_type': 'ENVI Standard',
                    'interleave': 'bsq',
                    'sensor_type': 'Unknown',
                    'byte_order': 1,
                    'wavelength_units': 'Unknown'}
            if parfile != 'None':
                par = ISPPar(parfile)
                self.samples = getattr(par, union(['width', 'range_samples', 'samples'], par.__dict__.keys())[0])
                self.lines = getattr(par, union(['nlines', 'azimuth_lines', 'lines'], par.__dict__.keys())[0])
                for arg in args:
                    setattr(self, arg, args[arg])
                # todo: is this all really correct?
                # 2x16 bit complex data (SCOMPLEX) is not supported by ENVI, thus float (1x32 bit) is written to the header file
                dtypes = {'FCOMPLEX': 6, 'SCOMPLEX': 4, 'FLOAT': 4, 'REAL*4': 4, 'INTEGER*2': 2, 'SHORT': 12}
                self.data_type = dtypes[getattr(par, union(['data_format', 'image_format'], par.__dict__.keys())[0])]
                if self.data_type == 6:
                    self.complex_function = 'Power'
                # projections = ['AEAC', 'EQA', 'LCC', 'LCC2', 'OMCH', 'PC', 'PS', 'SCH', 'TM', 'UTM']
                if hasattr(par, 'DEM_projection'):
                    if par.DEM_projection == 'UTM':
                        hem = 'North' if float(par.false_northing) == 0 else 'South'
                        self.map_info = ['UTM', '1.0000', '1.0000',
                                         par.corner_east, par.corner_north,
                                         str(abs(float(par.post_east))), str(abs(float(par.post_north))),
                                         par.projection_zone, hem, 'WGS-84', 'units=Meters']
                    elif par.DEM_projection == 'EQA':
                        self.map_info = ['Geographic Lat/Lon', '1.0000', '1.0000',
                                         par.corner_lon, par.corner_lat,
                                         str(abs(float(par.post_lon))), str(abs(float(par.post_lat))),
                                         'WGS-84', 'units=Degrees']
                    else:
                        raise IOError('unsupported projection')
            else:
                self.samples = 0
                self.lines = 0
                for arg in args:
                    setattr(self, arg, args[arg])

    def write(self, filename='same'):
        """
        write object to an ENVI header file
        """
        filename = os.path.splitext(self.filename)[0]+'.hdr' if filename is 'same' else filename
        with open(filename, 'w') as out:
            out.write('ENVI\n')
            for item in ['description', 'samples', 'lines', 'bands', 'header_offset', 'file_type', 'data_type', 'interleave', 'sensor_type', 'byte_order', 'map_info',
                         'coordinate_system_string', 'wavelength_units', 'band_names']:
                if hasattr(self, item):
                    value = getattr(self, item)
                    if isinstance(value, list):
                        out.write(item.replace('_', ' ') + ' = {' + ', '.join([str(x) for x in value]) + '}\n')
                    elif item in ['description', 'band_names', 'coordinate_system_string']:
                        out.write(item.replace('_', ' ') + ' = {' + value + '}\n')
                    else:
                        out.write(item.replace('_', ' ') + ' = ' + str(value) + '\n')


class ISPPar(object):
    """Reader for ISP parameter files of the GAMMA software package.

    This class allows to read all information from filed in GAMMA's parameter file format. Each key-value pair is parsed
    and added as attributes. For instance if the parameter file contains the pair 'sensor:    TSX-1' an attribute named
    'sensor' with the value 'TSX-1' will be available.

    The values are converted to native Python types, while unit identifiers like 'dB' or 'Hz' are removed. Please see
    the GAMMA reference manual for further information on the actual file format.
    """

    _re_kv_pair = re.compile(r'^(\w+):\s*(.+)\s*')
    _re_float_literal = re.compile(r'^[+-]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?')

    def __init__(self, filename):
        """Parses a ISP parameter file from disk.

        Args:
            filename: The filename or file object representing the ISP parameter file.
        """
        if isinstance(filename, basestring):
            par_file = open(filename, 'r')
        else:
            par_file = filename
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
                setattr(self, key, value)
        finally:
            par_file.close()
