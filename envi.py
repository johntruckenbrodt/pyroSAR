##############################################################
# ENVI header management
# John Truckenbrodt 2015-2016
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
# todo: remove gamma dependency
import os
import re
from gamma import ISPPar
from ancillary import union


# write ENVI header files
def hdr(data, filename='same'):
    hdrobj = data if isinstance(data, HDRobject) else HDRobject(data)
    hdrobj.write(filename)


# todo: check whether this function is of benefit
# http://gis.stackexchange.com/questions/48618/how-to-read-write-envi-metadata-using-gdal
def get_envi_header_dict(hdr):
    # Get all 'key = {val}' type matches
    regex = re.compile(r'^(.+?)\s*=\s*({\s*.*?\n*.*?})$', re.M | re.I)
    matches = regex.findall(hdr)

    # Remove them from the header
    subhdr = regex.sub('', hdr)

    # Get all 'key = val' type matches
    regex = re.compile(r'^(.+?)\s*=\s*(.*?)$', re.M | re.I)
    matches.extend(regex.findall(subhdr))

    return dict(matches)


# create ENVI hdr file object from existing .par or .hdr file
# for creating new headers from .par files currently only EQA and UTM projections with WGS-84 ellipsoid are supported
class HDRobject(object):
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
                dtypes = {'FCOMPLEX': 6, 'FLOAT': 4, 'REAL*4': 4, 'INTEGER*2': 2, 'SHORT': 12}
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

    # write object to an ENVI header file
    def write(self, filename='same'):
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
