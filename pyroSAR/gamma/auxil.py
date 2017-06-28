
import math
import os
import re
import subprocess as sp

from ..ancillary import parse_literal


class ISPPar(object):
    """Reader for ISP parameter files of the GAMMA software package.

    This class allows to read all information from filed in GAMMA's parameter file format. Each key-value pair is parsed
    and added as attributes. For instance if the parameter file contains the pair 'sensor:    TSX-1' a attribute named
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


# todo: make pretty or remove gamma dependency
class UTM(object):
    """
    convert a gamma parameter file corner coordinate from EQA to UTM
    """
    def __init__(self, parfile):
        par = ISPPar(parfile)
        inlist = [str(x) for x in ["WGS84", 1, "EQA", par.corner_lon, par.corner_lat, "", "WGS84", 1, "UTM", ""]]
        proc = sp.Popen(["coord_trans"], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True, shell=False).communicate("".join([x + "\n" for x in inlist]))
        proc = [x for x in filter(None, proc[0].split("\n")) if ":" in x]
        self.index = []
        for item in proc:
            entry = item.split(": ")
            entry = [entry[0].replace(" ", "_"), entry[1].split()]
            if len(entry[1]) > 1:
                setattr(self, entry[0], entry[1])
            else:
                setattr(self, entry[0], entry[1][0])
            self.index.append(entry[0])
            if "UTM" in entry[0]:
                self.zone, self.northing, self.easting = entry[1]
                self.index = list(set(self.index+["zone", "northing", "easting"]))


class Spacing(object):
    def __init__(self, par, targetres="automatic"):
        """
        compute ground multilooking factors and pixel spacings from an ISPPar object for a defined target resolution
        """
        # compute ground range pixel spacing
        par = par if isinstance(par, ISPPar) else ISPPar(par)
        self.groundRangePS = par.range_pixel_spacing/(math.sin(math.radians(par.incidence_angle)))
        # compute initial multilooking factors
        if targetres == "automatic":
            if self.groundRangePS > par.azimuth_pixel_spacing:
                ratio = self.groundRangePS/par.azimuth_pixel_spacing
                self.rlks = 1
                self.azlks = int(round(ratio))
            else:
                ratio = par.azimuth_pixel_spacing/self.groundRangePS
                self.rlks = int(round(ratio))
                self.azlks = 1
        else:
            self.rlks = int(round(float(targetres)/self.groundRangePS))
            self.azlks = int(round(float(targetres)/par.azimuth_pixel_spacing))


class Namespace(object):
    def __init__(self, directory, basename):
        self.__base = basename
        self.__outdir = directory
        self.__reg = []

    def appreciate(self, keys):
        for key in keys:
            setattr(self, key.replace('.', '_'), os.path.join(self.__outdir, self.__base+'_'+key))
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
