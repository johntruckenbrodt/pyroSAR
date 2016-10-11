
import os
import re
import ssl
import sys
import time
import spatial
import zipfile as zf
from ancillary import finder, dissolve
from datetime import datetime
from urllib import urlopen, urlencode
from urlparse import urlparse, urlunparse
import xml.etree.ElementTree as ElementTree
import gamma

try:
    import argparse
except ImportError:
    try:
        os.remove(os.path.join(os.path.dirname(sys.argv[0]), "locale.pyc"))
    finally:
        import argparse


# scans zipped or unzipped S1 scene folder names for key acquisition parameters
# satellite, beam, acquisition start and stop, product and orbit
class Identifier(object):
    def __init__(self, scene):

        self.base = os.path.realpath(scene)

        pattern = r"^(?P<sat>S1[AB])_(?P<beam>S1|S2|S3|S4|S5|S6|IW|EW|WV|EN|N1|N2|N3|N4|N5|N6|IM)_(?P<prod>SLC|GRD|OCN)(?:F|H|M|_)_(?:1|2)(?P<class>S|A)(?P<pols>SH|SV|DH|DV|HH|HV|VV|VH)_(?P<start>[0-9]{8}T[0-9]{6})_(?P<stop>[0-9]{8}T[0-9]{6})_(?:[0-9]{6})_(?:[0-9A-F]{6})_(?:[0-9A-F]{4})\.SAFE$"

        self.zipped = zf.is_zipfile(scene)

        if self.zipped:
            with zf.ZipFile(scene, "r") as z:
                self.namelist = z.namelist()
                self.scene = os.path.commonprefix(self.namelist).strip("/")
        else:
            self.scene = os.path.basename(os.path.normpath(scene))

        match = re.match(pattern, self.scene)
        if not match:
            raise IOError("folder does not match S1 scene naming convention")
        for key in re.compile(pattern).groupindex:
            setattr(self, key, match.group(key))

        self.orbit = "D" if float(re.findall("[0-9]{6}", self.start)[1]) < 120000 else "A"

    def getCorners(self):
        if self.zipped:
            with zf.ZipFile(self.base, "r") as z:
                kml = z.open([x for x in z.namelist() if re.search("map-overlay\.kml", x)][0], "r").read()
        else:
            with open(finder(self.base, ["*map-overlay.kml"])[0], "r") as infile:
                kml = infile.read()
        elements = ElementTree.fromstring(kml).findall(".//coordinates")

        coordinates = [x.split(",") for x in elements[0].text.split()]
        lat = [float(x[1]) for x in coordinates]
        lon = [float(x[0]) for x in coordinates]
        return {"xmin": min(lon), "xmax": max(lon), "ymin": min(lat), "ymax": max(lat)}

    def bbox(self, outname=None, format="ESRI Shapefile", overwrite=True):
        geo = self.getCorners()
        crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
        if outname is None:
            return spatial.bbox(geo, crs)
        else:
            spatial.bbox(geo, crs, outname=outname, format=format, overwrite=overwrite)


# initialize argument parser for S1 processing utilities
def init_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--transform", action="store_true", help="transform the final DEM to UTM coordinates")
    parser.add_argument("-l", "--logfiles", action="store_true", help="create logfiles of the executed GAMMA commands")
    parser.add_argument("-i", "--intermediates", action="store_true", help="keep intermediate files")
    parser.add_argument("-q", "--quiet", action="store_true", help="suppress standard console prints")
    parser.add_argument("-tr", "--targetresolution", default=20, help="the target resolution in meters for x and y", type=int)
    parser.add_argument("-fg", "--func_geoback", default=2, help="backward geocoding interpolation function; "
                                                                 "0 - Nearest Neighbor, 1 - Bicubic Spline, 2 - Bicubic Spline-Log; "
                                                                 "method 1: negative values possible (e.g. in urban areas) - use method 2 to avoid this", type=int)
    parser.add_argument("-fi", "--func_interp", default=0, help="function for interpolation of layover/shadow/foreshortening/DEM gaps; "
                                                                "0 - set to 0, 1 - linear interpolation, 2 - actual value, 3 - nn-thinned", type=int)
    parser.add_argument("-poe", "--poedir", default=None, help="directory containing aux_poeorb (precise orbit ephemerides) orbit state vector files")
    parser.add_argument("-res", "--resdir", default=None, help="directory containing aux_resorb (restituted orbit) orbit state vector files")
    parser.add_argument("zipfile", help="S1 zipped scene archive to be used")
    parser.add_argument("tempdir", help="temporary directory for intermediate files")
    parser.add_argument("outdir", help="output directory")
    parser.add_argument("srtmdir", help="directory containing SRTM hgt tiles (subdirectories possible)")
    return parser


# debursting of S1 SLC imagery
# the procedure consists of two steps. First antenna pattern deramping and then mosaicing of the single deramped bursts
# for mosaicing, the burst boundaries are calculated from the number of looks in range (rlks) and azimuth (azlks), in this case 5 range looks and 1 azimuth looks.
# Alternately 10 range looks and 2 azimuth looks could be used.
# if replace is set to True, the original files will be deleted
def deburst(burst1, burst2, burst3, name_out, rlks=5, azlks=1, replace=False, path_log=None):
    for burst in [burst1, burst2, burst3]:
        if not os.path.isfile(burst) or not os.path.isfile(burst+".par") or not os.path.isfile(burst+".tops_par"):
            raise IOError("input files missing; parameter files must be named e.g. {burst1}.par and {burst1}.tops_par")
    outpath = os.path.dirname(name_out)
    if not os.path.isdir(outpath):
        os.makedirs(outpath)
    tab_in = os.path.join(outpath, "tab_deramp1")
    tab_out = os.path.join(outpath, "tab_deramp2")
    with open(tab_in, "w") as out1:
        with open(tab_out, "w") as out2:
            for item in [burst1, burst2, burst3]:
                out1.write(item+"\t"+item+".par\t"+item+".tops_par\n")
                out2.write(item+"_drp\t"+item+"_drp.par\t"+item+"_drp.tops_par\n")
    gamma.process(["SLC_deramp_S1_TOPS", tab_in, tab_out, 0, 0], logpath=path_log)
    gamma.process(["SLC_mosaic_S1_TOPS", tab_out, name_out, name_out + ".par", rlks, azlks], logpath=path_log)

    if replace:
        for item in [burst1, burst2, burst3]:
            for subitem in [item+x for x in ["", ".par", ".tops_par"]]:
                os.remove(subitem)
    for item in [burst1, burst2, burst3]:
        for subitem in [item+x for x in ["_drp", "_drp.par", "_drp.tops_par"]]:
            os.remove(subitem)
    os.remove(tab_in)
    os.remove(tab_out)


# interface for management of S1 Orbit State Vector (OSV) files
# input are two directories, one for Precise Orbit Ephemerides (POE) and one for Restituted Orbit (RES) files; these directories are created if they do not exist
# actions performed upon calling the main function 'update':
# -the ESA Quality Control (QC) server is checked for any POE files not in the local directory
# -POE  files on the server and not in the local directory are downloaded
# -RES files newer than the latest POE file are downloaded; POE files are approximately 18 days behind the actual date, thus RES files can be used instead
# -delete all RES files for whose date a POE file has become available
# using function 'match' the corresponding POE (priority) or RES file is returned for a timestamp
# timestamps are always handled in the format YYYYMMDDThhmmss
class OSV(object):
    def __init__(self, outdir_poe, outdir_res):
        self.remote_poe = "https://qc.sentinel1.eo.esa.int/aux_poeorb/"
        self.remote_res = "https://qc.sentinel1.eo.esa.int/aux_resorb/"
        if outdir_poe == outdir_res:
            raise IOError("POE and RES directories must be different")
        self.outdir_poe = outdir_poe
        self.outdir_res = outdir_res
        self.pattern = "S1[AB]_OPER_AUX_(?:POE|RES)ORB_OPOD_[0-9TV_]{48}\.EOF"
        self.pattern_fine = "S1[AB]_OPER_AUX_(?P<type>(?:POE|RES)ORB)_OPOD_(?P<publish>[0-9]{8}T[0-9]{6})_V(?P<start>[0-9]{8}T[0-9]{6})_(?P<stop>[0-9]{8}T[0-9]{6})\.EOF"
        self.sslcontext = ssl._create_unverified_context()

    # create directories if they don't exist yet
    def _init_dir(self):
        for dir in [self.outdir_poe, self.outdir_res]:
            if not os.path.isdir(dir):
                os.makedirs(dir)

    # evaluate the 'type' function argument and return the corresponding local directory
    def _typeEvaluate(self, type):
        if type not in ["POE", "RES"]:
            raise IOError("type must be either 'POE' or 'RES'")
        if type == "POE":
            return self.remote_poe, self.outdir_poe
        else:
            return self.remote_res, self.outdir_res

    # check a server for files
    def catch(self, type="POE", start=None):
        address, outdir = self._typeEvaluate(type)
        address_parse = urlparse(address)
        query = {"page": 1}
        files = []
        if start is not None:
            date = datetime.strptime(start, "%Y%m%dT%H%M%S").strftime("%Y-%m-%d")
            query["validity_start_time"] = "{0}..{1}".format(date, time.strftime("%Y-%m-%d"))
        done = False
        while not done:
            subaddress = urlunparse(address_parse._replace(query=urlencode(query)))
            response = urlopen(subaddress, context=self.sslcontext).read()
            remotes = [os.path.join(address, x) for x in sorted(set(re.findall(self.pattern, response)))]
            if start is not None:
                remotes = [x for x in remotes if self.date(x, "stop") > start]
            counter = 0
            for item in remotes:
                if item not in files:
                    files.append(item)
                    counter += 1
            if counter == 0:
                done = True
            else:
                query["page"] += 1

        if type == "RES":
            files = [x for x in files if self.date(x, "stop") > self.maxdate("POE", "stop")]
        return files

    # extract a date from an OSV file name; types: 'publish', 'start', 'stop'
    def date(self, file, type):
        return re.match(self.pattern_fine, os.path.basename(file)).group(type)

    # delete all RES files for whose date a POE file exists
    def clean_res(self):
        depreceated = [x for x in self.getLocals("RES") if self.date(x, "stop") < self.maxdate("POE", "stop")]
        print "deleting {0} RES files".format(len(depreceated))
        for item in depreceated:
            os.remove(item)

    # get a list of local files
    def getLocals(self, type="POE"):
        address, directory = self._typeEvaluate(type)
        return finder(directory, [self.pattern], regex=True)

    # return the latest date of POE/RES files; datetypes: "publish", "start", "stop"
    def maxdate(self, type="POE", datetype="stop"):
        address, directory = self._typeEvaluate(type)
        files = finder(directory, [self.pattern], regex=True)
        return max([self.date(x, datetype) for x in files]) if len(files) > 0 else None

    # return the corresponding OSV file for the provided time stamp
    def match(self, timestamp):
        for item in dissolve([self.getLocals("POE"), self.getLocals("RES")]):
            if self.date(item, "start") <= timestamp <= self.date(item, "stop"):
                return item
        return None

    # download the newest files
    def retrieve(self, files, type="POE"):
        address, outdir = self._typeEvaluate(type)
        if not os.access(outdir, os.W_OK):
            raise RuntimeError("insufficient directory permissions")
        downloads = [x for x in files if not os.path.isfile(os.path.join(outdir, os.path.basename(x)))]
        print "downloading {0} {1} files".format(len(downloads), type)
        for item in downloads:
            infile = urlopen(item, context=self.sslcontext)
            with open(os.path.join(outdir, os.path.basename(item)), "wb") as outfile:
                outfile.write(infile.read())
            infile.close()

    # perform creating/updating operations for POE and RES files: download newest POE and RES files, delete RES files which can be replaced by newly downloaded POE files
    def update(self):
        self._init_dir()
        files_poe = self.catch("POE", start=self.maxdate("POE", "start"))
        self.retrieve(files_poe, "POE")
        files_res = self.catch("RES", start=self.maxdate("POE", "start"))
        self.retrieve(files_res, "RES")
        self.clean_res()


# unzip S1 scene archives
def unpack(zipfile, outdir):
    with zf.ZipFile(zipfile, "r") as z:
        scene = sorted(z.namelist())[0].strip("/")

        if not os.path.exists(os.path.join(outdir, scene)):
            if not z.testzip():
                print "unzipping data..."
                z.extractall(outdir)
            else:
                raise IOError("corrupt zip")
        else:
            raise IOError("scene already unpacked")
    return
