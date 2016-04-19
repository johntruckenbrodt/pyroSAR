##############################################################
# Reading and Organizing system for SAR images
# John Truckenbrodt 2016
# last update 2016-03-18
##############################################################
"""
this script is intended to contain several SAR scene identifier classes to read basic metadata from the scene folders/files, convert to GAMMA format and do simple pre-processing
"""

import os
import re
import abc
import math
from osgeo import gdal
from osgeo.gdalconst import GA_ReadOnly
import spatial
import zipfile as zf
import tarfile as tf
from ancillary import finder, parse_literal, run
from time import strptime, strftime
import xml.etree.ElementTree as ElementTree


def identify(scene):
    """Return a metadata handler of the given scene."""
    for handler in ID.__subclasses__():
        try:
            return handler(scene)
        except:
            pass
    raise IOError("data format not supported")


class ID(object):
    """Abstract class for SAR meta data handlers."""

    def bbox(self, outname=None, overwrite=True):
        """Return the bounding box."""
        if outname is None:
            return spatial.bbox(self.getCorners(), self.projection)
        else:
            spatial.bbox(self.getCorners(), self.projection, outname=outname, format="ESRI Shapefile", overwrite=overwrite)

    @abc.abstractmethod
    def calibrate(self, replace=False):
        return

    @property
    def compression(self):
        if os.path.isdir(self.scene):
            return None
        elif zf.is_zipfile(self.scene):
            return "zip"
        elif tf.is_tarfile(self.scene):
            return "tar"
        else:
            return None

    def findfiles(self, pattern):
        if os.path.isdir(self.scene):
            files = [self.scene] if re.search(pattern, os.path.basename(self.scene)) else finder(self.scene, [pattern], regex=True)
        elif zf.is_zipfile(self.scene):
            with zf.ZipFile(self.scene, "r") as zip:
                files = [os.path.join(self.scene, x.strip("/")) for x in zip.namelist() if re.search(pattern, x.strip("/"))]
        elif tf.is_tarfile(self.scene):
            tar = tf.open(self.scene)
            files = [os.path.join(self.scene, x) for x in tar.getnames() if re.search(pattern, x)]
            tar.close()
        else:
            files = [self.scene] if re.search(pattern, self.scene) else []
        return files

    @abc.abstractmethod
    def convert2gamma(self, directory):
        return

    def gdalinfo(self, scene):
        self.scene = os.path.realpath(scene)
        files = self.findfiles("(?:\.[NE][12]$|DAT_01\.001$|product\.xml$)")

        if len(files) == 1:
            prefix = {"zip": "/vsizip/", "tar": "/vsitar/", None: ""}[self.compression]
            self.file = files[0]
        elif len(files) > 1:
            raise IOError("file ambiguity detected")
        else:
            raise IOError("file type not supported")

        ext_lookup = {".N1": "ASAR", ".E1": "ERS1", ".E2": "ERS2"}
        extension = os.path.splitext(self.file)[1]
        if extension in ext_lookup:
            self.sensor = ext_lookup[extension]

        img = gdal.Open(prefix+self.file, GA_ReadOnly)
        meta = img.GetMetadata()
        self.cols, self.rows, self.bands = img.RasterXSize, img.RasterYSize, img.RasterCount
        self.projection = img.GetGCPProjection()
        self.gcps = [[[x.GCPPixel, x.GCPLine], [x.GCPX, x.GCPY, x.GCPZ]] for x in img.GetGCPs()]
        img = None

        for item in meta:
            entry = [item, parse_literal(meta[item].strip())]

            # todo: check module time for more general approaches
            for timeformat in ["%d-%b-%Y %H:%M:%S.%f", "%Y%m%d%H%M%S%f", "%Y-%m-%dT%H:%M:%S.%fZ"]:
                try:
                    entry[1] = strftime("%Y%m%dT%H%M%S", strptime(entry[1], timeformat))
                except (TypeError, ValueError):
                    pass

            if re.search("(?:LAT|LONG)", entry[0]):
                entry[1] /= 1000000.
            setattr(self, entry[0], entry[1])

    @abc.abstractmethod
    def getCorners(self):
        return

    def getGammaImages(self, directory):
        return [x for x in finder(directory, [self.outname_base], regex=True) if not x.endswith(".par")]

    def summary(self):
        for item in sorted(self.__dict__.keys()):
            print "{0}: {1}".format(item, getattr(self, item))

    @abc.abstractmethod
    def unpack(self, directory):
        return

    def _unpack(self, directory):
        if not os.path.isdir(directory):
            os.makedirs(directory)
        if tf.is_tarfile(self.scene):
            archive = tf.open(self.scene, "r")
            names = archive.getnames()
            header = os.path.commonprefix(names)
            if len(header) > 0:
                if archive.getmember(header).isdir():
                    for item in sorted(names):
                        if item != header:
                            member = archive.getmember(item)
                            outname = os.path.join(directory, item.replace(header+"/", ""))
                            if member.isdir():
                                os.makedirs(outname)
                            else:
                                with open(outname, "w") as outfile:
                                    outfile.write(member.tobuf())
                    archive.close()
                    return
            archive.extractall(directory)
            archive.close()
        elif zf.is_zipfile(self.scene):
            archive = zf.ZipFile(self.scene, "r")
            names = archive.namelist()
            header = os.path.commonprefix(names)
            if len(header) > 0:
                if header.endswith("/"):
                    for item in sorted(names):
                        if item != header:
                            outname = os.path.join(directory, item.replace(header, ""))
                            if item.endswith("/"):
                                os.makedirs(outname)
                            else:
                                with open(outname, "w") as outfile:
                                    outfile.write(archive.read(item))
                    archive.close()
                    return
            archive.extractall(directory)
            archive.close()
        self.scene = directory
        self.file = os.path.join(self.scene, os.path.basename(self.file))


class CEOS(ID):
    # todo: What sensors other than ERS1, ERS2 should be included?
    # todo: add a pattern to check if the scene could be handled by CEOS
    def __init__(self, scene):
        self.gdalinfo(scene)
        self.sensor = self.CEOS_MISSION_ID
        self.start = self.CEOS_ACQUISITION_TIME
        self.incidence = self.CEOS_INC_ANGLE
        self.spacing = [self.CEOS_PIXEL_SPACING_METERS, self.CEOS_LINE_SPACING_METERS]

        # todo: check whether this is correct:
        self.orbit = "D" if self.CEOS_PLATFORM_HEADING > 180 else "A"
        self.k_db = -10*math.log(self.CEOS_CALIBRATION_CONSTANT_K, 10)
        self.sc_db = {"ERS1": 59.61, "ERS2": 60}[self.sensor]
        self.outname_base = "{0}______{1}".format(*[self.sensor, self.start])

    #todo: should define a calibrate function
    def getCorners(self):
        lat = [x[1][1] for x in self.gcps]
        lon = [x[1][0] for x in self.gcps]
        return {"xmin": min(lon), "xmax": max(lon), "ymin": min(lat), "ymax": max(lat)}

    def convert2gamma(self, directory):
        if self.sensor in ["ERS1", "ERS2"]:
            outname = os.path.join(directory, self.outname_base+"_VV_slc")
            lea = os.path.join(self.scene, "LEA_01.001")
            title = os.path.basename(self.findfiles("\.PS$")[0]).replace(".PS", "")
            run(["par_ESA_ERS", lea, outname+".par", self.file, outname], inlist=[title])
        else:
            raise NotImplementedError("sensor {} not implemented yet".format(self.sensor))

    def unpack(self, directory):
        if self.sensor in ["ERS1", "ERS2"]:
            outdir = os.path.join(directory, re.sub("\.[EN][12]\.PS$", "", os.path.basename(self.findfiles("\.PS$")[0])))
            self._unpack(outdir)
        else:
            raise NotImplementedError("sensor {} not implemented yet".format(self.sensor))

# id = identify("/geonfs01_vol1/ve39vem/ERS/ERS1_0132_2529_20dec95")
# id = identify("/geonfs01_vol1/ve39vem/ERS/ERS1_0132_2529_20dec95.zip")


class ESA(ID):
    """Handle SAR data of the ESA format."""
    def __init__(self, scene):

        self.pattern = r"(?P<product_id>(?:SAR|ASA)_(?:IM(?:S|P|G|M|_)|AP(?:S|P|G|M|_)|WV(?:I|S|W|_))_[012B][CP])" \
                       r"(?P<processing_stage_flag>[A-Z])" \
                       r"(?P<originator_ID>[A-Z\-]{3})" \
                       r"(?P<start_day>[0-9]{8})_" \
                       r"(?P<start_time>[0-9]{6})_" \
                       r"(?P<duration>[0-9]{8})" \
                       r"(?P<phase>[0-9A-Z]{1})" \
                       r"(?P<cycle>[0-9]{3})_" \
                       r"(?P<relative_orbit>[0-9]{5})_" \
                       r"(?P<absolute_orbit>[0-9]{5})_" \
                       r"(?P<counter>[0-9]{4})\." \
                       r"(?P<satellite_ID>[EN][12])" \
                       r"(?P<extension>(?:\.zip|\.tar\.gz|))$"

        self.scene = os.path.realpath(scene)
        self.gdalinfo(self.scene)
        self.orbit = self.SPH_PASS[0]
        self.start = self.MPH_SENSING_START
        self.stop = self.MPH_SENSING_STOP
        self.spacing = [self.SPH_RANGE_SPACING, self.SPH_AZIMUTH_SPACING]
        self.looks = [self.SPH_RANGE_LOOKS, self.SPH_AZIMUTH_LOOKS]
        self.outname_base = "{0}______{1}".format(*[self.sensor, self.start])

    def getCorners(self):
        lon = [getattr(self, x) for x in self.__dict__.keys() if re.search("LONG", x)]
        lat = [getattr(self, x) for x in self.__dict__.keys() if re.search("LAT", x)]
        return {"xmin": min(lon), "xmax": max(lon), "ymin": min(lat), "ymax": max(lat)}

    def convert2gamma(self, directory):
        self.gammadir = directory
        outname = os.path.join(directory, self.outname_base)
        if len(self.getGammaImages(directory)) == 0:
            run(["par_ASAR", self.file, outname])
            os.remove(outname+".hdr")
            for item in finder(directory, [os.path.basename(outname)], regex=True):
                ext = ".par" if item.endswith(".par") else ""
                base = os.path.basename(item).strip(ext)
                base = base.replace(".", "_")
                base = base.replace("PRI", "mli")
                base = base.replace("SLC", "slc")
                newname = os.path.join(directory, base+ext)
                os.rename(item, newname)
        else:
            raise IOError("scene already processed")

    def calibrate(self, replace=False):
        k_db = {"ASAR": 55., "ERS1": 58.24, "ERS2": 59.75}[self.sensor]
        inc_ref = 90. if self. sensor == "ASAR" else 23.
        candidates = [x for x in self.getGammaImages(self.gammadir) if not x.endswith("_cal") and not os.path.isfile(x+"_cal")]
        for image in candidates:
            run(["radcal_PRI", image, image+".par", image+"_cal", image+"_cal.par", k_db, inc_ref])
            if replace:
                os.remove(image)
                os.remove(image+".par")

    def unpack(self, directory):
        outdir = os.path.join(directory, os.path.splitext(os.path.basename(self.file))[0])
        self._unpack(outdir)
# id = identify("/geonfs01_vol1/ve39vem/swos/ASA_APP_1PTDPA20040102_102928_000000162023_00051_09624_0240.N1")
# id = identify("/geonfs01_vol1/ve39vem/swos/SAR_IMP_1PXASI19920419_110159_00000017C083_00323_03975_8482.E1")
# id = identify("/geonfs01_vol1/ve39vem/swos/ER01_SAR_IMP_1P_19920419T110159_19920419T110216_IPA_03975_0000.ESA.tar.gz")

# scenes = finder("/geonfs01_vol1/ve39vem/swos", ["*"])
# counter = 0
# for scene in scenes:
#     counter += 1
#     try:
#         x = ESA(scene)
#         if len(x.findfiles(x.pattern)) == 0:
#             print scene
#             print x.findfiles("[EN][12]$")[0]
#             print "---------------------------------"
#     except RuntimeError as rue:
#         print scene
#         print rue
#         print"---------------------------------"
#     # progress = float(counter)/len(scenes)*100
#     # print progress
#     # if progress % 10 == 0:
#     #     print progress

# scenes = finder("/geonfs01_vol1/ve39vem/swos", ["*.tar.gz"])
# for scene in scenes:
#     print scene


class RS2(ID):
    def __init__(self, scene):

        self.pattern = r'^(?:RS2|RSAT2)_(?:OK[0-9]+)_(?:PK[0-9]+)_(?:DK[0-9]+)_' \
                       r'(?P<beam>[0-9A-Z]+)_' \
                       r'(?P<date>[0-9]{8})_' \
                       r'(?P<time>[0-9]{6})_' \
                       r'(?P<pols>[HV]{2}_' \
                       r'(?P<level>SLC|SGX|SGF|SCN|SCW|SSG|SPG)$'

        self.sensor = "RS2"
        self.scene = os.path.realpath(scene)
        self.gdalinfo(self.scene)
        self.start = self.ACQUISITION_START_TIME
        self.incidence = (self.FAR_RANGE_INCIDENCE_ANGLE + self.NEAR_RANGE_INCIDENCE_ANGLE)/2
        self.spacing = [self.PIXEL_SPACING, self.LINE_SPACING]
        self.orbit = self.ORBIT_DIRECTION[0]

    def getCorners(self):
        lat = [x[1][1] for x in self.gcps]
        lon = [x[1][0] for x in self.gcps]
        return {"xmin": min(lon), "xmax": max(lon), "ymin": min(lat), "ymax": max(lat)}

# id = identify("/geonfs01_vol1/ve39vem/RS2/RS2_OK53107_PK504800_DK448361_FQ1_20140606_055403_HH_VV_HV_VH_SLC.zip")


# todo: check self.file and self.scene assignment after unpacking
class SAFE(ID):
    def __init__(self, scene):

        self.scene = os.path.realpath(scene)

        self.pattern = r"^(?P<sensor>S1[AB])_" \
                       r"(?P<beam>S1|S2|S3|S4|S5|S6|IW|EW|WV|EN|N1|N2|N3|N4|N5|N6|IM)_" \
                       r"(?P<product>SLC|GRD|OCN)(?:F|H|M|_)_" \
                       r"(?:1|2)" \
                       r"(?P<category>S|A)" \
                       r"(?P<pols>SH|SV|DH|DV|HH|HV|VV|VH)_" \
                       r"(?P<start>[0-9]{8}T[0-9]{6})_" \
                       r"(?P<stop>[0-9]{8}T[0-9]{6})_" \
                       r"(?:[0-9]{6})_" \
                       r"(?:[0-9A-F]{6})_" \
                       r"(?:[0-9A-F]{4})" \
                       r"\.SAFE$"

        self.pattern_ds = r"^s1[ab]-" \
                          r"(?P<swath>s[1-6]|iw[1-3]?|ew[1-5]?|wv[1-2]|n[1-6])-" \
                          r"(?P<product>slc|grd|ocn)-" \
                          r"(?P<pol>hh|hv|vv|vh)-" \
                          r"(?P<start>[0-9]{8}t[0-9]{6})-" \
                          r"(?P<stop>[0-9]{8}t[0-9]{6})-" \
                          r"(?:[0-9]{6})-(?:[0-9a-f]{6})-" \
                          r"(?P<id>[0-9]{3})" \
                          r"\.xml$"

        files = self.findfiles(self.pattern)
        if len(files) == 1:
            self.file = files[0]
        elif len(files) == 0:
            raise IOError("folder does not match S1 scene naming convention")
        else:
            raise IOError("file ambiguity detected")

        match = re.match(re.compile(self.pattern), os.path.basename(self.file))

        if not match:
            raise IOError("folder does not match S1 scene naming convention")
        for key in re.compile(self.pattern).groupindex:
            setattr(self, key, match.group(key))

        self.orbit = "D" if float(re.findall("[0-9]{6}", self.start)[1]) < 120000 else "A"
        self.projection = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

    def calibrate(self, replace=False):
        print "calibration already performed during import"

    def convert2gamma(self, directory):
        if self.compression is not None:
            raise RuntimeError("scene is not yet unpacked")
        if self.product == "OCN":
            raise IOError("Sentinel-1 OCN products are not supported")
        if self.category == "A":
            raise IOError("Sentinel-1 annotation-only products are not supported")

        for xml_ann in finder(os.path.join(self.scene, "annotation"), [self.pattern_ds]):
            base = os.path.basename(xml_ann)
            match = re.compile(self.pattern_ds).match(os.path.basename(base))

            tiff = os.path.join(self.scene, "measurement", base.replace(".xml", ".tiff"))
            xml_cal = os.path.join(self.scene, "annotation", "calibration", "calibration-" + base)
            # the use of the noise xml file has been found to occasionally cause severe image artifacts of manifold nature and is thus excluded
            # the reason (GAMMA command error vs. bad ESA xml file entry) is yet to be discovered
            # xml_noise = os.path.join(self.scene, "annotation", "calibration", "noise-" + base)
            xml_noise = "-"
            fields = (self.sensor, match.group("swath"), self.start, match.group("pol").upper())
            if match.group("prod") == "SLC":
                name = os.path.join(directory, "{0}_{1}__{2}_{3}_slc".format(*fields))
                # print ["par_S1_SLC", tiff, xml_ann, xml_cal, xml_noise, name + ".par", name, name + ".tops_par"]
                run(["par_S1_SLC", tiff, xml_ann, xml_cal, xml_noise, name + ".par", name, name + ".tops_par"])
            else:
                name = os.path.join(directory, "{0}______{2}_{3}_mli".format(*fields))
                # print ["par_S1_GRD", tiff, xml_ann, xml_cal, xml_noise, name + ".par", name]
                run(["par_S1_GRD", tiff, xml_ann, xml_cal, xml_noise, name + ".par", name])

    def getCorners(self):
        if self.compression == "zip":
            with zf.ZipFile(self.scene, "r") as z:
                kml = z.open([x for x in z.namelist() if re.search("map-overlay\.kml", x)][0], "r").read()
        elif self.compression == "tar":
            tar = tf.open(self.scene, "r")
            kml = tar.extractfile().read()
            tar.close()
        else:
            with open(finder(self.scene, ["*map-overlay.kml"])[0], "r") as infile:
                kml = infile.read()
        elements = ElementTree.fromstring(kml).findall(".//coordinates")

        coordinates = [x.split(",") for x in elements[0].text.split()]
        lat = [float(x[1]) for x in coordinates]
        lon = [float(x[0]) for x in coordinates]
        return {"xmin": min(lon), "xmax": max(lon), "ymin": min(lat), "ymax": max(lat)}

    def unpack(self, directory):
        outdir = os.path.join(directory, os.path.basename(self.file))
        self._unpack(outdir)

# id = identify("/geonfs01_vol1/ve39vem/S1/archive/S1A_EW_GRDM_1SDH_20150408T053103_20150408T053203_005388_006D8D_5FAC.zip")


class ERS(object):
    # todo: add a pattern to check if the scene could be handled
    def __init__(self, scene):
        try:
            lea = finder(scene, ["LEA_01.001"])[0]
        except IndexError:
            raise IOError("wrong input format; no leader file found")
        with open(lea, "r") as infile:
            text = infile.read()
        # extract frame id
        frame_index = re.search("FRAME=", text).end()
        self.frame = text[frame_index:frame_index+4]
        # extract calibration meta information
        stripper = " \t\r\n\0"
        self.sensor = text[(720+395):(720+411)].strip(stripper)
        self.date = int(text[(720+67):(720+99)].strip(stripper)[:8])
        self.proc_fac = text[(720+1045):(720+1061)].strip(stripper)
        self.proc_sys = text[(720+1061):(720+1069)].strip(stripper)
        self.proc_vrs = text[(720+1069):(720+1077)].strip(stripper)
        text_subset = text[re.search("FACILITY RELATED DATA RECORD \[ESA GENERAL TYPE\]", text).start()-13:]
        self.cal = -10*math.log(float(text_subset[663:679].strip(stripper)), 10)
        self.antenna_flag = text_subset[659:663].strip(stripper)

        # the following section is only relevant for PRI products and can be considered future work
        # select antenna gain correction lookup file from extracted meta information
        # the lookup files are stored in a subfolder CAL which is included in the pythonland software package
        # if sensor == "ERS1":
        #     if date < 19950717:
        #         antenna = "antenna_ERS1_x_x_19950716"
        #     else:
        #         if proc_sys == "VMP":
        #             antenna = "antenna_ERS2_VMP_v68_x" if proc_vrs >= 6.8 else "antenna_ERS2_VMP_x_v67"
        #         elif proc_fac == "UKPAF" and date < 19970121:
        #             antenna = "antenna_ERS1_UKPAF_19950717_19970120"
        #         else:
        #             antenna = "antenna_ERS1"
        # else:
        #     if proc_sys == "VMP":
        #         antenna = "antenna_ERS2_VMP_v68_x" if proc_vrs >= 6.8 else "antenna_ERS2_VMP_x_v67"
        #     elif proc_fac == "UKPAF" and date < 19970121:
        #         antenna = "antenna_ERS2_UKPAF_x_19970120"
        #     else:
        #         antenna = "antenna_ERS2"
