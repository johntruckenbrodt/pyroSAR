#!/usr/bin/env python2.7
##############################################################
# geocoding of Sentinel-1 GRD data
# John Truckenbrodt 2015
# last update 2016-02-18
##############################################################

"""
The following tasks are performed by executing this script:
-unzipping of a S1 scene directory
-conversion into GAMMA format
-multilooking to retrieve approximately square pixels (20 m resolution)
-mosaicing of all SRTM HGT tiles overlapping with the SAR image (will be automatically downloaded if necessary)
-interpolation of SRTM data gaps
-forward geocoding and pixel area computation for one polarization
-backward geocoding, topographic normalization, conversion to dB and export to GeoTiff for all polarizations
the script supports the following additional options, which are to be added to its console call:
-l create logfiles?
-i keep intermediate files (which are stored in the temporary directory)?
-t transform the supplied DEM to UTM?
-q suppress standard console prints?
further information is given upon calling "python S1_geoGRD.py -h"
the output product will consist of a GeoTiff and the original annotation xml file for each polarization, a quicklook kml and png and,
in case the logfile option was set, a folder containing logfiles for each of the processing steps

change log
2016-01-05: changed the internal quicklook kml document descriptor from "Sentinel-1 Map Overlay" to the image basename (e.g. S1A_IW_20141003T054122_D_38) for disambiguation in Google Earth
2016-02-15: incidence angle is now rounded to next integer instead of just cutting off decimal digits
2016-02-15: changed shebang to explicitely use python2.7 interpreter
2016-02-17: moved the target resolution, func_geoback and func_interp arguments to the argparser
2016-02-17: outsourced the argument parser to script S1_auxil
2016-02-18: integrated orbit state vector correction
2016-02-18: adjusted the multilook factor computation to account for pixel spacing lower than the target resolution and also change factor rounding to always round to the smaller number (so that
            approximation of the target resolution is done but not overdone and the multilooked spacing is still under the target resolution) e.g. spacing 13.2 m and target resolution 20 m: closest
            approximation would be a factor of 2 to get 26.4 m spacing. Yet, the spacing is kept at 13.2 m since the resolution would be unnecessarily downgraded; final resolution resampling is
            performed by the geocoding routine at a later point
"""

import os
import re
import math
import srtm
import shutil
import subprocess as sp
from time import asctime
from gamma.util import correlate, ovs, ISPPar, gamma
from S1_auxil import Identifier, unpack, init_parser, OSV
from ancillary import finder, blockPrint, enablePrint
from pyroSAR import identify


def main(zipfile, tempdir, outdir, srtmdir, transform, logfiles=True, intermediates=False, res_target=20, func_geoback=2, func_interp=0, poedir=None, resdir=None):

    print "##########################################\n{0}\n-------------\nprocessing started: {1}\n".format(zipfile[:-4], asctime())

    id = identify(zipfile)
    outname_base = "_".join([os.path.join(outdir, id.sensor), id.beam, id.start, id.orbit])

    id.unpack(tempdir)


    # # scan the scene folder name for key acquisition attributes and concatenate the base output file name
    # id = Identifier(zipfile)
    # outname_base = "_".join([os.path.join(outdir, id.sat), id.beam, id.start, id.orbit])
    #
    # # unpack the zipped scene archive (if no processed files exist)
    # if len(finder(outdir, [os.path.basename(outname_base)], regex=True)) == 0:
    #     try:
    #         unpack(zipfile, tempdir)
    #         tempdir = os.path.join(tempdir, id.scene)
    #     except IOError as e:
    #         print "{}: {}".format(str(e), id.scene)
    #         return
    # else:
    #     print "scene already processed"
    #     return

    # create logfile folder if this option was selected
    if logfiles:
        path_log = outname_base+"_log"
        if os.path.exists(path_log):
            shutil.rmtree(path_log)
        os.makedirs(path_log)
    else:
        path_log = None

    ######################################################################
    print "converting to GAMMA format..."

    id.convert2gamma(id.scene)
    # try:
    #     run(["reader_old.py", tempdir], outdir=tempdir, logpath=path_log)
    # except ImportWarning:
    #     pass
    # except ImportError:
    #     print "...failed"
    #     return

    # collect all imported files
    # files_mli = finder(tempdir, ["*_mli"])
    files_mli = finder(tempdir, ["*_grd"])

    # correcting orbit state vectors
    if (poedir is not None and resdir is None) or (poedir is None and resdir is not None):
        print "both poedir and resdir must be set in order to utilize orbit state vector correction"
    if poedir is not None and resdir is not None:
        print "correcting orbit state vectors..."
        osv = OSV(poedir, resdir)
        osvfile = osv.match(id.start)
        if osvfile is None:
            print "...no appropriate file found"
        else:
            for item in files_mli:
                gamma(["S1_OPOD_vec", item + ".par", osvfile], logpath=path_log)

    # compute multilooking factors
    par = ISPPar(files_mli[0]+".par")
    rlks = int(res_target//par.range_pixel_spacing)
    azlks = int(res_target//par.azimuth_pixel_spacing)
    rlks = 1 if rlks == 0 else rlks
    azlks = 1 if azlks == 0 else azlks

    # perform multilooking
    if rlks > 1 or azlks > 1:
        print "multilooking..."
        for item in files_mli:
            gamma(["multi_look_MLI", item, item + ".par", item + "2", item + "2.par", rlks, azlks], logpath=path_log)

        # collect all newly created MLIs
        files_mli = [x+"2" for x in files_mli]

    # select master image
    master = files_mli[0]

    base = "_".join(master.split("_")[:-1])+"_"
    dem_seg = base+"dem"
    lut = base+"lut"
    lut_fine = base+"lut_fine"
    sim_sar = base+"sim_sar"
    u = base+"u"
    v = base+"v"
    inc = base+"inc"
    psi = base+"psi"
    pix = base+"pix"
    ls_map = base+"ls_map"
    pixel_area = base+"pixel_area"
    pixel_area2 = base+"pixel_area2"
    offs = base+"offs"
    coffs = base+"coffs"
    coffsets = base+"coffsets"
    snr = base+"snr"
    ellipse_pixel_area = base+"ellipse_pixel_area"
    ratio_sigma0 = base+"ratio_sigma0"

    # read image parameter file for meta information
    par = ISPPar(master+".par")

    # retrieve incidence angle from parameter file
    incidence = str(int(round(par.incidence_angle)))

    # append incidence angle to outname
    outname_base = outname_base+"_"+incidence

    ######################################################################
    # collect srtm files and mosaic them

    # define a name for the output mosaic
    name_srtm = os.path.join(tempdir, "srtm")

    # collect SRTM tiles (if tiles are not found in the defined SRTM directory, they are automatically downloaded to the temporary directory)
    targets = srtm.hgt_collect([x+".par" for x in files_mli], tempdir, demdir=srtmdir)

    if len(targets) == 0:
        print "no SRTM data found"
        return

    print "preparing SRTM data..."

    # mosaic SRTM tiles
    srtm.mosaic(targets, name_srtm)

    # interpolate data gaps
    srtm.fill(name_srtm, name_srtm+"_fill", path_log, replace=True)
    name_srtm += "_fill"

    # project DEM to UTM
    if transform:
        srtm.transform(name_srtm, name_srtm+"_utm", res_target)
        name_srtm += "_utm"
    ######################################################################
    # create DEM products

    # compute DEM oversampling factors
    ovs_lat, ovs_lon = ovs(name_srtm+".par", res_target)

    print "sar image simulation..."
    try:
        gamma(["gc_map", master + ".par", "-", name_srtm + ".par", name_srtm, dem_seg + ".par", dem_seg, lut, ovs_lat, ovs_lon, sim_sar, u, v, inc, psi, pix, ls_map, 8, func_interp], logpath=path_log)
    except sp.CalledProcessError:
        print "...failed"
        return

    ######################################################################
    print "initial pixel area estimation..."
    gamma(["pixel_area", master + ".par", dem_seg + ".par", dem_seg, lut, ls_map, inc, pixel_area], logpath=path_log)

    ######################################################################
    print "cross correlation..."
    try:
        correlate(master, pixel_area, master+"_diff.par", offs, snr, coffs=coffs, coffsets=coffsets, offsets=offs+".txt", path_log=path_log, maxwin=256)
    except RuntimeError:
        print "...failed"
        return
    ######################################################################
    print "supplementing lookup table with offset polynomials..."
    try:
        sim_width = ISPPar(dem_seg+".par").width
        gamma(["gc_map_fine", lut, sim_width, master + "_diff.par", lut_fine, 0], logpath=path_log)
    except sp.CalledProcessError:
        print "...failed"
        return

    ######################################################################
    print "refined pixel area estimation..."
    try:
        gamma(["pixel_area", master + ".par", dem_seg + ".par", dem_seg, lut_fine, ls_map, inc, pixel_area2], logpath=path_log)
    except sp.CalledProcessError:
        print "...failed"
        return

    ######################################################################
    print "radiometric calibration and normalization..."
    try:

        gamma(["radcal_MLI", master, master + ".par", "-", master + "_cal", "-", 0, 0, 1, 0.0, "-", ellipse_pixel_area], logpath=path_log)
        gamma(["ratio", ellipse_pixel_area, pixel_area2, ratio_sigma0, par.range_samples, 1, 1], logpath=path_log)
        for item in files_mli:
            gamma(["product", item, ratio_sigma0, item + "_pixcal", par.range_samples, 1, 1], logpath=path_log)
    except sp.CalledProcessError:
        print "...failed"
        return
    ######################################################################
    print "backward geocoding, incidence angle normalization and conversion to gamma backscatter..."
    for item in files_mli:
        gamma(["geocode_back", item + "_pixcal", par.range_samples, lut_fine, item + "_geo", sim_width, 0, func_geoback], logpath=path_log)

        gamma(["lin_comb", 1, item + "_geo", 0, math.cos(math.radians(par.incidence_angle)), item + "_geo_flat", sim_width], logpath=path_log)
        gamma(["sigma2gamma", item + "_geo_flat", inc, item + "_geo_norm", sim_width], logpath=path_log)

    ######################################################################

    print "creating final tiff files..."
    for item in finder(tempdir, ["*_geo_norm"]):
        polarization = re.findall("[HV]{2}", os.path.basename(item))[0].lower()
        outname = outname_base+"_"+polarization
        try:
            gamma(["data2geotiff", dem_seg + ".par", item, 2, outname + "_geocoded_norm.tif"], logpath=path_log)
        except ImportWarning:
            pass
        except sp.CalledProcessError:
            continue
        annotation_dir = os.path.join(tempdir, "annotation")
        annotation = os.path.join(annotation_dir, [x for x in os.listdir(annotation_dir) if polarization in os.path.basename(x)][0])
        os.rename(annotation, outname+"_annotation.xml")

    ######################################################################
    print "cleaning up..."
    # copy, rename and edit quicklook kml and png
    shutil.copyfile(os.path.join(tempdir, "preview", "map-overlay.kml"), outname_base+"_quicklook.kml")
    shutil.copyfile(os.path.join(tempdir, "preview", "quick-look.png"), outname_base+"_quicklook.png")
    with open(outname_base+"_quicklook.kml", "r") as infile:
        kml = infile.read()
    kml.replace("Sentinel-1 Map Overlay", outname_base)
    kml.replace("quick-look.png", os.path.basename(outname_base+"_quicklook.png"))
    with open(outname_base+"_quicklook.kml", "w") as outfile:
        outfile.write(kml)

    if not intermediates:
        shutil.rmtree(tempdir)

    if logfiles:
        os.rename(path_log, outname_base+"_log")

    print "...done:", asctime()
    print "##########################################"

if __name__ == "__main__":
    parser = init_parser()
    args = parser.parse_args()

    if args.quiet:
        blockPrint()

    main(args.zipfile, args.tempdir, args.outdir, args.srtmdir, args.transform, args.logfiles, args.intermediates, args.targetresolution, args.func_geoback, args.func_interp, args.poedir, args.resdir)

    if args.quiet:
        enablePrint()
