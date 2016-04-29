#!/usr/bin/env python
##############################################################
# universal core routines for processing SAR images is GAMMA
# John Truckenbrodt 2014-2015
##############################################################

"""
This script is intended as a set of generalized processing routines for modularized GAMMA work flows.
The function parametrization is intended to be applicable to any kind of situation and input data set. Thus, instead of choosing a specific parametrization for the data at hand,
core parameters are iterated over a set of values in order to find the one best suited for the task.
The approach of the single routines is likely to still have drawbacks and might fail in certain situations. Testing and suggestions on improvements are very welcome.
"""
import os
import math
from envi import hdr
import subprocess as sp
from auxil import ISPPar
from ancillary import run, Stack, haversine, union


# iterated cross-correlation offset and polynomial estimation
# this function is suited for SLCs (i.e. coregistration) as well as MLIs (i.e. geocoding); the procedure is selected based on the data type of the input data (complex for SLCs but not MLIs)
# starting from a user-defined size, the square offset search window is reduced until a sufficient number of offsets is found or a minimum size is reached
# the number of image kernels for offset search is computed based on the defined percentage of kernel overlap for range and azimuth respectively
# note: the implemented procedure is likely to be overly accurate for many applications; the procedure was chosen based on experience in extremely flat terrain where geocoding is
# problematic due to a lack of image contrast along topographic features and has proven to succeed even in those situations
# INPUT FILES:
# -master: master image
# -slave: image coregistered to the master
# -off: diffpar parameter file
# OUTPUT FILES:
# -offs: offset estimates (fcomplex)
# -offsets: range and azimuth offsets and SNR data (text format)
# -snr: offset estimation SNR (float)
# -coffs: culled range and azimuth offset estimates (fcomplex)
# -coffsets: culled offset estimates and SNR values (text format)
# FURTHER INPUT:
# -path_log: directory for created logfiles (e.g. offset_pwr will create file path_log/offset_pwr.log)
# -maxwin: maximum (initial) window size for offset search (will be iteratively divided by 2 until the minwin is reached)
# -minwin: minimum (final) window size
# -overlap: percentage overlap of the search windows (the number of search windows is computed based on their windows size and the overlap)
# -poly: the polynomial order for range and azimuth offset fitting
# -ovs: image oversampling factor for offset estimation
# -thres: offset estimation SNR quality threshold
# the default value "-" will result in no output file written
def correlate(master, slave, off, offs, snr, offsets="-", coffs="-", coffsets="-", path_log=None, maxwin=2048, minwin=128, overlap=.3, poly=4, ovs=2, thres=7.0):
    path_out = os.path.dirname(off)
    par = ISPPar(master+".par")

    if not os.path.isfile(off):
        run(["create_diff_par", master+".par", "-", off, 1, 0], logpath=path_log)

    if path_log is not None:
        if not os.path.isdir(path_log):
            os.makedirs(path_log)

    if par.image_format in ["FCOMPLEX", "SCOMPLEX"]:
        commands = ["offset_pwr", "offset_fit"]
        mode = "SLC"
    else:
        commands = ["offset_pwrm", "offset_fitm"]
        mode = "MLI"

    # compute the number of estimation windows in azimuth from the defined number of range windows
    dim_ratio = float(par.azimuth_lines)/float(par.range_samples)

    # iteratively reduce the size of the search windows until a sufficient number of offsets was found and/or a final window size is reached
    passed = False
    winsize = maxwin
    while winsize >= minwin:
        try:
            # compute the number of windows needed in range and azimuth based on the number of image pixels in both directions
            nr = int(round((float(par.range_samples)/winsize)*(1+overlap)))
            naz = str(int(int(nr)*dim_ratio))
            if mode == "SLC":
                run([commands[0], master, slave, master+".par", slave+".par", off, offs, snr, winsize, winsize, offsets, ovs, nr, naz, thres], path_out, path_log)
            else:
                run([commands[0], master, slave, off, offs, snr, winsize, winsize, offsets, ovs, nr, naz, thres], path_out, path_log)
            passed = True
            try:
                run([commands[1], offs, snr, off, coffs, coffsets, thres, poly], path_out, path_log)
            except RuntimeError:
                passed = False
        except ValueError:
            continue
        finally:
            winsize /= 2

    if not passed:
        for file in [offs, offsets, snr, coffs, coffsets]:
            if os.path.isfile(file):
                os.remove(file)
        raise RuntimeError("cross-correlation failed; consider verifying scene overlap or choice of polarization")


def geocode(scene, dem, tempdir, outdir, targetres, scaling="linear", func_geoback=2, func_interp=0):
    """
    scaling can be either 'linear', 'db' or a list of both (i.e. ['linear', 'db'])
    """

    scaling = [scaling] if isinstance(scaling, str) else scaling if isinstance(scaling, list) else []
    scaling = union(scaling, ["db", "linear"])
    if len(scaling) == 0:
        raise IOError("wrong input type for parameter scaling")

    scene.unpack(tempdir)

    scene.convert2gamma(scene.scene)

    scene.calibrate()

    images = [x for x in scene.getGammaImages(scene.scene) if x.endswith("_grd")]

    rlks = int(targetres//scene.spacing[0])
    azlks = int(targetres//scene.spacing[1])
    rlks = 1 if rlks == 0 else rlks
    azlks = 1 if azlks == 0 else azlks

    for image in images:
        run(["multi_look_MLI", image, image+".par", image+"_mli", image+"_mli.par", rlks, azlks])

    images = [x+"_mli" for x in images]

    master = images[0]

    dem_seg = master+"_dem"
    lut = master+"_lut"
    lut_fine = master+"_lut_fine"
    sim_sar = master+"_sim_sar"
    u = master+"_u"
    v = master+"_v"
    inc = master+"_inc"
    psi = master+"_psi"
    pix = master+"_pix"
    ls_map = master+"_ls_map"
    offs = master+"_offs"
    coffs = master+"_coffs"
    coffsets = master+"_coffsets"
    snr = master+"_snr"

    ovs_lat, ovs_lon = ovs(dem+".par", targetres)

    path_log = os.path.join(scene.scene, "logfiles")
    if not os.path.isdir(path_log):
        os.makedirs(path_log)

    master_par = ISPPar(master+".par")

    if master_par.image_geometry == "GROUND_RANGE":
        run(["gc_map_grd", master+".par", dem+".par", dem, dem_seg+".par", dem_seg, lut, ovs_lat, ovs_lon, sim_sar, u, v, inc, psi, pix, ls_map, 8, func_interp], logpath=path_log)
    else:
        run(["gc_map", master+".par", "-", dem+".par", dem, dem_seg+".par", dem_seg, lut, ovs_lat, ovs_lon, sim_sar, u, v, inc, psi, pix, ls_map, 8, func_interp], logpath=path_log)

    for item in [dem_seg, sim_sar, u, v, psi, pix, inc]:
        hdr(dem_seg+".par", item+".hdr")

    correlate(master, sim_sar, master+"_diff.par", offs, snr, coffs=coffs, coffsets=coffsets, offsets=offs+".txt", path_log=path_log, maxwin=256)

    try:
        sim_width = ISPPar(dem_seg+".par").width
        run(["gc_map_fine", lut, sim_width, master+"_diff.par", lut_fine, 0], logpath=path_log)
    except sp.CalledProcessError:
        print "...failed"
        return

    for image in images:
        run(["geocode_back", image, master_par.range_samples, lut_fine, image+"_geo", sim_width, 0, func_geoback], logpath=path_log)
        run(["product", image+"_geo", pix, image+"_geo_pix", sim_width, 1, 1, 0], logpath=path_log)
        run(["lin_comb", 1, image+"_geo_pix", 0, math.cos(math.radians(master_par.incidence_angle)), image+"_geo_pix_flat", sim_width], logpath=path_log)
        run(["sigma2gamma", image+"_geo_pix_flat", inc, image+"_geo_norm", sim_width], logpath=path_log)
        hdr(dem_seg+".par", image+"_geo_norm.hdr")

        for scale in scaling:
            suffix = {"linear": "", "db": "_db"}[scale]

            if scale == "db":
                run(["linear_to_dB", image+"_geo_norm", image+"_geo_norm_db", sim_width, 0, -99], logpath=path_log)
                hdr(dem_seg+".par", image+"_geo_norm_db.hdr")

            try:
                infile = image+"_geo_norm{}".format(suffix)
                outfile = os.path.join(outdir, os.path.basename(image)+"_geo_norm{}.tif".format(suffix))
                run(["data2geotiff", dem_seg+".par", infile, 2, outfile], logpath=path_log)
            except ImportWarning:
                pass


# wrapper for iterated initial offset estimation between two SLCs or MLIs
# from a starting search window size of 128x128, initial offsets are search for; if this fails, the window size is increased up to a maximum size of 1024x1024
# once a sufficient number of offsets is found, the window size is again decreased down to 128 and offset search repeated to refine the results on finer levels
# in case too few/no offsets were found on any window size, the routine will throw an error
# additionally, for each window size, the offset search is performed on three different levels of multilooking (factor 3, 2, 1) as this is reported to refine the results;
# if the routine fails on any level of multilooking, the routine will proceed with the next window size
def init_offset(master, slave, off, path_log, thres=7.0):
    path_out = os.path.dirname(off)

    if not os.path.isfile(off):
        run(["create_diff_par", master+".par", "-", off, 1, 0], logpath=path_log)

    par = ISPPar(master + ".par")
    mode = "SLC" if par.image_format in ["FCOMPLEX", "SCOMPLEX"] else "MLI"

    # first SLC offset estimation using orbit data (most important in case of very large offsets)
    if mode == "SLC":
        run(["init_offset_orbit", master+".par", slave+".par", off], path_out, path_log)

    mlk = Spacing(par)
    passed = False
    patchsizes = [128, 256, 512, 1024]
    patchstack = Stack(patchsizes[0])

    while not passed and not patchstack.empty():
        win = patchstack.pop()
        for factor in [3, 2, 1]:
            try:
                if mode == "SLC":
                    run(["init_offset", master, slave, master+".par", slave+".par", off, int(mlk.rlks)*factor, int(mlk.azlks)*factor, "-", "-", "-", "-", thres, win, win], path_out, path_log)
                else:
                    run(["init_offsetm", master, slave, off, int(mlk.rlks)*factor, int(mlk.azlks)*factor, "-", "-", "-", "-", thres, win], path_out, path_log)
                passed = True
            except ValueError:
                passed = False
                # break factor looping as soon as offset finding fails on any factor level
                break

        if not passed:
            # empty stack to finish loop if maximum window size did not result in success
            if win == max(patchsizes):
                patchstack.flush()
            # otherwise, move current window size and its successor onto the stack
            else:
                patchstack.push([win, patchsizes[patchsizes.index(win)+1]])

    if not passed:
        raise RuntimeError("no initial offset found; consider verifying scene overlap or choice of polarization")


# compute DEM oversampling factors for a target resolution in meters
def ovs(parfile, targetres):
    # read DEM parameter file
    dempar = ISPPar(parfile)

    # extract coordinates and pixel posting of the DEM
    if hasattr(dempar, "post_north"):
        post_north, post_east = [abs(float(x)) for x in [dempar.post_north, dempar.post_east]]
    else:
        res_lat, res_lon = [abs(float(x)) for x in [dempar.post_lat, dempar.post_lon]]

        # compute center coordinate
        lat = float(dempar.corner_lat)-(res_lat*(dempar.nlines//2))
        lon = float(dempar.corner_lon)+(res_lon*(dempar.width//2))

        # convert DEM resolution to meters
        post_north = haversine(lat, lon, lat+res_lat, lon)
        post_east = haversine(lat, lon, lat, lon+res_lon)

    # compute resampling factors for the DEM
    ovs_lat = post_north/targetres
    ovs_lon = post_east/targetres
    return ovs_lat, ovs_lon


# compute ground multilooking factors and pixel spacings from an ISPPar object for a defined target resolution
class Spacing(object):
    def __init__(self, par, targetres="automatic"):
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
