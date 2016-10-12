#!/usr/bin/env python
##############################################################
# universal core routines for processing SAR images is GAMMA
# John Truckenbrodt 2014-2016
##############################################################

"""
This script is intended as a set of generalized processing routines for modularized GAMMA work flows.
The function parametrization is intended to be applicable to any kind of situation and input data set. Thus, instead of choosing a specific parametrization for the data at hand,
core parameters are iterated over a set of values in order to find the one best suited for the task.
The approach of the single routines is likely to still have drawbacks and might fail in certain situations. Testing and suggestions on improvements are very welcome.
"""
import os
import re
import math
import envi
import subprocess as sp
from . import ISPPar, Spacing
from .error import gammaErrorHandler
from ancillary import run, Stack, union
from spatial import haversine
import pyroSAR


def process(cmd, outdir=None, logpath=None, inlist=None, void=True):
    log = os.path.join(logpath, cmd[0]+'.log') if logpath else None
    out, err = run(cmd, outdir=outdir, logfile=log, inlist=inlist, void=void)
    gammaErrorHandler(err)
    if not void:
        return out, err


def slc_corners(parfile):
    out, err = process(['SLC_corners', parfile], void=False)
    pts = {}
    for line in out.split('\n'):
        if line.startswith('min. latitude'):
            pts['ymin'], pts['ymax'] = [float(x) for x in re.findall('[0-9]+\.[0-9]+', line)]
        elif line.startswith('min. longitude'):
            pts['xmin'], pts['xmax'] = [float(x) for x in re.findall('[0-9]+\.[0-9]+', line)]
    return pts


# INPUT FILES:
# -master: master image
# -slave: image coregistered to the master
# -off: diffpar parameter file
# OUTPUT FILES:
# -offs: offset estimates (fcomplex)
# -offsets: range and azimuth offsets and cross-correlation data (text format)
# -ccp: cross-correlation of each patch (0.0->1.0) (float)
# -coffs: culled range and azimuth offset estimates (fcomplex)
# -coffsets: culled offset estimates and cross-correlation values (text format)
# FURTHER INPUT:
# -path_log: directory for created logfiles (e.g. offset_pwr will create file path_log/offset_pwr.log)
# -maxwin: maximum (initial) window size for offset search (will be iteratively divided by 2 until the minwin is reached)
# -minwin: minimum (final) window size
# -overlap: percentage overlap of the search windows (the number of search windows is computed based on their windows size and the overlap)
# -poly: the polynomial order for range and azimuth offset fitting
# -ovs: image oversampling factor for offset estimation
# -thres: cross-correlation threshold
# the default value "-" will result in no output file written
def correlate(master, slave, off, offs, ccp, offsets="-", coffs="-", coffsets="-", path_log=None, maxwin=2048, minwin=16, overlap=.3, poly=4, ovs=2, thres=0.15):
    """
    iterated cross-correlation offset and polynomial estimation
    this function is suited for SLCs (i.e. coregistration) as well as MLIs (i.e. geocoding); the procedure is selected based on the data type of the input data (complex for SLCs but not MLIs)
    starting from a user-defined size, the square offset search window is reduced until a sufficient number of offsets is found or a minimum size is reached
    the number of image kernels for offset search is computed based on the defined percentage of kernel overlap for range and azimuth respectively
    note: the implemented procedure is likely to be overly accurate for many applications; the procedure was chosen based on experience in extremely flat terrain where geocoding is
    problematic due to a lack of image contrast along topographic features and has proven to succeed even in those situations
    Args:
        (input) master: master image
        (input) slave: image coregistered to the master
        off: diffpar parameter file
        (output) offs: offset estimates (fcomplex)
        (output) ccp: cross-correlation of each patch (0.0->1.0) (float)
        (output) offsets: range and azimuth offsets and cross-correlation data (text format)
        (output) coffs: culled range and azimuth offset estimates (fcomplex)
        (output) coffsets: culled offset estimates and cross-correlation values (text format)
        path_log: directory for created logfiles (e.g. offset_pwr will create file path_log/offset_pwr.log)
        maxwin: maximum (initial) window size for offset search (will be iteratively divided by 2 until the minwin is reached)
        minwin: minimum (final) window size
        overlap: percentage overlap of the search windows (the number of search windows is computed based on their windows size and the overlap)
        poly: the polynomial order for range and azimuth offset fitting
        ovs: image oversampling factor for offset estimation
        thres: cross-correlation threshold
    """

    path_out = os.path.dirname(off)
    par = ISPPar(master+".par")

    if not os.path.isfile(off):
        process(["create_diff_par", master + ".par", "-", off, 1, 0], logpath=path_log)

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
        # compute the number of windows needed in range and azimuth based on the number of image pixels in both directions
        nr = int(round((float(par.range_samples)/winsize)*(1+overlap)))
        naz = str(int(int(nr)*dim_ratio))

        args = [off, offs, ccp, winsize, winsize, offsets, ovs, nr, naz, thres]
        try:
            if mode == "SLC":
                process([commands[0], master, slave, master + ".par", slave + ".par"] + args, path_out, path_log)
            else:
                process([commands[0], master, slave] + args, path_out, path_log)
            passed = True
            try:
                process([commands[1], offs, ccp, off, coffs, coffsets, thres, poly], path_out, path_log)
            except RuntimeError:
                passed = False
        except ValueError:
            continue
        finally:
            winsize /= 2

    if not passed:
        for file in [offs, offsets, ccp, coffs, coffsets]:
            if os.path.isfile(file):
                os.remove(file)
        raise RuntimeError("cross-correlation failed; consider verifying scene overlap or choice of polarization")


def cc_find(master, slave, diffpar, offs, ccp, winsize, mode, offsets="-", n_ovr=2, nr=400, naz=400, thres=.15, c_ovr=4):
    """
    *** Offset estimation between MLI or SLC images using intensity cross-correlation ***
    Args:
        master    (input) real valued intensity image 1 (reference)
        slave     (input) real valued intensity image 2
        diffpar   DIFF/GEO parameter file
        offs      (output) offset estimates (fcomplex)
        ccp       (output) cross-correlation of each patch (0.0->1.0) (float)
        winsize   patch size in range and azimuth (pixels)
        offsets   (output) range and azimuth offsets and cross-correlation data in text format, enter - for no output
        n_ovr     MLI oversampling factor (integer 2**N (1,2,4,8), enter - for default: 2)
        nr        number of offset estimates in range direction (enter - for default from offset parameter file)
        naz       number of offset estimates in azimuth direction (enter - for default from offset parameter file)
        thres     cross-correlation threshold (enter - for default from offset parameter file)
        c_ovr     correlation function oversampling factor (integer 2**N (1,2,4,8) default: 4)
    """
    args = [diffpar, offs, ccp, str(winsize), str(winsize), offsets, str(n_ovr), str(nr), str(naz), str(thres), str(c_ovr)]

    if mode == "MLI":
        cmd = ["offset_pwrm", master, slave]+args
    elif mode == "SLC":
        cmd = ["offset_pwr", master, slave, master+".par", slave+".par"]+args
    else:
        raise ValueError("mode must be either 'MLI' or 'SLC'")

    proc = sp.Popen(cmd, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = proc.communicate()
    for line in out.split("\n"):
        if line.startswith("number of offsets above correlation threshold"):
            vals = map(float, re.findall("[0-9]+", line))
            return vals[0]/vals[1]


def cc_fit(offs, ccp, diffpar, mode, coffs="-", coffsets="-", thres=.15, npoly=4):
    """
    *** Range and azimuth offset polynomial estimation ***
    Args:
          offs          (input) range and azimuth offset estimates for each patch (fcomplex)
          ccp           (input) cross-correlation or SNR of each patch (float)
          diffpar       (input) ISP offset/interferogram parameter file
          coffs         (output) culled range and azimuth offset estimates (fcomplex, enter - for none)
          coffsets      (output) culled offset estimates and cross-correlation values (text format, enter - for none)
          thres         cross-correlation threshold (enter - for default from diffpar)
          npoly         number of model polynomial parameters (enter - for default, 1, 3, 4, 6, default: 4)
    """
    args = [offs, ccp, diffpar, coffs, coffsets, str(thres), str(npoly)]

    if mode == "MLI":
        cmd = ["offset_fitm"]+args
    elif mode == "SLC":
        cmd = ["offset_fit"]+args
    else:
        raise ValueError("mode must be either 'MLI' or 'SLC'")

    proc = sp.Popen(cmd, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = proc.communicate()
    fit = None
    for line in out.split("\n"):
        if line.startswith("final model fit std. dev."):
            fit = tuple(map(float, re.findall("[0-9]+.[0-9]+", line)))
    print "final model fit std. dev. (range, azimuth):", fit
    if not fit:
        raise RuntimeError("failed at fitting offset estimation; command:\n{}".format(" ".join(cmd)))


def cc_fit2(offs, ccp, diffpar, mode, coffs="-", coffsets="-", thres=.15, npoly=4):
    """
    *** Range and azimuth offset polynomial estimation ***
    Args:
          offs          (input) range and azimuth offset estimates for each patch (fcomplex)
          ccp           (input) cross-correlation or SNR of each patch (float)
          diffpar       (input) ISP offset/interferogram parameter file
          coffs         (output) culled range and azimuth offset estimates (fcomplex, enter - for none)
          coffsets      (output) culled offset estimates and cross-correlation values (text format, enter - for none)
          thres         cross-correlation threshold (enter - for default from diffpar)
          npoly         number of model polynomial parameters (enter - for default, 1, 3, 4, 6, default: 4)
    """
    args = [offs, ccp, diffpar, coffs, coffsets, str(thres), str(npoly)]

    if mode == "MLI":
        cmd = ["offset_fitm"]+args
    elif mode == "SLC":
        cmd = ["offset_fit"]+args
    else:
        raise ValueError("mode must be either 'MLI' or 'SLC'")

    proc = sp.Popen(cmd, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = proc.communicate()
    errormessage = 'ERROR: insufficient offset points left after culling to determine offset model parameters'
    if errormessage in err:
        raise RuntimeError(errormessage)
    for line in out.split("\n"):
        if line.startswith("final model fit std. dev."):
            return tuple(map(float, re.findall("[0-9]+.[0-9]+", line)))
    raise RuntimeError(err)


def cc(master, slave, diffpar, offs, ccp, mode, offsets="-", nr=400, naz=400, thres=.15, c_ovr=4, minwin=8, maxwin=1024):
    best = 0
    winsize = minwin
    bestwin = 8
    while winsize <= maxwin:
        ratio = cc_find(master, slave, diffpar, offs, ccp, winsize, mode, offsets, 1, nr, naz, thres, c_ovr)
        if ratio > best:
            best = ratio
            bestwin = winsize
        else:
            break
        winsize *= 2
    print "optimal window size:", bestwin
    ratio = cc_find(master, slave, diffpar, offs, ccp, bestwin, mode, n_ovr=8)
    print "offset acceptance ratio:", ratio
    return bestwin


def cc_find2(master, slave, diffpar, offs, ccp, rwin, azwin, mode, offsets="-", n_ovr=2, nr=400, naz=400, thres=.15, c_ovr=4):
    """
    *** Offset estimation between MLI or SLC images using intensity cross-correlation ***
    Args:
        master    (input) real valued intensity image 1 (reference)
        slave     (input) real valued intensity image 2
        diffpar   DIFF/GEO parameter file
        offs      (output) offset estimates (fcomplex)
        ccp       (output) cross-correlation of each patch (0.0->1.0) (float)
        winsize   patch size in range and azimuth (pixels)
        offsets   (output) range and azimuth offsets and cross-correlation data in text format, enter - for no output
        n_ovr     MLI oversampling factor (integer 2**N (1,2,4,8), enter - for default: 2)
        nr        number of offset estimates in range direction (enter - for default from offset parameter file)
        naz       number of offset estimates in azimuth direction (enter - for default from offset parameter file)
        thres     cross-correlation threshold (enter - for default from offset parameter file)
        c_ovr     correlation function oversampling factor (integer 2**N (1,2,4,8) default: 4)
    """
    args = [diffpar, offs, ccp, str(rwin), str(azwin), offsets, str(n_ovr), str(nr), str(naz), str(thres), str(c_ovr)]
    if mode == "MLI":
        cmd = ["offset_pwrm", master, slave]+args
    elif mode == "SLC":
        cmd = ["offset_pwr", master, slave, master+".par", slave+".par"]+args
    else:
        raise ValueError("mode must be either 'MLI' or 'SLC'")
    proc = sp.Popen(cmd, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = proc.communicate()
    errormessage = 'ERROR: no offsets found above correlation threshold'
    if errormessage in err:
        raise RuntimeError(errormessage)
    for line in out.split("\n"):
        if line.startswith("number of offsets above correlation threshold"):
            vals = map(float, re.findall("[0-9]+", line))
            return int(vals[0])
    raise RuntimeError(err)


def correlate2(master, slave, diffpar, offs, ccp, offsets="-", coffs="-", coffsets="-", npoly=4, thres=0.15, minwin=8, maxwin=2048, niter=2, path_log=None):
    if not os.path.isfile(diffpar):
        process(["create_diff_par", master + ".par", "-", diffpar, 1, 0], logpath=path_log)
    if path_log is not None:
        if not os.path.isdir(path_log):
            os.makedirs(path_log)
    mode = "SLC" if ISPPar(master+".par").image_format in ["FCOMPLEX", "SCOMPLEX"] else "MLI"
    print "#####################################\ncross-correlation offset estimation"
    for i in range(0, niter):
        print "#############\niteration {}".format(i+1)
        cc(master, slave, diffpar, offs, ccp, mode, offsets, thres=thres, minwin=minwin, maxwin=maxwin)
        cc_fit(offs, ccp, diffpar, mode, coffs, coffsets, thres, npoly)
    print "#############\n...done\n#####################################"


def find_poly(offs, ccp, diffpar, mode, coffs, coffsets, thres=.15):
    bestpoly = 1
    bestfit = 999
    for npoly in [1, 3, 4]:
        try:
            fit = cc_fit2(offs, ccp, diffpar, mode, coffs, coffsets, thres, npoly)
            if sum(fit) < bestfit:
                bestpoly = npoly
                bestfit = fit
        except RuntimeError:
            continue
    print "best polynomial order:", bestpoly
    return cc_fit2(offs, ccp, diffpar, mode, coffs, coffsets, thres, bestpoly)


# def correlate3(master, slave, diffpar, offs, ccp, offsets="-", coffs="-", coffsets="-"):
#     if not os.path.isfile(diffpar):
#         process(["create_diff_par", master + ".par", "-", diffpar, 1, 0], logpath=path_log)
#     par = ISPPar(master+".par")
#     mode = "SLC" if par.image_format in ["FCOMPLEX", "SCOMPLEX"] else "MLI"
#     rwin = 2**(math.log(par.range_pixel_spacing*10, 2)//1)
#     azwin = 2**(math.log(par.azimuth_pixel_spacing*10, 2)//1)
#     nr = 3*(par.range_samples//rwin)
#     naz = 3*(par.azimuth_lines//azwin)
#     n_ovr = 2
#     thres = .15
#     # cc_find2(master, slave, diffpar, offs, ccp, rwin, azwin, mode, offsets, n_ovr, nr, naz, thres)
#     # print find_poly(offs, ccp, diffpar, mode, coffs, coffsets, thres)
#     # # cc_find2(master, slave, diffpar, offs, ccp, rwin/2, azwin/2, mode, offsets, n_ovr, nr*2, naz*2, thres)
#     # # fit2 = find_poly(offs, ccp, diffpar, mode, coffs, coffsets, thres)
#     #
#     cc_find2(master, slave, diffpar, offs, ccp, 8, 8, mode, offsets, n_ovr, 1000, 1000, thres)
#     print find_poly(offs, ccp, diffpar, mode, coffs, coffsets, thres)


def correlate4(master, slave, diffpar, offs, ccp, offsets="-", coffs="-", coffsets="-"):
    if not os.path.isfile(diffpar):
        try:
            process(["create_diff_par", master + ".par", "-", diffpar, "1", "0"])
        except IOError:
            process(["create_diff_par", slave + ".par", "-", diffpar, "1", "0"])
    try:
        par = ISPPar(master+".par")
    except IOError:
        par = ISPPar(slave + ".par")
    mode = "SLC" if par.image_format in ["FCOMPLEX", "SCOMPLEX"] else "MLI"
    # determine best window size
    bestwinsize = 0
    bestnoffs = 0
    for factor in [1, 2, 3, 4, 8]:
        try:
            noffs = cc_find2(master, slave, diffpar, offs, ccp, 8, 8, mode, offsets, 2, 1000/factor, 1000/factor)
            if noffs > bestnoffs:
                bestnoffs = noffs
                bestwinsize = factor*8
            else:
                break
        except RuntimeError:
            continue
    if bestwinsize == 0:
        bestwinsize = 8
    print "optimal window size:", bestwinsize
    if bestnoffs < 1000:
        print "initial number of offsets:", int(bestnoffs)
        # offset search refinement
        nr = (par.range_samples // 1000 * 1000)/(bestwinsize/8)
        naz = (par.azimuth_lines // 1000 * 1000)/(bestwinsize/8)
        nr = nr if nr <= 4000 else 4000
        naz = naz if naz <= 4000 else 4000
        print "new number of windows:", nr, naz
        bestnoffs = cc_find2(master, slave, diffpar, offs, ccp, bestwinsize, bestwinsize, mode, offsets, 8, nr, naz)
        print "final number of offsets:", int(bestnoffs)
    else:
        print "number of offsets:", int(bestnoffs)
    fit = find_poly(offs, ccp, diffpar, mode, coffs, coffsets)
    print "achieved fit:", fit


def correlate5(master, slave, diffpar, offs, ccp, offsets="-", coffs="-", coffsets="-"):
    if os.path.isfile(master+".par"):
        parfile = master+".par"
    elif os.path.isfile(slave+".par"):
        parfile = slave+".par"
    else:
        raise IOError("no appropriate parameter file found")
    process(["create_diff_par", parfile, "-", diffpar, 1, 0])
    par = ISPPar(parfile)
    mode = "SLC" if par.image_format in ["FCOMPLEX", "SCOMPLEX"] else "MLI"
    winsize = int(2**(math.log(par.range_pixel_spacing*10, 2)//1))
    dim_ratio = float(par.azimuth_lines)/float(par.range_samples)
    overlap = .5
    noffsets = 0
    while winsize >= 8:
        print "window size:", winsize
        nr = int(round((float(par.range_samples) / winsize) * (1 + overlap)))
        naz = int(int(nr) * dim_ratio)
        print "number of search windows (range, azimuth):", nr, naz
        try:
            noffsets = cc_find2(master, slave, diffpar, offs, ccp, winsize, winsize, mode, offsets, 2, nr, naz)
            print "number of found offsets:", noffsets
            if noffsets > 50:
                fit = find_poly(offs, ccp, diffpar, mode, coffs, coffsets)
                print "polynomial fit (range, azimuth):", fit
                if sum(fit) > 1.5:
                    raise RuntimeError
                else:
                    break
            winsize /= 2
        except RuntimeError:
            winsize /= 2
            continue
    if winsize == 4 and noffsets < 50:
        print "window size:", 8
        nr = (par.range_samples // 1000 * 1000)
        naz = (par.azimuth_lines // 1000 * 1000)
        print "number of search windows (range, azimuth):", nr, naz
        noffsets = cc_find2(master, slave, diffpar, offs, ccp, 8, 8, mode, offsets, 8, nr, naz)
        print "number of found offsets:", noffsets
        fit = find_poly(offs, ccp, diffpar, mode, coffs, coffsets)
        print "polynomial fit (range, azimuth):", fit


def geocode(scene, dem, tempdir, outdir, targetres, scaling="linear", func_geoback=2, func_interp=0, sarsimulation=True):
    """
    scaling can be either 'linear', 'db' or a list of both (i.e. ['linear', 'db'])

    intermediate output files (named <master_MLI>_<suffix>):
    dem: dem subsetted to the extent of the SAR image
    lut: rough geocoding lookup table
    lut_fine: fine geocoding lookup table
    sim_map: simulated SAR backscatter image in DEM geometry
    sim_sar: simulated SAR backscatter image in SAR geometry
    u: zenith angle of surface normal vector n (angle between z and n)
    v: orientation angle of n (between x and projection of n in xy plane)
    inc: local incidence angle (between surface normal and look vector)
    psi: projection angle (between surface normal and image plane normal)
    pix: pixel area normalization factor
    ls_map: layover and shadow map (in map projection)
    diffpar: ISP offset/interferogram parameter file
    offs: offset estimates (fcomplex)
    coffs: culled range and azimuth offset estimates (fcomplex)
    coffsets: culled offset estimates and cross correlation values (text format)
    ccp: cross-correlation of each patch (0.0->1.0) (float)
    """

    scene = scene if isinstance(scene, pyroSAR.ID) else pyroSAR.identify(scene)

    scaling = [scaling] if isinstance(scaling, str) else scaling if isinstance(scaling, list) else []
    scaling = union(scaling, ["db", "linear"])
    if len(scaling) == 0:
        raise IOError("wrong input type for parameter scaling")

    scene.unpack(tempdir)

    scene.convert2gamma(scene.scene)

    try:
        scene.correctOSV()
    except:
        pass

    scene.calibrate()

    images = [x for x in scene.getGammaImages(scene.scene) if x.endswith("_grd")]

    # rlks = int(targetres//scene.spacing[0])
    # azlks = int(targetres//scene.spacing[1])
    # rlks = 1 if rlks == 0 else rlks
    # azlks = 1 if azlks == 0 else azlks

    # rlks = int(targetres/round(scene.spacing[0]))
    # azlks = int(targetres/round(scene.spacing[1]))
    #
    # for image in images:
    #     gamma(["multi_look_MLI", image, image+".par", image+"_mli", image+"_mli.par", rlks, azlks])

    for image in images:
        multilook(image, image+"_mli", targetres)

    images = [x+"_mli" for x in images]

    master = images[0]

    dem_seg = master+"_dem"
    lut_coarse = master+"_lut_coarse"
    lut_fine = master+"_lut_fine"
    sim_map = master+"_sim_map"
    sim_sar = master+"_sim_sar"
    u = master+"_u"
    v = master+"_v"
    inc = master+"_inc"
    psi = master+"_psi"
    pix = master+"_pix"
    ls_map = master+"_ls_map"
    pixel_area_coarse = master+"_pixel_area_coarse"
    pixel_area_fine = master+"_pixel_area_fine"
    diffpar = master+"_off.par"
    offs = master+"_offs"
    offsets = offs + ".txt"
    coffs = master+"_coffs"
    coffsets = master+"_coffsets"
    ccp = master+"_ccp"
    ellipse_pixel_area = master+"_ellipse_pixel_area"
    ratio_sigma0 = master+"_ratio_sigma0"

    ovs_lat, ovs_lon = ovs(dem+".par", targetres)

    path_log = os.path.join(scene.scene, "logfiles")
    if not os.path.isdir(path_log):
        os.makedirs(path_log)

    master_par = ISPPar(master+".par")

    gc_map_opt = [dem+".par", dem, dem_seg+".par", dem_seg, lut_coarse, ovs_lat, ovs_lon, sim_map, u, v, inc, psi, pix, ls_map, 8, func_interp]

    if master_par.image_geometry == "GROUND_RANGE":
        process(["gc_map_grd", master + ".par"] + gc_map_opt, logpath=path_log)
    else:
        process(["gc_map", master + ".par", "-"] + gc_map_opt, logpath=path_log)

    for item in [dem_seg, sim_map, u, v, psi, pix, inc]:
        envi.hdr(dem_seg+".par", item+".hdr")

    sim_width = ISPPar(dem_seg+".par").width

    if sarsimulation is True:
        ######################################################################
        # cross correlation approach 1 #######################################
        # ####################################################################
        # dempar = ISPPar(dem_seg + ".par")
        # samples_dem = dempar.width
        # samples_mli = master_par.range_samples
        # lines_mli = master_par.azimuth_lines
        # gamma(["geocode", lut, sim_map, samples_dem, sim_sar, samples_mli, lines_mli, func_interp], logpath=path_log)
        # init_offset(master, sim_sar, diffpar, path_log)
        # correlate(master, sim_sar, diffpar, offs, ccp, coffs=coffs, coffsets=coffsets, offsets=offs+".txt", path_log=path_log, maxwin=256)
        ######################################################################
        # cross correlation approach 2 #######################################
        ######################################################################
        process(["pixel_area", master + ".par", dem_seg + ".par", dem_seg, lut_coarse, ls_map, inc, pixel_area_coarse], logpath=path_log)
        init_offset(master, pixel_area_coarse, diffpar, path_log)
        # correlate(master, pixel_area_coarse, diffpar, offs, ccp, coffs=coffs, coffsets=coffsets, offsets=offsets, path_log=path_log, maxwin=256)
        # correlate2(master, pixel_area_coarse, diffpar, offs, ccp, offsets=offsets, coffs=coffs, coffsets=coffsets, path_log=path_log)
        # correlate3(master, pixel_area_coarse, diffpar, offs, ccp, offsets=offsets, coffs=coffs, coffsets=coffsets)
        print "#####################################"
        print "cross correlation step 1"
        correlate5(pixel_area_coarse, master, diffpar, offs, ccp, offsets, coffs, coffsets)
        print "#####################################"
        print "cross correlation step 2"
        correlate5(pixel_area_coarse, master, diffpar, offs, ccp, offsets, coffs, coffsets)
        print "#####################################"
        # print "cross correlation step 3"
        # correlate4(pixel_area_coarse, master, diffpar, offs, ccp, offsets, coffs, coffsets)
        # print "#####################################"
        ######################################################################
        try:
            process(["gc_map_fine", lut_coarse, sim_width, diffpar, lut_fine, 0], logpath=path_log)
        except sp.CalledProcessError:
            print "...failed"
            return
    else:

        lut_fine = lut_coarse

    ######################################################################
    # normalization and backward geocoding approach 1 ####################
    ######################################################################
    for image in images:
        process(["geocode_back", image, master_par.range_samples, lut_fine, image + "_geo", sim_width, "-", func_geoback], logpath=path_log)
        process(["product", image + "_geo", pix, image + "_geo_pan", sim_width, 1, 1, 0], logpath=path_log)
        process(["lin_comb", 1, image + "_geo_pan", 0, math.cos(math.radians(master_par.incidence_angle)), image + "_geo_pan_flat", sim_width], logpath=path_log)
        process(["sigma2gamma", image + "_geo_pan_flat", inc, image + "_geo_norm", sim_width], logpath=path_log)
        envi.hdr(dem_seg+".par", image+"_geo_norm.hdr")
    ######################################################################
    # normalization and backward geocoding approach 2 ####################
    ######################################################################
    # try:
    #     gamma(["pixel_area", master+".par", dem_seg+".par", dem_seg, lut_fine, ls_map, inc, pixel_area_fine], logpath=path_log)
    # except sp.CalledProcessError:
    #     print "...failed"
    #     return
    # gamma(["radcal_MLI", master, master+".par", "-", master+"_cal", "-", 0, 0, 1, 0.0, "-", ellipse_pixel_area], logpath=path_log)
    # gamma(["ratio", ellipse_pixel_area, pixel_area_fine, ratio_sigma0, master_par.range_samples, 1, 1], logpath=path_log)
    #
    # for image in images:
    #     gamma(["product", image, ratio_sigma0, image+"_pan", master_par.range_samples, 1, 1], logpath=path_log)
    #     gamma(["geocode_back", image+"_pan", master_par.range_samples, lut_fine, image+"_pan_geo", sim_width, 0, func_geoback], logpath=path_log)
    #     gamma(["lin_comb", 1, image+"_pan_geo", 0, math.cos(math.radians(master_par.incidence_angle)), image+"_pan_geo_flat", sim_width], logpath=path_log)
    #     gamma(["sigma2gamma", image+"_pan_geo_flat", inc, image+"_geo_norm", sim_width], logpath=path_log)
    #     envi.hdr(dem_seg+".par", image+"_geo_norm.hdr")
    ######################################################################
    # conversion to (dB and) geotiff
    for image in images:
        for scale in scaling:
            if scale == "db":
                process(["linear_to_dB", image + "_geo_norm", image + "_geo_norm_db", sim_width, 0, -99], logpath=path_log)
                envi.hdr(dem_seg+".par", image+"_geo_norm_db.hdr")
            suffix = {"linear": "", "db": "_db"}[scale]
            infile = image+"_geo_norm{}".format(suffix)
            outfile = os.path.join(outdir, os.path.basename(image)+"_geo_norm{}.tif".format(suffix))
            try:
                process(["data2geotiff", dem_seg + ".par", infile, 2, outfile], logpath=path_log)
            except ImportWarning:
                pass


def init_offset(master, slave, off, path_log, thres=0.15):
    """
    wrapper for iterated initial offset estimation between two SLCs or MLIs
    from a starting search window size of 128x128, initial offsets are search for; if this fails, the window size is increased up to a maximum size of 1024x1024
    once a sufficient number of offsets is found, the window size is again decreased down to 128 and offset search repeated to refine the results on finer levels
    in case too few/no offsets were found on any window size, the routine will throw an error
    additionally, for each window size, the offset search is performed on three different levels of multilooking (factor 3, 2, 1) as this is reported to refine the results;
    if the routine fails on any level of multilooking, the routine will proceed with the next window size
    """
    path_out = os.path.dirname(off)

    if not os.path.isfile(off):
        process(["create_diff_par", master + ".par", "-", off, 1, 0], logpath=path_log)

    par = ISPPar(master + ".par")
    mode = "SLC" if par.image_format in ["FCOMPLEX", "SCOMPLEX"] else "MLI"

    # first offset estimation using orbit state-vectors and image parameters (most important in case of very large offsets)
    # this is only to be performed while coregistering two scenes not while geocoding
    # todo: query input file type in order to decide whether command is executed
    cmd = "init_offset_orbit" if mode == "SLC" else "init_offset_orbitm"
    if mode == "SLC":
        process([cmd, master + ".par", slave + ".par", off], path_out, path_log)

    mlk = Spacing(par)
    passed = False
    patchsizes = [128, 256, 512, 1024]
    patchstack = Stack(patchsizes[0])

    while not passed and not patchstack.empty():
        win = patchstack.pop()
        for factor in [3, 2, 1]:
            try:
                args = [int(mlk.rlks)*factor, int(mlk.azlks)*factor, "-", "-", "-", "-", thres, win, win]
                if mode == "SLC":
                    process(["init_offset", master, slave, master + ".par", slave + ".par", off] + args, path_out, path_log)
                else:
                    process(["init_offsetm", master, slave, off] + args, path_out, path_log)
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


def ovs(parfile, targetres):
    """
    compute DEM oversampling factors for a target resolution in meters
    """
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


def multilook(infile, outfile, targetres):
    par = ISPPar(infile+".par")

    if par.image_geometry == "SLANT_RANGE":
        groundRangePS = par.range_pixel_spacing/(math.sin(math.radians(par.incidence_angle)))
        rlks = int(round(float(targetres)/groundRangePS))
    else:
        rlks = int(round(float(targetres)/par.range_pixel_spacing))
    azlks = int(round(float(targetres)/par.azimuth_pixel_spacing))

    rlks = rlks if rlks > 0 else 1
    azlks = azlks if azlks > 0 else 1

    if rlks+azlks == 2:
        print "nothing to be done"
        return
    else:
        cmd = ["multi_look_MLI", infile, infile+".par", outfile, outfile+".par", str(rlks), str(azlks)]
        proc = sp.Popen(cmd, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
        out, err = proc.communicate()
        envi.hdr(outfile+".par")
