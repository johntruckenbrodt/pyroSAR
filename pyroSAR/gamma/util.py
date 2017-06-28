#!/usr/bin/env python
##############################################################
# universal core routines for processing SAR images with GAMMA
# John Truckenbrodt 2014-2017
##############################################################

"""
This script is intended as a set of generalized processing routines for modularized GAMMA work flows.
The function parametrization is intended to be applicable to any kind of situation and input data set. Thus, instead of choosing a specific parametrization for the data at hand,
core parameters are iterated over a set of values in order to find the one best suited for the task.
The approach of the single routines is likely to still have drawbacks and might fail in certain situations. Testing and suggestions on improvements are very welcome.
"""
import math
import os
import re
import shutil
import subprocess as sp
from collections import OrderedDict

from osgeo import ogr

from .. import envi
from ..spatial import haversine

from ..ancillary import run, Stack, union, finder
from . import ISPPar, Spacing, Namespace
from .error import gammaErrorHandler

ogr.UseExceptions()


def process(cmd, outdir=None, logpath=None, inlist=None, void=True):
    log = os.path.join(logpath, cmd[0] + '.log') if logpath else None
    out, err = run(cmd, outdir=outdir, logfile=log, inlist=inlist, void=False, errorpass=True)
    gammaErrorHandler(out, err)
    if not void:
        return out, err


def slc_corners(parfile):
    """
    extract the corner soordinates of a SAR scene
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


def slc_burst_corners(parfile, tops_parfile):
    proc = sp.Popen(['SLC_burst_corners', parfile, tops_parfile], stdin=sp.PIPE,
                    stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(err + '\nreading of burst corners failed')
    matches = [re.findall('[0-9.]+', x) for x in re.findall('Burst:[0-9 .]*', out)]
    bursts = {}
    for line in matches:
        lat = map(float, line[1::2])
        lon = map(float, line[2::2])
        bursts[int(line[0])] = zip(lat, lon)
    return OrderedDict(sorted(bursts.items()))


def coord2geom(coordinatelist):
    coordstrings = [' '.join(map(str, x)) for x in coordinatelist]
    coordstrings.append(coordstrings[0])
    return ogr.CreateGeometryFromWkt('POLYGON(({}))'.format(', '.join(coordstrings)))


def burstoverlap(slc1_parfile, slc2_parfile, slc1_tops_parfile, slc2_tops_parfile):
    b1 = slc_burst_corners(slc1_parfile, slc1_tops_parfile)
    b2 = slc_burst_corners(slc2_parfile, slc2_tops_parfile)
    master_select = []
    slave_select = []
    for mburst in b1.keys():
        for sburst in b2.keys():
            poly1 = coord2geom(b1[mburst])
            poly2 = coord2geom(b2[sburst])
            intersection = poly1.Intersection(poly2)
            area_ratio = intersection.GetArea() / poly1.GetArea()
            if area_ratio > 0.9:
                master_select.append(mburst)
                slave_select.append(sburst)
    return master_select, slave_select


def correlate(master, slave, off, offs, ccp, offsets='-', coffs='-', coffsets='-',
              path_log=None, maxwin=2048, minwin=16, overlap=.3, poly=4, ovs=2,
              thres=0.15):
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
    the default value '-' will result in no output file written
    """

    path_out = os.path.dirname(off)
    par = ISPPar(master + '.par')

    if not os.path.isfile(off):
        process(['create_diff_par', master + '.par', '-', off, 1, 0], logpath=path_log)

    if path_log is not None:
        if not os.path.isdir(path_log):
            os.makedirs(path_log)

    if par.image_format in ['FCOMPLEX', 'SCOMPLEX']:
        commands = ['offset_pwr', 'offset_fit']
        mode = 'SLC'
    else:
        commands = ['offset_pwrm', 'offset_fitm']
        mode = 'MLI'

    # compute the number of estimation windows in azimuth from the defined number of range windows
    dim_ratio = float(par.azimuth_lines) / float(par.range_samples)

    # iteratively reduce the size of the search windows until a sufficient number of offsets was found and/or a final window size is reached
    passed = False
    winsize = maxwin
    while winsize >= minwin:
        # compute the number of windows needed in range and azimuth based on the number of image pixels in both directions
        nr = int(round((float(par.range_samples) / winsize) * (1 + overlap)))
        naz = str(int(int(nr) * dim_ratio))

        args = [off, offs, ccp, winsize, winsize, offsets, ovs, nr, naz, thres]
        try:
            if mode == 'SLC':
                process(
                    [commands[0], master, slave, master + '.par', slave + '.par'] + args,
                    path_out, path_log)
            else:
                process([commands[0], master, slave] + args, path_out, path_log)
            passed = True
            try:
                process([commands[1], offs, ccp, off, coffs, coffsets, thres, poly],
                        path_out, path_log)
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
        raise RuntimeError(
            'cross-correlation failed; consider verifying scene overlap or choice of polarization')


def cc_find(master, slave, diffpar, offs, ccp, winsize, mode, offsets='-', n_ovr=2,
            nr=400, naz=400, thres=.15, c_ovr=4):
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
    args = [diffpar, offs, ccp, str(winsize), str(winsize), offsets, str(n_ovr), str(nr),
            str(naz), str(thres), str(c_ovr)]

    if mode == 'MLI':
        cmd = ['offset_pwrm', master, slave] + args
    elif mode == 'SLC':
        cmd = ['offset_pwr', master, slave, master + '.par', slave + '.par'] + args
    else:
        raise ValueError('mode must be either "MLI" or "SLC"')

    proc = sp.Popen(cmd, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = proc.communicate()
    for line in out.split('\n'):
        if line.startswith('number of offsets above correlation threshold'):
            vals = map(float, re.findall('[0-9]+', line))
            return vals[0] / vals[1]


def cc_fit(offs, ccp, diffpar, mode, coffs='-', coffsets='-', thres=.15, npoly=4):
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

    if mode == 'MLI':
        cmd = ['offset_fitm'] + args
    elif mode == 'SLC':
        cmd = ['offset_fit'] + args
    else:
        raise ValueError('mode must be either "MLI" or "SLC"')

    proc = sp.Popen(cmd, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = proc.communicate()
    fit = None
    for line in out.split('\n'):
        if line.startswith('final model fit std. dev.'):
            fit = tuple(map(float, re.findall('[0-9]+.[0-9]+', line)))
    print 'final model fit std. dev. (range, azimuth):', fit
    if not fit:
        raise RuntimeError(
            'failed at fitting offset estimation; command:\n{}'.format(' '.join(cmd)))


def cc_fit2(offs, ccp, diffpar, mode, coffs='-', coffsets='-', thres=.15, npoly=4):
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

    if mode == 'MLI':
        cmd = ['offset_fitm'] + args
    elif mode == 'SLC':
        cmd = ['offset_fit'] + args
    else:
        raise ValueError('mode must be either "MLI" or "SLC"')

    proc = sp.Popen(cmd, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = proc.communicate()
    errormessage = 'ERROR: insufficient offset points left after culling to determine offset model parameters'
    if errormessage in err:
        raise RuntimeError(errormessage)
    for line in out.split('\n'):
        if line.startswith('final model fit std. dev.'):
            return tuple(map(float, re.findall('[0-9]+.[0-9]+', line)))
    raise RuntimeError(err)


def cc(master, slave, diffpar, offs, ccp, mode, offsets='-', nr=400, naz=400, thres=.15,
       c_ovr=4, minwin=8, maxwin=1024):
    best = 0
    winsize = minwin
    bestwin = 8
    while winsize <= maxwin:
        ratio = cc_find(master, slave, diffpar, offs, ccp, winsize, mode, offsets, 1, nr,
                        naz, thres, c_ovr)
        if ratio > best:
            best = ratio
            bestwin = winsize
        else:
            break
        winsize *= 2
    print 'optimal window size:', bestwin
    ratio = cc_find(master, slave, diffpar, offs, ccp, bestwin, mode, n_ovr=8)
    print 'offset acceptance ratio:', ratio
    return bestwin


def cc_find2(master, slave, diffpar, offs, ccp, rwin, azwin, mode, offsets='-', n_ovr=2,
             nr=400, naz=400, thres=.15, c_ovr=4):
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
    args = [diffpar, offs, ccp, str(rwin), str(azwin), offsets, str(n_ovr), str(nr),
            str(naz), str(thres), str(c_ovr)]
    if mode == 'MLI':
        cmd = ['offset_pwrm', master, slave] + args
    elif mode == 'SLC':
        cmd = ['offset_pwr', master, slave, master + '.par', slave + '.par'] + args
    else:
        raise ValueError('mode must be either "MLI" or "SLC"')
    proc = sp.Popen(cmd, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = proc.communicate()
    errormessage = 'ERROR: no offsets found above correlation threshold'
    if errormessage in err:
        raise RuntimeError(errormessage)
    for line in out.split('\n'):
        if line.startswith('number of offsets above correlation threshold'):
            vals = map(float, re.findall('[0-9]+', line))
            return int(vals[0])
    raise RuntimeError(err)


def correlate2(master, slave, diffpar, offs, ccp, offsets='-', coffs='-', coffsets='-',
               npoly=4, thres=0.15, minwin=8, maxwin=2048, niter=2, path_log=None):
    if not os.path.isfile(diffpar):
        process(['create_diff_par', master + '.par', '-', diffpar, 1, 0],
                logpath=path_log)
    if path_log is not None:
        if not os.path.isdir(path_log):
            os.makedirs(path_log)
    mode = 'SLC' if ISPPar(master + '.par').image_format in ['FCOMPLEX',
                                                             'SCOMPLEX'] else 'MLI'
    print '#####################################\ncross-correlation offset estimation'
    for i in range(0, niter):
        print '#############\niteration {}'.format(i + 1)
        cc(master, slave, diffpar, offs, ccp, mode, offsets, thres=thres, minwin=minwin,
           maxwin=maxwin)
        cc_fit(offs, ccp, diffpar, mode, coffs, coffsets, thres, npoly)
    print '#############\n...done\n#####################################'


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
    print 'best polynomial order:', bestpoly
    return cc_fit2(offs, ccp, diffpar, mode, coffs, coffsets, thres, bestpoly)


# def correlate3(master, slave, diffpar, offs, ccp, offsets='-', coffs='-', coffsets='-'):
#     if not os.path.isfile(diffpar):
#         process(['create_diff_par', master + '.par', '-', diffpar, 1, 0], logpath=path_log)
#     par = ISPPar(master+'.par')
#     mode = 'SLC' if par.image_format in ['FCOMPLEX', 'SCOMPLEX'] else 'MLI'
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


def correlate4(master, slave, diffpar, offs, ccp, offsets='-', coffs='-', coffsets='-'):
    if not os.path.isfile(diffpar):
        try:
            process(['create_diff_par', master + '.par', '-', diffpar, '1', '0'])
        except IOError:
            process(['create_diff_par', slave + '.par', '-', diffpar, '1', '0'])
    try:
        par = ISPPar(master + '.par')
    except IOError:
        par = ISPPar(slave + '.par')
    mode = 'SLC' if par.image_format in ['FCOMPLEX', 'SCOMPLEX'] else 'MLI'
    # determine best window size
    bestwinsize = 0
    bestnoffs = 0
    for factor in [1, 2, 3, 4, 8]:
        try:
            noffs = cc_find2(master, slave, diffpar, offs, ccp, 8, 8, mode, offsets, 2,
                             1000 / factor, 1000 / factor)
            if noffs > bestnoffs:
                bestnoffs = noffs
                bestwinsize = factor * 8
            else:
                break
        except RuntimeError:
            continue
    if bestwinsize == 0:
        bestwinsize = 8
    print 'optimal window size:', bestwinsize
    if bestnoffs < 1000:
        print 'initial number of offsets:', int(bestnoffs)
        # offset search refinement
        nr = (par.range_samples // 1000 * 1000) / (bestwinsize / 8)
        naz = (par.azimuth_lines // 1000 * 1000) / (bestwinsize / 8)
        nr = nr if nr <= 4000 else 4000
        naz = naz if naz <= 4000 else 4000
        print 'new number of windows:', nr, naz
        bestnoffs = cc_find2(master, slave, diffpar, offs, ccp, bestwinsize, bestwinsize,
                             mode, offsets, 8, nr, naz)
        print 'final number of offsets:', int(bestnoffs)
    else:
        print 'number of offsets:', int(bestnoffs)
    fit = find_poly(offs, ccp, diffpar, mode, coffs, coffsets)
    print 'achieved fit:', fit


def correlate5(master, slave, diffpar, offs, ccp, offsets='-', coffs='-', coffsets='-'):
    if os.path.isfile(master + '.par'):
        parfile = master + '.par'
    elif os.path.isfile(slave + '.par'):
        parfile = slave + '.par'
    else:
        raise IOError('no appropriate parameter file found')
    process(['create_diff_par', parfile, '-', diffpar, 1, 0])
    par = ISPPar(parfile)
    mode = 'SLC' if par.image_format in ['FCOMPLEX', 'SCOMPLEX'] else 'MLI'
    winsize = int(2 ** (math.log(par.range_pixel_spacing * 10, 2) // 1))
    dim_ratio = float(par.azimuth_lines) / float(par.range_samples)
    overlap = .5
    noffsets = 0
    while winsize >= 8:
        print 'window size:', winsize
        nr = int(round((float(par.range_samples) / winsize) * (1 + overlap)))
        naz = int(int(nr) * dim_ratio)
        print 'number of search windows (range, azimuth):', nr, naz
        try:
            noffsets = cc_find2(master, slave, diffpar, offs, ccp, winsize, winsize, mode,
                                offsets, 2, nr, naz)
            print 'number of found offsets:', noffsets
            if noffsets > 50:
                fit = find_poly(offs, ccp, diffpar, mode, coffs, coffsets)
                print 'polynomial fit (range, azimuth):', fit
                if sum(fit) > 1.5:
                    raise RuntimeError
                else:
                    break
            winsize /= 2
        except RuntimeError:
            winsize /= 2
            continue
    if winsize == 4 and noffsets < 50:
        print 'window size:', 8
        nr = (par.range_samples // 1000 * 1000)
        naz = (par.azimuth_lines // 1000 * 1000)
        print 'number of search windows (range, azimuth):', nr, naz
        noffsets = cc_find2(master, slave, diffpar, offs, ccp, 8, 8, mode, offsets, 8, nr,
                            naz)
        print 'number of found offsets:', noffsets
        fit = find_poly(offs, ccp, diffpar, mode, coffs, coffsets)
        print 'polynomial fit (range, azimuth):', fit


def geocode(scene, dem, tempdir, outdir, targetres, scaling='linear', func_geoback=2,
            func_interp=0, nodata=(0, -99), sarsimulation=True, osvdir=None, cleanup=True):
    """
    scene: the SAR scene to be processed (can be the name of the file or a pyroSAR object)
    dem: the reference DEM in GAMMA format
    tempdir: a temporary directory for writing intermediate files
    outdir: the directory for the final GeoTiff output files
    targetres: the target resolution in meters
    scaling: can be either 'linear', 'db' or a list of both (i.e. ['linear', 'db'])
    func_geoback: backward geocoding interpolation mode (see GAMMA command geocode_back)
        0: nearest-neighbor
        1: bicubic spline
        2: bicubic-log spline, interpolates log(data)
        3: bicubic-sqrt spline, interpolates sqrt(data)
        NOTE: bicubic-log spline and bicubic-sqrt spline modes should only be used with non-negative data!
    func_interp: output lookup table values in regions of layover, shadow, or DEM gaps (see GAMMA command gc_map)
        0: set to (0.,0.)
        1: linear interpolation across these regions
        2: actual value
        3: nn-thinned
    nodata: the nodata values for the output files; defined as a tuple with two values, the first for linear, the second for logarithmic scaling, per default (0, -99)
    sarsimulation: perform geocoding with SAR simulation cross correlation? If False, geocoding is performed with the Range-Doppler approach using orbit state vectors
    cleanup: should all files written to the temporary directory during function execution be deleted after processing?

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
    if sarsimulation is True:
        raise RuntimeError('geocoding with cross correlation offset refinement is currently disabled')

    scene = scene if isinstance(scene, pyroSAR.ID) else pyroSAR.identify(scene)

    for dir in [tempdir, outdir]:
        if not os.path.isdir(dir):
            os.makedirs(dir)

    if scene.is_processed(outdir):
        print('scene {} already processed'.format(scene.outname_base()))
        return

    scaling = [scaling] if isinstance(scaling, str) else scaling if isinstance(scaling, list) else []
    scaling = union(scaling, ['db', 'linear'])
    if len(scaling) == 0:
        raise IOError('wrong input type for parameter scaling')

    if scene.compression is not None:
        print 'unpacking scene..'
        scene.unpack(tempdir)
    else:
        scene.scene = os.path.join(tempdir, os.path.basename(scene.file))
        os.makedirs(scene.scene)

    if scene.sensor in ['S1A', 'S1B']:
        print 'removing border noise..'
        scene.removeGRDBorderNoise()

    print 'converting scene to GAMMA format..'
    scene.convert2gamma(scene.scene)

    if scene.sensor in ['S1A', 'S1B']:
        print 'updating orbit state vectors..'
        try:
            scene.correctOSV(osvdir)
        except RuntimeError:
            return

    scene.calibrate()

    images = [x for x in scene.getGammaImages(scene.scene) if x.endswith('_grd') or x.endswith('_slc_cal')]

    for image in images:
        multilook(image, image + '_mli', targetres)

    images = [x + '_mli' for x in images]

    master = images[0]

    n = Namespace(scene.scene, scene.outname_base())
    n.appreciate(['dem_seg', 'lut_coarse', 'lut_fine', 'pix', 'ccp', 'inc', 'ls_map'])
    n.depreciate(['sim_map', 'u', 'v', 'psi'])

    if sarsimulation:
        n.appreciate(['ls_map'])

    ovs_lat, ovs_lon = ovs(dem + '.par', targetres)

    path_log = os.path.join(scene.scene, 'logfiles')
    if not os.path.isdir(path_log):
        os.makedirs(path_log)

    master_par = ISPPar(master + '.par')

    gc_map_args = [dem + '.par', dem, n.dem_seg + '.par', n.dem_seg, n.lut_coarse,
                   ovs_lat, ovs_lon, n.sim_map, n.u, n.v, n.inc, n.psi, n.pix, n.ls_map,
                   8, func_interp]

    if master_par.image_geometry == 'GROUND_RANGE':
        process(['gc_map_grd', master + '.par'] + gc_map_args, logpath=path_log)
    else:
        process(['gc_map', master + '.par', '-'] + gc_map_args, logpath=path_log)

    for item in ['dem_seg', 'sim_map', 'u', 'v', 'psi', 'pix', 'inc']:
        if n.isappreciated(item):
            envi.hdr(n.dem_seg + '.par', n.get(item) + '.hdr')

    sim_width = ISPPar(n.dem_seg + '.par').width

    if sarsimulation is True:

        n.appreciate(['lut_fine', 'diff.par', 'offs', 'ccp'])
        n.depreciate(['offsets', 'coffs', 'coffsets'])
        ######################################################################
        # cross correlation approach 1 #######################################
        # ####################################################################
        # dempar = ISPPar(dem_seg + '.par')
        # samples_dem = dempar.width
        # samples_mli = master_par.range_samples
        # lines_mli = master_par.azimuth_lines
        # gamma(['geocode', lut, sim_map, samples_dem, sim_sar, samples_mli, lines_mli, func_interp], logpath=path_log)
        # init_offset(master, sim_sar, diffpar, path_log)
        # correlate(master, sim_sar, diffpar, offs, ccp, coffs=coffs, coffsets=coffsets, offsets=offs+'.txt', path_log=path_log, maxwin=256)
        ######################################################################
        # cross correlation approach 2 #######################################
        ######################################################################
        n.appreciate(['pixel_area_coarse'])
        process(
            ['pixel_area', master + '.par', n.dem_seg + '.par', n.dem_seg, n.lut_coarse, n.ls_map, n.inc, n.pixel_area_coarse], logpath=path_log)
        init_offset(master, n.pixel_area_coarse, n.diff_par, path_log)
        # correlate(master, pixel_area_coarse, diffpar, offs, ccp, coffs=coffs, coffsets=coffsets, offsets=offsets, path_log=path_log, maxwin=256)
        # correlate2(master, pixel_area_coarse, diffpar, offs, ccp, offsets=offsets, coffs=coffs, coffsets=coffsets, path_log=path_log)
        # correlate3(master, pixel_area_coarse, diffpar, offs, ccp, offsets=offsets, coffs=coffs, coffsets=coffsets)
        print '#####################################'
        print 'cross correlation step 1'
        correlate5(n.pixel_area_coarse, master, n.diff_par, n.offs, n.ccp, n.offsets, n.coffs, n.coffsets)
        print '#####################################'
        print 'cross correlation step 2'
        correlate5(n.pixel_area_coarse, master, n.diff_par, n.offs, n.ccp, n.offsets, n.coffs, n.coffsets)
        print '#####################################'
        # print 'cross correlation step 3'
        # correlate4(pixel_area_coarse, master, diffpar, offs, ccp, offsets, coffs, coffsets)
        # print '#####################################'
        ######################################################################
        try:
            process(['gc_map_fine', n.lut_coarse, sim_width, n.diff_par, n.lut_fine, 0], logpath=path_log)
        except sp.CalledProcessError:
            print '...failed'
            return
        lut_final = n.lut_fine
    else:
        lut_final = n.lut_coarse

    ######################################################################
    # normalization and backward geocoding approach 1 ####################
    ######################################################################
    for image in images:
        process(['geocode_back', image, master_par.range_samples, lut_final, image + '_geo', sim_width, '-', func_geoback], logpath=path_log)
        process(['product', image + '_geo', n.pix, image + '_geo_pan', sim_width, 1, 1, 0], logpath=path_log)
        process(['lin_comb', 1, image + '_geo_pan', 0, math.cos(math.radians(master_par.incidence_angle)), image + '_geo_pan_flat', sim_width], logpath=path_log)
        process(['sigma2gamma', image + '_geo_pan_flat', n.inc, image + '_geo_norm', sim_width], logpath=path_log)
        envi.hdr(n.dem_seg + '.par', image + '_geo_norm.hdr')
    ######################################################################
    # normalization and backward geocoding approach 2 ####################
    ######################################################################
    # process(['pixel_area', master+'.par', dem_seg+'.par', dem_seg, lut_fine, ls_map, inc, pixel_area_fine], logpath=path_log)
    # process(['radcal_MLI', master, master+'.par', '-', master+'_cal', '-', 0, 0, 1, 0.0, '-', ellipse_pixel_area], logpath=path_log)
    # process(['ratio', ellipse_pixel_area, pixel_area_fine, ratio_sigma0, master_par.range_samples, 1, 1], logpath=path_log)
    #
    # for image in images:
    #     process(['product', image, ratio_sigma0, image+'_pan', master_par.range_samples, 1, 1], logpath=path_log)
    #     process(['geocode_back', image+'_pan', master_par.range_samples, lut_fine, image+'_pan_geo', sim_width, 0, func_geoback], logpath=path_log)
    #     process(['lin_comb', 1, image+'_pan_geo', 0, math.cos(math.radians(master_par.incidence_angle)), image+'_pan_geo_flat', sim_width], logpath=path_log)
    #     process(['sigma2gamma', image+'_pan_geo_flat', inc, image+'_geo_norm', sim_width], logpath=path_log)
    #     envi.hdr(dem_seg+'.par', image+'_geo_norm.hdr')
    ######################################################################
    # conversion to (dB and) geotiff
    for image in images:
        for scale in scaling:
            if scale == 'db':
                process(['linear_to_dB', image + '_geo_norm', image + '_geo_norm_db', sim_width, 0, -99], logpath=path_log)
                envi.hdr(n.dem_seg + '.par', image + '_geo_norm_db.hdr')
                nodata_out = nodata[1]
            else:
                nodata_out = nodata[0]
            suffix = {'linear': '', 'db': '_db'}[scale]
            infile = image + '_geo_norm{}'.format(suffix)
            outfile = os.path.join(outdir, os.path.basename(image) + '_geo_norm{}.tif'.format(suffix))

            process(['data2geotiff', n.dem_seg + '.par', infile, 2, outfile, nodata_out], logpath=path_log)

    if scene.sensor in ['S1A', 'S1B']:
        shutil.copyfile(os.path.join(scene.scene, 'manifest.safe'),
                        os.path.join(outdir, scene.outname_base() + '_manifest.safe'))
    if cleanup:
        shutil.rmtree(scene.scene)


def geocode2(scene, dem, tempdir, outdir, targetres, scaling='linear', func_geoback=2,
             func_interp=0, sarsimulation=True, cleanup=True):
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

    func_geoback: backward geocoding interpolation mode (see GAMMA command geocode_back)
        0: nearest-neighbor
        1: bicubic spline
        2: bicubic-log spline, interpolates log(data)
        3: bicubic-sqrt spline, interpolates sqrt(data)
        NOTE: bicubic-log spline and bicubic-sqrt spline modes should only be used with non-negative data!

    func_interp: output lookup table values in regions of layover, shadow, or DEM gaps (enter '-' for default) (see GAMMA command gc_map)
        0: set to (0.,0.)
        1: linear interpolation across these regions (default)
        2: actual value
        3: nn-thinned
    """
    ######################################################################################
    # general setup and checks
    if sarsimulation is True:
        raise RuntimeError(
            'geocoding with cross correlation offset refinement is currently disabled')

    scaling = [scaling] if isinstance(scaling, str) \
        else scaling if isinstance(scaling, list) \
        else []
    scaling = union(scaling, ['db', 'linear'])
    if len(scaling) == 0:
        raise IOError('wrong input type for parameter scaling')

    path_log = os.path.join(scene.scene, 'logfiles')
    if not os.path.isdir(path_log):
        os.makedirs(path_log)

    ######################################################################################
    # load scene into pyroSAR and check whether it has already been processed
    scene = scene if isinstance(scene, pyroSAR.ID) else pyroSAR.identify(scene)

    if len(finder(outdir, [scene.outname_base()], regex=True, recursive=False)) != 0:
        print('scene {} already processed'.format(scene.outname_base()))
        return
    ######################################################################################
    # general preparation of the scene for processing:
    # unpacking, conversion to GAMMA format, update of orbit state vectors, calibration,
    # multilooking to approximate target resolution
    scene.unpack(tempdir)

    scene.convert2gamma(scene.scene)

    if scene.sensor in ['S1A', 'S1B']:
        scene.correctOSV()

    scene.calibrate()

    images = [x for x in scene.getGammaImages(scene.scene) if x.endswith('_grd')]

    for image in images:
        multilook(image, image + '_mli', targetres)

    images = [x + '_mli' for x in images]
    ######################################################################################
    # set up some general parameters and environments for the processing
    master = images[0]

    n = Namespace(scene.scene, scene.outname_base())
    n.appreciate(['dem_seg', 'lut_coarse', 'lut_fine', 'ccp', 'inc', 'ls_map',
                  'pix_dem', 'pix_ell'])
    n.depreciate(['sim_map', 'u', 'v', 'psi', 'pix'])

    ovs_lat, ovs_lon = ovs(dem + '.par', targetres)
    master_par = ISPPar(master + '.par')
    ######################################################################################
    # Create geocoding lookup table for reference image
    gc_map_args = [dem + '.par', dem, n.dem_seg + '.par', n.dem_seg, n.lut_coarse,
                   ovs_lat, ovs_lon, n.sim_map, n.u, n.v, n.inc, n.psi, n.pix,
                   n.ls_map, 8, func_interp]

    if master_par.image_geometry == 'GROUND_RANGE':
        process(['gc_map_grd', master + '.par'] + gc_map_args, logpath=path_log)
    else:
        process(['gc_map', master + '.par', '-'] + gc_map_args, logpath=path_log)

    for item in ['dem_seg', 'sim_map', 'u', 'v', 'psi', 'pix', 'inc']:
        if n.isappreciated(item):
            envi.hdr(n.dem_seg + '.par', n.get(item) + '.hdr')

    sim_width = ISPPar(n.dem_seg + '.par').width

    ######################################################################################
    # Estimate pixel scattering area based on DEM and ellipsoid
    process(['pixel_area', master + '.par', n.dem_seg + '.par', n.dem_seg, n.lut_coarse,
             n.ls_map, n.inc, n.pix_dem], logpath=path_log)
    process(['radcal_MLI', master, master + '.par', '-', master + '_cal', '-', 0, 0, 1,
             0.0, '-', n.pix_ell], logpath=path_log)
    os.remove(master + '_cal')

    if sarsimulation is True:
        n.appreciate(['lut_fine', 'diff.par', 'offs', 'ccp'])
        n.depreciate(['offsets', 'coffs', 'coffsets'])
        ##################################################################################
        # cross correlation approach 1 ###################################################
        # ################################################################################
        # dempar = ISPPar(dem_seg + '.par')
        # samples_dem = dempar.width
        # samples_mli = master_par.range_samples
        # lines_mli = master_par.azimuth_lines
        # gamma(['geocode', lut, sim_map, samples_dem, sim_sar, samples_mli, lines_mli, func_interp], logpath=path_log)
        # init_offset(master, sim_sar, diffpar, path_log)
        # correlate(master, sim_sar, diffpar, offs, ccp, coffs=coffs, coffsets=coffsets, offsets=offs+'.txt', path_log=path_log, maxwin=256)
        ##################################################################################
        # cross correlation approach 2 ###################################################
        ##################################################################################
        init_offset(master, n.pix_dem, n.diff_par, path_log)
        # correlate(master, pixel_area_coarse, diffpar, offs, ccp, coffs=coffs, coffsets=coffsets, offsets=offsets, path_log=path_log, maxwin=256)
        # correlate2(master, pixel_area_coarse, diffpar, offs, ccp, offsets=offsets, coffs=coffs, coffsets=coffsets, path_log=path_log)
        # correlate3(master, pixel_area_coarse, diffpar, offs, ccp, offsets=offsets, coffs=coffs, coffsets=coffsets)
        print '#####################################'
        print 'cross correlation step 1'
        correlate5(n.pix_dem, master, n.diff_par, n.offs, n.ccp, n.offsets, n.coffs,
                   n.coffsets)
        print '#####################################'
        print 'cross correlation step 2'
        correlate5(n.pix_dem, master, n.diff_par, n.offs, n.ccp, n.offsets, n.coffs,
                   n.coffsets)
        print '#####################################'
        # print 'cross correlation step 3'
        # correlate4(pix_dem, master, diffpar, offs, ccp, offsets, coffs, coffsets)
        # print '#####################################'
        ##################################################################################

        # lookup table and pixel area refinement
        process(['gc_map_fine', n.lut_coarse, sim_width, n.diff_par, n.lut_fine, 0],
                logpath=path_log)
        process(['pixel_area', master + '.par', n.dem_seg + '.par', n.dem_seg,
                 n.lut_fine, n.ls_map, n.inc, n.pix_dem], logpath=path_log)

        lut_final = n.lut_fine
    else:
        lut_final = n.lut_coarse

    ######################################################################################
    # normalization and backward geocoding approach 1 ####################################
    ######################################################################################
    # for image in images:
    #     process(['geocode_back', image, master_par.range_samples, lut_final, image + '_geo', sim_width, '-', func_geoback], logpath=path_log)
    #     process(['product', image + '_geo', n.pix, image + '_geo_pan', sim_width, 1, 1, 0], logpath=path_log)
    #     process(['lin_comb', 1, image + '_geo_pan', 0, math.cos(math.radians(master_par.incidence_angle)), image + '_geo_pan_flat', sim_width], logpath=path_log)
    #     process(['sigma2gamma', image + '_geo_pan_flat', n.inc, image + '_geo_norm', sim_width], logpath=path_log)
    #     envi.hdr(n.dem_seg + '.par', image + '_geo_norm.hdr')
    ######################################################################################
    # normalization and backward geocoding approach 2 ####################################
    ######################################################################################
    n.appreciate(['pix_ratio', 'offs', 'ccp'])

    # Calculate pixel area normalization factor (ie., the ratio of the pixel area assuming
    # flat terrain to the pixel area estimated based on DEM)
    process(['ratio', n.pix_ell, n.pix_dem, n.pix_ratio, master_par.range_samples, 1, 1],
            logpath=path_log)

    for image in images:
        # pixel area compensation (i.e. topographic normalization)
        process(
            ['product', image, n.pix_ratio, image + '_pan', master_par.range_samples, 1,
             1], logpath=path_log)
        # conversion from range-doppler to map coordinates (i.e. geocoding)
        process(['geocode_back', image + '_pan', master_par.range_samples, lut_final,
                 image + '_pan_geo', sim_width, 0, func_geoback], logpath=path_log)
        #
        process(['lin_comb', 1, image + '_pan_geo', 0,
                 math.cos(math.radians(master_par.incidence_angle)),
                 image + '_pan_geo_flat', sim_width], logpath=path_log)
        # correction of incidence angle influence on backscatter
        process(['sigma2gamma', image + '_pan_geo_flat', n.inc, image + '_geo_norm',
                 sim_width], logpath=path_log)
        envi.hdr(n.dem_seg + '.par', image + '_geo_norm.hdr')
    ######################################################################################
    # conversion to (dB and) geotiff
    for image in images:
        for scale in scaling:
            if scale == 'db':
                process(['linear_to_dB', image + '_geo_norm', image + '_geo_norm_db',
                         sim_width, 0, -99], logpath=path_log)
                envi.hdr(n.dem_seg + '.par', image + '_geo_norm_db.hdr')
            suffix = {'linear': '', 'db': '_db'}[scale]
            infile = image + '_geo_norm{}'.format(suffix)
            outname = '{}_geo_norm{}.tif'.format(os.path.basename(image), suffix)
            outfile = os.path.join(outdir, outname)
            process(['data2geotiff', n.dem_seg + '.par', infile, 2, outfile],
                    logpath=path_log)

    if scene.sensor in ['S1A', 'S1B']:
        shutil.copyfile(os.path.join(scene.scene, 'manifest.safe'),
                        os.path.join(outdir, scene.outname_base() + '_manifest.safe'))
    if cleanup:
        shutil.rmtree(scene.scene)


def init_offset(master, slave, off, path_log, thres=0.2):
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
        process(['create_diff_par', master + '.par', '-', off, 1, 0], logpath=path_log)

    par = ISPPar(master + '.par')
    mode = 'SLC' if par.image_format in ['FCOMPLEX', 'SCOMPLEX'] else 'MLI'

    # first offset estimation using orbit state-vectors and image parameters (most important in case of very large offsets)
    # this is only to be performed while coregistering two scenes not while geocoding
    # todo: query input file type in order to decide whether command is executed
    cmd = 'init_offset_orbit' if mode == 'SLC' else 'init_offset_orbitm'
    if mode == 'SLC':
        process([cmd, master + '.par', slave + '.par', off], path_out, path_log)

    mlk = Spacing(par)
    passed = False
    patchsizes = [128, 256, 512, 1024]
    patchstack = Stack(patchsizes[0])

    while not passed and not patchstack.empty():
        win = patchstack.pop()
        for factor in [3, 2, 1]:
            try:
                args = [int(mlk.rlks) * factor, int(mlk.azlks) * factor, '-', '-', '-',
                        '-', thres, win, win]
                if mode == 'SLC':
                    process(
                        ['init_offset', master, slave, master + '.par', slave + '.par',
                         off] + args, path_out, path_log)
                else:
                    process(['init_offsetm', master, slave, off] + args, path_out,
                            path_log)
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
                patchstack.push([win, patchsizes[patchsizes.index(win) + 1]])

    if not passed:
        raise RuntimeError(
            'no initial offset found; consider verifying scene overlap or choice of polarization')


def ovs(parfile, targetres):
    """
    compute DEM oversampling factors for a target resolution in meters
    """
    # read DEM parameter file
    dempar = ISPPar(parfile)

    # extract coordinates and pixel posting of the DEM
    if hasattr(dempar, 'post_north'):
        post_north, post_east = [abs(float(x)) for x in
                                 [dempar.post_north, dempar.post_east]]
    else:
        res_lat, res_lon = [abs(float(x)) for x in [dempar.post_lat, dempar.post_lon]]

        # compute center coordinate
        lat = float(dempar.corner_lat) - (res_lat * (dempar.nlines // 2))
        lon = float(dempar.corner_lon) + (res_lon * (dempar.width // 2))

        # convert DEM resolution to meters
        post_north = haversine(lat, lon, lat + res_lat, lon)
        post_east = haversine(lat, lon, lat, lon + res_lon)

    # compute resampling factors for the DEM
    ovs_lat = post_north / targetres
    ovs_lon = post_east / targetres
    return ovs_lat, ovs_lon


def multilook(infile, outfile, targetres):
    """
    multilooking of SLC and MLI images

    targetres: the target resolution in ground range

    if the image is in slant range the ground range resolution is computed by dividing the range pixel spacing by
    the sine of the incidence angle

    the looks in range and azimuth are chosen to approximate the target resolution by rounding the ratio between
    target resolution and ground range/azimuth pixel spacing to the nearest integer
    """
    # read the input parameter file
    par = ISPPar(infile + '.par')

    # compute the range looks
    if par.image_geometry == 'SLANT_RANGE':
        # compute the ground range resolution
        groundRangePS = par.range_pixel_spacing / (math.sin(math.radians(par.incidence_angle)))
        rlks = int(round(float(targetres) / groundRangePS))
    else:
        rlks = int(round(float(targetres) / par.range_pixel_spacing))
    # compute the azimuth looks
    azlks = int(round(float(targetres) / par.azimuth_pixel_spacing))

    # set the look factors to 1 if they were computed to be 0
    rlks = rlks if rlks > 0 else 1
    azlks = azlks if azlks > 0 else 1

    if par.image_format in ['SCOMPLEX', 'FCOMPLEX']:
        # multilooking for SLC images
        process(['multi_look', infile, infile + '.par', outfile, outfile + '.par', rlks, azlks])
    else:
        # multilooking for MLI images
        process(['multi_look_MLI', infile, infile + '.par', outfile, outfile + '.par', rlks, azlks])
    envi.hdr(outfile + '.par')
