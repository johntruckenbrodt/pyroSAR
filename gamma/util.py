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
import re
import math
import subprocess as sp
from ancillary import run, Stack, parse_literal, haversine


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


# convert a gamma parameter file corner coordinate from EQA to UTM
class UTM(object):
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
