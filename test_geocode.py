
import os
import re
import sys
import shutil
from time import asctime
import subprocess as sp
from ancillary import dissolve
from gamma.util import ISPPar
from gamma.error import gammaErrorHandler


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
        print "windows size:", winsize
        try:
            # compute the number of windows needed in range and azimuth based on the number of image pixels in both directions
            nr = int(round((float(par.range_samples)/winsize)*(1+overlap)))
            naz = str(int(int(nr)*dim_ratio))
            if mode == "SLC":
                # run([commands[0], master, slave, master+".par", slave+".par", off, offs, snr, winsize, winsize, offsets, ovs, nr, naz, thres], path_out, path_log)
                proc = Cmd([commands[0], master, slave, master+".par", slave+".par", off, offs, snr, winsize, winsize, offsets, ovs, nr, naz, thres])
            else:
                # run([commands[0], master, slave, off, offs, snr, winsize, winsize, offsets, ovs, nr, naz, thres], path_out, path_log)
                proc = Cmd([commands[0], master, slave, off, offs, snr, winsize, winsize, offsets, ovs, nr, naz, thres])
            proc.process()
            for line in proc.stdout.split("\n"):
                if line.startswith("number of offsets above SNR threshold:"):
                    offstat = [int(x) for x in re.findall("[0-9]+", line)]
                    print "offsets: {0} of {1}".format(offstat[0], offstat[1])
            passed = True
            try:
                # run([commands[1], offs, snr, off, coffs, coffsets, thres, poly], path_out, path_log)
                proc = Cmd([commands[1], offs, snr, off, coffs, coffsets, thres, poly])
                proc.process()
                for line in proc.stdout.split("\n"):
                    if line.startswith("final model fit"):
                        modelstat = [float(x) for x in re.findall("[0-9]+\.[0-9]+", line)]
                        print "model fit (range, azimuth): {0}, {1}".format(modelstat[0], modelstat[1])
                        print "--------"
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


def correlate2(master, slave, off, offs, snr, offsets="-", coffs="-", coffsets="-", path_log=None, maxwin=2048, minwin=128, overlap=.3, poly=4, ovs=2, thres=7.0):
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
    winsize = minwin
    while winsize <= maxwin:
        print "windows size:", winsize
        try:
            # compute the number of windows needed in range and azimuth based on the number of image pixels in both directions
            nr = int(round((float(par.range_samples)/winsize)*(1+overlap)))
            naz = str(int(int(nr)*dim_ratio))
            if mode == "SLC":
                # run([commands[0], master, slave, master+".par", slave+".par", off, offs, snr, winsize, winsize, offsets, ovs, nr, naz, thres], path_out, path_log)
                proc = Cmd([commands[0], master, slave, master+".par", slave+".par", off, offs, snr, winsize, winsize, offsets, ovs, nr, naz, thres])
            else:
                # run([commands[0], master, slave, off, offs, snr, winsize, winsize, offsets, ovs, nr, naz, thres], path_out, path_log)
                proc = Cmd([commands[0], master, slave, off, offs, snr, winsize, winsize, offsets, ovs, nr, naz, thres])
            proc.process()
            for line in proc.stdout.split("\n"):
                if line.startswith("number of offsets above SNR threshold:"):
                    offstat = [int(x) for x in re.findall("[0-9]+", line)]
                    print "offsets: {0} of {1}".format(offstat[0], offstat[1])
            passed = True
            try:
                # run([commands[1], offs, snr, off, coffs, coffsets, thres, poly], path_out, path_log)
                proc = Cmd([commands[1], offs, snr, off, coffs, coffsets, thres, poly])
                proc.process()
                for line in proc.stdout.split("\n"):
                    if line.startswith("final model fit"):
                        modelstat = [float(x) for x in re.findall("[0-9]+\.[0-9]+", line)]
                        print "model fit (range, azimuth): {0}, {1}".format(modelstat[0], modelstat[1])
                        print "--------"
            except RuntimeError:
                passed = False
        except ValueError:
            continue
        finally:
            winsize *= 2

    if not passed:
        for file in [offs, offsets, snr, coffs, coffsets]:
            if os.path.isfile(file):
                os.remove(file)
        raise RuntimeError("cross-correlation failed; consider verifying scene overlap or choice of polarization")


# wrapper for subprocess execution including logfile writing and command prompt piping
def run(cmd, outdir=None, logpath=None, inlist=None):
    cmd = [str(x) for x in dissolve(cmd)]
    if outdir is None:
        outdir = os.getcwd()
    if logpath is None:
        log = sp.PIPE
    else:
        index = 1 if cmd[0] in [sys.executable, "Rscript"] else 0
        logfile = os.path.join(logpath, os.path.splitext(cmd[index])[0]+".log")
        log = open(logfile, "a")

    if inlist is None:
        # proc = sp.Popen(cmd, stdin=sp.PIPE, stdout=log, stderr=err, cwd=outdir)
        proc = sp.Popen(cmd, stdin=sp.PIPE, stdout=log, stderr=sp.PIPE, cwd=outdir)
        gammaErrorHandler(proc)
    else:
        out, err = sp.Popen(cmd, stdin=sp.PIPE, stdout=log, stderr=sp.PIPE, cwd=outdir, universal_newlines=True, shell=False).communicate("".join([str(x)+"\n" for x in inlist]))
    # add line for separating log entries of repeated function calls
    if logpath is not None:
        log.write("#####################################################################\n")
        log.close()


class Cmd(object):
    def __init__(self, cmd, outdir=os.getcwd(), inlist=None):
        self.cmd = [str(x) for x in dissolve(cmd)]
        self.inlist = inlist if inlist is None else [str(x) for x in dissolve(inlist)]
        self.outdir = outdir
    def process(self):
        stdin = self.inlist if self.inlist is None else "".join([x+"\n" for x in self.inlist])
        self.stdout, self.stderr = sp.Popen(self.cmd, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, cwd=os.getcwd(), universal_newlines=True, shell=False).communicate(stdin)
    def write(self, type, filename, overwrite=False):
        if not hasattr(self, type):
            raise IOError("object has no attribute {0}".format(type))
        mode = "w" if overwrite else "a"
        with open(filename, mode) as out:
            out.write(getattr(self, type))
            out.write("#####################################################################\n")


main = "/geonfs02_vol1/ve39vem/S1/test_in/S1A_IW_GRDH_1SDV_20141220T155633_20141220T155658_003805_0048BB_CE9B.SAFE/S1A______20141220T155633_"
master = main+"VH_mli3"
pixel_area = main+"VH_pixel_area"
off = master+"_diff.par"
offs = main+"VH_offs"
offsets = offs+".txt"
snr = main+"VH_snr"
coffs = main+"VH_coffs"
coffsets = main+"VH_coffsets"
path_log = "/geonfs02_vol1/ve39vem/S1/test_out/testlog"

products = [off, offs, offsets, snr, coffsets, coffs, coffsets]

# print asctime()
# for item in products:
#     if os.path.isfile(item):
#         os.remove(item)
# correlate(master, pixel_area, master+"_diff.par", offs, snr, coffs=coffs, coffsets=coffsets, offsets=offsets, path_log=path_log, minwin=256, maxwin=256, overlap=0)
# print "-----------"
# print asctime()
#
# for item in products:
#     shutil.copy(item, item+"1")


for size in [16, 32, 64, 128, 256]:
    for item in products:
        if os.path.isfile(item):
            os.remove(item)
    # for item in products:
    #     shutil.copy(item+"1", item)
    correlate(master, pixel_area, master+"_diff.par", offs, snr, coffs=coffs, coffsets=coffsets, offsets=offsets, path_log=path_log, minwin=size, maxwin=size, overlap=0.5, poly=3)
