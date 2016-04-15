#!/usr/bin/env python2.7
##############################################################
# multicore wrapper for batched geocoding of Sentinel-1 GRD data (script S1_geoGRD.py)
# John Truckenbrodt 2015
# last update 2015-02-17
##############################################################

"""
This script is intended as a wrapper for batched execution of script S1_geoGRD.py.
Instead of a single zipped S1 scene archive a folder is to be supplied, which contains a number of zipped archives.
Furthermore the number of cores can be defined, by default three cores are used.
This directory is scanned for files matching the pattern S1*GRDH*.zip (i.e. starting with "S1", having "GRDH" somewhere in the name and ending with ".zip".
Concluding, S1_geoGRD.py is called on different processing cores in parallel until all scenes are processed
"""
import S1_geoGRD
from time import asctime
from S1_auxil import init_parser
from ancillary import finder, multicore


def main(zipdir, tempdir, outdir, srtmdir, transform, logfiles, intermediates, cores, targetresolution, func_geoback, func_interp):

    files = finder(zipdir, ["S1*GRDH*.zip"])

    print "##########################################"
    print "started:", asctime()
    print "processing a total of {} scenes".format(len(files))
    multicore(S1_geoGRD.main, cores=cores, multiargs={"zipfile": files}, tempdir=tempdir, outdir=outdir, srtmdir=srtmdir,
              transform=transform, logfiles=logfiles, intermediates=intermediates, targetresolution=targetresolution, func_geoback=func_geoback, func_interp=func_interp)
    print "finished:", asctime()
    print "##########################################"


if __name__ == "__main__":
    parser = init_parser()
    parser.add_argument("-c", "--cores", default=3, help="number of cores to be used")
    args = parser.parse_args()

    main(args.zipdir, args.tempdir, args.outdir, args.srtmdir, args.transform, args.logfiles, args.intermediates, args.cores, args.targetresolution, args.func_geoback, args.func_interp)
