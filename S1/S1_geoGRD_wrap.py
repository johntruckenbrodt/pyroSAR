#!/usr/bin/env python
##############################################################
# wrapper for batched geocoding of Sentinel-1 GRD data (script S1_geoGRD.py)
# John Truckenbrodt 2015
# last update 2015-02-17
##############################################################

"""
This script is intended as a wrapper for batched execution of script S1_geoGRD.py.
Instead of a single zipped S1 scene archive a folder is to be supplied, which contains a number of zipped archives.
This directory is scanned for files matching the pattern S1*GRDH*.zip (i.e. starting with "S1", having "GRDH" somewhere in the name and ending with ".zip".
Concluding, S1_geoGRD.py is called sequentially for each scene until all scenes are processed
"""

import S1_geoGRD
from time import asctime
from S1_auxil import init_parser
from ancillary import finder, blockPrint, enablePrint


def main(zipdir, tempdir, outdir, srtmdir, transform, logfiles, intermediates, targetresolution, func_geoback, func_interp):

    files = finder(zipdir, ["S1*GRDH*.zip"])

    for file in files:
        S1_geoGRD.main(zipfile=file, tempdir=tempdir, outdir=outdir, srtmdir=srtmdir, transform=transform, logfiles=logfiles, intermediates=intermediates, res_target=targetresolution, func_geoback=func_geoback, func_interp=func_interp)

if __name__ == "__main__":
    parser = init_parser()
    args = parser.parse_args()

    print "##########################################"
    print "started:", asctime()

    if args.quiet:
        blockPrint()

    main(args.zipdir, args.tempdir, args.outdir, args.srtmdir, args.transform, args.logfiles, args.intermediates, args.targetresolution, args.func_geoback, args.func_interp)

    if args.quiet:
        enablePrint()

    print "finished:", asctime()
    print "##########################################"
