
import os
from ancillary import finder, multicore
from archivist import tar2zip, scan, compress
from pyroSAR import ESA

dir = "/geonfs01_vol1/ve39vem/swos"

# infiles = finder(dir, ["\.tar\.gz$"], regex=True)

# def operator(infile):
#     outname = os.path.join(os.path.dirname(infile), os.path.basename(scan(infile, "[NE][12]$")[0])+".zip")
#     tar2zip(infile, outname)
#     os.remove(infile)

# infiles = finder(dir, ["*N1"])
#
#
# def operator(infile):
#     compress(infile, infile+".zip")
#     os.remove(infile)
#
# print "processing {} scenes".format(len(infiles))
#
# multicore(operator, 20, {"infile": infiles})

files = finder(dir, ["*"])

fails = []

for i in range(len(files)):
    print format(i, "04d"), os.path.basename(files[i])
    try:
        id = ESA(files[i])
    except (AttributeError, RuntimeError):
        fails.append(files[i])
