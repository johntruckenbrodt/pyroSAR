import gdal
import math
import numpy as np
from gdalconst import *
from time import asctime

import raster


def bandmean(infile, outfile, maxmem=2000):
    print asctime()
    print "computing multitemporal mean"
    test = raster.Raster(infile)
    maxlines = maxmem//(test.cols*test.bands*4/1048576.)
    lines = range(0, test.rows)

    divisor = sorted([int(x % math.ceil(len(lines)/maxlines)) for x in lines])

    outarray = np.zeros((test.rows, test.cols), dtype=float)

    for x in set(divisor):
        indices = [y for y in lines if divisor[lines.index(y)] == x]
        arr = test.raster.ReadAsArray(0, indices[0], test.cols, len(indices))
        outarray[indices[0]:indices[-1]+1, :] = np.mean(arr, axis=0)

    driver = gdal.GetDriverByName("ENVI")
    out = driver.Create(outfile, test.cols, test.rows, 1, GDT_Float32)
    maskout = out.GetRasterBand(1)
    maskout.WriteArray(outarray, 0, 0)
    maskout.FlushCache()
    out.SetGeoTransform(test.raster.GetGeoTransform())
    out.SetProjection(test.raster.GetProjection())
    print asctime()
