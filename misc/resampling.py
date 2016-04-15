import gdal
import raster
from math import floor
import numpy as np


def seq(start, stop, step=1):
    n = int(round((stop - start)/float(step)))
    if n > 1:
        return [start + step*i for i in range(n+1)]
    else:
        return []


def bilinear(robj, spacing):

    # raster dimensions
    nx = robj.cols
    ny = robj.rows
    xres, yres, xmin, xmax, ymin, ymax = [robj.geotransform[x] for x in ["xres", "yres", "xmin", "xmax", "ymin", "ymax"]]

    # half raster cell widths
    hx = xres/2.0
    hy = yres/2.0

    # compute coordinates of the new grid
    xlist = seq(xmin + xmin % spacing, xmax - xmax % spacing, spacing)
    ylist = seq(ymin + ymin % spacing, ymax - ymax % spacing, spacing)
    zlist = []

    for x in xlist:
        for y in ylist:
            # calculate raster lower bound indices from point
            fx = (x - (xmin + hx)) / xres
            fy = (y - (ymax + hy)) / yres
            ix1 = int(floor(fx))
            iy1 = int(floor(fy))

            # special case where point is on upper bounds
            if fx == float(nx - 1):
                ix1 -= 1
            if fy == float(ny - 1):
                iy1 -= 1

            # upper bound indices on raster
            ix2 = ix1 + 1
            iy2 = iy1 + 1

            # test array bounds to ensure point is within raster midpoints
            if (ix1 < 0) or (iy1 < 0) or (ix2 > nx - 1) or (iy2 > ny - 1):
                out = robj.nodata
            else:
                # calculate differences from point to bounding raster midpoints
                dx1 = x - (xmin + ix1 * xres + hx)
                dy1 = y - (ymax + iy1 * yres + hy)
                dx2 = (xmin + ix2 * xres + hx) - x
                dy2 = (ymax + iy2 * yres + hy) - y

                # use the differences to weigh the four raster values
                div = xres * yres
                out = robj.matrix(1)[iy1, ix1] * dx2 * dy2 / div \
                      + robj.matrix(1)[iy1, ix2] * dx1 * dy2 / div \
                      + robj.matrix(1)[iy2, ix1] * dx2 * dy1 / div \
                      + robj.matrix(1)[iy2, ix2] * dx1 * dy1 / div

test = raster.Raster("/pvdata2/john/RADAR/ERS/test/workdir/S10W038.hgt")
print bilinear(test, .0016)


def bilinear(px, py, no_data=np.NAN):
    # Bilinear interpolated point at (px, py) on band_array
    # example: bilinear(2790501.920, 6338905.159)
    ny, nx = band_array.shape
    # Half raster cell widths
    hx = gt[1] / 2.0
    hy = gt[5] / 2.0
    # Calculate raster lower bound indices from point
    fx = (px - (gt[0] + hx)) / gt[1]
    fy = (py - (gt[3] + hy)) / gt[5]
    ix1 = int(np.floor(fx))
    iy1 = int(np.floor(fy))
    # Special case where point is on upper bounds
    if fx == float(nx - 1):
        ix1 -= 1
    if fy == float(ny - 1):
        iy1 -= 1
    # Upper bound indices on raster
    ix2 = ix1 + 1
    iy2 = iy1 + 1
    # Test array bounds to ensure point is within raster midpoints
    if (ix1 < 0) or (iy1 < 0) or (ix2 > nx - 1) or (iy2 > ny - 1):
        return no_data
    # Calculate differences from point to bounding raster midpoints
    dx1 = px - (gt[0] + ix1 * gt[1] + hx)
    dy1 = py - (gt[3] + iy1 * gt[5] + hy)
    dx2 = (gt[0] + ix2 * gt[1] + hx) - px
    dy2 = (gt[3] + iy2 * gt[5] + hy) - py
    # Use the differences to weigh the four raster values
    div = gt[1] * gt[5]
    return (band_array[iy1, ix1] * dx2 * dy2 / div +
            band_array[iy1, ix2] * dx1 * dy2 / div +
            band_array[iy2, ix1] * dx2 * dy1 / div +
            band_array[iy2, ix2] * dx1 * dy1 / div)


    # .astype(band_array.dtype)
