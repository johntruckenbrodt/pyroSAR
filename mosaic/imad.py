import argparse

import gdal
import numpy as np
from gdalconst import GDT_Float32
from scipy import linalg, stats

from mosaic_aux import Cpm, geneiv
from spatial import raster


def imad(file1, file2, outfile, lam=0.0):

    # load meta information
    data1 = raster.Raster(file1)
    data2 = raster.Raster(file2)

    # match dimensions
    if (data1.rows != data2.rows) or (data1.cols != data2.cols) or (data1.bands != data2.bands):
        return "Size mismatch"
    rows, cols, bands = data1.dim

    print 'delta    [canonical correlations]'
    # iteration of MAD
    cpm = Cpm(2 * bands)
    delta = 1.0
    oldrho = np.zeros(bands)
    itr = 0
    tile = np.zeros((cols, 2 * bands))
    sigMADs = 0
    means1 = 0
    means2 = 0
    A = 0
    B = 0
    rasterBands1 = data1.layers()
    rasterBands2 = data2.layers()

    while (delta > 0.001) and (itr < 100):
        # spectral tiling for statistics
        for row in range(rows):
            for k in range(bands):
                tile[:, k] = rasterBands1[k].ReadAsArray(0, row, cols, 1)
                tile[:, bands + k] = rasterBands2[k].ReadAsArray(0, row, cols, 1)
            # eliminate no-data pixels (assuming all zeroes)
            tst1 = np.sum(tile[:, 0:bands], axis=1)
            tst2 = np.sum(tile[:, bands::], axis=1)
            idx1 = set(np.where((tst1 > 0))[0])
            idx2 = set(np.where((tst2 > 0))[0])
            idx = list(idx1.intersection(idx2))
            if itr > 0:
                mads = np.asarray((tile[:, 0:bands] - means1) * A - (tile[:, bands::] - means2) * B)
                chisqr = np.sum((mads / sigMADs) ** 2, axis=1)
                wts = 1 - stats.chi2.cdf(chisqr, [bands])
                cpm.update(tile[idx, :], wts[idx])
            else:
                cpm.update(tile[idx, :])
            # weighted covariance matrices and means
        S = cpm.covariance()
        means = cpm.means()
        # reset prov means object
        cpm.__init__(2 * bands)
        s11 = S[0:bands, 0:bands]
        s11 = (1 - lam) * s11 + lam * np.eye(bands)
        s22 = S[bands:, bands:]
        s22 = (1 - lam) * s22 + lam * np.eye(bands)
        s12 = S[0:bands, bands:]
        s21 = S[bands:, 0:bands]
        c1 = s12 * linalg.inv(s22) * s21
        b1 = s11
        c2 = s21 * linalg.inv(s11) * s12
        b2 = s22
        # solution of generalized eigenproblems
        if bands > 1:
            mu2a, A = geneiv(c1, b1)
            mu2b, B = geneiv(c2, b2)
            # sort a
            idx = np.argsort(mu2a)
            A = A[:, idx]
            # sort b
            idx = np.argsort(mu2b)
            B = B[:, idx]
            mu2 = mu2b[idx]
        else:
            mu2 = c1 / b1
            A = 1 / np.sqrt(b1)
            B = 1 / np.sqrt(b2)
        # canonical correlations
        mu = np.sqrt(mu2)
        a2 = np.diag(A.T * A)
        b2 = np.diag(B.T * B)
        sigma = np.sqrt((2 - lam * (a2 + b2)) / (1 - lam) - 2 * mu)
        rho = mu * (1 - lam) / np.sqrt((1 - lam * a2) * (1 - lam * b2))
        # stopping criterion
        delta = max(abs(rho - oldrho))
        print delta, rho
        oldrho = rho
        # tile the sigmas and means
        sigMADs = np.tile(sigma, (cols, 1))
        means1 = np.tile(means[0:bands], (cols, 1))
        means2 = np.tile(means[bands::], (cols, 1))
        # ensure sum of positive correlations between X and U is positive
        D = np.diag(1 / np.sqrt(np.diag(s11)))
        s = np.ravel(np.sum(D * s11 * A, axis=0))
        A = A * np.diag(s / np.abs(s))
        # ensure positive correlation between each pair of canonical variates
        cov = np.diag(A.T * s12 * B)
        B = B * np.diag(cov / np.abs(cov))
        itr += 1
    # write results to disk
    driver = gdal.GetDriverByName("ENVI")
    outDataset = driver.Create(outfile, cols, rows, bands + 1, GDT_Float32)
    if data1.geotransform is not None:
        outDataset.SetGeoTransform(data1.raster.GetGeoTransform())
    if data1.projection is not None:
        outDataset.SetProjection(data1.projection)

    outBands = [outDataset.GetRasterBand(k + 1) for k in range(bands + 1)]
    for row in range(rows):
        for k in range(bands):
            tile[:, k] = rasterBands1[k].ReadAsArray(0, row, cols, 1)
            tile[:, bands + k] = rasterBands2[k].ReadAsArray(0, row, cols, 1)
        mads = np.asarray((tile[:, 0:bands] - means1) * A - (tile[:, bands::] - means2) * B)
        chisqr = np.sum((mads / sigMADs) ** 2, axis=1)
        for k in range(bands):
            outBands[k].WriteArray(np.reshape(mads[:, k], (1, cols)), 0, row)
        outBands[bands].WriteArray(np.reshape(chisqr, (1, cols)), 0, row)
    for outBand in outBands:
        outBand.FlushCache()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("file1", help="reference image", type=str)
    parser.add_argument("file2", help="image to be normalized", type=str)
    parser.add_argument("outfile", help="iMAD change detection image", type=str)
    parser.add_argument("-lam", "--no change threshold", type=float, default=.0, metavar="LAM", dest="lam", help=" (default: 0.0)")
    args = parser.parse_args()
    imad(args.file1, args.file2, args.outfile, args.lam)
