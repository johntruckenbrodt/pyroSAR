# !/usr/bin/env python
# ******************************************************************************
#  Name:     iMad.py
#  Purpose:  Perfrom IR-MAD change detection on bitemporal, multispectral
#            imagery
#  Usage:
#    python iMad.py
#
#  Copyright (c) 2013, Mort Canty
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

import canty_auxil.auxil as auxil
import numpy as np
from scipy import linalg, stats
from osgeo import gdal
from osgeo.gdalconst import GA_ReadOnly, GDT_Float32
import sets, os, sys, time


def main():
    gdal.AllRegister()
    path = auxil.select_directory('Choose working directory')
    if path:
        os.chdir(path)
    #  first image
    file1 = auxil.select_infile(title='Choose first image')
    if file1:
        inDataset1 = gdal.Open(file1, GA_ReadOnly)
        cols = inDataset1.RasterXSize
        rows = inDataset1.RasterYSize
        bands = inDataset1.RasterCount
    else:
        return
    pos1 = auxil.select_pos(bands)
    if not pos1:
        return
    dims = auxil.select_dims([0, 0, cols, rows])
    if dims:
        x10, y10, cols1, rows1 = dims
    else:
        return
    #  second image
    file2 = auxil.select_infile(title='Choose second image')
    if file2:
        inDataset2 = gdal.Open(file2, GA_ReadOnly)
        cols = inDataset2.RasterXSize
        rows = inDataset2.RasterYSize
        bands = inDataset2.RasterCount
    else:
        return
    pos2 = auxil.select_pos(bands)
    if not pos2:
        return
    dims = auxil.select_dims([0, 0, cols, rows])
    if dims:
        x20, y20, cols, rows = dims
    else:
        return
    #  penalization
    lam = auxil.select_penal(0.0)
    if lam is None:
        return
    #  outfile
    outfile, fmt = auxil.select_outfilefmt()
    if not outfile:
        return
    #  match dimensions
    bands = len(pos2)
    if (rows1 != rows) or (cols1 != cols) or (len(pos1) != bands):
        sys.stderr.write("Size mismatch")
        sys.exit(1)
    print '========================='
    print '       iMAD'
    print '========================='
    print time.asctime()
    print 'time1: ' + file1
    print 'time2: ' + file2
    print 'Delta    [canonical correlations]'
    #  iteration of MAD
    cpm = auxil.Cpm(2 * bands)
    delta = 1.0
    oldrho = np.zeros(bands)
    itr = 0
    tile = np.zeros((cols, 2 * bands))
    sigMADs = 0
    means1 = 0
    means2 = 0
    A = 0
    B = 0
    rasterBands1 = []
    rasterBands2 = []
    for b in pos1:
        rasterBands1.append(inDataset1.GetRasterBand(b))
    for b in pos2:
        rasterBands2.append(inDataset2.GetRasterBand(b))
    while (delta > 0.001) and (itr < 100):
        #      spectral tiling for statistics
        for row in range(rows):
            for k in range(bands):
                tile[:, k] = rasterBands1[k].ReadAsArray(x10, y10 + row, cols, 1)
                tile[:, bands + k] = rasterBands2[k].ReadAsArray(x20, y20 + row, cols, 1)
            #          eliminate no-data pixels (assuming all zeroes)
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
            #     weighted covariance matrices and means
        S = cpm.covariance()
        means = cpm.means()
        #     reset prov means object
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
        #     solution of generalized eigenproblems
        if bands > 1:
            mu2a, A = auxil.geneiv(c1, b1)
            mu2b, B = auxil.geneiv(c2, b2)
            #          sort a
            idx = np.argsort(mu2a)
            A = A[:, idx]
            #          sort b
            idx = np.argsort(mu2b)
            B = B[:, idx]
            mu2 = mu2b[idx]
        else:
            mu2 = c1 / b1
            A = 1 / np.sqrt(b1)
            B = 1 / np.sqrt(b2)
        #      canonical correlations
        mu = np.sqrt(mu2)
        a2 = np.diag(A.T * A)
        b2 = np.diag(B.T * B)
        sigma = np.sqrt((2 - lam * (a2 + b2)) / (1 - lam) - 2 * mu)
        rho = mu * (1 - lam) / np.sqrt((1 - lam * a2) * (1 - lam * b2))
        #      stopping criterion
        delta = max(abs(rho - oldrho))
        print delta, rho
        oldrho = rho
        #      tile the sigmas and means
        sigMADs = np.tile(sigma, (cols, 1))
        means1 = np.tile(means[0:bands], (cols, 1))
        means2 = np.tile(means[bands::], (cols, 1))
        #      ensure sum of positive correlations between X and U is positive
        D = np.diag(1 / np.sqrt(np.diag(s11)))
        s = np.ravel(np.sum(D * s11 * A, axis=0))
        A = A * np.diag(s / np.abs(s))
        #      ensure positive correlation between each pair of canonical variates
        cov = np.diag(A.T * s12 * B)
        B = B * np.diag(cov / np.abs(cov))
        itr += 1
    # write results to disk
    driver = gdal.GetDriverByName(fmt)
    outDataset = driver.Create(outfile, cols, rows, bands + 1, GDT_Float32)
    projection = inDataset1.GetProjection()
    geotransform = inDataset1.GetGeoTransform()
    if geotransform is not None:
        gt = list(geotransform)
        gt[0] = gt[0] + x10 * gt[1]
        gt[3] = gt[3] + y10 * gt[5]
        outDataset.SetGeoTransform(tuple(gt))
    if projection is not None:
        outDataset.SetProjection(projection)
    outBands = []
    for k in range(bands + 1):
        outBands.append(outDataset.GetRasterBand(k + 1))
    for row in range(rows):
        for k in range(bands):
            tile[:, k] = rasterBands1[k].ReadAsArray(x10, y10 + row, cols, 1)
            tile[:, bands + k] = rasterBands2[k].ReadAsArray(x20, y20 + row, cols, 1)
        mads = np.asarray((tile[:, 0:bands] - means1) * A - (tile[:, bands::] - means2) * B)
        chisqr = np.sum((mads / sigMADs) ** 2, axis=1)
        for k in range(bands):
            outBands[k].WriteArray(np.reshape(mads[:, k], (1, cols)), 0, row)
        outBands[bands].WriteArray(np.reshape(chisqr, (1, cols)), 0, row)
    for outBand in outBands:
        outBand.FlushCache()
    outDataset = None
    inDataset1 = None
    inDataset2 = None
    print 'result written to: ' + outfile
    print '--------done---------------------'


if __name__ == '__main__':
    main()
