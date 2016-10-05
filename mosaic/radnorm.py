import argparse

import numpy as np
from scipy import stats

from mosaic_aux import orthoregress
from spatial import raster


def radnorm(reference, target, imad, fs_in, fs_out, ncp_thres=0.95):

    # load meta information
    data1 = raster.Raster(reference)
    data2 = raster.Raster(target)
    data3 = raster.Raster(imad)

    # match dimensions
    if (len({data1.rows, data2.rows, data3.rows}) != 1) or (len({data1.cols, data2.cols, data3.cols}) != 1) or (len({data1.bands, data2.bands, data3.bands-1}) != 1):
        print "Size mismatch"
        return
    rows, cols, bands = data1.dim

    chisqr = data3.matrix(data3.bands).ravel()
    ncp = 1 - stats.chi2.cdf(chisqr, [data3.bands - 1])
    idx = np.where(ncp > ncp_thres)[0]
    # split train/test in ratio 2:1
    tmp = np.asarray(range(len(idx)))
    tst = idx[np.where(np.mod(tmp, 3) == 0)]
    trn = idx[np.where(np.mod(tmp, 3) > 0)]

    # print "no-change pixels (train):", len(trn)
    # print "no-change pixels (test):", len(tst)
    aa = []
    bb = []
    rr = []
    print "band [slope, intercept, correlation]"
    for k in range(1, bands+1):
        x = data1.matrix(k).astype(float).ravel()
        y = data2.matrix(k).astype(float).ravel()
        b, a, r = orthoregress(y[trn], x[trn])
        print k, [b, a, r]
        # print 'means(tgt,ref,nrm): ', np.mean(y[tst]), np.mean(x[tst]), np.mean(a + b * y[tst])
        # print 't-test, p-value:    ', stats.ttest_rel(x[tst], a + b * y[tst])
        # print 'vars(tgt,ref,nrm)   ', np.var(y[tst]), np.var(x[tst]), np.var(a + b * y[tst])
        # print 'F-test, p-value:    ', fv_test(x[tst], a + b * y[tst])
        aa.append(a)
        bb.append(b)
        rr.append(r)
    if min([abs(x) for x in rr]) < .6:
        raise ValueError("correlation too weak")
    if fs_in is not None:
        # print "normalizing full scene"

        fsdata = raster.Raster(fs_in)
        outDataset = raster.init(fsdata, fs_out)

        for k in range(1, bands+1):
            outBand = outDataset.GetRasterBand(k)
            outBand.SetNoDataValue(fsdata.nodata)
            mat = fsdata.matrix(k).astype(float)
            mat[mat != fsdata.nodata] = aa[k-1] + bb[k-1] * mat[mat != fsdata.nodata]
            outBand.WriteArray(mat)
            outBand.FlushCache()
        del outDataset


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", help="reference image", type=str)
    parser.add_argument("target", help="image to be normalized", type=str)
    parser.add_argument("imad", help="iMAD change detection image", type=str)
    parser.add_argument("fs_in", help="full dataset to be normalized", type=str)
    parser.add_argument("fs_out", help="name of full normalized image", type=str)
    parser.add_argument("-nct", "--no change threshold", type=float, default=.95, metavar="NCPT", dest="ncp_thres", help="no change pixel threshold (default: 0.95)")
    args = parser.parse_args()
    radnorm(args.reference, args.target, args.imad, args.fs_in, args.fs_out, args.ncp_thres)
