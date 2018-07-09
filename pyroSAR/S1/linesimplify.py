##############################################################
# Utilities for simplification of lines used by pyroSAR for border noise removal
# John Truckenbrodt 2017
##############################################################

from osgeo import ogr
import numpy as np
from pyroSAR.ancillary import rescale
from .polysimplify import VWSimplifier
# import matplotlib.pyplot as plt


def simplify(x, y, maxpoints=20):
    x = map(float, x)
    y = map(float, y)
    pts = np.array(zip(x, y))
    simplifier = VWSimplifier(pts)
    sqd = []
    iter_range = range(2, maxpoints + 1)
    for i in iter_range:
        VWpts = simplifier.from_number(i)
        xn, yn = zip(*VWpts)
        out = np.sum((y - np.interp(x, xn, yn)) ** 2)
        sqd.append(out)
    # sqd /= max(sqd)
    if min(sqd) == max(sqd):
        VWpts = simplifier.from_number(2)
        return VWpts
    else:
        sqd = rescale(sqd)
        # plt.plot(sqd)
        # plt.show()
        # iter = (np.array(iter_range) - 2) / (maxpoints - 2.)
        # plt.plot(iter_range, sqd, label='residual')
        # plt.plot(iter_range, iter, color='r', label='iteration')
        # plt.plot(iter_range, iter + sqd, color='g', label='residual+iteration')
        # plt.legend(loc='upper center', shadow=True)
        # plt.show()
        # npoints = np.argmin(iter + sqd) + 2
        npoints = np.argmax(np.array(sqd) < 0.01)+2
        VWpts = simplifier.from_number(npoints)
        return VWpts


def createPoly(xn, yn, xmax, ymax):
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(0, 0)
    for item in zip(xn, yn):
        item = map(int, item)
        if item != [0, 0] and item != [xmax, ymax]:
            ring.AddPoint(item[0], item[1])
    ring.AddPoint(xmax, ymax)
    ring.AddPoint(xmax, 0)
    ring.CloseRings()
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    return poly


def reduce(seq, maxpoints=20, straighten=False):
    if min(seq) == max(seq):
        return np.array(seq)
    x = map(float, range(0, len(seq)))
    # plt.plot(seq)
    VWpts = simplify(x, seq, maxpoints)
    xn, yn = map(list, zip(*VWpts))
    # plt.plot(xn, yn, linewidth=2, color='r')
    simple = np.interp(x, xn, yn)
    seq2 = np.copy(seq)
    points = []
    for xi, yi in enumerate(seq2):
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(xi, yi)
        points.append(point)
    points = np.array(points)
    while True:
        poly = createPoly(xn, yn, len(seq2), max(seq2))
        line = ogr.Geometry(ogr.wkbLineString)
        for xi, yi in zip(xn, yn):
            line.AddPoint(xi, yi)
        dists = np.array([line.Distance(point) for point in points])
        contain = np.array([point.Within(poly) for point in points])
        dists[~contain] = 0
        points = points[(dists > 0)]
        dists = dists[(dists > 0)]
        if len(dists) == 0:
            break
        candidate = points[np.argmax(dists)]
        cp = candidate.GetPoint()
        index = np.argmin(np.array(xn) < cp[0])
        xn.insert(index, cp[0])
        yn.insert(index, cp[1])
    if straighten:
        indices = [i for i in range(0, len(xn)) if (xn[i], yn[i]) in VWpts]
        for i, j in enumerate(indices):
            if i < (len(indices)-1):
                if indices[i+1] > j+1:
                    dx = abs(xn[j] - xn[indices[i + 1]])
                    dy = abs(yn[j] - yn[indices[i + 1]])
                    if dx > dy:
                        seg_y = yn[j:indices[i + 1]+1]
                        for k in range(j, indices[i + 1]+1):
                            yn[k] = min(seg_y)
    # plt.plot(xn, yn, linewidth=2, color='g')
    # plt.show()
    return np.interp(x, xn, yn).astype(int)
