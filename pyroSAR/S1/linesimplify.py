###############################################################################
# Utilities for simplification of lines used by pyroSAR for border noise removal

# Copyright (c) 2017-2020, the pyroSAR Developers.

# This file is part of the pyroSAR Project. It is subject to the
# license terms in the LICENSE.txt file found in the top-level
# directory of this distribution and at
# https://github.com/johntruckenbrodt/pyroSAR/blob/master/LICENSE.txt.
# No part of the pyroSAR project, including this file, may be
# copied, modified, propagated, or distributed except according
# to the terms contained in the LICENSE.txt file.
###############################################################################

from osgeo import ogr
import numpy as np
from spatialist.ancillary import rescale
from .polysimplify import VWSimplifier


import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
matplotlib.rcParams['font.size'] = 12


def simplify(x, y, maxpoints=20):
    x = list(map(float, x))
    y = list(map(float, y))
    pts = np.array(list(zip(x, y)))
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
        npoints = np.argmax(np.array(sqd) < 0.01) + 2
        VWpts = simplifier.from_number(npoints)
        return VWpts


def createPoly(xn, yn, xmax, ymax, plot=False):
    """
    create an OGR geometry from a sequence of indices
    
    Parameters
    ----------
    xn: numpy.ndarray
        the x indices of the points
    yn: numpy.ndarray
        the y  indices of the points
    xmax: int or float
        the maximum x index value
    ymax: int or float
        the maximum y index value

    Returns
    -------
    ogr.Geometry
    """
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint_2D(0, 0)
    for item in zip(xn, yn):
        item = list(map(int, item))
        if item != [0, 0] and item != [xmax, ymax]:
            ring.AddPoint_2D(item[0], item[1])
    ring.AddPoint_2D(xmax, ymax)
    ring.AddPoint_2D(xmax, 0)
    ring.CloseRings()
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    if plot:
        fig, ax = plt.subplots()
        pts = ring.GetPoints()
        arr = np.array(pts)
        polygon = Polygon(arr, True)
        p = PatchCollection([polygon], cmap=matplotlib.cm.jet, alpha=0.4)
        ax.add_collection(p)
        ax.autoscale_view()
        plt.scatter(arr[:, 0], arr[:, 1], s=10, color='red')
        plt.show()
    return poly


def reduce(seq, maxpoints=20, straighten=False, plot=False):
    """
    reduce the complexity of a line; the following steps are performed:
     - simplify the line using the Visvalingam-Whyatt method
     - iteratively add points on the original line back to the simplified line
       until the polygon spanned by the simplified line and (xmin, ymin) does not
       contain any further points of the original line; the polygon area is
       expected to only cover valid pixels of the image
     - optionally further straighten the result for smoother edges
    
    Parameters
    ----------
    seq: numpy.ndarray
        the 1D line sequence to be simplified
    maxpoints: int
        the maximum number points in the simplified sequence
    straighten: bool
        perform additional straightening on the simplified line?
    plot: bool
        plot the results?
    
    Returns
    -------
    numpy.ndarray
        the simplified line sequence
    """
    if min(seq) == max(seq):
        return np.array(seq)
    x = list(range(0, len(seq)))
    if plot:
        plt.plot(seq, label='ESA-corrected')
    # simplify the sequence using the Visvalingam-Whyatt algorithm
    VWpts = simplify(x, seq, maxpoints)
    xn, yn = [list(x) for x in zip(*VWpts)]
    if plot:
        plt.plot(xn, yn, linewidth=2, color='r', label='VW-simplified')
    simple = np.interp(x, xn, yn)
    # create a list of OGR points for the original border
    points = []
    for xi, yi in enumerate(seq):
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(int(xi), int(yi))
        points.append(point)
    points = np.array(points)
    while True:
        # create a polygon containing all pixels inside the simplified border
        # i.e., containing the area considered valid
        poly = createPoly(xn, yn, seq.size, int(max(seq)))
        # create an OGR line from the simplified border points
        line = ogr.Geometry(ogr.wkbLineString)
        for xi, yi in zip(xn, yn):
            line.AddPoint(xi, yi)
        # compute the distance of each original point to the simplified line
        dists = np.array([line.Distance(point) for point in points])
        # check which points are inside of the polygon
        contain = np.array([point.Within(poly) for point in points])
        # remove points outside the polygon and stop if
        # no further points outside the polygon exist
        dists[~contain] = 0
        points = points[(dists > 0)]
        dists = dists[(dists > 0)]
        if len(dists) == 0:
            break
        # select the point with the largest distance to the simplified
        # line and add it to the list of simplified points
        # this reduces the size of the polygon an thus the area considered valid
        candidate = points[np.argmax(dists)]
        cp = candidate.GetPoint()
        index = np.argmin(np.array(xn) < cp[0])
        xn.insert(index, cp[0])
        yn.insert(index, cp[1])
    if plot:
        plt.plot(xn, yn, linewidth=2, color='limegreen', label='corrected')
    
    # further straighten the line segments
    # def straight(xn, yn, VWpts):
    #     indices = [i for i in range(0, len(xn)) if (xn[i], yn[i]) in VWpts]
    #     print(indices)
    #     for i, j in enumerate(indices):
    #         if i < (len(indices) - 1):
    #             if indices[i + 1] > j + 1:
    #                 dx = abs(xn[j] - xn[indices[i + 1]])
    #                 dy = abs(yn[j] - yn[indices[i + 1]])
    #                 if dx > dy:
    #                     seg_y = yn[j:indices[i + 1] + 1]
    #                     for k in range(j, indices[i + 1] + 1):
    #                         yn[k] = min(seg_y)
    #     return yn
    
    def straight(xn, yn, VWpts):
        indices = [i for i in range(0, len(xn)) if (xn[i], yn[i]) in VWpts]
        xn_new = []
        yn_new = []
        # make all line segments horizontal or vertical
        for index in range(len(indices) - 1):
            i = indices[index]
            j = indices[index + 1]
            ymin = min(yn[i:j + 1])
            xn_new.extend([xn[i], xn[j]])
            yn_new.extend([ymin, ymin])
        # shift horizontal lines down if the preceding horizontal line has a lower y value
        # but only if the shift is less than the tolerance
        tolerance = 15
        for i in range(len(xn_new) - 2):
            if yn_new[i] == yn_new[i + 1]:
                if yn_new[i] < yn_new[i + 2] and abs(yn_new[i] - yn_new[i + 2]) < tolerance:
                    yn_new[i + 2] = yn_new[i]
                    yn_new[i + 3] = yn_new[i]
                elif (yn_new[i] > yn_new[i + 2]) \
                        and (yn_new[i + 2] == yn_new[i + 3]) \
                        and abs(yn_new[i] - yn_new[i + 2]) < tolerance:
                    yn_new[i] = yn_new[i + 2]
                    yn_new[i + 1] = yn_new[i + 2]
        return xn_new, yn_new
    
    if straighten:
        xn, yn = straight(xn, yn, VWpts)
        if plot:
            plt.plot(xn, yn, linewidth=2, color='m', label='straightened')
    if plot:
        plt.legend()
        plt.xlabel('row')
        plt.ylabel('column')
        plt.show()
    return np.interp(x, xn, yn).astype(int)
