"""
Visvalingam-Whyatt method of poly-line vertex reduction
Visvalingam, M and Whyatt J D (1993)
"Line Generalisation by Repeated Elimination of Points", Cartographic J., 30 (1), 46 - 51
Described here:
https://web.archive.org/web/20100428020453/http://www2.dcs.hull.ac.uk/CISRG/publications/DPs/DP10/DP10.html
=========================================
The MIT License (MIT)
Copyright (c) 2014 Elliot Hallmark
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
================================

code was obtained from https://github.com/Permafacture/Py-Visvalingam-Whyatt/blob/master/polysimplify.py
minor edits for Python3 compatibility by John Truckenbrodt 2019
"""

from numpy import array, argmin
import numpy as np


def triangle_area(p1, p2, p3):
    """
    calculates the area of a triangle given its vertices
    """
    return abs(p1[0] * (p2[1] - p3[1]) + p2[0] * (p3[1] - p1[1]) + p3[0] * (p1[1] - p2[1])) / 2.


def triangle_areas_from_array(arr):
    """
    take an (N,2) array of points and return an (N,1)
    array of the areas of those triangles, where the first
    and last areas are np.inf
    see triangle_area for algorithm
    """

    result = np.empty((len(arr),), arr.dtype)
    result[0] = np.inf
    result[-1] = np.inf

    p1 = arr[:-2]
    p2 = arr[1:-1]
    p3 = arr[2:]

    # an accumulators to avoid unnecessary intermediate arrays
    accr = result[1:-1]  # Accumulate directly into result
    acc1 = np.empty_like(accr)

    np.subtract(p2[:, 1], p3[:, 1], out=accr)
    np.multiply(p1[:, 0], accr, out=accr)
    np.subtract(p3[:, 1], p1[:, 1], out=acc1)
    np.multiply(p2[:, 0], acc1, out=acc1)
    np.add(acc1, accr, out=accr)
    np.subtract(p1[:, 1], p2[:, 1], out=acc1)
    np.multiply(p3[:, 0], acc1, out=acc1)
    np.add(acc1, accr, out=accr)
    np.abs(accr, out=accr)
    accr /= 2.
    # Notice: accr was writing into result, so the answer is in there
    return result


# the final value in thresholds is np.inf, which will never be
# the min value.  So, I am safe in "deleting" an index by
# just shifting the array over on top of it
def remove(s, i):
    """
    Quick trick to remove an item from a numpy array without
    creating a new object.  Rather than the array shape changing,
    the final value just gets repeated to fill the space.
    ~3.5x faster than numpy.delete
    """
    s[i:-1] = s[i + 1:]


class VWSimplifier(object):
    def __init__(self, pts):
        """
        Initialize with points. takes some time to build
        the thresholds but then all threshold filtering later
        is ultra fast
        """
        self.pts = np.array(pts)
        self.thresholds = self.build_thresholds()
        self.ordered_thresholds = sorted(self.thresholds, reverse=True)

    def build_thresholds(self):
        """
        compute the area value of each vertex, which one would
        use to mask an array of points for any threshold value.
        returns a numpy.array (length of pts)  of the areas.
        """
        pts = self.pts
        nmax = len(pts)
        real_areas = triangle_areas_from_array(pts)
        real_indices = list(range(nmax))

        # destructable copies
        # ARG! areas=real_areas[:] doesn't make a copy!
        areas = np.copy(real_areas)
        i = real_indices[:]

        # pick first point and set up for loop
        min_vert = int(argmin(areas))
        this_area = areas[min_vert]
        #  areas and i are modified for each point finished
        remove(areas, min_vert)  # faster
        # areas = np.delete(areas,min_vert) #slower
        real_idx = i.pop(min_vert)

        # cntr = 3
        while this_area < np.inf:
            '''min_vert was removed from areas and i.  Now,
            adjust the adjacent areas and remove the new
            min_vert.
            Now that min_vert was filtered out, min_vert points
            to the point after the deleted point.'''

            skip = False  # modified area may be the next minvert

            try:
                right_area = triangle_area(pts[i[min_vert - 1]],
                                           pts[i[min_vert]], pts[i[min_vert + 1]])
            except IndexError:
                # trying to update area of endpoint. Don't do it
                pass
            else:
                right_idx = i[min_vert]
                if right_area <= this_area:
                    # even if the point now has a smaller area,
                    # it ultimately is not more significant than
                    # the last point, which needs to be removed
                    # first to justify removing this point.
                    # Though this point is the next most significant
                    right_area = this_area

                    # min_vert refers to the point to the right of
                    # the previous min_vert, so we can leave it
                    # unchanged if it is still the min_vert
                    skip = min_vert

                # update both collections of areas
                real_areas[right_idx] = right_area
                areas[min_vert] = right_area

            if min_vert > 1:
                # cant try/except because 0-1=-1 is a valid index
                left_area = triangle_area(pts[i[min_vert - 2]],
                                          pts[i[min_vert - 1]], pts[i[min_vert]])
                if left_area <= this_area:
                    # same justification as above
                    left_area = this_area
                    skip = min_vert - 1
                real_areas[i[min_vert - 1]] = left_area
                areas[min_vert - 1] = left_area

            # only argmin if we have too.
            min_vert = skip or argmin(areas)
            real_idx = i.pop(min_vert)
            this_area = areas[min_vert]
            # areas = np.delete(areas,min_vert) #slower
            remove(areas, min_vert)  # faster
            '''if sum(np.where(areas==np.inf)[0]) != sum(list(reversed(range(len(areas))))[:cntr]):
              print "broke:",np.where(areas==np.inf)[0],cntr
              break
            cntr+=1
            #if real_areas[0]<np.inf or real_areas[-1]<np.inf:
            #  print "NO!", real_areas[0], real_areas[-1]
            '''
        return real_areas

    def from_threshold(self, threshold):
        return self.pts[self.thresholds >= threshold]

    def from_number(self, n):
        thresholds = self.ordered_thresholds
        try:
            threshold = thresholds[int(n)]
        except IndexError:
            return self.pts
        return self.pts[self.thresholds > threshold]

    def from_ratio(self, r):
        if r <= 0 or r > 1:
            raise ValueError("Ratio must be 0<r<=1")
        else:
            return self.from_number(r * len(self.thresholds))


class WKTSimplifier(VWSimplifier):
    """
    VWSimplifier that returns strings suitable for WKT creation
    """

    def __init__(self, *args, **kwargs):
        if 'precision' in kwargs:
            p = kwargs.pop('precision')
        else:
            p = None
        VWSimplifier.__init__(self, *args, **kwargs)
        self.set_precision(p)

    def set_precision(self, precision):
        if precision:
            self.pts_as_strs = self.pts.astype('S%s' % precision)
        else:
            self.pts_as_strs = self.pts.astype(str)

    '''slow
    def from_threshold(self,threshold,precision=None):
        arr = np.array2string(self.pts[self.thresholds > threshold],precision=precision)
        return arr.replace('[[ ','(').replace(']]',')').replace(']\n [ ',',')
    '''

    def wkt_from_threshold(self, threshold, precision=None):
        if precision:
            self.set_precision(precision)
        pts = self.pts_as_strs[self.thresholds >= threshold]
        return '(%s)' % ','.join(['%s %s' % (x, y) for x, y in pts])

    def wkt_from_number(self, n, precision=None):
        thresholds = self.ordered_thresholds
        if n < 3: n = 3  # For polygons. TODO something better
        try:
            threshold = thresholds[int(n)]
        except IndexError:
            threshold = 0

        return self.wkt_from_threshold(threshold, precision=precision)

    def wkt_from_ratio(self, r, precision=None):
        if r <= 0 or r > 1:
            raise ValueError("Ratio must be 0<r<=1")
        else:
            return self.wkt_from_number(r * len(self.thresholds))


try:
    from django.contrib.gis.gdal import OGRGeometry, OGRException
    from django.contrib.gis.geos import GEOSGeometry, fromstr
except ImportError:
    class GDALSimplifier(object):
        """
        Dummy object that would be replaced by a real one if correct module exists
        """

        def __init__(*args, **kwargs):
            print("""
                  django.contrib.gis.gdal not found.
                  GDALSimplifier not available.
                  """)
else:
    from json import loads
    import re

    p = re.compile('([ 0123456789.]+) ([0123456789.]+)')


    class GDALSimplifier(object):
        """
        Warning, there is a slight loss of precision just in the
        conversion from geometry object to numpy.array even if no
        threshold is applied.  ie:

        originalpolygeom.area   ->   413962.65495176613
        gdalsimplifierpoly.area ->   413962.65495339036
        """

        def __init__(self, geom, precision=None, return_GDAL=True):
            """
            accepts a gdal.OGRGeometry or geos.GEOSGeometry
            object and wraps multiple
            VWSimplifiers.  set return_GDAL to False for faster
            filtering with arrays of floats returned instead of
            geometry objects.
            """
            global p
            self.return_GDAL = return_GDAL
            if isinstance(geom, OGRGeometry):
                name = geom.geom_name
                self.Geometry = lambda w: OGRGeometry(w, srs=geom.srs)
                self.pts = np.array(geom.tuple)
            elif isinstance(geom, GEOSGeometry):
                name = geom.geom_type.upper()
                self.Geometry = lambda w: fromstr(w)
                self.pts = np.array(geom.tuple)
            elif isinstance(geom, unicode) or isinstance(geom, str):
                # assume wkt
                # for WKT
                def str2tuple(q):
                    return '(%s,%s)' % (q.group(1), q.group(2))

                self.return_GDAL = False  # don't even try
                self.Geometry = lambda w: w  # this will never be used
                name, pts = geom.split(' ', 1)
                self.pts = loads(p.sub(str2tuple, pts). \
                                 replace('(', '[').replace(')', ']'))
            self.precision = precision
            if name == 'LINESTRING':
                self.maskfunc = self.linemask
                self.buildfunc = self.linebuild
                self.fromnumfunc = self.notimplemented
            elif name == "POLYGON":
                self.maskfunc = self.polymask
                self.buildfunc = self.polybuild
                self.fromnumfunc = self.notimplemented
            elif name == "MULTIPOLYGON":
                self.maskfunc = self.multimask
                self.buildfunc = self.multibuild
                self.fromnumfunc = self.notimplemented
            else:
                raise RuntimeError("""
             Only types LINESTRING, POLYGON and MULTIPOLYGON
             supported, but got %s""" % name)
            # sets self.simplifiers to a list of VWSimplifiers
            self.buildfunc()

        # rather than concise, I'd rather be explicit and clear.

        def pt2str(self, pt):
            """make length 2 numpy.array.__str__() fit for wkt"""
            return ' '.join(pt)

        def linebuild(self):
            self.simplifiers = [WKTSimplifier(self.pts)]

        def line2wkt(self, pts):
            return u'LINESTRING %s' % pts

        def linemask(self, threshold):
            get_pts = self.get_pts
            pts = get_pts(self.simplifiers[0], threshold)
            if self.return_GDAL:
                return self.Geometry(self.line2wkt(pts))
            else:
                return pts

        def polybuild(self):
            list_of_pts = self.pts
            result = []
            for pts in list_of_pts:
                result.append(WKTSimplifier(pts))
            self.simplifiers = result

        def poly2wkt(self, list_of_pts):
            return u'POLYGON (%s)' % ','.join(list_of_pts)

        def polymask(self, threshold):
            get_pts = self.get_pts
            sims = self.simplifiers
            list_of_pts = [get_pts(sim, threshold) for sim in sims]
            if self.return_GDAL:
                return self.Geometry(self.poly2wkt(list_of_pts))
            else:
                return array(list_of_pts)

        def multibuild(self):
            list_of_list_of_pts = self.pts
            result = []
            for list_of_pts in list_of_list_of_pts:
                subresult = []
                for pts in list_of_pts:
                    subresult.append(WKTSimplifier(pts))
                result.append(subresult)
            self.simplifiers = result

        def multi2wkt(self, list_of_list_of_pts):
            outerlist = []
            for list_of_pts in list_of_list_of_pts:
                outerlist.append('(%s)' % ','.join(list_of_pts))
            return u'MULTIPOLYGON (%s)' % ','.join(outerlist)

        def multimask(self, threshold):
            loflofsims = self.simplifiers
            result = []
            get_pts = self.get_pts
            if self.return_GDAL:
                ret_func = lambda r: self.Geometry(self.multi2wkt(r))
            else:
                ret_func = lambda r: r
            for list_of_simplifiers in loflofsims:
                subresult = []
                for simplifier in list_of_simplifiers:
                    subresult.append(get_pts(simplifier, threshold))
                result.append(subresult)
            return ret_func(result)

        def notimplemented(self, n):
            print('This function is not yet implemented')

        def from_threshold(self, threshold):
            precision = self.precision
            if self.return_GDAL:
                self.get_pts = lambda obj, t: obj.wkt_from_threshold(t, precision)
            else:
                self.get_pts = lambda obj, t: obj.from_threshold(t)
            return self.maskfunc(threshold)

        def from_number(self, n):
            precision = self.precision
            if self.return_GDAL:
                self.get_pts = lambda obj, t: obj.wkt_from_number(t, precision)
            else:
                self.get_pts = lambda obj, t: obj.from_number(t)
            return self.maskfunc(n)

        def from_ratio(self, r):
            precision = self.precision
            if self.return_GDAL:
                self.get_pts = lambda obj, t: obj.wkt_from_ratio(t, precision)
            else:
                self.get_pts = lambda obj, t: obj.from_ratio(t)
            return self.maskfunc(r)


def fancy_parametric(k):
    """
    good k's: .33,.5,.65,.7,1.3,1.4,1.9,3,4,5
    """
    cos = np.cos
    sin = np.sin
    xt = lambda t: (k - 1) * cos(t) + cos(t * (k - 1))
    yt = lambda t: (k - 1) * sin(t) - sin(t * (k - 1))
    return xt, yt


if __name__ == "__main__":
    from time import time

    n = 5000
    thetas = np.linspace(0, 16 * np.pi, n)
    xt, yt = fancy_parametric(1.4)
    pts = np.array([[xt(t), yt(t)] for t in thetas])
    start = time()
    simplifier = VWSimplifier(pts)
    pts = simplifier.from_number(1000)
    end = time()
    print("%s vertices removed in %02f seconds" % (n - len(pts), end - start))

    import matplotlib

    matplotlib.use('AGG')
    import matplotlib.pyplot as plot

    plot.plot(pts[:, 0], pts[:, 1], color='r')
    plot.savefig('visvalingam.png')
    print("saved visvalingam.png")
    # plot.show()
