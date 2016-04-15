#!/usr/bin/env python
# ******************************************************************************
#  Name:     S1_auxil.py
#  Purpose:  auxiliary functions for processing multispectral imagery
#  Usage:
#    import auxil
#
#  Copyright (c) 2012, Mort Canty
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

import numpy as np
import math, platform, StringIO
# comment out for GAE deployment---------
from numpy.ctypeslib import ndpointer
import Tkinter, tkFileDialog, tkSimpleDialog
from scipy.special import betainc
import ctypes
from numpy.fft import fft2, ifft2, fftshift

if platform.system() == 'Windows':
    lib = ctypes.cdll.LoadLibrary('prov_means.dll')
elif platform.system() == 'Linux':
    lib = ctypes.cdll.LoadLibrary('prov_means.so')
provmeans = lib.provmeans
provmeans.restype = None
c_double_p = ctypes.POINTER(ctypes.c_double)
provmeans.argtypes = [ndpointer(np.float64),
                      ndpointer(np.float64),
                      ctypes.c_int,
                      ctypes.c_int,
                      c_double_p,
                      ndpointer(np.float64),
                      ndpointer(np.float64)]
# ----------------------------------------

# --------------------
# contrast enhancement
# -------------------
def linstr(x):
    # linear stretch
    x = np.asarray(x, dtype=np.float32)
    mx = np.max(x)
    mn = np.min(x)
    x = (x - mn) * 255.0 / (mx - mn)
    return np.asarray(x, dtype=np.uint8)


def histeqstr(x):
    #  histogram equalization stretch of input ndarray
    hist, bin_edges = np.histogram(x, 256, (0, 256))
    cdf = hist.cumsum()
    lut = 255 * cdf / float(cdf[-1])
    return np.interp(x, bin_edges[:-1], lut)


def lin2pcstr(x):
    #  2% linear stretch
    hist, bin_edges = np.histogram(x, 256, (0, 256))
    cdf = hist.cumsum()
    lower = 0
    i = 0
    while cdf[i] < 0.02 * cdf[-1]:
        lower += 1
        i += 1
    upper = 255
    i = 255
    while cdf[i] > 0.98 * cdf[-1]:
        upper -= 1
        i -= 1
    fp = (bin_edges - lower) * 255 / (upper - lower)
    fp = np.where(bin_edges <= lower, 0, fp)
    fp = np.where(bin_edges >= upper, 255, fp)
    return np.interp(x, bin_edges, fp)


# ---------------------
# orthogonal regression
# ---------------------

def orthoregress(x, y):
    Xm = np.mean(x)
    Ym = np.mean(y)
    s = np.cov(x, y)
    R = s[0, 1] / math.sqrt(s[1, 1] * s[0, 0])
    lam, vs = np.linalg.eig(s)
    idx = np.argsort(lam)
    vs = vs[:, idx]  # increasing order, so
    b = vs[1, 1] / vs[0, 1]  # first pc is second column
    return [b, Ym - b * Xm, R]


# -------------------------
# F-test for equal variance
# -------------------------

def fv_test(x0, x1):
    # taken from IDL library
    nx0 = len(x0)
    nx1 = len(x1)
    v0 = np.var(x0)
    v1 = np.var(x1)
    if v0 > v1:
        f = v0 / v1
        df0 = nx1 - 1
        df1 = nx0 - 1
    else:
        f = v1 / v0
        df0 = nx1 - 1
        df1 = nx0 - 1
    prob = 2.0 * betainc(0.5 * df1, 0.5 * df0, df1 / (df1 + df0 * f))
    if prob > 1:
        return (f, 2.0 - prob)
    else:
        return (f, prob)

    # ------------


# Gauss filter
# ------------

def dist(n, m):
    # n cols x m rows distance array
    result = []
    for i in range(m):
        for j in range(n):
            y = float(min(i, m - i))
            x = float(min(j, n - j))
            result.append(math.sqrt(x ** 2 + y ** 2))
    return result


def gaussfilter(sigma, n, m):
    dst = dist(n, m)
    result = []
    for d in dst:
        result.append(math.exp(-d ** 2 / (2 * sigma ** 2)))
    return result


# -----------------
# provisional means
# -----------------

class Cpm(object):
    '''Provisional means algorithm'''

    def __init__(self, N):
        self.mn = np.zeros(N)
        self.cov = np.zeros((N, N))
        self.sw = 0.0000001

    def update(self, Xs, Ws=None):
        n, N = np.shape(Xs)
        if Ws is None:
            Ws = np.ones(n)
        sw = ctypes.c_double(self.sw)
        mn = self.mn
        cov = self.cov
        provmeans(Xs, Ws, N, n, ctypes.byref(sw), mn, cov)
        self.sw = sw.value
        self.mn = mn
        self.cov = cov

    def covariance(self):
        c = np.mat(self.cov / (self.sw - 1.0))
        d = np.diag(np.diag(c))
        return c + c.T - d

    def means(self):
        return self.mn

    # --------------------------


# data array (design matrix)
# --------------------------

class DataArray(object):
    '''Representation of an image or image blob (string)
    as a samples*lines by bands float32 numpy array.
    Image must have BIP or BSQ interleave and uint8 or float32 format'''

    def __init__(self, image, samples, lines, bands, interleave, dtype):
        m = samples * lines
        if isinstance(image, str):
            if dtype == 1:
                data = np.fromstring(image, dtype=np.uint8)
            else:
                data = np.fromstring(image, dtype=np.float32)
        else:
            data = np.asarray(image, np.float32)
        if interleave == 'bsq':
            data = [data[i::m] for i in range(m)]
        self.data = np.reshape(data, (m, bands))
        self.samples = samples
        self.lines = lines
        self.pixels = samples * lines
        self.bands = bands

    def covw(self, da=None, w=None):
        #  return (weighted) means of self as row vector and (weighted) variance-covariance matrix of self
        #  if da is supplied, returns (weighted) covariance matrix of self with da
        try:
            if w is None:
                w = np.ones(self.pixels, dtype=np.float32)
            if da is not None:
                b = np.transpose(da.data)
                if self.pixels != da.pixels:
                    raise
            a = np.transpose(self.data)
            sumw = np.sum(w)
            ws = np.tile(w, (self.bands, 1))
            meansa = np.mat(np.sum(np.multiply(ws, a), 1) / sumw, dtype=np.float32)
            mns = np.tile(meansa.T, (1, self.pixels))
            ac = a - mns
            ac = np.mat(np.multiply(ac, np.sqrt(ws)))
            if da is None:
                bc = ac
            else:
                meansb = np.mat(np.sum(np.multiply(ws, b), 1) / sumw, dtype=np.float32)
                mns = np.tile(meansb.T, (1, self.pixels))
                bc = b - mns
                bc = np.mat(np.multiply(bc, np.sqrt(ws)))
            covmat = ac * bc.T / sumw
            return (meansa, covmat)
        except:
            return None


def make_png_rgb(samples, lines, band1, band2, band3):
    # return an rgb png image from 3 bytestrings
    #  make a boxed row flat pixel array for png.Writer
    size = samples * lines
    RGB = np.fromstring(band1 + band2 + band3, dtype=np.uint8)
    RGB = np.reshape(RGB, (3, size)).transpose()
    RGB = np.reshape(RGB, (lines, 3 * samples)).tolist()
    #  file object
    f = StringIO.StringIO()
    #  create PNG
    w = png.Writer(samples, lines)
    w.write(f, RGB)
    return f.getvalue()


# histogram stretches for bytestring image representations
def lin(band):
    band = np.fromstring(band, dtype=np.uint8)
    band = linstr(band)
    return np.asarray(band, np.uint8).tostring()


def lin2pc(band):
    #  2% linear stretch
    band = np.fromstring(band, dtype=np.uint8)
    band = lin2pcstr(band)
    return np.asarray(band, np.uint8).tostring()


def histeq(band):
    #  histogram equalization
    band = np.fromstring(band, dtype=np.uint8)
    band = histeqstr(band)
    return np.asarray(band, np.uint8).tostring()


def stretch(redband, greenband, blueband, enhance):
    if enhance == 'linear2pc':
        return (lin2pc(redband), lin2pc(greenband), lin2pc(blueband))
    elif enhance == "equalization":
        return (histeq(redband), histeq(greenband), histeq(blueband))
    elif enhance == 'linear':
        return (lin(redband), lin(greenband), lin(blueband))
    else:  # do nothing
        return (redband, greenband, blueband)


def byte_stretch(band, dtype=1, rng=None):
    #  byte stretch an image band coded as string
    if dtype == 1:
        tmp = np.fromstring(band, dtype=np.uint8)
    elif dtype == 2:
        tmp = np.fromstring(band, dtype=np.uint16)
    elif dtype == 4:
        tmp = np.fromstring(band, dtype=np.float32)
    else:
        tmp = np.fromstring(band, dtype=np.float64)
    if rng is None:
        rng = [np.min(tmp), np.max(tmp)]
    tmp = (tmp - rng[0]) * 255.0 / (rng[1] - rng[0])
    tmp = np.where(tmp < 0, 0, tmp)
    tmp = np.where(tmp > 255, 255, tmp)
    return np.asarray(tmp, np.uint8).tostring()


def normalize(da, coeffs):
    # return normalized bsq image string from input data array
    result = ''
    for k in range(da.bands):
        b = coeffs[k, 0]
        a = coeffs[k, 1]
        g = da.data[:, k] * b + a
        g = np.where(g < 0, 0, g)
        g = np.where(g > 255, 255, g)
        g = np.asarray(g, np.uint8)
        result += g.tostring()
    return result


def byteStretch(arr, rng=None):
    #  byte stretch image numpy array
    shp = arr.shape
    arr = arr.ravel()
    if rng is None:
        rng = [np.min(arr), np.max(arr)]
    tmp = (arr - rng[0]) * 255.0 / (rng[1] - rng[0])
    tmp = np.where(tmp < 0, 0, tmp)
    tmp = np.where(tmp > 255, 255, tmp)
    return np.asarray(np.reshape(tmp, shp), np.uint8)


# --------------
# simple dialogs
# --------------

def select_directory(title=None):
    root = Tkinter.Tk()
    root.withdraw()
    d = tkFileDialog.askdirectory(title=title)
    root.destroy()
    if d == '':
        return None
    else:
        return d


def select_infile(filt=None, title=None, mask=None):
    root = Tkinter.Tk()
    root.withdraw()
    if filt is None:
        filetypes = [('anyfile', '*_*')]
    else:
        filetypes = [('filtered', filt)]
    filename = tkFileDialog.Open(filetypes=filetypes, title=title).show()
    root.destroy()
    if filename == '':
        return None
    if mask:
        root = Tkinter.Tk()
        root.withdraw()
        filetypes = [('anyfile', '*.*')]
        maskname = tkFileDialog.Open(filetypes=filetypes, title='associated mask').show()
        root.destroy()
        if maskname:
            return (filename, maskname)
        else:
            return (filename, None)
    else:
        return filename


def select_outfilefmt(title=''):
    root = Tkinter.Tk()
    root.withdraw()
    fmt = tkSimpleDialog.askstring(title + ' Output Format',
                                   'Enter one of GTiff, PCIDSK, HFA, ENVI',
                                   initialvalue='GTiff')
    if fmt == 'GTiff':
        filetypes = [('GeoTiff', '*.tif')]
        defaultextension = '.tif'
    elif fmt == 'PCIDSK':
        filetypes = [('PCI', '*.pix')]
        defaultextension = '.pix'
    elif fmt == 'HFA':
        filetypes = [('ERDAS Imagine', '*.img')]
        defaultextension = '.img'
    elif fmt == 'ENVI':
        filetypes = [('ENVI', '*.*')]
        defaultextension = None
    else:
        root.destroy()
        return (None, None)
    filename = tkFileDialog.SaveAs(filetypes=filetypes, defaultextension=defaultextension, title=title).show()
    root.destroy()
    if filename:
        return (filename, fmt)
    else:
        return (None, None)


def select_outfile(filt='*', title=''):
    root = Tkinter.Tk()
    root.withdraw()
    filetypes = [('Output', filt)]
    filename = tkFileDialog.SaveAs(filetypes=filetypes, defaultextension=filt, title=title).show()
    root.destroy()
    if filename:
        return filename
    else:
        return None


def select_pos(bands, onlyone=None):
    root = Tkinter.Tk()
    root.withdraw()
    if onlyone:
        pos = tkSimpleDialog.askstring('Pos',
                                       'Single band position, 1 - %i' % bands,
                                       initialvalue='1')
    else:
        pos = tkSimpleDialog.askstring('Pos',
                                       'Band positions as list',
                                       initialvalue=str(range(1, bands + 1)))
    root.destroy()
    if pos:
        return eval(pos)
    else:
        return None


def select_dims(dimensions):
    root = Tkinter.Tk()
    root.withdraw()
    dims = tkSimpleDialog.askstring('Dims',
                                    'Dimensions as list [x0,y0,cols,rows]',
                                    initialvalue=str(dimensions))
    root.destroy()
    if dims:
        return tuple(eval(dims))
    else:
        return None


def select_penal(lam):
    root = Tkinter.Tk()
    root.withdraw()
    lam = tkSimpleDialog.askstring('Pen',
                                   'Penalization',
                                   initialvalue=str(lam))
    root.destroy()
    if lam:
        return eval(lam)
    else:
        return None


def select_ncp(ncp):
    root = Tkinter.Tk()
    root.withdraw()
    ncp = tkSimpleDialog.askstring('NCP',
                                   'Probability threshold',
                                   initialvalue=str(ncp))
    root.destroy()
    if ncp:
        return eval(ncp)
    else:
        return None


def select_rgb(bands):
    root = Tkinter.Tk()
    root.withdraw()
    if bands == 1:
        iv = '[1,1,1]'
    elif bands == 2:
        iv = '[1,1,2]'
    else:
        iv = '[1,2,3]'
    rgb = tkSimpleDialog.askstring('RGB',
                                   'RGB bands as sublist from [1 ... %i]' % bands,
                                   initialvalue=iv)
    root.destroy()
    if rgb:
        return tuple(eval(rgb))
    else:
        return None


def select_enhance(enh):
    root = Tkinter.Tk()
    root.withdraw()
    enh = tkSimpleDialog.askstring('Enhance',
                                   '1=linear255, 2=linear, 3=linear2pc, 4=equalization',
                                   initialvalue=enh)
    root.destroy()
    if enh:
        return enh
    else:
        return None


def select_integer(L, msg='Enter a number'):
    root = Tkinter.Tk()
    root.withdraw()
    L = tkSimpleDialog.askstring('Enter a number', msg, initialvalue=str(L))
    root.destroy()
    if L:
        return eval(L)
    else:
        return None

    # ------------------------


# generalized eigenproblem
# ------------------------

def choldc(A):
    # Cholesky-Banachiewicz algorithm,
    # A is a numpy matrix
    L = A - A
    for i in range(len(L)):
        for j in range(i):
            sm = 0.0
            for k in range(j):
                sm += L[i, k] * L[j, k]
            L[i, j] = (A[i, j] - sm) / L[j, j]
        sm = 0.0
        for k in range(i):
            sm += L[i, k] * L[i, k]
        L[i, i] = math.sqrt(A[i, i] - sm)
    return L


def geneiv(A, B):
    # solves A*x = lambda*B*x for numpy matrices A and B,
    # returns eigenvectors in columns
    Li = np.linalg.inv(choldc(B))
    C = Li * A * (Li.transpose())
    C = np.asmatrix((C + C.transpose()) * 0.5, np.float32)
    eivs, V = np.linalg.eig(C)
    return eivs, Li.transpose() * V


# ------------------------------------------------
# spectral transformations, use DataArray objects,
# return results as float32 image blobs
# ------------------------------------------------

def pca(da):
    try:
        means, covmat = da.covw()
        lams, eivs = np.linalg.eigh(covmat)
        #      sort eivs in decreasing order of variance
        idx = (np.argsort(lams))[::-1]
        lams = lams[idx]
        eivs = np.mat(eivs[:, idx], np.float32)
        #      project in bsq format, single precision, flattened
        mns = np.tile(means, (da.pixels, 1))
        ac = np.mat(da.data - mns)
        mns = None
        #      PCs in bsq format
        pcs = (ac * eivs).T
        pcs = pcs.ravel().tostring()
        return (lams, pcs)
    except:
        return None


def mnf(da, samples, lines, bands):
    try:
        #      center
        means, S = da.covw()
        mns = np.asarray(np.tile(means, (da.pixels, 1)))
        img = da.data - mns
        mns = None
        #      shifted image
        imgs = np.reshape(img, (samples, lines, bands))
        imgn = imgs - (np.roll(imgs, 1, axis=0) + np.roll(imgs, 1, axis=1)) / 2
        #      noise covariance
        dan = DataArray(imgn.tostring(), samples, lines, bands, 'bip', 4)
        Sn = dan.covw()[1] / 2
        #      mnf, eigenvalues sorted in increasing order
        lams, eivs = geneiv(Sn, S)
        idx = (np.argsort(lams))
        lams = lams[idx]
        eivs = (eivs[:, idx]).T
        #      MNFs in bsq format
        mnfs = eivs * np.mat(img).T
        mnfs = mnfs.ravel().tostring()
        return (lams, mnfs)
    except:
        return None


def similarity(bn0, bn1):
    """Register bn1 to bn0 ,  M. Canty 2012
bn0, bn1 and returned result are image bands
Modified from Imreg.py, see http://www.lfd.uci.edu/~gohlke/:
 Copyright (c) 2011-2012, Christoph Gohlke
 Copyright (c) 2011-2012, The Regents of the University of California
 Produced at the Laboratory for Fluorescence Dynamics
 All rights reserved.
    """

    def highpass(shape):
        """Return highpass filter to be multiplied with fourier transform."""
        x = np.outer(
            np.cos(np.linspace(-math.pi / 2., math.pi / 2., shape[0])),
            np.cos(np.linspace(-math.pi / 2., math.pi / 2., shape[1])))
        return (1.0 - x) * (2.0 - x)

    def logpolar(image, angles=None, radii=None):
        """Return log-polar transformed image and log base."""
        shape = image.shape
        center = shape[0] / 2, shape[1] / 2
        if angles is None:
            angles = shape[0]
            if radii is None:
                radii = shape[1]
        theta = np.empty((angles, radii), dtype=np.float64)
        theta.T[:] = -np.linspace(0, np.pi, angles, endpoint=False)
        #      d = radii
        d = np.hypot(shape[0] - center[0], shape[1] - center[1])
        log_base = 10.0 ** (math.log10(d) / (radii))
        radius = np.empty_like(theta)
        radius[:] = np.power(log_base, np.arange(radii, dtype = np.float64)) - 1.0
        x = radius * np.sin(theta) + center[0]
        y = radius * np.cos(theta) + center[1]
        output = np.empty_like(x)
        ndii.map_coordinates(image, [x, y], output=output)
        return output, log_base

    lines0, samples0 = bn0.shape
    #  make reference and warp bands same shape
    bn1 = bn1[0:lines0, 0:samples0]
    #  get scale, angle
    f0 = fftshift(abs(fft2(bn0)))
    f1 = fftshift(abs(fft2(bn1)))
    h = highpass(f0.shape)
    f0 *= h
    f1 *= h
    del h
    f0, log_base = logpolar(f0)
    f1, log_base = logpolar(f1)
    f0 = fft2(f0)
    f1 = fft2(f1)
    r0 = abs(f0) * abs(f1)
    ir = abs(ifft2((f0 * f1.conjugate()) / r0))
    i0, i1 = np.unravel_index(np.argmax(ir), ir.shape)
    angle = 180.0 * i0 / ir.shape[0]
    scale = log_base ** i1
    if scale > 1.8:
        ir = abs(ifft2((f1 * f0.conjugate()) / r0))
        i0, i1 = np.unravel_index(np.argmax(ir), ir.shape)
        angle = -180.0 * i0 / ir.shape[0]
        scale = 1.0 / (log_base ** i1)
        if scale > 1.8:
            raise ValueError("Images are not compatible. Scale change > 1.8")
    if angle < -90.0:
        angle += 180.0
    elif angle > 90.0:
        angle -= 180.0
    #  re-scale and rotate and then get shift
    bn2 = ndii.zoom(bn1, 1.0 / scale)
    bn2 = ndii.rotate(bn2, angle)
    if bn2.shape < bn0.shape:
        t = np.zeros_like(bn0)
        t[:bn2.shape[0], :bn2.shape[1]] = bn2
        bn2 = t
    elif bn2.shape > bn0.shape:
        bn2 = bn2[:bn0.shape[0], :bn0.shape[1]]
    f0 = fft2(bn0)
    f1 = fft2(bn2)
    ir = abs(ifft2((f0 * f1.conjugate()) / (abs(f0) * abs(f1))))
    t0, t1 = np.unravel_index(np.argmax(ir), ir.shape)
    if t0 > f0.shape[0] // 2:
        t0 -= f0.shape[0]
    if t1 > f0.shape[1] // 2:
        t1 -= f0.shape[1]
    #  return result
    return (scale, angle, [t0, t1])


# ---------------------------
# discrete wavelet transform
# ---------------------------

class DWTArray(object):
    #Partial DWT representation of image band which is input as 2-D uint8 array

    def __init__(self, band, samples, lines, itr=0):
        # Daubechies D4 wavelet
        self.H = np.asarray([(1 - math.sqrt(3)) / 8, (3 - math.sqrt(3)) / 8, (3 + math.sqrt(3)) / 8, (1 + math.sqrt(3)) / 8])
        self.G = np.asarray([-(1 + math.sqrt(3)) / 8, (3 + math.sqrt(3)) / 8, -(3 - math.sqrt(3)) / 8, (1 - math.sqrt(3)) / 8])
        self.num_iter = itr
        self.max_iter = 3
        # ignore edges if band dimension is not divisible by 2^max_iter
        r = 2 ** self.max_iter
        self.samples = r * (samples // r)
        self.lines = r * (lines // r)
        self.data = np.asarray(band[:self.lines, :self.samples], np.float32)

    def get_quadrant(self, quadrant):
        if self.num_iter == 0:
            m = 2 * self.lines
            n = 2 * self.samples
        else:
            m = self.lines / 2 ** (self.num_iter - 1)
            n = self.samples / 2 ** (self.num_iter - 1)
        if quadrant == 0:
            f = self.data[:m / 2, :n / 2]
        elif quadrant == 1:
            f = self.data[:m / 2, n / 2:n]
        elif quadrant == 2:
            f = self.data[m / 2:m, :n / 2]
        else:
            f = self.data[m / 2:m, n / 2:n]
        f = np.where(f < 0, 0, f)
        f = np.where(f > 255, 255, f)
        return np.asarray(f, np.uint8)

    def put_quadrant(self, f1, quadrant):
        if not (quadrant in range(4)) or (self.num_iter == 0):
            return 0
        m = self.lines / 2 ** (self.num_iter - 1)
        n = self.samples / 2 ** (self.num_iter - 1)
        f0 = self.data
        if quadrant == 0:
            f0[:m / 2, :n / 2] = f1
        elif quadrant == 1:
            f0[:m / 2, n / 2:n] = f1
        elif quadrant == 2:
            f0[m / 2:m, :n / 2] = f1
        else:
            f0[m / 2:m, n / 2:m] = f1
        return 1

    def normalize(self, a, b):
        #      normalize wavelet coefficients at all levels
        for c in range(1, self.num_iter + 1):
            m = self.lines / (2 ** c)
            n = self.samples / (2 ** c)
            self.data[:m, n:2 * n] = a[0] * self.data[:m, n:2 * n] + b[0]
            self.data[m:2 * m, :n] = a[1] * self.data[m:2 * n, :n] + b[1]
            self.data[m:2 * m, n:2 * n] = a[2] * self.data[m:2 * n, n:2 * n] + b[2]

    def filter(self):
        #      single application of filter bank
        if self.num_iter == self.max_iter:
            return 0
        #      get upper left quadrant
        m = self.lines / 2 ** self.num_iter
        n = self.samples / 2 ** self.num_iter
        f0 = self.data[:m, :n]
        #      temporary arrays
        f1 = np.zeros((m / 2, n))
        g1 = np.zeros((m / 2, n))
        ff1 = np.zeros((m / 2, n / 2))
        fg1 = np.zeros((m / 2, n / 2))
        gf1 = np.zeros((m / 2, n / 2))
        gg1 = np.zeros((m / 2, n / 2))
        #      filter columns and downsample
        ds = np.asarray(range(m / 2)) * 2 + 1
        for i in range(n):
            temp = np.convolve(f0[:, i].ravel(), \
                               self.H, 'same')
            f1[:, i] = temp[ds]
            temp = np.convolve(f0[:, i].ravel(), \
                               self.G, 'same')
            g1[:, i] = temp[ds]
        #      filter rows and downsample
        ds = np.asarray(range(n / 2)) * 2 + 1
        for i in range(m / 2):
            temp = np.convolve(f1[i, :], self.H, 'same')
            ff1[i, :] = temp[ds]
            temp = np.convolve(f1[i, :], self.G, 'same')
            fg1[i, :] = temp[ds]
            temp = np.convolve(g1[i, :], self.H, 'same')
            gf1[i, :] = temp[ds]
            temp = np.convolve(g1[i, :], self.G, 'same')
            gg1[i, :] = temp[ds]
        f0[:m / 2, :n / 2] = ff1
        f0[:m / 2, n / 2:] = fg1
        f0[m / 2:, :n / 2] = gf1
        f0[m / 2:, n / 2:] = gg1
        self.data[:m, :n] = f0
        self.num_iter = self.num_iter + 1

    def invert(self):
        H = self.H[::-1]
        G = self.G[::-1]
        m = self.lines / 2 ** (self.num_iter - 1)
        n = self.samples / 2 ** (self.num_iter - 1)
        #      get upper left quadrant
        f0 = self.data[:m, :n]
        ff1 = f0[:m / 2, :n / 2]
        fg1 = f0[:m / 2, n / 2:]
        gf1 = f0[m / 2:, :n / 2]
        gg1 = f0[m / 2:, n / 2:]
        f1 = np.zeros((m / 2, n))
        g1 = np.zeros((m / 2, n))
        #      upsample and filter rows
        for i in range(m / 2):
            a = np.ravel(np.transpose(np.vstack((ff1[i, :], np.zeros(n / 2)))))
            b = np.ravel(np.transpose(np.vstack((fg1[i, :], np.zeros(n / 2)))))
            f1[i, :] = np.convolve(a, H, 'same') + np.convolve(b, G, 'same')
            a = np.ravel(np.transpose(np.vstack((gf1[i, :], np.zeros(n / 2)))))
            b = np.ravel(np.transpose(np.vstack((gg1[i, :], np.zeros(n / 2)))))
            g1[i, :] = np.convolve(a, H, 'same') + np.convolve(b, G, 'same')
        #      upsample and filter columns
        for i in range(n):
            a = np.ravel(np.transpose(np.vstack((f1[:, i], np.zeros(m / 2)))))
            b = np.ravel(np.transpose(np.vstack((g1[:, i], np.zeros(m / 2)))))
            f0[:, i] = 4 * (np.convolve(a, H, 'same') + np.convolve(b, G, 'same'))
        self.data[:m, :n] = f0
        self.num_iter = self.num_iter - 1


class ATWTArray(object):
    '''A trous wavelet transform'''

    def __init__(self, band):
        self.num_iter = 0
        #      cubic spline filter
        self.H = np.array([1.0 / 16, 1.0 / 4, 3.0 / 8, 1.0 / 4, 1.0 / 16])
        #      data arrays
        self.lines, self.samples = band.shape
        self.bands = np.zeros((4, self.lines, self.samples), np.float32)
        self.bands[0, :, :] = np.asarray(band, np.float32)

    def inject(self, band):
        m = self.lines
        n = self.samples
        self.bands[0, :, :] = band[0:m, 0:n]

    def get_band(self, i):
        return self.bands[i, :, :]

    def normalize(self, a, b):
        if self.num_iter > 0:
            for i in range(1, self.num_iter + 1):
                self.bands[i, :, :] = a * self.bands[i, :, :] + b

    def filter(self):
        if self.num_iter < 3:
            self.num_iter += 1
            #          a trous filter
            n = 2 ** (self.num_iter - 1)
            H = np.vstack((self.H, np.zeros((2 ** (n - 1), 5))))
            H = np.transpose(H).ravel()
            H = H[0:-n]
            #          temporary arrays
            f1 = np.zeros((self.lines, self.samples))
            ff1 = f1 * 0.0
            #          filter columns
            f0 = self.bands[0, :, :]
            #          filter columns
            for i in range(self.samples):
                f1[:, i] = np.convolve(f0[:, i].ravel(), H, 'same')
            #          filter rows
            for j in range(self.lines):
                ff1[j, :] = np.convolve(f1[j, :], H, 'same')
            self.bands[self.num_iter, :, :] = self.bands[0, :, :] - ff1
            self.bands[0, :, :] = ff1

    def invert(self):
        if self.num_iter > 0:
            self.bands[0, :, :] += self.bands[self.num_iter, :, :]
            self.num_iter -= 1


if __name__ == '__main__':
    pass
