
import os
import math
import ctypes
import platform
import numpy as np
from glob import glob
from scipy.special import betainc
from numpy.ctypeslib import ndpointer as ndp

# provide C-implemented provisional means algorithm
provmeans = ctypes.cdll.LoadLibrary(os.path.join(os.getcwd(), "prov_means.dll" if platform.system() == "Windows" else "libprov_means.so")).provmeans
provmeans.argtypes = [ndp(np.float64), ndp(np.float64), ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_double), ndp(np.float64), ndp(np.float64)]
provmeans.restype = None


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


class Cpm(object):
    # Provisional means algorithm

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


def finder(folder, matchlist):
    out = list(set([f for files in [glob(os.path.join(item[0], pattern)) for item in os.walk(folder) for pattern in matchlist] for f in files]))
    return sorted(out)


# F-test for equal variance
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
        return f, 2.0 - prob
    else:
        return f, prob


def geneiv(A, B):
    # solves A*x = lambda*B*x for numpy matrices A and B,
    # returns eigenvectors in columns
    Li = np.linalg.inv(choldc(B))
    C = Li * A * (Li.transpose())
    C = np.asmatrix((C + C.transpose()) * 0.5, np.float32)
    eivs, V = np.linalg.eig(C)
    return eivs, Li.transpose() * V


# orthogonal regression
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
