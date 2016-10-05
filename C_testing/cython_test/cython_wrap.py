from time import asctime

import numpy as np
from C_testing.cython_test import fancymath

from spatial import raster

filename = "/media/john/Data/DATA/Viterbo_Italy_test_layerstack_subset"

inarr = raster.Raster(filename).raster.ReadAsArray()

print asctime()
# outarr1 = np.zeros([inarr.shape[1], inarr.shape[2]], dtype=np.float32)
# for x in range(inarr.shape[1]):
#     for y in range(inarr.shape[2]):
#         outarr1[x, y] = max(inarr[:, x, y]) - min(inarr[:, x, y])
outarr1 = np.mean(inarr, axis=0)
print asctime()

# print inarr.shape
# print inarr[:, 50, 50]
# print len(inarr[:, 50, 50])
# print max(inarr[:, 50, 50])

# print inarr.shape
# print inarr.dtype
#
outarr2 = fancymath(inarr)

print asctime()
