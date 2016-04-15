
import os
from mosaic import mosaic
from mosaic_aux import finder

dir_in = "/homes4/geoinf/ve39vem/Landsat"
dir_out = "/homes4/geoinf/ve39vem/Landsat/outdir"

mosaik = os.path.join(dir_out, "mosaic_test")
# master = "/homes4/geoinf/ve39vem/Landsat/LT51130262010246.envi"
master = "/homes4/geoinf/ve39vem/Landsat/LT51180272005219.envi"

master = mosaik if os.path.isfile(mosaik) else master

slaves = [x for x in finder(dir_in, ["*.envi"]) if x != master and "2005" in x]

mosaic(master, slaves, mosaik)
