import os
import re
import spatial
import vector
import srtm
from datetime import datetime
# from time import asctime, mktime, strptime
from ancillary import finder, multicore
from gamma.util import geocode
from pyroSAR import identify

sitename = "Azraq"
maindir = "/geonfs01_vol1/ve39vem/swos_process"
srtmdir = "/geonfs02_vol1/SRTM_1_HGT"

targetres = 20
func_geoback = 2
func_interp = 0
scaling = ["linear", "db"]

sitedir = os.path.join(maindir, sitename)
tempdir = os.path.join(sitedir, "proc_in")
outdir = os.path.join(sitedir, "proc_out")

for dir in [tempdir, outdir]:
    if not os.path.isdir(dir):
        os.makedirs(dir)

site = vector.Vector("/geonfs01_vol1/ve39vem/swos_testsites/Test_Sites_project_v29_{}.shp".format(sitename))
site_area = site.getArea()

scenes = finder("/geonfs02_vol1/swos_fsu/sentinel1/GRD", ["^S1[AB]"], regex=True)
# scenes = finder("/geonfs01_vol1/ve39vem/swos", ["*zip"])

# selection = []
# a = datetime.now()
# for scene in scenes:
#     try:
#         id = identify(scene)
#         if len(finder(outdir, [id.outname_base()], regex=True)) == 0:
#             intersect = spatial.intersect(id.bbox(), site)
#             if intersect is not None:
#                 if "VV" in id.polarizations and id.acquisition_mode == "IW":
#                     # if id.start > "20160302T054342":
#                     # if intersect.GetArea()/site_area == 1:
#                     # if id.beam == "IW" and id.pols in ["SV", "DV"]:
#                     print scene
#                     selection.append(id)
#     except IOError:
#         continue
#
# b = datetime.now()
# print b - a
#
# with open(os.path.join(sitedir, "scenelist"), "w") as outfile:
#     outfile.write("\n".join([x.scene for x in selection])+"\n")
#
with open(os.path.join(sitedir, "scenelist"), "r") as infile:
    selection = [identify(x) for x in infile.read().strip().split()]

# scenes = [identify(x) for x in selection]
#
#
# def seconds(scene): return mktime(strptime(scene.start, "%Y%m%dT%H%M%S"))
#
# scenes = sorted(scenes, key=seconds)
#
# overlaps = []
# for item in scenes:
#     intersect = spatial.intersect(item.bbox(), camarque)
#     overlaps.append(intersect.GetArea()/camarque_area)
#
# for i in range(1, len(scenes)):
#     print scenes[i].file, seconds(scenes[i])-seconds(scenes[i-1])


srtm_mosaic = os.path.join(maindir, sitename, "srtm_{}".format(sitename))
srtm_mosaic_utm = srtm_mosaic + "_utm"

if not os.path.isfile(srtm_mosaic):
    srtm.makeSRTM(selection, srtmdir, srtm_mosaic)

if not os.path.isfile(srtm_mosaic_utm):
    srtm.transform(srtm_mosaic, srtm_mosaic_utm, targetres)

# scene = selection[1]
# geocode(scene, srtm_mosaic_utm, tempdir, outdir, targetres, scaling, func_geoback, func_interp)

print "start geocoding"
print datetime.now()
multicore(geocode, cores=10, multiargs={"scene": selection}, dem=srtm_mosaic_utm,
          tempdir=tempdir, outdir=outdir,
          targetres=targetres, scaling=scaling,
          func_geoback=func_geoback, func_interp=func_interp)
print datetime.now()
print "done"
