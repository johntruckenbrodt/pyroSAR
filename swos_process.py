
import re
import spatial
import vector
from datetime import datetime
from time import asctime, mktime, strptime
from ancillary import finder, multicore

from S1 import S1_geoGRD
from pyroSAR import identify

camarque = vector.Vector("/geonfs01_vol1/ve39vem/S1/test_camarque/boundary/Test_Sites_project_v11_Camarque.shp")
camarque_area = camarque.getArea()

# scenes = finder("/geonfs02_vol1/swos_fsu/sentinel1/GRD", ["^S1A"], regex=True)
# scenes = finder("/geonfs01_vol1/ve39vem/swos", ["*zip"])
#
# selection = []
# a = datetime.now()
# for scene in scenes:
#     try:
#         id = identify(scene)
#         intersect = spatial.intersect(id.bbox(), camarque)
#         if intersect is not None:
#             # if intersect.GetArea()/camarque_area == 1 and id.beam == "IW" and id.pols in ["SV", "DV"]:
#             # if intersect.GetArea()/camarque_area == 1:
#             # if id.beam == "IW" and id.pols in ["SV", "DV"]:
#             print scene
#             selection.append(scene)
#     except IOError:
#         continue
# b = datetime.now()
# print b-a
# with open("/geonfs01_vol1/ve39vem/scenelist", "w") as outfile:
#     outfile.write("\n".join(selection)+"\n")

with open("/geonfs01_vol1/ve39vem/scenelist", "r") as infile:
    selection = infile.read().strip().split()

scenes = [identify(x) for x in selection]


def seconds(scene): return mktime(strptime(scene.start, "%Y%m%dT%H%M%S"))

scenes = sorted(scenes, key=seconds)

overlaps = []
for item in scenes:
    intersect = spatial.intersect(item.bbox(), camarque)
    overlaps.append(intersect.GetArea()/camarque_area)

for i in range(1, len(scenes)):
    print scenes[i].file, seconds(scenes[i])-seconds(scenes[i-1])


# print "started:", asctime()
# print "processing a total of {} scenes".format(len(selection))
# tempdir = "/geonfs01_vol1/ve39vem/S1/test_camarque/test_in"
# outdir = "/geonfs01_vol1/ve39vem/S1/test_camarque/test_out"
# srtmdir = "/geonfs02_vol1/SRTM_1_HGT"

# multicore(S1_geoGRD.main, cores=10, multiargs={"zipfile": selection}, tempdir=tempdir, outdir=outdir, srtmdir=srtmdir, transform=True)
# print "finished:", asctime()

# testcase = "/geonfs02_vol1/swos_fsu/sentinel1/GRD/S1A_IW_GRDH_1SDV_20150314T173821_20150314T173846_005031_0064FC_DBE3.zip"
# S1_geoGRD.main(testcase, tempdir=tempdir, outdir=outdir, srtmdir=srtmdir, transform=True)
