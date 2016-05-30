
import os
import spatial
import vector
import srtm
from datetime import datetime
from time import asctime, mktime, strptime
from ancillary import finder, multicore
from gamma.util import geocode
from pyroSAR import identify

sitename = "Camargue"
maindir = "/geonfs01_vol1/ve39vem/swos_process"
srtmdir = "/geonfs02_vol1/SRTM_1_HGT"

targetres = 20
func_geoback = 2
func_interp = 0
scaling = ["linear", "db"]

tempdir = os.path.join(maindir, sitename, "proc_in")
outdir = os.path.join(maindir, sitename, "proc_out")

for dir in [tempdir, outdir]:
    if not os.path.isdir(dir):
        os.makedirs(dir)

site = vector.Vector("/geonfs01_vol1/ve39vem/swos_testsites/Test_Sites_project_v29_{}.shp".format(sitename))
site_area = site.getArea()

scenes = finder("/geonfs02_vol1/swos_fsu/sentinel1/GRD", ["^S1A"], regex=True)
# scenes = finder("/geonfs01_vol1/ve39vem/swos", ["*zip"])
#
selection = []
a = datetime.now()
for scene in scenes:
    try:
        id = identify(scene)
        intersect = spatial.intersect(id.bbox(), site)
        if intersect is not None:
            if "VV" in id.polarisations and id.beam == "IW" and id.start > "20160302T054342":
            # if intersect.GetArea()/camarque_area == 1 and id.beam == "IW" and id.pols in ["SV", "DV"]:
            # if intersect.GetArea()/camarque_area == 1:
            # if id.beam == "IW" and id.pols in ["SV", "DV"]:
                print scene
                selection.append(scene)
    except IOError:
        continue

b = datetime.now()
print b-a
# with open("/geonfs01_vol1/ve39vem/scenelist", "w") as outfile:
#     outfile.write("\n".join(selection)+"\n")

# with open("/geonfs01_vol1/ve39vem/scenelist", "r") as infile:
#     selection = infile.read().strip().split()

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

# scene = selection[0]
# geocode(scene, srtm_mosaic_utm, tempdir, outdir, targetres, scaling, func_geoback, func_interp)

print "start geocoding"
multicore(geocode, cores=10, multiargs={"scene": selection}, dem=srtm_mosaic_utm,
          tempdir=tempdir, outdir=outdir,
          targetres=targetres, scaling=scaling,
          func_geoback=func_geoback, func_interp=func_interp)
print datetime.now()
print "done"
