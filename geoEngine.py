
import os
import srtm
import vector
import spatial
from datetime import datetime
from pyroSAR import identify
from gamma.util import geocode
from ancillary import finder, multicore


srtmdir = "/geonfs02_vol1/SRTM_1_HGT"
scenelist = "/geonfs01_vol1/ve39vem/scenelist"

tempdir = "/geonfs01_vol1/ve39vem/temp"
outdir = "/geonfs01_vol1/ve39vem/out"

targetres = 20
func_geoback = 2
func_interp = 0
scaling = ["linear", "db"]


camarque = vector.Vector("/geonfs01_vol1/ve39vem/S1/test_camarque/boundary/Test_Sites_project_v11_Camarque.shp")
camarque_area = camarque.getArea()

scenes = finder("/geonfs02_vol1/swos_fsu/sentinel1/GRD", ["^S1[AB]"], regex=True)
# scenes = finder("/geonfs01_vol1/ve39vem/swos", ["*zip"])

selection = []
print datetime.now()
print "collecting files"
for scene in scenes:
    try:
        id = identify(scene)
        intersect = spatial.intersect(id.bbox(), camarque)
        if intersect is not None:
            if "VV" in id.polarisations and id.beam == "IW" and id.start > "20160302T054342":
            # if intersect.GetArea()/camarque_area == 1 and id.beam == "IW" and id.pols in ["SV", "DV"]:
            # if intersect.GetArea()/camarque_area == 1:
            # if id.beam == "IW" and id.pols in ["SV", "DV"]:
                print scene
                selection.append(id)
    except IOError:
        continue
print datetime.now()
# with open(scenelist, "r") as infile:
#     files = infile.read().strip().split()


srtm_mosaic = os.path.join(tempdir, "srtm_camarque")
srtm_mosaic_utm = srtm_mosaic + "_utm"

# srtm.makeSRTM(selection, srtmdir, srtm_mosaic)
# srtm.transform(srtm_mosaic, srtm_mosaic_utm, targetres)

# scene = selection[0]
# geocode(scene, srtm_mosaic_utm, tempdir, outdir, targetres, scaling, func_geoback, func_interp)

print "start geocoding"
multicore(geocode, cores=10, multiargs={"scene": selection}, dem=srtm_mosaic_utm,
          tempdir=tempdir, outdir=outdir,
          targetres=targetres, scaling=scaling,
          func_geoback=func_geoback, func_interp=func_interp)
print datetime.now()
print "done"
