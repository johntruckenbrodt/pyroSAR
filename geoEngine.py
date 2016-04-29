
from pyroSAR import identify
from gamma.util import geocode

srtmdir = "/geonfs02_vol1/SRTM_1_HGT"
scenelist = "/geonfs01_vol1/ve39vem/scenelist"

tempdir = "/geonfs01_vol1/ve39vem/temp"
outdir = "/geonfs01_vol1/ve39vem/out"

targetres = 20
func_geoback = 2
func_interp = 0
scaling = "db"

with open(scenelist, "r") as infile:
    files = infile.read().strip().split()

scenes = [identify(x) for x in files]

scenes = [x for x in scenes if "VV" in x.polarisations]

# srtm_mosaic = os.path.join(tempdir, "srtm")
#
# srtm.makeSRTM(scenes, srtmdir, srtm_mosaic)
#
# srtm.transform(srtm_mosaic, srtm_mosaic+"_utm", targetres)
#
# srtm_mosaic += "_utm"

srtm_mosaic = "/geonfs01_vol1/ve39vem/temp/srtm_swap_fill_utm"

scene = scenes[0]

geocode(scene, srtm_mosaic, tempdir, outdir, targetres, scaling, func_geoback, func_interp)
