import os
import re
from time import asctime

from ancillary import finder, run
from spatial import haversine, raster

# path = "/media/john/Data/DATA/Sentinel"
# path_tiling = "/media/john/Data/DATA/Sentinel/tiling_test"
# path = "/geonfs02_vol1/baci_fsu/02_Sentinel/03_FTS/02_Viterbo_Italy/04_output"
path = "/homes4/geoinf/ve39vem/RADAR/Sentinel/test_out"
path_tiling = "/homes4/geoinf/ve39vem/RADAR/Sentinel/tiling"
files = finder(path, ["*vv*.tif"])

for filename in files:
    print os.path.basename(filename), asctime()

    data = raster.Raster(filename)

    lat = [data.geo["ymin"], data.geo["ymax"]]
    lon = [data.geo["xmin"], data.geo["xmax"]]

    lat = [int(x//1) for x in lat]
    lon = [int(x//1) for x in lon]

    lat = range(min(lat), max(lat)+1)
    lon = range(min(lon), max(lon)+1)

    resampling = "bilinear"
    interleave = "BAND"

    # targetres = data.res
    warpoptions = ["gdalwarp", "-overwrite", "--config", "GDAL_CACHEMAX", 2000, "-wm", 6000, "-multi",
                   "-r", resampling, "-of", "GTiff", "-co", "INTERLEAVE="+interleave]

    for x in lon:
        for y in lat:
            tilename = "{0}{1}{2}{3}".format("N" if y > 0 else "S", str(y).zfill(2), "W" if x < 0 else "E", str(x).zfill(3))

            outpath = os.path.join(path_tiling, tilename)

            lat_center = y+.5
            lon_center = x+.5
            post_north = haversine(lat_center, lon_center, lat_center+data.res[1], lon_center)
            post_east = haversine(lat_center, lon_center, lat_center, lon_center+data.res[0])

            targetres = ["-tr", data.res[0]*(post_north/post_east), data.res[1]]

            if not os.path.isdir(outpath):
                os.makedirs(outpath)

            timestamp = re.search("[0-9]{8}_[0-9]{6}", filename).group().replace("_", "T")
            # dtime = datetime.strptime(timestamp, "%Y%m%dT%H%M%S")

            outname = os.path.join(outpath, tilename+"_"+timestamp+".tif")
            extent_arg = ["-te", x, y, x+1, y+1]

            # print "resample:", asctime()
            run([warpoptions, extent_arg, targetres, filename, outname])

            stats = raster.Raster(outname).allstats[0]

            # check whether resulting file contains valid pixels; if not delete it otherwise reduce it to lines and columns containing any number of valid data entries
            if stats is None:
                for item in finder(os.path.dirname(outname), [os.path.basename(outname)+"*"]):
                    os.remove(item)
            # else:
            #     # print "load raster:", asctime()
            #     rast = raster.Raster(outname)
            #     # print "reduce:", asctime()
            #     rast.reduce()
            #
            #     def fun(x):
            #         return x*1000
            #
            #     # print "rescale:", asctime()
            #     rast.rescale(fun)
            #     # print "write:", asctime()
            #     rast.write(outname=outname)
                # rast.write(outname=outname, dtype="Int16")
                # print "done:", asctime()
                # print "--------"

            if os.listdir(outpath) == []:
                os.rmdir(outpath)


def mosaic(directory, timestamp, outname, format="GTiff"):
    if format == "GTiff" and not outname.endswith(".tif"):
        outname += ".tif"
    files = finder(directory, [timestamp+"$"], regex=True)
    run(["gdalwarp", "-overwrite", "--config", "GDAL_CACHEMAX", 2000, "-wm", 6000, "-multi", "-of", format, files, outname])
