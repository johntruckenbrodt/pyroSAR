import raster
#
test = raster.Raster("/homes4/geoinf/ve39vem/RADAR/Sentinel/tiling/N62E014/N62E014_20150504T163006.tif")

# print test.geo
#
# print test.extract(test.geo["xmax"], test.geo["ymax"], 1)
#
# print test.extract(test.geo["xmax"], test.geo["ymin"], 1)
#
# print test.extract(test.geo["xmin"], test.geo["ymin"], 1)
#
# print test.extract(test.geo["xmin"], test.geo["ymax"], 1)

print test.geo

counter = 0
i = test.geo["xmin"]
while i <= test.geo["xmax"]:
    print i
    j = test.geo["ymin"]
    while j <= test.geo["ymax"]:
        x = test.extract(i, j, 1)
        j += test.res[1]
        counter += 1
    i += test.res[0]
print counter

