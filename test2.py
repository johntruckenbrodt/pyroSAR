
import vector

test = vector.Vector("/homes4/geoinf/ve39vem/RADAR/Sentinel/projections/test4.shp")

newsrs = '+proj=aeqd +lat_0=53 +lon_0=24 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'

test.reproject(newsrs)

test.write("/homes4/geoinf/ve39vem/RADAR/Sentinel/projections/test5.shp")
