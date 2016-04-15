import raster
import vector
import spatial


class Continent(object):
    def __init__(self, name, acronym, center_lon, center_lat):
        self.name = name
        self.acronym = acronym
        self.center_lat = center_lat
        self.center_lon = center_lon
        self.intersect = None
        self.shape = vector.Vector("/homes4/geoinf/ve39vem/RADAR/Sentinel/projections/{}_t3.shp".format(self.acronym))


continents = [Continent("Africa", "AF", 21.5, 8.5),
              Continent("Asia", "AS", 94.0, 47.0),
              Continent("Europe", "EU", 24.0, 53.0),
              Continent("NorthAmerica", "NA", -97.5, 52.0),
              Continent("Oceania", "OC", 131.5, -19.5),
              Continent("SouthAmerica", "SA", -60.5, -14.0),
              Continent("Antarctica", "AN", 0.0, -90.0)]

# continent = "Europe"
# srs = 'PROJCS["World_Azimuthal_Equidistant", ' \
#       'GEOGCS["GCS_WGS_1984", DATUM["WGS_1984", SPHEROID["WGS_1984",6378137,298.257223563]], PRIMEM["Greenwich",0], UNIT["Degree",0.017453292519943295]], ' \
#       'PROJECTION["Azimuthal_Equidistant"], ' \
#       'PARAMETER["False_Easting",0], ' \
#       'PARAMETER["False_Northing",0], ' \
#       'PARAMETER["longitude_of_center",{0}], ' \
#       'PARAMETER["latitude_of_center",{1}], ' \
#       'UNIT["Meter",1], ' \
#       'AUTHORITY["EPSG","54032"]]'.format(center_lon[continent], center_lat[continent])
# print srs

test = raster.Raster("/homes4/geoinf/ve39vem/RADAR/Sentinel/tiling/N62E014/N62E014_20150504T163006.tif")


relations = [(region, spatial.intersect(test, region.shape), spatial.centerdist(test, region.shape)) for region in continents]

candidates = [entry for entry in relations if entry[1] is not None]

dists = [entry[2] for entry in candidates]
winner = candidates[dists.index(min(dists))][0]
