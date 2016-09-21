import math
import os

import raster
import vector
from osgeo import ogr, osr

osr.UseExceptions()


def bbox(coordinates, crs, outname=None, format="ESRI Shapefile", overwrite=True):
    """
    create a bounding box vector object or shapefile from coordinates and coordinate reference system
    coordinates must be provided in a dictionary containing numerical variables with names 'xmin', 'xmax', 'ymin' and 'ymax'
    the coordinate reference system can be in either WKT, EPSG or PROJ4 format
    """
    srs = crsConvert(crs, "osr")

    ring = ogr.Geometry(ogr.wkbLinearRing)

    ring.AddPoint(coordinates["xmin"], coordinates["ymin"])
    ring.AddPoint(coordinates["xmin"], coordinates["ymax"])
    ring.AddPoint(coordinates["xmax"], coordinates["ymax"])
    ring.AddPoint(coordinates["xmax"], coordinates["ymin"])
    ring.CloseRings()

    geom = ogr.Geometry(ogr.wkbPolygon)
    geom.AddGeometry(ring)

    geom.FlattenTo2D()

    bbox = vector.Vector(driver="Memory")
    bbox.addlayer("bbox", srs, ogr.wkbPolygon)
    bbox.addfield("id", width=4)
    bbox.addfeature(geom, "id", 1)
    geom.Destroy()
    if outname is None:
        return bbox
    else:
        bbox.write(outname, format, overwrite)


def init_vector(obj):
    if isinstance(obj, raster.Raster):
        shape = obj.bbox()
    elif isinstance(obj, vector.Vector):
        shape = obj
    else:
        raise IOError("object must be of type Raster or Vector")
    return shape


def centerdist(obj1, obj2):
    shape1 = init_vector(obj1)
    shape2 = init_vector(obj2)

    feature1 = shape1[0]
    geometry1 = feature1.GetGeometryRef()
    center1 = geometry1.Centroid()

    feature2 = shape2[0]
    geometry2 = feature2.GetGeometryRef()
    center2 = geometry2.Centroid()

    return center1.Distance(center2)


def intersect(obj1, obj2):
    shape1 = init_vector(obj1)
    shape2 = init_vector(obj2)

    shape1.reproject(shape2.srs)

    feature1 = shape1[0]
    geometry1 = feature1.GetGeometryRef()

    feature2 = shape2[0]
    geometry2 = feature2.GetGeometryRef()

    intersect = geometry2.Intersection(geometry1)

    return intersect if intersect.GetArea() > 0 else None


def crsConvert(crsIn, crsOut):
    """
    convert between epsg, wkt, proj4 and opengis spatial references
    crsText must be a osr.SpatialReference object, an opengis URL (e.g. "http://www.opengis.net/def/crs/EPSG/0/4326") or a string of type WKT, PROJ4 or EPSG
    crsOut must be either "wkt", "proj4", "epsg", "osr" or "opengis"
    if type "osr" is selected the function will return a spatial reference object of type osr.SpatialReference()
    """
    if isinstance(crsIn, osr.SpatialReference):
        srs = crsIn.Clone()
    else:
        srs = osr.SpatialReference()
        try:
            if "opengis.net/def/crs/EPSG/0/" in crsIn:
                crsIn = int(os.path.basename(crsIn.strip("/")))
            srs.ImportFromEPSG(crsIn)
        except (TypeError, RuntimeError):
            try:
                srs.ImportFromWkt(crsIn)
            except (TypeError, RuntimeError):
                try:
                    srs.ImportFromProj4(crsIn)
                except (TypeError, RuntimeError):
                    raise TypeError("crsText not recognized; must be of type WKT, PROJ4 or EPSG")
    if crsOut == "wkt":
        return srs.ExportToWkt()
    elif crsOut == "proj4":
        return srs.ExportToProj4()
    elif crsOut == "epsg":
        srs.AutoIdentifyEPSG()
        return int(srs.GetAuthorityCode(None))
    elif crsOut == "opengis":
        srs.AutoIdentifyEPSG()
        return "http://www.opengis.net/def/crs/EPSG/0/{}".format(srs.GetAuthorityCode(None))
    elif crsOut == "osr":
        return srs
    else:
        raise ValueError("crsOut not recognized; must be either wkt, proj4 or epsg")


def haversine(lat1, lon1, lat2, lon2):
    """
    compute distance in meters between two points in latlon
    """
    radius = 6371000
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    a = math.sin((lat2-lat1)/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin((lon2-lon1)/2)**2
    c = 2 * math.asin(math.sqrt(a))
    return radius * c