from osgeo import ogr, osr

import raster
import vector
from .auxil import crsConvert

osr.UseExceptions()


def bbox(coordinates, crs, outname=None, format='ESRI Shapefile', overwrite=True):
    """
    create a bounding box vector object or shapefile from coordinates and coordinate reference system
    coordinates must be provided in a dictionary containing numerical variables with names 'xmin', 'xmax', 'ymin' and 'ymax'
    the coordinate reference system can be in either WKT, EPSG or PROJ4 format
    """
    srs = crsConvert(crs, 'osr')

    ring = ogr.Geometry(ogr.wkbLinearRing)

    ring.AddPoint(coordinates['xmin'], coordinates['ymin'])
    ring.AddPoint(coordinates['xmin'], coordinates['ymax'])
    ring.AddPoint(coordinates['xmax'], coordinates['ymax'])
    ring.AddPoint(coordinates['xmax'], coordinates['ymin'])
    ring.CloseRings()

    geom = ogr.Geometry(ogr.wkbPolygon)
    geom.AddGeometry(ring)

    geom.FlattenTo2D()

    bbox = vector.Vector(driver='Memory')
    bbox.addlayer('bbox', srs, ogr.wkbPolygon)
    bbox.addfield('id', width=4)
    bbox.addfeature(geom, 'id', 1)
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
        raise IOError('object must be of type Raster or Vector')
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
