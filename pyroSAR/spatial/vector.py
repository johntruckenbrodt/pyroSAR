# -*- coding: utf-8 -*-
##############################################################
# OGR wrapper for convenient vector data handling and processing
# John Truckenbrodt 2015-2018
##############################################################

"""
This is intended as a vector meta information handler with options for reading and writing vector data in a convenient
manner by simplifying the numerous options provided by the OGR python binding
"""

import os

from osgeo import ogr, osr

from pyroSAR.spatial.auxil import crsConvert
from pyroSAR.ancillary import parse_literal

ogr.UseExceptions()
osr.UseExceptions()


class Vector(object):
    def __init__(self, filename=None, driver='ESRI Shapefile'):

        if driver not in ['ESRI Shapefile', 'Memory']:
            raise IOError('driver not supported')

        if filename is None:
            driver = 'Memory'
        else:
            self.filename = filename

        self.driver = ogr.GetDriverByName(driver)

        self.vector = self.driver.CreateDataSource('out') if driver == 'Memory' else self.driver.Open(filename)

        nlayers = self.vector.GetLayerCount()
        if nlayers > 1:
            raise IOError('multiple layers are currently not supported')
        elif nlayers == 1:
            self.init_layer()

    def __getitem__(self, expression):
        """
        subset the vector object by index or attribute.
        See `ogr.Layer.SetAttributeFilter <http://gdal.org/python/osgeo.ogr.Layer-class.html#SetAttributeFilter>`_
        for details on the expression syntax.
        Parameters
        ----------
        expression: int or str
            the key or expression to be used for subsetting
        Returns
        -------
        Vector
            a vector object matching the specified criteria
        """
        if not isinstance(expression, (int, str)):
            raise RuntimeError('expression must be of type int or str')
        expression = parse_literal(expression) if isinstance(expression, str) else expression
        if isinstance(expression, int):
            feat = self.getFeatureByIndex(expression)
        else:
            self.layer.SetAttributeFilter(expression)
            feat = self.getfeatures()
            feat = feat if len(feat) > 0 else None
            self.layer.SetAttributeFilter('')
        if feat is None:
            return None
        else:
            return feature2vector(feat, ref=self)

    def init_layer(self):
        self.layer = self.vector.GetLayer()
        self.__features = [None]*self.nfeatures

    def init_features(self):
        del self.__features
        self.__features = [None]*self.nfeatures

    @property
    def extent(self):
        return dict(zip(["xmin", "xmax", "ymin", "ymax"], self.layer.GetExtent()))

    @property
    def fieldDefs(self):
        return [self.layerdef.GetFieldDefn(x) for x in range(0, self.nfields)]

    @property
    def fieldnames(self):
        return [field.GetName() for field in self.fieldDefs]

    @property
    def geomType(self):
        return self.layerdef.GetGeomType()

    @property
    def layerdef(self):
        return self.layer.GetLayerDefn()

    @property
    def layername(self):
        return self.layer.GetName()

    @property
    def nlayers(self):
        return self.vector.GetLayerCount()

    @property
    def nfeatures(self):
        return len(self.layer)

    @property
    def nfields(self):
        return self.layerdef.GetFieldCount()

    def getProjection(self, type):
        """
        type can be either "epsg", "wkt", "proj4" or "osr"
        """
        return crsConvert(self.layer.GetSpatialRef(), type)

    @property
    def proj4(self):
        return self.srs.ExportToProj4().strip()

    @property
    def srs(self):
        return self.layer.GetSpatialRef()

    # todo Should return the wkt of the object, not of the projection
    @property
    def wkt(self):
        return self.srs.ExportToWkt()

    def addfield(self, name, type=ogr.OFTString, width=10):
        fieldDefn = ogr.FieldDefn(name, type)
        if type == ogr.OFTString:
            fieldDefn.SetWidth(width)
        self.layer.CreateField(fieldDefn)

    def addlayer(self, name, srs, geomType):
        self.vector.CreateLayer(name, srs, geomType)
        self.init_layer()

    def addfeature(self, geometry, fieldname, fieldvalue):

        if fieldname not in self.fieldnames:
            raise IOError("field does not exist")

        featureDefn = self.layerdef
        feature = ogr.Feature(featureDefn)
        feature.SetGeometry(geometry)
        feature.SetField(fieldname, fieldvalue)
        self.layer.CreateFeature(feature)
        feature.Destroy()
        self.init_features()

    def convert2wkt(self, set3D=True):
        """
        export the geometry of each feature as a wkt string
        """
        features = self.getfeatures()
        for feature in features:
            try:
                feature.geometry().Set3D(set3D)
            except AttributeError:
                dim = 3 if set3D else 2
                feature.geometry().SetCoordinateDimension(dim)

        return [feature.geometry().ExportToWkt() for feature in features]

    def getArea(self):
        return sum([x.GetGeometryRef().GetArea() for x in self.getfeatures()])

    def getFeatureByAttribute(self, fieldname, attribute):
        attr = attribute.strip() if isinstance(attribute, str) else attribute
        if fieldname not in self.fieldnames:
            raise KeyError('invalid field name')
        out = []
        self.layer.ResetReading()
        for feature in self.layer:
            field = feature.GetField(fieldname)
            field = field.strip() if isinstance(field, str) else field
            if field == attr:
                out.append(feature.Clone())
        self.layer.ResetReading()
        if len(out) == 0:
            return None
        elif len(out) == 1:
            return out[0]
        else:
            return out

    def getFeatureByIndex(self, index):
        feature = self.layer[index]
        if feature is None:
            feature = self.getfeatures()[index]
        return feature

    def getfeatures(self):
        self.layer.ResetReading()
        features = [x.Clone() for x in self.layer]
        self.layer.ResetReading()
        return features

    def load(self):
        self.layer.ResetReading()
        for i in range(self.nfeatures):
            if self.__features[i] is None:
                self.__features[i] = self.layer[i]

    def reproject(self, projection):

        srs_out = crsConvert(projection, "osr")

        # create the CoordinateTransformation
        coordTrans = osr.CoordinateTransformation(self.srs, srs_out)

        layername = self.layername
        geomType = self.geomType
        features = self.getfeatures()

        self.__init__()
        self.addlayer(layername, srs_out, geomType)

        for feature in features:
            geom = feature.GetGeometryRef()
            geom.Transform(coordTrans)
            newfeature = feature.Clone()
            newfeature.SetGeometry(geom)
            self.layer.CreateFeature(newfeature)
            newfeature.Destroy()
        self.init_features()

    def setCRS(self, crs):
        """
        directly reset the spatial reference system of the vector object
        Parameters
        ----------
        crs: int, str or osr.SpatialReference
            the input CRS

        Returns
        -------

        Example
        -------
        >>> site = Vector('shape.shp')
        >>> site.setCRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs ')

        """
        # try to convert the input crs to osr.SpatialReference
        srs_out = crsConvert(crs, 'osr')

        # save all relevant info from the existing vector object
        layername = self.layername
        geomType = self.geomType
        layer_definition = ogr.Feature(self.layer.GetLayerDefn())
        fields = [layer_definition.GetFieldDefnRef(x) for x in range(layer_definition.GetFieldCount())]
        features = self.getfeatures()

        # initialize a new vector object and create a layer
        self.__init__()
        self.addlayer(layername, srs_out, geomType)

        # add the fields to new layer
        self.layer.CreateFields(fields)

        # add the features to the newly created layer
        for feat in features:
            self.layer.CreateFeature(feat)
        self.init_features()

    def write(self, outfile, format='ESRI Shapefile', overwrite=True):
        (outfilepath, outfilename) = os.path.split(outfile)
        basename = os.path.splitext(outfilename)[0]

        driver = ogr.GetDriverByName(format)

        if os.path.exists(outfile):
            if overwrite:
                driver.DeleteDataSource(outfile)
            else:
                raise IOError("file already exists")

        outdataset = driver.CreateDataSource(outfile)
        outlayer = outdataset.CreateLayer(self.layername, geom_type=self.geomType)
        outlayerdef = outlayer.GetLayerDefn()

        for fieldDef in self.fieldDefs:
            outlayer.CreateField(fieldDef)

        for feature in self.layer:
            outFeature = ogr.Feature(outlayerdef)
            outFeature.SetGeometry(feature.GetGeometryRef())
            for j in range(0, self.nfields):
                outFeature.SetField(self.fieldnames[j], feature.GetField(j))
            # add the feature to the shapefile
            outlayer.CreateFeature(outFeature)
            outFeature.Destroy()

        srs_out = self.srs.Clone()
        srs_out.MorphToESRI()
        with open(os.path.join(outfilepath, basename+".prj"), "w") as prj:
            prj.write(srs_out.ExportToWkt())

        outdataset.Destroy()


def feature2vector(feature, ref, layername=None):
    """
    create a Vector object from ogr features
    Parameters
    ----------
    feature: ogr.Feature or list
        a single feature or a list of features
    ref: Vector
        a reference Vector object to retrieve geo information
    layername: str or None
        the name of the output layer; retrieved from ref if None

    Returns
    -------

    """
    features = feature if isinstance(feature, list) else [feature]
    layername = layername if layername is not None else ref.layername
    vec = Vector(driver='Memory')
    vec.addlayer(layername, ref.srs, ref.geomType)
    feat_def = features[0].GetDefnRef()
    fields = [feat_def.GetFieldDefn(x) for x in range(0, feat_def.GetFieldCount())]
    vec.layer.CreateFields(fields)
    for feat in features:
        vec.layer.CreateFeature(feat)
    vec.init_features()
    return vec


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

    bbox = Vector(driver='Memory')
    bbox.addlayer('bbox', srs, ogr.wkbPolygon)
    bbox.addfield('id', width=4)
    bbox.addfeature(geom, 'id', 1)
    geom.Destroy()
    if outname is None:
        return bbox
    else:
        bbox.write(outname, format, overwrite)


def centerdist(obj1, obj2):
    if not isinstance(obj1, Vector) or isinstance(obj2, Vector):
        raise IOError('both objects must be of type Vector')

    feature1 = obj1[0]
    geometry1 = feature1.GetGeometryRef()
    center1 = geometry1.Centroid()

    feature2 = obj2[0]
    geometry2 = feature2.GetGeometryRef()
    center2 = geometry2.Centroid()

    return center1.Distance(center2)


def intersect(obj1, obj2):
    if not isinstance(obj1, Vector) or not isinstance(obj2, Vector):
        raise RuntimeError('both objects must be of type Vector')

    if obj1.nfeatures > 1 or obj2.nfeatures > 1:
        raise RuntimeError('only objects with one feature are currently supported')

    obj1.reproject(obj2.srs)

    feature1 = obj1.getFeatureByIndex(0)
    geometry1 = feature1.GetGeometryRef()

    feature2 = obj2.getFeatureByIndex(0)
    geometry2 = feature2.GetGeometryRef()

    intersect = geometry2.Intersection(geometry1)

    return intersect if intersect.GetArea() > 0 else None
