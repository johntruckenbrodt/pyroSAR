##############################################################
# OGR wrapper for convenient vector data handling and processing
# John Truckenbrodt 2015
# last update 2015-12-09
##############################################################

"""
This is intended as a vector meta information handler with options for reading and writing vector data in a convenient manner by simplifying the numerous options provided
by the OGR python binding
"""

import os
from osgeo import ogr, osr
from ancillary import crsConvert

ogr.UseExceptions()


class Vector(object):
    #todo Define get_projection which returns the projection in a given format
    def __init__(self, filename=None, driver="ESRI Shapefile"):

        if driver not in ["ESRI Shapefile", "Memory"]:
            raise IOError("driver not supported")

        if filename is None:
            driver = "Memory"
        else:
            self.filename = filename

        self.driver = ogr.GetDriverByName(driver)

        self.vector = self.driver.CreateDataSource("out") if driver == "Memory" else self.driver.Open(filename)

        self.nlayers = self.vector.GetLayerCount()
        if self.nlayers > 1:
            raise IOError("multiple layers are currently not supported")
        elif self.nlayers == 1:
            self.init_layer()

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
    def nfeatures(self):
        return self.layer.GetFeatureCount()

    @property
    def nfields(self):
        return self.layerdef.GetFieldCount()

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

    def getArea(self):
        return sum([x.GetGeometryRef().GetArea() for x in self.getfeatures()])

    def getfeatures(self):
        features = [None]*self.nfeatures
        self.layer.SetNextByIndex(0)
        for i in range(self.nfeatures):
            if self.__features[i] is None:
                features[i] = self.layer.GetNextFeature()
            else:
                features[i] = self.__features[i].Clone()
                self.layer.SetNextByIndex(i+1)
        self.layer.SetNextByIndex(0)
        return features

    def load(self):
        self.layer.SetNextByIndex(0)
        for i in range(self.nfeatures):
            if self.__features[i] is None:
                self.__features[i] = self.layer.GetNextFeature()
            else:
                self.layer.SetNextByIndex(i+1)
        self.layer.SetNextByIndex(0)

    def reproject(self, projection):

        if isinstance(projection, osr.SpatialReference):
            srs_out = projection.Clone()
        else:
            srs_out = osr.SpatialReference()
            srs_out.ImportFromWkt(crsConvert(projection, "wkt"))

        # create the CoordinateTransformation
        coordTrans = osr.CoordinateTransformation(self.srs, srs_out)

        layername = self.layername
        geomType = self.geomType
        features = self.getfeatures()

        self.__init__()
        self.addlayer(layername, srs_out, geomType)

        for i in range(len(features)):
            geom = features[i].GetGeometryRef()
            geom.Transform(coordTrans)

            newfeature = features[i].Clone()
            newfeature.SetGeometry(geom)
            self.layer.CreateFeature(newfeature)
            newfeature.Destroy()
        self.init_features()

    def write(self, outfile, format="ESRI Shapefile", overwrite=True):
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

        inFeatures = self.getfeatures()

        for i in range(len(inFeatures)):
            outFeature = ogr.Feature(outlayerdef)
            outFeature.SetGeometry(inFeatures[i].GetGeometryRef())
            for j in range(0, self.nfields):
                outFeature.SetField(self.fieldnames[j], inFeatures[i].GetField(j))

            # add the feature to the shapefile
            outlayer.CreateFeature(outFeature)

            # destroy the features and get the next input feature
            outFeature.Destroy()
            inFeatures[i].Destroy()

        srs_out = self.srs.Clone()
        srs_out.MorphToESRI()
        with open(os.path.join(outfilepath, basename+".prj"), "w") as prj:
            prj.write(srs_out.ExportToWkt())

        outdataset.Destroy()
