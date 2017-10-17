##############################################################
# Convenience functions for general spatial applications
# John Truckenbrodt, 2016-2017
##############################################################
import math
import os
from osgeo import osr, gdal

osr.UseExceptions()


def crsConvert(crsIn, crsOut):
    """
    convert between epsg, wkt, proj4 and opengis spatial references
    crsText must be a osr.SpatialReference object, an opengis URL (e.g. 'http://www.opengis.net/def/crs/EPSG/0/4326') or a string of type WKT, PROJ4 or EPSG
    crsOut must be either 'wkt', 'proj4', 'epsg', 'osr', 'opengis' or 'prettyWkt' (a wkt string formatted for readability)
    if type 'osr' is selected the function will return a spatial reference object of type osr.SpatialReference()
    """
    if isinstance(crsIn, osr.SpatialReference):
        srs = crsIn.Clone()
    else:
        srs = osr.SpatialReference()
        try:
            if 'opengis.net/def/crs/EPSG/0/' in crsIn:
                crsIn = int(os.path.basename(crsIn.strip('/')))
            srs.ImportFromEPSG(crsIn)
        except (TypeError, RuntimeError):
            try:
                srs.ImportFromWkt(crsIn)
            except (TypeError, RuntimeError):
                try:
                    srs.ImportFromProj4(crsIn)
                except (TypeError, RuntimeError):
                    raise TypeError('crsText not recognized; must be of type WKT, PROJ4 or EPSG')
    if crsOut == 'wkt':
        return srs.ExportToWkt()
    elif crsOut == 'prettyWkt':
        return srs.ExportToPrettyWkt()
    elif crsOut == 'proj4':
        return srs.ExportToProj4()
    elif crsOut == 'epsg':
        srs.AutoIdentifyEPSG()
        return int(srs.GetAuthorityCode(None))
    elif crsOut == 'opengis':
        srs.AutoIdentifyEPSG()
        return 'http://www.opengis.net/def/crs/EPSG/0/{}'.format(srs.GetAuthorityCode(None))
    elif crsOut == 'osr':
        return srs
    else:
        raise ValueError('crsOut not recognized; must be either wkt, proj4, opengis or epsg')


def haversine(lat1, lon1, lat2, lon2):
    """
    compute distance in meters between two points in latlon
    """
    radius = 6371000
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    a = math.sin((lat2-lat1)/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin((lon2-lon1)/2)**2
    c = 2 * math.asin(math.sqrt(a))
    return radius * c


# a simple wrapper for gdal.Warp
def gdalwarp(src, dst, options):
    out = gdal.Warp(dst, src, options=gdal.WarpOptions(**options))
    out = None


# a simple wrapper for gdal.BuildVRT
def gdalbuildvrt(src, dst, options):
    out = gdal.BuildVRT(dst, src, options=gdal.BuildVRTOptions(**options))
    out = None


# a simple wrapper for gdal.Translate
def gdal_translate(src, dst, options):
    out = gdal.Translate(dst, src, options=gdal.TranslateOptions(**options))
    out = None


# a simple wrapper for gdal.VectorTranslate aka ogr2ogr
def ogr2ogr(src, dst, options):
    out = gdal.VectorTranslate(dst, src, options=gdal.VectorTranslateOptions(**options))
    out = None


# a simple wrapper for gdal.Rasterize
def gdal_rasterize(src, dst, options):
    out = gdal.Rasterize(dst, src, options=gdal.RasterizeOptions(**options))
    out = None
