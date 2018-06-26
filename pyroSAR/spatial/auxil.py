##############################################################
# Convenience functions for general spatial applications
# John Truckenbrodt, 2016-2018
##############################################################
import math
import os
from osgeo import osr, gdal

osr.UseExceptions()


def crsConvert(crsIn, crsOut):
    """
    convert between different types of spatial references

    Parameters
    ----------
    crsIn: int, str or osr.SpatialReference
        the input CRS
    crsOut: {'wkt', 'proj4', 'epsg', 'osr', 'opengis' or 'prettyWkt'}
        the output CRS type

    Returns
    -------
    int, str or osr.SpatialReference
        the output CRS

    Examples
    --------
    convert an integer EPSG code to PROJ4:

    >>> crsConvert(4326, 'proj4')
    '+proj=longlat +datum=WGS84 +no_defs '

    convert a PROJ4 string to an opengis URL:

    >>> crsConvert('+proj=longlat +datum=WGS84 +no_defs ', 'opengis')
    'http://www.opengis.net/def/crs/EPSG/0/4326'

    convert the opengis URL back to EPSG:

    >>> crsConvert('http://www.opengis.net/def/crs/EPSG/0/4326', 'epsg')
    4326

    """
    if isinstance(crsIn, osr.SpatialReference):
        srs = crsIn.Clone()
    else:
        srs = osr.SpatialReference()
        try:
            if 'opengis.net/def/crs/EPSG/0/' in str(crsIn):
                crsIn = int(os.path.basename(crsIn.strip('/')))
            srs.ImportFromEPSG(crsIn)
        except (TypeError, RuntimeError):
            try:
                srs.ImportFromWkt(crsIn)
            except (TypeError, RuntimeError):
                try:
                    srs.ImportFromProj4(crsIn)
                except (TypeError, RuntimeError):
                    raise TypeError('crsIn not recognized; must be of type WKT, PROJ4 or EPSG')
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
    compute the distance in meters between two points in latlon

    Parameters
    ----------
    lat1: int or float
        the latitude of point 1
    lon1: int or float
        the longitude of point 1
    lat2: int or float
        the latitude of point 2
    lon2: int or float
        the longitude of point 2

    Returns
    -------
    float
        the distance between point 1 and point2 in meters

    """
    radius = 6371000
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    a = math.sin((lat2-lat1)/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin((lon2-lon1)/2)**2
    c = 2 * math.asin(math.sqrt(a))
    return radius * c


def gdalwarp(src, dst, options):
    """
    a simple wrapper for gdal.Warp

    Parameters
    ----------
    src: str, ogr.DataSource or gdal.DataSource
        the input data set
    dst: str
        the output data set
    options: dict
        additional parameters passed to gdal.Warp;
        see http://gdal.org/python/osgeo.gdal-module.html#WarpOptions

    Returns
    -------

    """
    out = gdal.Warp(dst, src, options=gdal.WarpOptions(**options))
    out = None


def gdalbuildvrt(src, dst, options=None):
    """
    a simple wrapper for gdal.BuildVRT

    Parameters
    ----------
    src: str, ogr.DataSource or gdal.DataSource
        the input data set
    dst: str
        the output data set
    options: dict
        additional parameters passed to gdal.BuildVRT;
        see http://gdal.org/python/osgeo.gdal-module.html#BuildVRTOptions

    Returns
    -------

    """
    options = {} if options is None else options
    out = gdal.BuildVRT(dst, src, options=gdal.BuildVRTOptions(**options))
    out = None


def gdal_translate(src, dst, options):
    """
    a simple wrapper for gdal.Translate

    Parameters
    ----------
    src: str, ogr.DataSource or gdal.DataSource
        the input data set
    dst: str
        the output data set
    options: dict
        additional parameters passed to gdal.Translate;
        see http://gdal.org/python/osgeo.gdal-module.html#TranslateOptions

    Returns
    -------

    """
    out = gdal.Translate(dst, src, options=gdal.TranslateOptions(**options))
    out = None


def ogr2ogr(src, dst, options):
    """
    a simple wrapper for gdal.VectorTranslate aka ogr2ogr

    Parameters
    ----------
    src: str or ogr.DataSource
        the input data set
    dst: str
        the output data set
    options: dict
        additional parameters passed to gdal.VectorTranslate;
        see http://gdal.org/python/osgeo.gdal-module.html#VectorTranslateOptions

    Returns
    -------

    """
    out = gdal.VectorTranslate(dst, src, options=gdal.VectorTranslateOptions(**options))
    out = None


def gdal_rasterize(src, dst, options):
    """
    a simple wrapper for gdal.Rasterize

    Parameters
    ----------
    src: str or ogr.DataSource
        the input data set
    dst: str
        the output data set
    options: dict
        additional parameters passed to gdal.Rasterize;
        see http://gdal.org/python/osgeo.gdal-module.html#RasterizeOptions

    Returns
    -------

    """
    out = gdal.Rasterize(dst, src, options=gdal.RasterizeOptions(**options))
    out = None
