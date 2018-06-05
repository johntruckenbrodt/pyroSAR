##############################################################
# GDAL wrapper for convenient raster data handling and processing
# John Truckenbrodt 2015-2018
##############################################################


# todo: function to write data with the same metadata as a given file
# todo: documentation

from __future__ import division
import os
import re
import shutil
import subprocess as sp
from math import sqrt
from time import gmtime, strftime
import numpy as np

# from pyroSAR._common import (SpatialResults, nan_to_num, subset, test_import_type,
#                             test_string, test_data)
from pyroSAR import envi
from pyroSAR.spatial.auxil import (gdalwarp, gdalbuildvrt)
from pyroSAR.spatial.vector import (Vector, bbox, crsConvert, intersect)
from pyroSAR.ancillary import (dissolve, multicore)

from osgeo import (gdal, gdal_array, osr)
from osgeo.gdalconst import (GA_ReadOnly, GDT_Byte, GDT_Int16, GDT_UInt16,
                             GDT_Int32, GDT_UInt32, GDT_Float32, GDT_Float64)

os.environ['GDAL_PAM_PROXY_DIR'] = '/tmp'

gdal.UseExceptions()


class Raster(object):
    """
    This is intended as a raster meta information handler with options for reading and writing raster data in a
    convenient manner by simplifying the numerous options provided by the GDAL python binding.
    Several functions are provided along with this module to directly modify the raster object in memory or directly
    write a newly created file to disk (without modifying the rasterobject itself).
    Upon initializing a Raster object only metadata is loaded, the actual data can be, for example,
    loaded to memory by calling functions matrix or load.
    """

    # todo: init a Raster object from array data not only from a filename
    def __init__(self, filename):
        if os.path.isfile(filename):
            self.filename = filename if os.path.isabs(filename) else os.path.join(os.getcwd(), filename)
            self.raster = gdal.Open(filename, GA_ReadOnly)
        else:
            raise OSError('file does not exist')

        # a list to contain arrays
        self.__data = [None] * self.bands

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        self.raster = None

    @property
    def cols(self):
        return self.raster.RasterXSize

    @property
    def rows(self):
        return self.raster.RasterYSize

    @property
    def bands(self):
        return self.raster.RasterCount

    @property
    def dim(self):
        return [self.rows, self.cols, self.bands]

    @property
    def driver(self):
        return self.raster.GetDriver()

    @property
    def format(self):
        return self.driver.ShortName

    @property
    def dtype(self):
        return gdal.GetDataTypeName(self.raster.GetRasterBand(1).DataType)

    @property
    def projection(self):
        return self.raster.GetProjection()

    @property
    def srs(self):
        return osr.SpatialReference(wkt=self.projection)

    @property
    def proj4(self):
        return self.srs.ExportToProj4()

    @property
    def epsg(self):
        try:
            return crsConvert(self.srs, 'epsg')
        except RuntimeError:
            return None

    @property
    def geogcs(self):
        return self.srs.GetAttrValue('geogcs')

    @property
    def projcs(self):
        return self.srs.GetAttrValue('projcs') if self.srs.IsProjected() else None

    @property
    def geo(self):
        out = dict(zip(['xmin', 'xres', 'rotation_x', 'ymax', 'rotation_y', 'yres'],
                       self.raster.GetGeoTransform()))

        # note: yres is negative!
        out['xmax'] = out['xmin'] + out['xres'] * self.cols
        out['ymin'] = out['ymax'] + out['yres'] * self.rows
        return out

    @property
    def res(self):
        return [abs(float(self.geo['xres'])), abs(float(self.geo['yres']))]

    @property
    def nodata(self):
        return self.raster.GetRasterBand(1).GetNoDataValue()

    @nodata.setter
    def nodata(self, value):
        for i in range(1, self.bands + 1):
            self.raster.GetRasterBand(i).SetNoDataValue(value)

    @property
    def proj4args(self):
        args = [x.split('=') for x in re.split('[+ ]+', self.proj4) if len(x) > 0]
        return dict([(x[0], None) if len(x) == 1 else tuple(x) for x in args])

    @property
    def allstats(self):
        statcollect = []
        for x in self.layers():
            try:
                stats = x.ComputeStatistics(False)
            except RuntimeError:
                stats = None
            statcollect.append(stats)
        return statcollect

    def assign(self, array, dim='full'):
        """
        assign an array to an existing Raster object
        """
        self.__data = [array]
        if dim != 'full':
            shape = array.shape
            if len(shape) == 2:
                self.bands = 1
                self.rows, self.cols = shape
            else:
                self.bands, self.rows, self.cols = shape

            # print shape
            # print self.cols, self.rows
            # print self.raster.RasterXSize, self.raster.RasterYSize

            self.dim = [self.rows, self.cols, self.bands]
            self.geo['xmin'] += dim[0] * self.geo['xres']
            self.geo['ymax'] += dim[1] * self.geo['yres']
            self.geo['xmax'] = self.geo['xmin'] + self.geo['xres'] * self.cols
            self.geo['ymin'] = self.geo['ymax'] + self.geo['yres'] * self.rows
            self.raster.SetGeoTransform(
                [self.geo[x] for x in ['xmin', 'xres', 'rotation_x', 'ymax', 'rotation_y', 'yres']])

    def bbox(self, outname=None, format='ESRI Shapefile', overwrite=True):
        if outname is None:
            return bbox(self.geo, self.proj4)
        else:
            bbox(self.geo, self.proj4, outname=outname, format=format, overwrite=overwrite)

    def extract(self, px, py, radius=1, no_data=0):
        """
        extract weighted average of pixels intersecting with a defined radius to a point
        radius is a multiple of the pixel resolution

        Args:
            px:
            py:
            radius:
            no_data:

        Returns:

        """
        xres, yres = self.res

        hx = xres / 2.0
        hy = yres / 2.0

        xlim = float(xres * radius)
        ylim = float(yres * radius)

        # compute minimum x and y pixel coordinates
        xmin = int((px - self.geo['xmin'] - xlim) // xres)
        ymin = int((self.geo['ymax'] - py - xlim) // yres)

        xmin = xmin if xmin >= 0 else 0
        ymin = ymin if ymin >= 0 else 0

        # compute maximum x and y pixel coordinates
        xmax = int((px - self.geo['xmin'] + xlim) // xres) + 2
        ymax = int((self.geo['ymax'] - py + ylim) // yres) + 2

        xmax = xmax if xmax <= self.cols else self.cols
        ymax = ymax if ymax <= self.rows else self.rows

        # load array subset
        array = self.raster.GetRasterBand(1).ReadAsArray(xmin, ymin, xmax - xmin, ymax - ymin)

        sum = 0
        counter = 0
        weightsum = 0
        for x in range(xmin, xmax):
            for y in range(ymin, ymax):
                # check whether point is a valid image index
                val = array[y - ymin, x - xmin]
                if val != no_data:
                    # compute distances of pixel center coordinate to requested point

                    xc = x * xres + hx + self.geo['xmin']
                    yc = self.geo['ymax'] - y * yres + hy

                    dx = abs(xc - px)
                    dy = abs(yc - py)

                    # check whether point lies within ellipse: if ((dx ** 2) / xlim ** 2) + ((dy ** 2) / ylim ** 2) <= 1
                    weight = sqrt(dx ** 2 + dy ** 2)
                    sum += val * weight
                    weightsum += weight
                    counter += 1

        if sum > 0:
            return sum / weightsum
        else:
            if counter > 0:
                return 0
            else:
                return no_data

    def is_valid(self):
        """
        check image integrity
        tries to compute the checksum for each raster layer and returns False if this fails
        https://lists.osgeo.org/pipermail/gdal-dev/2013-November/037520.html

        :return: (logical) is the file valid?
        """
        for i in range(self.raster.RasterCount):
            try:
                checksum = self.raster.GetRasterBand(i + 1).Checksum()
            except RuntimeError:
                return False
        return True

    def layers(self):
        """
        get specific raster layer information objects

        Returns:

        """
        return [self.raster.GetRasterBand(band) for band in range(1, self.bands + 1)]

    def load(self, dim='full'):
        """
        load all raster data to arrays

        Args:
            dim:

        Returns:

        """
        dim = [0, 0, self.cols, self.rows] if dim == 'full' else dim
        for i in range(1, self.bands + 1):
            self.__data[i - 1] = self.matrix(i, dim)

    def matrix(self, band=1, dim='full'):
        """
        returns an array of a raster band

        Args:
            band:
            dim:

        Returns:

        """
        dim = [0, 0, self.cols, self.rows] if dim == 'full' else dim
        if self.__data[band - 1] is not None:
            return self.__data[band - 1][dim[1]:dim[3], dim[0]:dim[2]]
        else:
            return self.raster.GetRasterBand(band).ReadAsArray(*dim)

    # compute basic statistic measures from selected bands (provided by either single integer keys or a list of integers)
    # def getstat(self, statistic, bands='all'):
    #     statistics = {'min': 0, 'max': 1, 'mean': 2, 'sdev': 3}
    #     if statistic not in statistics:
    #         raise IOError('statistic not supported')
    #     if type(bands) == int:
    #         return self.allstats[bands-1][statistics[statistic]]
    #     elif bands == 'all':
    #         return [self.allstats[x-1][statistics[statistic]] for x in range(1, self.bands+1)]
    #     elif type(bands) == list:
    #         return [self.allstats[x-1][statistics[statistic]] for x in bands]

    # crop a raster object using another raster or extent object
    # if no name for an output file is provided, a list of pixel coordinates for cropping is returned
    # def crop(self, clipobject, outname=None):
    #     ext = Extent(self)
    #     inter = util.intersect(self, clipobject)
    #     if inter is None:
    #         raise IOError('no extent overlap')
    #     clip = [int(ceil((inter.left-ext.left)/self.res[0])), int(ceil((ext.top-inter.top)/self.res[1])),
    #             int(floor((inter.right-inter.left)/self.res[0])), int(floor((inter.top-inter.bottom)/self.res[1]))]
    #     if outname is not None:
    #         self.write(outname, dim=clip)
    #     else:
    #         return clip

    def reduce(self, outname=None, format='ENVI'):
        """
        remove all lines and columns containing only no data values

        Args:
            outname:
            format:

        Returns:

        """
        if self.bands != 1:
            raise ValueError('only single band images supported')

        stats = self.allstats[0]

        if stats[0] == stats[1]:
            raise ValueError('file does not contain valid pixels')

        # load raster layer into an array
        mat = self.matrix()

        mask1 = ~np.all(mat == self.nodata, axis=0)
        mask2 = ~np.all(mat == self.nodata, axis=1)
        mask1_l = mask1.tolist()
        mask2_l = mask2.tolist()

        left = mask1_l.index(True)
        cols = len(mask1_l) - mask1_l[::-1].index(True) - left
        top = mask2_l.index(True)
        rows = len(mask2_l) - mask2_l[::-1].index(True) - top

        mat = mat[mask2, :]
        mat = mat[:, mask1]

        if outname is None:
            self.assign(mat, dim=[left, top, cols, rows])
        else:
            self.write(outname, dim=[left, top, cols, rows], format=format)

    def rescale(self, function):
        """
        perform raster computations with custom functions and assign them to the existing raster object in memory

        Args:
            function:

        Returns:

        """
        if self.bands != 1:
            raise ValueError('only single band images are currently supported')

        # load array
        mat = self.matrix()

        # scale values
        scaled = function(mat)

        # round to nearest integer
        rounded = np.rint(scaled)

        # assign newly computed array to raster object
        self.assign(rounded)

    def write(self, outname, dtype='default', format='ENVI', dim='full', nodata='default', compress_tif=False):
        """
        write the raster object to a file
        if the data itself has been loaded to self.data (by function load), the in-memory data will be written to the file, otherwise the data is copied from the source file
        the parameter dim gives the opportunity to write a cropped version of the raster file; a dim-formatted list is, for example, returned by function crop

        Args:
            outname: str
                the name of the file to be written
            dtype: str
            format: str
            dim:
            nodata: int, float or str
                the nodata value to set to the output image;
                if set to 'default' the value is read from the currently opened file
            compress_tif: bool
                compress the created GeoTiff?

        Returns:

        """
        if os.path.isfile(outname):
            raise RuntimeError('target file already exists')

        if format == 'GTiff' and not re.search('\.tif[f]*$', outname):
            outname += '.tif'

        dtype = dtypes(self.dtype if dtype == 'default' else dtype)
        nodata = self.nodata if nodata == 'default' else nodata

        geo = list(self.raster.GetGeoTransform())

        options = []
        if format == 'GTiff' and compress_tif:
            options += ['COMPRESS=DEFLATE', 'PREDICTOR=2']

        if dim != 'full':
            geo[0] += dim[0] * geo[1]
            geo[3] += dim[1] * geo[5]

        dim = [0, 0, self.cols, self.rows] if dim == 'full' else dim
        driver = gdal.GetDriverByName(format)
        outDataset = driver.Create(outname, dim[2], dim[3], self.bands, dtype, options if len(options) > 0 else None)
        outDataset.SetMetadata(self.raster.GetMetadata())
        if self.geo is not None:
            outDataset.SetGeoTransform(geo)
        if self.projection is not None:
            outDataset.SetProjection(self.projection)
        for i in range(1, self.bands + 1):
            outband = outDataset.GetRasterBand(i)
            outband.SetNoDataValue(nodata)
            if self.__data[i - 1] is None:
                mat = self.raster.GetRasterBand(i).ReadAsArray(*dim)
            else:
                mat = self.__data[i - 1]
            outband.WriteArray(mat)
            outband.FlushCache()
        if format == 'GTiff':
            outDataset.SetMetadataItem('TIFFTAG_DATETIME', strftime('%Y:%m:%d %H:%M:%S', gmtime()))
        outDataset = None

        # write a png image of three raster bands (provided in a list of 1-based integers); percent controls the size ratio of input and output
        # def png(self, bands, outname, percent=10):
        #     if len(bands) != 3 or max(bands) not in range(1, self.bands+1) or min(bands) not in range(1, self.bands+1):
        #         print 'band indices out of range'
        #         return
        #     if not outname.endswith('.png'):
        #         outname += '.png'
        #     exp_bands = ' '.join(['-b '+str(x) for x in bands]).split()
        #     exp_scale = [['-scale', self.getstat('min', x), self.getstat('max', x), 0, 255] for x in bands]
        #     exp_size = ['-outsize', str(percent)+'%', str(percent)+'%']
        #     cmd = dissolve([['gdal_translate', '-q', '-of', 'PNG', '-ot', 'Byte'], exp_size, exp_bands, exp_scale, self.filename, outname])
        #     sp.check_call([str(x) for x in cmd])


def reproject(rasterobject, reference, outname, resampling='bilinear', format='ENVI'):
    """
    reproject a raster file
    """
    rasterobject = rasterobject if type(rasterobject) == Raster else Raster(rasterobject)
    projection = reference.projection if type(reference) == Raster else reference
    sp.check_call(['gdalwarp', '-overwrite', '-q', '-r', resampling, '-of', format,
                   '-tr', str(rasterobject.res[0]), str(rasterobject.res[1]),
                   '-srcnodata', str(rasterobject.nodata), '-dstnodata', str(rasterobject.nodata),
                   '-t_srs', projection, rasterobject.filename, outname])


# todo improve speed until aborting when all target files already exist
def stack(srcfiles, dstfile, resampling, targetres, srcnodata, dstnodata, shapefile=None, layernames=None, sortfun=None,
          separate=False, overwrite=False, compress=True, cores=4):
    """
    function for mosaicking, resampling and stacking of multiple raster files into a 3D data cube

    Parameters
    ----------
    srcfiles: list
        a list of file names or a list of lists; each sub-list is treated as an order to mosaic its containing files
    dstfile: str
        the destination file (if sesparate) or a directory
    resampling: {near, bilinear, cubic, cubicspline, lanczos, average, mode, max, min, med, Q1, Q3}
        the resampling method; see documentation of gdalwarp
    targetres: tuple
        a list with two entries for x and y spatial resolution
    srcnodata: int or float
        the nodata value of the source files
    dstnodata: int or float
        the nodata value of the destination file(s)
    shapefile: str or spatial.vector.Vector
        a shapefile for defining the area of the destination files
    layernames: list
        the names of the output layers; if None, the basenames of the input files are used
    sortfun: function
        a function for sorting the input files; this is needed for defining the mosaicking order
    separate: bool
        should the files be written to a single raster block or separate files? If separate, each tile is written to geotiff.
    overwrite: bool
        overwrite the file if it already exists?
    compress: bool
        compress the geotiff files?
    cores: int
        the number of CPU threads to use; this is only relevant if separate = True

    Returns
    -------
    """
    if len(dissolve(srcfiles)) == 0:
        raise IOError('no input files provided to function raster.stack')

    if layernames is not None:
        if len(layernames) != len(srcfiles):
            raise IOError('mismatch between number of source file groups and layernames')

    if not isinstance(targetres, (list, tuple)) or len(targetres) != 2:
        raise RuntimeError('targetres must be a list or tuple with two entries for x and y resolution')

    if len(srcfiles) == 1 and not isinstance(srcfiles[0], list):
        raise IOError('only one file specified; nothing to be done')

    if resampling not in ['near', 'bilinear', 'cubic', 'cubicspline', 'lanczos', 'average', 'mode', 'max', 'min', 'med',
                          'Q1', 'Q3']:
        raise IOError('resampling method not supported')

    projections = list()
    for x in dissolve(srcfiles):
        try:
            projection = Raster(x).projection
        except OSError as e:
            print('cannot read file: {}'.format(x))
            raise e
        projections.append(projection)

    projections = list(set(projections))
    if len(projections) > 1:
        raise IOError('raster projection mismatch')
    elif len(projections) == 0:
        raise RuntimeError('could not retrieve the projection from any of the {} input images'.format(len(srcfiles)))
    else:
        srs = projections[0]

    # read shapefile bounding coordinates and reduce list of rasters to those overlapping with the shapefile
    if shapefile is not None:
        shp = shapefile if isinstance(shapefile, Vector) else Vector(shapefile)
        shp.reproject(srs)
        ext = shp.extent
        arg_ext = (ext['xmin'], ext['ymin'], ext['xmax'], ext['ymax'])
        for i in range(len(srcfiles)):
            group = sorted(srcfiles[i], key=sortfun) if isinstance(srcfiles[i], list) else [srcfiles[i]]
            group = [x for x in group if intersect(shp, Raster(x).bbox())]
            if len(group) > 1:
                srcfiles[i] = group
            elif len(group) == 1:
                srcfiles[i] = group[0]
            else:
                srcfiles[i] = None
        srcfiles = filter(None, srcfiles)
    else:
        arg_ext = None

    # create temporary directory for writing intermediate files
    dst_base = os.path.splitext(dstfile)[0]
    tmpdir = dst_base + '__tmp'
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)

    options_warp = {'options': ['-q'],
                    'format': 'GTiff' if separate else 'ENVI',
                    'outputBounds': arg_ext, 'multithread': True,
                    'srcNodata': srcnodata, 'dstNodata': dstnodata,
                    'xRes': targetres[0], 'yRes': targetres[1],
                    'resampleAlg': resampling}

    if overwrite:
        options_warp['options'] += ['-overwrite']

    if separate and compress:
        options_warp['options'] += ['-co', 'COMPRESS=DEFLATE', '-co', 'PREDICTOR=2']

    options_buildvrt = {'outputBounds': arg_ext, 'srcNodata': srcnodata}

    # create VRT files for mosaicing
    for i in range(len(srcfiles)):
        base = srcfiles[i][0] if isinstance(srcfiles[i], list) else srcfiles[i]
        vrt = os.path.join(tmpdir, os.path.splitext(os.path.basename(base))[0] + '.vrt')
        gdalbuildvrt(srcfiles[i], vrt, options_buildvrt)
        srcfiles[i] = vrt

    # if no specific layernames are defined and sortfun is not set to None,
    # sort files by custom function or, by default, the basename of the raster/VRT file
    if layernames is None and sortfun is not None:
        srcfiles = sorted(srcfiles, key=sortfun if sortfun else os.path.basename)

    bandnames = [os.path.splitext(os.path.basename(x))[0] for x in srcfiles] if layernames is None else layernames

    if separate or len(srcfiles) == 1:
        if not os.path.isdir(dstfile):
            os.makedirs(dstfile)
        dstfiles = [os.path.join(dstfile, x) + '.tif' for x in bandnames]
        if overwrite:
            files = [x for x in zip(srcfiles, dstfiles)]
        else:
            files = [x for x in zip(srcfiles, dstfiles) if not os.path.isfile(x[1])]
            if len(files) == 0:
                print('all target tiff files already exist, nothing to be done')
                shutil.rmtree(tmpdir)
                return
        srcfiles, dstfiles = map(list, zip(*files))

        multicore(gdalwarp, cores=cores, multiargs={'src': srcfiles, 'dst': dstfiles}, options=options_warp)
    else:
        # create VRT for stacking
        vrt = os.path.join(tmpdir, os.path.basename(dst_base) + '.vrt')
        options_buildvrt['options'] = ['-separate']
        gdalbuildvrt(srcfiles, vrt, options_buildvrt)

        # warp files
        gdalwarp(vrt, dstfile, options_warp)

        # edit ENVI HDR files to contain specific layer names
        par = envi.HDRobject(dstfile + '.hdr')
        par.band_names = bandnames
        envi.hdr(par)

    # remove temporary directory and files
    shutil.rmtree(tmpdir)


def rasterize(vectorobject, outname, reference, burn_values=1, expressions=None, nodata=0):
    """
    rasterize a vector object

    Parameters
    ----------
    vectorobject: Vector
        the vector object to be rasterized
    outname: str
        the name of the GeoTiff output file
    reference: Raster
        a reference Raster object to retrieve geo information and extent from
    burn_values: int, or list
        the values to be written to the raster file
    expressions: list
        SQL expressions to filter the vector object by attributes
    nodata: int
        the nodata value of the target raster file

    Returns
    -------

    Example
    -------
    >>> from pyroSAR.spatial import Vector, Raster, rasterize
    >>> vec = Vector('source.shp')
    >>> ref = Raster('reference.tif')
    >>> outname = 'target.tif'
    >>> expressions = ['ATTRIBUTE=1', 'ATTRIBUTE=2']
    >>> burn_values = [1, 2]
    >>> rasterize(vec, outname, reference, burn_values, expressions)
    """
    if expressions is None:
        expressions = ['']
    if isinstance(burn_values, (int, float)):
        burn_values = [burn_values]
    if len(expressions) != len(burn_values):
        raise RuntimeError('expressions and burn_values of different length')
    target_ds = gdal.GetDriverByName('GTiff').Create(outname, reference.cols, reference.rows, 1, gdal.GDT_Byte)
    target_ds.SetGeoTransform(reference.raster.GetGeoTransform())
    target_ds.SetProjection(reference.raster.GetProjection())
    band = target_ds.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    band.FlushCache()
    for expression, value in zip(expressions, burn_values):
        vectorobject.layer.SetAttributeFilter(expression)
        gdal.RasterizeLayer(target_ds, [1], vectorobject.layer, burn_values=[value])
    vectorobject.layer.SetAttributeFilter('')
    target_ds = None


def dtypes(typestring):
    """
    translate raster data type descriptions to GDAl data type codes

    Args:
        typestring: str
            the data type string to be converted

    Returns:
    int
        the GDAL data type code
    """
    # create dictionary with GDAL style descriptions
    dictionary = {'Byte': GDT_Byte, 'Int16': GDT_Int16, 'UInt16': GDT_UInt16, 'Int32': GDT_Int32,
                  'UInt32': GDT_UInt32, 'Float32': GDT_Float32, 'Float64': GDT_Float64}

    # add numpy style descriptions
    dictionary.update(typemap())

    if typestring not in dictionary.keys():
        raise ValueError("unknown data type; use one of the following: ['{}']".format("', '".join(dictionary.keys())))

    return dictionary[typestring]


def typemap():
    """
    create a dictionary for mapping numpy data types to GDAL data type codes

    Returns
    -------
    dict
        the type map
    """
    TYPEMAP = {}
    for name in dir(np):
        obj = getattr(np, name)
        if hasattr(obj, 'dtype'):
            try:
                npn = obj(0)
                code = gdal_array.NumericTypeCodeToGDALTypeCode(npn.dtype.type)
                if code:
                    TYPEMAP[npn.dtype.name] = code
            except (TypeError, ValueError, AttributeError, OSError):
                pass
    return TYPEMAP
