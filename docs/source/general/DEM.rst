###############
DEM Preparation
###############

SAR processing requires a high resolution Digital Elevation Model for ortho-rectification and normalization of
terrain-specific imaging effects.

In SNAP, the DEM is usually auto-downloaded by the software itself and the user only specifies the DEM source to be
used, e.g. SRTM. pyroSAR's convenience function :func:`pyroSAR.snap.util.geocode` can additionally pass SNAP's option to use an
external DEM file via parameters `externalDEMFile`, `externalDEMNoDataValue` and `externalDEMApplyEGM`.

GAMMA does not provide ways to automatically download DEMs for processing and the user thus also needs to provide an
external DEM file in GAMMA's own format. However, several commands are available to prepare these DEMs including
conversion from geoid heights to WGS84 ellipsoid heights.

pyroSAR offers several convenience functions to automatically prepare DEM mosaics from different
sources to use them in either SNAP or GAMMA.

Download of DEM Tiles
=====================

The function :func:`pyroSAR.auxdata.dem_autoload` offers convenient download of tiles from different sources
overlapping with user-defined geometries. Optionally, a buffer in degrees can be defined.
This function internally makes use of the function :func:`spatialist.auxil.gdalbuildvrt`.

.. code-block:: python

    from pyroSAR.auxdata import dem_autoload
    from spatialist import Vector

    site = 'mysite.shp'
    vrt = 'mosaic.vrt'

    with Vector(site) as vec:
        vrt = dem_autoload(geometries=[vec],
                           demType='SRTM 1Sec HGT',
                           vrt=vrt,
                           buffer=0.1)

The tiles, which are delivered in compressed archives, are directly connected to a virtual mosaic using GDAL's VRT
format, making it easier to work with them by treating them as a single file.
For downloading tiles of some DEM types, e.g. `TDX90m`, an account needs to be created and the user credentials be passed to
function :func:`~pyroSAR.auxdata.dem_autoload`. See the function's documentation for further details.

The files are stored in SNAP's location for auxiliary data, which per default is `$HOME/.snap/auxdata/dem`.
The function :func:`~pyroSAR.auxdata.dem_autoload` has proven beneficial in server environments where not each node has internet access and the tiles thus
need to be downloaded prior to processing on these nodes.

DEM Mosaicing
=============

In a next step we create a mosaic GeoTIFF cropped to the boundaries defined in the VRT using the function
:func:`pyroSAR.auxdata.dem_create`.
The spatial reference system, WGS84 UTM 32N in this case, is defined by its EPSG code but also several other options
are available. Since for SAR processing we are interested in ellipsoid heights, we call the function with the according
parameter `geoid_convert` set to `True`.
This function makes use of :func:`spatialist.auxil.gdalwarp`.
Conversion of vertical reference systems, e.g. from geoid to ellipsoid, requires GDAL version >=2.2.

.. code-block:: python

    from pyroSAR.auxdata import dem_create

    outname = 'mysite_srtm.tif'

    dem_create(src=vrt, dst=outname,
               t_srs=32632, tr=(20, 20),
               resampling_method='bilinear',
               geoid_convert=True, geoid='EGM96')

GAMMA Import
============

For convenience, pyroSAR's :mod:`~pyroSAR.gamma` submodule contains a function :func:`pyroSAR.gamma.dem.dem_autocreate`, which is a
combination of functions :func:`~pyroSAR.auxdata.dem_autoload` and :func:`~pyroSAR.auxdata.dem_create` and further
executes GAMMA commands for format conversion.
It offers the same parameters as these two functions and a user can additionally decide whether geoid-ellipsoid
conversion is done in GDAL or in GAMMA via parameter `geoid_mode`. The output is a file in GAMMA format, which can
directly be used for processing by e.g. function :func:`pyroSAR.gamma.geocode`.
