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
Users can either mosaic the result in GDAL's VRT format, return a file list with GDAL-readable paths
(e.g. pointing to a GeoTIFF inside a zip using the /vsizip/ directive), or just download products without any return.

When writing a VRT and ``crop=True`` (the default), the resulting mosaic is cropped to the extent of the (buffered) input geometry.
The ``crop`` argument does not have an effect when returning a file list.

.. note::

    VRTs do not support mosaics crossing the antimeridian. Use a file list in this case.

.. tab-set::
    :sync-group: dem_prep

    .. tab-item:: VRT
        :sync: vrt

        .. code-block:: python

            from pyroSAR.auxdata import dem_autoload
            from spatialist import bbox

            extent = {'xmin': 11.5, 'xmax': 12, 'ymin': 50.5, 'ymax': 51}
            vrt = 'mosaic.vrt'

            with bbox(extent, crs=4326) as vec:
                dem_autoload(
                    geometry=vec,
                    demType='Copernicus 30m Global DEM',
                    buffer=0.1,
                    vrt=vrt,
                )

    .. tab-item:: file list
        :sync: list

        .. code-block:: python

            from pyroSAR.auxdata import dem_autoload
            from spatialist import Vector

            extent = {'xmin': 11.5, 'xmax': 12, 'ymin': 50.5, 'ymax': 51}


            with bbox(extent, crs=4326) as vec:
                tiles = dem_autoload(
                    geometry=vec,
                    demType='Copernicus 30m Global DEM',
                    buffer=0.1,
                    return_fname=True
                )

For downloading tiles of some DEM types, an account needs to be created and the user credentials be passed to
function :func:`~pyroSAR.auxdata.dem_autoload`. See the function's documentation for further details.

The files are stored in SNAP's location for auxiliary data, which per default is `$HOME/.snap/auxdata/dem`.
This path can be modified using :attr:`pyroSAR.examine.ExamineSnap.auxdatapath`.

DEM Mosaicing
=============

In a next step we create a mosaic GeoTIFF using the function :func:`pyroSAR.auxdata.dem_create`.
The spatial reference system, WGS84 UTM 32N in this case, is defined by its EPSG code but also several other options
are available (see function :func:`spatialist.auxil.crsConvert` for options).
Since for SAR processing we are interested in WGS84 ellipsoid heights (and not geoid heights as for most DEMs),
the function defaults to `geoid_convert=True`. The correct geoid model is inferred automatically from the input DEM type.
This function makes use of :func:`spatialist.auxil.gdalwarp`.
Conversion of vertical reference systems, e.g. from geoid to ellipsoid, requires GDAL version >=2.2.

Since ``crop`` does not have an effect when returning a file list from
:func:`~pyroSAR.auxdata.dem_autoload`,
``geometry`` and ``buffer`` need also to be passed to
:func:`~pyroSAR.auxdata.dem_create` to achieve the same result as with the VRT.
If omitting them, the result will be the same as if creating a VRT in
:func:`~pyroSAR.auxdata.dem_autoload` with ``crop=False``,
i.e. creating a mosaic covering the extent of all input DEM tiles.

.. tab-set::
    :sync-group: dem_prep

    .. tab-item:: VRT
        :sync: vrt

        .. code-block:: python

            from pyroSAR.auxdata import dem_create


            outname = 'cop-dem.tif'


            dem_create(
                src=vrt,
                dst=outname,
                t_srs=32632,
                tr=(20, 20)


            )

    .. tab-item:: file list
        :sync: list

        .. code-block:: python

            from pyroSAR.auxdata import dem_create

            extent = {'xmin': 11.5, 'xmax': 12, 'ymin': 50.5, 'ymax': 51}
            outname = 'cop-dem.tif'

            with bbox(extent, crs=4326) as vec:
                dem_create(
                    geometry=vec,
                    buffer=0.1,
                    src=tiles,
                    dst=outname,
                    t_srs=32632,
                    tr=(20, 20)
                )

Next to the advantage of supporting the antimeridian case, the file list approach might also be preferred for the coverage of the output mosaic.
The mosaic created from a VRT covers exactly the EPSG:4326 extent of the input geometry plus buffer.
The mosaic created from the file list will cover the bounding box of the extent of the input geometry plus buffer projected to ``t_srs``:

.. tab-set::
    :sync-group: dem_prep

    .. tab-item:: VRT
        :sync: vrt

        .. figure:: figures/dem_from_vrt.png

    .. tab-item:: file list
        :sync: list

        .. figure:: figures/dem_from_list.png

GAMMA Import
============

For convenience, pyroSAR's :mod:`~pyroSAR.gamma` submodule contains a function :func:`pyroSAR.gamma.dem.dem_autocreate`, which is a
combination of functions :func:`~pyroSAR.auxdata.dem_autoload` and :func:`~pyroSAR.auxdata.dem_create` and further
executes GAMMA commands for format conversion.
It offers the same parameters as these two functions and a user can additionally decide whether geoid-ellipsoid
conversion is done in GDAL or in GAMMA via parameter `geoid_mode`. The output is a file in GAMMA format, which can
directly be used for processing by e.g. function :func:`pyroSAR.gamma.util.geocode`.
