Handling of Orbit State Vector Files
------------------------------------
SAR products require additional orbit state vector (OSV) information to improve their spatial location accuracy.
This information is found in externally hosted files, which need to be downloaded separately and are then used by SAR
processing software to update the product's metadata. Currently, pyroSAR only supports handling of Sentinel-1 OSV files.

In SNAP, the corresponding processing node is called `Apply-Orbit-File`, which automatically downloads the OSV file and
updates the scene's metadata. The files are stored in SNAP's location for auxiliary data,
which per default is `$HOME/.snap/auxdata/Orbits`.

In GAMMA, on the other hand, the downloading has to be done manually after which the command `isp.S1_OPOD_vec` can be
used for updating the metadata. pyroSAR offers several approaches for automatically downloading these
files. The central tool for managing existing files and downloading new ones is the class :class:`pyroSAR.S1.OSV`, which
is used for all approaches.

**approach 1: direct download by time span**

In case a large number of scenes is to be processed and/or no internet access is available during processing, the files
can be downloaded by time span to a central directory. This is the most basic approach using the central class
:class:`~pyroSAR.S1.OSV` mentioned above, making use of its methods :meth:`~pyroSAR.S1.OSV.catch` and
:meth:`~pyroSAR.S1.OSV.retrieve`.

.. code-block:: python

    from pyroSAR.S1 import OSV

    osvdir = '/path/to/osvdir'

    with OSV(osvdir) as osv:
        files = osv.catch(sensor='S1A', osvtype='POE',
                          start='20170101T000000', stop='20180101T000000')
        osv.retrieve(files)

Two sub-directories `POEORB` and `RESORB` will be created in `osvdir` containing the downloaded files. `POEORB` will
contain the `Precise Orbit Ephemerides` files, which are the most accurate but are first available about two weeks after
the scene's acquisition. `RESORB` describes the `Restituted Orbit` files, which are less accurate but available
directly after acquisition. The files can be manually downloaded here: https://qc.sentinel1.eo.esa.int/.

**approach 2: manual download per scene**

The method :meth:`pyroSAR.drivers.SAFE.getOSV` can be used to directly retrieve the files relevant for the scene.
This method internally uses the methods described above with a time span limited to that of the scene acquisition.

.. code-block:: python

    from pyroSAR import identify
    scene = 'S1A_IW_GRDH_1SDV_20180101T170648_20180101T170713_019964_021FFD_DA78.zip'
    id = identify(scene)
    id.getOSV(outdir='/path/to/osvdir', osvType='POE')

**approach 3: direct download and scene metadata update**

The convenience function :func:`pyroSAR.gamma.correctOSV` internally makes use of approach 2 and additionally directly
executes the GAMMA command `isp.S1_OPOD_vec` for updating the scene's metadata with the information of the OSV file.
The scene has to be unpacked first.

.. code-block:: python

    from pyroSAR import identify
    from pyroSAR.gamma import correctOSV
    scene = 'S1A_IW_GRDH_1SDV_20180101T170648_20180101T170713_019964_021FFD_DA78.zip'
    id = identify(scene)
    id.unpack('tmpdir')
    correctOSV(id=id, osvdir='/path/to/osvdir', osvType='POE')

**approach 4: automatic download and use during processing**

The processing function :func:`pyroSAR.gamma.geocode` automatically downloads OSV files needed for processing and
updates the scene's metadata using function :func:`~pyroSAR.gamma.correctOSV`.
It is thus the most convenient way to handle these files and related processing steps.
The parameter `allow_RES_OSV` can be used to allow processing with `RES` files if no `POE` file is available yet.

.. code-block:: python

    from pyroSAR.gamma import geocode
    scene = 'S1A_IW_GRDH_1SDV_20180101T170648_20180101T170713_019964_021FFD_DA78.zip'
    geocode(scene=scene,
            dem='/path/to/demfile',
            tempdir='tmpdir',
            outdir='outdir',
            targetres=20,
            osvdir='/path/to/osvdir',
            allow_RES_OSV=False)

DEM Preparation
---------------

SAR processing requires a high resolution Digital Elevation Model for ortho-rectification and normalization of
terrain-specific imaging effects.

In SNAP, the DEM is usually auto-downloaded by the software itself and the user only specifies the DEM source to be
used, e.g. SRTM. pyroSAR's convenience function :func:`pyroSAR.snap.util.geocode` additional passes SNAP's option to use an
external DEM file via parameters `externalDEMFile`, `externalDEMNoDataValue` and `externalDEMApplyEGM`.

GAMMA does not provide ways to automatically download DEMs for processing and the user thus also needs to provide an
external DEM file in GAMMA's own format. However, several commands are available to prepare these DEMs including
conversion from EGM96 geoid heights to WGS84 ellipsoid heights.

pyroSAR offers several convenience functions to automatically prepare DEM mosaics from different openly available
sources to use them in either SNAP or GAMMA.

**download of DEM tiles**

The function :func:`pyroSAR.auxdata.dem_autoload` offers convenient download of tiles from four different sources
overlapping with user-defined geometries. Optionally, a buffer in degrees can be defined.
This function internally makes use of function :func:`spatialist.auxil.gdalbuildvrt`.

.. code-block:: python

    from pyroSAR.auxdata import dem_autoload
    from spatialist import Vector

    site = 'mysite.shp'
    vrt = 'mosaic.vrt'

    with Vector(site) as shp:
        vrt = dem_autoload(geometries=[vec],
                           demType='SRTM 1Sec HGT'
                           vrt = vrt,
                           buffer=0.1)

The tiles, which are delivered in compressed archives, are directly connected to a virtual mosaic using GDAL's VRT
format, making it easier to work with them by treating them as a single file.
For downloading TanDEM-X tiles (DEM type `TDX90m`), an account needs to be created and the user credentials be passed to
function :func:`~pyroSAR.auxdata.dem_autoload`. See the function's documentation for further details.

The files are stored in SNAP's location for auxiliary data, which per default is `$HOME/.snap/auxdata/dem`.
This function has proven beneficial in server environments where not each node has internet access and the tiles thus
need to be downloaded prior to processing on these nodes.

**DEM Mosaicing**

In a next step we create a mosaic GeoTiff cropped to the boundaries defined in the VRT using function
:func:`pyroSAR.auxdata.dem_create`.
The spatial reference system, WGS84 UTM 32N in this case, is defined by its EPSG code but also several other options
are available. Since for SAR processing we are interested in ellipsoid heights, we call the function with the according
parameter `geoid_convert` set to `True`.
This function makes use of :func:`spatialist.auxil.gdalwarp`.
Conversion of vertical reference systems, e.g. from geoid to ellipsoid, require GDAL version 2.2 at least.

.. code-block:: python

    from pyroSAR.auxdata import dem_create

    outname = 'mysite_srtm.tif'

    def dem_create(src=vrt, dst=outname,
                   t_srs=32632, tr=20,
                   resampling_method='bilinear',
                   geoid_convert=True, geoid='EGM96')

**GAMMA Import**

For convenience, pyroSAR's `gamma` submodule contains a function :func:`pyroSAR.gamma.dem.dem_autocreate`, which is a
combination of functions :func:`~pyroSAR.auxdata.dem_autoload` and :func:`~pyroSAR.auxdata.dem_create` and further
executes GAMMA commands for format conversion.
It offers the same parameters as these two functions and a user can additionally decide whether geoid-ellipsoid
conversion is done in GDAL or in GAMMA via parameter `geoid_mode`. The output is a file in GAMMA format, which can
directly be used for processing by e.g. function :func:`pyroSAR.gamma.geocode`.
