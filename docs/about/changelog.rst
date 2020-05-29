#########
Changelog
#########

0.6 / 2018-11-20
================

SAR metadata
------------
- new standardized  metadata fields `orbitNumber_abs`, `orbitNumber_rel`, `cycleNumber` and `frameNumber` for all SAR
  formats
- customization of output file names with additional metadata fields (e.g. orbit numbers)

software configuration
----------------------
- pyroSAR configuration file handling: the paths to the SNAP and Gamma installation as well as relevant metadata
  directories are now registered in a configuration file `config.ini`, which is stored in a directory `.pyrosar` in the
  user home directory
- improved SNAP installation verification: pyroSAR now performs a deeper check of the SNAP installation to make sure
  it is not mistaken with e.g. the Ubuntu package manager snap; relevant installation executables and directories are
  stored in the configuration file

general functionality
---------------------
- deeper integration of package `spatialist <https://github.com/johntruckenbrodt/spatialist>`_: all the spatial file
  handling functionality that was part of pyroSAR is now part of package spatialist; now all the functionality is imported
  from spatialist and removed from pyroSAR
- improved search for datasets processed by pyroSAR: new helper functions exist, which make it easier to search for
  datasets by metadata fields, which are internally searched for in the respective file names
- introduced gamma function parser: these new tools search for a GAMMA_HOME environment variable and, if found, parse
  Python functions from the docstring of respective command line tools; for this, new Python scripts are created, which
  are stored alongside the configuration file in the user home directory; this way users can easily use Python functions
  with named parameters instead of the positional arguments of the Gamma command line tools
- improved documentation

Open Data Cube Export
---------------------
functionality to export processed datasets directly to an Open Data Cube:
it is now possible to create Open Data Cube product YML files as well as YML files for data indexing and ingestion
into this product; pyroSAR also internally checks for compatibility of a particular dataset with the target product;
this way, the resulting files can easily be passed to the Open Data Cube command line tools
several bug fixes

SNAP API
--------
improved SNAP processing workflow node linking: it is now possible to add a node also before an existing one, instead
of just after it

Python package integrity
------------------------
- add trove classifiers for supported operating systems and MIT license for easier online search
- exchange http with https for all URLs that support it

0.7 / 2019-01-03
================

several changes to the functioning of the Gamma command API

GAMMA API
---------

processing
++++++++++
- :func:`pyroSAR.gamma.geocode`:

  * optionally write all Gamma commands to shellscript
  * newly introduced choice of normalization method
  * changed normalization default approach

- :func:`pyroSAR.gamma.process`:

  * new parameter `logfile` to specify a logfile instead of just a directory with automated file naming
  * new parameter `shellscript` to write the executed command to a shell script protocol

command parser
++++++++++++++
- add parameters `outdir` and `shellscript` to parsed functions
- extensive improvement to accurately parse more commands
- add parameter `inlist` to some commands, which require interactive input via `stdin`

general
+++++++
- several bug fixes
- extended documentation
- make use of parsed command functions internally
- enable passing `logpath`, `outdir` and `shellscript` to all parsed functions via additional parameters for other
  convenience functions

0.8 / 2019-02-11
================

Auxiliary Data Handling
-----------------------

- new module auxdata with function :func:`pyroSAR.auxdata.dem_autoload` to automatically download tiles of
  different DEM types overlapping with given geometries
- class :class:`pyroSAR.S1.OSV`: reduced search time for new RES orbit state vector files;
  included more meaningful status messages

GAMMA API
---------

- new function :func:`pyroSAR.gamma.srtm.dem_autocreate` to automatically create DEMs in Gamma format from the output
  of function :func:`pyroSAR.auxdata.dem_autoload`
- improved writing of ENVI HDR files from class :class:`pyroSAR.gamma.ISPPar`
- class :class:`pyroSAR.gamma.UTM`: improved to work with newer Gamma versions
- function :func:`pyroSAR.gamma.geocode`:

  + improved documentation
  + clarified code for better readability
  + more consistent naming scheme for all temporarily written files
  + export temporarily written files (e.g. local incidence angle) via new parameter `export_extra`
  + additional parametrization tests to ensure best processing result
  + changed default of parameter `func_interp` to 2 to work best with default of parameter `normalization_method`
    (see documentation of Gamma command pixel_area)

SNAP API
--------

- function :func:`pyroSAR.snap.util.geocode`:

  + export temporarily written files (e.g. local incidence angle) via new parameter `export_extra`

0.9 / 2019-06-15
================

Drivers
-------

- :class:`pyroSAR.drivers.SAFE`: read heading angle, incident angle and image geometry (e.g. Ground Range) from metadata
- :class:`pyroSAR.drivers.Archive`: improved cross-compatibility with Python2 and Python3


SNAP API
--------

- function :func:`pyroSAR.snap.util.geocode`:

  + option to export `DEM` via parameter `export_extra`
  + added Sentinel-1 `ThermalNoiseRemoval` node via new parameter `removeS1ThermalNoise`
  + added `Multilook` node which is executed to approximate the target resolution if necessary
    (currently only for Sentinel-1 since metadata entries `incidence` and `image_geometry` are required)
  + new parameter `groupsize` to split workflows into several groups, which are executed separately with
    intermediate products written to disk. This increases processing speed
  + simplified internal node parametrization for easier use in future functions
  + fail if no POE orbit state vector file is found
  + `Terrain-Flattening`:

    * added additional parameters `additionalOverlap` and `oversamplingMultiple`
    * use bilinear instead of bicubic interpolation
  + `Remove-GRD-Border-Noise`: decrease `borderLimit` from 1000 to 500 (SNAP default)
  + new parameter `gpt_exceptions` to execute workflows containing specific nodes with different GPT versions than
    the default one
  + automatically remove node parameters on GPT fail and re-run the modified workflow; this is relevant if a node is
    executed in an older GPT version (e.g. via parameter `gpt_exceptions`), which does not accept parameters which were
    introduced in later GPT versions (e.g. those described above for node `Terrain-Flattening`)
  + disable/enable terrain flattening via new parameter `terrainFlattening`
  + optionally return workflow filename with new parameter `returnWF`
  + execute custom pyroSAR S1 GRD border noise removal (see :func:`pyroSAR.S1.removeGRDBorderNoise`)
  + new parameters `demResamplingMethod` and `imgResamplingMethod`

GAMMA API
---------

- SRTM Tools renamed to DEM Tools

  + function :func:`pyroSAR.gamma.dem.dem_autocreate`:

    * define arbitrary output CRS and resolution via new parameters `t_srs` and `tr`
    * optionally perform geoid to ellipsoid conversion in either GDAL or GAMMA via new parameter `geoid_mode`

- function :func:`pyroSAR.gamma.geocode`:

  + removed multiplication of backscatter with cosine of incident angle via command `lin_comb`
  + fixed bug in writing correct nodata values to ancillary products defined via parameter `export_extra`
  + changed default of parameter `func_geoback` from 2 to 1 (GAMMA default)

- function :func:`pyroSAR.gamma.correctOSV`:

  + fixed bug in using the first OSV file in a directory for correcting an image, which resulted in S1B files being
    corrected with S1A OSV files. This occasionally resulted in errors of no DEM overlap while processing S1B scenes

- fixed bug in treating GAMMA image pixel coordinates as top left instead of pixel center. This is relevant for writing
  ENVI HDR files for GAMMA images via function :func:`pyroSAR.gamma.par2hdr` resulting in the image to be shifted
  by 1/2 pixel to Southeast

Command Parser
++++++++++++++
- compatibility with GAMMA version released in November 2018
- delete parsed modules if environment variable `GAMMA_HOME` was reset causing them to be re-parsed with the new version
  on module import

general functionality
---------------------

- new function :func:`pyroSAR.ancillary.multilook_factors` to compute factors depending on image geometry and target resolution
- :func:`pyroSAR.S1.removeGRDBorderNoise`: reached Python3 compatibility

Auxiliary Data Handling
-----------------------

- new function :func:`pyroSAR.auxdata.dem_create` for convenient creation of DEM mosaics as downloaded by
  :func:`pyroSAR.auxdata.dem_autoload`

- function :func:`pyroSAR.auxdata.dem_autoload`: download 1 degree tiles instead of 5 degree tiles

- class :class:`pyroSAR.S1.OSV`:

  + download files specific to the Sentinel-1 sensor (S1A/S1B) instead of all matching the acquisition time
  + improved time span search, which occasionally resulted in missing OSV files

0.9.1 / 2019-07-05
==================

Auxiliary Data Handling
-----------------------

- function :func:`pyroSAR.auxdata.dem_create`: new parameter `resampling_method`

GAMMA API
---------

- function :func:`pyroSAR.gamma.dem.dem_autocreate`: new parameter `resampling_method`

SNAP API
--------

- function :func:`pyroSAR.snap.util.geocode`: fixed typo of parameter `removeS1BorderNoise`

0.10 / 2019-12-06
=================

Drivers
-------

- method :meth:`~pyroSAR.drivers.ID.bbox`: choose the output vector file format via new parameter `driver` or by
  using one of spatialist's supported file name extensions (see :meth:`spatialist.vector.Vector.write`)

- :class:`pyroSAR.drivers.SAFE`

  + new method :meth:`~pyroSAR.drivers.SAFE.quicklook` for writing KMZ quicklooks
  + method :meth:`~pyroSAR.drivers.SAFE.getOSV`: renamed parameter `outdir` to `osvdir`

- :class:`pyroSAR.drivers.Archive`: remove scenes from the database if they cannot be found at their file location.
  This is performed at each initialization of an `Archive` object.

GAMMA API
---------

- new parameter `basename_extensions` for adding extra metadata fields to output image names; affects:

  + :func:`pyroSAR.gamma.convert2gamma`
  + :func:`pyroSAR.gamma.geocode`

- :func:`pyroSAR.gamma.correctOSV`: make use of OSV files in SNAP's auxdata structure
- :func:`pyroSAR.gamma.geocode`: made border nose removal optional with new parameter `removeS1BorderNoise`

SNAP API
--------
- workflow parsing

  + improved output XML for better display in SNAP GUI
  + support for nodes with multiple input scenes, e.g. `SliceAssembly`

- SAR processor (function :func:`~pyroSAR.snap.auxil.gpt`)

  + write Sentinel-1 manifest.safe with processing results
  + two methods for border noise removal: `ESA` and `pyroSAR` via new parameter `removeS1BorderNoiseMethod`

- function :func:`pyroSAR.snap.util.geocode`

  + optional speckle filtering with new parameter `speckleFilter`
  + choose the output backscatter reference area (`beta0`/`gamma0`/`sigma0`) with new parameter `refarea`
  + default of parameter `groupsize` changed to 1
  + internally download S1 OSV files
  + internally download SNAP's `EGM96` geoid to `WGS84` ellipsoid DEM conversion lookup table via new function
    :func:`pyroSAR.snap.auxil.get_egm96_lookup`
  + support for multi-scene `SliceAssembly`; can be invoke by passing a list of scenes to parameter `infile`
  + new parameter `removeS1BorderNoiseMethod`
  + new parameter `gpt_args` to pass additional arguments to the GPT call

Datacube Tools
--------------

- :meth:`pyroSAR.datacube_util.Product.export_ingestion_yml`: new parameter `chunking`

Auxiliary Data Handling
-----------------------

- OSV download functionality (class :class:`pyroSAR.S1.OSV`)

  + made definition of OSV download directory optional; default is SNAP's auxdata directory
  + organization of downloaded files into SNAP's auxdata structure:

    * compression to zip
    * sort files into subdirs for sensor, year, month

  + removed method :meth:`~pyroSAR.S1.OSV.update`

Ancillary Tools
---------------
- :func:`pyroSAR.ancillary.parse_datasetname`

  + support for datasets in NetCDF format
  + enable parsing of ancillary products like local incidence angle (\*inc_geo.tif)

- :func:`pyroSAR.ancillary.find_datasets`:  new parameters `start` and `stop` for time filtering

general
-------
- bug fixes and documentation improvements

0.10.1 / 2019-12-12
===================

GAMMA API
---------

- :ref:`Command API <gamma-command-api>` compatibility with GAMMA version 20191203

0.11 / 2020-05-29
=================

Drivers
-------

- :class:`pyroSAR.drivers.Archive`: completely restructured to use the `SQLAlchemy <https://www.sqlalchemy.org/>`_
  Object Relational Mapper (ORM). This makes it possible to switch between SQLite+Spatialite and PostgreSQL+PostGIS
  database backends.

- :meth:`pyroSAR.drivers.SAFE.getOSV`: new argument `returnMatch` to also return the name of an OSV file instead of just
  downloading it.

SNAP API
--------

- arbitrary nodes can now be parsed. Before, only a small selection of nodes (those used by function
  :func:`~pyroSAR.snap.util.geocode`) were available. Now, any node and its default parametrization can be parsed to XML
  from the GPT documentation by internally calling e.g.:

    ::

        gpt Terrain-Flattening -h

  The parsed XML representation is saved for faster future reuse. See function :func:`~pyroSAR.snap.auxil.parse_node`
  for details. In all cases the standard SNAP file suffix is used for output products, e.g. `_TF` for
  `Terrain-Flattening`.

- multi-source nodes like `SliceAssembly` now take any number of sources, not just two.
  See class :class:`~pyroSAR.snap.auxil.Node`.

- function :func:`pyroSAR.snap.util.geocode`:

  + new argument `nodataValueAtSea` to decide whether sea areas are masked out.
    Depends on the quality of the sea mask in the input DEM.
  + automatically download required Sentinel-1 Orbit State Vector (OSV) files.
  + new argument `allow_RES_OSV` to decide whether to allow usage of the less accurate Sentinel-1 RES OSV files in
    case the POE file is not available yet.
  + new argument `demName` to choose the type of the auto-downloaded DEM.

Auxiliary Data Handling
-----------------------

- class :class:`pyroSAR.S1.OSV`:

  + removed progressbar from method :meth:`~pyroSAR.S1.OSV.catch` and made it optional in method
    :meth:`~pyroSAR.S1.OSV.retrieve` with new argument `pbar`

general
-------
- bug fixes, new automated tests, documentation improvements
