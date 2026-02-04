#########
Changelog
#########

0.6 | 2018-11-20
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

0.7 | 2019-01-03
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

0.8 | 2019-02-11
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

0.9 | 2019-06-15
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

0.9.1 | 2019-07-05
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

0.10 | 2019-12-06
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

0.10.1 | 2019-12-12
===================

GAMMA API
---------

- :ref:`Command API <gamma-command-api>` compatibility with GAMMA version 20191203

0.11 | 2020-05-29
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

0.11.1 | 2020-07-17
===================

- bug fixes

GAMMA API
---------

- :ref:`Command API <gamma-command-api>` compatibility with GAMMA version 20200713

0.12 | 2021-02-19
=================

Drivers
-------

- :class:`pyroSAR.drivers.Archive`:

  + new argument `cleanup` to automatically remove missing scenes from database on initialization
  + method :meth:`~pyroSAR.drivers.Archive.insert`: improved insertion speed
  + method :meth:`~pyroSAR.drivers.Archive.select_duplicates`: new argument `value`
  + method :meth:`~pyroSAR.drivers.Archive.get_colnames`: new argument `table` to get column names from arbitrary
    tables, not just the main `data` table
  + method :meth:`~pyroSAR.drivers.Archive.drop_element`: option to remove scene from `data` and `duplicates` tables
    simultaneously by removing argument `table` and adding argument `with_duplicates`
  + method :meth:`~pyroSAR.drivers.Archive.drop_table`:

    * new argument `verbose`
    * remove arbitrary tables, not just `data` and `duplicates`

  + method :meth:`~pyroSAR.drivers.Archive.drop_database`: replaced by new function :func:`pyroSAR.drivers.drop_archive`
  + new method :meth:`~pyroSAR.drivers.Archive.add_tables` to add custom tables to a database
  + bug fixes

- :class:`pyroSAR.drivers.CEOS_PSR`:

  + added support for ALOS-1 PALSAR
  + added basic support for Level 1.0 data

- :class:`pyroSAR.drivers.SAFE`:

  + method :meth:`~pyroSAR.drivers.SAFE.getOSV`: new argument `useLocal` to not search online if local matching
    files are found

GAMMA API
---------

- :ref:`Command API <gamma-command-api>` compatibility with GAMMA version 20201216

- function :func:`pyroSAR.gamma.convert2gamma`:

  + renamed argument `S1_noiseremoval` to `S1_tnr` (thermal noise removal)
  + new argument `S1_bnr` (border noise removal)

- function :func:`pyroSAR.gamma.geocode`:

  + new default ``removeS1BorderNoiseMethod='gamma'``
  + renamed argument `tempdir` to `tmpdir`

SNAP API
--------

- function :func:`pyroSAR.snap.util.geocode`:

  + enable grid alignment with new arguments `alignToStandardGrid`, `standardGridOriginX` and `standardGridOriginY`
  + new argument `tmpdir` to choose the location of temporarily created files
  + bug fixes

- function :func:`pyroSAR.snap.auxil.gpt`:

  + perform custom pyroSAR S1 GRD border noise removal only if IPF<2.9

Auxiliary Data Handling
-----------------------

- function :func:`pyroSAR.auxdata.dem_autoload`: return `None` if a VRT was defined

0.12.1 | 2021-03-09
===================

SNAP API
--------

- function :func:`pyroSAR.snap.util.geocode`:

  + output both sigma0 and gamma0 via argument `refarea`
  + new `export_extra` option 'layoverShadowMask'

- numerous bug fixes and API improvements

Auxiliary Data Handling
-----------------------

- class :class:`pyroSAR.S1.OSV`:

  + download files from https://scihub.copernicus.eu/gnss

0.13 | 2021-09-10
=================

Drivers
-------

- new class :class:`pyroSAR.drivers.EORC_PSR`
- new argument `exist_ok` for ID object unpack methods to enable reuse of already unpacked scenes
- :meth:`pyroSAR.drivers.SAFE.getOSV`: new argument `url_option` to choose between different download URLs
- :class:`pyroSAR.drivers.SAFE` align coordinate sorting of attribute `meta['coordinates']` with CRS description
- :func:`pyroSAR.drivers.identify_many`: disable progressbar by default

GAMMA API
---------

- adaptations to enable processing of :class:`~pyroSAR.drivers.EORC_PSR` data:

  + :func:`pyroSAR.gamma.calibrate`
  + :func:`pyroSAR.gamma.convert2gamma`
  + :func:`pyroSAR.gamma.geocode`

- :func:`pyroSAR.gamma.geocode`:

  + experimental optional refinement of the geocoding lookup table with new argument `refine_lut`
  + removed arguments `normalization_method`, `func_interp`, `removeS1BorderNoise`, `sarSimCC`
  + limit radiometric normalization to RTC correction method
  + simplify and improve computation of RTC contribution area
  + file suffices `pan` and `norm` have been replaced with `gamma0-rtc`
  + argument `export_extra` options:

    * removed `pix_geo`
    * renamed `pix_fine` to `pix_ratio`
    * added `pix_area_sigma0`, `pix_area_sigma0_geo`, `pix_area_gamma0_geo`, `gs_ratio` , `gs_ratio_geo`, `pix_ratio_geo`

  + use a dedicated temporary directory to unpack the scene and write GAMMA files so that they are separated (the GAMMA
    files used to be written to the unpacked scene's directory)
  + enable multiple scenes as input so that they can be mosaiced in SAR geometry before geocoding

- :func:`pyroSAR.gamma.correctOSV`: new argument `directory`

- :func:`pyroSAR.gamma.multilook`: new argument `exist_ok`

- :func:`pyroSAR.gamma.convert2gamma`: new argument `exist_ok`

- function :func:`pyroSAR.gamma.dem.dem_autocreate`:

  + do not apply an extent buffer by default
  + allow geometry in arbitrary CRS

SNAP API
--------

- function :func:`pyroSAR.snap.util.geocode`:

  + new `export_extra` option `scatteringArea`

- extended support for `BandMaths` operator

Auxiliary Data Handling
-----------------------

- method :meth:`pyroSAR.S1.OSV.catch`: new argument `url_option` with two download URLs to choose from

- function :func:`pyroSAR.auxdata.dem_autoload`:

  + added new DEM option `GETASSE30`
  + align pixels of subsetted VRT with original tiles

- function :func:`pyroSAR.auxdata.dem_create`:

  + new argument `outputBounds`

general
-------
- replaced print messages with logging. This made the `verbose` argument that was used by several functions and
  methods obsolete; affects the following:

  + :func:`pyroSAR.drivers.identify_many`: replaced by argument `pbar`
  + :meth:`pyroSAR.drivers.Archive.add_tables`: removed
  + :meth:`pyroSAR.drivers.Archive.drop_table`: removed
  + :meth:`pyroSAR.drivers.Archive.insert`: replaced by argument `pbar`
  + :meth:`pyroSAR.drivers.Archive.import_outdated`: removed
  + :meth:`pyroSAR.drivers.Archive.move`: replaced by argument `pbar`
  + :meth:`pyroSAR.drivers.Archive.select`: removed
  + :func:`pyroSAR.snap.auxil.execute`: removed

  See section :doc:`Logging </general/logging>` for details.

0.14.0 | 2021-10-12
===================

Drivers
-------
- raise more appropriate errors (`c430c59 <https://github.com/johntruckenbrodt/pyroSAR/commit/c430c59289016b5fe2e0f3044225dc5166c39e80>`_)

- :func:`pyroSAR.drivers.findfiles`: removed (functionality contained in :meth:`pyroSAR.drivers.ID.findfiles`,
  now making use of :func:`spatialist.ancillary.finder`)

- :meth:`pyroSAR.drivers.Archive.select`:

  + show progressbar for scene identification if ``pbar=True``
  + enabled input of :obj:`~datetime.datetime` objects for arguments ``mindate`` and ``maxdate``

- :func:`pyroSAR.drivers.identify_many`: issue a warning when a file cannot be accessed
  (instead of raising a :obj:`PermissionError`)

GAMMA API
---------
- :func:`pyroSAR.gamma.dem.dem_autocreate`: support for new DEM options provided by :func:`pyroSAR.auxdata.dem_autoload`

SNAP API
--------
- :func:`pyroSAR.snap.auxil.get_egm96_lookup` removed in favor of new function :func:`pyroSAR.auxdata.get_egm_lookup`

Auxiliary Data Handling
-----------------------
- method :meth:`pyroSAR.S1.OSV.retrieve`: thread-safe writing of orbit files

- new function :func:`pyroSAR.auxdata.get_egm_lookup`

- function :func:`pyroSAR.auxdata.dem_create`

  + new geoid option 'EGM2008'
  + make use of :func:`~pyroSAR.auxdata.get_egm_lookup` for auto-download of EGM lookup files
  + several bug fixes related to vertical CRS transformation
  + bug fix for target pixel alignment

- function :func:`pyroSAR.auxdata.dem_autoload`: new DEM options:

  + 'Copernicus 10m EEA DEM'
  + 'Copernicus 30m Global DEM'
  + 'Copernicus 90m Global DEM'

general
-------
- replaced http URLs with https where applicable
- improved documentation

0.15.0 | 2022-01-04
===================

Drivers
-------
- :meth:`pyroSAR.drivers.ID.geometry`: new method

GAMMA API
---------
- :ref:`Command API <gamma-command-api>` compatibility with GAMMA version 20211208

- renamed argument `resolution` to `spacing`; affects:

  + :func:`pyroSAR.gamma.geocode`
  + :func:`pyroSAR.gamma.ovs`
  + :func:`pyroSAR.gamma.multilook`

- function :func:`pyroSAR.gamma.calibrate`

  + removed argument `replace`
  + added argument `return_fnames`

- function :func:`pyroSAR.gamma.convert2gamma`

  + added argument `return_fnames`

- function :func:`pyroSAR.gamma.multilook`

  + pass multiple Sentinel-1 sub-swaths to argument `infile` which are then
    combined into a single MLI using GAMMA command `isp.multi_look_ScanSAR`

- class :class:`pyroSAR.gamma.ISPPar`:

  + new object attribute `filetype` with possible values 'isp' and 'dem'

SNAP API
--------
- function :func:`pyroSAR.snap.util.geocode`:

  + enabled SLC processing
  + enable processing of sigma nought RTC
  + new `export_extra` argument `gammaSigmaRatio`
  + simplified workflow by writing layover-shadow mask directly from `Terrain-Correction`
  + changed processing node sequence:

    * was: Read->ThermalNoiseRemoval->SliceAssembly->Remove-GRD-Border-Noise->Calibration
    * is:  Read->Remove-GRD-Border-Noise->Calibration->ThermalNoiseRemoval->SliceAssembly

  + new output image naming scheme, e.g.

    * S1A__IW___A_20210914T191350_VV_gamma0-rtc.tif
    * S1A__IW___A_20210914T191350_VH_sigma0-elp.tif

- function :func:`pyroSAR.snap.auxil.gpt`:

  + removed argument `multisource`
  + added argument `tmpdir`

Auxiliary Data Handling
-----------------------
- function :func:`pyroSAR.auxdata.dem_autoload`:

  + updated version of 'Copernicus 10m EEA DEM' from '2020_1' to '2021_1'
  + new DEM options:

    * 'Copernicus 30m Global DEM II'
    * 'Copernicus 90m Global DEM II'

general
-------
- compatibility with sqlalchemy>=1.4

0.15.1 | 2022-01-07
===================
general
-------
- bug fixes

0.16.0 | 2022-03-03
===================

Drivers
-------
- :class:`pyroSAR.drivers.BEAM_DIMAP`: new driver supporting SNAP's BEAM-DIMAP format
- :class:`pyroSAR.drivers.SAFE`:

  + corrected SLC metadata (was read from first sub-swath, now from center sub-swath or as sum of all sub-swaths):
    center: spacing, heading, incidence; sum: samples, lines
  + new property :attr:`pyroSAR.drivers.SAFE.resolution`

Auxiliary Data Handling
-----------------------
- create water body mask mosaics from ancillary DEM products. Affects the following:

  + function :func:`pyroSAR.auxdata.dem_autoload`: new arguments `nodata` and `hide_nodata`

- function :func:`pyroSAR.auxdata.dem_create`:

  + new arguments `pbar` and `threads`

SNAP API
--------
- new method :meth:`pyroSAR.snap.auxil.Par_BandMath.add_equation`
- new function :func:`pyroSAR.snap.util.noise_power`
- new function :func:`pyroSAR.snap.auxil.erode_edges`
- function :func:`pyroSAR.snap.auxil.writer`:

  + new arguments `clean_edges` and `clean_edges_npixels`
    (to make use of function :func:`~pyroSAR.snap.auxil.erode_edges`)
  + enabled conversion of BEAM-DIMAP files

- function :func:`pyroSAR.snap.util.geocode`:

  + new arguments `clean_edges` and `clean_edges_npixels` (see function :func:`~pyroSAR.snap.auxil.writer`)
  + renamed argument `tr` to `spacing`
  + new arguments `rlks` and `azlks` to manually set the number of looks

GAMMA API
---------
- function :func:`pyroSAR.gamma.geocode`:

  + new arguments `rlks` and `azlks`

- function :func:`pyroSAR.gamma.multilook`:

  + new arguments `rlks` and `azlks`

general
-------
- correction of multi-look factor computation. Before: approximate target pixel spacing but never exceed it.
  Now: first best approximate the azimuth spacing as close as possible (even if this means exceeding the target spacing)
  and then choose the range looks to approximate a square pixel as close as possible. API changes:

  + function :func:`pyroSAR.ancillary.multilook_factors`:

    * renamed argument `sp_rg` to `source_rg`
    * renamed argument `sp_az` to `source_az`
    * replaced arguments `tr_rg` and `tr_az` with unified `target`

0.16.1 | 2022-03-07
===================

Auxiliary Data Handling
-----------------------
- function :func:`pyroSAR.auxdata.get_egm_lookup`:

  + changed URL for PROJ geoid models, which results in better performance for
    function :func:`pyroSAR.auxdata.dem_create`
    (See `pyroSAR#200 <https://github.com/johntruckenbrodt/pyroSAR/issues/200>`_).

0.16.2 | 2022-03-14
===================

SNAP API
--------
- function :func:`pyroSAR.snap.util.noise_power`: added missing orbit state vector refinement

0.16.3 | 2022-03-23
===================

SNAP API
--------
- function :func:`pyroSAR.snap.util.noise_power`: pass argument `cleanup` to :func:`~pyroSAR.snap.auxil.gpt` call
- function :func:`~pyroSAR.snap.auxil.gpt`: shortened names of temporary directories
- function :func:`~pyroSAR.snap.auxil.erode_edges`: fixed bug in polygon selection
- function :func:`~pyroSAR.snap.auxil.writer`: do not erode edges of layover-shadow mask

0.17.0 | 2022-05-30
===================

SNAP API
--------

- function :func:`pyroSAR.snap.erode_edges`: reuse mask for all images

GAMMA API
---------

- new function :func:`pyroSAR.gamma.dem.dem_import`

- function :func:`pyroSAR.gamma.geocode`:

  + new argument `update_osv`

general
-------

- full support for Sentinel-1 stripmap mode; renamed `SM` naming pattern to `S1..S6` to differentiate different beams
- bug fixes

0.17.2 | 2022-06-23
===================

Auxiliary Data Handling
-----------------------
- function :func:`pyroSAR.auxdata.dem_create`:

  + use maximum possible value of `dtype` (e.g. 255 for unit8) instead of -32767.0 if the nodata value cannot be read from the source file
  + always use the same value for source and destination nodata

0.17.3 | 2022-07-03
===================

Auxiliary Data Handling
-----------------------
- function :func:`pyroSAR.auxdata.dem_create`:

  + In case the nodata value could not be read from the source file, the function used to define a value itself, which is prone to errors. This value now needs to be set by a user via new argument `nodata` if it cannot be read from the source file.
  + bug fix: no longer try to download 'Copernicus 30m Global DEM' or 'Copernicus 90m Global DEM' tiles that don't exist.

- function :func:`pyroSAR.auxdata.dem_autoload`:

  + new argument `dst_nodata`. This can be used to temporarily override the native nodata value for extrapolation of ocean areas (in combination with ``hide_nodata=True``).

0.18.0 | 2022-08-24
===================

Drivers
-------
- method :meth:`pyroSAR.drivers.SAFE.quicklook`: new argument `na_transparent`
- new class :class:`~pyroSAR.drivers.TDM`
- method :meth:`pyroSAR.drivers.TSX.getCorners`: fixed bug in longitude computation
- class :class:`~pyroSAR.drivers.ESA`: improved support for ERS and ASAR


GAMMA API
---------
- :ref:`Command API <gamma-command-api>` compatibility with GAMMA version 20220629

SNAP API
--------
- compatibility with SNAP version 9
- function :func:`~pyroSAR.snap.util.geocode`: improved support for ERS and ASAR

0.19.0 | 2022-09-28
===================

Drivers
-------
- class :class:`pyroSAR.drivers.ESA`: added support for ASAR WSM

SNAP API
--------
- new convenience functions:

  + :func:`pyroSAR.snap.auxil.geo_parametrize`
  + :func:`pyroSAR.snap.auxil.sub_parametrize`
  + :func:`pyroSAR.snap.auxil.mli_parametrize`
  + :func:`pyroSAR.snap.auxil.dem_parametrize`

- function :func:`pyroSAR.snap.auxil.orb_parametrize`: removed args `workflow`, `before`, `continueOnFail`; added `kwargs`
- function :func:`pyroSAR.snap.auxil.erode_edges`: extended to also take a BEAM-DIMAP product as input or a folder of multiple ENVI files (and not just and individual ENVI file)
- function :func:`pyroSAR.snap.auxil.Workflow.insert_node`: option to insert multiple nodes at once

Auxiliary Data Handling
-----------------------
- function :func:`pyroSAR.auxdata.dem_autoload`:

  + new argument `crop` to optionally return the full extent of all overlapping DEM tiles
  + added download status print messages
  + download and modify a Copernicus DEM index file for future reuse; this removes the need to search the FTP server for files and thus greatly accelerates the process of collecting all files overlapping with the AOI

0.20.0 | 2022-12-27
===================

Drivers
-------
- class :class:`pyroSAR.drivers.ESA`: changed ASAR orbit type from DELFT to DORIS
- class :class:`pyroSAR.drivers.BEAM_DIMAP`: new attributes `meta['incidence']` and `meta['image_geometry']`
- class :class:`pyroSAR.drivers.Archive`: new argument `date_strict` for method :meth:`~pyroSAR.drivers.Archive.select`

SNAP API
--------
- function :func:`pyroSAR.snap.util.geocode`: force multi-looking for ERS1, ERS2, ASAR even if range and azimuth factor are both 1

Auxiliary Data Handling
-----------------------
- function :func:`pyroSAR.auxdata.dem_autoload`:

  + no longer require DEM tiles for creating a mosaic to address ocean cases
  + simplified handling and removed arguments `nodata`, `dst_nodata` and `hide_nodata`
  + the DEM option 'Copernicus 30m Global DEM' now also includes several auxiliary layers that can be downloaded automatically
  + the URLs for DEM options 'SRTM 3Sec' and 'TDX90m' have been updated

- function :func:`pyroSAR.auxdata.dem_create`:

  + option to customize the output DEM via additional keyword arguments to be passed to :func:`spatialist.auxil.gdalwarp`
  + no longer require a nodata value

0.21.0 | 2023-05-11
===================

Drivers
-------
- class :class:`pyroSAR.drivers.Archive`:

  + improved PostgreSQL connection stability
  + method :meth:`~pyroSAR.drivers.Archive.select`: the `vectorobject` geometry is now cloned before being reprojected to EPSG:4326 so that the source geometry remains unaltered

GAMMA API
---------
- the `LAT` module is no longer needed: new pyroSAR-internal implementations can be used if the module is missing (concerns commands `product`, `ratio` and `linear_to_dB`)
- improved backwards compatibility:

  + use `multi_look_ScanSAR` if present and `multi_S1_TOPS` otherwise
  + use `gc_map2` if possible (present and with all needed arguments) and `gc_map` otherwise
  + addressed the case where `gc_map` does not have an argument `OFF_par`

- function `gamma.pixel_area_wrap`: new argument `exist_ok` (this function will be made more visible in the documentation once matured)
- bug fixes:

  + :func:`pyroSAR.gamma.convert2gamma`: raise an error if `S1_bnr=True` but the GAMMA command does not support border noise removal
  + :func:`pyroSAR.gamma.geocode`: removed unneeded underscore in HDR file naming
  + `gamma.pixel_area_wrap`: fixed some issues with occasionally missing intermediate files, e.g. for computing ratios

SNAP API
--------
- function :func:`pyroSAR.snap.util.geocode`: new argument `dem_oversampling_multiple` with default 2 to increase the DEM oversampling factor for terrain flattening
- function :func:`pyroSAR.snap.auxil.erode_edges`:

  + do not attempt to perform erosion if the image only contains nodata (this might happen if only parts of the image were geocoded)
  + make sure that a backscatter image is used for erosion (auxiliary data like the local incidence angle often has a larger valid data extent and using such image for erosion would thus not properly erode edges of the backscatter images; additionally this has the effect that all images will have the same valid data extent after erosion)
  + the written mask files (delineating valid data and nodata after erosion of the backscatter image and used for masking all other images) are now compressed (deflate) so that data volume is decreased significantly

Auxiliary Data Handling
-----------------------
- function :func:`pyroSAR.auxdata.dem_create`:

  + new argument `resampleAlg` to change the resampling algorithm

0.22.0 | 2023-09-21
===================

Drivers
-------
- class :class:`pyroSAR.drivers.Archive`:

  + allow multiple products with same `outname_base`, e.g. Sentinel-1 GRD and SLC; this required the introduction of a second primary key in the database
  + method :meth:`~pyroSAR.drivers.Archive.import_outdated`: option to import data from an old database with only one primary key; this requires the old
    database to be opened in legacy mode (new argument `legacy=True`)

- class :class:`pyroSAR.drivers.SAFE`: support for handling Sentinel-1 OCN products (metadata reading and database handling)

Auxiliary Data Handling
-----------------------
- class :class:`pyroSAR.auxdata.DEMHandler`: enabled handling of southern hemisphere geometries.

0.22.1 | 2023-10-11
===================

Drivers
-------
- class :class:`pyroSAR.drivers.BEAM_DIMAP`: enable calling inherited method :meth:`~pyroSAR.drivers.ID.geometry`

0.22.2 | 2023-11-16
===================

SNAP API
--------
- function :func:`pyroSAR.snap.auxil.writer`: fixed bug in ignoring `erode_edges` argument
- function :func:`pyroSAR.snap.auxil.erode_edges`: enable handling of polarimetric matrices

Drivers
-------
- function :func:`pyroSAR.drivers.identify`: enable reading of :class:`~pyroSAR.drivers.TDM` products

Misc
----
- class :class:`pyroSAR.examine.ExamineGamma`: enhanced flexibility in finding GAMMA installation

0.23.0 | 2023-11-23
===================

Drivers
-------
- class :class:`pyroSAR.drivers.Archive`: fixed bug in loading spatialite on Darwin-based systems

Auxiliary Data Handling
-----------------------

changes to Sentinel-1 OSV data handling:

- method :meth:`pyroSAR.S1.OSV.catch`:

  + removed `url_option` 1 (https://scihub.copernicus.eu/gnss)
  + made option 2 the new default option 1 (https://step.esa.int/auxdata/orbits/Sentinel-1)

- added new arguments to the following functions:

  + :func:`pyroSAR.gamma.correctOSV`: `url_option`
  + :func:`pyroSAR.gamma.geocode`: `s1_osv_url_option`
  + :func:`pyroSAR.snap.auxil.orb_parametrize`: `url_option`
  + :func:`pyroSAR.snap.util.geocode`: `s1_osv_url_option`
  + :func:`pyroSAR.snap.util.noise_power`: `osv_url_option`

0.24.0 | 2024-01-10
===================

Drivers
-------
- new base attribute `coordinates`
- enable method :meth:`~pyroSAR.drivers.ID.geometry` for all driver classes
- classes :class:`~pyroSAR.drivers.ESA` and :class:`~pyroSAR.drivers.CEOS_ERS`: removed call to `gdalinfo`
  (for increased test capability and speed)
- outsourced regular expressions for product identification into separate module `patterns`

Auxiliary Data Handling
-----------------------
- method :meth:`pyroSAR.S1.OSV.catch`: fixed bug in finding files starting in previous month

0.25.0 | 2024-04-16
===================

Drivers
-------
- class :class:`pyroSAR.drivers.Archive`:

  + replaced column `bbox` with `geometry`; requires database migration
  + method :meth:`~pyroSAR.drivers.Archive.export2shp`: improved column name laundering

SNAP API
--------
- function :func:`pyroSAR.snap.auxil.gpt`: fixed bug that occurred during removal of BNR node

Ancillary Tools
---------------
- new classes :class:`pyroSAR.ancillary.Lock` and :class:`pyroSAR.ancillary.LockCollection`
  for custom file/folder locking

Auxiliary Data Handling
-----------------------
- function :func:`pyroSAR.auxdata.dem_create`:

  + make use of new classes :class:`~pyroSAR.ancillary.Lock` and :class:`~pyroSAR.ancillary.LockCollection`
    for DEM download and mosaic creation (new argument `lock_timeout`)
  + check whether all VRT source files exist

0.26.0 | 2024-05-15
===================

SNAP API
--------
- compatibility with SNAP 10.
- completely revised configuration mechanisms. See

  + :doc:`/general/configuration`
  + :class:`pyroSAR.examine.ExamineSnap`
  + :class:`pyroSAR.examine.SnapProperties`

0.26.1 | 2024-10-01
===================

Drivers
-------
- method :meth:`pyroSAR.drivers.Archive.select`: do not accept multi-feature vectorobjects

SNAP API
--------
- fixed bug in writing SNAP properties configuration

Auxiliary Data Handling
-----------------------
- class :class:`pyroSAR.auxdata.DEMHandler`: lock created VRT files

0.27.0 | 2024-12-19
===================

Auxiliary Data Handling
-----------------------
- class :class:`pyroSAR.S1.OSV`: fixed bug in searching STEP OSV repository
- function :func:`pyroSAR.auxdata.dem_create`: removed argument `lock_timeout`, do no longer lock the target file.
- function :func:`pyroSAR.auxdata.dem_autoload`: the target VRT file is no longer locked. However, the individual downloaded DEM tiles now are.

Ancillary Tools
---------------
- classes :class:`~pyroSAR.ancillary.Lock` and :class:`~pyroSAR.ancillary.LockCollection`:
  enable nested locking

Misc
----
- removed upper Python dependency limit

0.28.0 | 2025-02-20
===================

General
-------
- support for SNAP 11 (tested, no modifications necessary)
- support for Sentinel-1C and D

Drivers
-------
- function :func:`pyroSAR.drivers.identify_many`: new argument `cores` for parallel scene identification
- class :class:`pyroSAR.drivers.SAFE`: enable unzipping of products from CDSE

Auxiliary Data Handling
-----------------------
- removed option for `TDX90m` DEM download because the FTP server has been shut down
  (perhaps reactivated in the future if HTTPS authentication can be implemented)

0.29.0 | 2025-04-09
===================

General
-------
- extended support for Sentinel-1C and D

Drivers
-------
- :meth:`pyroSAR.drivers.SAFE.geo_grid`: new method

0.29.1 | 2025-05-12
===================

SNAP API
--------
- support for SNAP 12

0.30.0 | 2025-05-14
===================

Drivers
-------
- changed polygon coordinate order to counter-clockwise for methods

  - :meth:`pyroSAR.drivers.ID.bbox`
  - :meth:`pyroSAR.drivers.ID.geometry`
- method :meth:`pyroSAR.drivers.Archive.select`: new argument `return_value`

0.30.1 | 2025-08-22
===================

Drivers
-------
- :class:`~pyroSAR.drivers.ESA`: read all `GEOLOCATION GRID ADS` segments to obtain GCPs, not just the first one (bugfix)

GAMMA API
---------
- support for GAMMA version 20250625
- support for polar stereographic projections (via :meth:`~pyroSAR.gamma.auxil.ISPPar.envidict`)
- class :class:`~pyroSAR.gamma.auxil.ISPPar`: raise error if file type is unknown
  (instead of setting the `filetype` attribute to `unknown`)
- :func:`~pyroSAR.gamma.util.pixel_area_wrap`:

  + create ENVI HDR files for inputs to :func:`~pyroSAR.gamma.util.lat_ratio` (bugfix)
  + fixed bug in ignoring conditions for writing ENVI HDR files of `pix*` and `gs_ratio` products

- improved readability of tests

0.31.0 | 2025-09-23
===================

Drivers
-------
- :meth:`pyroSAR.drivers.ID.bbox`: new argument `buffer`
- :class:`~pyroSAR.drivers.SAFE`, :class:`~pyroSAR.drivers.BEAM_DIMAP`: new argument `looks`
- :class:`~pyroSAR.drivers.Archive`: context-manage all database handles (code improvement)

GAMMA API
---------
- :func:`~pyroSAR.gamma.util.convert2gamma`, :func:`~pyroSAR.gamma.util.correctOSV`: add file locking
- fixed argument names of `isp.MLI_cat`

0.32.0 | 2025-10-29
===================

SNAP API
--------
- :func:`~pyroSAR.snap.auxil.orb_parametrize`: improved ERS/ASAR orbit handling (more work necessary to always select the best available file, because all options are limited in time (e.g. use option 1 if possible, fall back to option 2 otherwise, etc.); needs a download functionality like :class:`pyroSAR.S1.auxil.OSV` to know which ones are available)
- :func:`~pyroSAR.snap.util.geocode`:

  + explicitly use 'Latest Auxiliary File' for Envisat calibration (just for readability, this is already the default value of the parsed node; other options: 'Product Auxiliary File', 'External Auxiliary File')
  + leave calibration node polarizations field empty when processing all polarizations (otherwise processing may finish without errors but no product is being written; looks like a SNAP bug, also reported in `step-44830 <https://forum.step.esa.int/t/naming-of-source-bands/44830>`_)
  + `Calibration` in/out band handling improvements

    * select source bands based on sensor and acquisition mode (also described in `step-44830 <https://forum.step.esa.int/t/naming-of-source-bands/44830>`_)
    * more explicit handling of output bands: all that are not needed set to `False`
    * commented out output bands that are apparently not needed

  + fixed sarsim-cc geocoding:

    * old: `SAR-Simulation->Cross-Correlation->Terrain-Flattening->SARSim-Terrain-Correction` (does not work because `Terrain-Flattening` does not pass through any source layers)
    * new: `SAR-Simulation->Cross-Correlation->Warp->Terrain-Flattening->Terrain-Correction`
    * this reveals a flaw in current SNAP processing: the additional `Warp` step introduces unnecessary resampling, the created lookup table is not passed between operators and thus makes the process inefficient, the whole procedure only works with EPSG:4326 as map geometry thus, by the looks of it, requiring three forward geocoding steps (for `SAR-Simulation`, `Terrain-Flattening` and `Terrain-Correction`, respectively)

- :func:`~pyroSAR.snap.auxil.groupbyWorkers`: add `Warp` operator to the group of its source node, because it cannot be executed alone (just like `ThermalNoiseRemoval`)
- ancillary layer writing fix: a layover-shadow-mask can also be created by `SAR-Simulation`, but the output layer is named differently ('layover_shadow_mask' instead of 'layoverShadowMask' by `Terrain-Correction`); this must be handled correctly in :func:`pyroSAR.snap.auxil.writer`

Drivers
-------
- :class:`~pyroSAR.drivers.ESA`:

  + :meth:`~pyroSAR.drivers.ESA.scanMetadata`:

    * read out all MPH, SPH, DSD and GEOLOCATION_GRID_ADS metadata and expose it via `meta['origin']`
    * use absolute orbit number as `frameNumber` instead of product counter (which often seems to be 0)
    * convert original metadata to Python types (int, float, datetime)
    * renamed several meta attributes:

      - `incidenceAngleMin` -> `incidence_nr`
      - `incidenceAngleMax` -> `incidence_fr`
      - `rangeResolution`, `azimuthResolution` -> `resolution` (tuple)
      - `neszNear`, `neszFar` -> `nesz` (tuple)

  + new method :meth:`~pyroSAR.drivers.ESA.geo_grid` (like for `SAFE`)
  + corrected `acquisition_mode` for ASAR WSM, WSS
  + added MR product type

- :class:`~pyroSAR.drivers.BEAM_DIMAP`

  + improved metadata parsing

    * `incidenceAngleMidSwath` not always present, use `incidence_near` and `incidence_far` alternatively
    * the cycle number may be named `orbit_cycle` or `CYCLE`
    * for pyroSAR `frameNumber`, use `ABS_ORBIT`, not `data_take_id` as for Sentinel-1
    * added further `meta` attributes: `swath`, `looks`
    * always four `Polarizations` fields present, some may be set to None -> filtered out
    * for Sentinel-1 the product and acquisition_mode attributes can be obtained from `ACQUISITION_MODE` and `PRODUCT_TYPE` respectively; for ASAR/ERS `ACQUISITION_MODE` is missing and `PRODUCT_TYPE` contains the original values, e.g. 'ASA_APP_1P' -> must be abstracted

  + added MR product type

- :class:`~pyroSAR.drivers.ID`

  + added methods `start_dt` and `stop_dt` returning timezone-aware datetime objects

Ancillary Tools
---------------

- :meth:`~pyroSAR.ancillary.multilook_factors`: fixed bug in returning 0 as range factor

0.32.1 | 2025-11-06
===================

Auxiliary Data Handling
-----------------------
- class :class:`pyroSAR.S1.OSV`: lock local target files for download (to avoid multi-download and conflicts in parallel processes)

0.33.0 | 2025-12-17
===================

Drivers
-------
- :class:`~pyroSAR.drivers.ESA`:

    + convert coordinates in `meta['origin']` to floats
    + read incident angles directly from metadata, not from custom mapping `ANGLES_RESOLUTION` (from which they have been removed)
    + `ERS.mapping` renaming:

          * `ANGLES_RESOLUTION` -> `RESOLUTION_NESZ`
          * `get_angles_resolution` -> `get_resolution_nesz`
          * `range` -> `res_rg`
          * `azimuth` -> `res_az`
          * `nesz_near` -> `nesz_nr`
          * `nesz_far` -> `nesz_fr`

    + made code more robust by reading SPH and DSD sizes from MPH
    + added WSS mode to `RESOLUTION_NESZ` (although all values are just `None` because they could not be found yet)
    + simplified code and added typing

- :class:`~pyroSAR.drivers.BEAM_DIMAP`:

    + more robust incident angle reading

SNAP API
--------
- support for SNAP 13

Ancillary Tools
---------------

- :meth:`~pyroSAR.ancillary.multilook_factors`: complete reimplementation for more robustness

Auxiliary Data Handling
-----------------------
- class :class:`pyroSAR.auxdata.DEMHandler`: handle ocean areas without DEM coverage using a dummy DEM spanning the target extent instead of the whole globe. The latter is no longer supported by GDAL.

0.33.1 | 2026-01-19
===================

Drivers
-------
- :meth:`pyroSAR.drivers.SAFE.geo_grid`: fixed datetime handling bug by requiring spatialist>=0.16.2

0.33.2 | 2026-01-21
===================

Auxiliary Data Handling
-----------------------
- :meth:`S1.OSV.__catch_step_auxdata` do not stop if no file was found on first URL

0.33.3 | 2026-01-30
===================

GAMMA API
---------
- :class:`pyroSAR.gamma.auxil.ISPPar`: fixed `date` attribute handling