Changelog
=========

0.6 / 2018-11-20
----------------

SAR metadata
************
- new standardized  metadata fields `orbitNumber_abs`, `orbitNumber_rel`, `cycleNumber` and `frameNumber` for all SAR
  formats
- customization of output file names with additional metadata fields (e.g. orbit numbers)

software configuration
**********************
- pyroSAR configuration file handling: the paths to the SNAP and Gamma installation as well as relevant metadata
  directories are now registered in a configuration file `config.ini`, which is stored in a directory `.pyrosar` in the
  user home directory
- improved SNAP installation verification: pyroSAR now performs a deeper check of the SNAP installation to make sure
  it is not mistaken with e.g. the Ubuntu package manager snap; relevant installation executables and directories are
  stored in the configuration file

general functionality
*********************
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
*********************
functionality to export processed datasets directly to an Open Data Cube:
it is now possible to create Open Data Cube product YML files as well as YML files for data indexing and ingestion
into this product; pyroSAR also internally checks for compatibility of a particular dataset with the target product;
this way, the resulting files can easily be passed to the Open Data Cube command line tools
several bug fixes

SNAP API
********
improved SNAP processing workflow node linking: it is now possible to add a node also before an existing one, instead
of just after it

Python package integrity
************************
- add trove classifiers for supported operating systems and MIT license for easier online search
- exchange http with https for all URLs that support it

0.7 / 2019-01-03
----------------

several changes to the functioning of the Gamma command API

GAMMA API
*********

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
----------------

Auxiliary Data Handling
***********************

- new module auxdata with function :func:`pyroSAR.auxdata.dem_autoload` to automatically download tiles of
  different DEM types overlapping with given geometries
- class :class:`pyroSAR.S1.OSV`: reduced search time for new RES orbit state vector files;
  included more meaningful status messages

GAMMA API
*********

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
********

- function :func:`pyroSAR.snap.util.geocode`:

  + export temporarily written files (e.g. local incidence angle) via new parameter `export_extra`
