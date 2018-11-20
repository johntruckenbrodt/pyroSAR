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

changes related to Python package integrity
*******************************************
- add trove classifiers for supported operating systems and MIT license for easier online search
- exchange http with https for all URLs that support it
