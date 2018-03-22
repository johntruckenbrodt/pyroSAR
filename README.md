# pyroSAR
[![Build Status](https://www.travis-ci.org/johntruckenbrodt/pyroSAR.svg?branch=master)](https://www.travis-ci.org/johntruckenbrodt/pyroSAR)      [![Coverage Status](https://coveralls.io/repos/github/johntruckenbrodt/pyroSAR/badge.svg?branch=master)](https://coveralls.io/github/johntruckenbrodt/pyroSAR?branch=master)

### a Python framework for large-scale SAR satellite data processing

The pyroSAR package aims at providing a complete solution for the scalable organization and processing of SAR satellite data:
* Reading of data from various past and present satellite missions
* Handling of acquisition metadata
* User-friendly access to processing utilities in SNAP and GAMMA Remote Sensing software
* Formatting of the preprocessed data for further analysis

This software used to run on a local server and is currently being restructured into a Python package.
Not everything is working properly, stay tuned...


### Installation of dependencies

##### GDAL
pyroSAR requires GDAL version 2.1 with GEOS and PROJ4 as dependencies as well as the GDAL Python binding. Alternatively, one can use <a href="https://github.com/nextgis/pygdal">pygdal</a>, a virtualenv and setuptools friendly version of standard GDAL python bindings.
###### Ubuntu
Currently Ubuntu comes with GDAL 1.11. By adding the ubuntugis repository to e.g. apt you can install 
version >2.1:
```sh
sudo add-apt-repository ppa:ubuntugis/ppa
sudo apt-get update
sudo apt-get install python-gdal python3-gdal gdal-bin
```
This way the required dependencies (GEOS and PROJ4 in particular) are also installed.
You can check the version by typing:
```sh
gdalinfo --version
```
###### Debian
Starting with Debian 9 (Stretch) GDAL is available in version >2.1 in the official repository.
###### Building from source
Alternatively, you can build GDAL and the dependencies from source. The script `pyroSAR/install/install_deps.sh` 
gives specific instructions on how to do it. It is not yet intended to run this script via shell, but rather to 
follow the instructions step by step.
##### SQLite + SpatiaLite
While sqlite3 and its Python binding are usually already installed, the spatialite extension needs to be 
added. Two packages exist, libspatialite and mod_spatialite. Both can be used by pyroSAR.
mod_spatialite has been found to be easier to setup with sqlite and can be installed via apt:
```sh
sudo apt-get install libsqlite3-mod-spatialite
```

The following can be run in Python to test the needed functionality:
```Python
import sqlite3
# setup an in-memory database
con=sqlite3.connect(':memory:')
# enable loading extensions and load spatialite
con.enable_load_extension(True)
try:
    con.load_extension('mod_spatialite')
except sqlite3.OperationalError:
    con.load_extension('libspatialite')
```
In case loading extensions is not permitted you might need to install the package `pysqlite2`. 
See the script `pyroSAR/install/install_deps.sh` for instructions. 
There you can also find instructions on how to install spatialite from source.
To test `pysqlite2` you can import it as follows and then run the test above:
```Python
from pysqlite2 import dbapi2 as sqlite3
```
Installing this package is likely to cause problems with the sqlite3 library installed on the system. 
Thus, it is safer to build a static sqlite3 library for it (see installation script).
### Installation of pyroSAR
Once everything is set up, pyroSAR is ready to be installed:
```sh
sudo pip install git+https://github.com/johntruckenbrodt/pyroSAR.git
```
You might need to install pip and git for this to work:
```sh
sudo apt-get install python-pip
sudo apt-get install git
```


### Long Description

The launch of recent satellite missions, the Sentinel fleet of ESA’s Copernicus programme in particular, has led to a
tremendous increase in available earth observation data provided at no cost. The increase in data availability opens up
new opportunities for analysing data not only in the spatial but also temporal domain by observing time series and thus
the possibility to visualise processes on the earth’s surface. Although this is not entirely new to optical satellite
data, Synthetic Aperture Radar (SAR) data was only infrequently acquired in the past. With the new ESA SAR satellites
Sentinel-1A and Sentinel-1B there is now the possibility to observe the earth’s surface with a repeat rate of up to six
days and a spatial resolution of 20 m independent of atmospheric effects and sun illumination.
Together with the increase in data availability comes the challenge of organizing data and preparing it for scientific
analysis. While traditional software aimed at analysing single images, the need arises for fully automated organization
of large image archives with thousands of images together with a highly capable processing framework to make full use of
available hardware resources.
The pyroSAR environment aims at providing a complete solution for the organization and processing of SAR satellite data
for applications scalable from personal computers to large server infrastructures using various open source tools and
libraries. Its purpose is to provide complex functionality for reading various data formats from past and present
satellite missions, handling metadata about acquisition characteristics in a database, and providing homogenised
user-friendly access to processing utilities in ESA’s Sentinel Application Platform (SNAP) as well as GAMMA Remote
Sensing software.
The data reader uses the Geo Data Abstraction Library (GDAL) where possible and own implementations otherwise.
The metadata attributes are homogenised to enable database access of specific acquisition characteristics across
different sensor platforms. The metadata is ingested into a SpatiaLite database from which the original imagery can be
selected for processing.
The SAR processor provides functionality to distribute the tasks on different computing cores as well as different
server nodes. By following a stringent naming scheme of processed images as well as annotated metadata XMLs, processing
can be organized to be performed by several operators in a server network. This way, redundant usage of disk space and
processing resources can be reduced.
Once the images are processed, further functionalities are available for mosaicking and resampling images to common
pixel boundaries suited for time series analysis. Thus, the scientist is provided with data stacks cropped to the study
area and directly formatted for analysis without spending time with SAR-specific processing and general data management
issues.
This software is currently being developed within EU Horizon-2020 project ‘Satellite-based Wetland Observation Service
(SWOS)’. In an effort to better monitor wetlands from space with both optical and radar data, the dense time series of
the Sentinel-1 satellites is exploited in order to derive high temporal resolution surface water dynamics maps. This is
realized by applying clustering techniques in the temporal image domain of all available datasets.
