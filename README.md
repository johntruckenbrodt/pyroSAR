<h1 align="center">
  <br>
  <a>pyroSAR</a>
</h1>
<h3 align="center">A Python Framework for Large-Scale SAR Satellite Data Processing</h3>

<p align="center">
  <a href="#description">Description</a> •
  <a href="#installation">Installation</a> •
  <a href="#dependencies">Dependencies</a> •
  <a href="#example">Example</a> •
  <a href="#notes">Notes</a> •
  <a href="#authors">Authors</a>
</p>

<p align="center">
  <a href="https://www.travis-ci.org/johntruckenbrodt/pyroSAR"><img src="https://www.travis-ci.org/johntruckenbrodt/pyroSAR.svg?branch=master", alt="Travis Status"></a>
  <a href='https://ci.appveyor.com/project/johntruckenbrodt/pyrosar'><img src='https://ci.appveyor.com/api/projects/status/won0layps8mkss9h?svg=true' alt='Appveyor Status' /></a>
  <a href='https://coveralls.io/github/johntruckenbrodt/pyroSAR?branch=master'>
    <img src='https://coveralls.io/repos/github/johntruckenbrodt/pyroSAR/badge.svg?branch=master' alt='Coveralls Status' />   </a>
  <a href='https://pyrosar.readthedocs.io/en/latest/?badge=latest'>
    <img src='https://readthedocs.org/projects/pyrosar/badge/?version=latest' alt='Documentation Status' /></a>
  <a href='https://badge.fury.io/py/pyroSAR'>
    <img src='https://badge.fury.io/py/pyroSAR.svg' alt='PIP Status' /></a>
</p>

# Description

The pyroSAR package aims at providing a complete solution for the scalable organization and processing of SAR satellite data:
* Reading of data from various past and present satellite missions
* Handling of acquisition metadata
* User-friendly access to processing utilities in SNAP and GAMMA Remote Sensing software
* Formatting of the preprocessed data for further analysis

[Here](https://pyrosar.readthedocs.io/en/latest/?badge=latest) you can find the documentation.

pyroSAR is/was used in these projects:
* [BACI](http://www.baci-h2020.eu/index.php/Main/HomePage)
* [CCI Biomass](http://cci.esa.int/biomass)
* [GlobBiomass](http://globbiomass.org/)
* [SALDI](https://www.saldi.uni-jena.de/)
* [SenThIS](https://eos-jena.com/en/projects/)
* [Sentinel4REDD](https://www.dlr.de/rd/en/Portaldata/28/Resources/dokumente/re/Projektblatt_Sentinel4REDD_engl.pdf)
* [SWOS](https://www.swos-service.eu/)

You know of other projects? We'd be happy to know.

# Installation
First we need to install pip. On Windows pip is installed together with Anaconda, on Linux you can easily install it
 via command line:
```sh
sudo apt-get install python-pip
```
The latest stable release of pyroSAR can then be installed:
```sh
python -m pip install pyroSAR
```

For installation of the latest master branch on GitHub, we need the version control system git. On Windows, git can be 
downloaded from [here](https://git-scm.com/downloads). On Linux you can install it via command line:
```sh
sudo apt-get install git
```
Once everything is set up, pyroSAR is ready to be installed:
```sh
python -m pip install git+https://github.com/johntruckenbrodt/pyroSAR.git
```

# Dependencies
If you are using Windows, the easiest way to work with pyroSAR and Python in general is by using 
[Anaconda](https://www.anaconda.com/download/). It comes with all basic requirements of pyroSAR.
The more specific instructions below are intended for Linux users.
### GDAL
pyroSAR requires GDAL version 2.1 with GEOS and PROJ4 as dependencies as well as the GDAL Python binding.
#### Ubuntu
Starting with release Yakkety (16.10), Ubuntu comes with GDAL >2.1. 
See <a href="https://launchpad.net/ubuntu/yakkety/amd64/gdal-bin">here</a>. 
You can install it like this:
```bash
sudo apt-get install python-gdal python3-gdal gdal-bin
```
For older Ubuntu releases you can add the ubuntugis repository to apt prior to installation to install version >2.1:
```sh
sudo add-apt-repository ppa:ubuntugis/ppa
sudo apt-get update
```
This way the required dependencies (GEOS and PROJ4 in particular) are also installed.
You can check the version by typing:
```sh
gdalinfo --version
```
#### Debian
Starting with Debian 9 (Stretch) GDAL is available in version >2.1 in the official repository.
#### Building from source
Alternatively, you can build GDAL and the dependencies from source. The script `pyroSAR/install/install_deps.sh` 
gives specific instructions on how to do it. It is not yet intended to run this script via shell, but rather to 
follow the instructions step by step.
### SQLite + SpatiaLite
While `sqlite3` and its Python binding are usually already installed, the `spatialite` extension needs to be 
added. Two packages exist, `libspatialite` and `mod_spatialite`. Both can be used by pyroSAR.
On Ubuntu, `mod_spatialite` has been found to be easier to setup with `sqlite` and can be installed via `apt`:
```sh
sudo apt-get install libsqlite3-mod-spatialite
```
On CentOS, `libspatialite` including shared objects for extension loading can be installed via `yum`:
```sh
sudo yum install libspatialite-devel
```
The following can be run in Python to test the needed functionality:
```Python
import sqlite3
# setup an in-memory database
con=sqlite3.connect(':memory:')
# enable loading extensions and load spatialite
con.enable_load_extension(True)
try:
    con.load_extension('mod_spatialite.so')
except sqlite3.OperationalError:
    con.load_extension('libspatialite.so')
```
In case loading extensions is not permitted you might need to install the package `pysqlite2` 
together with a static build of `sqlite3`. 
See the script `pyroSAR/install/install_deps.sh` for instructions. 
There you can also find instructions on how to install `spatialite` from source.
To test `pysqlite2` you can import it as follows and then run the test above:
```Python
from pysqlite2 import dbapi2 as sqlite3
```
Installing this package is likely to cause problems with the `sqlite3` library installed on the system. 
Thus, it is safer to build a static `sqlite3` library for it (see installation script).


#  Example
Now that everything is installed, we can start working with our satellite data.
Let's assume you have a Sentinel-1 scene in a local directory. 
At first we load the scene into pyroSAR for analysis of the metadata:
```python
from pyroSAR import identify
name = 'S1A_IW_GRDH_1SDV_20150222T170750_20150222T170815_004739_005DD8_3768.zip'
scene = identify(name)
print(scene)
```
This will automatically identify the scene, scan it for metadata and print a summary of selected metadata entries.
The names of the attributes (e.g. sensor and acquisition_mode) are standardized for all SAR scenes.
Further entries, whose names are not standardized, can be found in a dictionary `scene.meta`.

Now that we have made ourselves familiar with the scene, we import it into a sqlite database:
```python
from pyroSAR import Archive
dbfile = 'scenes.db'
with Archive(dbfile) as archive:
    archive.insert(scene)
```
`dbfile` is a file either containing an already existing database or one to be created.

Let's assume our database contains a number of scenes and we want to select some for processing.  
We have a shapefile, which contains a geometry delimiting our test site for which we want to 
process some Sentinel-1 scenes.  
We already processed some scenes in the past and the results are stored in a directory
`outdir`. We only want to select scenes which have not been processed to this directory before.  
Furthermore, we are only interested in scenes acquired in Ground Range Detected (GRD) Interferometric Wide 
Swath mode (IW), which contain a VV band.

```python
from spatialist import Vector
archive = Archive('scenes.db')
site = Vector('site.shp')
outdir = '/path/to/processed/results'
maxdate = '20171231T235959'
selection_proc = archive.select(vectorobject=site,
                                processdir=outdir,
                                maxdate=maxdate,
                                sensor=('S1A', 'S1B'),
                                product='GRD',
                                acquisition_mode='IW',
                                vv=1)
archive.close()
```
Here we use pyroSAR's own vector geometry driver for loading the shapefile and pass it, together with the other parameters,
to the method `Archive.select`. You can also use the `with` statement like in the code block above.
The returned `selection_proc` is a list of file names for the scenes we selected from the database, which we can now 
pass to a processing function:
```python
from pyroSAR.snap import geocode

# the target pixel resolution in meters
resolution = 20

for scene in selection_proc:
    geocode(infile=scene, outdir=outdir, tr=resolution, scaling='db', shapefile=site)
```
The function `geocode` is a basic utility for SNAP. It will perform all necessary steps to subset, resample, orthorectify,
topographically normalize and scale the input image and write GeoTiff files to the selected output directory.  
All necessary files like orbit state vectors and SRTM DEM tiles are downloaded automatically in the background by SNAP.  
SNAP is most conveniently used with workflow XMLs. The function geocode parses a workflow for the particular scene,
parametrizes it depending on the scene type and selected processing parameters and writes it to the output directory.  
It then calls the command `gpt`, which is SNAP's command line interface, on the workflow to execute the processing steps. 

# Notes
### A Word on File Naming
pyroSAR internally uses a fixed naming scheme to keep track of processed results. For each scene an identifier is created,
which contains the sensor, acquisition mode, orbit (ascending or dsescending) and the time stamp of the acquisition start.
For the example above it is `S1A__IW___A_20150222T170750`, which can be created by calling `scene.outname_base()`. For each
attribute a fixed number of digits is reserved. In case the attribute is shorter than this number, 
the rest of the digits is filled with underscores. I.e., the sensor field is four digits long, but 'S1A' only three.
Thus, `S1A_` is the sensor slot. In the same way, `IW__` is the acquisition mode slot, which is also four digits long.  
`A` denotes ascending orbit, the time stamp is in format YYYYmmddTHHMMSS.
  
Processing functions like `geocode` add suffixes to this identifier to further keep track of individual processing
steps performed on the dataset.  
This core concept is used by many pyroSAR functions internally to keep track of which scenes have been processed before.

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

# Authors
* **John Truckenbrodt** (john.truckenbrodt@uni-jena.de)
* **Felix Cremer** (felix.cremer@uni-jena.de)
* **Ismail Baris** (i.baris@outlook.de)
