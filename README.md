# pyroSAR #

### a Python framework for large-scale SAR satellite data processing ###

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


Just ask John.