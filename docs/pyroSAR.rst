Drivers
=======

This is the core module of package pyroSAR. It contains the drivers for the different SAR image formats and offers
functionality for retrieving metadata, unpacking images, downloading ancillary files like DEMs and
Orbit State Vector files as well as archiving scenes in a database. The ID class and its subclasses allow easy and
standardized access to the metadata of images from different SAR sensors.

.. automodule:: pyroSAR.drivers

    .. rubric:: classes

    .. autosummary::
        :nosignatures:

        ID
        CEOS_PSR
        CEOS_ERS
        ESA
        SAFE
        TSX
        Archive

    .. rubric:: functions

    .. autosummary::
        :nosignatures:

        identify
        identify_many
        filter_processed
        findfiles
        getFileObj
        parse_date

.. automodule:: pyroSAR.drivers
    :members:
    :undoc-members:
    :show-inheritance:

SNAP Processing
===============

.. automodule:: pyroSAR.snap.util
    :members:
    :undoc-members:
    :show-inheritance:

GAMMA Processing
================

.. automodule:: pyroSAR.gamma.util
    :members: geocode, convert2gamma
    :undoc-members:
    :show-inheritance:

Sentinel-1 Tools
================

.. automodule:: pyroSAR.S1.auxil
    :members: OSV
    :undoc-members:
    :show-inheritance:
