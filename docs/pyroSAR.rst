Drivers
=======

This is the core module of package pyroSAR. It contains the drivers for the different SAR image formats and offers
functionality for retrieving metadata, unpacking images, downloading ancillary files like DEMs and
Orbit State Vector files as well as archiving scenes in a database. The ID class and its subclasses allow easy and
standardized access to the metadata of SAR scenes of different sensors.

.. autosummary::
    :nosignatures:

    ~pyroSAR.drivers.ID
    ~pyroSAR.drivers.Archive
    ~pyroSAR.drivers.CEOS_PSR
    ~pyroSAR.drivers.CEOS_ERS
    ~pyroSAR.drivers.ESA
    ~pyroSAR.drivers.SAFE
    ~pyroSAR.drivers.TSX
    ~pyroSAR.drivers.identify
    ~pyroSAR.drivers.identify_many
    ~pyroSAR.drivers.filter_processed

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
