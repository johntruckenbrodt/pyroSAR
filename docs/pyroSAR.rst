Drivers
=======

.. automodule:: pyroSAR.drivers
    :members:
    :undoc-members:
    :show-inheritance:

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

SNAP Processing
===============

.. automodule:: pyroSAR.snap.util
    :members:
    :undoc-members:
    :show-inheritance:

GAMMA Processing
================

.. automodule:: pyroSAR.gamma
    :members: geocode, convert2gamma, ISPPar, process, ovs, S1_deburst, correctOSV
    :undoc-members:
    :show-inheritance:

GAMMA Command API
-----------------

Parser Documentation
********************

.. automodule:: pyroSAR.gamma.parser
    :members:
    :undoc-members:
    :show-inheritance:

API Demo
********

This is a demonstration of an output script as generated automatically by function
:func:`~pyroSAR.gamma.parser.parse_module` for the Gamma module `ISP`. The command name and all parameters are passed
to function :func:`~pyroSAR.gamma.process`, which converts all input to str and then calls the command via the
:mod:`subprocess` module.

.. automodule:: pyroSAR.gamma.parser_demo
    :members:
    :undoc-members:
    :show-inheritance:

Sentinel-1 Tools
================

.. automodule:: pyroSAR.S1.auxil
    :members: OSV
    :undoc-members:
    :show-inheritance:

Datacube Tools
==============
.. automodule:: pyroSAR.datacube_util
    :members:
    :undoc-members:
    :show-inheritance:

Ancillary Functions
===================

.. automodule:: pyroSAR.ancillary
    :members:
    :undoc-members:
    :show-inheritance:
