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
        BEAM_DIMAP
        CEOS_PSR
        CEOS_ERS
        EORC_PSR
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
        getFileObj
        parse_date
        drop_archive

SNAP Processing
===============

.. automodule:: pyroSAR.snap.util
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        geocode
        noise_power

Workflow Parsing and Execution
------------------------------

.. automodule:: pyroSAR.snap.auxil
    :members: gpt, execute, parse_node, parse_recipe, split, groupbyWorkers, Workflow, Node, Par, Par_BandMath
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        gpt
        execute
        parse_node
        parse_recipe
        split
        groupbyWorkers
        Workflow
        Node
        Par
        Par_BandMath

General Utilities
-----------------

.. automodule:: pyroSAR.snap.auxil
    :members: erode_edges, writer
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        erode_edges
        writer

GAMMA Processing
================

.. automodule:: pyroSAR.gamma
    :members: geocode, convert2gamma, ISPPar, process, ovs, S1_deburst, correctOSV, multilook, par2hdr, UTM, calibrate
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        calibrate
        convert2gamma
        correctOSV
        geocode
        ISPPar
        multilook
        ovs
        par2hdr
        process
        S1_deburst
        UTM

DEM tools
---------

.. automodule:: pyroSAR.gamma.dem
    :members: dem_autocreate, dem_import, dempar, fill, hgt, hgt_collect, makeSRTM, mosaic, swap
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        dem_autocreate
        dem_import
        dempar
        fill
        hgt
        hgt_collect
        makeSRTM
        mosaic
        swap

.. _gamma-command-api:

GAMMA Command API
-----------------

This is an attempt to make it easier to execute GAMMA commands by offering automatically parsed Python functions.
Thus, instead of executing the command via shell:

.. code-block:: shell

    offset_fit offs ccp off.par coffs - 0.15 3 0 > offset_fit.log

one can wrap it in a Python script:

.. code-block:: python

    import os
    from pyroSAR.gamma.api import isp

    workdir = '/data/gamma_workdir'

    parameters = {'offs': os.path.join(workdir, 'offs'),
                  'ccp': os.path.join(workdir, 'ccp'),
                  'OFF_par': os.path.join(workdir, 'off.par'),
                  'coffs': os.path.join(workdir, 'coffs'),
                  'thres': 0.15,
                  'npoly': 3,
                  'interact_flag': 0,
                  'logpath': workdir}

    isp.offset_fit(**parameters)

A file `offset_fit.log` containing the output of the command is written in both cases. Any parameters, which should
not be written and need to be set to - in the shell can be omitted in the Python call since all optional parameters
of the functions are already defined with '-' as a default.
The documentation can be called like with any Python function:

.. code-block:: python

    from pyroSAR.gamma.api import isp
    help(isp.offset_fit)

Parser Documentation
********************

.. automodule:: pyroSAR.gamma.parser
    :members:
    :undoc-members:
    :show-inheritance:

API Demo
********

This is a demonstration of an output script as generated automatically by function
:func:`~pyroSAR.gamma.parser.parse_module` for the GAMMA module `ISP`.
Within each function, the command name and all parameters are passed to function
:func:`~pyroSAR.gamma.process`, which converts all input to :py:obj:`str` and then calls the command via the
:mod:`subprocess` module.

.. automodule:: pyroSAR.gamma.parser_demo
    :members:
    :undoc-members:
    :show-inheritance:

Sentinel-1 Tools
================

.. automodule:: pyroSAR.S1
    :members: OSV, removeGRDBorderNoise
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        OSV
        removeGRDBorderNoise

Auxiliary Data Tools
====================

.. automodule:: pyroSAR.auxdata
    :members: dem_autoload, dem_create, get_egm_lookup
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        dem_autoload
        dem_create
        get_egm_lookup

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

    .. autosummary::
        :nosignatures:

        find_datasets
        getargs
        groupby
        groupbyTime
        hasarg
        multilook_factors
        parse_datasetname
        seconds
