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