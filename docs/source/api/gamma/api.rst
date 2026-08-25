.. _gamma-command-api:

GAMMA Command API
-----------------

This is an attempt to make it easier to execute GAMMA commands by offering automatically parsed Python functions.
Thus, instead of executing the command via the shell:

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

The parameter names found in the docstring of the command are now function arguments, which makes writing scripts much more user-friendly.
Additionally, three arguments are added to each command: ``logpath``, ``outdir`` and ``shellscript``.

``logpath`` defines the folder where logfiles will be stored.
The file names follow the pattern ``<command name>.log``.
Hence, in both examples above a file ``offset_fit.log`` containing the output of the command call is written.
The outputs of repeated calls are written to the same file and are delimited by a multiple-hashtag line.

``outdir`` controls the folder where the command is executed.
This might be helpful if users want to define input and output file names relative to a certain folder.
Moreover, some GAMMA commands output files with fixed or auto-generated names into the execution folder.

With ``shellscript``, a user can define the name of a bash script in which all executed GAMMA commands are stored.
Two environment variables are defined for better readability: ``$GAMMA_HOME`` and ``$OUTDIR``.
The latter corresponds to the user-defined ``outdir`` and is redefined in the script every time it is changed by the user.
Hence, setting ``outdir`` may increase the script's readability.

.. caution::

    While this script is generally executable, it is no full equivalent of what is done in pyroSAR.
    Some steps are done directly in Python and there is thus no GAMMA equivalent.
    The script is currently more intended for increased transparency on what pyroSAR does and to facilitate bug fixing and knowledge exchange.

Any parameters, which should not be written and need to be set to ``-`` in the shell, can be omitted in the Python call since all optional parameters of the functions are already defined with ``'-'`` as a default.
The documentation can be called like with any Python function:

.. code-block:: python

    from pyroSAR.gamma.api import isp
    help(isp.offset_fit)

Common pitfalls
***************

The parser relies on a correct GAMMA installation so that all commands can be executed and their docstrings extracted.
It is recommended to turn on logging on first import of pyroSAR's gamma functionality to obtain parsing details.
Generally, the parser does not raise errors when parsing of individual commands fails.
Hence, reading through the logs is recommended for checking whether all commands were parsed successfully.
Arguably, this is can be improved.
For re-parsing modules, delete them from their storage location under ``$HOME/.pyrosar``.

In general, the parser is kept up to date with the latest GAMMA version.
General support for older versions is too much effort.
Exceptions can be discussed. Feel free to reach out.

Parser API Documentation
************************

.. automodule:: pyroSAR.gamma.parser
    :members:
    :undoc-members:
    :show-inheritance:

API Demo
********

This is a demonstration of an output script as generated automatically by function
:func:`~pyroSAR.gamma.parser.parse_module` for the GAMMA module `ISP`.
Within each function, the command name and all parameters are passed to function
:func:`~pyroSAR.gamma.auxil.process`, which converts all input to :py:obj:`str` and then calls the command via the
:mod:`subprocess` module.

.. automodule:: pyroSAR.gamma.parser_demo
    :members:
    :undoc-members:
    :show-inheritance: