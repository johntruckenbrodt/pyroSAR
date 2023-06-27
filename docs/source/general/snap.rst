########
SNAP API
########

pyroSAR offers a collection of tools to parse SNAP XML workflows and execute them with SNAP's Graph Processing Tool
(`GPT <https://senbox.atlassian.net/wiki/spaces/SNAP/pages/70503475/Bulk+Processing+with+GPT>`_). All functionality is
purely performed in Python and only the command line calls to GPT interact with SNAP. SNAP's Python API
`snappy <https://senbox.atlassian.net/wiki/spaces/SNAP/pages/19300362/How+to+use+the+SNAP+API+from+Python>`_ is not used
due to installation limitations and processing performance.

The following serves as a minimal example to showcase the core API functionality. A more complex example is given with
function :func:`pyroSAR.snap.util.geocode`.

.. code-block:: python

    from pyroSAR.snap.auxil import parse_recipe, parse_node

    workflow = parse_recipe('blank')

    read = parse_node('Read')
    read.parameters['file'] = 'S1A_IW_GRDH_1SDV_20150222T170750_20150222T170815_004739_005DD8_3768.zip'
    read.parameters['formatName'] = 'SENTINEL-1'
    workflow.insert_node(read)

    tnr = parse_node('ThermalNoiseRemoval')
    workflow.insert_node(tnr, before=read.id)

    bnr = parse_node('Remove-GRD-Border-Noise')
    bnr.parameters['selectedPolarisations'] = ['VV']
    workflow.insert_node(bnr, before=tnr.id)

    write = parse_node('Write')
    write.parameters['file'] = 'outname'
    write.parameters['formatName'] = 'BEAM-DIMAP'
    workflow.insert_node(write, before=bnr.id)

    workflow.write('outname_proc')

Here, the function :func:`~pyroSAR.snap.auxil.parse_recipe` is first used to create an empty workflow object of type
:class:`~pyroSAR.snap.auxil.Workflow`.
Using the function :func:`~pyroSAR.snap.auxil.parse_node`, individual processing nodes can be loaded as
:class:`~pyroSAR.snap.auxil.Node` objects and parameterized using a :class:`~pyroSAR.snap.auxil.Par` object via
``<node>.parameters``.
The method :meth:`~pyroSAR.snap.auxil.Workflow.insert_node` is then used to insert the nodes into the workflow including
linking of the nodes by modifying the source node entries. E.g. `Read` is set as source of the newly inserted
`Remove-GRD-Border-Noise` node. As a last step, the workflow is written to an XML file with method
:meth:`~pyroSAR.snap.auxil.Workflow.write`.

This XML file can then be passed to function :func:`~pyroSAR.snap.auxil.gpt` to process the workflow by internally
calling the GPT command line tool:

.. code-block:: python

    from pyroSAR.snap.auxil import gpt

    gpt('outname_proc.xml', tmpdir='.')

workflow splitting
==================

Simple workflows like the one shown above take only a few seconds to process, but the more processing nodes are added,
the more time it obviously takes to execute them. However, it was observed that executing long workflows takes longer
and consumes more memory than executing each node individually. pyroSAR offers functionality to split long workflows
into smaller groups and execute them in sequence with intermediate files being written in a temporary directory.
First, the workflow nodes are grouped to contain a defined number of processing nodes, i.e. everything but `Read` and
`Write`, using function :func:`~pyroSAR.snap.auxil.groupbyWorkers`:

.. code-block:: python

    from pyroSAR.snap.auxil import groupbyWorkers

    groupbyWorkers('outname_proc.xml', n=1)

This will return

.. code-block:: python

    [['Read', 'ThermalNoiseRemoval'], ['Remove-GRD-Border-Noise', 'Write']]

These groups can directly be passed passed to function :func:`~pyroSAR.snap.auxil.gpt` via parameter ``groups``.
Internally the workflow is then split based on the groups and written to new XML files in a temporary directory using
function :func:`~pyroSAR.snap.auxil.split`. In this case, two workflows would be created:

- `Read` -> `ThermalNoiseRemoval` -> `Write`
- `Read` -> `Remove-GRD-Border-Noise` -> `Write`

These new files are then executed in sequence with intermediate `BEAM-DIMAP`
files written in the same directory as the sub-workflow XML files. After processing this directory is deleted unless
parameter ``cleanup`` of function :func:`~pyroSAR.snap.auxil.gpt` is set to ``False``.

backwards compatibility
=======================

With new versions of SNAP, new parameters are introduced and others removed. If a new parameter is not listed in the
node's XML description its default is used by SNAP during processing. If, however, a parameter is contained in the
workflow that is no longer supported by SNAP, the processing will be terminated. This can easily happen if the workflow
was created by an older version of SNAP. pyroSAR reads the error messages and, if an unknown parameter is mentioned,
deletes this parameter from the workflow, saves it to a new file and executes it instead.

troubleshooting
===============

SNAP as well as pyroSAR's SNAP API are constantly being developed and bugs are unfortunately inevitable.
This section is intended to guide users to better interpret errors and unexpected behaviour.

*The process is running but seems inactive without any progress.*

This might be related to SNAP's inability to download needed DEM tiles.
SNAP will be stuck in a loop infinitely trying to download the missing tiles.
This can be identified by directly running gpt in the command line.
However, by operating gpt through a Python subprocess, it is not possible to see those command line messages.
Only after a process has terminated, all messages can be retrieved and be written to log or error files.

A simple approach to interpret such a behaviour is to first create a workflow XML file with
:func:`~pyroSAR.snap.util.geocode`'s parameter ``test=True`` (so that only the XML is written but it is not executed):

.. code-block:: python

    from pyroSAR.snap import geocode
    geocode(scene='S1A_IW_GRDH_1SDV_20200720T023849_20200720T023914_033532_03E2B5_2952.zip',
            outdir='/test', test=True)


and then run gpt on it directly in the shell (i.e. outside of Python):

::

    gpt /test/S1A__IW___D_20200720T023849_VV_Orb_ML_TC_proc.xml

This way one can directly see gpt's status, which in this case might be

::

    SEVERE: org.esa.snap.core.dataop.dem.ElevationFile: java.lang.reflect.InvocationTargetException