####################################
Handling of Orbit State Vector Files
####################################
SAR products require additional orbit state vector (OSV) information to improve their spatial location accuracy.
This information is found in externally hosted files, which need to be downloaded separately and are then used by SAR
processing software to update the product's metadata. Currently, pyroSAR only supports handling of Sentinel-1 OSV files.

In SNAP, the corresponding processing node is called `Apply-Orbit-File`, which automatically downloads the OSV file and
updates the scene's metadata. The files are stored in SNAP's location for auxiliary data,
which per default is `$HOME/.snap/auxdata/Orbits`.

In GAMMA, on the other hand, the downloading has to be done manually after which the command `isp.S1_OPOD_vec` can be
used for updating the metadata. pyroSAR offers several approaches for automatically downloading these
files. The central tool for managing existing files and downloading new ones is the class :class:`pyroSAR.S1.OSV`, which
is used for all approaches.

.. note::

    in the following a dedicated directory is defined into which the files will be downloaded. If this directory is
    not defined (default is `None`), the files will be downloaded to SNAP's auxiliary data location (see above). This is
    recommended as the files are kept in a central location that is accessible both by SNAP and by pyroSAR's GAMMA
    functionality.

approach 1: direct download by time span
========================================

In case a large number of scenes is to be processed and/or no internet access is available during processing, the files
can be downloaded by time span to a central directory. This is the most basic approach using the central class
:class:`~pyroSAR.S1.OSV` mentioned above, making use of its methods :meth:`~pyroSAR.S1.OSV.catch` and
:meth:`~pyroSAR.S1.OSV.retrieve`.

.. code-block:: python

    from pyroSAR.S1 import OSV

    osvdir = '/path/to/osvdir'

    with OSV(osvdir) as osv:
        files = osv.catch(sensor='S1A', osvtype='POE',
                          start='20170101T000000', stop='20180101T000000',
                          url_option=1)
        osv.retrieve(files)

Two sub-directories `POEORB` and `RESORB` will be created in `osvdir` containing the downloaded files. `POEORB` will
contain the `Precise Orbit Ephemerides` files, which are the most accurate but are first available about two weeks after
the scene's acquisition. `RESORB` describes the `Restituted Orbit` files, which are less accurate but available
directly after acquisition. See method :meth:`~pyroSAR.S1.OSV.catch` for download URL options.

approach 2: manual download per scene
=====================================

The method :meth:`pyroSAR.drivers.SAFE.getOSV` can be used to directly retrieve the files relevant for the scene.
This method internally uses the methods described above with a time span limited to that of the scene acquisition.

.. code-block:: python

    from pyroSAR import identify
    scene = 'S1A_IW_GRDH_1SDV_20180101T170648_20180101T170713_019964_021FFD_DA78.zip'
    id = identify(scene)
    match = id.getOSV(osvdir='/path/to/osvdir', osvType='POE', returnMatch=True)
    print(match)

approach 3: direct download and scene metadata update (GAMMA only)
==================================================================

The convenience function :func:`pyroSAR.gamma.correctOSV` internally makes use of approach 2 and additionally directly
executes the GAMMA command `isp.S1_OPOD_vec` for updating the scene's metadata with the information of the OSV file.
The scene has to be unpacked first (see :meth:`pyroSAR.drivers.SAFE.unpack`).

.. code-block:: python

    from pyroSAR import identify
    from pyroSAR.gamma import correctOSV
    scene = 'S1A_IW_GRDH_1SDV_20180101T170648_20180101T170713_019964_021FFD_DA78.zip'
    id = identify(scene)
    id.unpack('tmpdir')
    correctOSV(id=id, osvdir='/path/to/osvdir', osvType='POE')

approach 4: automatic download and use during processing
========================================================

The processing function :func:`pyroSAR.gamma.geocode` automatically downloads OSV files needed for processing and
updates the scene's metadata using function :func:`~pyroSAR.gamma.correctOSV`.
It is thus the most convenient way to handle these files and related processing steps.
The parameter `allow_RES_OSV` can be used to allow processing with `RES` files if no `POE` file is available yet.

.. code-block:: python

    from pyroSAR.gamma import geocode
    scene = 'S1A_IW_GRDH_1SDV_20180101T170648_20180101T170713_019964_021FFD_DA78.zip'
    geocode(scene=scene,
            dem='/path/to/demfile',
            tmpdir='tmpdir',
            outdir='outdir',
            targetres=20,
            osvdir='/path/to/osvdir',
            allow_RES_OSV=False)

Similarly, the function :func:`pyroSAR.snap.util.geocode` also automatically downloads OSV files and chooses the best
matching OSV type for processing.

.. code-block:: python

    from pyroSAR.snap import geocode
    scene = 'S1A_IW_GRDH_1SDV_20180101T170648_20180101T170713_019964_021FFD_DA78.zip'
    geocode(infile=scene,
            outdir='outdir',
            allow_RES_OSV=True)

In contrast to the GAMMA function, the OSV download directory cannot be set because of the fixed SNAP auxiliary data
location. The type of the available OSV file is written to the workflow XML file for processing:

.. code-block:: xml

    <node id="Apply-Orbit-File">
        <operator>Apply-Orbit-File</operator>
        <sources>
            <sourceProduct refid="Read"/>
        </sources>
        <parameters class="com.bc.ceres.binding.dom.XppDomElement">
            <orbitType>Sentinel Restituted (Auto Download)</orbitType>
            <polyDegree>3</polyDegree>
            <continueOnFail>false</continueOnFail>
        </parameters>
    </node>
