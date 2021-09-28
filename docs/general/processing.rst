#################################
SAR Image Handling and Processing
#################################

Image Metadata
==============

Let's start working with our actual satellite data.
At first we load the scene into pyroSAR for analysis of the metadata:

.. code-block:: python

    from pyroSAR import identify
    name = 'S1A_IW_GRDH_1SDV_20150222T170750_20150222T170815_004739_005DD8_3768.zip'
    scene = identify(name)
    print(scene)

This will automatically identify the scene, scan it for metadata and print a summary of selected metadata entries.
Several attribute names (e.g. `sensor` and `acquisition_mode`) are standardized for all SAR scenes.
Further entries, whose names are not standardized, can be found in a dictionary `scene.meta`.
The function :func:`~pyroSAR.drivers.identify` will loop through all SAR images classes (:mod:`pyroSAR.drivers`) and return an
object of the class that was successful in identifying the scene (:class:`~pyroSAR.drivers.SAFE` in this case).

.. _database-handling:

Database Handling
=================

Now that we have made ourselves familiar with the scene, we can import its metadata into an SQLite database using class
:class:`~pyroSAR.drivers.Archive`:

.. code-block:: python

    from pyroSAR import Archive
    dbfile = 'scenes.db'
    with Archive(dbfile) as archive:
        archive.insert(scene)

`dbfile` is a file either containing an already existing database or one to be created.
In this case an SQLite database with SpatiaLite extension is created.
Alternatively, PostgreSQL + PostGIS can be used.

Let's assume our database contains a number of scenes and we want to select some for processing.
We have a shapefile, which contains a geometry delimiting our test site for which we want to
process some Sentinel-1 scenes.
We already processed some scenes in the past and the results are stored in a directory
`outdir`. We only want to select scenes which have not been processed to this directory before.
Furthermore, we are only interested in scenes acquired in Ground Range Detected (GRD) Interferometric Wide
Swath mode (IW), which contain a VV band.

.. code-block:: python

    from spatialist import Vector
    archive = Archive('scenes.db')
    outdir = '/path/to/processed/results'
    maxdate = '20171231T235959'
    with Vector('site.shp') as site:
        selection_proc = archive.select(vectorobject=site,
                                        processdir=outdir,
                                        maxdate=maxdate,
                                        sensor=('S1A', 'S1B'),
                                        product='GRD',
                                        acquisition_mode='IW',
                                        vv=1)
    archive.close()

Here we use the vector geometry driver of package :doc:`spatialist <spatialist:index>`, which is developed alongside of pyroSAR.
The :class:`spatialist.Vector <spatialist.vector.Vector>` object is then passed to method
:meth:`Archive.select <pyroSAR.drivers.Archive.select>`.

.. _processing:

Processing
==========

The returned `selection_proc` is a list of file names for the scenes we selected from the database, which we can now
pass to a processing function:

.. code-block:: python

    from pyroSAR.snap import geocode

    # the target pixel spacing in meters
    spacing = 20

    for scene in selection_proc:
        geocode(infile=scene, outdir=outdir, tr=spacing, scaling='db', shapefile=site)

The function :func:`snap.geocode <pyroSAR.snap.util.geocode>` is a basic utility for SNAP.
It will perform all necessary steps to subset, resample, topographically normalize, geocode and scale the input
image and write GeoTIFF files to the selected output directory.
All necessary files like orbit state vectors and SRTM DEM tiles are downloaded automatically in the background by SNAP.
SNAP is most conveniently used with workflow XMLs. The function geocode parses a workflow for the particular scene,
parametrizes it (depending on the scene type and selected processing parameters) and writes it to the output directory.
It then calls the command `gpt`, which is SNAP's command line interface, on the workflow to execute the processing steps.

