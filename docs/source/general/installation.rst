############
Installation
############

conda
=====

Starting with version 0.11, pyroSAR is distributed via `conda-forge <https://anaconda.org/conda-forge/pyrosar>`_
and can easily be installed with

::

    conda install --channel conda-forge pyrosar

This is by far the easiest way to work with pyroSAR on any operating system.

pip
===

Installation with pip is also supported and offers the advantage to install intermediate development stages directly
from the GitHub repository. Mind however that several dependencies like GDAL cannot fully be installed this way.
See further below for detailed Linux dependency installation instructions.

Installation of pip (Linux):

::

    sudo apt-get install python-pip

The latest stable release of pyroSAR can then be installed:

::

    python -m pip install pyroSAR

For installation of the latest master branch on GitHub, we need the version control system git. On Windows, git can be
downloaded from `git-scm.com <https://git-scm.com/downloads>`_. On Linux you can install it via command line:

::

    sudo apt-get install git

Once everything is set up, pyroSAR is ready to be installed:

::

    python -m pip install git+https://github.com/johntruckenbrodt/pyroSAR.git

Dependencies
============
The more specific instructions below are intended for Linux users who like to work outside of the Anaconda environment.

GDAL
----
pyroSAR requires GDAL version 2.1 with GEOS and PROJ4 as dependencies as well as the GDAL Python binding.

Ubuntu
++++++
Starting with release Yakkety (16.10), Ubuntu comes with GDAL >2.1.
You can install it like this:

::

    sudo apt-get install python-gdal python3-gdal gdal-bin

For older Ubuntu releases you can add the ubuntugis repository to apt prior to installation to install version >2.1:

::

    sudo add-apt-repository ppa:ubuntugis/ppa
    sudo apt-get update

This way the required dependencies (GEOS and PROJ4 in particular) are also installed.
You can check the version by typing:

::

    gdalinfo --version

Debian
++++++
Starting with Debian 9 (Stretch) GDAL is available in version >2.1 in the official repository.

Building from source
++++++++++++++++++++
Alternatively, you can build GDAL and the dependencies from source. The script `pyroSAR/install/install_deps.sh`
gives specific instructions on how to do it. It is not yet intended to run this script via shell, but rather to
follow the instructions step by step.

SQLite + SpatiaLite
-------------------
While `sqlite3` and its Python binding are usually already installed, the `spatialite` extension needs to be
added. Two packages exist, `libspatialite` and `mod_spatialite`. Both can be used by pyroSAR.
On Ubuntu, `mod_spatialite` has been found to be easier to setup with `sqlite` and can be installed via `apt`:

::

    sudo apt-get install libsqlite3-mod-spatialite

On CentOS, `libspatialite` including shared objects for extension loading can be installed via `yum`:

::

    sudo yum install libspatialite-devel

The following can be run in Python to test the needed functionality:

.. code-block:: python

    import sqlite3

    # setup an in-memory database
    con=sqlite3.connect(':memory:')

    # enable loading extensions and load spatialite
    con.enable_load_extension(True)
    try:
        con.load_extension('mod_spatialite.so')
    except sqlite3.OperationalError:
        con.load_extension('libspatialite.so')

In case loading extensions is not permitted you might need to install the package `pysqlite2`
together with a static build of `sqlite3`. See the script `pyroSAR/install/install_deps.sh` for instructions.
There you can also find instructions on how to install `spatialite` from source.
To test `pysqlite2` you can import it as follows and then run the test above:

.. code-block:: python

    from pysqlite2 import dbapi2 as sqlite3

Installing this package is likely to cause problems with the `sqlite3` library installed on the system.
Thus, it is safer to build a static `sqlite3` library for it (see installation script).

GAMMA
-----
GAMMA's home directory as environment variable 'GAMMA_HOME' is expected to end either as GAMMA_SOFTWARE-<VERSIONNUMBER> or GAMMA_SOFTWARE/<VERSIONNUMBER>. 
If this differs in your install and cannot be changed, a workaround is adjusting the expected pattern in :class:`~pyroSAR.examine.ExamineGamma`.
