#############
Configuration
#############

pyroSAR stores configuration under `$HOME/.pyrosar`.
It contains a file `config.ini` which stores installation paths of SNAP and GAMMA.
The installations are first identified by running the respective `Examine*` class (e.g. :class:`~pyroSAR.examine.ExamineSnap`):

.. code-block:: python

    from pyroSAR.examine import ExamineSnap
    config = ExamineSnap()

SNAP configuration can also be modified with this class, either by the object properties `userpath` and `auxdatapath` or by the underlying :class:`~pyroSAR.examine.SnapProperties` object:

.. code-block:: python

    config.userpath = '/path/to/snap/data'
    config.snap_properties['snap.userdir'] = '/path/to/snap/data'

The values are directly written to either `snap.auxdata.properties` or `snap.properties` under `$HOME/.snap/etc`.
The content of these files will override that in the files found under `etc` in the SNAP installation folder.
Setting a parameter to `None` will comment out the value in the respective file.
