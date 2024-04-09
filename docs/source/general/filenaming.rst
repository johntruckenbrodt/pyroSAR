###########
File Naming
###########

pyroSAR internally uses a fixed naming scheme to keep track of processed results. For each scene an identifier is created,
which contains the sensor, acquisition mode, orbit (ascending or descending) and the time stamp of the acquisition start.
For example `S1A__IW___A_20150222T170750`, which is created by calling method :meth:`~pyroSAR.drivers.ID.outname_base`:

.. code-block:: python

    from pyroSAR import identify
    id = identify('S1A_IW_GRDH_1SDV_20150222T170750_20150222T170815_004739_005DD8_3768.zip')
    print(id.outname_base())

For each attribute a fixed number of digits is reserved. In case the attribute is shorter than this number,
the rest of the digits is filled with underscores. I.e., the sensor field is four digits long, but 'S1A' only three.
Thus, `S1A_` is the sensor slot. In the same way, `IW__` is the acquisition mode slot, which is also four digits long.
`A` denotes ascending orbit, the time stamp is in format YYYYmmddTHHMMSS.

Processing functions like :func:`~pyroSAR.gamma.geocode` add suffixes to this identifier to further keep track of
individual processing steps performed on the dataset.
This core concept is used by many pyroSAR functions internally to keep track of which scenes have been processed before.
