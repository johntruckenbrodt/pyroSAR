#######
Logging
#######

pyroSAR makes use of the :mod:`logging` module to display status messages for running processes.
See `Logging HOWTO <https://docs.python.org/3/howto/logging.html>`_ for a basic tutorial.
To display log messages you may add one of the following examples to your script:

.. code-block:: python

  import logging

  # basic info
  logging.basicConfig(level=logging.INFO)

  # basic info with some message filtering
  logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

  # detailed debug info
  logging.basicConfig(level=logging.DEBUG)
