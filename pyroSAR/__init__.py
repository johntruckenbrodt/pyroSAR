from .drivers import *
from .archive import Archive, drop_archive
from . import ancillary, drivers

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version(__name__)
except PackageNotFoundError:
    # package is not installed
    pass
