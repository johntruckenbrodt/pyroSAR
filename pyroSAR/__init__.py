from .drivers import *
from ._dev_config import ConfigHandler

ConfigHandler = ConfigHandler()

from . import ancillary, drivers

from pkg_resources import get_distribution, DistributionNotFound
try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    pass
