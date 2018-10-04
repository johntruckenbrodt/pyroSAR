import os
import sys
import warnings

from .parser import autoparse

try:
    autoparse()

except RuntimeError:
    warnings.warn('could not find Gamma installation; module autoparsing failed')

sys.path.insert(0, os.path.join(os.path.expanduser('~'), '.pyrosar'))

try:
    from gammaparse import *
except ImportError:
    pass
