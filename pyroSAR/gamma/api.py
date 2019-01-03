import os
import sys
import warnings

from .parser import autoparse

try:
    autoparse()
    
    sys.path.insert(0, os.path.join(os.path.expanduser('~'), '.pyrosar'))
    
    try:
        from gammaparse import *
    except ImportError:
        warnings.warn('found a Gamma installation directory, but module parsing failed')

except RuntimeError:
    warnings.warn('could not find Gamma installation directory; please set the GAMMA_HOME environment variable')
