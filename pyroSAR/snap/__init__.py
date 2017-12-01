from .util import geocode
from .auxil import getAuxdata, gpt, check_executable
import warnings

SNAP_EXECUTABLE = ['snap32.exe', 'snap64.exe', 'snap.exe', 'snap']

if check_executable(SNAP_EXECUTABLE):
    pass
else:
    warnings.warn("Lorem ipsum dolor sit amet, consectetuer adipiscing elit: http://step.esa.int/main/download/", 
                  Warning)
