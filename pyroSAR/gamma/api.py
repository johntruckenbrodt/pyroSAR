###############################################################################
# import wrapper for the pyroSAR GAMMA API

# Copyright (c) 2018-2026, the pyroSAR Developers.

# This file is part of the pyroSAR Project. It is subject to the
# license terms in the LICENSE.txt file found in the top-level
# directory of this distribution and at
# https://github.com/johntruckenbrodt/pyroSAR/blob/master/LICENSE.txt.
# No part of the pyroSAR project, including this file, may be
# copied, modified, propagated, or distributed except according
# to the terms contained in the LICENSE.txt file.
###############################################################################
import os
import sys

from .parser import autoparse


# These placeholders allow pyroSAR.gamma modules to be imported when GAMMA
# is unavailable. They are replaced by the generated modules when parsing
# succeeds.
class _UnavailableGammaModule:
    def __init__(self, module: str) -> None:
        self.module = module
    
    def __getattr__(self, command: str):
        raise AttributeError(
            f"The command '{command}' is not available. "
            f"Please install GAMMA module '{self.module.upper()}'."
        )


diff = _UnavailableGammaModule('diff')
disp = _UnavailableGammaModule('disp')
isp = _UnavailableGammaModule('isp')
lat = _UnavailableGammaModule('lat')

if 'GAMMA_HOME' in os.environ:
    autoparse()
    
    sys.path.insert(0, os.path.join(os.path.expanduser('~'), '.pyrosar'))
    
    from gammaparse import *
