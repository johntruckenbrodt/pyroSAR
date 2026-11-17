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
from typing import Self

from pyroSAR.examine import ExamineGamma


# These placeholders allow pyroSAR.gamma modules to be imported when GAMMA
# is unavailable. They are replaced by the generated modules when parsing
# succeeds.
class _UnavailableGammaModule:
    def __init__(self: Self, module: str) -> None:
        self.module = module
    
    def __getattr__(self: Self, command: str):
        raise AttributeError(
            f"The command '{command}' is not available. "
            f"Please install GAMMA module '{self.module.upper()}'."
        )


diff = _UnavailableGammaModule('diff')
disp = _UnavailableGammaModule('disp')
isp = _UnavailableGammaModule('isp')
lat = _UnavailableGammaModule('lat')

if 'GAMMA_HOME' in os.environ:
    
    config = ExamineGamma()
    gammaVersion = config.version
    
    refDir = f"{sys.prefix}/share/pyrosar/gamma/{gammaVersion}"
    if not os.path.isdir(f"{refDir}/parse") : 
        log.warning(f"Environment variable GAMMA_HOME is set but gamma helper dir {refDir}/parse not found. Parsing after GAMMA installation not correct")
    else :

        sys.path.insert(0,refDir)

        from  parse import *
