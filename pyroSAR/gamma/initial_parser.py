#!/usr/bin/env python3

"""
Creates python help functions for GAMMA programs
helpers are stored in {sys.prefix}/share/pyrosar/gamma/{gammaVersion}/parse

""" 

import sys
import os
import importlib

from  pyroSAR.gamma import parser , api

def main():
    # check needed environment 

    neededEnvVars = ('GAMMA_HOME','GAMMA_RASTER','LD_LIBRARY_PATH') 
    for actVar in neededEnvVars :
        if actVar not in os.environ :
            raise ValueError(f"needed environment variable {actVar} not set")
    
    parser.initialparse()

    # check if parsing found GAMMA programs

    importlib.reload(sys.modules['pyroSAR.gamma.api'])

    if hasattr(api,'diff') :
        print(f"GAMMA helper successfully installed")

if __name__ == '__main__':
    main()

