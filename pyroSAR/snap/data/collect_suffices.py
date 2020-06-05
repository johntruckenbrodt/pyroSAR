##############################################################
# SNAP source code scan for retrieving operator suffices

# Copyright (c) 2020, the pyroSAR Developers.

# This file is part of the pyroSAR Project. It is subject to the
# license terms in the LICENSE.txt file found in the top-level
# directory of this distribution and at
# https://github.com/johntruckenbrodt/pyroSAR/blob/master/LICENSE.txt.
# No part of the pyroSAR project, including this file, may be
# copied, modified, propagated, or distributed except according
# to the terms contained in the LICENSE.txt file.
##############################################################
import os
import re
import subprocess as sp
from spatialist.ancillary import finder

"""
This script clones the SNAP source code from GitHub and reads the suffices for SNAP operators.
E.g. The operator Terrain-Flattening has a suffix TF. If Terrain-Flattening is added to a workflow
in SNAP's graph builder, this suffix is appended to the automatically created output file name.
As pyroSAR also automatically creates file names with processing step suffices, it is convenient to just
use those defined by SNAP.
Currently I am not aware of any way to retrieve them directly from a SNAP installation.
A question in the STEP forum is asked: https://forum.step.esa.int/t/snappy-get-operator-product-suffix/22885

Feel free to contact me if you have ideas on how to improve this!
"""


def main():
    # some arbitrary directory for the source code
    workdir = os.path.join(os.path.expanduser('~'), '.pyrosar', 'snap_code')
    os.makedirs(workdir, exist_ok=True)
    
    # the name of the Java properties file containing the operator-suffix lookup
    outfile = 'snap.suffices.properties'
    
    # clone all relevant toolboxes
    for tbx in ['snap-engine', 'snap-desktop', 's1tbx']:
        print(tbx)
        target = os.path.join(workdir, tbx)
        if not os.path.isdir(target):
            url = 'https://github.com/senbox-org/{}'.format(tbx)
            sp.check_call(['git', 'clone', url], cwd=workdir)
        else:
            sp.check_call(['git', 'pull'], cwd=target)
    
    # search patterns for relevant files
    # Usually files containing operator classes are named <operator>Op.java but with out dashes
    # e.g. TerrainFlatteningOp.java for the Terrain-Flattening operator
    # One exception is Calibration for which there is a sub-class for each SAR sensor
    operators = finder(workdir, ['*Op.java', 'BaseCalibrator.java'])
    
    # a list for collection the suffices
    collect = []
    
    for op in operators:
        with open(op) as infile:
            content = infile.read()
        
        # the suffix is defined as a class attribute PRODUCT_SUFFIX
        pattern = 'String PRODUCT_SUFFIX = \"_([a-zA-Z]*)\"'
        match = re.search(pattern, content)
        if match:
            suffix = match.groups()[0]
        else:
            suffix = ''
        
        # the name of the operator as available in the UI
        pattern = 'alias = \"([a-zA-Z-]*)\"'
        match = re.search(pattern, content)
        if match:
            alias = match.groups()[0]
        else:
            alias = None
        
        # only collect operators for which an alias exists, i.e. which are exposed in the UI,
        # and for which a suffix is defined. In the UI, all operators for which no suffix exists
        # will just get no suffix in any written file.
        if alias is not None and suffix != '':
            print(alias, suffix)
            collect.append('{0}={1}'.format(alias, suffix))
    
    print('found {} matching operators'.format(len(collect)))
    
    with open(outfile, 'w') as out:
        out.write('\n'.join(sorted(collect, key=str.lower)))


if __name__ == '__main__':
    main()
