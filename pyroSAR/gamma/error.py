###############################################################################
# interface for translating GAMMA errors messages into Python error types

# Copyright (c) 2015-2019, the pyroSAR Developers.

# This file is part of the pyroSAR Project. It is subject to the
# license terms in the LICENSE.txt file found in the top-level
# directory of this distribution and at
# https://github.com/johntruckenbrodt/pyroSAR/blob/master/LICENSE.txt.
# No part of the pyroSAR project, including this file, may be
# copied, modified, propagated, or distributed except according
# to the terms contained in the LICENSE.txt file.
###############################################################################

import re


def gammaErrorHandler(out, err):
    """
    Function to raise errors in Python. This function is not intended for direct use, but as part of function gamma.util.process
    Args:
        out: the stdout message returned by a subprocess call of a gamma command
        err: the stderr message returned by a subprocess call of a gamma command

    Raises: IOError | ValueError | RuntimeError | None

    """
    
    # scan stdout and stdin messages for lines starting with 'ERROR'
    messages = out.split('\n') if out else []
    messages.extend(err.strip().split('\n'))
    errormessages = [x for x in messages if x.startswith('ERROR')]
    
    # registry of known gamma error messages and corresponding Python error types
    # do not change the Python error types of specific messages! This will change the behavior of several functions
    # in case no error is to be thrown define None as error type
    knownErrors = {'image data formats differ': IOError,
                   'cannot open': IOError,
                   r'no coverage of SAR image by DEM(?: \(in (?:latitude/northing|longitude/easting)\)|)': RuntimeError,
                   'libgdal.so.1: no version information available': None,
                   'line outside of image': ValueError,
                   'no offsets found above SNR threshold': ValueError,
                   'window size < 4': ValueError,
                   'MLI oversampling factor must be 1, 2, 4, 8': ValueError,
                   'no points available for determining average intensity': ValueError,
                   'p_interp(): time outside of range': RuntimeError,
                   'no overlap with lookup table': RuntimeError,
                   'insufficient offset points to determine offset model parameters': RuntimeError,
                   'insufficient offset points left after culling to determine offset model parameters': RuntimeError,
                   'calloc_1d: number of elements <= 0': ValueError,
                   'multi-look output line:': RuntimeError,
                   'no OPOD state vector found with the required start time!': RuntimeError,
                   'gc_map operates only with slant range geometry, image geometry in SLC_par: GROUND_RANGE': RuntimeError,
                   'OPOD state vector data ends before start of the state vector time window': RuntimeError,
                   'non-zero exit status': RuntimeError,
                   'unsupported DEM projection': RuntimeError,
                   'tiffWriteProc:No space left on device': RuntimeError}
    
    # check if the error message is known and throw the mapped error from knownErrors accordingly.
    # Otherwise throw an GammaUnknownError.
    # The actual message is passed to the error and thus visible for backtracing
    if len(errormessages) > 0:
        errormessage = errormessages[-1]
        err_out = '\n\n'.join([re.sub('ERROR[: ]*', '', x) for x in errormessages])
        for error in knownErrors:
            if re.search(error, errormessage):
                errortype = knownErrors[error]
                if errortype:
                    raise errortype(err_out)
                else:
                    return
        raise GammaUnknownError(err_out)


class GammaUnknownError(Exception):
    """
    This is a general error, which is raised if the error message is not yet integrated
    into the known errors of function gammaErrorHandler.
    If this error occurs the message should be included in function gammaErrorHandler.
    """
    
    def __init__(self, errormessage):
        Exception.__init__(self, errormessage)
