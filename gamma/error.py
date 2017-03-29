##############################################################
# interface for translating GAMMA errors into Python error types
# John Truckenbrodt 2015-2017
##############################################################

import re


def gammaErrorHandler(out, err):
    messages = out.split('\n') if out else []
    messages.extend(err.strip().split('\n'))
    errormessages = [x for x in messages if x.startswith('ERROR')]
    knownErrors = {'image data formats differ': IOError,
                   'cannot open': IOError,
                   'no coverage of SAR image by DEM \(in (?:latitude/northing|longitude/easting)\)': IOError,
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
                   'gc_map operates only with slant range geometry, image geometry in SLC_par: GROUND_RANGE': RuntimeError}
    if len(errormessages) > 0:
        errormessage = errormessages[-1]
        for error in knownErrors:
            if re.search(error, errormessage):
                errortype = knownErrors[error]
                if errortype:
                    raise knownErrors[error](errormessage)
                else:
                    return
        raise GammaUnknownError(errormessage)


class GammaUnknownError(Exception):
    def __init__(self, errormessage):
        Exception.__init__(self, errormessage)
