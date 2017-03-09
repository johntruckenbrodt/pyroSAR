##############################################################
# interface for translating GAMMA errors into Python error types
# John Truckenbrodt 2015-2017
##############################################################

import re


def gammaErrorHandler(message):
    errormessage = message.strip().replace('ERROR: ', '')
    knownErrors = {'image data formats differ': IOError,
                   'cannot open': IOError,
                   'no coverage of SAR image by DEM \(in (?:latitude/northing|longitude/easting)\)': IOError,
                   'libgdal.so.1: no version information available': None,
                   'line outside of image': ValueError,
                   'no offsets found above SNR threshold': ValueError,
                   'window size < 4': ValueError,
                   'MLI oversampling factor must be 1, 2, 4, 8': ValueError,
                   'no points available for determining average intensity': ValueError,
                   'no overlap with lookup table': RuntimeError,
                   'nominal terrain altitude (m):      0.000': None,
                   'insufficient offset points to determine offset model parameters': RuntimeError,
                   'insufficient offset points left after culling to determine offset model parameters': RuntimeError,
                   'calloc_1d: number of elements <= 0': ValueError,
                   'multi-look output line:': RuntimeError,
                   'no OPOD state vector found with the required start time!': RuntimeError,
                   'gc_map operates only with slant range geometry, image geometry in SLC_par: GROUND_RANGE': RuntimeError}
    if errormessage != '':
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
