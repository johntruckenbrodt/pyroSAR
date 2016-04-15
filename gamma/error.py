
import re


def gammaErrorHandler(process):
    errormessage = process.stderr.read().strip().replace("ERROR: ", "")
    knownErrors = {"image data formats differ": IOError,
                   "cannot open": IOError,
                   "no coverage of SAR image by DEM \(in (?:latitude/northing|longitude/easting)\)": IOError,
                   "libgdal.so.1: no version information available": ImportWarning,
                   "line outside of image": ValueError,
                   "no offsets found above SNR threshold": ValueError,
                   "window size < 4": ValueError,
                   "MLI oversampling factor must be 1, 2, 4, 8": ValueError,
                   "no points available for determining average intensity": ValueError,
                   "insufficient offset points to determine offset model parameters": RuntimeError,
                   "insufficient offset points left after culling to determine offset model parameters": RuntimeError,
                   "calloc_1d: number of elements <= 0": ValueError,
                   "multi-look output line:": RuntimeError}
    if errormessage != "":
        for error in knownErrors:
            if re.search(error, errormessage):
                raise knownErrors[error](errormessage)
        raise GammaUnknownError(errormessage)


class GammaUnknownError(Exception):
    def __init__(self, errormessage):
        Exception.__init__(self, errormessage)
