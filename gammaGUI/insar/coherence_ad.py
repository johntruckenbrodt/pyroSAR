##############################################################
# Coherence Estimation using phase slope and texture
# module of software gammaGUI
# Axel Schmidt, John Truckenbrodt 2014-2015
##############################################################

"""
README

# SYNOPSIS: cc_ad <interf> <pwr1> <pwr2> <slope> <texture> <cc_ad> <width> [box_min] [box_max] [wgt_flag] [loff] [nl]

cc_ad takes the following arguments:

* <interf>   - (input) complex interferogram
* <pwr1>     - (input) intensity image of the first scene
* <pwr2>     - (input) intensity image of the second scene
* <slope>    - (input) phase slope data (or -)
    * <texture>  - (input) backscatter texture data (or -)
    * <cc_ad>    - (output) adaptive coherence estimate file
    * <width>    - number of samples per row
    * [box_min]  - smallest correlation average box size (default = 3.0)
    * [box_max]  - largest correlation average box size (default = 9.0)
    * [wgt_flag] - flag for weighting function selection (default = 0):
                    0: constant weights
                    1: gaussian weighting function
    * [loff]     - offset in lines to starting line (default = 0)
    * [nl]       - number of lines to process (default = 0, to end of File)

    Script description:

        The Script batch processes the coherence (cc_wave) for the SAR files given by the user.

        NOTE: For cc_wave the arguments [xmin], [xmax], [ymin] and <ymax> are omitted in the
              script. For cc_ad the arguments [loff] and [nl] are omitted in the script.

        The coherence images are stored in a unique folder named by the images for which the
        coherence was estimated.
        Furthermore a logfile is created in the home directory containing the processed products
        and the filepaths where they were stored.

        The output filenames for coherence are assembled from the filenames of the two input
        (.mli, .rmli) files and the extensions '_cc' and 'cc_ad' are added, respectively.

        ANNOTATIONS:
        (*1) 'The size of the estimation window is a crucial factor determining the coherence
              estimate. For increasing window size the estimation bias and the estimation un-
              certainty decrease while the spatial resolution of the coherence image decreases.'
              (Gamma Remote Sensing (2007): Docu-InSAR Processing)

            Attention: The resulting window size also dependends on the previously used
                   multilooking factors, i.e. the size of the coherence estimation
                   window size results from the multilooking factors times the factors
                   defined for coherence window size within this script

        (*2) 'To decrease the effect of resolution loss due to the windowing operation, weight-
              ing functions can be applied within the window for CC_WAVE (0: constant (default) |
              1: triangular | 2: gaussian | 3: none (phase only) and CC_AD (0: constant weights |
              1: gaussian weighting function), respectively.

        (*3) 'The user can select the minimum and maximum window size' in order to compute co-
              herence 'using an adaptive window size [that] has the advantage to decrease the
              estimation bias and improve the accuracy of the estimate, while attempting to pre-
              serve spatial resolution'. The coherence image obtained with the adaptive estimation
              method presents a stronger contrast due to the fact that low coherence areas have
              smaller bias. The level of noise has also decreased because of the smaller uncertain-
              ty of the coherence estimates, in particular in case of low coherence (Gamma Remote
              Sensing (2013): LAT Users Guide)
"""

import os
import re

from ancillary import finder, run, ReadPar
from gamma.util import ISPPar


# define (and create) directories for processing results and logfile
path_log = os.path.join(os.getcwd(), "LOG/ISP/")
path_out = os.path.join(os.getcwd(), "ISP/")
for path in [path_log, path_out]:
    if not os.path.exists(path):
        os.makedirs(path)

par = ReadPar(os.path.join(os.getcwd(), "PAR", "coherence_ad.par"), type="exe")

# retrieve additional arguments from script call
differential = True if par.differential == "True" else False

# find flattened interferograms
list_flt = finder(path_out, ["*_int_diff"]) if differential else finder(path_out, ["*_flt"])

if len(list_flt) > 0:
    print "#############################################"
    print "estimation started..."

    for name_flt in list_flt:
        # extract timestamps from flt name
        id_pwr = re.findall("[A-Z0-9_]{9,10}[0-9T]{15}_[HV]{2}_slc(?:_cal|)", name_flt)
        # find mli/rmli images matching the extracted timestamps
        try:
            name_mli = finder(os.getcwd(), ["*"+id_pwr[0]+"_mli"])[0]
            # name_rmli = finder(os.getcwd(), ["*"+id_pwr[1]+"_reg_mli"])[0]
            name_rmli = name_flt[:-8]+"reg_mli" if differential else name_flt[:-3]+"reg_mli"
        except:
            raise IOError("multilooked images missing")

        # concatenate coherence image name
        name_cc = name_flt+"_cc_ad"

        if not os.path.isfile(name_cc):
            print os.path.basename(name_cc)
            # read image samples
            samples = str(ISPPar(name_mli + ".par").range_samples)
            # run gamma command
            run(['cc_ad', name_flt, name_mli, name_rmli, "-", "-", name_cc, samples, par.box_min, par.box_max, par.wgt_ad], path_out, path_log)
        else:
            print "coherence image", name_cc, "already exists"

    print "...estimation finished"
else:
    print "#############################################"
    print "no {0} interferograms found".format("differential" if differential else "flattened")
