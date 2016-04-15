##############################################################
# Interferogram Flattening
# module of software gammaGUI
# Axel Schmidt, John Truckenbrodt 2015
##############################################################

"""
README

    # SYNOPSIS: ph_slope_base <int_in> <SLC-1.par> <OFF_par> <base> <int_out> [int_type} [inverse]

    ph_slope_base takes the following arguments:

    * <int_in>   - (input) interferogram (FCOMPLEX) or unwrapped phase (FLOAT) (unflattened)
    * <SLC_par>  - (input) ISP parameter file for the reference SLC
    * <OFF_par>  - (input)ISP offset/interferogram parameter file
    * <base>     - (input) baseline file
    * <int_out>  - (output) interferogram (fcomplex) or unwrapped interferogram (float) with phase
               trend subtracted/added
    * <int_type> - interferogram type: 0=unwwrapped phase, 1=complex interf. (default=1)
    * <inverse>  - subtract/add inversion flag (0=subtract phase ramp, 1=add phase ramp (default=0))

    Script description:

        The Script batch processes the flattening of the interferogram (ph_slope_base) for the
        SAR files given by the user.
        Thereby the files to process are given by the file 'files_ph_slope_base'. The file
        contains the paths and filenames separated by tabs. The function 'iterate_tabfile'
        iterates through the elements in the file. One iteration is referred to one line in
        that file. This specific file is expected to be created by the GUI.

FLATTENING

    ph_slope_base reads the input complex (or real valued unwrapped) interferogram and subtracts/adds
    the flat-Earth phase trend calculated from the initial baseline (first two lines of provided \
    baseline file).

    The most common use is to remove the flat-Earth phase ramp of a complex valued interferogram (*.int),
    a process often named interferogram "flattening", and write out the flattened complex valued interferogram (*flt).

    The subtracted/added phase ramp is calculated by determining the component of the baseline parallel
    to the look vector for each point across the range swath.

    Removal of the phase ramp is necessary for accurate estimation of the correlation function, and phase
    unwrapping. This is especially true for larger interferometric baselines with a large range phase trend.
    This is because high fringe rates violate the assumption that the phase is constant over the region used
    to estimate the correlation.

    For the baseline a time offset (in seconds) can be specified. The time of a specific image time
    corresponds to the time as specified by the SLC and ISP/offset interferogram parameter files\ plus
    the indicated time offset. (Gamma Remote Sensing (2003): ISP Reference Manual)
    """

import re
import os.path

from ancillary import finder, run

# define (and create) directories for processing results and logfile
path_log = os.path.join(os.getcwd(), "LOG/ISP/")
path_out = os.path.join(os.getcwd(), "ISP/")
for path in [path_log, path_out]:
    if not os.path.exists(path):
        os.mkdir(path)

list_int = finder(path_out, ["*_int"])

if len(list_int) > 0:
    print "#############################################"
    print "interferogram flattening started..."

    for name_int in list_int:
        slc_par = finder(os.getcwd(), [re.findall("[A-Z0-9_]{9,10}[0-9T]{15}_[HV]{2}_slc(?:_cal|)", name_int)[0]])[0]+".par"
        name_off = name_int[:-3]+"off"
        if not os.path.isfile(name_off):
            raise IOError("offset file missing")

        if os.path.isfile(name_int[:-3]+"base_refine"):
            name_base = name_int[:-3]+"base_refine"
        elif os.path.isfile(name_int[:-3]+"base_init"):
            name_base = name_int[:-3]+"base_init"
        else:
            raise IOError("baseline file missing")
        name_flt = name_int[:-3]+"flt"
        print os.path.basename(name_flt)
        run(["ph_slope_base", name_int, slc_par, name_off, name_base, name_flt], path_out, path_log)

    print "...done"
else:
    print "#############################################"
    print "no interferograms found"
