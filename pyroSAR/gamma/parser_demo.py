from pyroSAR.gamma.auxil import process


def adapt_filt(int, sm, width, low_SNR_thr='-', filt_width='-', xmin='-', xmax='-', ymin='-', ymax='-', logpath=None):
    """
    | Adaptive bandpass filtering of interferograms
    | Copyright 2016, Gamma Remote Sensing, v3.5 clw 17-Feb-2016

    Parameters
    ----------
    int:
        (input) complex interferogram image filename
    sm:
        (output) smoothed interferogram filename
    width:
        number of samples/row
    low_SNR_thr:
        low SNR threshold (default = .25);
    filt_width:
        filter width in pixels (default = 1.0)
    xmin:
        offset to starting range pixel(default = 0)
    xmax:
        offset last range pixel (default = width-1)
    ymin:
        offset to starting azimuth row (default = 0)
    ymax:
        offset to last azimuth row (default = nlines-1)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/adapt_filt', int, sm,
             width, low_SNR_thr, filt_width, xmin, xmax, ymin, ymax], logpath=logpath)


def adf(interf, sm, cc, width, alpha='-', nfft='-', cc_win='-', step='-', loff='-', nlines='-', wfrac='-', logpath=None):
    """
    | Adaptive spectral filtering for complex interferograms
    | Copyright 2016, Gamma Remote Sensing, v3.5 15-Feb-2016 clw

    Parameters
    ----------
    interf:
        (input) interferogram (fcomplex)
    sm:
        (output) filtered interferogram (fcomplex)
    cc:
        (output) coherence derived from filtered interferogram (float)
    width:
        number of samples/line
    alpha:
        exponent for non-linear filtering (enter - for default: 0.40)
    nfft:
        filtering FFT window size, 2\*\*N, 8 --> 512, (enter - for default: 32)
    cc_win:
        coherence parameter estimation window size odd, max: 15 (enter - for default: 5)
    step:
        processing step (enter - for default: nfft/8)
    loff:
        offset to starting line to process (enter - for default: 0)
    nlines:
        number of lines to process (enter - for default: to end of file)
    wfrac:
        minimum fraction of points required to be non-zero in the filter window (enter - for default: 0.200)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/adf', interf, sm, cc,
             width, alpha, nfft, cc_win, step, loff, nlines, wfrac], logpath=logpath)


def af_SLC(SLC_par, SLC, rwin='-', azwin='-', dr='-', daz='-', thres='-', a1_flg='-', b0_flg='-', offsets='-', n_ovr='-', roff='-', azoff='-', logpath=None):
    """
    | Focus testing for SLC data using autofocus estimation of effective velocity
    | Copyright 2016, Gamma Remote Sensing, v1.4 16-Feb-2016 clw/uw

    Parameters
    ----------
    SLC_par:
        (input) ISP SLC image parameter file
    SLC:
        (input) single-look complex image
    rwin:
        range window size (enter - for default: 1024)
    azwin:
        azimuth window size (enter - for default: 4096)
    dr:
        range sample increment (enter - for default: 1024,  enter 0 for single patch)
    daz:
        azimuth line increment (enter - for default: 8192,  enter 0 for single patch)
    thres:
        offset estimation SNR threshold (enter - for default: 10.000)
    a1_flg:
        fit a1 for first derivative of the effective velocity w.r.t.range
            * 0: no (default)
            * 1: yes

    b0_flg:
        fit b0 for first derivative of the effective velocity w.r.t. along-track time
            * 0: no (default)
            * 1: yes

    offsets:
        (output) range and azimuth offsets and SNR data in text format, enter - for no output
    n_ovr:
        SLC oversampling factor (1,2,4: enter - for default: 1)
    roff:
        range offset for single patch center
    azoff:
        azimuth offset for single patch center
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/af_SLC', SLC_par, SLC, rwin, azwin,
             dr, daz, thres, a1_flg, b0_flg, offsets, n_ovr, roff, azoff], logpath=logpath)


def ASAR_LO_phase_drift(SLC1_par, SLC2_par, OFF_par, ph_drift, logpath=None):
    """
    | Calculate interferometric phase correction due to drift of the ASAR local oscillator
    | Copyright 2015, Gamma Remote Sensing, v1.1 3-Dec-2015 clw

    Parameters
    ----------
    SLC1_par:
        (input) SLC-1 ISP image parameter file
    SLC2_par:
        (input) SLC-2 ISP image parameter file
    OFF_par:
        (input) ISP offset/interferogram parameter file 
    ph_drift:
        (output) interferometric phase correction due to drift of the ASAR LO (radians)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/ASAR_LO_phase_drift',
             SLC1_par, SLC2_par, OFF_par, ph_drift], logpath=logpath)


def ASAR_XCA(ASA_XCA, antenna, swath='-', pol='-', logpath=None):
    """
    | Interpretation of ASAR external calibration data file (ASA_XCA)
    | Copyright 2006, Gamma Remote Sensing, v1.1 7-June-2006 awi/uw/clw

    Parameters
    ----------
    ASA_XCA:
        (input) ASAR external calibration data file (binary)
    antenna:
        (output) 1-way antenna gain pattern file or '-' (if not provided)
            or 'all' to generate all ASAR antenna diagrams
    swath:
        ASAR swath (IS1,IS2,...IS7;SS1,SS2,...SS5)
    pol:
        polarization (HH,VV,HV,VH)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/ASAR_XCA',
             ASA_XCA, antenna, swath, pol], logpath=logpath)


def ave_image(im_list, width, ave, start='-', nlines='-', pixav_x='-', pixav_y='-', zflag='-', nmin='-', logpath=None):
    """
    | Calculate average of a stack of images (float format)
    | Copyright 2015, Gamma Remote Sensing, v1.9 20-Nov-2015 clw

    Parameters
    ----------
    im_list:
        (input) text file with names of co-registered images in column 1 (float)
    width:
        number of samples/line
    ave:
        (output) average of input image data files (float)
    start:
        starting line (default: 1)
    nlines:
        number of lines to process (enter -  for default: entire file)
    pixav_x:
        number of pixels to average in width  (default: 1)
    pixav_y:
        number of pixels to average in height (default: 1)
    zflag:
        zero flag:
            * 0: interpret 0.0 as missing data value (default)
            * 1: interpret 0.0 as valid data

    nmin:
        minimum number of images required to calculate the average if zflag = 0 (default: 3/4\*nfiles)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/ave_image', im_list,
             width, ave, start, nlines, pixav_x, pixav_y, zflag, nmin], logpath=logpath)


def az_integrate(data, width, azi, cflg, scale='-', lz='-', logpath=None):
    """
    | Calculate azimuth integral of float data (unwrapped phase or azimuth offsets)
    | Copyright 2012, Gamma Remote Sensing, v1.2 6-Feb-2012

    Parameters
    ----------
    data:
        (input) input data (example: SBI dtrapped phase) (float)
    width:
        (input) number of range samples/line 
    azi:
        (output) input data integrated along azimuth (float)
    cflg:
        integration constant flag:
            * 0: set azimuth integral value to 0.0 at specified line
            * 1: set average of the azimuth integral to 0.0

    scale:
        scale factor to apply to the data (enter - for default, default: 1.0)
    lz:
        line offset where the azimuth integral is set to 0.0 (cflg = 0, enter - for default, default: 0)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/az_integrate',
             data, width, azi, cflg, scale, lz], logpath=logpath)


def az_spec_SLC(SLC, SLC_par, spectrum, roff='-', namb='-', pltflg='-', logpath=None):
    """
    | Doppler centroid estimate from SLC images
    | Copyright 2016, Gamma Remote Sensing, v2.9 clw 15-Feb-2016

    Parameters
    ----------
    SLC:
        (input) SAR image data file (fcomplex or scomplex format)
    SLC_par:
        (input) ISP SLC image parameter file
    spectrum:
        (output) Doppler spectrum (text format)
    roff:
        range sample offset to center of estimation window (enter - for default=center_swath)
    namb:
        number of multiples of the PRF to add to the estimated centroid (default=0)
    pltflg:
        azimuth spectrum plotting flag:
            * 0: none (default)
            * 1: output plot in PNG format

    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/az_spec_SLC',
             SLC, SLC_par, spectrum, roff, namb, pltflg], logpath=logpath)


def base_copy(SLC1_par, baseline_1, SLC2_par, baseline_2, time_rev='-', logpath=None):
    """
    | Calculate baseline file for a subsection of a reference SLC
    | Copyright 2003, Gamma Remote Sensing, v1.1 6-Jan-2003 ts/clw/uw

    Parameters
    ----------
    SLC1_par:
        (input) ISP image parameter file of the reference SLC
    baseline-1:
        (input) baseline file derived using the reference SLC geometry
    SLC2_par:
        (input) ISP image parameter file corresponding to the subsecton of the reference SLC
    baseline-2:
        (output) baseline file derived using the geometry and timing of the SLC subsection
    time_rev:
        SLC image normal=1,  time-reversed = -1 (default=1)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/base_copy', SLC1_par,
             baseline_1, SLC2_par, baseline_2, time_rev], logpath=logpath)


def base_est_fft(interf, SLC1_par, OFF_par, baseline, nazfft='-', r_samp='-', az_line='-', logpath=None):
    """
    | Estimate baseline from interferogram fringe spectrum
    | Copyright 2016, Gamma Remote Sensing, v2.1 clw/uw 20-Feb-2016

    Parameters
    ----------
    interf:
        (input) multi-look interferogram with range phase 
    SLC1_par:
        (input) SLC-1 ISP image parameter file
    OFF_par:
        (input) ISP offset/interferogram parameter file
    baseline:
        (output) baseline file
    nazfft:
        size of azimuth FFT (lines read from file, 2\*\*N) (default: 512)
    r_samp:
        range pixel offset to center of the FFT window (default: center)
    az_line:
        line offset from start of the interf. for the  FFT window (default=center)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/base_est_fft', interf,
             SLC1_par, OFF_par, baseline, nazfft, r_samp, az_line], logpath=logpath)


def base_ls(SLC_par, OFF_par, gcp_ph, baseline, ph_flag='-', bc_flag='-', bn_flag='-', bcdot_flag='-', bndot_flag='-', bperp_min='-', SLC2R_par='-', logpath=None):
    """
    | Least squares baseline estimation using terrain heights
    | Copyright 2005, Gamma Remote Sensing, v2.2 5-Sep-2005 clw/uw

    Parameters
    ----------
    SLC_par:
        (input) ISP parameter file of the reference SLC
    OFF_par:
        (input) ISP interferogram/offset parameter file
    gcp_ph:
        (input) ground control point heights + extracted unwrapped phase (text format)
    baseline:
        (input) baseline parameter file
    ph_flag:
        restore range phase ramp (default=0: do not restore  1: restore)
    bc_flag:
        cross-track baseline component estimate (0:orbit derived  1:estimate from data, default=1)
    bn_flag:
        normal baseline component estimate      (0:orbit derived  1:estimate from data, default=1)
    bcdot_flag:
        cross-track baseline rate estimate      (0:orbit derived  1:estimate from data, default=1)
    bndot_flag:
        normal baseline rate estimate           (0:orbit derived  1:estimate from data, default=0)
    bperp_min:
        minimum perpendicular baseline required for L.S estimation (m, default=  10.0)
    SLC2R_par:
        (input) parameter file of resampled SLC, required if SLC-2 frequency differs from SLC-1
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/base_ls', SLC_par, OFF_par, gcp_ph, baseline,
             ph_flag, bc_flag, bn_flag, bcdot_flag, bndot_flag, bperp_min, SLC2R_par], logpath=logpath)


def base_orbit(SLC1_par, SLC2_par, baseline, logpath=None):
    """
    | Estimate baseline from orbit state vectors
    | Copyright 2015, Gamma Remote Sensing, v4.1 clw 18-Apr-2015

    Parameters
    ----------
    SLC1_par:
        (input) SLC-1 ISP image parameter file
    SLC2_par:
        (input) SLC-2 ISP image parameter file
    baseline:
        (output) baseline file (text format, enter - for none)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/base_orbit',
             SLC1_par, SLC2_par, baseline], logpath=logpath)


def base_perp(baseline, SLC1_par, OFF_par, time_rev='-', logpath=None):
    """
    | Calculate baseline components perpendicular and parallel to look vector
    | Copyright 2005, Gamma Remote Sensing, v3.5 10-May-2005 clw/uw

    Parameters
    ----------
    baseline:
        (input) baseline file (text)
    SLC1_par:
        (input) ISP parameter file of SLC-1 (reference SLC)
    OFF_par:
        (input) ISP interferogram/offset parameter file
    time_rev:
        SLC image normal=1 (default), image time-reversed = -1
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/base_perp',
             baseline, SLC1_par, OFF_par, time_rev], logpath=logpath)


def bpf(data_in, data_out, width, fc_x, bw_x, fc_y, bw_y, roff='-', azoff='-', nr='-', naz='-', data_type='-', f_mode='-', beta='-', fir_len='-', logpath=None):
    """
    | Interferometric SAR Processor (ISP): Program /cluster/GAMMA_SOFTWARE-20161207/ISP/bin/bpf.c
    | Copyright 2016, Gamma Remote Sensing, v1.7 clw 3-Mar-2016
    | Bandpass filter for 2-dimensional complex image data (FCOMPLEX or SCOMPLEX format)

    Parameters
    ----------
    data_in:
        (input) input data file (fcomplex, scomplex, float)
    data_out:
        (output) output data file (fcomplex, scomplex, float)
    width:
        number of samples/line
    fc_x:
        normalized x-coord. (across) filter center frequency (range: -0.5 --> 0.5)
    bw_x:
        normalized x-coord. bandwidth (range: 0 --> 1.0)
    fc_y:
        normalized y-coord. (down) filter center frequency (range: -0.5 --> 0.5)
    bw_y:
        normalized y-coord. bandwidth (range: 0 --> 1.0)
    roff:
        offset to starting range to filter   (default: 0)
    azoff:
        offset to starting azimuth to filter (default: 0)
    nr:
        number of range pixels to filter  (default - : width - roff)
    naz:
        number of azimuth lines to filter (default - : nlines - azoff)
    data_type:
        data type (default 0:fcomplex,  1:scomplex,  2:float)
    f_mode:
        fill mode (default 0:force filtered value to 0.0 for input value 0.0, 1:no forcing)
    beta:
        Kaiser window beta parameter (default - :    1.000)
    fir_len:
        finite impulse reponse filter length (default - : 128)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/bpf', data_in, data_out, width, fc_x,
             bw_x, fc_y, bw_y, roff, azoff, nr, naz, data_type, f_mode, beta, fir_len], logpath=logpath)


def bridge(int, flag, unw, bridge, width, xmin='-', xmax='-', ymin='-', ymax='-', logpath=None):
    """
    | Phase unwrap new regions with bridges to regions already unwrapped
    | Copyright 2010, Gamma Remote Sensing, v1.5 clw 4-Nov-2010

    Parameters
    ----------
    int:
        (input) interferogram (fcomplex)
    flag:
        (input) unwrapping flag file
    unw:
        (input/output) unwrapped phase (float) 
    bridge:
        (input) bridge data file (text format)
    width:
        number of samples/row
    xmin:
        starting range pixel offset to unwrap (default = 0)
    xmax:
        last range pixel offset to unwrap (default=width-1)
    ymin:
        starting azimuth row offset to unwrap, relative to start (default = 0)
    ymax:
        last azimuth row offset to unwrap, relative to start (default = nlines-1)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/bridge', int,
             flag, unw, bridge, width, xmin, xmax, ymin, ymax], logpath=logpath)


def cc_wave(interf, pwr1, pwr2, corr, width, bx, by, wflg, xmin, xmax, ymin, ymax, logpath=None):
    """
    | Estimate interferometric coherence
    | Copyright 2015, Gamma Remote Sensing, v5.8 27-Jan-2015 clw/uw

    Parameters
    ----------
    interf:
        (input) normalized complex interferogram
    pwr1:
        (input) intensity image of the first scene (enter - for none)
    pwr2:
        (input) intensity image of the second scene (enter - for none)
    corr:
        (output) estimated degree of coherence filename
    width:
        number of samples/row
    bx:
        coherence window size (columns) (default: 5.0)
    by:
        coherence window size (rows) (default: 5.0)
    wflg:
        magnitude weighting function:
            * 0: constant (default)
            * 1: triangular
            * 2: gaussian
            * 3: none (phase only)

    xmin:
        starting range pixel offset (default = 0)
    xmax:
        last range pixel offset (default = width-1)
    ymin:
        starting azimuth row offset, relative to start (default = 0)
    ymax:
        last azimuth row offset, relative to start (default = nlines-1)
            * NOTE: omitting pwr1 and pwr2 or setting wflg = 3 selects a coherence estimate algorithm that only
              uses the complex interferogram values. In the case of wflg = 3, only the interferogram phase is used.
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/cc_wave', interf, pwr1,
             pwr2, corr, width, bx, by, wflg, xmin, xmax, ymin, ymax], logpath=logpath)


def clear_flag(flag, width, flag_bits, xmin, xmax, ymin, ymax, logpath=None):
    """
    | Clear phase unwrapping flag bits
    | Copyright 2005, Gamma Remote Sensing, v1.6 clw 17-Oct-2005

    Parameters
    ----------
    flag:
        (input)phase unwrapping flag filename 
    width:
        number of samples/row
    flag_bits:
        byte with value of flag(s) to be cleared: 
            Charges = 3	Guides = 4	Low SNR = 8	Visited = 16
            BRANCH PT. = 32	Cuts   = 64	Lawn    = 128
    xmin:
        starting range pixel offset (default = 0)
    xmax:
        last range pixel offset (default = width-1)
    ymin:
        starting azimuth row offset, relative to start (default = 0)
    ymax:
        last azimuth row offset, relative to start (default = nlines-1)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/clear_flag',
             flag, width, flag_bits, xmin, xmax, ymin, ymax], logpath=logpath)


def corr_flag(corr, flag, width, corr_thr, xmin='-', xmax='-', ymin='-', ymax='-', border='-', logpath=None):
    """
    | Low correlation region detection for phase unwrapping
    | Copyright 2005, Gamma Remote Sensing, v2.4 1-Mar-2005 clw/uw

    Parameters
    ----------
    corr:
        (input)interferometric correlation file
    flag:
        (input/output) phase unwrapping flag filename 
    width:
        number of samples/row
    corr_thr:
        corrrelation threshold (0 --> 1.0)
    xmin:
        starting range pixel offset (default = 0)
    xmax:
        last range pixel offset (default = width-1)
    ymin:
        starting azimuth row offset, relative to start (default = 0)
    ymax:
        last azimuth row offset, relative to start (default = nlines-1)
    border:
        effective range of low coherence pixels to set low coherence flag (default=2)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/corr_flag', corr, flag,
             width, corr_thr, xmin, xmax, ymin, ymax, border], logpath=logpath)


def create_offset(SLC1_par, SLC2_par, OFF_par, algorithm='-', rlks='-', azlks='-', iflg='-', logpath=None):
    """
    | Create and update ISP offset and interferogram parameter files
    | Copyright 2015 Gamma Remote Sensing v5.3 clw/uw 10-Nov-2015

    Parameters
    ----------
    SLC1_par:
        (input) SLC-1/MLI-1 ISP image parameter filename (reference)
    SLC2_par:
        (input) SLC-2/MLI-2 ISP image parameter filename
    OFF_par:
        (input/output) ISP offset/interferogram parameter file
    algorithm:
        offset estimation algorithm
            * 1: intensity cross-correlation (default)
            * 2: fringe visibility

    rlks:
        number of interferogram range looks (enter -  for default: 1)
    azlks:
        number of interferogram azimuth looks (enter - for default: 1)
    iflg:
        interactive mode flag (enter -  for default)
            * 0: non-interactive
            * 1: interactive (default)

    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/create_offset', SLC1_par,
             SLC2_par, OFF_par, algorithm, rlks, azlks, iflg], logpath=logpath)


def dcomp_sirc(infile, outfile, samples, loff='-', nlines='-', logpath=None):
    """
    | Extract SIR-C SLC compressed single-pol data
    | Copyright 2009, Gamma Remote Sensing, v1.4 16-Oct-2009 clw

    Parameters
    ----------
    infile:
        (input) SIR-C single-pol SLC compressed data
    outfile:
        (output) complex floating point data
    samples:
        number of polarimetric samples per input line (4 bytes/sample)
    loff:
        offset to starting line (default: 0)
    nlines:
        number of lines to copy(default: entire file, 0 = entire file)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/dcomp_sirc',
             infile, outfile, samples, loff, nlines], logpath=logpath)


def dcomp_sirc_quad(infile, outfile, samples, parameter, loff='-', nlines='-', logpath=None):
    """
    | Extract SIR-C MLC or SLC compressed quad-pol data
    | Copyright 2009, Gamma Remote Sensing, v1.4 16-Oct-2009 uw/clw

    Parameters
    ----------
    infile:
        (input) SIR-C SLC or MLC quad-pol compressed data
    outfile:
        (output) complex floating point data
    samples:
        number of polarimetric samples per input line (10 bytes/sample)
    parameter:
        polarimetric parameter to extract from SLC or MLC product:
            * 0:  SLC total power
            * 1:  SLC-HH
            * 2:  SLC-HV
            * 3:  SLC-VH
            * 4:  SLC-VV
            * 5:  MLC total power
            * 6:  MLC-HVHV\*
            * 7:  MLC-VVVV\*
            * 8:  MLC-HHHH\*
            * 9:  MLC-HHHV\*
            * 10: MLC-HHVV\*
            * 11: MLC-HVVV\*

    loff:
        offset to starting line (default: 0)
    nlines:
        number of lines to copy(default: entire file, 0 = entire file)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/dcomp_sirc_quad',
             infile, outfile, samples, parameter, loff, nlines], logpath=logpath)


def DELFT_vec2(SLC_par, DELFT_dir, nstate='-', interval='-', ODR='-', logpath=None):
    """
    | Extract and interpolate DELFT ERS-1, ERS-2, and ENVISAT state vectors
    | Copyright 2012, Gamma Remote Sensing, v2.6 clw 24-Oct-2012

    Parameters
    ----------
    SLC_par:
        (input) ISP image parameter file
    DELFT_dir:
        directory containing Delft orbit arclist and ODR files for ERS-1, ERS-2 or ENVISAT
            * NOTE: enter . for current directory

    nstate:
        number of state vectors to generate (enter - for default (>= 15)
    interval:
        time interval between state vectors in the ISP image parameter file (s) (default: 10.0)
    ODR:
        ODR file to use (include path) rather than ODR file determined from the Delft orbit arclist
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/DELFT_vec2',
             SLC_par, DELFT_dir, nstate, interval, ODR], logpath=logpath)


def DORIS_vec(SLC_PAR, DOR, nstate='-', logpath=None):
    """
    | Extract ENVISAT DORIS state vectors and write to an ISP image parameter file
    | Copyright 2008, Gamma Remote Sensing, v1.4 11-Jun-2008 clw

    Parameters
    ----------
    SLC_PAR:
        (input/output)ISP SLC/MLI image parameter file
    DOR:
        (input) ASAR DORIS data file (DOR_VOR_AXVF)
    nstate:
        number of state vectors to extract (enter - for default: 11)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/DORIS_vec',
             SLC_PAR, DOR, nstate], logpath=logpath)


def fspf(data_in, data_out, width, dtype='-', r_max='-', spf_type='-', MLI_par='-', logpath=None):
    """
    | ISP Program /cluster/GAMMA_SOFTWARE-20161207/ISP/bin/fspf.c
    | Copyright 2014, Gamma Remote Sensing, v1.2 28-May-2014 of/clw/uw
    | Fast spatial filter for 2D data

    Parameters
    ----------
    data_in:
        (input) input image data
    data_out:
        (output) spatially filtered image data
    width:
        number of samples/row
    dtype:
        data type (enter - for default):
            * 0: FCOMPLEXn                1: SCOMPLEX
            * 2: FLOAT (default)

    r_max:
        maximum filter radius (range samples) (enter - for default: 64)
    spf_type:
        spatial filter type (enter - for default):
            * 0: uniform average (default for fcomplex and scomplex)
            * 1: triangular weighted average: 1 - (r/r_max)
            * 2: quadratic weighted average: 1 - (r/r_max)\*\*2
            * 3: Gaussian weighted average: exp(-2.\*(r\*\*2/r_max\*\*2))
            * 4: linear least-squares (default for float data)

    MLI_par:
        MLI or SLC parameter file with the same number of looks as the input image, required for GPRI data
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/fspf', data_in,
             data_out, width, dtype, r_max, spf_type, MLI_par], logpath=logpath)


def gcp_phase(unw, OFF_par, gcp, gcp_ph, win_sz='-', logpath=None):
    """
    | Extract unwrapped phase at GCP locations
    | Copyright 2006, Gamma Remote Sensing, v1.5 8-Mar-2006 clw

    Parameters
    ----------
    unw:
        (input) unwrapped interferometric phase
    OFF_par:
        (input) ISP interferogram/offset parameter file
    gcp:
        (input) ground control point data (text format)
    gcp_ph:
        (output) ground control point data + extracted unwrapped phase (text)
    win_sz:
        window size for averaging phase for each gcp, must be odd (default: 1)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/gcp_phase',
             unw, OFF_par, gcp, gcp_ph, win_sz], logpath=logpath)


def grasses(int, flag, unw, width, xmin='-', xmax='-', ymin='-', ymax='-', xinit='-', yinit='-', init_ph='-', logpath=None):
    """
    | Phase unwrapping by region growing
    | Copyright 2005, Gamma Remote Sensing, v4.2 1-Mar-2005 clw/uw

    Parameters
    ----------
    int:
        (input) interferogram filename
    flag:
        (input) unwrapping flag filename
    unw:
        (output) unwrapped phase filename
    width:
        number of samples/row
    xmin:
        starting range pixel offset (default = 0)
    xmax:
        last range pixel offset (default=width-1)
    ymin:
        starting azimuth row offset, relative to start (default = 0)
    ymax:
        last azimuth row offset, relative to start (default = nlines-1)
    xinit:
        starting range pixel for unwrapping (default = width/2)
    yinit:
        starting row to unwrap (default = height/2)
    init_ph:
        flag to set phase at starting point to 0.0 (default 0: not set to 0.0, 1: set to 0.0)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/grasses', int, flag, unw,
             width, xmin, xmax, ymin, ymax, xinit, yinit, init_ph], logpath=logpath)


def hgt_map(unw, SLC_par, OFF_par, baseline, hgt, gr, ph_flag='-', loff='-', nlines='-', SLC2R_par='-', logpath=None):
    """
    | Interferometric height/ground range estimation vs. slant range
    | Copyright 2005, Gamma Remote Sensing, v5.1 clw/uw 9-Sep-2005

    Parameters
    ----------
    unw:
        (input) unwrapped interferometric phase
    SLC_par:
        (input) ISP parameter file for the reference SLC
    OFF_par:
        (input) ISP offset/interferogram processing parameters
    baseline:
        (input) baseline parameter file
    hgt:
        (output) height file (in slant range geometry) relative to the WGS-84 ellipsoid
    gr:
        (output) cross-track ground ranges on the WGS-84 ellipsoid (in slant range geometry)
    ph_flag:
        restore phase slope flag (0:no phase change  default=1:add back phase ramp)
    loff:
        offset to starting line (default = 0)
    nlines:
        number of lines to calculate (enter - for default: to end of file)
    SLC2R_par:
        (input) parameter file of resampled SLC, required if SLC-2 frequency differs from SLC-1
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/hgt_map', unw, SLC_par,
             OFF_par, baseline, hgt, gr, ph_flag, loff, nlines, SLC2R_par], logpath=logpath)


def image_stat(image, width, roff, loff, nr, nl, report, logpath=None):
    """
    | Calculate mean, standard deviation and number of non-zero values for a rectangular image region (float format)
    | Copyright 2016, Gamma Remote Sensing, v1.3 3-Nov-2016

    Parameters
    ----------
    image:
        (input) image data file (float)
    width:
        image line width (samples/line)
    roff:
        sample offset to region start (enter - for default: 0)
    loff:
        line offset to region start (enter - for default: 0)
    nr:
        region width (samples, enter - for default: width - roff)
    nl:
        number of lines in the region (enter - for default: image_lines - loff)
    report:
        output text file (keyword:value format)
            keywords: file, mean, stdev, total_samples, non_zero_samples, fraction_valid)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/image_stat',
             image, width, roff, loff, nr, nl, report], logpath=logpath)


def init_offset(SLC_1, SLC_2, SLC1_par, SLC2_par, OFF_par, rlks='-', azlks='-', rpos='-', azpos='-', offr='-', offaz='-', thres='-', rwin='-', azwin='-', cflag='-', logpath=None):
    """
    | Determine initial offset between SLC images using correlation of image intensity
    | Copyright 2016, Gamma Remote Sensing, v3.1 clw 12-Apr-2016

    Parameters
    ----------
    SLC-1:
        (input) single-look complex image 1 (reference)
    SLC-2:
        (input) single-look complex image 2
    SLC1_par:
        (input) SLC-1 ISP image parameter file
    SLC2_par:
        (input) SLC-2 ISP image parameter file
    OFF_par:
        (input) ISP offset/interferogram parameter file
    rlks:
        number of range looks (default: 1)
    azlks:
        number of azimuth looks (default: 1)
    rpos:
        center of patch in range (samples) (enter - for default: image center)
    azpos:
        center of patch in azimuth (lines) (enter - for default: image center)
    offr:
        initial range offset (samples) (enter - for default: 0)
    offaz:
        initial azimuth offset (lines) (enter - for default: 0)
    thres:
        cross-correlation threshold (enter - for default: 0.150)
    rwin:
        range window size (default: 512)
    azwin:
        azimuth window size (default: 512)
    cflag:
        copy offsets to the range and azimuth offset polynomials in the OFF_par
            * 0: do not copy
            * 1: copy constant range and azimuth offset (default)

    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/init_offset', SLC_1, SLC_2, SLC1_par, SLC2_par,
             OFF_par, rlks, azlks, rpos, azpos, offr, offaz, thres, rwin, azwin, cflag], logpath=logpath)


def init_offset_orbit(SLC1_par, SLC2_par, OFF_par, rpos='-', azpos='-', cflag='-', logpath=None):
    """
    | Initial SLC image offset estimation from orbit state-vectors and image parameters
    | Copyright 2016, Gamma Remote Sensing, v1.7 21-Apr-2016 clw/uw

    Parameters
    ----------
    SLC1_par:
        (input) SLC-1 parameter file
    SLC2_par:
        (input) SLC-2 parameter file
    OFF_par:
        (input/output) ISP/offset parameter file
    rpos:
        range position for offset estimation (enter - for default: center of SLC-1)
    azpos:
        azimuth position for offset estimation (enter - for default: center of SLC-1)
    cflag:
        copy offsets to the range and azimuth offset polynomials in the OFF_par
            * 0: do not copy
            * 1: copy constant range and azimuth offset (default)

    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/init_offset_orbit',
             SLC1_par, SLC2_par, OFF_par, rpos, azpos, cflag], logpath=logpath)


def interp_ad(data_in, data_out, width, r_max='-', np_min='-', np_max='-', w_mode='-', type='-', cp_data='-', logpath=None):
    """
    | Weighted interpolation of gaps in 2D data using an adaptive smoothing window
    | Copyright 2016, Gamma Remote Sensing, v2.1 23-Nov-2016 clw/uw

    Parameters
    ----------
    data_in:
        (input) data with gaps
    data_out:
        (output) data with gaps filled by interpolation
    width:
        number of samples/row
    r_max:
        maximum interpolation window radius (default(-): 16)
    np_min:
        minimum number of points used for the interpolation (default(-): 16)
    np_max:
        maximum number of points used for the interpolation (default(-): 16)
    w_mode:
        data weighting mode (enter - for default):
            * 0: constant
            * 1: 1 - (r/r_max)
            * 2: 1 - (r/r_max)\*\*2  (default)
            * 3: exp(-2.\*(r\*\*2/r_max\*\*2))

    type:
        input and output data type:
            * 0: FCOMPLEX
            * 1: SCOMPLEX
            * 2: FLOAT (default)
            * 3: INT
            * 4: SHORT

    cp_data:
        copy data flag:
            * 0: do not copy input data values to output
            * 1: copy input data values to output (default)

    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/interp_ad', data_in, data_out,
             width, r_max, np_min, np_max, w_mode, type, cp_data], logpath=logpath)


def mask_data(data_in, width, data_out, mask, format_flag='-', logpath=None):
    """
    | Mask float or fcomplex data using an 8-bit SUN/BMP/TIFF format raster image
    | Copyright 2015, Gamma Remote Sensing, v1.3 3-Dec-2015 clw

    Parameters
    ----------
    data_in:
        (input) data file (FLOAT or FCOMPLEX format)
    width:
        width of input data file
    data_out:
        (output) data file, same data format as input
    mask:
        (input) mask file, SUN/BMP/TIFF raster format, 8-bits/pixel
            output data values are set to 0.0 at all locations where the mask is black (0,0,0) or dn = 0
            * NOTE: mask file must have the same width as the input data file

    format_flag:
        data format:
            * 0: FLOAT (default)
            * 1: FCOMPLEX

    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/mask_data',
             data_in, width, data_out, mask, format_flag], logpath=logpath)


def mcf(interf, wgt, mask, unw, width, tri_mode='-', roff='-', loff='-', nr='-', nlines='-', npat_r='-', npat_az='-', ovrlap='-', r_init='-', az_init='-', init_flag='-', logpath=None):
    """
    | Phase unwrapping using Minimum Cost Flow (MCF) and triangulation
    | Copyright 2016, Gamma Remote Sensing, v2.2 clw/uw 30-Nov-2016

    Parameters
    ----------
    interf:
        (input) interferogram (\*.int,\*.flt)(fcomplex)
    wgt:
        (input) weight factors (0 -> 1.0) file (float)(enter - for uniform weight)
    mask:
        (input) validity mask (SUN/BMP/TIFF raster format, value 0 -> pixel not used)(enter - if no mask)
    unw:
        (output) unwrapped phase image (\*.unw)(float)
    width:
        number of samples/row
    tri_mode:
        triangulation mode
            * 0: filled triangular mesh (default)
            * 1: Delaunay triangulation

    roff:
        offset to starting range of section to unwrap (default: 0)
    loff:
        offset to starting line of section to unwrap (default: 0)
    nr:
        number of range samples of section to unwrap (default(-): width - roff)
    nlines:
        number of lines of section to unwrap (default(-): total number of lines - loff)
    npat_r:
        number of patches in range
    npat_az:
        number of patches in azimuth
    ovrlap:
        overlap between patches in pixels (overlap >= 7, default(-): 512)
    r_init:
        phase reference point range offset (default(-): roff)
    az_init:
        phase reference point azimuth offset (default(-): loff)
    init_flag:
        flag to set phase at reference point
            * 0: use initial point phase value (default)
            * 1: set phase to 0.0 at initial point

    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/mcf', interf, wgt, mask, unw, width, tri_mode,
             roff, loff, nr, nlines, npat_r, npat_az, ovrlap, r_init, az_init, init_flag], logpath=logpath)


def MLI_cat(MLI_1, MLI_2, MLI1_par, MLI2_par, MLI_3, MLI3_par, logpath=None):
    """
    | Concatenate two MLI images using bicubic spline interpolation
    | Copyright 2015, Gamma Remote Sensing, v1.0 23-Jul-2015 awi

    Parameters
    ----------
    MLI-1:
        (input) MLI-1 image (single-look)
    MLI-2:
        (input) MLI-2 image to be appended to MLI-1
    MLI1_par:
        (input) MLI-1 ISP image parameter file
    MLI2_par:
        (input) MLI-2 ISP image parameter file
    MLI-3:
        (output) concatenated MLI image
    MLI3_par:
        (output) ISP image parameter file for concatenated image
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/MLI_cat', MLI_1,
             MLI_2, MLI1_par, MLI2_par, MLI_3, MLI3_par], logpath=logpath)


def MLI_copy(MLI_in, MLI_in_par, MLI_out, MLI_out_par, roff='-', nr='-', loff='-', nl='-', logpath=None):
    """
    | Copy MLI data file with options for segment extraction
    | Copyright 2013, Gamma Remote Sensing, v4.4 10-Jan-2013 uw/clw

    Parameters
    ----------
    MLI_in:
        (input) multi-look intensity image (float format)
    MLI_in_par:
        (input) ISP image parameter file for input MLI
    MLI_out:
        (output) selected MLI section (float format)
    MLI_out_par:
        (output) ISP image parameter file for output MLI
    roff:
        offset to starting range sample (enter - for default: 0)
    nr:
        number of range samples (enter - for default: to end of line
    loff:
        offset to starting line (enter - for default: 0)
    nl:
        number of lines to copy (enter - for default: to end of file)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/MLI_copy', MLI_in,
             MLI_in_par, MLI_out, MLI_out_par, roff, nr, loff, nl], logpath=logpath)


def mosaic_WB(data_tab, dtype, data_out, data_par_out, sc_flg='-', logpath=None):
    """
    | ISP: Program /cluster/GAMMA_SOFTWARE-20161207/ISP/bin/mosaic_WB.c
    | Copyright 2011, Gamma Remote Sensing, v1.2 6-Apr-2011 clw
    | Mosaic Wide-Beam ScanSAR data processed by the MSP

    Parameters
    ----------
    data_tab:
        (input) 2 column list of data  and ISP image parameter files for the beams in the mosaic (text)
    dtype:
        (input) input data type:
            * 0: FLOAT
            * 1: FCOMPLEX

    data_out:
        (output) output image mosaic
    data_par_out:
        (output) ISP image parameter file for output image mosaic
    sc_flg:
        intensity scaling flag:
            * 0: do not scale different beam data values
            * 1: use overlap regions to scale beam intensities (default)

    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/mosaic_WB',
             data_tab, dtype, data_out, data_par_out, sc_flg], logpath=logpath)


def multi_cpx(data_in, OFF_par_in, data_out, OFF_par_out, rlks='-', azlks='-', loff='-', nlines='-', roff='-', nsamp='-', logpath=None):
    """
    | Calculate multi-look averaged or interpolated 2D image (fcomplex data)
    | Copyright 2013, Gamma Remote Sensing, v2.5 28-Mar-2013 clw/uw

    Parameters
    ----------
    data_in:
        (input) input fcomplex image file
    OFF_par_in:
        (input) offset parameter file for input image
    data_out:
        (output) output multi-look or interpolated fcomplex data file
    OFF_par_out:
        (input/output) offset parameter file for output, if already exists, then used as input
    rlks:
        number of range looks, values < -1, interpreted as an image oversampling factor (default: 1)
    azlks:
        number azimuth looks,  values < -1, interpreted as an image oversampling factor (default: 1)
    loff:
        line offset to starting line (default: 0)
    nlines:
        number of lines (default: 0, to end of file)
    roff:
        offset to starting range sample (default:0)
    nsamp:
        number of range samples to extract (default: 0, to end of line)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/multi_cpx', data_in, OFF_par_in,
             data_out, OFF_par_out, rlks, azlks, loff, nlines, roff, nsamp], logpath=logpath)


def multi_look(SLC, SLC_par, MLI, MLI_par, rlks, azlks, loff='-', nlines='-', scale='-', exp='-', logpath=None):
    """
    | Calculate a multi-look intensity (MLI) image from an SLC image
    | Copyright 2016, Gamma Remote Sensing, v4.1 18-Nov-2016 clw/uw

    Parameters
    ----------
    SLC:
        (input) single-look complex image
    SLC_par:
        (input) SLC ISP image parameter file
    MLI:
        (output) multi-look intensity image
    MLI_par:
        (output) MLI ISP image parameter file
    rlks:
        number of range looks
    azlks:
        number of azimuth looks
    loff:
        offset to starting line (default: 0)
    nlines:
        number of SLC lines to process (enter - for default: entire file)
    scale:
        scale factor for output MLI (default: 1.0)
    exp:
        exponent for the output MLI (default: 1.0)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/multi_look', SLC, SLC_par,
             MLI, MLI_par, rlks, azlks, loff, nlines, scale, exp], logpath=logpath)


def multi_real(data_in, OFF_par_in, data_out, OFF_par_out, rlks='-', azlks='-', loff='-', nlines='-', roff='-', nsamp='-', logpath=None):
    """
    | Calculate multi-look averaged or interpolated 2D image (float data)
    | Copyright 2012, Gamma Remote Sensing, v2.5 16-Jul-2013 clw/uw

    Parameters
    ----------
    data_in:
        (input) input float image file
    OFF_par_in:
        (input) interferogram/offset parameter file for input image
    data_out:
        (output) output multi-look or interpolated float data file
    OFF_par_out:
        (input/output) interferogram/offset parameter file for output, if already existing, used as input
    rlks:
        number of range looks, values < -1, interpreted as an image oversampling factor (default: 1)
    azlks:
        number azimuth looks,  values < -1, interpreted as an image oversampling factor (default: 1)
    loff:
        line offset to starting line (default:0)
    nlines:
        number of lines (default:0, to end of file)
    roff:
        offset to starting range sample (default:0)
    nsamp:
        number of range samples to extract (default:0, to end of line)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/multi_real', data_in, OFF_par_in,
             data_out, OFF_par_out, rlks, azlks, loff, nlines, roff, nsamp], logpath=logpath)


def multi_S1_TOPS(SLC_tab, MLI, MLI_par, rlks, azlks, wflg='-', SLCR_tab='-', logpath=None):
    """
    | Calculate MLI mosaic from Sentinel-1 TOPS SLC burst data (FCOMPLEX and SCOMPLEX)
    | Copyright 2016, Gamma Remote Sensing v3.3 23-Aug-2016 awi/clw/uw

    Parameters
    ----------
    SLC_tab:
        (input) 3 column list of SLC, SLC_par, Sentinel-1 TOPS_par, rows sorted in the order IW1, IW2, IW3
    MLI:
        (output) multi-look intensity image
    MLI_par:
        (output) MLI image parameter file
    rlks:
        number of range looks
    azlks:
        number of azimuth looks
    wflg:
        burst window calculation flag:
            * 0: use existing burst window parameters if they exist, otherwise calculate burst window parameters (default)
            * 1: calculate burst window parameters from burst parameters and the number of range and azimuth looks

    SLCR_tab:
        (input) SLC_tab of the reference scene, 3 column list of SLC, SLC_par, TOPS_par sorted in the order IW1, IW2, IW3
            * NOTE: When generating an MLI mosaic of a resampled SLC, the SLC_tab of the reference scene is required

    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/multi_S1_TOPS',
             SLC_tab, MLI, MLI_par, rlks, azlks, wflg, SLCR_tab], logpath=logpath)


def multi_SLC_WSS(SLC, SLC_par, MLI, MLI_par, logpath=None):
    """
    | Calculate multi-look intensity image (MLI) from a ASAR Wide-Swath SLC
    | Copyright 2008, Gamma Remote Sensing v1.2 08-Jan-2008 clw/awi

    Parameters
    ----------
    SLC:
        (input) ASAR Wide-Swath SLC image
    SLC_par:
        (input) ASAR Wide-Swath SLC image parameter file
    MLI:
        (output) multi-look intensity image
    MLI_par:
        (output) MLI image parameter file
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/multi_SLC_WSS',
             SLC, SLC_par, MLI, MLI_par], logpath=logpath)


def neutron(intensity, flag, width, n_thres, ymin='-', ymax='-', logpath=None):
    """
    | Generate phase unwrapping neutrons using image intensity
    | Copyright 2014, Gamma Remote Sensing, v2.3 20-Jan-2014 clw/uw

    Parameters
    ----------
    intensity:
        (input) image intensity 
    flag:
        (input) phase unwrapping flag file
    width:
        number of samples/row
    n_thres:
        neutron threshold, multiples of the average intensity (default=6.0)
    ymin:
        offset to starting azimuth row (default = 0)
    ymax:
        offset to last azimuth row (default = nlines-1)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/neutron',
             intensity, flag, width, n_thres, ymin, ymax], logpath=logpath)


def offset_add(OFF_par1, OFF_par2, OFF_par3, logpath=None):
    """
    | Add range and azimuth offset polynomial coefficients
    | Copyright 2008, Gamma Remote Sensing, v1.1 12-Feb-2008 clw

    Parameters
    ----------
    OFF_par1:
        (input) ISP offset/interferogram parameter file
    OFF_par2:
        (input) ISP offset/interferogram parameter file
    OFF_par3:
        (output) ISP offset/interferogram parameter file with sums of the
            range and azimuth offset polynomials in OFF_par1 and OFF_par2
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/offset_add',
             OFF_par1, OFF_par2, OFF_par3], logpath=logpath)


def offset_pwr(SLC1, SLC2, SLC1_par, SLC2_par, OFF_par, offs, ccp, rwin='-', azwin='-', offsets='-', n_ovr='-', nr='-', naz='-', thres='-', c_ovr='-', pflag='-', pltflg='-', ccs='-', logpath=None):
    """
    | Offset tracking between SLC images using intensity cross-correlation
    | Copyright 2016, Gamma Remote Sensing, v5.1 clw 22-Oct-2016

    Parameters
    ----------
    SLC1:
        (input) single-look complex image 1 (reference)
    SLC2:
        (input) single-look complex image 2
    SLC1_par:
        (input) SLC-1 ISP image parameter file
    SLC2_par:
        (input) SLC-2 ISP image parameter file
    OFF_par:
        (input) ISP offset/interferogram parameter file
    offs:
        (output) offset estimates in range and azimuth (fcomplex)
    ccp:
        (output) cross-correlation of each patch (0.0->1.0) (float)
    rwin:
        range patch size (range pixels, (enter - for default from offset parameter file)
    azwin:
        azimuth patch size (azimuth lines, (enter - for default from offset parameter file)
    offsets:
        (output) range and azimuth offsets and cross-correlation data in text format, enter - for no output
    n_ovr:
        SLC oversampling factor (integer 2\*\*N (1,2,4,8), enter - for default: 2)
    nr:
        number of offset estimates in range direction (enter - for default from offset parameter file)
    naz:
        number of offset estimates in azimuth direction (enter - for default from offset parameter file)
    thres:
        cross-correlation threshold (enter - for default from offset parameter file)
    c_ovr:
        correlation function oversampling factor (integer 2\*\*N (1,2,4,8,16) default: 4)
    pflag:
        print flag (enter - for default)
            * 0: print offset summary
            * 1: print all offset data

    pltflg:
        plotting flag (enter - for default)
            * 0: none (default)
            * 1: screen output
            * 2: screen output and PNG format plots
            * 3: output plots in PDF format

    ccs:
        (output) cross-correlation standard deviation of each patch (float)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/offset_pwr', SLC1, SLC2, SLC1_par, SLC2_par, OFF_par,
             offs, ccp, rwin, azwin, offsets, n_ovr, nr, naz, thres, c_ovr, pflag, pltflg, ccs], logpath=logpath)


def offset_pwr_tracking(SLC1, SLC2, SLC1_par, SLC2_par, OFF_par, offs, ccp, rwin='-', azwin='-', offsets='-', n_ovr='-', thres='-', rstep='-', azstep='-', rstart='-', rstop='-', azstart='-', azstop='-', c_ovr='-', pflag='-', pltflg='-', ccs='-', logpath=None):
    """
    | Offset tracking between SLC images using intensity cross-correlation
    | Copyright 2016, Gamma Remote Sensing, v5.1 clw 22-Oct-2016

    Parameters
    ----------
    SLC1:
        (input) single-look complex image 1 (reference)
    SLC2:
        (input) single-look complex image 2
    SLC1_par:
        (input) SLC-1 ISP image parameter file
    SLC2_par:
        (input) SLC-2 ISP image parameter file
    OFF_par:
        (input) ISP offset/interferogram parameter file
    offs:
        (output) offset estimates in range and azimuth (fcomplex)
    ccp:
        (output) cross-correlation of each patch (0.0->1.0) (float)
    rwin:
        range patch size (range pixels, (enter - for default from offset parameter file)
    azwin:
        azimuth patch size (azimuth lines, (enter - for default from offset parameter file)
    offsets:
        (output) range and azimuth offsets and cross-correlation data in text format, enter - for no output
    n_ovr:
        SLC oversampling factor (integer 2\*\*N (1,2,4,8), enter - for default: 2)
    thres:
        cross-correlation threshold (0--> 1.)(enter - for default from offset parameter file)
    rstep:
        step in range pixels (enter - for default: rwin/2)
    azstep:
        step in azimuth pixels (enter - for default: azwin/2)
    rstart:
        offset to starting range pixel (enter - for default: 0)
    rstop:
        offset to ending range pixel (enter - for default: nr-1)
    azstart:
        offset to starting azimuth line (enter - for default: 0)
    azstop:
        offset to ending azimuth line  (enter - for default: nlines-1)
    c_ovr:
        correlation function oversampling factor (integer 2\*\*N (1,2,4,8,16) default: 4)
    pflag:
        print flag (enter - for default)
            * 0: print offset summary
            * 1: print all offset data)

    pltflg:
        plotting flag (enter - for default)
            * 0: none (default)
            * 1: screen output
            * 2: screen output and PNG format plots
            * 3: output plots in PDF format

    ccs:
        (output) cross-correlation standard deviation of each patch (float)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/offset_pwr_tracking', SLC1, SLC2, SLC1_par, SLC2_par, OFF_par, offs, ccp,
             rwin, azwin, offsets, n_ovr, thres, rstep, azstep, rstart, rstop, azstart, azstop, c_ovr, pflag, pltflg, ccs], logpath=logpath)


def offset_pwr_tracking2(SLC1, SLC2, SLC1_par, SLC2_par, OFF_par, offs, ccp, OFF_par2='-', offs2='-', rwin='-', azwin='-', offsets='-', n_ovr='-', thres='-', rstep='-', azstep='-', rstart='-', rstop='-', azstart='-', azstop='-', c_ovr='-', pflag='-', pltflg='-', ccs='-', logpath=None):
    """
    | Intensity cross-correlation offset tracking with the initial offset for each patch determined from input offset map
    | Copyright 2016, Gamma Remote Sensing, v1.3 clw 22-Oct-2016

    Parameters
    ----------
    SLC1:
        (input) single-look complex image 1 (reference)
    SLC2:
        (input) single-look complex image 2
    SLC1_par:
        (input) SLC-1 ISP image parameter file
    SLC2_par:
        (input) SLC-2 ISP image parameter file
    OFF_par:
        (input) ISP offset/interferogram parameter file
    offs:
        (output) offset estimates in range and azimuth (fcomplex)
    ccp:
        (output) cross-correlation of each patch (0.0->1.0) (float)
    OFF_par2:
        (input) ISP offset/interferogram parameter file of the offset map to determine initial offsets (enter - for none)
    offs2:
        (input) input range and azimuth offset map to determine initial offsets (enter - for none)
    rwin:
        range patch size (range pixels, (enter - for default from offset parameter file)
    azwin:
        azimuth patch size (azimuth lines, (enter - for default from offset parameter file)
    offsets:
        (output) range and azimuth offsets and cross-correlation data in text format, enter - for no output
    n_ovr:
        SLC oversampling factor (integer 2\*\*N (1,2,4,8), enter - for default: 2)
    thres:
        cross-correlation threshold (0--> 1.)(enter - for default from offset parameter file)
    rstep:
        step in range pixels (enter - for default: rwin/2)
    azstep:
        step in azimuth pixels (enter - for default: azwin/2)
    rstart:
        offset to starting range pixel (enter - for default: 0)
    rstop:
        offset to ending range pixel (enter - for default: nr-1)
    azstart:
        offset to starting azimuth line (enter - for default: 0)
    azstop:
        offset to ending azimuth line  (enter - for default: nlines-1)
    c_ovr:
        correlation function oversampling factor (integer 2\*\*N (1,2,4,8,16) default: 4)
    pflag:
        print flag (enter - for default)
            * 0: print offset summary
            * 1: print all offset data

    pltflg:
        plotting flag (enter - for default)
            * 0: none (default)
            * 1: screen output
            * 2: screen output and PNG format plots
            * 3: output plots in PDF format

    ccs:
        (output) cross-correlation standard deviation of each patch (float)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/offset_pwr_tracking2', SLC1, SLC2, SLC1_par, SLC2_par, OFF_par, offs, ccp, OFF_par2,
             offs2, rwin, azwin, offsets, n_ovr, thres, rstep, azstep, rstart, rstop, azstart, azstop, c_ovr, pflag, pltflg, ccs], logpath=logpath)


def offset_SLC(SLC_1, SLC_2, SLC1_par, SLC2_par, OFF_par, offs, snr, rwin='-', azwin='-', offsets='-', n_ovr='-', nr='-', naz='-', thres='-', ISZ='-', pflag='-', logpath=None):
    """
    | Offsets between SLC images using fringe visibility
    | Copyright 2016, Gamma Remote Sensing, v2.9 clw 4-Mar-2016

    Parameters
    ----------
    SLC-1:
        (input) single-look complex image 1 (reference)
    SLC-2:
        (input) single-look complex image 2
    SLC1_par:
        (input) SLC-1 ISP image parameter file
    SLC2_par:
        (input) SLC-2 ISP image parameter file
    OFF_par:
        (input) ISP offset/interferogram parameter file
    offs:
        (output) offset estimates (fcomplex)
    snr:
        (output) offset estimation snr (float)
    rwin:
        search window size (range pixels, (enter - for default from offset parameter file))
    azwin:
        search window size (azimuth lines, (enter - for default from offset parameter file))
    offsets:
        (output) range and azimuth offsets and snr data in text format, enter - for no output
    n_ovr:
        SLC oversampling factor (integer 2\*\*N (1,2,4) default = 2)
    nr:
        number of offset estimates in range direction (enter - for default from offset parameter file)
    naz:
        number of offset estimates in azimuth direction (enter - for default from offset parameter file)
    thres:
        offset estimation quality threshold (enter - for default from offset parameter file)
    ISZ:
        search chip interferogram size (in non-oversampled pixels, default=16)
    pflag:
        print flag (0:print offset summary  default=1:print all offset data)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/offset_SLC', SLC_1, SLC_2, SLC1_par, SLC2_par,
             OFF_par, offs, snr, rwin, azwin, offsets, n_ovr, nr, naz, thres, ISZ, pflag], logpath=logpath)


def offset_SLC_tracking(SLC_1, SLC_2, SLC1_par, SLC2_par, OFF_par, offs, snr, rsw='-', azsw='-', offsets='-', n_ovr='-', thres='-', rstep='-', azstep='-', rstart='-', rstop='-', azstart='-', azstop='-', ISZ='-', pflag='-', logpath=None):
    """
    | Offset tracking between SLC images using fringe visibility
    | Copyright 2016, Gamma Remote Sensing, v3.6 clw 4-Mar-2016

    Parameters
    ----------
    SLC-1:
        (input) single-look complex image 1 (reference)
    SLC-2:
        (input) single-look complex image 2
    SLC1_par:
        (input) SLC-1 ISP image parameter file
    SLC2_par:
        (input) SLC-2 ISP image parameter file
    OFF_par:
        (input) ISP offset/interferogram parameter file
    offs:
        (output) offset estimates (fcomplex)
    snr:
        (output) offset estimation snr (float)
    rsw:
        range search window size (range pixels) (enter - for default from offset parameter file)
    azsw:
        azimuth search window size (azimuth lines) (enter - for default from offset parameter file)
    offsets:
        (output) range and azimuth offsets and snr data in text format, enter - for no output
    n_ovr:
        SLC over-sampling factor (integer 2\*\*N (1,2,4) default: 2)
    thres:
        offset estimation quality threshold (enter - for default from offset parameter file)
    rstep:
        step in range pixels (enter - for default: rsw/2)
    azstep:
        step in azimuth pixels (enter - for default: azsw/2)
    rstart:
        starting range pixel (enter - for default: rsw/2)
    rstop:
        ending range pixel (enter - for default: nr - rsw/2)
    azstart:
        starting azimuth line (enter - for default: azsw/2)
    azstop:
        ending azimuth line  (enter - for default: nlines - azsw/2)
    ISZ:
        search chip interferogram size (in non-oversampled pixels, default: 16)
    pflag:
        print flag:
            * 0: print offset summary  (default)
            * 1: print all offset data

    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/offset_SLC_tracking', SLC_1, SLC_2, SLC1_par, SLC2_par, OFF_par, offs,
             snr, rsw, azsw, offsets, n_ovr, thres, rstep, azstep, rstart, rstop, azstart, azstop, ISZ, pflag], logpath=logpath)


def offset_tracking(offs, ccp, SLC_par, OFF_par, disp_map, disp_val='-', mode='-', thres='-', poly_flag='-', logpath=None):
    """
    | Conversion of range and azimuth offsets files to displacement map
    | Copyright 2015, Gamma Remote Sensing, v1.8 28-Nov-2015 ts/clw/uw

    Parameters
    ----------
    offs:
        (input) range and azimuth offset estimates (fcomplex)
    ccp:
        (input) cross-correlation of the offset estimates (float)
    SLC_par:
        (input) SLC parameter file of reference SLC
    OFF_par:
        (input) offset parameter file used in the offset tracking
    disp_map:
        (output) range and azimuth displacement estimates (fcomplex)
    disp_val:
        (output) range and azimuth displacement estimates and SNR values (enter - for none) (text)
    mode:
        flag indicating displacement mode:
            * 0: displacement in range and azimuth pixels
            * 1: displacement in meters in slant range and azimuth directions
            * 2: displacement in meters in ground range and azimuth directions (default)

    thres:
        SNR threshold to accept offset value (default from OFF_par)
    poly_flag:
        flag indicating if trend calculated using offset polynomials from OFF_par is subtracted:
            * 0: do not subtract polynomial trend from offset data
            * 1: subtract polynomial trend from offset data (default)

    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/offset_tracking', offs, ccp,
             SLC_par, OFF_par, disp_map, disp_val, mode, thres, poly_flag], logpath=logpath)


def ORB_prop_SLC(SLC_par, nstate='-', interval='-', extra='-', mode='-', logpath=None):
    """
    | Calculate state vectors using orbit propagation and interpolation
    | Copyright 2008, Gamma Remote Sensing, v1.8 11-Jun-2008 clw/awi

    Parameters
    ----------
    SLC_par:
        (input) ISP image parameter file with at least 1 state vector
    nstate:
        number of state vectors to calculate (enter - for default: nstate from image duration + extra)
    interval:
        time interval between state vectors (enter - for default: state vector time interval in SLC_par)
    extra:
        extra time for state vectors at start and end of image (sec.) (enter - for default: 30.0)
    mode:
        orbit propagation mode:
            * 0: polynomial interpolation (default, if 3 or more state vectors available)
            * 1: integration of the equations of motion (default, if less than 3 state vectors available)
            * 2: interpolate between state vectors, minimum of 3 state vectors;
              interpolation of the equations of motion outside of the time span of the existing state vectors
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/ORB_prop_SLC',
             SLC_par, nstate, interval, extra, mode], logpath=logpath)


def ORRM_vec(SLC_par, ORRM, nstate='-', logpath=None):
    """
    | Calculate state vectors extraction from ORRM file
    | Copyright 2008, Gamma Remote Sensing, v1.4 15-Nov-2004 clw

    Parameters
    ----------
    SLC_par:
        (input/output) ISP SLC/MLI image parameter file
    ORRM:
        (input) ORRM state vector file
    nstate:
        number of state vectors (default=5, maximum=64)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/ORRM_vec',
             SLC_par, ORRM, nstate], logpath=logpath)


def par_ACS_ERS(CEOS_SAR_leader, SLC_par, logpath=None):
    """
    | ISP parameter file generation for ERS SLC data from the ACS processor
    | Copyright 2005, Gamma Remote Sensing, v1.3 17-Oct-2005 clw/uw

    Parameters
    ----------
    CEOS_SAR_leader:
        (input) ERS CEOS SAR leader file
    SLC_par:
        (output) ISP SLC parameter file (example <orbit>.slc.par)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_ACS_ERS',
             CEOS_SAR_leader, SLC_par], logpath=logpath)


def par_ASAR(ASAR_file, output_name, K_dB='-', logpath=None):
    """
    | Extract SLC/MLI image parameters and images from ENVISAT ASAR SLC, WSS, APP, and PRI products
    | Copyright 2014, Gamma Remote Sensing, v2.7 20-Aug-2014 clw/uw/awi

    Parameters
    ----------
    ASAR_file:
        (input)ASAR data file including header and image as provided by ESA
    output_name:
        (output)common part of output file names (e.g. orbit number)
    K_dB:
        Calibration factor in dB (nominal value for all ASAR modes = 55.0)
            * NOTE: Use - for the calibration coefficient provided in the header of the ASAR_file
            * NOTE: In the case that a calibration factor is provided, PRI images are converted
              to radiometrically calibrated ground-range intensity images in float format
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_ASAR',
             ASAR_file, output_name, K_dB], logpath=logpath)


def par_ASF_91(CEOS_leader, CEOS_trailer, SLC_par, logpath=None):
    """
    | SLC parameter file for data data from theAlaska SAR Facility (1991-1996)
    | Copyright 2008, Gamma Remote Sensing, v3.3 25-Mar-2008 clw/uw

    Parameters
    ----------
    CEOS_leader:
        (input) ASF CEOS leader file
    CEOS_trailer:
        (input) ASF CEOS trailer file
    SLC_par:
        (output) ISP SLC image parameter file
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_ASF_91',
             CEOS_leader, CEOS_trailer, SLC_par], logpath=logpath)


def par_ASF_96(CEOS_SAR_leader, SLC_par, logpath=None):
    """
    | ISP parameter file for ASF data 1996-->present v1.1
    | Copyright 2003, Gamma Remote Sensing, v1.4 4-Aug-2003 clw/uw

    Parameters
    ----------
    CEOS_SAR_leader:
        (input) CEOS SAR leader file
    SLC_par:
        (output) ISP SLC parameter file (example <orbit>.slc.par)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_ASF_96',
             CEOS_SAR_leader, SLC_par], logpath=logpath)


def par_ASF_PRI(CEOS_leader, CEOS_data, GRD_par, GRD, logpath=None):
    """
    | ISP parameter file for ASF detected ground range images (L1) Sep 1996 --> present
    | Copyright 2014, Gamma Remote Sensing, v1.3 3-Apr-2014 clw/uw

    Parameters
    ----------
    CEOS_leader:
        (input) CEOS leader file
    CEOS_data:
        (input) CEOS data file binary)
    GRD_par:
        (output) ISP ground range image parameter file
    GRD:
        (output) ISP ground range image (enter -  for none, float intensity)
            * NOTE: The input data converted to intensity using the expression:  (dn/1000.)\*\*2

    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_ASF_PRI',
             CEOS_leader, CEOS_data, GRD_par, GRD], logpath=logpath)


def par_ASF_RSAT_SS(CEOS_leader, CEOS_data, GRD_par, GRD, logpath=None):
    """
    | ISP parameter file for ASF Radarsat-1 SCANSAR images
    | Copyright 2004, Gamma Remote Sensing, v1.0 27-Aug-2004 clw/uw

    Parameters
    ----------
    CEOS_leader:
        (input) CEOS leader file (Radarsat-1 SCANSAR)
    CEOS_data:
        (input) CEOS data file (Radarsat-1 SCANSAR)
    GRD_par:
        (output) ISP image parameter file (example <orbit>.mli.par)
    GRD:
        (output) ISP image (example <orbit>.mli) (enter -  for none, short integer)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_ASF_RSAT_SS',
             CEOS_leader, CEOS_data, GRD_par, GRD], logpath=logpath)


def par_ATLSCI_ERS(CEOS_SAR_leader, CEOS_Image, SLC_par, logpath=None):
    """
    | ISP parameter file for ATL-SCI ERS SLC data
    | Copyright 2003, Gamma Remote Sensing, v2.8 24-Nov-2003 clw

    Parameters
    ----------
    CEOS_SAR_leader:
        (input) CEOS SAR leader file (LEA_01.001)
    CEOS_Image:
        (input) CEOS image data segment (DAT_01.001)
    SLC_par:
        (output) ISP SLC parameter file (example <orbit>.slc.par)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_ATLSCI_ERS',
             CEOS_SAR_leader, CEOS_Image, SLC_par], logpath=logpath)


def par_CS_SLC(HDF5, trunk, logpath=None):
    """
    | Generate ISP SLC parameter and image files for Cosmo-Skymed SCS data
    | Copyright 2015, Gamma Remote Sensing, v1.7 21-Aug-2015 awi/ms/cw

    Parameters
    ----------
    HDF5:
        (input) SCS data file in HDF5 format
    trunk:
        (output) output file name trunk used for output filenames 
            (example: yyyymmdd -> yyyymmdd_pol_beamid.slc yyyymmdd_pol_beamid.slc.par)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_CS_SLC',
             HDF5, trunk], logpath=logpath)


def par_CS_SLC_TIF(GeoTIFF, XML, trunk, logpath=None):
    """
    | Generate ISP SLC parameter and image files for Cosmo Skymed SCS data in GeoTIFF format
    | Copyright 2015, Gamma Remote Sensing, v1.4 12-Aug-2015 awi/ms/clw

    Parameters
    ----------
    GeoTIFF:
        (input) SCS data file in GeoTIFF format
    XML:
        (input) SCS meta data file in XML format
    trunk:
        (output) output file name trunk used for output filenames 
            (example: yyyymmdd -> yyyymmdd_pol_beamid.slc yyyymmdd_pol_beamid.slc.par)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_CS_SLC_TIF',
             GeoTIFF, XML, trunk], logpath=logpath)


def par_EORC_PALSAR(CEOS_leader, SLC_par, CEOS_data, SLC='-', logpath=None):
    """
    | Reformat EORC PALSAR + PALSAR2 level 1.1 CEOS format SLC data and generate the ISP parameter file
    | Copyright 2016, Gamma Remote Sensing, v2.7 27-Apr-2016 clw

    Parameters
    ----------
    CEOS_leader:
        (input) CEOS leader file for PALSAR or PALSAR-2 Level 1.1 SLC data (LED...)
    SLC_par:
        (output) ISP image parameter file (example: yyyymmdd.SLC.par)
    CEOS_data:
        (input)  PALSAR CEOS format Level 1.1 SLC (IMG...)
    SLC:
        (output) reformatted PALSAR SLC (example: yyyymmdd.SLC, enter - for none)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_EORC_PALSAR',
             CEOS_leader, SLC_par, CEOS_data, SLC], logpath=logpath)


def par_ESA_ERS(CEOS_SAR_leader, SLC_par, CEOS_DAT='-', SLC='-', logpath=None):
    """
    | ISP parameter file generation for ERS SLC data from the PGS, VMP, and SPF processors
    | Copyright 2012, Gamma Remote Sensing, v1.4 12-Jan-2012 clw/uw

    Parameters
    ----------
    CEOS_SAR_leader:
        (input) ERS CEOS SAR leader file
    SLC_par:
        (output) ISP SLC parameter file (example: <date>.SLC.par)
    CEOS_DAT:
        (input) CEOS data file (example: DAT_01.001)
    SLC:
        (output) SLC data with file and line headers removed (example: <date>.SLC)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_ESA_ERS',
             CEOS_SAR_leader, SLC_par, CEOS_DAT, SLC], logpath=logpath)


def par_KC_PALSAR_slr(facter_m, CEOS_leader, SLC_par, pol, pls_mode, KC_data, pwr, fdtab='-', logpath=None):
    """
    | Generate ISP parameter file, Doppler table, and images for PALSAR KC Slant-Range data
    | Copyright 2013, Gamma Remote Sensing, v1.9.1 20-Aug-2013 ms,awi,clw

    Parameters
    ----------
    facter_m:
        (input) PALSAR Kyoto-Carbon parameter file
    CEOS_leader:
        (input) PALSAR Kyoto-Carbon leader file (LED)
    SLC_par:
        (output) ISP image parameter file (example: yyyymmdd.mli.par)
    pol:
        polarization e.g. HH or HV
    pls_mode:
        PALSAR acquisition mode:
            * 1: Fine Beam Single
            * 2: Fine Beam Double
            * 3: Wide Beam

    KC_data:
        (input) PALSAR Kyoto-Carbon data (short, little endian, amplitude)
    pwr:
        (output) PALSAR intensity (float, GAMMA Software endianness)
    fdtab:
        (output)table of output polynomials, one polynomial/block used as input to gc_map_fd
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_KC_PALSAR_slr', facter_m,
             CEOS_leader, SLC_par, pol, pls_mode, KC_data, pwr, fdtab], logpath=logpath)


def par_KS_DGM(HDF5, trunk, logpath=None):
    """
    | Generate ISP SLC parameter and PRI image files for Kompsat DGM data
    | Copyright 2014, Gamma Remote Sensing, v1.0 5-May-2014 awi

    Parameters
    ----------
    HDF5:
        (input) DGM data file in HDF5 format
    trunk:
        (output) output file name trunk used for output filenames 
            (example: yyyymmdd -> yyyymmdd_pol_beamid.slc yyyymmdd_pol_beamid.pri.par)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_KS_DGM',
             HDF5, trunk], logpath=logpath)


def par_KS_SLC(HDF5, trunk, logpath=None):
    """
    | Generate ISP SLC parameter and image files for Kompsat SCS data
    | Copyright 2016, Gamma Remote Sensing, v1.4 11-Feb-2016 awi/clw

    Parameters
    ----------
    HDF5:
        (input) SCS data file in HDF5 format
    trunk:
        (output) output file name trunk used for output filenames 
            (example: yyyymmdd -> yyyymmdd_pol_beamid.slc yyyymmdd_pol_beamid.slc.par)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_KS_SLC',
             HDF5, trunk], logpath=logpath)


def par_PRI(CEOS_SAR_leader, PRI_par, CEOS_DAT, PRI, logpath=None):
    """
    | ISP parameter file generation for ERS PRI data from the PGS and VMP processors
    | Copyright 2012, Gamma Remote Sensing, v1.6 12-Jan-2012 clw

    Parameters
    ----------
    CEOS_SAR_leader:
        (input) ERS CEOS SAR leader file for PRI product
    PRI_par:
        (output) ISP image parameter file (example: <yyyymmdd>.PRI.par)
    CEOS_DAT:
        (input) CEOS data file (example: DAT_01.001)
    PRI:
        (output) PRI data with file and line headers removed (example: <yyyymmdd>.PRI)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_PRI',
             CEOS_SAR_leader, PRI_par, CEOS_DAT, PRI], logpath=logpath)


def par_PRI_ESRIN_JERS(CEOS_SAR_leader, PRI_par, CEOS_DAT, PRI, logpath=None):
    """
    | ISP GRD parameter file for ESRIN processed JERS PRI data
    | Copyright 2008, Gamma Remote Sensing, v1.8 16-May-2008 clw/uw

    Parameters
    ----------
    CEOS_SAR_leader:
        (input) ERS CEOS SAR leader file for PRI product
    PRI_par:
        (output) ISP image parameter file (example: <yyyymmdd>.PRI.par)
    CEOS_DAT:
        (input) CEOS data file (example: DAT_01.001)
    PRI:
        (output) PRI data with file and line headers removed (example: <yyyymmdd>.PRI)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_PRI_ESRIN_JERS',
             CEOS_SAR_leader, PRI_par, CEOS_DAT, PRI], logpath=logpath)


def par_PulSAR(CEOS_SAR_leader, SLC_par, logpath=None):
    """
    | ISP parameter file generation for ERS SLC data from the PULSAR SAR processor
    | Copyright 2003, Gamma Remote Sensing, v1.2 4-Aug-2003 clw/uw

    Parameters
    ----------
    CEOS_SAR_leader:
        (input) ERS CEOS SAR leader file
    SLC_par:
        (output) ISP SLC parameter file (example <orbit>.slc.par)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_PulSAR',
             CEOS_SAR_leader, SLC_par], logpath=logpath)


def par_RISAT_GRD(CEOS_leader, BAND_META, GRD_par, CEOS_image, GRD='-', line_dir='-', pix_dir='-', cal_flg='-', KdB='-', logpath=None):
    """
    | Read RISAT-1 Ground-Range data from a CEOS data set and perform radiometric calibration
    | Copyright 2015, Gamma Remote Sensing, v1.2 24-Feb-2015 clw

    Parameters
    ----------
    CEOS_leader:
        (input) CEOS SAR leader file (example: lea_01.001)
    BAND_META:
        (input) BAND_META.txt, additional RISAT system parameters for the scene (format keywork=value)
    GRD_par:
        (output) ISP GRD parameter file (example: YYYYMMDD.GRD.par)
    CEOS_image:
        (input) CEOS Ground-Range image file (example: dat_01.001)
    GRD:
        (output) Ground-Range data with file and line headers removed (enter - for none: example: YYYYMMDD.GRD)
    line_dir:
        set output image line direction (enter - for default):
            * 0: used value derived from CEOS leader file
            * 1: retain input data line direction  (default)
            * -1: reverse input data line direction

    pix_dir:
        set output pixel direction (enter - for default):
            * 0: used value derived from CEOS leader file
            * 1: retain input data pixel direction (default)
            * -1: reverse input data pixel direction

    cal_flg:
        calibration flag (enter - for default):
            * 0: do not apply radiometric calibration
            * 1: apply radiometric calibration including KdB and incidence angle correction (default)

    KdB:
        calibration constant (dB) (enter - to use value in the CEOS leader)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_RISAT_GRD', CEOS_leader,
             BAND_META, GRD_par, CEOS_image, GRD, line_dir, pix_dir, cal_flg, KdB], logpath=logpath)


def par_RISAT_SLC(CEOS_leader, BAND_META, SLC_par, CEOS_image, SLC='-', line_dir='-', pix_dir='-', cal_flg='-', KdB='-', logpath=None):
    """
    | Read RISAT-1 CEOS format SLC data and perform radiometric calibration
    | Copyright 2013, Gamma Remote Sensing, v1.1 3-Jun-2013 clw

    Parameters
    ----------
    CEOS_leader:
        (input) CEOS SAR leader file (example: lea_01.001)
    BAND_META:
        (input) BAND_META.txt, additional RISAT system parameters for the scene (format keywork=value)
    SLC_par:
        (output) ISP SLC image parameter file (example: YYYYMMDD.grd.par)
    CEOS_image:
        (input) CEOS SLC image file (example: dat_01.001)
    SLC:
        (output) SLC data with file and line headers removed (enter - for none: example: YYYYMMDD.grd)
    line_dir:
        set output image line direction (enter - for default):
            * 0: used value derived from CEOS leader file
            * 1: retain input data line direction  (default)
            * -1: reverse input data line direction

    pix_dir:
        set output pixel direction (enter - for default):
            * 0: used value derived from CEOS leader file
            * 1: retain input data pixel direction (default)
            * -1: reverse input data pixel direction

    cal_flg:
        calibration flag (enter - for default):
            * 0: do not apply radiometric calibration
            * 1: apply radiometric calibration including KdB and incidence angle correction (default)

    KdB:
        calibration constant (dB) (enter - to use value in the CEOS leader)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_RISAT_SLC', CEOS_leader,
             BAND_META, SLC_par, CEOS_image, SLC, line_dir, pix_dir, cal_flg, KdB], logpath=logpath)


def par_RSAT2_SG(product_XML, lut_XML, GeoTIFF, polarization, GRD_par, GRD, logpath=None):
    """
    | Generate SLC parameter and ground range image files for Radarsat 2 SGF/SGX data
    | Copyright 2015, Gamma Remote Sensing, v1.8 13-Aug-2015 awi/cw

    Parameters
    ----------
    product_XML:
        (input) Radarsat-2 product annotation XML file (product.xml)
    lut_XML:
        (input) Radarsat-2 calibration XML file (lutSigma.xml), use - for no calibration
    GeoTIFF:
        (input) image data file in GeoTIFF format (imagery_PP.tif)
    polarization:
        (input) image polarization: HH, VV, HV, VH
    GRD_par:
        (output) ISP GRD parameter file (example: yyyymmdd_PP.GRD.par)
    GRD:
        (output) float GRD data file (example: yyyymmdd_pp.GRD)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_RSAT2_SG', product_XML,
             lut_XML, GeoTIFF, polarization, GRD_par, GRD], logpath=logpath)


def par_RSAT2_SLC(product_XML, lut_XML, GeoTIFF, polarization, SLC_par, SLC, logpath=None):
    """
    | Generate SLC parameter and image files for Radarsat 2 SLC data from GeoTIFF
    | Copyright 2015, Gamma Remote Sensing, v2.5 13-Aug-2015 awi/clw

    Parameters
    ----------
    product_XML:
        (input) Radarsat-2 product annotation XML file (product.xml)
    lut_XML:
        (input) Radarsat-2 calibration XML file (lutSigma.xml), use - for no calibration
    GeoTIFF:
        (input) image data file in GeoTIFF format (imagery_PP.tif)
    polarization:
        (input) image polarization: HH, VV, HV, VH
    SLC_par:
        (output) ISP SLC parameter file (example: yyyymmdd_pp.SLC.par)
    SLC:
        (output) SLC data file (example: yyyymmdd_pp.SLC)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_RSAT2_SLC',
             product_XML, lut_XML, GeoTIFF, polarization, SLC_par, SLC], logpath=logpath)


def par_RSAT_SCW(CEOS_leader, CEOS_trailer, CEOS_data, GRD_par, GRD, sc_dB='-', dt='-', logpath=None):
    """
    | ISP parameter file for SCANSAR Wide Swath Data
    | Copyright 2012, Gamma Remote Sensing, v2.0 14-Feb-2012 clw

    Parameters
    ----------
    CEOS_leader:
        (input) CEOS SAR leader file
    CEOS_trailer:
        (input) CEOS SAR trailer file
    CEOS_data:
        (input) CEOS data file binary)
    GRD_par:
        (output) ISP ground range image parameter file (example <orbit>.mli.par)
    GRD:
        (output) ISP ground range image (example <orbit>.mli) (enter -  for none, float)
    sc_dB:
        intensity scale factor in dB (enter - for default:   0.00)
    dt:
        azimuth image time offset (s) (enter - for default = 0.0)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_RSAT_SCW', CEOS_leader,
             CEOS_trailer, CEOS_data, GRD_par, GRD, sc_dB, dt], logpath=logpath)


def par_RSAT_SGF(CEOS_leader, CEOS_data, GRD_par, GRD, sc_dB='-', dt='-', logpath=None):
    """
    | ISP parameter file for RSI/Atlantis Radarsat SGF (ground range) and SCANSAR SCW16 data
    | Copyright 2012, Gamma Remote Sensing, v2.2 14-Feb-2012 clw

    Parameters
    ----------
    CEOS_leader:
        (input) CEOS leader file (RSI SGF or SCW16 products, LEA_01.001)
    CEOS_data:
        (input) CEOS data file (RSI SGF or SCW16 products, DAT_01.001) 
    GRD_par:
        (output) ISP ground range image parameter file (example <orbit>.mli.par)
    GRD:
        (output) ISP ground range image (example <orbit>.GRD.par) (enter -  for none, float)
    sc_dB:
        intensity scale factor in dB (enter - for default:   0.00)
    dt:
        azimuth image time offset (s) (enter - for default = 0.0)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_RSAT_SGF',
             CEOS_leader, CEOS_data, GRD_par, GRD, sc_dB, dt], logpath=logpath)


def par_RSAT_SLC(CEOS_leader, SLC_par, CEOS_data, SLC='-', sc_dB='-', dt='-', logpath=None):
    """
    | ISP parameter file for RSI/Atlantis/ASF processed Radarsat SLC data
    | Copyright 2012, Gamma Remote Sensing, v4.0 5-Sep-2012 clw

    Parameters
    ----------
    CEOS_leader:
        (input) CEOS SAR leader file (example: lea_01.001)
    SLC_par:
        (output) ISP SLC parameter file (example: <date>.SLC.par)
    CEOS_data:
        (input) CEOS data file (example: dat_01.001)
    SLC:
        (output) SLC data with file and line headers removed (example: <date>.SLC)
    sc_dB:
        intensity scale factor in dB (enter - for default:  60.00)
    dt:
        azimuth image time offset (s) (enter - for default = 0.0)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_RSAT_SLC',
             CEOS_leader, SLC_par, CEOS_data, SLC, sc_dB, dt], logpath=logpath)


def par_RSI_ERS(CEOS_SAR_leader, SLC_par, logpath=None):
    """
    | ISP parameter file for RSI processed ERS SLC data
    | Copyright 2003, Gamma Remote Sensing, v1.7 4-Aug-2003 clw/uw

    Parameters
    ----------
    CEOS_SAR_leader:
        (input) ERS CEOS SAR leader file
    SLC_par:
        (output) ISP SLC parameter file (example <orbit>.slc.par)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_RSI_ERS',
             CEOS_SAR_leader, SLC_par], logpath=logpath)


def par_S1_GRD(GeoTIFF, annotation_XML, calibration_XML, noise_XML, MLI_par, MLI, GRD_par='-', GRD='-', eflg='-', rps='-', noise_pwr='-', logpath=None):
    """
    | Generate MLI and GRD images and parameter files from a Sentinel-1 GRD product
    | Copyright 2016, Gamma Remote Sensing, v2.8 17-Aug-2016 awi/clw/ts

    Parameters
    ----------
    GeoTIFF:
        (input) image data file in GeoTIFF format (\*.tiff)
    annotation_XML:
        (input) Sentinel-1 L1 XML annotation file
    calibration_XML:
        (input) Sentinel-1 L1 radiometric calibration XML file (enter - for no radiometric calibration)
    noise_XML:
        (input) Sentinel-1 L1 noise XML file (enter - to not add back thermal noise)
            * NOTE: The L1 GRD product has thermal noise subtracted, enter noise_XML to add back thermal noise

    MLI_par:
        (output) MLI parameter file (example: yyyymmdd_pp.MLI.par)
    MLI:
        (output) MLI data file in slant range geometry (example: yyyymmdd_pp.MLI)
    GRD_par:
        (output) GRD parameter file (example: yyyymmdd_pp.GRD.par, enter - for none)
    GRD:
        (output) GRD data file (example: yyyymmdd_pp.GRD, enter - for none)
    eflg:
        GR-SR grid extrapolation flag:
            * 0: no extrapolation of the GR-SR grid beyond the grid boundaries
            * 1: permit extrapolation of the GR-SR grid to cover the entire image (default)
            * NOTE: extrapolation of the GR-SR grid may introduce geocoding errors

    rps:
        slant range pixel spacing (m) (enter - for default: calculated from ground-range parameters)
    noise_pwr:
        noise intensity for each MLI sample in slant range using data from noise_XML
            * NOTE:  when the noise_pwr file is specified, noise power correction will NOT be applied to the MLI data values

    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_S1_GRD', GeoTIFF, annotation_XML,
             calibration_XML, noise_XML, MLI_par, MLI, GRD_par, GRD, eflg, rps, noise_pwr], logpath=logpath)


def par_S1_SLC(GeoTIFF, annotation_XML, calibration_XML, noise_XML, SLC_par, SLC, TOPS_par='-', dtype='-', sc_dB='-', noise_pwr='-', logpath=None):
    """
    | Generate SLC parameter and image files for Sentinel-1 SLC data
    | Copyright 2016, Gamma Remote Sensing, v3.3 17-Aug-2016 awi/clw

    Parameters
    ----------
    GeoTIFF:
        (input) image data file in GeoTIFF format (\*.tiff)
    annotation_XML:
        (input) Sentinel-1 L1 XML annotation file
    calibration_XML:
        (input) Sentinel-1 L1 radiometric calibration XML file (enter - for no radiometric calibration)
    noise_XML:
        (input) Sentinel-1 L1 noise XML file (enter - to not subtract thermal noise power level)
    SLC_par:
        (output) ISP SLC parameter file (example: yyyymmdd_iw1_vv.SLC.par)
    SLC:
        (output) SLC data file (example: yyyymmdd_iw1_vv.SLC)
    TOPS_par:
        (output) SLC burst annotation file, TOPS and EW SLC data only (enter - for none, example: yyyymmdd_iw1_vv.TOPS_par)
    dtype:
        output data type:
            * 0: FCOMPLEX (default)
            * 1: SCOMPLEX

    sc_dB:
        scale factor for FCOMPLEX -> SCOMPLEX, (enter - for default: HH,VV (dB): 60.0000,  VH,HV: 70.0000)
    noise_pwr:
        noise intensity for each SLC sample in slant range using data from noise_XML
            * NOTE:  when the noise_pwr file is specified, noise power will NOT be subtracted from the image data values

    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_S1_SLC', GeoTIFF, annotation_XML,
             calibration_XML, noise_XML, SLC_par, SLC, TOPS_par, dtype, sc_dB, noise_pwr], logpath=logpath)


def par_TX_SLC(annotation_XML, COSAR, SLC_par, SLC, pol='-', logpath=None):
    """
    | Generate SLC parameter file and SLC image from a Terrasar-X SSC data set
    | Copyright 2016, Gamma Remote Sensing, v2.3 26-May-2016 awi/clw

    Parameters
    ----------
    annotation_XML:
        (input) TerraSAR-X product annotation XML file
    COSAR:
        (input) COSAR SSC strip-mode SLC data file
    SLC_par:
        (output) ISP SLC parameter file (example: yyyymmdd.SLC.par)
    SLC:
        (output) SLC data file, example: yyyymmdd.SLC (enter - for none, SLC output will not be produced)
    pol:
        polarisation HH, HV, VH, VV (default: first polarisation found in the annotation_XML)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/par_TX_SLC',
             annotation_XML, COSAR, SLC_par, SLC, pol], logpath=logpath)


def ph_slope_base(int_in, SLC_par, OFF_par, base, int_out, int_type='-', inverse='-', logpath=None):
    """
    | Subtract/add interferogram flat-Earth phase trend as estimated from initial baseline
    | Copyright 2006, Gamma Remote Sensing, v4.4 3-Nov-2006 clw

    Parameters
    ----------
    int_in:
        (input) interferogram (FCOMPLEX) or unwrapped phase (FLOAT) (unflattened)
    SLC_par:
        (input) ISP parameter file for the reference SLC
    OFF_par:
        (input) ISP offset/interferogram parameter file
    base:
        (input) baseline file
    int_out:
        (output) interferogram (FCOMPLEX) or unwrapped phase (FLOAT) with phase trend subtracted/added
    int_type:
        interferogram type: 0=unwrapped phase, 1=complex interf. (default=1)
    inverse:
        subtract/add inversion flag (0=subtract phase ramp, 1=add phase ramp (default=0)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/ph_slope_base', int_in,
             SLC_par, OFF_par, base, int_out, int_type, inverse], logpath=logpath)


def phase_slope(interf, slopes, width, win_sz='-', thres='-', xmin='-', xmax='-', ymin='-', ymax='-', logpath=None):
    """
    | Calculate interferogram phase slopes in range and azimuth
    | Copyright 2011, Gamma Remote Sensing, v1.3 19-Apr-2011 clw/uw

    Parameters
    ----------
    interf:
        (input) interferogram (fcomplex)
    slopes:
        (output) range and azimuth phase slopes (fcomplex)
    width:
        number of samples/row
    win_sz:
        size of region used for slopes determination (default = 5)
    thres:
        correlation threshold for accepting slope estimates 0.0 -> 1.0 (default=.4)
    xmin:
        starting range pixel offset (default = 0)
    xmax:
        last range pixel offset (default = width-1)
    ymin:
        starting azimuth row offset (default = 0)
    ymax:
        last azimuth row offset (default = nlines-1)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/phase_slope', interf,
             slopes, width, win_sz, thres, xmin, xmax, ymin, ymax], logpath=logpath)


def PRC_vec(SLC_par, PRC, nstate='-', logpath=None):
    """
    | State vectors from ERS PRC orbit data for ISP processing clw/uw
    | Copyright 2008, Gamma Remote Sensing, v1.7 clw 11-Jun-2008

    Parameters
    ----------
    SLC_par:
        (input/output) ISP SLC/MLI image parameter file
    PRC:
        (input) PRC state vector file
    nstate:
        number of state vectors (default=5, maximum=64)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/PRC_vec',
             SLC_par, PRC, nstate], logpath=logpath)


def ptarg_cal_MLI(MLI_par, MLI, r_samp, az_samp, psigma, c_r_samp, c_az_samp, ptr_image, r_plot, az_plot, pcal, osf='-', win='-', pltflg='-', psz='-', csz='-', theta_inc='-', logpath=None):
    """
    | Point target analysis and radiometric calibration of slant-range and ground-range (GRD) images
    | Copyright 2016, Gamma Remote Sensing, v2.6 19-Feb-2016 clw

    Parameters
    ----------
    MLI_par:
        (input) slant-range or ground-range image parameter file for detected intensity data
    MLI:
        (input) ground-range or slant range detected image in FLOAT format
    r_samp:
        point target range sample number, target region size is 16x16
    az_samp:
        point target azimuth line number, target region size is 16x16
    psigma:
        radar cross-section of the calibration target in m\*\*2
    c_r_samp:
        clutter region center range sample number, clutter region size is 16x16
    c_az_samp:
        clutter region center azimuth line number, clutter region size is 16x16
    ptr_image:
        (output) oversampled point target image, with and without phase gradient, nominal width: 256
    r_plot:
        (output) range point target response plot data (text format)
    az_plot:
        (output) azimuth point target response plot data (text format)
    pcal:
        (output) measured point target parameters and radiometric calibration factor (text format)
    osf:
        image over-sampling factor, 2, 4, 8, 16, 32, 64 (enter - for default: 16)
    win:
        maximum search window offset (samples) (enter - for default: 1)
    pltflg:
        plotting mode flag:
            * 0: none
            * 1: output plots in PNG format (default)
            * 2: screen output
            * 3: output plots in PDF format

    psz:
        point target region size (samples) (enter - for default: 16)
    csz:
        clutter region size (samples) (enter - for default: 16)
    theta_inc:
        incidence angle required for calibration of terrain corrrected RISAT-1 images
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/ptarg_cal_MLI', MLI_par, MLI, r_samp, az_samp, psigma,
             c_r_samp, c_az_samp, ptr_image, r_plot, az_plot, pcal, osf, win, pltflg, psz, csz, theta_inc], logpath=logpath)


def ptarg_cal_SLC(SLC_par, SLC, r_samp, az_samp, psigma, c_r_samp, c_az_samp, ptr_image, r_plot, az_plot, pcal, osf='-', win='-', pltflg='-', psz='-', csz='-', c_image='-', logpath=None):
    """
    | Point target analysis and radiometric calibration of SLC images
    | Copyright 2016, Gamma Remote Sensing, v2.4 19-Feb-2016 clw

    Parameters
    ----------
    SLC_par:
        (input) SLC image parameter file
    SLC:
        (input) SLC image in FCOMPLEX or SCOMPLEX format
    r_samp:
        point target range sample number, target region size is 16x16
    az_samp:
        point target azimuth line number, target region size is 16x16
    psigma:
        radar cross-section of the calibration target in m\*\*2
    c_r_samp:
        clutter region center range sample number, clutter region size is 16x16
    c_az_samp:
        clutter region center azimuth line number, clutter region size is 16x16
    ptr_image:
        (output) oversampled point target image, with and without phase gradient, nominal width: 256
    r_plot:
        (output) range point target response plot data (text format)
    az_plot:
        (output) azimuth point target response plot data (text format)
    pcal:
        (output) measured point target parameters and radiometric calibration factor (text format)
    osf:
        image over-sampling factor, 2, 4, 8, 16, 32, 64 (enter - for default: 16)
    win:
        maximum search window offset (samples) (enter - for default: 1)
    pltflg:
        plotting mode flag:
            * 0: none
            * 1: output plots in PNG format (default)
            * 2: screen output
            * 3: output plots in PDF format

    psz:
        point target region size (samples) (enter - for default: 16)
    csz:
        clutter region size (samples) (enter - for default: 16)
    c_image:
        (output) clutter region image (FCOMPLEX format)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/ptarg_cal_SLC', SLC_par, SLC, r_samp, az_samp, psigma,
             c_r_samp, c_az_samp, ptr_image, r_plot, az_plot, pcal, osf, win, pltflg, psz, csz, c_image], logpath=logpath)


def ptarg_SLC(SLC_par, SLC, r_samp, az_samp, ptr_image, r_plot, az_plot, ptr_par='-', osf='-', win='-', pltflg='-', logpath=None):
    """
    | Point target response analysis and interpolation for SLC images
    | Copyright 2016, Gamma Remote Sensing, v1.9 19-Feb-2016 clw

    Parameters
    ----------
    SLC_par:
        (input) SLC image parameter file
    SLC:
        (input) SLC image in FCOMPLEX or SCOMPLEX format
    r_samp:
        point target range sample number
    az_samp:
        point target azimuth line number
    ptr_image:
        (output) oversampled point target image (fcomplex, 1024x1024 samples), with and without phase gradient
    r_plot:
        (output) range point target response plot data (text format)
    az_plot:
        (output) azimuth point target response plot data (text format)
    ptr_par:
        (output) measured point target parameters (text format)
    osf:
        image over-sampling factor, 2, 4, 8, 16, 32, 64 (enter - for default: 16)
    win:
        maximum search window offset (samples) (enter - for default: 1)
    pltflg:
        plotting mode flag:
            * 0: none
            * 1: output plots in PNG format (default)
            * 2: screen output
            * 3: output plots in PDF format

    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/ptarg_SLC', SLC_par, SLC, r_samp,
             az_samp, ptr_image, r_plot, az_plot, ptr_par, osf, win, pltflg], logpath=logpath)


def radcal_MLI(MLI, MLI_PAR, OFF_par, CMLI, antenna='-', rloss_flag='-', ant_flag='-', refarea_flag='-', sc_dB='-', K_dB='-', pix_area='-', logpath=None):
    """
    | Radiometric calibration for multi-look intensity (MLI) data
    | Copyright 2016, Gamma Remote Sensing, v2.0 9-Nov-2016 uw/clw/of

    Parameters
    ----------
    MLI:
        (input) MLI image (float)
    MLI_PAR:
        (input) SLC parameter file of input MLI image
    OFF_par:
        (input) ISP offset/interferogram parameter file (enter - for images in MLI geometry)
    CMLI:
        (output) radiometrically calibrated output MLI (float)
    antenna:
        (input) 1-way antenna gain pattern file or - if not provided
    rloss_flag:
        range spreading loss correction:
            * 0: no correction (default)
            * 1: apply r^3 correction  (all modes except ASAR APS)
            * 2: apply r^4 correction (used only for ASAR APS mode)
            * -1: undo r^3 correction
            * -2: undo r^4 correction)

    ant_flag:
        antenna pattern correction:
            * 0: no correction (default)
            * 1: apply antenna pattern correction
            * -1: undo antenna pattern correction)

    refarea_flag:
        reference pixel area correction:
            * 0: no pixel area correction (default)
            * 1: calculate sigma0, scale area by sin(inc_ang)/sin(ref_inc_ang)
            * 2: calculate gamma0, scale area by sin(inc_ang)/(cos(inc_ang)\*sin(ref_inc_ang)
            * -1: undo sigma0 area scaling factor
            * -2: undo gamma0 area scaling factor

    sc_dB:
        scale factor in dB (default: 0.0)
    K_dB:
        calibration factor in dB (default: -(value from MLI_PAR))
    pix_area:
        (output) ellipsoid-based ground range sigma0 or gamma0 pixel reference area (float)
            refarea_flag 1 or -1: sigma0 ref. area
            refarea_flag 2 or -2: gamma0 ref. area
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/radcal_MLI', MLI, MLI_PAR, OFF_par, CMLI,
             antenna, rloss_flag, ant_flag, refarea_flag, sc_dB, K_dB, pix_area], logpath=logpath)


def radcal_PRI(PRI, PRI_PAR, GRD, GRD_PAR, K_dB='-', inc_ref='-', roff='-', nr='-', loff='-', nl='-', logpath=None):
    """
    | Convert ESA processed short integer format PRI to radiometrically calibrated GRD image (float)
    | Copyright 2016, Gamma Remote Sensing, v1.5 5-Mar-2016 uw/clw

    Parameters
    ----------
    PRI:
        (input) PRI ground-range image (short integer, sqrt(backscat. intensity)
    PRI_PAR:
        (input) SLC parameter file of input PRI ground-range image (yyyymmdd.PRI.par)
    GRD:
        (output) calibrated ground-range image (float, backscat. intensity)
    GRD_PAR:
        (output) ISP image parameter file of output calibrated ground-range image (yyyymmdd.GRD.par)
    K_dB:
        calibration factor in decibels (default: 59.75 dB)
            ERS1 (D-Paf,ESRIN): 58.24 dB, ERS2 (D-Paf,ESRIN,I-Paf,UK-Paf after 1997): 59.75 dB
            ENVISAT ASAR: 55.0 dB (all modes)
            for details see product specifications and ESA publications.
    inc_ref:
        reference incidence angle in deg. (default: 23.0 deg.)
            ENVISAT ASAR: 90.0 deg. (all modes)
    roff:
        offset to starting range sample (default: 0)
    nr:
        number of range samples (default: 0, to end of line)
    loff:
        offset to starting line (default: 0, 1 header line in the input file is assumed for ERS)
    nl:
        number of lines to copy (default: 0, to end of file)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/radcal_PRI', PRI, PRI_PAR,
             GRD, GRD_PAR, K_dB, inc_ref, roff, nr, loff, nl], logpath=logpath)


def radcal_pwr_stat(SLC_tab, SLC_tab_cal, plist, MSR_cal, PWR_cal, roff='-', loff='-', nr='-', nl='-', plist_out='-', logpath=None):
    """
    | Generate calibrated SLC image files using point targets determined from the Mean/Sigma Ratio and Intensity
    | Copyright 2012, Gamma Remote Sensing, v1.3 11-May-2012 clw/uw

    Parameters
    ----------
    SLC_tab:
        (input) two column list of the SLC filenames and SLC parameter filenames of the uncalibrated SLC images
    SLC_tab_cal:
        (input) two column list of the SLC filenames and SLC parameter filenames of the calibrated SLC images (enter - for none)
    plist:
        (input) point list for the point to use for calibraton (int, enter - to use the data to determine the calibration points)
    MSR_cal:
        mean/sigma ratio for point target selection for relative calibration between scenes:    1.500
    PWR_cal:
        intensity threshold ratio for point target selection for relative calibration between scenes:    1.000
    roff:
        offset to starting range of section to analyze (default -: 0)
    loff:
        offset to starting line of section to analyze (default -: 0)
    nr:
        number of range pixels to analyze (default -: to end of line)
    nl:
        number of azimuth lines to analyze (default -: to end of file)
    plist_out:
        point list of points used to determine calibration using MSR_cal and PWR_cal thresholds
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/radcal_pwr_stat', SLC_tab,
             SLC_tab_cal, plist, MSR_cal, PWR_cal, roff, loff, nr, nl, plist_out], logpath=logpath)


def radcal_SLC(SLC, SLC_PAR, CSLC, CSLC_PAR, fcase='-', antenna='-', rloss_flag='-', ant_flag='-', refarea_flag='-', sc_dB='-', K_dB='-', pix_area='-', logpath=None):
    """
    | Radiometric calibration of SLC data
    | Copyright 2016, Gamma Remote Sensing, v2.3 9-Nov-2016 uw/clw/of

    Parameters
    ----------
    SLC:
        (input) SLC (fcomplex or scomplex)
    SLC_PAR:
        (input) SLC parameter file of input SLC
    CSLC:
        (output) radiometrically calibrated SLC (fcomplex or scomplex)
    CSLC_PAR:
        (output) SLC parameter file of output calibrated SLC
    fcase:
        format case (default = 1)
            * 1: fcomplex --> fcomplex (pairs of float)
            * 2: fcomplex --> scomplex (pairs of short integer)
            * 3: scomplex --> fcomplex
            * 4: scomplex --> scomplex

    antenna:
        1-way antenna gain pattern file or - (if not provided)
    rloss_flag:
        range spreading loss correction:
            * 0: no correction (default)
            * 1: apply r^3 correction  (all modes except ASAR APS)
            * 2: apply r^4 correction (used only for ASAR APS mode)
            * -1: undo r^3 correction
            * -2: undo r^4 correction)

    ant_flag:
        antenna pattern correction:
            * 0: no correction (default)
            * 1: apply antenna pattern correction
            * -1: undo antenna pattern correction)

    refarea_flag:
        reference pixel area correction:
            * 0: no pixel area correction (default)
            * 1: calculate sigma0, scale area by sin(inc_ang)/sin(ref_inc_ang)
            * 2: calculate gamma0, scale area by sin(inc_ang)/(cos(inc_ang)\*sin(ref_inc_ang)
            * -1: undo sigma0 area scaling factor
            * -2: undo gamma0 area scaling factor

    sc_dB:
        scale factor in dB (default: 0.0)
    K_dB:
        calibration factor in dB (default: -(value from SLC_PAR) )
    pix_area:
        (output) ellipsoid-based ground range sigma0 or gamma0 pixel reference area (float)
            refarea_flag 1 or -1: sigma0 ref. area
            refarea_flag 2 or -2: gamma0 ref. area
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/radcal_SLC', SLC, SLC_PAR, CSLC, CSLC_PAR,
             fcase, antenna, rloss_flag, ant_flag, refarea_flag, sc_dB, K_dB, pix_area], logpath=logpath)


def rascc_mask(cc, pwr, width, start_cc='-', start_pwr='-', nlines='-', pixavr='-', pixavaz='-', cc_thres='-', pwr_thres='-', cc_min='-', cc_max='-', scale='-', exp='-', LR='-', rasf='-', logpath=None):
    """
    | Generate phase unwrapping validity mask using correlation and intensity
    | Copyright 2016, Gamma Remote Sensing, v2.0 12-Sep-2016 clw/uw

    Parameters
    ----------
    cc:
        (input)interferometric correlation image (float)
    pwr:
        (input)intensity image (float, enter - if not available)
    width:
        number of samples/row
    start_cc:
        starting line of coherence image (default: 1)
    start_pwr:
        starting line of intensity image (default: 1)
    nlines:
        number of lines to display (default=0: to end of file)
    pixavr:
        number of pixels to average in range (default: 1)
    pixavaz:
        number of pixels to average in azimuth (default: 1)
    cc_thres:
        coherence threshold for masking, pixels with cc < cc_thres are set to 0 (default: 0.0)
    pwr_thres:
        relative intensity threshold for masking, pixels with intensity < pwr_thres \* average intensity are set to 0 (default: 0)
    cc_min:
        minimum coherence value used for color display (default: 0.1)
    cc_max:
        maximum coherence value used for color display (default: 0.9)
    scale:
        intensity display scale factor (default: 1.)
    exp:
        intensity display exponent (default: .35)
    LR:
        left/right mirror image flag, (1: normal (default), -1: mirror image)
    rasf:
        (output) image filename, extension determines the format, enter - for default: \*.ras
            \*.bmp BMP format
            \*.ras Sun raster format
            \*.tif TIFF format
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/rascc_mask', cc, pwr, width, start_cc, start_pwr,
             nlines, pixavr, pixavaz, cc_thres, pwr_thres, cc_min, cc_max, scale, exp, LR, rasf], logpath=logpath)


def rascc_mask_thinning(ras_in, in_file, width, ras_out, nmax='-', thresh_1='-', thresh_nmax='-', logpath=None):
    """
    | Adaptive sampling reduction for phase unwrapping validity mask
    | Copyright 2015, Gamma Remote Sensing, v1.5 5-Dec-2015 uw/clw

    Parameters
    ----------
    ras_in:
        (input) validity mask (SUN/BMP/TIFF raster format 8-bit image)
    in_file:
        (input) file used for adaptive sampling reduction, e.g. correlation coefficient (float)
    width:
        number of samples/row of in_file
    ras_out:
        (output) validity mask with reduced sampling (8-bit SUN rasterfile or BMP format image)
    nmax:
        number of sampling reduction runs (default: 3)
    thresh_1:
        first threshold (used for smallest scale sampling reduction)
            ...          further thresholds
    thresh_nmax:
        threshold nmax (used for largest scale sampling reduction)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/rascc_mask_thinning', ras_in,
             in_file, width, ras_out, nmax, thresh_1, thresh_nmax], logpath=logpath)


def residue(int, flag, width, xmin='-', xmax='-', ymin='-', ymax='-', logpath=None):
    """
    | Determine interferometric phase unwrapping residues
    | Copyright 2014, Gamma Remote Sensing, v2.6 14-Jan-2014 clw/uw

    Parameters
    ----------
    int:
        (input) interferogram (fcomplex)
    flag:
        (input) flag file (unsigned char)
    width:
        number of samples/row
    xmin:
        offset to starting range pixel(default = 0)
    xmax:
        offset last range pixel (default = width-1)
    ymin:
        offset to starting azimuth row (default = 0)
    ymax:
        offset to last azimuth row (default = nlines-1)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/residue', int,
             flag, width, xmin, xmax, ymin, ymax], logpath=logpath)


def residue_cc(int, flag, width, xmin='-', xmax='-', ymin='-', ymax='-', logpath=None):
    """
    | Determine interferometric phase unwrapping residues considering low coherence regions
    | Copyright 2014, Gamma Remote Sensing, v2.6  20-Jan-2014 clw/uw/ts

    Parameters
    ----------
    int:
        (input) interferogram (fcomplex)
    flag:
        (input) flag file (unsigned char)
    width:
        number of samples/row
    xmin:
        offset to starting range pixel(default = 0)
    xmax:
        offset last range pixel (default = width-1)
    ymin:
        offset to starting azimuth row (default = 0)
    ymax:
        offset to last azimuth row (default = nlines-1)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/residue_cc',
             int, flag, width, xmin, xmax, ymin, ymax], logpath=logpath)


def RSAT2_vec(SLC_par, RSAT2_orb, nstate='-', logpath=None):
    """
    | Extract Radarsat-2 state vectors from a definitive orbit file
    | Copyright 2010, Gamma Remote Sensing, v1.0 clw 13-May-2010

    Parameters
    ----------
    SLC_par:
        (input) ISP image parameter file
    RSAT2_orb:
        Radarsat-2 definitive orbit data file available from MDA. (orbit_number_def.orb)
    nstate:
        number of state vectors to extract (enter - for default: 9)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/RSAT2_vec',
             SLC_par, RSAT2_orb, nstate], logpath=logpath)


def S1_burstloc(annotation_XML, logpath=None):
    """
    | Print Burst information found in the Sentinel-1 annotation file
    | Copyright 2016, Gamma Remote Sensing, v1.0 22-Jan-2016 awi

    Parameters
    ----------
    annotation_XML:
        (input) Sentinel-1 L1 XML annotation file
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/S1_burstloc',
             annotation_XML], logpath=logpath)


def S1_OPOD_vec(SLC_PAR, OPOD, nstate='-', logpath=None):
    """
    | Extract Sentinel-1 OPOD state vectors and copy into the ISP image parameter file
    | Copyright 2015, Gamma Remote Sensing, v1.3 17-Aug-2016 awi/clw

    Parameters
    ----------
    SLC_PAR:
        (input/output)ISP SLC/MLI image parameter file
    OPOD:
        (input) Sentinel-1 OPOD orbit data file (AUX_POEORB or AUX_RESORB)
            https://qc.sentinel1.eo.esa.int/aux_resorb/
    nstate:
        number of state vectors to extract (default: include 60 sec extention at the start and end of the SLC data)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/S1_OPOD_vec',
             SLC_PAR, OPOD, nstate], logpath=logpath)


def sbi_filt(SLC_1, SLC1_par, SLC2R_par, SLCf, SLCf_par, SLCb, SLCb_par, norm_sq, iwflg='-', logpath=None):
    """
    | Azimuth filtering of SLC data to support split-beam interferometry to measure azimuth offsets
    | Copyright 2016, Gamma Remote Sensing, v1.2 clw 5-Mar-2016

    Parameters
    ----------
    SLC-1:
        (input) SLC image (SCOMPLEX or FCOMPLEX format)
    SLC1_par:
        (input) SLC image parameter file
    SLC2R_par:
        (input) SLC2 ISP image parameter file for the co-registered image of the interferometric pair,
            used to determine azimuth common-band for each output SLC (enter - for none)
    SLCf:
        (output) SLC image (forward-looking, FCOMPLEX format)
    SLCf_par:
        (output) SLC parameter file (forward-looking)
    SLCb:
        (output) SLC image (backward-looking, FCOMPLEX format)
    SLCb_par:
        (output) SLC parameter file (backward-looking)
    norm_sq:
        squint between beams as a fraction of the azimuth spectrum width (default: 0.5)
    iwflg:
        inverse weighting flag:
            * 0: no compensation for azimuth spectrum weighting
            * 1: compensate for the azimuth spectrum weighting (default)

    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/sbi_filt', SLC_1, SLC1_par,
             SLC2R_par, SLCf, SLCf_par, SLCb, SLCb_par, norm_sq, iwflg], logpath=logpath)


def sbi_offset(sbi_unw, SLCf_par, SLCb_par, OFF_par, az_offset, logpath=None):
    """
    | Calculate azimuth offsets from unwrapped split-beam interferogram
    | Copyright 2011, Gamma Remote Sensing, v1.0 25-Nov-2011

    Parameters
    ----------
    sbi_unw:
        (input) unwrapped phase of split-beam interferogram (float)
    SLCf_par:
        (input) reference SLC parameter file (forward-looking)
    SLCb_par:
        (input) reference SLC parameter file (backward-looking)
    OFF_par:
        (input) offset parameter file
    az_offset:
        (output) azimuth offsets (m)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/sbi_offset',
             sbi_unw, SLCf_par, SLCb_par, OFF_par, az_offset], logpath=logpath)


def slant_range(SLC_par, slr, logpath=None):
    """
    | Calculate slant range for every range sample
    | Copyright 2013, Gamma Remote Sensing v1.1 28-Aug-2013

    Parameters
    ----------
    SLC_par:
        (input) SLC or MLI image parameter file
    slr:
        (output) slant range for every sample in the image (float)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/slant_range',
             SLC_par, slr], logpath=logpath)


def SLC_burst_copy(SLC, SLC_par, TOPS_par, SLC_out, SLC_out_par, burst_num, drflg='-', SLC_par2='-', logpath=None):
    """
    | Copy selected burst from Sentinel-1 TOPS SLC to a file
    | Copyright 2014, Gamma Remote Sensing, v1.3 21-Oct-2014 awi/clw

    Parameters
    ----------
    SLC:
        (input) Sentinel-1 TOPS mode burst SLC
    SLC_par:
        (input) SLC parameter file for the TOPS burst SLC
    TOPS_par:
        (input) TOPS parameter file for the TOPS burst SLC
    SLC_out:
        (output) SLC file containing a single burst
    SLC_out_par:
        (output) SLC parameter file for the single burst SLC
    burst_num:
        burst number of selected burst (1->number of bursts in the SLC)
    drflg:
        deramp phase flag:
            * 0: no modification of the burst SLC phase (default)
            * 1: subtract TOPS Doppler phase ramp (deramp)

    SLC_par2:
        (output) SLC parameter file for the single burst SLC with deramped phase (drflg: 1)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/SLC_burst_copy', SLC, SLC_par,
             TOPS_par, SLC_out, SLC_out_par, burst_num, drflg, SLC_par2], logpath=logpath)


def SLC_burst_corners(SLC_par, TOPS_par, logpath=None):
    """
    | Calculate corner geographic coordinates of Sentinel-1 TOPS SLC bursts
    | Copyright 2016, Gamma Remote Sensing, v1.1 14-Apr-2016 awi/rc/cw

    Parameters
    ----------
    SLC_par:
        (input) SLC parameter file for the TOPS burst SLC
    TOPS_par:
        (input) TOPS parameter file for the TOPS burst SLC
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/SLC_burst_corners',
             SLC_par, TOPS_par], logpath=logpath)


def SLC_cat(SLC_1, SLC_2, SLC1_par, SLC2_par, OFF_par, SLC_3, SLC3_par, dopflg='-', iflg='-', phflg='-', logpath=None):
    """
    | Concatenate two SLC images using 2-D SINC interpolation
    | Copyright 2015, Gamma Remote Sensing, v1.6 11-Nov-2015 clw

    Parameters
    ----------
    SLC-1:
        (input) SLC-1 image
    SLC-2:
        (input) SLC-2 image to be appended to SLC-1
    SLC1_par:
        (input) SLC-1 ISP image parameter file
    SLC2_par:
        (input) SLC-2 ISP image parameter file
    OFF_par:
        (input) ISP offset parameter file containing offset polynomials between SLC-1 and SLC-2
    SLC-3:
        (output) concatenated SLC
    SLC3_par:
        (output) ISP image parameter file for concatenated image
    dopflg:
        Doppler flag:
            * 0: ignore Doppler centroid information, assume 0 Doppler centroid
            * 1: use Doppler centroid information for interpolation (default)

    iflg:
        input data type flag: 
            * 0: input data are SLC images, use data type specified in SLC_par files (SCOMPLEX or FCOMPLEX) (default)
            * 1: input scenes are interferograms, force FCOMPLEX data type

    phflg:
        phase offset correction flag:
            * 0: no phase offset correction for SLC-2
            * 1: apply phase offset correction to SLC-2 (default)

    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/SLC_cat', SLC_1, SLC_2, SLC1_par,
             SLC2_par, OFF_par, SLC_3, SLC3_par, dopflg, iflg, phflg], logpath=logpath)


def SLC_cat_S1_TOPS(SLC_tab1, SLC_tab2, SLC_tab3, logpath=None):
    """
    | Concatenate adjacent Sentinel-1 TOPS SLC images
    | Copyright 2016, Gamma Remote Sensing v1.9 4-Feb-2016

    Parameters
    ----------
    SLC_tab1:
        (input) 3 column list of the reference TOPS SLC swaths in row order IW1, IW2, IW3... (earlier time)
            SLC_tab line entries:   SLC    SLC_par   TOPS_par
    SLC_tab2:
        (input) 3 column list of TOPS SLC-2 swaths in the same order as the SLC_tab1 IW1, IW2, IW3... (later time)
    SLC_tab3:
        (input) 3 column list of the output concatenated TOPS swaths in the order IW1, IW2, IW3...
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/SLC_cat_S1_TOPS',
             SLC_tab1, SLC_tab2, SLC_tab3], logpath=logpath)


def SLC_copy(SLC_in, SLC_par_in, SLC_out, SLC_par_out, fcase='-', sc='-', roff='-', nr='-', loff='-', nl='-', swap='-', header_lines='-', logpath=None):
    """
    | Copy SLC with options for data format conversion, segment extraction, and byte swapping
    | Copyright 2015, Gamma Remote Sensing, v5.1 13-Aug-2015 uw/clw

    Parameters
    ----------
    SLC_in:
        (input) SLC (FCOMPLEX or scOMPLEX format)
    SLC_par_in:
        (input) ISP SLC parameter file for input SLC
    SLC_out:
        (output) selected SLC section (FCOMPLEX or scOMPLEX format)
    SLC_par_out:
        (output) ISP SLC parameter file of output SLC
    fcase:
        data format conversion (enter - for default: output format = input format)
            * 1: FCOMPLEX --> FCOMPLEX (default sc = 1.0)
            * 2: FCOMPLEX --> scOMPLEX (default sc = 10000.0)
            * 3: scOMPLEX --> FCOMPLEX (default sc = 0.0001)
            * 4: scOMPLEX --> scOMPLEX (default sc = 1.0)

    sc:
        scale factor for input SLC data (enter - for default)
    roff:
        offset to starting range sample (enter - for default: 0)
    nr:
        number of range samples (enter - for default: to end of line)
    loff:
        offset to starting line (enter - for default: 0)
    nl:
        number of lines to copy (enter - for default: to end of file)
    swap:
        swap data (enter - for default)
            * 0: normal (default)
            * 1: swap real/imaginary part of complex data
            * 2: swap left/right (near/far range)

    header_lines:
        number of input file header lines (enter - for default: 0)
            * NOTE: CEOS format SLC data have 1 header line
            * NOTE: file offset pointer size (bytes): 8

    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/SLC_copy', SLC_in, SLC_par_in, SLC_out,
             SLC_par_out, fcase, sc, roff, nr, loff, nl, swap, header_lines], logpath=logpath)


def SLC_copy_S1_TOPS(SLC1_tab, SLC2_tab, BURST_tab, dtype='-', logpath=None):
    """
    | Copy multiple bursts from a Sentinel-1 TOPS SLC to an output TOPS SLC
    | Copyright 2016, Gamma Remote Sensing v1.9 16-Sep-2016 clw

    Parameters
    ----------
    SLC1_tab:
        (input) 3 column list of TOPS SLC-1 swaths to be copied in row order IW1, IW2, IW3:
            SLC_tab line entries:   SLC    SLC_par   TOPS_par
    SLC2_tab:
        (input) 3 column list of the output copied SLC-1 TOPS swaths in the order IW1, IW2, IW3
    BURST_tab:
        (input) 2 column list of the first and last burst to copy from each swath, one line for each swath
            BURST_tab line entries: first_burst  last_burst    Note: first burst is 1,  enter - to select last physical burst
            Note: if first_burst <= 0, then blank bursts are generated at the start of the output swath
            if last_burst exceeds the number of bursts in the input data swath, then blank bursts are appended to the end of the output swath
    dtype:
        output data type (default: same as input data):
            * 0: FCOMPLEX
            * 1: SCOMPLEX

    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/SLC_copy_S1_TOPS',
             SLC1_tab, SLC2_tab, BURST_tab, dtype], logpath=logpath)


def SLC_corners(SLC_par, terra_alt='-', logpath=None):
    """
    | Calculate SLC/MLI image corners in geodetic latitude and longitude (deg.)
    | Copyright 2014, Gamma Remote Sensing, v1.6 21-Aug-2014 clw

    Parameters
    ----------
    SLC_par:
        (input) ISP SLC/MLI image parameter file
    terra_alt:
        (input) average terrain altitude (default: 300.000 meters)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/SLC_corners',
             SLC_par, terra_alt], logpath=logpath)


def SLC_deramp(SLC_1, SLC_par1, SLC_2, SLC_par2, mode, dop_ph='-', logpath=None):
    """
    | Calculate and subtract Doppler phase from an SLC image
    | Copyright 2016, Gamma Remote Sensing, v1.5 4-Feb-2016 clw

    Parameters
    ----------
    SLC-1:
        (input) SLC data file (fcomplex or scomplex format)
    SLC_par1:
        (input) SLC parameter file with Doppler information
    SLC-2:
        (output) SLC with Doppler phase removed (or added)
    SLC_par2:
        (output) SLC parameter file for the output SLC
    mode:
        mode of operation:
            * 0: subtract Doppler phase ramp (deramp)
            * 1: add Doppler phase ramp (reramp)

    dop_ph:
        (output) Doppler phase (FLOAT)
            Note: SLC_par1 contains the Doppler polynomial that is used to calculate the Doppler phase ramp
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/SLC_deramp',
             SLC_1, SLC_par1, SLC_2, SLC_par2, mode, dop_ph], logpath=logpath)


def SLC_deramp_S1_TOPS(SLC1_tab, SLC2_tab, mode, phflg, logpath=None):
    """
    | Calculate and subtract S1 TOPS Doppler phase from burst SLC data
    | Copyright 2015, Gamma Remote Sensing v1.4 18-Jun-2015

    Parameters
    ----------
    SLC1_tab:
        (input) 3 column list of TOPS SLC-1 swaths to be deramped in row order IW1, IW2, IW3:
            SLC_tab line entries:   SLC    SLC_par   TOPS_par
    SLC2_tab:
        (input) 3 column list of the output deramped SLC-1 TOPS swaths in the order IW1, IW2, IW3
    mode:
        mode of operation:
            * 0: subtract TOPS Doppler phase (deramp)
            * 1: add Doppler phase ramp (reramp)

    phflg:
        deramp phase flag:
            * 0: do not save TOPS Doppler phase (default)
            * 1: save TOPS Doppler phase, output filename is the same as the deramped SLC with extension .dph

    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/SLC_deramp_S1_TOPS',
             SLC1_tab, SLC2_tab, mode, phflg], logpath=logpath)


def SLC_interp(SLC_2, SLC1_par, SLC2_par, OFF_par, SLC_2R, SLC2R_par, loff='-', nlines='-', logpath=None):
    """
    | SLC complex image resampling using 2-D SINC interpolation
    | Copyright 2015, Gamma Remote Sensing, v4.3 11-Nov-2015 clw

    Parameters
    ----------
    SLC-2:
        (input) SLC-2 image to be resampled to the geometry of the SLC-1 reference image
    SLC1_par:
        (input) SLC-1 ISP image parameter file
    SLC2_par:
        (input) SLC-2 ISP image parameter file
    OFF_par:
        (input) ISP offset/interferogram parameter file
    SLC-2R:
        (output) single-look complex image 2 coregistered to SLC-1
    SLC2R_par:
        (output) SLC-2R ISP image parameter file for coregistered image
    loff:
        offset to first valid output line (in SLC-1 lines) (default: 0)
    nlines:
        number of valid output lines (default: 0, to end of file)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/SLC_interp', SLC_2, SLC1_par,
             SLC2_par, OFF_par, SLC_2R, SLC2R_par, loff, nlines], logpath=logpath)


def SLC_interp_S1_TOPS(SLC2_tab, SLC2_par, SLC1_tab, SLC1_par, OFF_par, SLC2R_tab, SLC_2R='-', SLC2R_par='-', logpath=None):
    """
    | Resample S1 TOPS (IW mode) SLC using global offset polynomial
    | Copyright 2015, Gamma Remote Sensing v1.9 4-Dec-2015

    Parameters
    ----------
    SLC2_tab:
        (input) 3 column list of TOPS SLC-2 swaths to be resampled to the geometry of the reference SLC1 in row order IW1, IW2, IW3:
            SLC_tab line entries:   SLC    SLC_par   TOPS_par
    SLC2_par:
        SLC parameter file of TOPS SLC-2 mosaic, SLC-2 is generated from the TOPS swaths listed in SLC2_tab
    SLC1_tab:
        (input) 3 column list of the reference TOPS SLC swaths in row order IW1, IW2, IW3
    SLC1_par:
        SLC parameter file of the reference TOPS SLC-1 mosaic, SLC-1 is generated from the TOPS swaths listed in SLC1_tab
    OFF_par:
        (input) global ISP offset and interferogram parameter file, the offset model is determined from the TOPS SLC mosaics
    SLC2R_tab:
        (input) 3 column list of the output resampled SLC-2 TOPS swaths in the order IW1, IW2, IW3
    SLC-2R:
        (output) resampled mosaic generated from the swaths listed in SLC2R_tab, coregisted to the TOPS SLC-1 mosaic (enter - for none)
    SLC2R_par:
        (output) SLC parameter file associated with the resampled TOPS SLC-2R mosaic
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/SLC_interp_S1_TOPS', SLC2_tab,
             SLC2_par, SLC1_tab, SLC1_par, OFF_par, SLC2R_tab, SLC_2R, SLC2R_par], logpath=logpath)


def SLC_mosaic_S1_TOPS(SLC_tab, SLC, SLC_par, rlks, azlks, wflg='-', SLCR_tab='-', logpath=None):
    """
    | Calculate SLC mosaic of Sentinel-1 TOPS burst SLC data
    | Copyright 2016, Gamma Remote Sensing v3.5 23-August-2016 clw/awi

    Parameters
    ----------
    SLC_tab:
        (input) 3 column list of SLC, SLC_par, Sentinel-1 TOPS_par sorted in the order IW1, IW2, IW3...
    SLC:
        (output) SLC mosaic image
    SLC_par:
        (output) SLC mosaic image parameter file
    rlks:
        number of range looks used to determine burst window boundaries for the mosaic
    azlks:
        number of azimuth looks used to determine burst window boundaries for the mosaic
    wflg:
        burst window calculation flag:
            * 0: use existing burst window parameters if they exist, otherwise calculate burst window parameters (default)
            * 1: calculate burst window parameters from burst parameters and the number of range and azimuth looks

    SLCR_tab:
        (input) SLC_tab of the reference scene, 3 column list of SLC, SLC_par, TOPS_par sorted sorted in the order IW1, IW2, IW3
            * NOTE: When generating a mosaic of a resampled SLC, the SLC_tab of the reference scene is required

    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/SLC_mosaic_S1_TOPS',
             SLC_tab, SLC, SLC_par, rlks, azlks, wflg, SLCR_tab], logpath=logpath)


def SLC_ovr(SLC, SLC_par, SLC_ovr, SLC_ovr_par, r_ovr, logpath=None):
    """
    | ISP Program /cluster/GAMMA_SOFTWARE-20161207/ISP/bin/SLC_ovr.c
    | Copyright 2016, Gamma Remote Sensing, v1.8 5-Mar-2016 clw
    | Oversample or subsample SLC data in slant-range

    Parameters
    ----------
    SLC:
        (input) SLC file (fcomplex or scomplex)
    SLC_par:
        (input) SLC parameter file of SLC file
    SLC_ovr:
        (output) range resampled SLC file (fcomplex or scomplex)
    SLC_ovr_par:
        (output) SLC parameter file of range resampled SLC data file
    r_ovr:
        integer range oversampling factor (2 --> 16)
            if r_ovr < 0, the SLC will be subsampled, integer range subsampling factor (-2 --> -16)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/SLC_ovr', SLC,
             SLC_par, SLC_ovr, SLC_ovr_par, r_ovr], logpath=logpath)


def SLC_phase_shift(SLC_1, SLC_par1, SLC_2, SLC_par2, ph_shift, logpath=None):
    """
    | Add a constant phase from an SLC image
    | Copyright 2015, Gamma Remote Sensing, v1.1 1-Dec-2015 clw

    Parameters
    ----------
    SLC-1:
        (input) SLC data file (fcomplex or scomplex format)
    SLC_par1:
        (input) SLC parameter file
    SLC-2:
        (output) SLC with phase shift
    SLC_par2:
        (output) SLC parameter file for the output SLC
    ph_shift:
        phase shift to add to SLC phase (radians)
            * NOTE: Used to apply a constant phase shift of -1.25 radians to Sentinel-1 TOPS SLC data
              from swath IW1 acquired up to 10-Mar-2015.
              Used to apply a constant phase shift of -3.83 radians to Sentinel-1 TOPS SLC data with
              H-POL on receive (e.g. VH) acquired up to 10-Mar-2015.
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/SLC_phase_shift',
             SLC_1, SLC_par1, SLC_2, SLC_par2, ph_shift], logpath=logpath)


def split_WB(data_in, data_par_in, data_tab, dtype, logpath=None):
    """
    | ISP: Program /cluster/GAMMA_SOFTWARE-20161207/ISP/bin/split_WB.c
    | Copyright 2011, Gamma Remote Sensing, v1.2 31-May-2011 clw
    | Split WB mosaic image into individual beams using ISP parameter files

    Parameters
    ----------
    data_in:
        (input) input mosaicked data in slant-range geometry (e.g. DEM data)
    data_par_in:
        (input) ISP image parameter file for data in the input mosaic
    data_tab:
        (input) 2 column list of output data filenames and ISP image parameter files for each beam in the mosaic (text)
    dtype:
        (input) input data type:
            * 0: FLOAT
            * 1: FCOMPLEX

    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/split_WB',
             data_in, data_par_in, data_tab, dtype], logpath=logpath)


def SR_to_GRD(MLI_par, OFF_par, GRD_par, in_file, out_file, rlks='-', azlks='-', interp_mode='-', grd_rsp='-', grd_azsp='-', logpath=None):
    """
    | Conversion to ground range for ISP MLI and INSAR data of type float
    | Copyright 2009, Gamma Remote Sensing, v1.9 7-May-2009 uw/clw

    Parameters
    ----------
    MLI_par:
        (input) MLI image parameter file of input slant range image (float)
    OFF_par:
        (input) ISP offset/interferogram parameter file of input image (enter - image in MLI geometry)
    GRD_par:
        (input/output) image parameter file of output ground range image
    in_file:
        (input) slant range image (float)
    out_file:
        (output) ground range image (float)
    rlks:
        multi-looking in range (prior to resampling, default=1)
    azlks:
        multi-looking in azimuth (prior to resampling, default=1)
    interp_mode:
        interpolation mode
            * 0: nearest neighbor (default)
            * 1: spline
            * 2: spline log

    grd_rsp:
        output image ground range sample spacing (m) (default = (input image azimuth spacing) \* azlks)
    grd_azsp:
        output image azimuth sample spacing (m) (default = (input image azimuth spacing) \* azlks)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/SR_to_GRD', MLI_par, OFF_par, GRD_par,
             in_file, out_file, rlks, azlks, interp_mode, grd_rsp, grd_azsp], logpath=logpath)


def subtract_phase(interf_in, phase_file, interf_out, width, factor='-', logpath=None):
    """
    | Land Application Tools: Program /cluster/GAMMA_SOFTWARE-20161207/ISP/bin/subtract_phase.c
    | Copyright 2001, Gamma Remote Sensing, v3.1 23-Jan-2001 uw/clw
    | subtract scaled phase image from a complex interferogram

    Parameters
    ----------
    interf_in:
        (input) input interferogram (fcomplex format)
    phase_file:
        (input) unwrapped interferometric phase (float)
    interf_out:
        (output) output interferogram (input interferogram - scaled phase) (fcomplex)
    width:
        number of samples/line
    factor:
        constant scale factor for input phase data [default=1.0]
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/subtract_phase',
             interf_in, phase_file, interf_out, width, factor], logpath=logpath)


def tree_cc(flag, width, mbl='-', xmin='-', xmax='-', ymin='-', ymax='-', logpath=None):
    """
    | Phase unwrapping tree generation with low correlation search (modified ARW algorithm)
    | Copyright 2014, Gamma Remote Sensing, v2.9 20-Jan-2014 clw/uw

    Parameters
    ----------
    flag:
        (input) phase unwrapping flag file
    width:
        number of samples/row
    mbl:
        maximum branch length (default=32, maximum=64) 
    xmin:
        starting range pixel offset (default = 0)
    xmax:
        last range pixel offset (default = width-1)
    ymin:
        starting azimuth row, relative to start (default = 0)
    ymax:
        last azimuth row, relative to start (default = nlines-1)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/tree_cc', flag,
             width, mbl, xmin, xmax, ymin, ymax], logpath=logpath)


def tree_gzw(flag, width, mbl='-', xmin='-', xmax='-', ymin='-', ymax='-', logpath=None):
    """
    | Phase unwrapping tree generation (GZW algorithm)
    | Copyright 2008, Gamma Remote Sensing, v3.6 5-Sep-2008 clw/uw

    Parameters
    ----------
    flag:
        (input) phase unwrapping flag file
    width:
        number of samples/row
    mbl:
        maximum branch length (default=32)
    xmin:
        starting range pixel offset (default = 0)
    xmax:
        last range pixel offset (default = width-1)
    ymin:
        starting azimuth row, relative to start (default = 0)
    ymax:
        last azimuth row, relative to start (default = nlines-1)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/tree_gzw',
             flag, width, mbl, xmin, xmax, ymin, ymax], logpath=logpath)


def unw_model(interf, unw_model, unw, width, xinit='-', yinit='-', ref_ph='-', width_model='-', logpath=None):
    """
    | Phase unwrapping using a model of the unwrapped phase
    | Copyright 2008, Gamma Remote Sensing, v1.6 5-Sep-2008 clw/uw

    Parameters
    ----------
    interf:
        (input) complex interferogram
    unw_model:
        (input) approximate unwrapped phase model (float)
    unw:
        (output) unwrapped phase (float)
    width:
        number of samples/row of the interferogram
    xinit:
        offset to phase reference location in range (col)
    yinit:
        offset to phase reference location in azimuth (row)
    ref_ph:
        reference point phase (radians) (enter - for phase at the reference point )
    width_model:
        number of samples/row of the unwrapped phase model (default: interferogram width)
    logpath: str or None
        a directory to write command logfiles to
    """
    process(['/cluster/GAMMA_SOFTWARE-20161207/ISP/bin/unw_model', interf,
             unw_model, unw, width, xinit, yinit, ref_ph, width_model], logpath=logpath)
