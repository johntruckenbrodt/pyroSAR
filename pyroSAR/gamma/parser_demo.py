from pyroSAR.gamma.auxil import process


def adapt_filt(int, sm, width, low_SNR_thr='-', filt_width='-', xmin='-', xmax='-', ymin='-', ymax='-', logpath=None, outdir=None, shellscript=None):
    """
    | Adaptive bandpass filtering of interferograms
    | Copyright 2023, Gamma Remote Sensing, v3.6 clw 18-Apr-2023
    
    Parameters
    ----------
    int:
        (input) complex interferogram image filename
    sm:
        (output) smoothed interferogram filename
    width:
        number of samples/row
    low_SNR_thr:
        low SNR threshold (enter - for default: .25);
    filt_width:
        filter width in pixels (enter - for default: 1.0)
    xmin:
        offset to starting range pixel(enter - for default: 0)
    xmax:
        offset last range pixel (enter - for default: width-1)
    ymin:
        offset to starting azimuth row (enter - for default: 0)
    ymax:
        offset to last azimuth row (enter - for default: nlines-1)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/adapt_filt', int, sm, width, low_SNR_thr, filt_width, xmin, xmax, ymin, ymax]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def adf(interf, sm, cc, width, alpha='-', nfft='-', cc_win='-', step='-', loff='-', nlines='-', wfrac='-', logpath=None, outdir=None, shellscript=None):
    """
    | Adaptive interferogram bandpass filter based on the power spectral density
    | Copyright 2024, Gamma Remote Sensing, v3.9 12-Mar-2024 clw/cm
    
    Parameters
    ----------
    interf:
        (input) interferogram (fcomplex)
    sm:
        (output) filtered interferogram (fcomplex)
    cc:
        (output) filtered interferogram correlation coefficient (float)
    width:
        number of samples/line
    alpha:
        exponent for non-linear filtering (enter - for default: 0.40)
    nfft:
        filtering FFT window size, 2\\*\\*N, 8 --> 512, (enter - for default: 32)
    cc_win:
        correlation parameter estimation window size odd, max: 15 (enter - for default: 5)
    step:
        processing step (enter - for default: nfft/8)
    loff:
        offset to starting line to process (enter - for default: 0)
    nlines:
        number of lines to process (enter - for default: to end of file)
    wfrac:
        minimum fraction of points required to be non-zero in the filter window (enter - for default: 0.500)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/adf', interf, sm, cc, width, alpha, nfft, cc_win, step, loff, nlines, wfrac]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def adf2(interf, cc_interf, sm, cc_filt, width, alpha_max='-', nfft='-', cc_win='-', step='-', loff='-', nlines='-', wfrac='-', logpath=None, outdir=None, shellscript=None):
    """
    | Adaptive interferogram filter based on the power spectral density and correlation coefficient
    | Copyright 2023, Gamma Remote Sensing, v1.2 18-Apr-2023 clw/cm
    
    Parameters
    ----------
    interf:
        (input) complex interferogram (fcomplex)
    cc_interf:
        (input) correlation coefficient of the input interferogram (float)
    sm:
        (output) filtered interferogram (fcomplex)
    cc_filt:
        (output) filtered interferogram correlation coefficient (float)
    width:
        number of samples/line
    alpha_max:
        maximum value for the adaptive filter exponent (enter - for default: 0.50)
    nfft:
        filter window FFT size, 2\\*\\*N, 8->512, (enter - for default: 32)
    cc_win:
        filtered interferogram correlation estimation window size odd, max: 21 (enter - for default: 9)
    step:
        processing step in range and azimuth (enter - for default: nfft/8)
    loff:
        offset to starting line to process (enter - for default: 0)
    nlines:
        number of lines to process (enter - for default: to end of file)
    wfrac:
        minimum fraction of points required to be non-zero in the filter window (enter - for default: 0.200)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/adf2', interf, cc_interf, sm, cc_filt, width, alpha_max, nfft, cc_win, step, loff, nlines, wfrac]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def af_SLC(SLC_par, SLC, rwin='-', azwin='-', dr='-', daz='-', thres='-', a1_flg='-', b0_flg='-', offsets='-', n_ovr='-', roff='-', azoff='-', logpath=None, outdir=None, shellscript=None):
    """
    | Focus testing for SLC data using autofocus estimation of effective velocity
    | Copyright 2023, Gamma Remote Sensing, v1.6 18-Apr-2023 clw/uw
    
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
        fit a1 for first derivative of the effective velocity w.r.t.range (enter - for default)
            * 0: no (default)
            * 1: yes
    
    b0_flg:
        fit b0 for first derivative of the effective velocity w.r.t. along-track time (enter - for default)
            * 0: no (default)
            * 1: yes
    
    offsets:
        (output) range and azimuth offsets and SNR data in text format (enter - for no output)
    n_ovr:
        SLC oversampling factor (1,2,4: enter - for default: 1)
    roff:
        range offset for single patch center (enter - for default: image center in range)
    azoff:
        azimuth offset for single patch center (enter - for default: image center in azimuth)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/af_SLC', SLC_par, SLC, rwin, azwin, dr, daz, thres, a1_flg, b0_flg, offsets, n_ovr, roff, azoff]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def ASAR_LO_phase_drift(SLC1_par, SLC2_par, OFF_par, ph_drift, logpath=None, outdir=None, shellscript=None):
    """
    | Calculate interferometric phase correction due to drift of the ASAR local oscillator
    | Copyright 2023, Gamma Remote Sensing, v1.2 18-Apr-2023 clw
    
    Parameters
    ----------
    SLC1_par:
        (input) SLC1 ISP image parameter file
    SLC2_par:
        (input) SLC2 ISP image parameter file
    OFF_par:
        (input) ISP offset/interferogram parameter file 
    ph_drift:
        (output) interferometric phase correction due to drift of the ASAR LO (radians)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/ASAR_LO_phase_drift', SLC1_par, SLC2_par, OFF_par, ph_drift]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def ASAR_XCA(ASA_XCA, antenna, swath='-', pol='-', logpath=None, outdir=None, shellscript=None):
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
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/ASAR_XCA', ASA_XCA, antenna, swath, pol]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def ave_cpx(cpx_list, width, cpx_ave, start='-', nlines='-', zflag='-', logpath=None, outdir=None, shellscript=None):
    """
    | Calculate average of a set of FCOMPLEX images
    | Copyright 2022, Gamma Remote Sensing, v2.1 17-Aug-2022 clw/cm
    
    Parameters
    ----------
    cpx_list:
        (input) list of coregistered images (FCOMPLEX)
    width:
        number of samples/line
    cpx_ave:
        (output) average of images listed in cpx_list (FCOMPLEX)
    start:
        starting line (enter - for default: 1)
    nlines:
        number of lines to process (enter - for default: entire file)
    zflag:
        zero flag (enter - for default)
            * 0: interpret 0.0 as missing data value (default)
            * 1: interpret 0.0 as valid data
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/ave_cpx', cpx_list, width, cpx_ave, start, nlines, zflag]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def ave_image(im_list, width, ave_image, start='-', nlines='-', pixav_x='-', pixav_y='-', zflag='-', nmin='-', logpath=None, outdir=None, shellscript=None):
    """
    | Calculate average of a set of FLOAT images
    | Copyright 2022, Gamma Remote Sensing, v2.6 17-Aug-2022 clw/cm
    
    Parameters
    ----------
    im_list:
        (input) list of coregistered images (FLOAT)
    width:
        number of samples/line
    ave_image:
        (output) average of images listed in im_list (FLOAT)
    start:
        starting line (enter - for default: 1)
    nlines:
        number of lines to process (enter - for default: entire file)
    pixav_x:
        number of pixels to average in width  (enter - for default: 1)
    pixav_y:
        number of pixels to average in height (enter - for default: 1)
    zflag:
        zero flag (enter - for default)
            * 0: interpret 0.0 as missing data value (default)
            * 1: interpret 0.0 as valid data
    
    nmin:
        minimum number of images required to calculate the average if zflag = 0 (enter - for default: 3/4\\*nfiles)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/ave_image', im_list, width, ave_image, start, nlines, pixav_x, pixav_y, zflag, nmin]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def az_integrate(data, width, azi, cflg, scale='-', lz='-', logpath=None, outdir=None, shellscript=None):
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
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/az_integrate', data, width, azi, cflg, scale, lz]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def az_spec_SLC(SLC, SLC_par, spectrum, roff='-', namb='-', pltflg='-', logpath=None, outdir=None, shellscript=None):
    """
    | Doppler centroid estimate from SLC images
    | Copyright 2023, Gamma Remote Sensing, v3.0 19-Apr-2023 clw
    
    Parameters
    ----------
    SLC:
        (input) SAR image data file (FCOMPLEX or SCOMPLEX format)
    SLC_par:
        (input) ISP SLC image parameter file
    spectrum:
        (output) Doppler spectrum (text format)
    roff:
        range sample offset to center of estimation window (enter - for default: center of swath)
    namb:
        number of multiples of the PRF to add to the estimated centroid (enter - for default: 0)
    pltflg:
        azimuth spectrum plotting flag (enter - for default)
            * 0: none (default)
            * 1: output plot in PNG format
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/az_spec_SLC', SLC, SLC_par, spectrum, roff, namb, pltflg]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def base_copy(SLC1_par, baseline1, SLC2_par, baseline2, time_rev='-', logpath=None, outdir=None, shellscript=None):
    """
    | Calculate baseline file for a subsection of a reference SLC
    | Copyright 2023, Gamma Remote Sensing, v1.2 24-Apr-2023 ts/clw/uw
    
    Parameters
    ----------
    SLC1_par:
        (input) ISP image parameter file of the reference SLC
    baseline1:
        (input) baseline file derived using the reference SLC geometry
    SLC2_par:
        (input) ISP image parameter file corresponding to the subsection of the reference SLC
    baseline2:
        (output) baseline file derived using the geometry and timing of the SLC subsection
    time_rev:
        SLC image time reversal flag (enter - for default)
            * 1: normal (default)
            * -1: time-reversed
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/base_copy', SLC1_par, baseline1, SLC2_par, baseline2, time_rev]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def base_est_fft(interf, SLC1_par, OFF_par, baseline, nazfft='-', r_samp='-', az_line='-', nrfft='-', logpath=None, outdir=None, shellscript=None):
    """
    | Estimate baseline from interferogram phase spectrum
    | Copyright 2023, Gamma Remote Sensing, v2.3 clw/uw 18-Apr-2023
    
    Parameters
    ----------
    interf:
        (input) multilook interferogram with residual range and azimuth fringes
    SLC1_par:
        (input) SLC1 ISP image parameter file
    OFF_par:
        (input) ISP offset/interferogram parameter file
    baseline:
        (output) baseline file
    nazfft:
        size of azimuth FFT (2\\*\\*N) (enter - for default: 512)
    r_samp:
        range pixel offset to center of the FFT window (enter - for default: center)
    az_line:
        line offset from start of the interf. for the  FFT window (enter - for default: center)
    nrfft:
        size of the range FFT (2\\*\\*N), minimum: 32 (enter - for default: 512)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/base_est_fft', interf, SLC1_par, OFF_par, baseline, nazfft, r_samp, az_line, nrfft]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def base_init(SLC1_par, SLC2_par, OFF_par, interf, baseline, mflag='-', nrfft='-', nazfft='-', r_samp='-', az_line='-', logpath=None, outdir=None, shellscript=None):
    """
    | Estimate initial baseline using orbit state vectors, offsets, and interferogram phase
    | Copyright 2023, Gamma Remote Sensing, v2.8 18-Apr-2023 clw/uw/cm
    
    Parameters
    ----------
    SLC1_par:
        (input) SLC1 ISP image parameter file
    SLC2_par:
        (input) SLC2 ISP image parameter file
    OFF_par:
        (input) ISP offset/interferogram parameter file (enter - for none)
    interf:
        (input) unflattened interferogram (enter - for none)
            base      (output) baseline parameter file
    baseline:
        not documented
    mflag:
        baseline estimation method flag (enter - for default)
            mflag    b_para    b_perp    input
            * 0:     orbits    orbits    p1,p2  (default)
            * 1:     offsets   offsets   p1,p2,off
            * 2:     orbits    fft       p1,p2,off,int
            * 3:     offsets   fft       p1,p2,off,int
            * 4:     fft       fft       p1,off,int   
    
    nrfft:
        size of range FFT   (512, 1024,...) (enter - for default determined from image width)
    nazfft:
        size of azimuth FFT (512, 1024,...) (enter - for default determined from image azimuth lines)
    r_samp:
        range pixel offset to center of the FFT window (enter - for default, default: range center)
    az_line:
        line offset from start of the interf. for the  FFT window (enter - for default, default: azimuth center)
            * NOTE: Not all input data files are required for the different methods
              enter - for files that are not provided
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/base_init', SLC1_par, SLC2_par, OFF_par, interf, baseline, mflag, nrfft, nazfft, r_samp, az_line]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def base_ls(SLC_par, OFF_par, gcp_ph, baseline, ph_flag='-', bc_flag='-', bn_flag='-', bcdot_flag='-', bndot_flag='-', bperp_min='-', SLC2R_par='-', logpath=None, outdir=None, shellscript=None):
    """
    | Least squares baseline estimation using terrain heights
    | Copyright 2023, Gamma Remote Sensing, v2.4 18-Apr-2023 clw/uw/cm
    
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
        restore range phase ramp (enter - for default)
            * 0: do not restore (default)
            * 1: restore
    
    bc_flag:
        cross-track baseline component estimate (enter - for default)
            * 0: orbit-derived
            * 1: estimate from data (default)
    
    bn_flag:
        normal baseline component estimate (enter - for default)
            * 0: orbit-derived
            * 1: estimate from data (default)
    
    bcdot_flag:
        cross-track baseline rate estimate (enter - for default)
            * 0: orbit-derived
            * 1: estimate from data (default)
    
    bndot_flag:
        normal baseline rate estimate (enter - for default)
            * 0: orbit-derived (default)
            * 1: estimate from data
    
    bperp_min:
        minimum perpendicular baseline required for L.S estimation (m, enter - for default:  10.0)
    SLC2R_par:
        (input) parameter file of resampled SLC, required if SLC2 frequency differs from SLC1 (enter - for none)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/base_ls', SLC_par, OFF_par, gcp_ph, baseline, ph_flag, bc_flag, bn_flag, bcdot_flag, bndot_flag, bperp_min, SLC2R_par]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def base_orbit(SLC1_par, SLC2_par, baseline, logpath=None, outdir=None, shellscript=None):
    """
    | Estimate baseline from orbit state vectors
    | Copyright 2023, Gamma Remote Sensing, v4.5 clw/cm 18-Apr-2023
    
    Parameters
    ----------
    SLC1_par:
        (input) SLC1 ISP image parameter file
    SLC2_par:
        (input) SLC2 ISP image parameter file
    baseline:
        (output) baseline file (text format, enter - for none)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/base_orbit', SLC1_par, SLC2_par, baseline]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def base_perp(baseline, SLC1_par, OFF_par, time_rev='-', logpath=None, outdir=None, shellscript=None):
    """
    | Calculate baseline components perpendicular and parallel to look vector
    | Copyright 2023, Gamma Remote Sensing, v3.6 18-Apr-2023 clw/uw
    
    Parameters
    ----------
    baseline:
        (input) baseline file (text)
    SLC1_par:
        (input) ISP parameter file of SLC1 (reference SLC)
    OFF_par:
        (input) ISP interferogram/offset parameter file
    time_rev:
        SLC image time reversal flag (enter - fo default)
            * 1: normal (default)
            * -1: time-reversed
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/base_perp', baseline, SLC1_par, OFF_par, time_rev]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def bpf(data_in, data_out, width, fc_x, bw_x, fc_y, bw_y, roff='-', azoff='-', nr='-', naz='-', dtype='-', zflag='-', beta='-', fir_len='-', logpath=None, outdir=None, shellscript=None):
    """
    | Interferometric SAR Processor (ISP): Program GAMMA_SOFTWARE-20250625/ISP/bin/bpf.c
    | Copyright 2023, Gamma Remote Sensing, v1.9 18-Apr-2023 clw
    | Bandpass filter for 2-dimensional image data (FCOMPLEX, SCOMPLEX, and FLOAT)
    
    Parameters
    ----------
    data_in:
        (input) input image data  file
    data_out:
        (output) bandpass filtered image data
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
        offset to starting range to filter (enter - for default: 0)
    azoff:
        offset to starting azimuth to filter (enter - for default: 0)
    nr:
        number of range pixels to filter  (enter - for default: width - roff)
    naz:
        number of azimuth lines to filter (enter - for default: nlines - azoff)
    dtype:
        data type (enter - for default)
            * 0: FCOMPLEX (default)
            * 1: SCOMPLEX
            * 2: FLOAT
    
    zflag:
        zero data flag (enter - for default)
            * 0: set output to 0.0 when the input data are 0.0 (no_data)(default)
            * 1: 0.0 values are considered as valid data
    
    beta:
        Kaiser window beta parameter (enter - for default:    4.538)
    fir_len:
        finite impulse response filter length (enter - for default: 64)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/bpf', data_in, data_out, width, fc_x, bw_x, fc_y, bw_y, roff, azoff, nr, naz, dtype, zflag, beta, fir_len]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def bridge_unw(int, flag, unw, bridge, width, xmin='-', xmax='-', ymin='-', ymax='-', logpath=None, outdir=None, shellscript=None):
    """
    | Phase unwrap new regions with bridges to regions already unwrapped
    | Copyright 2023, Gamma Remote Sensing, v1.5 19-Apr-2023 clw
    
    Parameters
    ----------
    int:
        (input) interferogram (FCOMPLEX)
    flag:
        (input) unwrapping flag file
    unw:
        (input/output) unwrapped phase (FLOAT) 
    bridge:
        (input) bridge data file (text format)
    width:
        number of samples/row
    xmin:
        starting range pixel offset to unwrap (enter - for default: 0)
    xmax:
        last range pixel offset to unwrap (enter - for default: width-1)
    ymin:
        starting azimuth row offset to unwrap, relative to start (enter - for default: 0)
    ymax:
        last azimuth row offset to unwrap, relative to start (enter - for default: nlines-1)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/bridge_unw', int, flag, unw, bridge, width, xmin, xmax, ymin, ymax]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def cc_wave(interf, MLI1, MLI2, cc, width, bx='-', by='-', wflg='-', xmin='-', xmax='-', ymin='-', ymax='-', logpath=None, outdir=None, shellscript=None):
    """
    | Estimate interferometric correlation coefficient
    | Copyright 2023, Gamma Remote Sensing, v6.4 6-Dec-2023 clw/uw/cm
    
    Parameters
    ----------
    interf:
        (input) normalized complex interferogram (FCOMPLEX)
    MLI1:
        (input) multilook intensity image of the first scene (FLOAT) (enter - for none)
    MLI2:
        (input) multilook intensity image of the second scene (FLOAT) (enter - for none)
    cc:
        (output) estimated correlation coefficient (FLOAT)
    width:
        number of samples/line
    bx:
        estimation window size in columns (enter - for default: 5.0)
    by:
        estimation window size in lines (enter - for default: 5.0)
    wflg:
        estimation window (enter - for default):
            * 0: rectangular (constant weighting) (default)
            * 1: circular triangular
            * 2: circular Gaussian
            * 3: normalized vector sum with rectangular window (constant weighting)
            * NOTE: This estimator does not use the MLI data
    
    xmin:
        starting range pixel offset (enter - for default: 0)
    xmax:
        last range pixel offset (enter - for default: width - 1)
    ymin:
        starting azimuth row offset, relative to start (enter -  for default: 0)
    ymax:
        last azimuth row offset, relative to start (enter - for default: nlines - 1)
            * NOTE:   The normalized vector sum (wflg = 3) is used as estimator when the MLI images are not provided.
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/cc_wave', interf, MLI1, MLI2, cc, width, bx, by, wflg, xmin, xmax, ymin, ymax]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def clear_flag(flag, width, flag_bits, xmin='-', xmax='-', ymin='-', ymax='-', logpath=None, outdir=None, shellscript=None):
    """
    | Clear phase unwrapping flag bits
    | Copyright 2023, Gamma Remote Sensing, v1.7 19-Apr-2023 clw
    
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
        starting range pixel offset (enter - for default: 0)
    xmax:
        last range pixel offset (enter - for default: width-1)
    ymin:
        starting azimuth row offset, relative to start (enter - for default: 0)
    ymax:
        last azimuth row offset, relative to start (enter - for default: nlines-1)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/clear_flag', flag, width, flag_bits, xmin, xmax, ymin, ymax]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def corr_flag(corr, flag, width, corr_thr, xmin='-', xmax='-', ymin='-', ymax='-', border='-', logpath=None, outdir=None, shellscript=None):
    """
    | Low correlation region detection for phase unwrapping
    | Copyright 2023, Gamma Remote Sensing, v2.6 19-Apr-2023 clw/uw
    
    Parameters
    ----------
    corr:
        (input)interferometric correlation file
    flag:
        (input/output) phase unwrapping flag filename 
    width:
        number of samples/row
    corr_thr:
        correlation threshold (0 --> 1.0)
    xmin:
        starting range pixel offset (enter - for default: 0)
    xmax:
        last range pixel offset (enter - for default: width-1)
    ymin:
        starting azimuth row offset, relative to start (enter - for default: 0)
    ymax:
        last azimuth row offset, relative to start (enter - for default: nlines-1)
    border:
        effective range of low coherence pixels to set low coherence flag (enter - for default: 2)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/corr_flag', corr, flag, width, corr_thr, xmin, xmax, ymin, ymax, border]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def create_offset(SLC1_par, SLC2_par, OFF_par, algorithm='-', rlks='-', azlks='-', iflg='-', logpath=None, outdir=None, shellscript=None):
    """
    | Create and update ISP offset and interferogram parameter files
    | Copyright 2023 Gamma Remote Sensing v5.6 18-Apr-2023 clw/uw/cm
    
    Parameters
    ----------
    SLC1_par:
        (input) SLC1/MLI1 ISP image parameter filename (reference)
    SLC2_par:
        (input) SLC2/MLI2 ISP image parameter filename
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
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/create_offset', SLC1_par, SLC2_par, OFF_par, algorithm, rlks, azlks, iflg]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def dcomp_sirc(infile, outfile, samples, loff='-', nlines='-', logpath=None, outdir=None, shellscript=None):
    """
    | Extract SIR-C SLC compressed single-pol data
    | Copyright 2023, Gamma Remote Sensing, v1.5 18-Apr-2023 clw
    
    Parameters
    ----------
    infile:
        (input) SIR-C single-pol SLC compressed data
    outfile:
        (output) complex floating point data
    samples:
        number of polarimetric samples per input line (4 bytes/sample)
    loff:
        offset to starting line (enter - for default: 0)
    nlines:
        number of lines to copy (enter - or 0 for default: entire file)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/dcomp_sirc', infile, outfile, samples, loff, nlines]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def dcomp_sirc_quad(infile, outfile, samples, parameter, loff='-', nlines='-', logpath=None, outdir=None, shellscript=None):
    """
    | Extract SIR-C MLC or SLC compressed quad-pol data
    | Copyright 2023, Gamma Remote Sensing, v1.5 18-Apr-2023 uw/clw
    
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
            * 6:  MLC-HVHV\\*
            * 7:  MLC-VVVV\\*
            * 8:  MLC-HHHH\\*
            * 9:  MLC-HHHV\\*
            * 10: MLC-HHVV\\*
            * 11: MLC-HVVV\\*
    
    loff:
        offset to starting line (enter - for default: 0)
    nlines:
        number of lines to copy (enter - or 0 for default: entire file)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/dcomp_sirc_quad', infile, outfile, samples, parameter, loff, nlines]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def DELFT_vec2(SLC_par, DELFT_dir, nstate='-', interval='-', ODR='-', logpath=None, outdir=None, shellscript=None):
    """
    | Extract and interpolate DELFT ERS-1, ERS-2, and ENVISAT state vectors
    | Copyright 2023, Gamma Remote Sensing, v2.7 19-Apr-2023 clw
    
    Parameters
    ----------
    SLC_par:
        (input) ISP image parameter file
    DELFT_dir:
        directory containing Delft orbit arclist and ODR files for ERS-1, ERS-2 or ENVISAT
            * NOTE: enter . for current directory
    
    nstate:
        number of state vectors to generate (enter - for default, >= 15)
    interval:
        time interval between state vectors in the ISP image parameter file (s) (enter - for default: 10.0)
    ODR:
        ODR file to use (include path) rather than ODR file determined from the Delft orbit arclist (enter - for none)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/DELFT_vec2', SLC_par, DELFT_dir, nstate, interval, ODR]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def doppler_2d_SLC(SLC, SLC_par, dop2d, loff='-', blsz='-', nbl='-', a2_flg='-', b0_flg='-', b1_flg='-', c0_flg='-', namb='-', logpath=None, outdir=None, shellscript=None):
    """
    | 2-D Doppler centroid trend estimation from SLC data
    | Copyright 2025, Gamma Remote Sensing, v1.3 14-May-2025 clw/cm
    
    Parameters
    ----------
    SLC:
        (input) SLC image (SCOMPLEX or FCOMPLEX format)
    SLC_par:
        (input) SLC parameter file
    dop2d:
        (output) estimated doppler centroid as a function of range for each block (text format) (enter - for none)
    loff:
        number of lines offset (enter - for default: 0)
    blsz:
        block size lines, minimum: 256 (enter - for default: 2048)
    nbl:
        number of blocks (enter - for default: calculated automatically)
    a2_flg:
        fit a2 for second derivative of the Doppler centroid w.r.t.range (Hz/m/m) (enter - for default)
            * 0: no (default)
            * 1: yes
    
    b0_flg:
        fit b0 for first derivative of the Doppler centroid w.r.t. along-track time (Hz/sec) (enter - for default)
            * 0: no
            * 1: yes (default)
    
    b1_flg:
        fit b1 for along-track rate of the change in slope of Doppler w.r.t. range (Hz/sec/m)(enter - for default)
            * 0: no
            * 1: yes (default)
    
    c0_flg:
        fit c0 for second derivative of the Doppler centroid w.r.t. along-track time (Hz/sec/sec) (enter - for default)
            * 0: no (default)
            * 1: yes
    
    namb:
        user defined number of Doppler ambiguities to add to the Doppler function (enter - for default: 0)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/doppler_2d_SLC', SLC, SLC_par, dop2d, loff, blsz, nbl, a2_flg, b0_flg, b1_flg, c0_flg, namb]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def DORIS_vec(SLC_par, DOR, nstate='-', logpath=None, outdir=None, shellscript=None):
    """
    | Extract ENVISAT DORIS state vectors and write to an ISP image parameter file
    | Copyright 2023, Gamma Remote Sensing, v1.5 18-Apr-2023 clw
    
    Parameters
    ----------
    SLC_par:
        (input/output)ISP SLC/MLI image parameter file
    DOR:
        (input) ASAR DORIS data file (DOR_VOR_AXVF)
    nstate:
        number of state vectors to extract (enter - for default: 11)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/DORIS_vec', SLC_par, DOR, nstate]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def error_stat(d1, d2, width, dtype, roff, loff, nr, nl, report, logpath=None, outdir=None, shellscript=None):
    """
    | Calculate statistics for two data files and their difference (FLOAT or FCOMPLEX)
    | Copyright 2017, Gamma Remote Sensing, v1.2 clw 7-Jan-2016
    
    Parameters
    ----------
    d1:
        (input) data file 1
    d2:
        (input) data file 2
    width:
        image line width (samples/line)
    dtype:
        data type for d1 and d2:
            * 0: FLOAT
            * 1: FCOMPLEX
    
    roff:
        sample offset to region start (enter - for default: 0)
    loff:
        line offset to region start (enter - for default: 0)
    nr:
        region width (samples, enter - for default: width - roff)
    nl:
        number of lines in the region (enter - for default: data_lines - loff)
    report:
        output text file (keyword:value format)
            keywords: data_1, data_2, d1_mean, d2_mean, d1_stddev, d2_stddev, root_mean_square_error, normalized_mean_square_error,
            cross_correlation_coefficient, cross_correlation_angle, total_samples, non_zero_samples, valid_fraction
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/error_stat', d1, d2, width, dtype, roff, loff, nr, nl, report]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def fill_gaps(data_in, width, data_out, dtype='-', method='-', max_dist='-', bp_flag='-', win='-', ds_method='-', ds_size='-', ds_data='-', logpath=None, outdir=None, shellscript=None):
    """
    | Fill gaps in 2D raster file
    | Copyright 2023, Gamma Remote Sensing, v2.4 18-Apr-2023 cm
    
    Parameters
    ----------
    data_in:
        (input) input data file (FLOAT / FCOMPLEX)
    width:
        width of input data
    data_out:
        (output) output data file (FLOAT / FCOMPLEX)
    dtype:
        input and output data type (enter - for default)
            * 0: FLOAT (default)
            * 1: FCOMPLEX
    
    method:
        method flag (enter - for default: 4)
            * 0: Laplace interpolation and linear extrapolation - least squares solution
            * 1: Laplace interpolation and linear extrapolation - smaller system of linear equations than in method #0 in case of few missing values - least squares solution
            * 2: Laplace interpolation and linear extrapolation - solves a direct linear system of equations for the missing values (not a least squares solution)
            * 3: biharmonic interpolation - implementation similar to method #1 - least squares solution
            * 4: spring analogy: assumes springs (with a nominal length of zero) connect each node with every neighbor - least squares solution (default)
            * 5: average of the 8 nearest neighbors - this method solves a direct linear system for the missing values (not a least squares solution)
            * NOTE: small gaps: use method #0, #1 or #3 - large gaps: use method #2, #4 or #5 - most demanding: method #3
    
    max_dist:
        maximum interpolation / extrapolation distance in pixels (enter - or 0 for default: unlimited)
    bp_flag:
        perform block processing (enter - for default: 0)
            * 0: no block processing (default)
            * 1: block processing (faster, avoid overflow, however might be slightly less accurate)
            * NOTE: when block processing is selected, a two-step process is carried out: 1: solving the downsampled array (coarse processing), 2: block processing
    
    win:
        block size (pixels, 10 < win < 1000, enter - for default: 100)
    ds_method:
        method flag (0 - 5, same choices as for [method] option) (enter - for default: same as [method])
            * NOTE: for an input containing large gaps, method #2, #4 or #5 may yield more appropriate results.
    
    ds_size:
        maximum size of downsampled data (for both width and height) (pixels, ds_size > 10, enter - for default: 400)
    ds_data:
        (output) write intermediate data after solving the downsampled array (FLOAT / FCOMPLEX)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/fill_gaps', data_in, width, data_out, dtype, method, max_dist, bp_flag, win, ds_method, ds_size, ds_data]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def fspf(data_in, data_out, width, dtype='-', r_max='-', spf_type='-', MLI_par='-', interp_mode='-', order='-', logpath=None, outdir=None, shellscript=None):
    """
    | ISP fspf: Fast spatial filter for 2D data
    | Copyright 2025, Gamma Remote Sensing, v2.0 9-Apr-2025 of/clw/uw/cm
    
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
            * 0: FCOMPLEX
            * 1: SCOMPLEX
            * 2: FLOAT (default)
    
    r_max:
        maximum filter radius (range samples) (enter - for default: 64)
    spf_type:
        spatial filter type (enter - for default):
            * 0: uniform average (default for FCOMPLEX and SCOMPLEX)
            * 1: triangular weighted average: 1 - (r/r_max)
            * 2: quadratic weighted average: 1 - (r/r_max)^2
            * 3: Gaussian weighted average: exp(-2.\\*(r^2/r_max^2))
            * 4: linear least-squares (default for FLOAT data)
            * 5: median
    
    MLI_par:
        MLI or SLC parameter file with the same number of looks as the input image, required for GPRI data (enter - for none)
    interp_mode:
        interpolation method for resampling the data to the original size after filtering
            * 0: bicubic spline (default)
            * 1: bicubic spline sqrt(x)
            * 2: B-spline interpolation (default B-spline degree: 3)
            * 3: B-spline interpolation sqrt(x) (default B-spline degree: 3)
    
    order:
        B-Spline interpolation degree (2->9) (enter - default: 3)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/fspf', data_in, data_out, width, dtype, r_max, spf_type, MLI_par, interp_mode, order]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def gcp_phase(unw, OFF_par, gcp, gcp_ph, win_sz='-', logpath=None, outdir=None, shellscript=None):
    """
    | Extract unwrapped phase at GCP locations
    | Copyright 2023, Gamma Remote Sensing, v1.6 19-Apr-2023 clw
    
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
        window size for averaging phase for each GCP, must be odd (enter - for default: 1)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/gcp_phase', unw, OFF_par, gcp, gcp_ph, win_sz]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def grasses(int, flag, unw, width, xmin='-', xmax='-', ymin='-', ymax='-', xinit='-', yinit='-', init_ph='-', logpath=None, outdir=None, shellscript=None):
    """
    | Phase unwrapping by region growing
    | Copyright 2023, Gamma Remote Sensing, v4.4 19-Apr-2023 clw/uw
    
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
        starting range pixel offset (enter - for default: 0)
    xmax:
        last range pixel offset (enter - for default: width-1)
    ymin:
        starting azimuth row offset, relative to start (enter - for default: 0)
    ymax:
        last azimuth row offset, relative to start (enter - for default: nlines-1)
    xinit:
        starting range pixel for unwrapping (enter - for default: width/2)
    yinit:
        starting row to unwrap (enter - for default: height/2)
    init_ph:
        flag to set phase at starting point to 0.0 (enter - for default)
            * 0: not set to 0.0 (default)
            * 1: set to 0.0
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/grasses', int, flag, unw, width, xmin, xmax, ymin, ymax, xinit, yinit, init_ph]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def GRD_to_SR(GRD_par, MLI_par, OFF_par, in_file, out_file, rlks='-', azlks='-', interp_mode='-', sr_rsp='-', sr_azsp='-', degree='-', logpath=None, outdir=None, shellscript=None):
    """
    | Conversion to slant range for ISP MLI and INSAR ground range data of type FLOAT
    | Copyright 2023, Gamma Remote Sensing, v2.5 18-Apr-2023 uw/clw/cm
    
    Parameters
    ----------
    GRD_par:
        (input) SLC parameter file of output ground range image
    MLI_par:
        (input/output) MLI ISP image parameter file for slant range image
            * NOTE: delete an existing MLI parameter file to recalculate the output MLI parameters
    
    OFF_par:
        (input) ISP offset/interferogram parameter file of input image (enter - image in MLI geometry)
    in_file:
        (input) ground range image (FLOAT)
    out_file:
        (output) slant range image (FLOAT)
    rlks:
        multi-looking in range (prior to resampling, enter - for default: 1)
    azlks:
        multi-looking in azimuth (prior to resampling, enter - for default: 1)
    interp_mode:
        interpolation mode (enter - for default)
            * 0: nearest-neighbor
            * 1: bicubic spline
            * 2: bicubic spline log(x)
            * 3: bicubic spline sqrt(x)
            * 4: B-spline interpolation (default B-spline degree: 3)
            * 5: B-spline interpolation sqrt(x) (default) (default B-spline degree: 3)
            * NOTE: log and sqrt interpolation modes should only be used with non-negative data!
    
    sr_rsp:
        output image slant range sample spacing (m) (enter - for default: c/(2\\*adc_sampling_rate)
    sr_azsp:
        output image azimuth sample spacing (m) (enter - for default: (input image azimuth spacing) \\* azlks)
    degree:
        B-spline degree (2->9) (enter - for default: 3)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/GRD_to_SR', GRD_par, MLI_par, OFF_par, in_file, out_file, rlks, azlks, interp_mode, sr_rsp, sr_azsp, degree]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def hgt_map(unw, SLC_par, OFF_par, baseline, hgt, gr, ph_flag='-', loff='-', nlines='-', SLC2R_par='-', logpath=None, outdir=None, shellscript=None):
    """
    | Interferometric height/ground range estimation vs. slant range
    | Copyright 2023, Gamma Remote Sensing, v5.3 clw/uw 18-Apr-2023
    
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
        restore phase slope flag (enter - for default)
            * 0: no phase change
            * 1: add back phase ramp (default)
    
    loff:
        offset to starting line (enter - for default: 0)
    nlines:
        number of lines to calculate (enter - for default: to end of file)
    SLC2R_par:
        (input) parameter file of resampled SLC, required if SLC2 frequency differs from SLC1 (enter - for none)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/hgt_map', unw, SLC_par, OFF_par, baseline, hgt, gr, ph_flag, loff, nlines, SLC2R_par]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def image_stat(image, width, roff='-', loff='-', nr='-', nl='-', report='-', median_flg='-', logpath=None, outdir=None, shellscript=None):
    """
    | Calculate mean, standard deviation, number of non-zero values, min, max and median for a rectangular image region (FLOAT format)
    | Copyright 2025, Gamma Remote Sensing, v1.6 27-May-2025 clw/cm
    
    Parameters
    ----------
    image:
        (input) image data file (FLOAT)
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
        output text file (keyword:value format, enter - for none)
            keywords: file, mean, stdev, total_samples, non_zero_samples, fraction_valid, min, max, median
    median_flg:
        median calculation flag (enter - for default)
            * 0: do not calculate median
            * 1: calculate median (default, memory use may be large)
            * NOTE: only the non-zero samples are considered in the statistical values
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/image_stat', image, width, roff, loff, nr, nl, report, median_flg]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def init_offset(SLC1, SLC2, SLC1_par, SLC2_par, OFF_par, rlks='-', azlks='-', rpos='-', azpos='-', offr='-', offaz='-', thres='-', rwin='-', azwin='-', cflag='-', deramp='-', logpath=None, outdir=None, shellscript=None):
    """
    | Determine initial offset between SLC images using correlation of image intensity
    | Copyright 2023, Gamma Remote Sensing, v3.3 clw/cm 18-Apr-2023
    
    Parameters
    ----------
    SLC1:
        (input) single-look complex image 1 (reference)
    SLC2:
        (input) single-look complex image 2
    SLC1_par:
        (input) SLC1 ISP image parameter file
    SLC2_par:
        (input) SLC2 ISP image parameter file
    OFF_par:
        (input) ISP offset/interferogram parameter file
    rlks:
        number of range looks (enter - for default: 1)
    azlks:
        number of azimuth looks (enter - for default: 1)
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
        range window size (enter - for default: 512)
    azwin:
        azimuth window size (enter - for default: 512)
    cflag:
        copy offsets to the range and azimuth offset polynomials in the OFF_par (enter - for default)
            * 0: do not copy
            * 1: copy constant range and azimuth offset (default)
    
    deramp:
        deramp SLC phase flag (enter - for default)
            * 0: no deramp (Doppler centroid close to 0) (default)
            * 1: deramp SLC phase
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/init_offset', SLC1, SLC2, SLC1_par, SLC2_par, OFF_par, rlks, azlks, rpos, azpos, offr, offaz, thres, rwin, azwin, cflag, deramp]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def init_offset_orbit(SLC1_par, SLC2_par, OFF_par, rpos='-', azpos='-', cflag='-', logpath=None, outdir=None, shellscript=None):
    """
    | Initial SLC image offset estimation from orbit state-vectors and image parameters
    | Copyright 2020, Gamma Remote Sensing, v1.9 18-Apr-2023 clw/uw/cm
    
    Parameters
    ----------
    SLC1_par:
        (input) SLC1 parameter file
    SLC2_par:
        (input) SLC2 parameter file
    OFF_par:
        (input/output) ISP/offset parameter file
    rpos:
        range position for offset estimation (enter - for default: center of SLC1)
    azpos:
        azimuth position for offset estimation (enter - for default: center of SLC1)
    cflag:
        copy offsets to the range and azimuth offset polynomials in the OFF_par (enter - for default)
            * 0: do not copy
            * 1: copy constant range and azimuth offset (default)
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/init_offset_orbit', SLC1_par, SLC2_par, OFF_par, rpos, azpos, cflag]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def interf_SLC(SLC1, SLC2, SLC1_par, SLC2_par, OFF_par, MLI1, MLI2, interf, rlks='-', azlks='-', loff='-', nltot='-', rfilt='-', azfilt='-', s_off='-', logpath=None, outdir=None, shellscript=None):
    """
    | Interferogram generation using a pair of SLC images
    | Copyright 2023, Gamma Remote Sensing, v5.0 clw/uw 18-Apr-2023
    
    Parameters
    ----------
    SLC1:
        (input) single-look complex image 1 (reference)
    SLC2:
        (input) single-look complex image 2
    SLC1_par:
        (input) SLC1 ISP image parameter file
    SLC2_par:
        (input) SLC2 ISP image parameter file
    OFF_par:
        (input) ISP offset/interferogram parameter file
    MLI1:
        (output) multi-look intensity image 1
    MLI2:
        (output) multi-look intensity image 2
    interf:
        interferogram from SLC1 and SLC2
    rlks:
        number of interferogram range looks (enter - for default: 2)
    azlks:
        number of interferogram azimuth looks (enter - for default: 10)
    loff:
        offset to starting line of interferogram (relative to start of SLC1) (enter - for default: 0)
    nltot:
        number of SLC lines to process (enter - or 0 for default: to end of file)
    rfilt:
        range common band filtering flag (enter - for default)
            * 0: OFF
            * 1: ON (default)
    
    azfilt:
        azimuth common band filtering flag (enter - for default)
            * 0: OFF
            * 1: ON (default)
    
    s_off:
        offset to the nominal range spectral shift (frac. of range sampling freq.) (enter - for default: 0.0)
            * NOTE: enter - as filename to avoid creation of corresponding output file
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/interf_SLC', SLC1, SLC2, SLC1_par, SLC2_par, OFF_par, MLI1, MLI2, interf, rlks, azlks, loff, nltot, rfilt, azfilt, s_off]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def interp_ad(data_in, data_out, width, r_max='-', np_min='-', np_max='-', w_mode='-', dtype='-', cp_data='-', logpath=None, outdir=None, shellscript=None):
    """
    | Weighted interpolation of gaps in 2D data using an adaptive smoothing window
    | Copyright 2018, Gamma Remote Sensing, v2.1 13-Jun-2018 clw/uw
    
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
            * 2: 1 - (r/r_max)\\*\\*2  (default)
            * 3: exp(-2.\\*(r\\*\\*2/r_max\\*\\*2))
    
    dtype:
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
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/interp_ad', data_in, data_out, width, r_max, np_min, np_max, w_mode, dtype, cp_data]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def mask_data(data_in, width, data_out, mask, dtype='-', logpath=None, outdir=None, shellscript=None):
    """
    | Mask float or fcomplex data using an 8-bit SUN/BMP/TIFF format raster image
    | Copyright 2022, Gamma Remote Sensing, v1.6 8-Nov-2022 clw/cm
    
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
    
    dtype:
        data format:
            * 0: FLOAT (default)
            * 1: FCOMPLEX
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/mask_data', data_in, width, data_out, mask, dtype]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def mcf(interf, wgt, mask, unw, width, tri_mode='-', roff='-', loff='-', nr='-', nlines='-', npat_r='-', npat_az='-', ovrlap='-', r_init='-', az_init='-', init_flag='-', logpath=None, outdir=None, shellscript=None):
    """
    | Phase unwrapping using Minimum Cost Flow (MCF) on a triangular mesh
    | Copyright 2024, Gamma Remote Sensing, v2.9 clw/uw/cm 4-Apr-2024
    
    Parameters
    ----------
    interf:
        (input) interferogram (\\*.int,\\*.diff,\\*.flt) (FCOMPLEX)
    wgt:
        (input) weight factors (0 -> 1.0, e.g. coherence map) file (FLOAT) (enter - for uniform weights)
    mask:
        (input) validity mask (SUN/BMP/TIFF raster format, value 0 -> pixel not used) (enter - if no mask)
    unw:
        (output) unwrapped phase image (\\*.unw) (FLOAT)
    width:
        number of samples/row
    tri_mode:
        triangulation mode (enter - for default)
            * 0: filled triangular mesh
            * 1: Delaunay triangulation
            * 2: filled triangular mesh, replacing gaps with noise (default)
            * 3: filled triangular mesh, replacing gaps and outside boundary with noise
    
    roff:
        offset to starting range of section to unwrap (enter - for default: 0)
    loff:
        offset to starting line of section to unwrap (enter - for default: 0)
    nr:
        number of range samples of section to unwrap (enter - for default: width - roff)
    nlines:
        number of lines of section to unwrap (enter - for default: total number of lines - loff)
    npat_r:
        number of patches in range (enter - for default: 1, enter 0 to automatically define number of patches)
    npat_az:
        number of patches in azimuth (enter - for default: 1, enter 0 to automatically define number of patches)
    ovrlap:
        overlap between patches in pixels (overlap >= 7, enter - for default: 1024)
    r_init:
        phase reference point range offset (enter - for default: center of valid data bounding box)
    az_init:
        phase reference point azimuth offset (enter - for default: center of valid data bounding box)
    init_flag:
        flag to set phase at reference point (enter - for default)
            * 0: use initial point phase value (default)
            * 1: set phase to 0.0 at initial point
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/mcf', interf, wgt, mask, unw, width, tri_mode, roff, loff, nr, nlines, npat_r, npat_az, ovrlap, r_init, az_init, init_flag]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def MLI_cat(MLI1, MLI2, MLI1_par, MLI2_par, MLI3, MLI3_par, dtype='-', mflg='-', overlap='-', interp_mode='-', degree='-', extrapol='-', logpath=None, outdir=None, shellscript=None):
    """
    | Concatenate two MLI images using B-spline interpolation
    | Copyright 2023, Gamma Remote Sensing, v2.0 18-Apr-2023 awi/cm/clw
    
    Parameters
    ----------
    MLI1:
        (input) MLI1 image (single-look)
    MLI2:
        (input) MLI2 image to be appended to MLI1
    MLI1_par:
        (input) MLI1 ISP image parameter file
    MLI2_par:
        (input) MLI2 ISP image parameter file
    MLI3:
        (output) concatenated MLI image
    MLI3_par:
        (output) ISP image parameter file for concatenated image
    dtype:
        input/output data type (enter - for default)
            * 0: FLOAT (default)
            * 1: FCOMPLEX
            * NOTE: FCOMPLEX is for differential interferograms
    
    mflg:
        mosaicking option flag (enter - for default)
            * 0: in overlapping areas, use MLI2 data to fill MLI1 empty areas (default)
            * 1: in overlapping areas, do not use MLI2 data to fill MLI1 empty areas
    
    overlap:
        number of pixels at the edge of MLI1 valid areas to replace by MLI2 data (only if mflg=0, enter - for default: 0)
    interp_mode:
        interpolation mode in case of different geometries (enter - for default)
            * 0: B-spline interpolation (default for FCOMPLEX)
            * 1: B-spline interpolation sqrt(x) (default for FLOAT)
            * NOTE: sqrt interpolation mode should only be used with non-negative data!
    
    degree:
        B-spline degree (2->9) (enter - default: 4)
    extrapol:
        extrapolation flag (enter - for default)
            * 0: do not extrapolate (default)
            * 1: extrapolate last line if needed
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/MLI_cat', MLI1, MLI2, MLI1_par, MLI2_par, MLI3, MLI3_par, dtype, mflg, overlap, interp_mode, degree, extrapol]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def MLI_copy(MLI_in, MLI_in_par, MLI_out, MLI_out_par, roff='-', nr='-', loff='-', nl='-', logpath=None, outdir=None, shellscript=None):
    """
    | Copy MLI data file with options for segment extraction
    | Copyright 2019, Gamma Remote Sensing, v4.9 15-Oct-2019 uw/clw/cm
    
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
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/MLI_copy', MLI_in, MLI_in_par, MLI_out, MLI_out_par, roff, nr, loff, nl]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def mosaic_WB(data_tab, dtype, data_out, data_par_out, sc_flg='-', logpath=None, outdir=None, shellscript=None):
    """
    | ISP: Program GAMMA_SOFTWARE-20250625/ISP/bin/mosaic_WB.c
    | Copyright 2018, Gamma Remote Sensing, v1.3 26-Apr-2018 clw/cm
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
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/mosaic_WB', data_tab, dtype, data_out, data_par_out, sc_flg]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def multi_look(SLC, SLC_par, MLI, MLI_par, rlks, azlks, loff='-', nlines='-', scale='-', exp='-', logpath=None, outdir=None, shellscript=None):
    """
    | Calculate a multi-look intensity (MLI) image from an SLC image
    | Copyright 2022, Gamma Remote Sensing, v4.7 8-Aug-2022 clw/uw/cm
    
    Parameters
    ----------
    SLC:
        (input) single-look complex image (SCOMPLEX or FCOMPLEX)
    SLC_par:
        (input) SLC ISP image parameter file
    MLI:
        (output) multi-look intensity image (FLOAT)
    MLI_par:
        (output) MLI ISP image parameter file
    rlks:
        number of range looks (INT)
    azlks:
        number of azimuth looks (INT)
    loff:
        offset to starting line (enter - for default: 0)
    nlines:
        number of SLC lines to process (enter - for default: entire file)
    scale:
        scale factor for output MLI (enter - for default: 1.0)
    exp:
        exponent for the output MLI (enter - for default: 1.0)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/multi_look', SLC, SLC_par, MLI, MLI_par, rlks, azlks, loff, nlines, scale, exp]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def multi_look2(SLC, SLC_par, MLI, MLI_par, r_dec, az_dec, rwin='-', azwin='-', wflg='-', n_ovr='-', lanczos='-', beta='-', scale='-', exp='-', logpath=None, outdir=None, shellscript=None):
    """
    | Calculate an MLI image from an SLC with optional oversampling and separate multilooking and decimation factors
    | Copyright 2024, Gamma Remote Sensing, v1.9 10-Jun-2024 clw/cm
    
    Parameters
    ----------
    SLC:
        (input) single-look complex image (SCOMPLEX or FCOMPLEX)
    SLC_par:
        (input) SLC image parameter file
    MLI:
        (output) multi-look intensity image (FLOAT)
    MLI_par:
        (output) MLI image parameter file
    r_dec:
        range decimation factor (int)
    az_dec:
        azimuth decimation factor (int)
    rwin:
        averaging window width (int)  (enter - for default: r_dec)
    azwin:
        averaging window height (int) (enter - for default: az_dec)
    wflg:
        window weighting function (enter - for default):
            * 0: rectangular (default)
            * 1: Kaiser
            * 2: circular Gaussian
    
    n_ovr:
        oversampling factor 1 -> 2 (enter - for default: 1)
    lanczos:
        Lanczos interpolator order 5 -> 9 (enter - for default: 7)
    beta:
        Gaussian or Kaiser window parameter (enter - for default: 2.0)
    scale:
        scale factor for output MLI (enter - for default: 1.0)
    exp:
        exponent for the output MLI (enter - for default: 1.0)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/multi_look2', SLC, SLC_par, MLI, MLI_par, r_dec, az_dec, rwin, azwin, wflg, n_ovr, lanczos, beta, scale, exp]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def multi_look_MLI(MLI_in, MLI_in_par, MLI_out, MLI_out_par, rlks, azlks, loff='-', nlines='-', scale='-', e_flag='-', logpath=None, outdir=None, shellscript=None):
    """
    | Multilooking (averaging and decimation) of MLI images
    | Copyright 2019, Gamma Remote Sensing, v1.9 29-Oct-2019 clw/cm
    
    Parameters
    ----------
    MLI_in:
        (input) multi-look intensity image (MLI) file (float)
    MLI_in_par:
        (input) MLI parameter file
    MLI_out:
        (output) multi-looked MLI image (float)
    MLI_out_par:
        (output) MLI parameter file for output MLI
    rlks:
        range looks for multi-looking
    azlks:
        azimuth looks for multi-looking
    loff:
        offset to starting line (enter - for default: 0)
    nlines:
        number of input MLI lines to process (enter - for default: entire file)
    scale:
        scale factor for output MLI (enter - for default: 1.0)
    e_flag:
        extent flag (enter - for default)
            * 0: only permit pixels with the full number of looks (default)
            * 1: permit pixels without the full number of looks
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/multi_look_MLI', MLI_in, MLI_in_par, MLI_out, MLI_out_par, rlks, azlks, loff, nlines, scale, e_flag]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def multi_look_ScanSAR(SLC_tab, MLI, MLI_par, rlks, azlks, bflg='-', SLCR_tab='-', scale='-', logpath=None, outdir=None, shellscript=None):
    """
    | Calculate MLI mosaic from ScanSAR SLC burst data (Sentinel-1, TerraSAR-X, RCM...)
    | Copyright 2023, Gamma Remote Sensing v4.6 30-Nov-2023 awi/clw/uw/cm
    
    Parameters
    ----------
    SLC_tab:
        (input) 3 column list of ScanSAR SLC, swaths are listed in order from near to far range
            SLC_tab line entries:   SLC   SLC_par  TOPS_par
    MLI:
        (output) mosaicked MLI image (non-overlapping burst windows)
    MLI_par:
        (output) MLI image parameter file
    rlks:
        number of range looks
    azlks:
        number of azimuth looks
    bflg:
        burst window calculation flag (enter - for default):
            * 0: use existing burst window parameters if they exist, otherwise calculate burst window parameters (default)
            * 1: calculate burst window parameters from burst parameters and the number of range and azimuth looks
    
    SLCR_tab:
        (input) 3 column list of the reference scene, swaths are listed in order from near to far range, (enter - for default: none)
            SLCR_tab line entries: SLC  SLC_par  TOPS_par
            When generating an MLI mosaic from resampled ScanSAR SLC data, the SLC_tab of the reference scene must be provided
    scale:
        scale factor for output MLI (enter - for default: calculate from calibration gain in SLC parameter file)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/multi_look_ScanSAR', SLC_tab, MLI, MLI_par, rlks, azlks, bflg, SLCR_tab, scale]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def multi_real(data_in, OFF_par_in, data_out, OFF_par_out, rlks='-', azlks='-', loff='-', nlines='-', roff='-', nsamp='-', logpath=None, outdir=None, shellscript=None):
    """
    | Calculate multi-look averaged or interpolated 2D image (float data)
    | Copyright 2023, Gamma Remote Sensing, v2.7 19-Apr-2023 clw/uw/cm
    
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
        number of range looks, values < -1, interpreted as an image oversampling factor (enter - for default: 1)
    azlks:
        number azimuth looks,  values < -1, interpreted as an image oversampling factor (enter - for default: 1)
    loff:
        line offset to starting line (enter - for default: 0)
    nlines:
        number of lines (enter - for default: to end of file)
    roff:
        offset to starting range sample (enter - for default: 0)
    nsamp:
        number of range samples to extract (enter - for default: to end of line)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/multi_real', data_in, OFF_par_in, data_out, OFF_par_out, rlks, azlks, loff, nlines, roff, nsamp]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def multi_SLC_WSS(SLC, SLC_par, MLI, MLI_par, logpath=None, outdir=None, shellscript=None):
    """
    | Calculate multi-look intensity image (MLI) from a ASAR Wide-Swath SLC
    | Copyright 2023, Gamma Remote Sensing v1.3 18-Apr-2023 clw/awi
    
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
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/multi_SLC_WSS', SLC, SLC_par, MLI, MLI_par]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def neutron(intensity, flag, width, n_thres, ymin='-', ymax='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate phase unwrapping neutrons using image intensity
    | Copyright 2023, Gamma Remote Sensing, v2.4 19-Apr-2023 clw/uw
    
    Parameters
    ----------
    intensity:
        (input) image intensity 
    flag:
        (input) phase unwrapping flag file
    width:
        number of samples/row
    n_thres:
        neutron threshold, multiples of the average intensity (enter - for default: 6.0)
    ymin:
        offset to starting azimuth row (enter - for default: 0)
    ymax:
        offset to last azimuth row (enter - for default: nlines-1)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/neutron', intensity, flag, width, n_thres, ymin, ymax]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def offset_add(OFF_par1, OFF_par2, OFF_par3, logpath=None, outdir=None, shellscript=None):
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
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/offset_add', OFF_par1, OFF_par2, OFF_par3]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def offset_fit(offs, ccp, OFF_par, coffs='-', coffsets='-', thres='-', npoly='-', interact_mode='-', logpath=None, outdir=None, shellscript=None):
    """
    | Range and azimuth offset polynomial estimation
    | Copyright 2023, Gamma Remote Sensing, v3.9 18-Apr-2023 clw/uw/cm
    
    Parameters
    ----------
    offs:
        (input) range and azimuth offset estimates for each patch (FCOMPLEX)
    ccp:
        (input) cross-correlation or SNR of each patch (FLOAT)
    OFF_par:
        (input) ISP offset/interferogram parameter file
    coffs:
        (output) culled range and azimuth offset estimates (FCOMPLEX, enter - for none)
    coffsets:
        (output) culled offset estimates and cross-correlation values (text format, enter - for none)
    thres:
        cross-correlation threshold (enter - for default from OFF_par)
    npoly:
        number of model polynomial parameters (enter - for default, 1, 3, 4, 6, default: 4)
    interact_mode:
        interactive culling of input data (enter - for default)
            * 0: off (default)
            * 1: on
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/offset_fit', offs, ccp, OFF_par, coffs, coffsets, thres, npoly, interact_mode]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def offset_pwr(SLC1, SLC2, SLC1_par, SLC2_par, OFF_par, offs, ccp, rwin='-', azwin='-', offsets='-', n_ovr='-', nr='-', naz='-', thres='-', lanczos='-', bw_frac='-', deramp='-', int_filt='-', pflag='-', pltflg='-', ccs='-', logpath=None, outdir=None, shellscript=None):
    """
    | Offset estimation between SLC images using intensity cross-correlation
    | Copyright 2023, Gamma Remote Sensing, v5.8 clw/cm 18-Apr-2023
    
    Parameters
    ----------
    SLC1:
        (input) single-look complex image 1 (reference)
    SLC2:
        (input) single-look complex image 2
    SLC1_par:
        (input) SLC1 ISP image parameter file
    SLC2_par:
        (input) SLC2 ISP image parameter file
    OFF_par:
        (input) ISP offset/interferogram parameter file
    offs:
        (output) offset estimates in range and azimuth (FCOMPLEX)
    ccp:
        (output) cross-correlation of each patch (0.0->1.0) (FLOAT)
    rwin:
        range patch size (range pixels, enter - for default from offset parameter file)
    azwin:
        azimuth patch size (azimuth lines, enter - for default from offset parameter file)
    offsets:
        (output) range and azimuth offsets and cross-correlation data in text format, enter - for no output
    n_ovr:
        SLC oversampling factor (integer 2\\*\\*N (1,2,4), enter - for default: 2)
    nr:
        number of offset estimates in range direction (enter - for default from offset parameter file)
    naz:
        number of offset estimates in azimuth direction (enter - for default from offset parameter file)
    thres:
        cross-correlation threshold (0.0->1.0) (enter - for default from offset parameter file)
    lanczos:
        Lanczos interpolator order 5 -> 9 (enter - for default: 5)
    bw_frac:
        bandwidth fraction of low-pass filter on complex data (0.0->1.0) (enter - for default: 1.0)
    deramp:
        deramp SLC phase flag (enter - for default)
            * 0: no deramp (Doppler centroid close to 0) (default)
            * 1: deramp SLC phase
    
    int_filt:
        intensity low-pass filter flag (enter - for default)
            * 0: no filter
            * 1: low-pass filter of intensity data, highly recommended when no oversampling used (default)
    
    pflag:
        print flag (enter - for default)
            * 0: print offset summary (default)
            * 1: print all offset data
    
    pltflg:
        plotting flag (enter - for default)
            * 0: none (default)
            * 1: screen output
            * 2: screen output and PNG format plots
            * 3: output plots in PDF format
    
    ccs:
        (output) cross-correlation standard deviation of each patch (FLOAT) (enter - for none)
            * NOTE: ScanSAR and TOPS data need to be previously deramped
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/offset_pwr', SLC1, SLC2, SLC1_par, SLC2_par, OFF_par, offs, ccp, rwin, azwin, offsets, n_ovr, nr, naz, thres, lanczos, bw_frac, deramp, int_filt, pflag, pltflg, ccs]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def offset_pwr_tracking(SLC1, SLC2, SLC1_par, SLC2_par, OFF_par, offs, ccp, rwin='-', azwin='-', offsets='-', n_ovr='-', thres='-', rstep='-', azstep='-', rstart='-', rstop='-', azstart='-', azstop='-', lanczos='-', bw_frac='-', deramp='-', int_filt='-', pflag='-', pltflg='-', ccs='-', logpath=None, outdir=None, shellscript=None):
    """
    | Offset tracking between SLC images using intensity cross-correlation
    | Copyright 2023, Gamma Remote Sensing, v6.4 clw/cm 18-Apr-2023
    
    Parameters
    ----------
    SLC1:
        (input) single-look complex image 1 (reference)
    SLC2:
        (input) single-look complex image 2
    SLC1_par:
        (input) SLC1 ISP image parameter file
    SLC2_par:
        (input) SLC2 ISP image parameter file
    OFF_par:
        (input) ISP offset/interferogram parameter file
    offs:
        (output) offset estimates in range and azimuth (FCOMPLEX)
    ccp:
        (output) cross-correlation of each patch (0.0->1.0) (FLOAT)
    rwin:
        range patch size (range pixels, enter - for default from offset parameter file)
    azwin:
        azimuth patch size (azimuth lines, enter - for default from offset parameter file)
    offsets:
        (output) range and azimuth offsets and cross-correlation data in text format, enter - for no output
    n_ovr:
        SLC oversampling factor (integer 2\\*\\*N (1,2,4), enter - for default: 2)
    thres:
        cross-correlation threshold (0.0->1.0) (enter - for default from offset parameter file)
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
        offset to ending azimuth line (enter - for default: nlines-1)
    lanczos:
        Lanczos interpolator order 5 -> 9 (enter - for default: 5)
    bw_frac:
        bandwidth fraction of low-pass filter on complex data (0.0->1.0) (enter - for default: 1.0)
    deramp:
        deramp SLC phase flag (enter - for default)
            * 0: no deramp (Doppler centroid close to 0) (default)
            * 1: deramp SLC phase
    
    int_filt:
        intensity low-pass filter flag (enter - for default)
            * 0: no filter
            * 1: low-pass filter of intensity data, highly recommended when no oversampling used (default)
    
    pflag:
        print flag (enter - for default)
            * 0: print offset summary (default)
            * 1: print all offset data
    
    pltflg:
        plotting flag (enter - for default)
            * 0: none (default)
            * 1: screen output
            * 2: screen output and PNG format plots
            * 3: output plots in PDF format
    
    ccs:
        (output) cross-correlation standard deviation of each patch (FLOAT) (enter - for none)
            * NOTE: ScanSAR and TOPS data need to be previously deramped
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/offset_pwr_tracking', SLC1, SLC2, SLC1_par, SLC2_par, OFF_par, offs, ccp, rwin, azwin, offsets, n_ovr, thres, rstep, azstep, rstart, rstop, azstart, azstop, lanczos, bw_frac, deramp, int_filt, pflag, pltflg, ccs]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def offset_pwr_tracking2(SLC1, SLC2, SLC1_par, SLC2_par, OFF_par, offs, ccp, OFF_par2='-', offs2='-', rwin='-', azwin='-', offsets='-', n_ovr='-', thres='-', rstep='-', azstep='-', rstart='-', rstop='-', azstart='-', azstop='-', bw_frac='-', deramp='-', int_filt='-', pflag='-', pltflg='-', ccs='-', logpath=None, outdir=None, shellscript=None):
    """
    | Intensity cross-correlation offset tracking with the initial offset for each patch determined from input offset map
    | Copyright 2023, Gamma Remote Sensing, v2.1 clw/cm 18-Apr-2023
    
    Parameters
    ----------
    SLC1:
        (input) single-look complex image 1 (reference)
    SLC2:
        (input) single-look complex image 2
    SLC1_par:
        (input) SLC1 ISP image parameter file
    SLC2_par:
        (input) SLC2 ISP image parameter file
    OFF_par:
        (input) ISP offset/interferogram parameter file
    offs:
        (output) offset estimates in range and azimuth (FCOMPLEX)
    ccp:
        (output) cross-correlation of each patch (0.0->1.0) (FLOAT)
    OFF_par2:
        (input) ISP offset/interferogram parameter file of the offset map to determine initial offsets (enter - for none)
    offs2:
        (input) input range and azimuth offset map to determine initial offsets (enter - for none)
    rwin:
        range patch size (range pixels, enter - for default from offset parameter file)
    azwin:
        azimuth patch size (azimuth lines, enter - for default from offset parameter file)
    offsets:
        (output) range and azimuth offsets and cross-correlation data in text format, enter - for no output
    n_ovr:
        SLC oversampling factor (integer 2\\*\\*N (1,2,4), enter - for default: 2)
    thres:
        cross-correlation threshold (0.0->1.0) (enter - for default from offset parameter file)
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
        offset to ending azimuth line (enter - for default: nlines-1)
    bw_frac:
        bandwidth fraction of low-pass filter on complex data (0.0->1.0) (enter - for default: 1.0)
    deramp:
        deramp SLC phase flag (enter - for default)
            * 0: no deramp (Doppler centroid close to 0) (default)
            * 1: deramp SLC phase
    
    int_filt:
        intensity low-pass filter flag (enter - for default)
            * 0: no filter
            * 1: low-pass filter of intensity data, highly recommended when no oversampling used (default)
    
    pflag:
        print flag (enter - for default)
            * 0: print offset summary (default)
            * 1: print all offset data
    
    pltflg:
        plotting flag (enter - for default)
            * 0: none (default)
            * 1: screen output
            * 2: screen output and PNG format plots
            * 3: output plots in PDF format
    
    ccs:
        (output) cross-correlation standard deviation of each patch (FLOAT) (enter - for none)
            * NOTE: ScanSAR and TOPS data need to be previously deramped
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/offset_pwr_tracking2', SLC1, SLC2, SLC1_par, SLC2_par, OFF_par, offs, ccp, OFF_par2, offs2, rwin, azwin, offsets, n_ovr, thres, rstep, azstep, rstart, rstop, azstart, azstop, bw_frac, deramp, int_filt, pflag, pltflg, ccs]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def offset_pwr_tracking_polygons(SLC_par, OFF_par, rlks, azlks, rwin, azwin, polygons, rstep='-', azstep='-', rstart='-', rstop='-', azstart='-', azstop='-', rb='-', azb='-', logpath=None, outdir=None, shellscript=None):
    """
    | Offset tracking polygon calculation in MLI coordinates
    | Copyright 2023, Gamma Remote Sensing, v1.2 18-Apr-2023 cw
    
    Parameters
    ----------
    SLC_par:
        (input) reference SLC ISP image parameter file
    OFF_par:
        (input/output) ISP offset/interferogram parameter file
    rlks:
        range decimation factor for MLI geometry  (enter - for default: 1)
    azlks:
        azimuth decimation factor for the MLI geometry (enter - for default: 1)
    rwin:
        range patch size (range pixels, enter - for default from offset parameter file)
    azwin:
        azimuth patch size (azimuth lines, enter - for default from offset parameter file)
    polygons:
        (output) polygon vertices in text format
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
        offset to ending azimuth line (enter - for default: nlines-1)
    rb:
        polygon range border in MLI samples: (enter - for default: 7)
    azb:
        polygon azimuth border in MLI lines: (enter - for default: 7)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/offset_pwr_tracking_polygons', SLC_par, OFF_par, rlks, azlks, rwin, azwin, polygons, rstep, azstep, rstart, rstop, azstart, azstop, rb, azb]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def offset_SLC(SLC1, SLC2, SLC1_par, SLC2_par, OFF_par, offs, snr, rwin='-', azwin='-', offsets='-', n_ovr='-', nr='-', naz='-', thres='-', ISZ='-', pflag='-', logpath=None, outdir=None, shellscript=None):
    """
    | Offsets between SLC images using fringe visibility
    | Copyright 2023, Gamma Remote Sensing, v3.1 18-Apr-2023 clw
    
    Parameters
    ----------
    SLC1:
        (input) single-look complex image 1 (reference)
    SLC2:
        (input) single-look complex image 2
    SLC1_par:
        (input) SLC1 ISP image parameter file
    SLC2_par:
        (input) SLC2 ISP image parameter file
    OFF_par:
        (input) ISP offset/interferogram parameter file
    offs:
        (output) offset estimates (FCOMPLEX)
    snr:
        (output) offset estimation SNR (FLOAT)
    rwin:
        search window size (range pixels) (enter - for default from offset parameter file)
    azwin:
        search window size (azimuth lines) (enter - for default from offset parameter file)
    offsets:
        (output) range and azimuth offsets and SNR data in text format, enter - for no output
    n_ovr:
        SLC oversampling factor (integer 2\\*\\*N (1,2,4) enter - for default: 2)
    nr:
        number of offset estimates in range direction (enter - for default from offset parameter file)
    naz:
        number of offset estimates in azimuth direction (enter - for default from offset parameter file)
    thres:
        offset estimation quality threshold (enter - for default from offset parameter file)
    ISZ:
        search chip interferogram size (in non-oversampled pixels, enter - for default: 16)
    pflag:
        print flag (enter - for default)
            * 0: print offset summary (default)
            * 1: print all offset data
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/offset_SLC', SLC1, SLC2, SLC1_par, SLC2_par, OFF_par, offs, snr, rwin, azwin, offsets, n_ovr, nr, naz, thres, ISZ, pflag]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def offset_SLC_tracking(SLC1, SLC2, SLC1_par, SLC2_par, OFF_par, offs, snr, rsw='-', azsw='-', offsets='-', n_ovr='-', thres='-', rstep='-', azstep='-', rstart='-', rstop='-', azstart='-', azstop='-', ISZ='-', pflag='-', logpath=None, outdir=None, shellscript=None):
    """
    | Offset tracking between SLC images using fringe visibility
    | Copyright 2023, Gamma Remote Sensing, v3.8 18-Apr-2023 clw
    
    Parameters
    ----------
    SLC1:
        (input) single-look complex image 1 (reference)
    SLC2:
        (input) single-look complex image 2
    SLC1_par:
        (input) SLC1 ISP image parameter file
    SLC2_par:
        (input) SLC2 ISP image parameter file
    OFF_par:
        (input) ISP offset/interferogram parameter file
    offs:
        (output) offset estimates (FCOMPLEX)
    snr:
        (output) offset estimation SNR (FLOAT)
    rsw:
        range search window size (range pixels) (enter - for default from offset parameter file)
    azsw:
        azimuth search window size (azimuth lines) (enter - for default from offset parameter file)
    offsets:
        (output) range and azimuth offsets and SNR data in text format, enter - for no output
    n_ovr:
        SLC over-sampling factor (integer 2\\*\\*N (1,2,4) enter - for default: 2)
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
        ending azimuth line (enter - for default: nlines - azsw/2)
    ISZ:
        search chip interferogram size (in non-oversampled pixels, enter - for default: 16)
    pflag:
        print flag (enter - for default)
            * 0: print offset summary (default)
            * 1: print all offset data
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/offset_SLC_tracking', SLC1, SLC2, SLC1_par, SLC2_par, OFF_par, offs, snr, rsw, azsw, offsets, n_ovr, thres, rstep, azstep, rstart, rstop, azstart, azstop, ISZ, pflag]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def offset_sub(offs, OFF_par, offs_sub, logpath=None, outdir=None, shellscript=None):
    """
    | Subtraction of polynomial from range and azimuth offset estimates
    | Copyright 2017, Gamma Remote Sensing, v1.0 27-Mar-2017 cm
    
    Parameters
    ----------
    offs:
        (input) range and azimuth offset estimates (fcomplex)
    OFF_par:
        (input) ISP offset/interferogram parameter file
    offs_sub:
        (output) range and azimuth offset estimates after polynomial subtraction (fcomplex)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/offset_sub', offs, OFF_par, offs_sub]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def offset_tracking(offs, ccp, SLC_par, OFF_par, disp_map, disp_val='-', mode='-', thres='-', poly_flag='-', logpath=None, outdir=None, shellscript=None):
    """
    | Conversion of range and azimuth offsets files to displacement map
    | Copyright 2017, Gamma Remote Sensing, v2.0 4-Apr-2017 ts/clw/uw
    
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
        (output) range and azimuth displacement estimates and cross-correlation values (enter - for none) (text)
    mode:
        flag indicating displacement mode:
            * 0: displacement in range and azimuth pixels
            * 1: displacement in meters in slant range and azimuth directions
            * 2: displacement in meters in ground range and azimuth directions (default)
    
    thres:
        cross-correlation threshold to accept offset value (default from OFF_par)
    poly_flag:
        flag indicating if trend calculated using offset polynomials from OFF_par is subtracted:
            * 0: do not subtract polynomial trend from offset data
            * 1: subtract polynomial trend from offset data (default)
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/offset_tracking', offs, ccp, SLC_par, OFF_par, disp_map, disp_val, mode, thres, poly_flag]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def ORB_filt(SLC_par_in, SLC_par_out, interval='-', extra='-', logpath=None, outdir=None, shellscript=None):
    """
    | Filter state vectors using a least-squares polynomial model
    | Copyright 2020, Gamma Remote Sensing, v1.3 20-May-2020 clw/cm
    
    Parameters
    ----------
    SLC_par_in:
        (input) ISP image parameter file at least 5 state vectors
    SLC_par_out:
        (output) ISP image parameter file with state vectors filtered using least-squares
    interval:
        time interval between state vectors (enter - for default: state vector time interval in SLC_par)
    extra:
        extra time for state vectors at start and end of image (sec.) (enter - for default: 5.0)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/ORB_filt', SLC_par_in, SLC_par_out, interval, extra]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def ORB_prop_SLC(SLC_par, nstate='-', interval='-', extra='-', mode='-', logpath=None, outdir=None, shellscript=None):
    """
    | Calculate state vectors using orbit propagation and interpolation
    | Copyright 2022, Gamma Remote Sensing, v2.0 1-Feb-2022 clw/awi/cm
    
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
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/ORB_prop_SLC', SLC_par, nstate, interval, extra, mode]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def ORRM_vec(SLC_par, ORRM, nstate='-', logpath=None, outdir=None, shellscript=None):
    """
    | Calculate state vectors extraction from ORRM file
    | Copyright 2023, Gamma Remote Sensing, v1.5 19-Apr-2023 clw
    
    Parameters
    ----------
    SLC_par:
        (input/output) ISP SLC/MLI image parameter file
    ORRM:
        (input) ORRM state vector file
    nstate:
        number of state vectors (enter - for default: 5, maximum: 1024)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/ORRM_vec', SLC_par, ORRM, nstate]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_ACS_ERS(CEOS_SAR_leader, SLC_par, logpath=None, outdir=None, shellscript=None):
    """
    | ISP parameter file generation for ERS SLC data from the ACS processor
    | Copyright 2020, Gamma Remote Sensing, v1.4 3-Sep-2020 clw/uw/cm
    
    Parameters
    ----------
    CEOS_SAR_leader:
        (input) ERS CEOS SAR leader file
    SLC_par:
        (output) ISP SLC parameter file (example <orbit>.slc.par)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_ACS_ERS', CEOS_SAR_leader, SLC_par]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_ASAR(ASAR_ERS_file, output_name, K_dB='-', logpath=None, outdir=None, shellscript=None):
    """
    | Extract SLC/MLI image parameters and images from ENVISAT ASAR SLC, WSS, APP, and PRI products
    | Copyright 2023, Gamma Remote Sensing, v2.9 20-Oct-2023 clw/uw/awi/cm
    
    Parameters
    ----------
    ASAR_ERS_file:
        (input) ASAR or ERS data in ASAR format (SAR_IMS_1P) including header and image as provided by ESA
    output_name:
        (output) common part of output file names (e.g. YYYMMDD date)
    K_dB:
        Calibration factor in dB (nominal value for all ASAR modes: 55.0)
            * NOTE: Use - to use the calibration factor provided in the ASAR file header
            * NOTE: In the case that a calibration factor is specified on the command line, PRI images are converted
              to radiometrically calibrated ground-range intensity images in float format
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_ASAR', ASAR_ERS_file, output_name, K_dB]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_ASF_91(CEOS_leader, CEOS_trailer, SLC_par, logpath=None, outdir=None, shellscript=None):
    """
    | SLC parameter file for data data from theAlaska SAR Facility (1991-1996)
    | Copyright 2020, Gamma Remote Sensing, v3.4 3-Sep-2020 clw/uw/cm
    
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
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_ASF_91', CEOS_leader, CEOS_trailer, SLC_par]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_ASF_96(CEOS_SAR_leader, SLC_par, logpath=None, outdir=None, shellscript=None):
    """
    | ISP parameter file for ASF data 1996-->present v1.1
    | Copyright 2020, Gamma Remote Sensing, v1.4 3-Sep-2020 clw/uw/cm
    
    Parameters
    ----------
    CEOS_SAR_leader:
        (input) CEOS SAR leader file
    SLC_par:
        (output) ISP SLC parameter file (example <orbit>.slc.par)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_ASF_96', CEOS_SAR_leader, SLC_par]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_ASF_PRI(CEOS_leader, CEOS_data, GRD_par, GRD, logpath=None, outdir=None, shellscript=None):
    """
    | ISP parameter file for ASF detected ground range images (L1) Sep 1996 --> present
    | Copyright 2021, Gamma Remote Sensing, v1.5 14-Jun-2021 clw/uw/cm
    
    Parameters
    ----------
    CEOS_leader:
        (input) CEOS leader file
    CEOS_data:
        (input) CEOS data file binary
    GRD_par:
        (output) ISP ground range image parameter file
    GRD:
        (output) ISP ground range image (enter - for none, FLOAT intensity)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_ASF_PRI', CEOS_leader, CEOS_data, GRD_par, GRD]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_ASF_RSAT_SS(CEOS_leader, CEOS_data, GRD_par, GRD, logpath=None, outdir=None, shellscript=None):
    """
    | ISP parameter file for ASF Radarsat-1 SCANSAR images
    | Copyright 2020, Gamma Remote Sensing, v1.1 3-Sep-2020 clw/uw/cm
    
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
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_ASF_RSAT_SS', CEOS_leader, CEOS_data, GRD_par, GRD]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_ASF_SLC(CEOS_leader, SLC_par, CEOS_data='-', SLC='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate SLC image parameter file and reformat data
    | Copyright 2023, Gamma Remote Sensing, v1.1 18-Apr-2023 clw/uw
    
    Parameters
    ----------
    CEOS_leader:
        (input) CEOS SAR leader file
    SLC_par:
        (output) ISP SLC parameter file (example <date>.slc.par)
    CEOS_data:
        (input) CEOS data file (example: dat_01.001) (enter - for none)
    SLC:
        (output) SLC data with file and line headers removed (example: <date>.slc) (enter - for none)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_ASF_SLC', CEOS_leader, SLC_par, CEOS_data, SLC]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_ASNARO2(CEOS_data, CEOS_leader, SLC_par, SLC='-', reramp='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate SLC parameter and image files for ASNARO-2 Spotlight, Stripmap and ScanSAR level 1.1 data
    | Copyright 2023, Gamma Remote Sensing, v1.4 15-Jun-2023 cm/uw
    
    Parameters
    ----------
    CEOS_data:
        (input) CEOS format SLC data (IMG-PP-AS2\\*)
    CEOS_leader:
        (input) CEOS SAR leader file for ASNARO-2 data (LED-AS2\\*)
    SLC_par:
        (output) ISP SLC parameter file (example: yyyymmdd_pp.slc.par)
    SLC:
        (output) SLC (Spotlight and Stripmap) or SLI (ScanSAR) data file (enter - for none, example: yyyymmdd_pp.slc)
    reramp:
        reramp SLC phase flag (enter - for default)
            * 0: no reramp
            * 1: reramp SLC phase (default)
            * NOTE: ASNARO2 geocoded and georeferenced data in GeoTIFF format (level 1.5) can be read using par_ASNARO2_geo program.
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_ASNARO2', CEOS_data, CEOS_leader, SLC_par, SLC, reramp]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_ATLSCI_ERS(CEOS_SAR_leader, CEOS_Image, SLC_par, logpath=None, outdir=None, shellscript=None):
    """
    | ISP parameter file for ATL-SCI ERS SLC data
    | Copyright 2020, Gamma Remote Sensing, v2.9 21-Sep-2020 clw/cm
    
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
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_ATLSCI_ERS', CEOS_SAR_leader, CEOS_Image, SLC_par]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_Capella_SLC(GeoTIFF, ext_JSON, SLC_par, SLC='-', radcal='-', noise='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate SLC parameter and image files for Capella SLC data
    | Copyright 2025, Gamma Remote Sensing, v2.0 28-Apr-2025 cm
    
    Parameters
    ----------
    GeoTIFF:
        (input) Capella image data file in GeoTIFF format (\\*.tif)
    ext_JSON:
        (input) Capella extended metadata file in JSON format (\\*_extended.json)
    SLC_par:
        (output) ISP SLC parameter file (example: yyyymmdd.slc.par)
    SLC:
        (output) SLC data file (enter - for none, example: yyyymmdd.slc)
    radcal:
        radiometric calibration flag (enter - for default)
            * 0: beta0 (default)
            * 1: sigma0
    
    noise:
        noise levels flag (enter - for default)
            * 0: do not use noise levels (default)
            * 1: use noise levels
            * NOTE: Capella terrain geocoded data in GeoTIFF format can be read using par_Capella_geo program
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_Capella_SLC', GeoTIFF, ext_JSON, SLC_par, SLC, radcal, noise]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_CS_DGM(HDF5, trunk, logpath=None, outdir=None, shellscript=None):
    """
    | Generate ISP MLI parameter and image files for COSMO-Skymed DGM data
    | Copyright 2024, Gamma Remote Sensing, v1.1 12-Sep-2024 cm/awi/ms/cw/uw
    
    Parameters
    ----------
    HDF5:
        (input) COSMO-Skymed DGM data file in HDF5 format
    trunk:
        (output) output file name trunk used for output filenames 
            (example: yyyymmdd -> yyyymmdd_pol.mli yyyymmdd_pol.mli.par)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_CS_DGM', HDF5, trunk]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_CS_SLC(HDF5, trunk, logpath=None, outdir=None, shellscript=None):
    """
    | Generate ISP SLC parameter and image files for Cosmo-Skymed SCS data
    | Copyright 2024, Gamma Remote Sensing, v2.2 12-Sep-2024 awi/ms/cw/uw
    
    Parameters
    ----------
    HDF5:
        (input) SCS data file in HDF5 format
    trunk:
        (output) output file name trunk used for output filenames 
            (example: yyyymmdd -> yyyymmdd_pol_beamid.slc yyyymmdd_pol_beamid.slc.par)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_CS_SLC', HDF5, trunk]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_CS_SLC_TIF(GeoTIFF, XML, trunk, logpath=None, outdir=None, shellscript=None):
    """
    | Generate ISP SLC parameter and image files for Cosmo Skymed SCS data in GeoTIFF format
    | Copyright 2023, Gamma Remote Sensing, v1.6 16-May-2023 awi/ms/clw/cm
    
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
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_CS_SLC_TIF', GeoTIFF, XML, trunk]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_CSG_DGM(HDF5, trunk, logpath=None, outdir=None, shellscript=None):
    """
    | Generate ISP MLI parameter and image files for COSMO-Skymed Second Generation DGM data
    | Copyright 2024, Gamma Remote Sensing, v1.1 12-Sep-2024 cm/awi/ms/cw/uw
    
    Parameters
    ----------
    HDF5:
        (input) COSMO-Skymed Second Generation DGM data file in HDF5 format
    trunk:
        (output) output file name trunk used for output filenames 
            (example: yyyymmdd -> yyyymmdd_pol.mli yyyymmdd_pol.mli.par)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_CSG_DGM', HDF5, trunk]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_CSG_SLC(HDF5, trunk, logpath=None, outdir=None, shellscript=None):
    """
    | Generate ISP SLC parameter and image files for COSMO-Skymed Second Generation SCS data
    | Copyright 2024, Gamma Remote Sensing, v1.3 12-Sep-2024 cm/awi/ms/cw/uw
    
    Parameters
    ----------
    HDF5:
        (input) SCS data file in HDF5 format
    trunk:
        (output) output file name trunk used for output filenames 
            (example: yyyymmdd -> yyyymmdd_pol_beamid.slc yyyymmdd_pol_beamid.slc.par)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_CSG_SLC', HDF5, trunk]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_EORC_JERS_SLC(CEOS_SAR_leader, SLC_par, CEOS_data='-', SLC='-', logpath=None, outdir=None, shellscript=None):
    """
    | Reformat EORC processed JERS-1 SLC and generate the ISP parameter file
    | Copyright 2023, Gamma Remote Sensing, v1.6 18-Apr-2023 clw/cm
    
    Parameters
    ----------
    CEOS_SAR_leader:
        (input) CEOS SAR leader file for JERS SLC processed by EORC
    SLC_par:
        (output) ISP image parameter file
    CEOS_data:
        (input) CEOS format SLC data (IMOP_01.DAT, enter - for none)
    SLC:
        (output) reformated JERS SLC (example: yyyymmdd.SLC, enter - for none)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_EORC_JERS_SLC', CEOS_SAR_leader, SLC_par, CEOS_data, SLC]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_EORC_PALSAR(CEOS_leader, SLC_par, CEOS_data, SLC='-', dtype='-', sc_dB='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate SLC image and parameter files for PALSAR, PALSAR-2, and PALSAR-3 level 1.1 SLC data produced by EORC/JAXA and ESA
    | Copyright 2025, Gamma Remote Sensing, v3.9 12-Jun-2025 clw/cm
    
    Parameters
    ----------
    CEOS_leader:
        (input) CEOS leader file for PALSAR, PALSAR-2, or PALSAR-3 Level 1.1 SLC data (LED...)
    SLC_par:
        (output) ISP image parameter file (example: yyyymmdd.slc.par)
    CEOS_data:
        (input) PALSAR CEOS format Level 1.1 SLC (IMG...)
    SLC:
        (output) reformatted PALSAR SLC (example: yyyymmdd.slc, enter - for none)
    dtype:
        output data type (enter - for default)
            * 0: FCOMPLEX (default)
            * 1: SCOMPLEX
    
    sc_dB:
        scale factor for FCOMPLEX -> SCOMPLEX, (enter - for default: HH,VV (dB): 60.0000, VH,HV: 70.0000)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_EORC_PALSAR', CEOS_leader, SLC_par, CEOS_data, SLC, dtype, sc_dB]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_EORC_PALSAR_ScanSAR(CEOS_data, CEOS_leader, SLC_par, SLC='-', TOPS_par='-', afmrate='-', shift='-', reramp='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate SLC parameter and image files from EORC PALSAR2 ScanSAR burst SLC data in CEOS format
    | Copyright 2023, Gamma Remote Sensing, v1.4 18-Apr-2023 cm/clw
    
    Parameters
    ----------
    CEOS_data:
        (input) CEOS image file for a PALSAR2 ScanSAR burst data subswath (IMG...)
    CEOS_leader:
        (input) CEOS leader file for PALSAR2 ScanSAR burst data (LED...)
    SLC_par:
        (output) ISP image parameter file (example: yyyymmdd_b1_hh.slc.par)
    SLC:
        (output) SLC data file (enter - for none, example: yyyymmdd_b1_hh.slc)
    TOPS_par:
        (output) SLC burst annotation file (enter - for none, example: yyyymmdd_b1_hh.slc.tops_par)
    afmrate:
        azimuth FM rate estimation method (enter - for default)
            * 0: beam velocity on the ground
            * 1: platform velocity (default)
    
    shift:
        shift azimuth spectrum by fs/2 (enter - for default)
            * 0: no
            * 1: yes (default)
    
    reramp:
        reramp data using Doppler centroid and azimuth FM rate estimate (enter - for default)
            * 0: no
            * 1: yes (default)
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_EORC_PALSAR_ScanSAR', CEOS_data, CEOS_leader, SLC_par, SLC, TOPS_par, afmrate, shift, reramp]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_ERSDAC_PALSAR(ERSDAC_SLC_par, SLC_par, logpath=None, outdir=None, shellscript=None):
    """
    | Generate the ISP image parameter file from ERSDAC PALSAR level 1.1 SLC data
    | Copyright 2023, Gamma Remote Sensing, v1.7 5-Jun-2023 clw
    
    Parameters
    ----------
    ERSDAC_SLC_par:
        (input) ERSDAC SLC parameter file Level 1.1 (PASL11\\*.SLC.par)
    SLC_par:
        (output) ISP image parameter file (example: yyyymmdd.slc.par)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_ERSDAC_PALSAR', ERSDAC_SLC_par, SLC_par]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_ESA_ERS(CEOS_SAR_leader, SLC_par, inlist, CEOS_DAT='-', SLC='-', logpath=None, outdir=None, shellscript=None):
    """
    | ISP parameter file generation for ERS SLC data from the PGS, VMP, and SPF processors
    | Copyright 2020, Gamma Remote Sensing, v1.5 21-Sep-2020 clw/uw/cm
    
    Parameters
    ----------
    CEOS_SAR_leader:
        (input) ERS CEOS SAR leader file
    SLC_par:
        (output) ISP SLC parameter file (example: <date>.slc.par)
    inlist:
        a list of arguments to be passed to stdin
    CEOS_DAT:
        (input) CEOS data file (example: DAT_01.001)
    SLC:
        (output) SLC data with file and line headers removed (example: <date>.slc)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_ESA_ERS', CEOS_SAR_leader, SLC_par, CEOS_DAT, SLC]
    process(cmd, logpath=logpath, outdir=outdir, inlist=inlist, shellscript=shellscript)


def par_ESA_JERS_SEASAT_SLC(CEOS_data, CEOS_leader, SLC_par, SLC='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate SLC parameter and image files for ESA-provided JERS and SEASAT SLC data
    | Copyright 2023, Gamma Remote Sensing, v1.4 15-Jun-2023 cm/clw/ts
    
    Parameters
    ----------
    CEOS_data:
        (input) CEOS format SLC data (DAT_01.001)
    CEOS_leader:
        (input) CEOS SAR leader file for JERS SLC processed by ESA (LEA_01.001)
    SLC_par:
        (output) ISP SLC parameter file (example: yyyymmdd.slc.par)
    SLC:
        (output) SLC data file (enter - for none, example: yyyymmdd.slc)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_ESA_JERS_SEASAT_SLC', CEOS_data, CEOS_leader, SLC_par, SLC]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_ESA_PALSAR_GDH(CEOS_data, CEOS_leader, MLI_par, MLI='-', GRD_par='-', GRD='-', rps='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate MLI and GRD image and parameter files for PALSAR + PALSAR2 level 1.5 GDH data provided by ESA
    | Copyright 2023, Gamma Remote Sensing, v1.4 5-Jun-2023 clw/cm
    
    Parameters
    ----------
    CEOS_data:
        (input) CEOS image file for PALSAR or PALSAR-2 Level 1.5 GDH data (IMG...)
    CEOS_leader:
        (input) CEOS leader file for PALSAR or PALSAR-2 Level 1.5 GDH data (LED...)
    MLI_par:
        (output) MLI parameter file (example: yyyymmdd_pp.mli.par)
    MLI:
        (output) MLI data file in slant range geometry (example: yyyymmdd_pp.mli, enter - for none)
    GRD_par:
        (output) GRD parameter file (example: yyyymmdd_pp.grd.par, enter - for none)
    GRD:
        (output) GRD data file (example: yyyymmdd_pp.grd, enter - for none)
    rps:
        slant range pixel spacing (m) (enter - for default: calculated from ground-range parameters)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_ESA_PALSAR_GDH', CEOS_data, CEOS_leader, MLI_par, MLI, GRD_par, GRD, rps]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_Fucheng_SLC(GeoTIFF, annotation_XML, calibration_XML, noise_XML, SLC_par, SLC='-', dtype='-', radcal='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate SLC parameter and image files for Spacety Fucheng SLC data
    | Copyright 2024, Gamma Remote Sensing, v1.1 7-Jun-2024 cm/clw/awi
    
    Parameters
    ----------
    GeoTIFF:
        (input) image data file in \\*.tiff GeoTIFF format (enter - for default: none)
    annotation_XML:
        (input) Fucheng XML annotation file
    calibration_XML:
        (input) Fucheng radiometric calibration XML file to generate output as sigma0
            (enter - for default: return uncalibrated digital numbers)
    noise_XML:
        (input) Fucheng noise XML file (enter - for default: no subtraction of thermal noise power)
    SLC_par:
        (output) SLC parameter file (e.g.: yyyymmdd_vv.slc.par)
    SLC:
        (output) SLC data file (enter - for default: none, e.g.: yyyymmdd_vv.slc)
    dtype:
        output data type (enter - for default)
            * 0: FCOMPLEX (default)
            * 1: SCOMPLEX
    
    radcal:
        radiometric calibration flag (enter - for default)
            * 0: none
            * 1: Beta Nought
            * 2: Sigma Nought (default)
            * 3: Gamma
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_Fucheng_SLC', GeoTIFF, annotation_XML, calibration_XML, noise_XML, SLC_par, SLC, dtype, radcal]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_GF3_SLC(GeoTIFF, annotation_XML, SLC_par, SLC='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate SLC parameter file and SLC image from a Gaofen-3 data set in GeoTIFF format
    | Copyright 2023, Gamma Remote Sensing, v1.3 14-Jun-2023 cm
    
    Parameters
    ----------
    GeoTIFF:
        (input) Gaofen-3 data file in GeoTIFF format (\\*.tiff) (enter - for none)
    annotation_XML:
        (input) Gaofen-3 annotation file in XML format (\\*.meta.xml)
    SLC_par:
        (output) ISP SLC parameter file (example: yyyymmdd.slc.par)
    SLC:
        (output) ISP SLC data file (example: yyyymmdd.slc) (enter - for none, SLC output will not be produced)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_GF3_SLC', GeoTIFF, annotation_XML, SLC_par, SLC]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_HISEA1_SLC(GeoTIFF, annotation_XML, calibration_XML, SLC_par, SLC='-', dtype='-', sc_dB='-', shift='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate SLC parameter and image files for Hisea-1 SLC data
    | Copyright 2023, Gamma Remote Sensing, v1.4 11-May-2023 awi/cm
    
    Parameters
    ----------
    GeoTIFF:
        (input) image data file in GeoTIFF format (enter - for none, \\*.tiff)
    annotation_XML:
        (input) Hisea-1 L1 XML annotation file
    calibration_XML:
        (input) Hisea-1 L1 radiometric calibration XML file (enter - for no radiometric calibration)
    SLC_par:
        (output) ISP SLC parameter file (example: yyyymmdd_vv.slc.par)
    SLC:
        (output) SLC data file (enter - for none, example: yyyymmdd_vv.slc)
    dtype:
        output data type (enter - for default)
            * 0: FCOMPLEX (default)
            * 1: SCOMPLEX
    
    sc_dB:
        scale factor for FCOMPLEX -> SCOMPLEX, (enter - for default: HH,VV (dB): 60.0000,  VH,HV: 70.0000)
    shift:
        shift azimuth spectrum by fs/2 (enter - for default)
            * 0: no
            * 1: yes (default)
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_HISEA1_SLC', GeoTIFF, annotation_XML, calibration_XML, SLC_par, SLC, dtype, sc_dB, shift]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_HT1_SLC(GeoTIFF, annotation_XML, SLC_par, SLC='-', dtype='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate SLC parameter and image files for HT1 / Hongtu-1 / PIESAT-1 SLC data
    | Copyright 2024, Gamma Remote Sensing, v1.0 5-Jun-2024 cm/clw/awi
    
    Parameters
    ----------
    GeoTIFF:
        (input) image data file in \\*.tiff GeoTIFF format (enter - for default: none)
    annotation_XML:
        (input) HT1 XML annotation file
    SLC_par:
        (output) SLC parameter file (e.g.: yyyymmdd_vv.slc.par)
    SLC:
        (output) SLC data file (enter - for default: none, e.g.: yyyymmdd_vv.slc)
    dtype:
        output data type (enter - for default)
            * 0: FCOMPLEX (default)
            * 1: SCOMPLEX
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_HT1_SLC', GeoTIFF, annotation_XML, SLC_par, SLC, dtype]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_ICEYE_GRD(GeoTIFF, XML, MLI_par, MLI='-', GRD_par='-', GRD='-', rps='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate MLI and GRD image and parameter files for ICEYE GRD data
    | Copyright 2024, Gamma Remote Sensing, v1.4 13-Jun-2024 cm
    
    Parameters
    ----------
    GeoTIFF:
        (input) ICEYE GRD data file in GeoTIFF format (enter - for none, \\*.tif)
    XML:
        (input) ICEYE XML annotation file
    MLI_par:
        (output) MLI parameter file (example: yyyymmdd.mli.par)
    MLI:
        (output) MLI data file in slant range geometry (example: yyyymmdd.mli, enter - for none)
    GRD_par:
        (output) GRD parameter file (example: yyyymmdd.grd.par, enter - for none)
    GRD:
        (output) GRD data file (example: yyyymmdd.grd, enter - for none)
    rps:
        slant range pixel spacing (m) (enter - for default: calculated from ground-range parameters)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_ICEYE_GRD', GeoTIFF, XML, MLI_par, MLI, GRD_par, GRD, rps]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_ICEYE_SLC(HDF5, SLC_par, SLC='-', dtype='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate ISP SLC parameter and binary files for ICEYE SLC data
    | Copyright 2024, Gamma Remote Sensing, v1.9 28-Oct-2024 cm
    
    Parameters
    ----------
    HDF5:
        (input) ICEYE SLC data file in HDF5 format
    SLC_par:
        (output) ISP SLC parameter file (example: yyyymmdd.slc.par)
    SLC:
        (output) SLC data file (enter - for none, example: yyyymmdd.slc)
    dtype:
        output data type (enter - for default: same as input)
            * 0: FCOMPLEX
            * 1: SCOMPLEX
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_ICEYE_SLC', HDF5, SLC_par, SLC, dtype]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_IECAS_SLC(aux_data, slc_Re, slc_Im, date, SLC_par, SLC, logpath=None, outdir=None, shellscript=None):
    """
    | Generate SLC parameter and image files for IECAS SLC data
    | Copyright 2023, Gamma Remote Sensing, v1.3 18-Apr-2023
    
    Parameters
    ----------
    aux_data:
        (input) IECAS SAR auxillary data (POS\\*.dat)
    slc_Re:
        (input) real part of complex SLC data
    slc_Im:
        (input) imaginary part of complex SLC data
    date:
        (input) acquistion date format: YYYYMMDD (example 20110121) from aux_data filename
    SLC_par:
        (output) ISP SLC parameter file
    SLC:
        (output) SLC image
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_IECAS_SLC', aux_data, slc_Re, slc_Im, date, SLC_par, SLC]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_KC_PALSAR_slr(facter_m, CEOS_leader, SLC_par, pol, pls_mode, KC_data, pwr='-', fdtab='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate ISP parameter file, Doppler table, and images for PALSAR KC Slant-Range data
    | Copyright 2023, Gamma Remote Sensing, v2.3 5-Jun-2023 ms/awi/clw/cm
    
    Parameters
    ----------
    facter_m:
        (input) PALSAR Kyoto-Carbon parameter file
    CEOS_leader:
        (input) PALSAR Kyoto-Carbon leader file (LED)
    SLC_par:
        (output) ISP image parameter file (example: yyyymmdd_pp.mli.par)
    pol:
        polarization e.g. HH or HV
    pls_mode:
        PALSAR acquisition mode:
            * 1: Fine Beam Single
            * 2: Fine Beam Double
            * 3: Wide Beam
    
    KC_data:
        (input) PALSAR Kyoto-Carbon data (named sar_Q\\*.dat_\\*)
    pwr:
        (output) PALSAR Kyoto-Carbon data strip expressed as SAR intensity (enter - for none, example: yyyymmdd_pp.mli)
    fdtab:
        (output) table of output polynomials, one polynomial/block used as input to gc_map_fd (enter - for none)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_KC_PALSAR_slr', facter_m, CEOS_leader, SLC_par, pol, pls_mode, KC_data, pwr, fdtab]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_KS_DGM(HDF5, trunk, logpath=None, outdir=None, shellscript=None):
    """
    | Generate ISP SLC parameter and PRI image files for Kompsat DGM data
    | Copyright 2023, Gamma Remote Sensing, v1.4 13-Jul-2023 awi/cm
    
    Parameters
    ----------
    HDF5:
        (input) DGM data file in HDF5 format
    trunk:
        (output) output file name trunk used for output filenames 
            (example: yyyymmdd -> yyyymmdd_pol_beamid.slc yyyymmdd_pol_beamid.pri.par)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_KS_DGM', HDF5, trunk]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_KS_SLC(HDF5, trunk, logpath=None, outdir=None, shellscript=None):
    """
    | Generate ISP SLC parameter and image files for Kompsat SCS data
    | Copyright 2023, Gamma Remote Sensing, v1.7 13-Jul-2023 awi/clw/cm
    
    Parameters
    ----------
    HDF5:
        (input) SCS data file in HDF5 format
    trunk:
        (output) output file name trunk used for output filenames 
            (example: yyyymmdd -> yyyymmdd_pol_beamid.slc yyyymmdd_pol_beamid.slc.par)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_KS_SLC', HDF5, trunk]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_LT1_SLC(GeoTIFF, annotation_XML, SLC_par, SLC='-', dtype='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate SLC parameter file and SLC image from a LT-1 data set
    | Copyright 2024, Gamma Remote Sensing, v1.3 17-Jul-2024 awi/cm
    
    Parameters
    ----------
    GeoTIFF:
        (input) image data file in GeoTIFF format (enter - for none, \\*.tiff)
    annotation_XML:
        (input) LT-1 product annotation XML file (\\*.meta.xml)
    SLC_par:
        (output) ISP SLC parameter file (example: yyyymmdd.slc.par)
    SLC:
        (output) SLC data file, example: yyyymmdd.slc (enter - for none, SLC output will not be produced)
    dtype:
        output data type (enter - for default)
            * 0: FCOMPLEX (default)
            * 1: SCOMPLEX
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_LT1_SLC', GeoTIFF, annotation_XML, SLC_par, SLC, dtype]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_MSP(SAR_par, PROC_par, SLC_par, image_format='-', logpath=None, outdir=None, shellscript=None):
    """
    | ISP image parameter file from MSP processing parameter and sensor files
    | Copyright 2024, Gamma Remote Sensing, v3.7 8-May-2024 clw/uw/of
    
    Parameters
    ----------
    SAR_par:
        (input) MSP SAR sensor parameter file
    PROC_par:
        (input) MSP processing parameter file
    SLC_par:
        (output) ISP SLC/MLI image parameter file
    image_format:
        image format flag (enter - for default: from MSP processing parameter file)
            * 0: FCOMPLEX (pairs of 4-byte float)
            * 1: SCOMPLEX (pairs of 2-byte short integer)
            * 2: FLOAT (4-bytes/value)
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_MSP', SAR_par, PROC_par, SLC_par, image_format]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_NISAR_RSLC(HDF5, root_name, radcal='-', noise='-', band='-', freq='-', pol='-', out_flag='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate ISP SLC parameter and image files for NISAR Level-1 RSLC data
    | Copyright 2025, Gamma Remote Sensing, v1.4 19-May-2025 cm
    
    Parameters
    ----------
    HDF5:
        (input) NISAR RSLC data file in HDF5 format (Level-1 Range Doppler Single Look Complex)
    root_name:
        (output) root name of the generated output files (example: yyyymmdd)
    radcal:
        radiometric calibration flag (enter - for default)
            * 0: none
            * 1: beta0
            * 2: sigma0 (default)
            * 3: gamma0
    
    noise:
        noise subtraction using noise equivalent backscatter look-up table (enter - for default)
            * 0: do not apply noise subtraction (default)
            * 1: apply noise subtraction
    
    band:
        radar band L or S (enter - for default: all available radar bands)
    freq:
        frequencies A or B in case of split imaging bands (enter - for default: all available frequencies)
    pol:
        polarization HH, HV, RH, RV, VH, or VV (enter - for default: all available polarizations)
    out_flag:
        output flag (enter - for default)
            * 0: write data and parameter files (default)
            * 1: only write parameter files
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_NISAR_RSLC', HDF5, root_name, radcal, noise, band, freq, pol, out_flag]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_NovaSAR_GRD(GeoTIFF, XML, polarization, MLI_par, MLI='-', GRD_par='-', GRD='-', rps='-', radcal='-', noise='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate MLI and GRD image and parameter files for NovaSAR GRD and SCD data
    | Copyright 2023, Gamma Remote Sensing, v1.8 3-Mar-2023 cm
    
    Parameters
    ----------
    GeoTIFF:
        (input) NovaSAR image data file in GeoTIFF format (enter - for none, \\*.tif)
    XML:
        (input) NovaSAR XML annotation file
    polarization:
        image polarization: HH, VV, HV, VH, CH, CV
    MLI_par:
        (output) MLI parameter file (example: yyyymmdd_pp.mli.par)
    MLI:
        (output) MLI data file in slant range geometry (example: yyyymmdd_pp.mli, enter - for none)
    GRD_par:
        (output) GRD parameter file (example: yyyymmdd_pp.grd.par, enter - for none)
    GRD:
        (output) GRD data file (example: yyyymmdd_pp.grd, enter - for none)
    rps:
        slant range pixel spacing (m) (enter - for default: calculated from ground-range parameters)
    radcal:
        radiometric calibration flag (enter - for default)
            * 0: beta0 (default)
            * 1: sigma0
    
    noise:
        noise levels flag (enter - for default)
            * 0: do not use noise levels (default)
            * 1: use noise levels
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_NovaSAR_GRD', GeoTIFF, XML, polarization, MLI_par, MLI, GRD_par, GRD, rps, radcal, noise]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_NovaSAR_SLC(GeoTIFF, XML, polarization, SLC_par, SLC='-', dtype='-', radcal='-', noise='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate SLC parameter and image files for NovaSAR SLC data
    | Copyright 2023, Gamma Remote Sensing, v1.6 3-Mar-2023 cm
    
    Parameters
    ----------
    GeoTIFF:
        (input) NovaSAR image data file in GeoTIFF format (enter - for none, \\*.tif)
    XML:
        (input) NovaSAR XML annotation file
    polarization:
        image polarization: HH, VV, HV, VH, CH, CV
    SLC_par:
        (output) ISP SLC parameter file (example: yyyymmdd_pp.slc.par)
    SLC:
        (output) SLC data file (enter - for none, example: yyyymmdd_pp.slc)
    dtype:
        output data type (enter - for default: same as input)
            * 0: FCOMPLEX
            * 1: SCOMPLEX
    
    radcal:
        radiometric calibration flag (enter - for default)
            * 0: beta0 (default)
            * 1: sigma0
    
    noise:
        noise levels flag (enter - for default)
            * 0: do not use noise levels (default)
            * 1: use noise levels
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_NovaSAR_SLC', GeoTIFF, XML, polarization, SLC_par, SLC, dtype, radcal, noise]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_NovaSAR_SRD(GeoTIFF, XML, polarization, MLI_par, MLI='-', radcal='-', noise='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate MLI image and parameter files for NovaSAR SRD data
    | Copyright 2023, Gamma Remote Sensing, v1.3 3-Mar-2023 cm
    
    Parameters
    ----------
    GeoTIFF:
        (input) NovaSAR image data file in GeoTIFF format (enter - for none, \\*.tif)
    XML:
        (input) NovaSAR XML annotation file
    polarization:
        image polarization: HH, VV, HV, VH, CH, CV
    MLI_par:
        (output) MLI parameter file (example: yyyymmdd_pp.mli.par)
    MLI:
        (output) MLI data file in slant range geometry (example: yyyymmdd_pp.mli, enter - for none)
    radcal:
        radiometric calibration flag (enter - for default)
            * 0: beta0 (default)
            * 1: sigma0
    
    noise:
        noise levels flag (enter - for default)
            * 0: do not use noise levels (default)
            * 1: use noise levels
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_NovaSAR_SRD', GeoTIFF, XML, polarization, MLI_par, MLI, radcal, noise]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_PRI(CEOS_SAR_leader, PRI_par, CEOS_DAT, PRI, logpath=None, outdir=None, shellscript=None):
    """
    | ISP parameter file generation for ERS PRI data from the PGS and VMP processors
    | Copyright 2020, Gamma Remote Sensing, v1.7 21-Sep-2020 clw/cm
    
    Parameters
    ----------
    CEOS_SAR_leader:
        (input) ERS CEOS SAR leader file for PRI product
    PRI_par:
        (output) ISP image parameter file (example: <yyyymmdd>.pri.par)
    CEOS_DAT:
        (input) CEOS data file (example: DAT_01.001)
    PRI:
        (output) PRI data with file and line headers removed (example: <yyyymmdd>.pri)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_PRI', CEOS_SAR_leader, PRI_par, CEOS_DAT, PRI]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_PRI_ESRIN_JERS(CEOS_SAR_leader, PRI_par, CEOS_DAT, PRI, logpath=None, outdir=None, shellscript=None):
    """
    | ISP GRD parameter file for ESRIN processed JERS PRI data
    | Copyright 2020, Gamma Remote Sensing, v1.9 21-Sep-2020 clw/uw/cm
    
    Parameters
    ----------
    CEOS_SAR_leader:
        (input) ERS CEOS SAR leader file for PRI product
    PRI_par:
        (output) ISP image parameter file (example: <yyyymmdd>.pri.par)
    CEOS_DAT:
        (input) CEOS data file (example: DAT_01.001)
    PRI:
        (output) PRI data with file and line headers removed (example: <yyyymmdd>.pri)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_PRI_ESRIN_JERS', CEOS_SAR_leader, PRI_par, CEOS_DAT, PRI]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_PulSAR(CEOS_SAR_leader, SLC_par, logpath=None, outdir=None, shellscript=None):
    """
    | ISP parameter file generation for ERS SLC data from the PULSAR SAR processor
    | Copyright 2020, Gamma Remote Sensing, v1.3 21-Sep-2020 clw/uw/cm
    
    Parameters
    ----------
    CEOS_SAR_leader:
        (input) ERS CEOS SAR leader file
    SLC_par:
        (output) ISP SLC parameter file (example <orbit>.slc.par)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_PulSAR', CEOS_SAR_leader, SLC_par]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_RCM_GRC(RCM_dir, polarization, radcal, noise, SLC_par='-', SLC='-', GRC_par='-', GRC='-', rps='-', noise_pwr='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate slant and ground range complex images and parameter files from a Radarsat Constellation GRC (Ground Range georeferenced Complex) product
    | Copyright 2024, Gamma Remote Sensing, v3.0 21-Oct-2024 cm
    
    Parameters
    ----------
    RCM_dir:
        (input) Radarsat Constellation main directory path (e.g.: RCM3_OK1001322_PK1001415_1_5M4_20160417_004803_VV_GRC)
    polarization:
        image polarization: HH, VV, HV, VH, CH, CV
    radcal:
        radiometric calibration flag (enter - for default)
            * 0: none (default)
            * 1: Beta Nought
            * 2: Sigma Nought
            * 3: Gamma
    
    noise:
        noise levels flag (enter - for default)
            * 0: do not use noise levels file (default)
            * 1: use noise levels file
            * NOTE: noise levels file can only be used for radiometrically calibrated data (radcal flag: 1, 2, or 3)
    
    SLC_par:
        (output) SLC parameter file (example: yyyymmdd_pp.slc.par, enter - for none)
    SLC:
        (output) SLC data file in slant range geometry (example: yyyymmdd_pp.slc, enter - for none)
    GRC_par:
        (output) GRC parameter file (example: yyyymmdd_pp.grc.par, enter - for none)
    GRC:
        (output) GRC data file (example: yyyymmdd_pp.grc, enter - for none)
    rps:
        slant range pixel spacing (m) (enter - for default: calculated from ground-range parameters)
    noise_pwr:
        (output) noise intensity for each SLC sample in slant range using data from noise levels file (enter - for none)
            * NOTE: when the noise_pwr file is specified, noise power correction will NOT be applied to the GRC / SLC data values
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_RCM_GRC', RCM_dir, polarization, radcal, noise, SLC_par, SLC, GRC_par, GRC, rps, noise_pwr]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_RCM_GRD(RCM_dir, polarization, radcal, noise, MLI_par='-', MLI='-', GRD_par='-', GRD='-', rps='-', noise_pwr='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate MLI and GRD images and parameter files from a Radarsat Constellation GRD (Ground Range georeferenced Detected) product
    | Copyright 2024, Gamma Remote Sensing, v2.9 21-Oct-2024 cm
    
    Parameters
    ----------
    RCM_dir:
        (input) Radarsat Constellation main directory path (e.g.: RCM1_OK1001327_PK1001418_1_3M28_20160417_013625_HH_GRD)
    polarization:
        image polarization: HH, VV, HV, VH, CH, CV
    radcal:
        radiometric calibration flag (enter - for default)
            * 0: none (default)
            * 1: Beta Nought
            * 2: Sigma Nought
            * 3: Gamma
    
    noise:
        noise levels flag (enter - for default)
            * 0: do not use noise levels file (default)
            * 1: use noise levels file
            * NOTE: noise levels file can only be used for radiometrically calibrated data (radcal flag: 1, 2, or 3)
    
    MLI_par:
        (output) MLI parameter file (example: yyyymmdd_pp.mli.par, enter - for none)
    MLI:
        (output) MLI data file in slant range geometry (example: yyyymmdd_pp.mli, enter - for none)
    GRD_par:
        (output) GRD parameter file (example: yyyymmdd_pp.grd.par, enter - for none)
    GRD:
        (output) GRD data file (example: yyyymmdd_pp.grd, enter - for none)
    rps:
        slant range pixel spacing (m) (enter - for default: calculated from ground-range parameters)
    noise_pwr:
        (output) noise intensity for each MLI sample in slant range using data from noise levels file (enter - for none)
            * NOTE: when the noise_pwr file is specified, noise power correction will NOT be applied to the GRD / MLI data values
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_RCM_GRD', RCM_dir, polarization, radcal, noise, MLI_par, MLI, GRD_par, GRD, rps, noise_pwr]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_RCM_MLC(RCM_dir, radcal, noise, root_name, logpath=None, outdir=None, shellscript=None):
    """
    | Generate parameter and image files for Radarsat Constellation MLC (Multi-Look Complex) data from GeoTIFF or NITF format
    | Copyright 2024, Gamma Remote Sensing, v1.4 21-Oct-2024 cm
    
    Parameters
    ----------
    RCM_dir:
        (input) Radarsat Constellation main directory path (e.g.: RCM2_OK1782060_PK1782073_2_SC30MCPC_20200504_105537_CH_CV_MLC)
    radcal:
        radiometric calibration flag (enter - for default)
            * 0: none (default)
            * 1: Beta Nought
            * 2: Sigma Nought
            * 3: Gamma
    
    noise:
        noise levels flag (enter - for default)
            * 0: do not use noise levels file (default)
            * 1: use noise levels file
            * NOTE: noise levels file can only be used for radiometrically calibrated data (radcal flag: 1, 2, or 3)
    
    root_name:
        (output) root name of the generated output files (example: yyyymmdd)
            * NOTE: the program will automatically complete the root_name and add extensions for each covariance matrix element
              for both data and parameter files, such as 20210927_CH.mlc, 20210927_CH.mlc.par, 20210927_XC.mlc, etc.
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_RCM_MLC', RCM_dir, radcal, noise, root_name]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_RCM_SLC(RCM_dir, polarization, radcal, noise, SLC_par, SLC, noise_pwr='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate SLC parameter and image files for Radarsat Constellation SLC data from GeoTIFF or NITF file
    | Copyright 2024, Gamma Remote Sensing, v2.7 21-Oct-2024 cm
    
    Parameters
    ----------
    RCM_dir:
        (input) Radarsat Constellation main directory path (e.g.: RCM2_OK1002260_PK1002436_3_SC50MB_20160417_002427_VH_VV_SLC)
    polarization:
        image polarization: HH, VV, HV, VH, CH, CV
    radcal:
        radiometric calibration flag (enter - for default)
            * 0: none (default)
            * 1: Beta Nought
            * 2: Sigma Nought
            * 3: Gamma
    
    noise:
        noise levels flag (enter - for default)
            * 0: do not use noise levels file (default)
            * 1: use noise levels file
            * NOTE: noise levels file can only be used for radiometrically calibrated data (radcal flag: 1, 2, or 3)
    
    SLC_par:
        (output) ISP SLC parameter file (example: yyyymmdd_pp.slc.par)
    SLC:
        (output) SLC data file (example: yyyymmdd_pp.slc)
    noise_pwr:
        (output) noise intensity for each SLC sample in slant range using data from noise levels file (enter - for none)
            * NOTE: when the noise_pwr file is specified, noise power correction will NOT be applied to the SLC data values
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_RCM_SLC', RCM_dir, polarization, radcal, noise, SLC_par, SLC, noise_pwr]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_RCM_SLC_ScanSAR(RCM_dir, polarization, radcal, noise_in, root_name, SLC_tab='-', beam='-', noise_out='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate SLC parameter and image files from Radarsat Constellation ScanSAR SLC data in GeoTIFF or NITF format
    | Copyright 2024, Gamma Remote Sensing, v3.3 21-Oct-2024 cm
    
    Parameters
    ----------
    RCM_dir:
        (input) Radarsat Constellation main directory path (e.g.: RCM2_OK1002260_PK1002436_3_SC50MB_20160417_002427_VH_VV_SLC)
    polarization:
        image polarization: HH, VV, HV, VH, CH, CV
    radcal:
        radiometric calibration flag (enter - for default)
            * 0: none (default)
            * 1: Beta Nought
            * 2: Sigma Nought
            * 3: Gamma
    
    noise_in:
        noise levels flag (enter - for default)
            * 0: do not use noise levels file (default)
            * 1: use noise levels file
            * NOTE: noise levels file can only be used for radiometrically calibrated data (radcal flag: 1, 2, or 3)
    
    root_name:
        (output) root name of the generated output files (example: yyyymmdd_pp)
            * NOTE: the program will automatically complete the root_name with beam numbers and extensions for the SLC, SLC_par, and TOPS_par files
    
    SLC_tab:
        (output) 3 column list of SLC, SLC_par, and TOPS_par files, with the beams sorted from near to far range (example: yyyymmdd_pp.SLC_tab)
    beam:
        number specifying the desired ScanSAR beam number (enter - for default: extract all beams)
            * NOTE: enter 0 to get the list of the available beams
    
    noise_out:
        output noise intensity for each SLC sample in slant range flag (enter - for default)
            * 0: do not write noise intensity files (default)
            * 1: write noise intensity files (file name(s) automatically defined)
            * NOTE: when noise intensity files are written, noise power correction will NOT be applied to the SLC data values
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_RCM_SLC_ScanSAR', RCM_dir, polarization, radcal, noise_in, root_name, SLC_tab, beam, noise_out]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_RISAT_GRD(CEOS_leader, BAND_META, GRD_par, CEOS_image, GRD='-', line_dir='-', pix_dir='-', cal_flg='-', KdB='-', logpath=None, outdir=None, shellscript=None):
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
        (output) ISP GRD parameter file (example: YYYYMMDD.grd.par)
    CEOS_image:
        (input) CEOS Ground-Range image file (example: dat_01.001)
    GRD:
        (output) Ground-Range data with file and line headers removed (enter - for none: example: YYYYMMDD.grd)
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
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_RISAT_GRD', CEOS_leader, BAND_META, GRD_par, CEOS_image, GRD, line_dir, pix_dir, cal_flg, KdB]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_RISAT_SLC(CEOS_leader, BAND_META, SLC_par, CEOS_image, SLC='-', line_dir='-', pix_dir='-', cal_flg='-', KdB='-', logpath=None, outdir=None, shellscript=None):
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
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_RISAT_SLC', CEOS_leader, BAND_META, SLC_par, CEOS_image, SLC, line_dir, pix_dir, cal_flg, KdB]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_RSAT2_SG(product_XML, lut_XML, GeoTIFF, polarization, MLI_par='-', MLI='-', GRD_par='-', GRD='-', rps='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate MLI and GRD images and parameter files from Radarsat 2 SGF/SGX/SCF data
    | Copyright 2023, Gamma Remote Sensing, v2.2 7-Jun-2023 awi/cw/cm
    
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
    MLI_par:
        (output) MLI parameter file (example: yyyymmdd_pp.mli.par, enter - for none)
    MLI:
        (output) MLI data file in slant range geometry (example: yyyymmdd_pp.mli, enter - for none)
    GRD_par:
        (output) GRD parameter file (example: yyyymmdd_pp.grd.par, enter - for none)
    GRD:
        (output) GRD data file (example: yyyymmdd_pp.grd, enter - for none)
    rps:
        slant range pixel spacing (m) (enter - for default: calculated from ground-range parameters)
            * NOTE: Ground range geometry is less accurate than slant range geometry and should be avoided
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_RSAT2_SG', product_XML, lut_XML, GeoTIFF, polarization, MLI_par, MLI, GRD_par, GRD, rps]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_RSAT2_SLC(product_XML, lut_XML, GeoTIFF, polarization, SLC_par, SLC, logpath=None, outdir=None, shellscript=None):
    """
    | Generate SLC parameter and image files for Radarsat 2 SLC data from GeoTIFF
    | Copyright 2023, Gamma Remote Sensing, v2.9 7-Jun-2023 awi/clw/cm
    
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
        (output) ISP SLC parameter file (example: yyyymmdd_pp.slc.par)
    SLC:
        (output) SLC data file (example: yyyymmdd_pp.slc)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_RSAT2_SLC', product_XML, lut_XML, GeoTIFF, polarization, SLC_par, SLC]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_RSAT_SCW(CEOS_leader, CEOS_trailer, CEOS_data, GRD_par, GRD, sc_dB='-', dt='-', logpath=None, outdir=None, shellscript=None):
    """
    | ISP parameter file for SCANSAR Wide Swath Data
    | Copyright 2020, Gamma Remote Sensing, v2.2 3-Sep-2020 clw/cm
    
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
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_RSAT_SCW', CEOS_leader, CEOS_trailer, CEOS_data, GRD_par, GRD, sc_dB, dt]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_RSAT_SGF(CEOS_leader, CEOS_data, GRD_par, GRD, sc_dB='-', dt='-', logpath=None, outdir=None, shellscript=None):
    """
    | ISP parameter file for RSI/Atlantis Radarsat SGF (ground range) and SCANSAR SCW16 data
    | Copyright 2020, Gamma Remote Sensing, v2.4 3-Sep-2020 clw/cm
    
    Parameters
    ----------
    CEOS_leader:
        (input) CEOS leader file (RSI SGF or SCW16 products, LEA_01.001)
    CEOS_data:
        (input) CEOS data file (RSI SGF or SCW16 products, DAT_01.001) 
    GRD_par:
        (output) ISP ground range image parameter file (example <orbit>.mli.par)
    GRD:
        (output) ISP ground range image (example <orbit>.grd.par) (enter -  for none, float)
    sc_dB:
        intensity scale factor in dB (enter - for default:   0.00)
    dt:
        azimuth image time offset (s) (enter - for default = 0.0)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_RSAT_SGF', CEOS_leader, CEOS_data, GRD_par, GRD, sc_dB, dt]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_RSAT_SLC(CEOS_leader, SLC_par, CEOS_data, SLC='-', sc_dB='-', dt='-', logpath=None, outdir=None, shellscript=None):
    """
    | ISP parameter file for RSI/Atlantis/ASF processed Radarsat SLC data
    | Copyright 2012, Gamma Remote Sensing, v4.0 5-Sep-2012 clw
    
    Parameters
    ----------
    CEOS_leader:
        (input) CEOS SAR leader file (example: lea_01.001)
    SLC_par:
        (output) ISP SLC parameter file (example: <date>.slc.par)
    CEOS_data:
        (input) CEOS data file (example: dat_01.001)
    SLC:
        (output) SLC data with file and line headers removed (example: <date>.slc)
    sc_dB:
        intensity scale factor in dB (enter - for default:  60.00)
    dt:
        azimuth image time offset (s) (enter - for default = 0.0)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_RSAT_SLC', CEOS_leader, SLC_par, CEOS_data, SLC, sc_dB, dt]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_RSI_ERS(CEOS_SAR_leader, SLC_par, logpath=None, outdir=None, shellscript=None):
    """
    | ISP parameter file for RSI processed ERS SLC data
    | Copyright 2020, Gamma Remote Sensing, v1.8 3-Sep-2020 clw/uw/cm
    
    Parameters
    ----------
    CEOS_SAR_leader:
        (input) ERS CEOS SAR leader file
    SLC_par:
        (output) ISP SLC parameter file (example <orbit>.slc.par)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_RSI_ERS', CEOS_SAR_leader, SLC_par]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_S1_GRD(GeoTIFF, annotation_XML, calibration_XML, noise_XML, MLI_par, MLI, GRD_par='-', GRD='-', eflg='-', rps='-', noise_pwr='-', edge_flag='-', loff='-', nl='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate MLI and GRD images and parameter files from a Sentinel-1 GRD product
    | Copyright 2023, Gamma Remote Sensing, v4.8 27-Apr-2023 awi/clw/ts/cm
    
    Parameters
    ----------
    GeoTIFF:
        (input) image data file in GeoTIFF format (enter - for none, \\*.tiff)
    annotation_XML:
        (input) Sentinel-1 L1 XML annotation file
    calibration_XML:
        (input) Sentinel-1 L1 radiometric calibration XML file (enter - for no radiometric calibration)
    noise_XML:
        (input) Sentinel-1 L1 noise XML file (enter - to not subtract thermal noise power level)
    MLI_par:
        (output) MLI parameter file (example: yyyymmdd_pp.mli.par)
    MLI:
        (output) MLI data file in slant range geometry (example: yyyymmdd_pp.mli, enter - for none)
    GRD_par:
        (output) GRD parameter file (example: yyyymmdd_pp.grd.par, enter - for none)
    GRD:
        (output) GRD data file (example: yyyymmdd_pp.grd, enter - for none)
    eflg:
        GR-SR grid extrapolation flag (enter - for default)
            * 0: no extrapolation of the GR-SR grid beyond the grid boundaries
            * 1: permit extrapolation of the GR-SR grid to cover the entire image (default)
            * NOTE: extrapolation of the GR-SR grid may introduce geocoding errors
    
    rps:
        slant range pixel spacing (m) (enter - for default: calculated from ground-range parameters)
    noise_pwr:
        noise intensity for each MLI sample in slant range using data from noise_XML (enter - for none)
            * NOTE: when the noise_pwr file is specified, noise power correction will NOT be applied to the MLI data values
    
    edge_flag:
        edge cleaning flag (enter - for default)
            * 0: do not clean edges (default for Sentinel-1 IPF version >= 2.90)
            * 1: basic method
            * 2: elaborate method based on Canny edge detection (default for Sentinel-1 IPF version < 2.90)
            * 3: force basic method when Sentinel-1 IPF version >= 2.90
            * 4: force elaborate method based on Canny edge detection when Sentinel-1 IPF version >= 2.90
            * NOTE: options 1 and 2 are changed to 0 when Sentinel-1 IPF version >= 2.90
    
    loff:
        offset to starting line of the input segment (enter - for default: 0)
    nl:
        number of lines to read from the file beginning at loff (enter - for default: to end of file)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_S1_GRD', GeoTIFF, annotation_XML, calibration_XML, noise_XML, MLI_par, MLI, GRD_par, GRD, eflg, rps, noise_pwr, edge_flag, loff, nl]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_S1_SLC(GeoTIFF, annotation_XML, calibration_XML, noise_XML, SLC_par, SLC, TOPS_par='-', dtype='-', sc_dB='-', noise_pwr='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate SLC parameter and image files for Sentinel-1 SLC data
    | Copyright 2025, Gamma Remote Sensing, v5.7 14-Apr-2025 awi/clw/cm
    
    Parameters
    ----------
    GeoTIFF:
        (input) image data file in \\*.tiff GeoTIFF format (enter - for default: none)
    annotation_XML:
        (input) Sentinel-1 L1 XML annotation file
    calibration_XML:
        (input) Sentinel-1 L1 radiometric calibration XML file to generate output as sigma0
            (enter - for default: return uncalibrated digital numbers)
    noise_XML:
        (input) Sentinel-1 L1 noise XML file (enter - for default: no subtraction of thermal noise power)
    SLC_par:
        (output) ISP SLC parameter file. Example: yyyymmdd_iw1_vv.slc.par
    SLC:
        (output) SLC data file (enter - for default: none). Example: yyyymmdd_iw1_vv.slc
    TOPS_par:
        (output) SLC burst annotation file; for TOPS and EW SLC data only (enter - for default: none). Example: yyyymmdd_iw1_vv.slc.tops_par
    dtype:
        output data type (enter - for default)
            * 0: FCOMPLEX (default)
            * 1: SCOMPLEX
    
    sc_dB:
        scale factor for FCOMPLEX -> SCOMPLEX, (enter - for default: HH,VV (dB): 60.0000,  VH,HV: 70.0000)
    noise_pwr:
        noise intensity for each SLC sample in slant range using data from noise_XML (enter - for none)
            * NOTE: when the noise_pwr file is specified, noise power will NOT be subtracted from the image data values
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_S1_SLC', GeoTIFF, annotation_XML, calibration_XML, noise_XML, SLC_par, SLC, TOPS_par, dtype, sc_dB, noise_pwr]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_SAOCOM_GRD(data, XML, MLI_par, MLI='-', GRD_par='-', GRD='-', rps='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate MLI parameter and image files for SAOCOM L1B Ground Range Detected Images
    | Copyright 2025, Gamma Remote Sensing, v1.0 13-Jan-2025 cm
    
    Parameters
    ----------
    data:
        (input) SAOCOM image data file in binary format (enter - for none, e.g. di--acqId0000729082-a-tw--2411281122-hh-m)
    XML:
        (input) SAOCOM XML annotation file (e.g. di--acqId0000729082-a-tw--2411281122-hh-m.xml)
    MLI_par:
        (output) MLI parameter file (example: yyyymmdd_pp.mli.par)
    MLI:
        (output) MLI data file (FCOMPLEX, enter - for none, example: yyyymmdd_pp.mli)
    GRD_par:
        (output) GRD parameter file (example: yyyymmdd_pp.grd.par, enter - for none)
    GRD:
        (output) GRD data file (example: yyyymmdd_pp.grd, enter - for none)
    rps:
        slant range pixel spacing (m) (enter - for default: calculated from ground-range parameters)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_SAOCOM_GRD', data, XML, MLI_par, MLI, GRD_par, GRD, rps]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_SAOCOM_SLC(data, XML, SLC_par, SLC='-', TOPS_par='-', RSLC_par='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate SLC parameter and image files for SAOCOM stripmap and TOPS SLC data
    | Copyright 2023, Gamma Remote Sensing, v1.6 21-Mar-2023 cm
    
    Parameters
    ----------
    data:
        (input) SAOCOM image data file in binary format (enter - for none, e.g. slc-acqId0000089010-a-tna-0000000000-s3qp-hh)
    XML:
        (input) SAOCOM XML annotation file (e.g. slc-acqId0000089010-a-tna-0000000000-s3qp-hh.xml)
    SLC_par:
        (output) SLC parameter file (example: yyyymmdd_s3_pp.slc.par)
    SLC:
        (output) SLC data file (FCOMPLEX, enter - for none, example: yyyymmdd_s3_pp.slc)
    TOPS_par:
        (output) SLC burst annotation file, TOPS data only (enter - for none, example: yyyymmdd_s3_vv.slc.tops_par)
    RSLC_par:
        (input) reference SLC parameter file to keep consistent range pixel spacing (example: yyyymmdd_s1_pp.slc.par)
            * NOTE: SAOCOM geocoded data in GeoTIFF format (GEC and GTC / level 1C and 1D data) can be read using par_SAOCOM_geo program
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_SAOCOM_SLC', data, XML, SLC_par, SLC, TOPS_par, RSLC_par]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_SICD_SLC(NITF, radcal, noise, SLC_par, SLC='-', XML='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate SLC parameter and image files for SICD SLC data
    | Copyright 2025, Gamma Remote Sensing, v2.0 28-Apr-2025 cm
    
    Parameters
    ----------
    NITF:
        (input) Sensor Independent Complex Data (SICD) file in NITF 2.1 container file (e.g.: CAPELLA_C03_SM_SICD_HH_20210512034455_20210512034459.ntf)
    radcal:
        radiometric calibration flag (enter - for default)
            * 0: none
            * 1: beta0 (default)
            * 2: sigma0
            * 3: gamma0
            * 4: RCS (target radar cross section in m^2)
    
    noise:
        noise levels flag (enter - for default)
            * 0: do not use noise levels (default)
            * 1: use noise levels
    
    SLC_par:
        (output) ISP SLC parameter file (example: yyyymmdd.slc.par)
    SLC:
        (output) SLC data file (enter - for none, example: yyyymmdd.slc)
    XML:
        (output) XML metadata file (enter - for none, example: CAPELLA_C03_SM_SICD_HH_20210512034455_20210512034459.xml)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_SICD_SLC', NITF, radcal, noise, SLC_par, SLC, XML]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_SIRC(CEOS_leader, SLC_par, UTC_MET='-', logpath=None, outdir=None, shellscript=None):
    """
    | ISP SLC parameter file from SIR-C CEOS leader file
    | Copyright 2025, Gamma Remote Sensing, v2.7 28-May-2025 clw/uw
    
    Parameters
    ----------
    CEOS_leader:
        (input) JPL SIR-C CEOS leader file
    SLC_par:
        (output) ISP SLC parameter file
    UTC_MET:
        time reference for state vectors: MET (Mission Elapsed Time) or UTC (enter - for default: UTC)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_SIRC', CEOS_leader, SLC_par, UTC_MET]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_STRIX(CEOS_leader, SLC_par, CEOS_data, SLC='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate SLC parameter and image files for Synspective StriX SLC data
    | Copyright 2023, Gamma Remote Sensing, v1.5 9-May-2023 awi/cm
    
    Parameters
    ----------
    CEOS_leader:
        (input) CEOS leader file for STRIX-alpha SLC data (LED-STRIXA...)
    SLC_par:
        (output) ISP image parameter file (example: yyyymmdd.slc.par)
    CEOS_data:
        (input) STRIX-alpha CEOS format SLC (IMG-pp-STRIXA...)
    SLC:
        (output) reformatted STRIX SLC (example: yyyymmdd.slc, enter - for none)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_STRIX', CEOS_leader, SLC_par, CEOS_data, SLC]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_SV2_SLC(GeoTIFF, annotation_XML, SLC_par, SLC='-', dtype='-', radcal='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate SLC parameter file and SLC image from a SuperView Neo-2 / SuperView-2 / Gaojing-2 data set
    | Copyright 2025, Gamma Remote Sensing, v1.3 12-May-2025 awi/cm
    
    Parameters
    ----------
    GeoTIFF:
        (input) image data file in GeoTIFF format (enter - for none, \\*.tiff)
    annotation_XML:
        (input) SV-2 product annotation XML file (\\*.meta.xml)
    SLC_par:
        (output) ISP SLC parameter file (example: yyyymmdd.slc.par)
    SLC:
        (output) SLC data file, example: yyyymmdd.slc (enter - for none, SLC output will not be produced)
    dtype:
        output data type (enter - for default)
            * 0: FCOMPLEX (default)
            * 1: SCOMPLEX
    
    radcal:
        output radiometric calibration flag (enter - for default)
            * 0: beta0
            * 1: sigma0 (default)
            * 2: gamma0
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_SV2_SLC', GeoTIFF, annotation_XML, SLC_par, SLC, dtype, radcal]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_SWOT_SLC(NETCDF, trunk, DEM='-', DEM_par='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate SLC parameter and image files for SWOT level 1B KaRIn SLC data
    | Copyright 2024, Gamma Remote Sensing, v1.2 30-Oct-2024 cm
    
    Parameters
    ----------
    NETCDF:
        (input) SWOT level 1B KaRIn SLC data file in NETCDF format (``SWOT_L1B_..._PIC0_01.nc``)
    trunk:
        (output) file name trunk used for output filenames
            (example: yyyymmdd -> yyyymmdd_L_minus_y.slc yyyymmdd_L_minus_y.slc.par)
    DEM:
        (output) DEM file in SCH coordinates (enter - for none)
    DEM_par:
        (output) DEM parameter file in SCH coordinates (enter - for none)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_SWOT_SLC', NETCDF, trunk, DEM, DEM_par]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_TX_GRD(annotation_XML, GeoTIFF, GRD_par, GRD='-', pol='-', MLI_par='-', MLI='-', rps='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate ground range image and image parameter file for Terrasar-X MGD data in GeoTIFF format
    | Copyright 2023, Gamma Remote Sensing, v1.5 8-May-2023 awi/clw/cm
    
    Parameters
    ----------
    annotation_XML:
        (input) Terrasar-X product annotation XML file
    GeoTIFF:
        (input) image data file in GeoTIFF format
            * NOTE: make sure the data set contains the selected polarization
    
    GRD_par:
        (output) ISP ground range image parameter file (example: yyyymmdd.grd.par, enter - for none)
    GRD:
        (output) calibrated ground range data file (example: yyyymmdd.grd, enter - for none)
    pol:
        polarization: HH, HV, VH, VV (enter - for default: first polarization found in the annotation_XML)
    MLI_par:
        (output) MLI parameter file (example: yyyymmdd.mli.par, enter - for none)
    MLI:
        (output) MLI data file in slant range geometry (example: yyyymmdd.mli, enter - for none)
    rps:
        slant range pixel spacing (m) (enter - for default: calculated from ground-range parameters)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_TX_GRD', annotation_XML, GeoTIFF, GRD_par, GRD, pol, MLI_par, MLI, rps]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_TX_ScanSAR(annotation_XML, swath, SLC_par, SLC, TOPS_par, bwflg='-', dtype='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate SLC, SLC_par and TOPS_par from a Terrasar-X ScanSAR data set
    | Copyright 2023, Gamma Remote Sensing, v2.4 18-Apr-2023 clw/cm/awi
    
    Parameters
    ----------
    annotation_XML:
        (input) TerraSAR-X ScanSAR product annotation XML file including path
            * NOTE: The path to the image products is determined from the path to the XML annotation
    
    swath:
        number specifying the desired ScanSAR swath (1 -> maximum number of swaths (4 or 6))
            * NOTE: The image product name is specified in the XML file
    
    SLC_par:
        (output) ISP SLC parameter file (example: yyyymmdd.slc.par)
    SLC:
        (output) SLC ScanSAR data file, example: yyyymmdd.slc
            (enter - for none, SLC output will not be produced)
    TOPS_par:
        (output) SLC ScanSAR burst annotation file (example: yyyymmdd_s1.slc.tops_par
    bwflg:
        burst window flag (enter - for default)
            * 0: use first and last annotation line values specified in the annotation_XML
            * 1: extend first and last valid line to include all data lines (default)
    
    dtype:
        output data type (enter - for default)
            * 0: same as input (default)
            * 1: FCOMPLEX
            * NOTE: While TSX ScanSAR data are not acquired in TOPS mode, the same data structure can be used for burst annotation
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_TX_ScanSAR', annotation_XML, swath, SLC_par, SLC, TOPS_par, bwflg, dtype]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_TX_SLC(annotation_XML, COSAR, SLC_par, SLC, pol='-', dtype='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate SLC parameter file and SLC image from a Terrasar-X SSC data set
    | Copyright 2023, Gamma Remote Sensing, v2.5 6-Mar-2023 awi/clw/cm
    
    Parameters
    ----------
    annotation_XML:
        (input) TerraSAR-X product annotation XML file
    COSAR:
        (input) COSAR SSC stripmap or spotlight mode SLC data file
    SLC_par:
        (output) ISP SLC parameter file (example: yyyymmdd.slc.par)
    SLC:
        (output) SLC data file, example: yyyymmdd.slc (enter - for none, SLC output will not be produced)
    pol:
        polarization HH, HV, VH, VV (enter - for default: first polarization found in the annotation_XML)
    dtype:
        output data type (enter - for default)
            * 0: same as input (default)
            * 1: FCOMPLEX
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_TX_SLC', annotation_XML, COSAR, SLC_par, SLC, pol, dtype]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def par_UAVSAR_SLC(ann, SLC_MLC_in, SLC_MLI_par, SLC_MLI_out='-', image_type='-', image_format='-', DOP='-', logpath=None, outdir=None, shellscript=None):
    """
    | ISP image parameter file from UAVSAR annotation file (ann) for SLC and MLC products
    | Copyright 2025, Gamma Remote Sensing, v2.0 31-Mar-2025 clw/cm
    
    Parameters
    ----------
    ann:
        (input) UAVSAR annotation file (\\*ann.txt or \\*.ann)
    SLC_MLC_in:
        (input) UAVSAR binary data file (required for annotation file version 1.2) (enter - for none)
    SLC_MLI_par:
        (output) ISP image parameter file
    SLC_MLI_out:
        (output) SLC data file (enter - for none)
    image_type:
        image type flag (enter - for default)
            * 0: SLC (slc) in slant range coordinates (default)
            * 1: MLC (mlc) in slant range coordinates
              HHHH\\*, VVVV\\*, HVHV\\* are FLOAT format
              HHHV\\*, HHVV\\*, HVVV\\* are FCOMPLEX format
    image_format:
        image data format flag (enter - for default)
            * 0: FCOMPLEX (pairs of 4-byte float (re,im)) (default)
            * 2: FLOAT  (4-bytes/value)
    
    DOP:
        (input) UAVSAR Doppler look-up table (if not zero-Doppler) (enter - for none)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/par_UAVSAR_SLC', ann, SLC_MLC_in, SLC_MLI_par, SLC_MLI_out, image_type, image_format, DOP]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def ph_slope_base(int_in, SLC_par, OFF_par, base, int_out, int_type='-', inverse='-', logpath=None, outdir=None, shellscript=None):
    """
    | Subtract/add interferogram flat-Earth phase trend as estimated from initial baseline
    | Copyright 2023, Gamma Remote Sensing, v4.5 19-Apr-2023 clw
    
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
        interferogram type (enter - for default)
            * 0: unwrapped phase
            * 1: complex interferogram (default)
    
    inverse:
        subtract/add inversion flag (enter - for default)
            * 0: subtract phase ramp (default)
            * 1: add phase ramp
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/ph_slope_base', int_in, SLC_par, OFF_par, base, int_out, int_type, inverse]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def phase_slope(interf, slopes, width, win_sz='-', thres='-', xmin='-', xmax='-', ymin='-', ymax='-', logpath=None, outdir=None, shellscript=None):
    """
    | Calculate interferogram phase slopes in range and azimuth
    | Copyright 2023, Gamma Remote Sensing, v1.4 18-Apr-2023 clw/uw
    
    Parameters
    ----------
    interf:
        (input) interferogram (fcomplex)
    slopes:
        (output) range and azimuth phase slopes (fcomplex)
    width:
        number of samples/row
    win_sz:
        size of region used for slopes determination (enter - for default: 5)
    thres:
        correlation threshold for accepting slope estimates 0.0 -> 1.0 (enter - for default: .4)
    xmin:
        starting range pixel offset (enter - for default: 0)
    xmax:
        last range pixel offset (enter - for default: width-1)
    ymin:
        starting azimuth row offset (enter - for default: 0)
    ymax:
        last azimuth row offset (enter - for default: nlines-1)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/phase_slope', interf, slopes, width, win_sz, thres, xmin, xmax, ymin, ymax]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def PRC_vec(SLC_par, PRC, nstate='-', logpath=None, outdir=None, shellscript=None):
    """
    | State vectors from ERS PRC orbit data for ISP processing clw/uw
    | Copyright 2023, Gamma Remote Sensing, v1.9 11-Oct-2023 clw
    
    Parameters
    ----------
    SLC_par:
        (input/output) ISP SLC/MLI image parameter file
    PRC:
        (input) PRC state vector file
    nstate:
        number of state vectors (enter - for default: 5, maximum: 1024)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/PRC_vec', SLC_par, PRC, nstate]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def ptarg_cal_MLI(MLI_par, MLI, r_samp, az_samp, psigma, c_r_samp, c_az_samp, ptr_image, r_plot, az_plot, pcal, osf='-', win='-', pltflg='-', psz='-', csz='-', theta_inc='-', logpath=None, outdir=None, shellscript=None):
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
        radar cross-section of the calibration target in m\\*\\*2
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
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/ptarg_cal_MLI', MLI_par, MLI, r_samp, az_samp, psigma, c_r_samp, c_az_samp, ptr_image, r_plot, az_plot, pcal, osf, win, pltflg, psz, csz, theta_inc]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def ptarg_cal_SLC(SLC_par, SLC, r_samp, az_samp, psigma, c_r_samp, c_az_samp, ptr_image, r_plot, az_plot, pcal, osf='-', win='-', pltflg='-', psz='-', csz='-', c_image='-', logpath=None, outdir=None, shellscript=None):
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
        radar cross-section of the calibration target in m\\*\\*2
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
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/ptarg_cal_SLC', SLC_par, SLC, r_samp, az_samp, psigma, c_r_samp, c_az_samp, ptr_image, r_plot, az_plot, pcal, osf, win, pltflg, psz, csz, c_image]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def ptarg_SLC(SLC_par, SLC, r_samp, az_samp, ptr_image, r_plot, az_plot, ptr_par='-', osf='-', win='-', pltflg='-', logpath=None, outdir=None, shellscript=None):
    """
    | Point target response analysis and interpolation for SLC images
    | Copyright 2024, Gamma Remote Sensing, v2.0 4-Oct-2024 clw
    
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
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/ptarg_SLC', SLC_par, SLC, r_samp, az_samp, ptr_image, r_plot, az_plot, ptr_par, osf, win, pltflg]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def radcal_MLI(MLI, MLI_par, OFF_par, CMLI, antenna='-', rloss_flag='-', ant_flag='-', refarea_flag='-', sc_dB='-', K_dB='-', pix_area='-', logpath=None, outdir=None, shellscript=None):
    """
    | Radiometric calibration for multi-look intensity (MLI) data
    | Copyright 2023, Gamma Remote Sensing, v2.4 6-Jul-2023 uw/clw/of
    
    Parameters
    ----------
    MLI:
        (input) MLI image (FLOAT)
    MLI_par:
        (input) SLC parameter file of input MLI image
    OFF_par:
        (input) ISP offset/interferogram parameter file (enter - for images in MLI geometry)
    CMLI:
        (output) radiometrically calibrated output MLI (FLOAT)
    antenna:
        (input) 1-way antenna gain pattern file (enter - for none)
    rloss_flag:
        range spreading loss correction (enter - for default)
            * 0: no correction (default)
            * 1: apply r^3 correction  (all modes except ASAR APS)
            * 2: apply r^4 correction (used only for ASAR APS mode)
            * -1: undo r^3 correction
            * -2: undo r^4 correction
    
    ant_flag:
        antenna pattern correction (enter - for default)
            * 0: no correction (default)
            * 1: apply antenna pattern correction
            * -1: undo antenna pattern correction
    
    refarea_flag:
        reference pixel area correction (enter - for default)
            * 0: no pixel area correction (default)
            * 1: calculate sigma0, scale area by sin(inc_ang)/sin(ref_inc_ang)
            * 2: calculate gamma0, scale area by sin(inc_ang)/(cos(inc_ang)\\*sin(ref_inc_ang)
            * -1: undo sigma0 area scaling factor
            * -2: undo gamma0 area scaling factor
    
    sc_dB:
        scale factor in dB (enter - for default: 0.0)
    K_dB:
        calibration factor in dB (enter - for default: value from MLI_par)
    pix_area:
        (output) ellipsoid-based ground range sigma0 or gamma0 pixel reference area (FLOAT) (enter - for none)
            refarea_flag 1 or -1: sigma0 ref. area
            refarea_flag 2 or -2: gamma0 ref. area
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/radcal_MLI', MLI, MLI_par, OFF_par, CMLI, antenna, rloss_flag, ant_flag, refarea_flag, sc_dB, K_dB, pix_area]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def radcal_PRI(PRI, PRI_par, GRD, GRD_par, K_dB='-', inc_ref='-', roff='-', nr='-', loff='-', nl='-', logpath=None, outdir=None, shellscript=None):
    """
    | Convert ESA processed short integer format PRI to radiometrically calibrated GRD image (float)
    | Copyright 2023, Gamma Remote Sensing, v1.7 19-Apr-2023 uw/clw
    
    Parameters
    ----------
    PRI:
        (input) PRI ground-range image (short integer, sqrt(backscat. intensity)
    PRI_par:
        (input) SLC parameter file of input PRI ground-range image (yyyymmdd.pri.par)
    GRD:
        (output) calibrated ground-range image (float, backscat. intensity)
    GRD_par:
        (output) ISP image parameter file of output calibrated ground-range image (yyyymmdd.grd.par)
    K_dB:
        calibration factor in decibels (enter - for default: 59.75 dB)
            ERS1 (D-Paf,ESRIN): 58.24 dB, ERS2 (D-Paf,ESRIN,I-Paf,UK-Paf after 1997): 59.75 dB
            ENVISAT ASAR: 55.0 dB (all modes)
            for details see product specifications and ESA publications.
    inc_ref:
        reference incidence angle in deg. (enter - for default: 23.0 deg.)
            ENVISAT ASAR: 90.0 deg. (all modes)
    roff:
        offset to starting range sample (enter - for default: 0)
    nr:
        number of range samples (enter - for default: to end of line)
    loff:
        offset to starting line (enter - for default: 0, 1 header line in the input file is assumed for ERS)
    nl:
        number of lines to copy (enter - for default: to end of file)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/radcal_PRI', PRI, PRI_par, GRD, GRD_par, K_dB, inc_ref, roff, nr, loff, nl]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def radcal_pwr_stat(SLC_tab, SLC_tab_cal, plist, MSR_cal, PWR_cal, roff='-', loff='-', nr='-', nl='-', plist_out='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate calibrated SLC image files using point targets determined from the Mean/Sigma Ratio and Intensity
    | Copyright 2022, Gamma Remote Sensing, v1.5 8-Nov-2022 clw/uw/cm
    
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
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/radcal_pwr_stat', SLC_tab, SLC_tab_cal, plist, MSR_cal, PWR_cal, roff, loff, nr, nl, plist_out]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def radcal_SLC(SLC, SLC_par, CSLC, CSLC_par, fcase='-', antenna='-', rloss_flag='-', ant_flag='-', refarea_flag='-', sc_dB='-', K_dB='-', pix_area='-', logpath=None, outdir=None, shellscript=None):
    """
    | Radiometric calibration of SLC data
    | Copyright 2023, Gamma Remote Sensing, v2.8 6-Jul-2023 uw/clw/of
    
    Parameters
    ----------
    SLC:
        (input) SLC (FCOMPLEX or SCOMPLEX)
    SLC_par:
        (input) SLC parameter file of input SLC
    CSLC:
        (output) radiometrically calibrated SLC (FCOMPLEX or SCOMPLEX)
    CSLC_par:
        (output) SLC parameter file of output calibrated SLC
    fcase:
        format case (enter - for default)
            * 1: FCOMPLEX --> FCOMPLEX (pairs of FLOAT) (default)
            * 2: FCOMPLEX --> SCOMPLEX (pairs of SHORT INTEGER)
            * 3: SCOMPLEX --> FCOMPLEX
            * 4: SCOMPLEX --> SCOMPLEX
    
    antenna:
        1-way antenna gain pattern file (enter - for none)
    rloss_flag:
        range spreading loss correction (enter - for default)
            * 0: no correction (default)
            * 1: apply r^3 correction  (all modes except ASAR APS)
            * 2: apply r^4 correction (used only for ASAR APS mode)
            * -1: undo r^3 correction
            * -2: undo r^4 correction
    
    ant_flag:
        antenna pattern correction (enter - for default)
            * 0: no correction (default)
            * 1: apply antenna pattern correction
            * -1: undo antenna pattern correction
    
    refarea_flag:
        reference pixel area correction (enter - for default)
            * 0: no pixel area correction (default)
            * 1: calculate sigma0, scale area by sin(inc_ang)/sin(ref_inc_ang)
            * 2: calculate gamma0, scale area by sin(inc_ang)/(cos(inc_ang)\\*sin(ref_inc_ang)
            * -1: undo sigma0 area scaling factor
            * -2: undo gamma0 area scaling factor
    
    sc_dB:
        scale factor in dB (enter - for default: 0.0)
    K_dB:
        calibration factor in dB (enter - for default: value from SLC_par)
    pix_area:
        (output) ellipsoid-based ground range sigma0 or gamma0 pixel reference area (FLOAT) (enter - for none)
            refarea_flag 1 or -1: sigma0 ref. area
            refarea_flag 2 or -2: gamma0 ref. area
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/radcal_SLC', SLC, SLC_par, CSLC, CSLC_par, fcase, antenna, rloss_flag, ant_flag, refarea_flag, sc_dB, K_dB, pix_area]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def rascc_mask(cc, pwr, width, start_cc='-', start_pwr='-', nlines='-', pixavr='-', pixavaz='-', cc_thres='-', pwr_thres='-', cc_min='-', cc_max='-', scale='-', exp='-', LR='-', rasf='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate phase unwrapping validity mask using correlation and intensity
    | Copyright 2023, Gamma Remote Sensing, v2.2 19-Apr-2023 clw/uw
    
    Parameters
    ----------
    cc:
        (input) interferometric correlation image (FLOAT)
    pwr:
        (input) intensity image (FLOAT, enter - if not available)
    width:
        number of samples/row
    start_cc:
        starting line of coherence image (enter - for default: 1)
    start_pwr:
        starting line of intensity image (enter - for default: 1)
    nlines:
        number of lines to display (enter - or 0 for default: to end of file)
    pixavr:
        number of pixels to average in range (enter - for default: 1)
    pixavaz:
        number of pixels to average in azimuth (enter - for default: 1)
    cc_thres:
        coherence threshold for masking, pixels with cc < cc_thres are set to 0 (enter - for default: 0.0)
    pwr_thres:
        relative intensity threshold for masking, pixels with intensity < pwr_thres \\* average intensity are set to 0 (enter - for default: 0)
    cc_min:
        minimum coherence value used for color display (enter - for default: 0.1)
    cc_max:
        maximum coherence value used for color display (enter - for default: 0.9)
    scale:
        intensity display scale factor (enter - for default: 1.0)
    exp:
        intensity display exponent (enter - for default: 0.35)
    LR:
        image mirror flag (enter - for default)
            * 1: normal (default)
            * -1: mirror image
    
    rasf:
        (output) image filename, extension determines the format, enter - for default: \\*.tif
            \\*.bmp BMP format
            \\*.ras Sun raster format
            \\*.tif TIFF format
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/rascc_mask', cc, pwr, width, start_cc, start_pwr, nlines, pixavr, pixavaz, cc_thres, pwr_thres, cc_min, cc_max, scale, exp, LR, rasf]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def rascc_mask_thinning(ras_in, in_file, width, ras_out, nmax='-', thresholds='-', logpath=None, outdir=None, shellscript=None):
    """
    | Adaptive sampling reduction for phase unwrapping validity mask
    | Copyright 2023, Gamma Remote Sensing, v1.7 19-Apr-2023 uw/clw
    
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
        number of sampling reduction runs (enter - for default: 3)
    thresholds:
        a list of thresholds sorted from smallest to largest scale sampling reduction
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/rascc_mask_thinning', ras_in, in_file, width, ras_out, nmax, thresholds]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def res_map(hgt, gr, data, SLC_par, OFF_par, res_hgt, res_data, nr='-', naz='-', azps_res='-', loff='-', nlines='-', logpath=None, outdir=None, shellscript=None):
    """
    | Slant range to ground range transformation based on interferometric ground-range
    | Copyright 2023, Gamma Remote Sensing, v2.6 18-Apr-2023 clw/uw
    
    Parameters
    ----------
    hgt:
        (input) height file in slant range geometry
    gr:
        (input) ground range file in slant range geometry
    data:
        (input) data file in slant range geometry (float) (intensity \\*.pwr or correlation \\*.cc)
    SLC_par:
        (input) ISP parameter file of reference SLC
    OFF_par:
        (input) offset/interferogram processing parameters
    res_hgt:
        (output) resampled height file in ground range geometry
    res_data:
        (output) resampled data file in ground range geometry
    nr:
        number of range samples for L.S. estimate (enter - for default: 7, must be odd)
    naz:
        number of azimuth samples for L.S. extimate (enter - for default: 7, must be odd)
    azps_res:
        azimuth output map sample spacing in meters (enter - for default: azimuth spacing)
    loff:
        offset to starting line for height calculations (enter - for default: 0)
    nlines:
        number of lines to calculate (enter - for default: to end of file)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/res_map', hgt, gr, data, SLC_par, OFF_par, res_hgt, res_data, nr, naz, azps_res, loff, nlines]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def residue(int, flag, width, xmin='-', xmax='-', ymin='-', ymax='-', logpath=None, outdir=None, shellscript=None):
    """
    | Determine interferometric phase unwrapping residues
    | Copyright 2023, Gamma Remote Sensing, v2.8 18-Apr-2023 clw/uw
    
    Parameters
    ----------
    int:
        (input) interferogram (fcomplex)
    flag:
        (input) flag file (unsigned char)
    width:
        number of samples/row
    xmin:
        offset to starting range pixel (enter - for default: 0)
    xmax:
        offset last range pixel (enter - for default: width-1)
    ymin:
        offset to starting azimuth row (enter - for default: 0)
    ymax:
        offset to last azimuth row (enter - for default: nlines-1)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/residue', int, flag, width, xmin, xmax, ymin, ymax]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def residue_cc(int, flag, width, xmin='-', xmax='-', ymin='-', ymax='-', logpath=None, outdir=None, shellscript=None):
    """
    | Determine interferometric phase unwrapping residues considering low coherence regions
    | Copyright 2023, Gamma Remote Sensing, v2.8 18-Apr-2023 clw/uw/ts
    
    Parameters
    ----------
    int:
        (input) interferogram (fcomplex)
    flag:
        (input) flag file (unsigned char)
    width:
        number of samples/row
    xmin:
        offset to starting range pixel (enter - for default: 0)
    xmax:
        offset last range pixel (enter - for default: width-1)
    ymin:
        offset to starting azimuth row (enter - for default: 0)
    ymax:
        offset to last azimuth row (enter - for default: nlines-1)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/residue_cc', int, flag, width, xmin, xmax, ymin, ymax]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def RSAT2_vec(SLC_par, RSAT2_orb, nstate='-', logpath=None, outdir=None, shellscript=None):
    """
    | Extract Radarsat-2 state vectors from a definitive orbit file
    | Copyright 2022, Gamma Remote Sensing, v1.1 clw/cm 7-Nov-2022
    
    Parameters
    ----------
    SLC_par:
        (input) ISP image parameter file
    RSAT2_orb:
        Radarsat-2 definitive orbit data file available from MDA (orbit_number_def.orb)
    nstate:
        number of state vectors to extract (enter - for default: 9)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/RSAT2_vec', SLC_par, RSAT2_orb, nstate]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def S1_burstloc(annotation_XML, logpath=None, outdir=None, shellscript=None):
    """
    | Print Burst information found in the Sentinel-1 annotation file
    | Copyright 2025, Gamma Remote Sensing, v1.4 3-Feb-2025 awi/cm
    
    Parameters
    ----------
    annotation_XML:
        (input) Sentinel-1 L1 XML annotation file
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/S1_burstloc', annotation_XML]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def S1_ETAD_SLC(ETAD, SLC1_tab, SLC2_tab, OPOD='-', corr='-', phase='-', tropo='-', iono='-', tides='-', bistatic='-', Doppler='-', FM_rate='-', mode='-', order='-', logpath=None, outdir=None, shellscript=None):
    """
    | Read and apply Sentinel-1 Extended Timing Annotation Dataset (ETAD) to correct range and azimuth timings of Sentinel-1 SLC images
    | Copyright 2025, Gamma Remote Sensing, v1.1 25-Jun-2025 cm
    
    Parameters
    ----------
    ETAD:
        (input) ETAD directory (e.g. S1A_IW_ETA__AXDV_20240807T172347_20240807T172414_055110_06B719_202E.SAFE)
            ETAD can be downloaded from https://dataspace.copernicus.eu/
    SLC1_tab:
        (input) SLC_tab of Sentinel-1 TOPS or Stripmap SLC (e.g. 20240807.SLC_tab)
    SLC2_tab:
        (output) SLC_tab of Sentinel-1 TOPS or Stripmap SLC with ETAD correction (e.g. 20240807.ETAD.SLC_tab)
    OPOD:
        replace state vectors by precision orbit data (OPOD) provided with ETAD data (enter - for default)
            * 0: no
            * 1: yes (default)
    
    corr:
        apply following timing corrections (enter - for default)
            * 0: no correction
            * 1: all corrections (default)
            * 2: all corrections in range only
            * 3: all corrections in azimuth only
            * 4: select individual corrections (defined in subsequent options)
    
    phase:
        apply phase corrections corresponding to the selected timing corrections in range (enter - for default)
            * 0: no
            * 1: yes (default)
            * 2: yes, experimental mode (phase corrections written to file(s))
    
    tropo:
        apply corrections for tropospheric delay in range (enter - for default)
            * 0: no
            * 1: yes (default)
    
    iono:
        apply corrections for ionospheric delay in range (enter - for default)
            * 0: no
            * 1: yes (default)
    
    tides:
        apply corrections for solid Earth tides (enter - for default)
            * 0: no
            * 1: yes, in range and azimuth (default)
            * 2: range only
            * 3: azimuth only
    
    bistatic:
        apply corrections for bistatic azimuth shifts (enter - for default)
            * 0: no
            * 1: yes (default)
    
    Doppler:
        apply corrections for Doppler-induced range shifts (enter - for default)
            * 0: no
            * 1: yes (default)
    
    FM_rate:
        apply corrections for FM-rate mismatch azimuth shifts (enter - for default)
            * 0: no
            * 1: yes (default)
    
    mode:
        complex data interpolation mode (enter - for default)
            * 0: Lanczos (default)
            * 1: B-spline
    
    order:
        Lanczos interpolator order / B-spline degree 4 -> 9 (enter - for default: 5)
            NOTES: - SLC1_tab and SLC2_tab or their contents can be the same files (the files will be overwritten in that case)
            - if SLC2_tab doesn't exist, it will be automatically created with file names derived from SLC1_tab contents
            - SLC_tab line entries:
            - TOPS mode:      SLC   SLC_par   TOPS_par
            - Stripmap mode:  SLC   SLC_par
            - with [phase] = 1, phase corrections only use ionospheric delays and solid Earth tides in range direction
            - with [phase] = 2, phase corrections also include compensation for tropospheric path delays
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/S1_ETAD_SLC', ETAD, SLC1_tab, SLC2_tab, OPOD, corr, phase, tropo, iono, tides, bistatic, Doppler, FM_rate, mode, order]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def S1_OPOD_vec(SLC_par, OPOD, nstate='-', logpath=None, outdir=None, shellscript=None):
    """
    | Extract Sentinel-1 OPOD state vectors and copy into the ISP image parameter file
    | Copyright 2025, Gamma Remote Sensing, v1.8 23-Jan-2024 awi/clw/cm
    
    Parameters
    ----------
    SLC_par:
        (input/output) ISP SLC/MLI image parameter file
    OPOD:
        (input) Sentinel-1 OPOD orbit data file (AUX_POEORB or AUX_RESORB)
            orbit files can be downloaded from https://s1qc.asf.alaska.edu/ or https://dataspace.copernicus.eu/
    nstate:
        number of state vectors to extract (enter - for default: include 60 sec extention at the start and end of the SLC data)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/S1_OPOD_vec', SLC_par, OPOD, nstate]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def sbi_filt(SLC1, SLC1_par, SLC2R_par, SLCf, SLCf_par, SLCb, SLCb_par, norm_sq, iwflg='-', logpath=None, outdir=None, shellscript=None):
    """
    | Azimuth filtering of SLC data to support split-beam interferometry to measure azimuth offsets
    | Copyright 2023, Gamma Remote Sensing, v1.6 clw/cm 18-Apr-2023
    
    Parameters
    ----------
    SLC1:
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
        inverse weighting flag (enter - for default)
            * 0: no compensation for azimuth spectrum weighting
            * 1: compensate for the azimuth spectrum weighting (default)
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/sbi_filt', SLC1, SLC1_par, SLC2R_par, SLCf, SLCf_par, SLCb, SLCb_par, norm_sq, iwflg]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def sbi_offset(sbi_unw, SLCf_par, SLCb_par, OFF_par, az_offset, logpath=None, outdir=None, shellscript=None):
    """
    | Calculate azimuth offsets from unwrapped split-beam interferogram
    | Copyright 2022, Gamma Remote Sensing, v1.1 8-Nov-2022
    
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
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/sbi_offset', sbi_unw, SLCf_par, SLCb_par, OFF_par, az_offset]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def ScanSAR_burst_copy(SLC, SLC_par, TOPS_par, SLC_out, SLC_out_par, burst_num, drflg='-', SLC_par2='-', dtype='-', logpath=None, outdir=None, shellscript=None):
    """
    | Copy selected burst from Sentinel-1 TOPS SLC to a file
    | Copyright 2023, Gamma Remote Sensing, v2.1 18-Apr-2023 awi/clw/cm
    
    Parameters
    ----------
    SLC:
        (input) ScanSAR mode burst SLC
    SLC_par:
        (input) SLC parameter file for the ScanSAR burst scene
    TOPS_par:
        (input) burst parameter file for the ScanSAR burst SLC
    SLC_out:
        (output) SLC file containing a single burst
    SLC_out_par:
        (output) SLC parameter file for the single burst
    burst_num:
        burst number of selected burst (1 -> number of bursts in the SLC)
    drflg:
        deramp phase flag (enter - for default)
            * 0: no modification of the burst SLC phase (default)
            * 1: subtract TOPS mode Doppler phase ramp for Sentinel-1 (deramp)
    
    SLC_par2:
        (output) SLC parameter file for the single burst SLC with deramped phase (drflg: 1, enter - for none)
    dtype:
        output data type (enter - for default: same as input data):
            * 0: FCOMPLEX
            * 1: SCOMPLEX
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/ScanSAR_burst_copy', SLC, SLC_par, TOPS_par, SLC_out, SLC_out_par, burst_num, drflg, SLC_par2, dtype]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def ScanSAR_burst_corners(SLC_par, TOPS_par, KML='-', logpath=None, outdir=None, shellscript=None):
    """
    | Calculate corner geographic coordinates of ScanSAR burst data and generate a KML with burst rectangles
    | Copyright 2025, Gamma Remote Sensing, v1.5 10-Feb-2025 awi/rc/cw
    
    Parameters
    ----------
    SLC_par:
        (input) SLC parameter file for the ScanSAR burst data
    TOPS_par:
        (input) ScanSAR burst parameter file
    KML:
        (output) KML output file
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/ScanSAR_burst_corners', SLC_par, TOPS_par, KML]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def ScanSAR_burst_MLI(SLC_tab, MLI_tab, rlks, azlks, bflg='-', SLCR_tab='-', MLI_dir='-', scale='-', logpath=None, outdir=None, shellscript=None):
    """
    | Generate MLI burst data from ScanSAR burst SLC data (Sentinel-1, RCM, and TSX)
    | Copyright 2024, Gamma Remote Sensing v2.5 25-Jun-2024 clw/cm
    
    Parameters
    ----------
    SLC_tab:
        (input) 3 column list of ScanSAR SLC, swaths are listed in order from near to far range
            SLC_tab line entries:   SLC   SLC_par  TOPS_par
    MLI_tab:
        (output) 3 column list of MLI swaths listed in order from near to far range
            MLI_tab line entries:   MLI   MLI_par  TOPS_par
            * NOTE: if the MLI_tab does not yet exist, the file entries will be created with names derived from the SLC_tab entries
    
    rlks:
        number of range looks  (1...80)
    azlks:
        number of azimuth look (1...20)
    bflg:
        burst window calculation flag (enter - for default)
            * 0: use existing burst window parameters if they exist, otherwise calculate burst window parameters (default)
            * 1: calculate burst window parameters from burst parameters and the number of range and azimuth looks
    
    SLCR_tab:
        (input) 3 column list of the reference scene with swaths, listed in order from near to far range (enter - for none)
            SLCR_tab line entries:   SLC    SLC_par   TOPS_par
    MLI_dir:
        directory for output burst MLI data, ignored if the MLI_tab already exists (enter - for default: current directory)
    scale:
        scale factor for output MLI (enter - for default: calculate from calibration gain in SLC parameter file)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/ScanSAR_burst_MLI', SLC_tab, MLI_tab, rlks, azlks, bflg, SLCR_tab, MLI_dir, scale]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def ScanSAR_burst_overlap(SLC_tab, root_name, rlks, azlks, mode='-', bflg='-', SLCR_tab='-', dburst='-', bound='-', logpath=None, outdir=None, shellscript=None):
    """
    | Extract and mosaic overlapping parts of ScanSAR / TOPS burst data
    | Copyright 2023, Gamma Remote Sensing v1.8 18-Apr-2023 cm/clw/uw
    
    Parameters
    ----------
    SLC_tab:
        (input) 3 column list of SLC, SLC_par, Sentinel-1 TOPS_par sorted in the order IW1, IW2, IW3...
    root_name:
        (output) output data root name (example: yyyymmdd_pp_overlap)
    rlks:
        number of range looks used to determine burst window boundaries
    azlks:
        number of azimuth looks used to determine burst window boundaries
    mode:
        output mode (enter - for default)
            * 0: output data are mosaics, non-overlapping parts are set to 0 (default)
            * 1: output data are mosaics, non-overlapping parts are written
            * 2: output data are burst data containing only overlapping parts
            * 3: output data is a polygon file with polygons encompassing overlapping areas in the SLC mosaic
            * 4: output data is a polygon file with polygons encompassing overlapping areas in the MLI mosaic
    
    bflg:
        burst window calculation flag (enter - for default)
            * 0: use existing burst window parameters if they exist, otherwise calculate burst window parameters (default)
            * 1: recalculate burst window parameters from burst parameters and the number of range and azimuth looks
    
    SLCR_tab:
        (input) SLC_tab of the reference scene, 3 column list of SLC, SLC_par, TOPS_par sorted sorted in the order IW1, IW2, IW3 (enter - for none)
            * NOTE: When generating a mosaic of a resampled SLC, the SLC_tab of the reference scene is required
    
    dburst:
        delta burst number (1=overlap of subsequent bursts, enter - for default: 1)
    bound:
        boundary pixels in polygon (enter - for default: 0)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/ScanSAR_burst_overlap', SLC_tab, root_name, rlks, azlks, mode, bflg, SLCR_tab, dburst, bound]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def ScanSAR_burst_to_mosaic(DATA_tab, mosaic, MLI_par, mflg='-', data_tab_ref='-', min_ovr='-', max_ovr='-', logpath=None, outdir=None, shellscript=None):
    """
    | Calculate mosaic of multilook ScanSAR burst data (FLOAT or FCOMPLEX)
    | Copyright 2023, Gamma Remote Sensing v2.5 18-Apr-2023 clw/cm
    
    Parameters
    ----------
    DATA_tab:
        (input) 3 column list of swaths in ML_DATA burst geometry listed in the order from near to far range
            DATA_tab line entries:   DATA   MLI_par  TOPS_par
            * NOTE: The data type (FLOAT or FCOMPLEX) is specified in the MLI_par and the burst parameters (TOPS_par) must agree
    
    mosaic:
        (output) mosaic image from bursts in multi-look geometry
    MLI_par:
        (output) mosaic image parameter file
    mflg:
        mosaicking option flag (enter - for default)
            * 0: no overlap between bursts or image swaths (default)
            * 1: average data in the overlap between bursts and in the overlap between image swaths
            * 2: average data in the overlap between bursts but not in the overlap between image swaths
    
    data_tab_ref:
        (input) reference scene DATA_tab, 3 column list of DATA, MLI_par, TOPS_par listed in order from near to far range (enter - for none)
            * NOTE: When generating a mosaic produced using data from a resampled scene, the MLI_tab of the reference scene is required
    
    min_ovr:
        minimum number of overlapping bursts (using mflg = 1 or 2, enter - for default: 1)
    max_ovr:
        maximum number of overlapping bursts (using mflg = 1 or 2, enter - for default: unlimited)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/ScanSAR_burst_to_mosaic', DATA_tab, mosaic, MLI_par, mflg, data_tab_ref, min_ovr, max_ovr]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def ScanSAR_full_aperture_SLC(SLC1_tab, SLC2_tab, SLCR_tab='-', SLC2_dir='-', vmode='-', wflg='-', imode='-', order='-', n_ovr='-', logpath=None, outdir=None, shellscript=None):
    """
    | Calculate continuous SLC data from ScanSAR burst data (Sentinel-1, RCM, and TSX)
    | Copyright 2023, Gamma Remote Sensing v1.9 18-Apr-2023 clw/cm
    
    Parameters
    ----------
    SLC1_tab:
        (input) 3 column list of ScanSAR SLC swaths listed in order from near to far range
            SLC1_tab line entries:   SLC   SLC_par  TOPS_par
    SLC2_tab:
        (input/output) 3 column list of oversampled continuous SLC swaths listed in order from near to far range
            SLC2_tab line entries:   SLC   SLC_par
            * NOTE: if the SLC2_tab does not yet exist, the file entries will be created with names derived from the SLC1_tab entries
    
    SLCR_tab:
        (input) 3 column list of the reference scene with swaths, listed in order from near to far range (enter - for none)
            SLCR_tab line entries:   SLC    SLC_par   TOPS_par
    SLC2_dir:
        directory for output oversampled continuous SLC, ignored if the SLC2_tab already exists (enter - or . for the current directory)
    vmode:
        sample validity mode (enter - for default):
            * 0: all data in the burst are considered valid (default)
            * 1: interpolate samples between the valid data bounds of the burst
    
    wflg:
        burst window calculation flag (enter - for default):
            * 0: use existing burst window parameters if they exist, otherwise calculate burst window parameters (default)
            * 1: calculate burst window parameters from burst parameters and the number of range and azimuth looks
    
    imode:
        interpolation mode (enter - for default):
            * 0: Lanczos (default)
            * 1: B-spline
    
    order:
        Lanczos interpolator order / B-spline degree 4 -> 9 (enter - for default: 5)
            dtype     output data type, (enter - for default: same as input data):
            * 0: FCOMPLEX
            * 1: SCOMPLEX
    
    n_ovr:
        SLC oversampling factor, must be in the range 2 --> 32 (enter - for default: automatically calculated)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/ScanSAR_full_aperture_SLC', SLC1_tab, SLC2_tab, SLCR_tab, SLC2_dir, vmode, wflg, imode, order, n_ovr]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def ScanSAR_mosaic_to_burst(DATA, MLI_par, DATA_tab, logpath=None, outdir=None, shellscript=None):
    """
    | Resample image data in the MLI mosaic geometry to burst MLI geometry (FLOAT or FCOMPLEX)
    | Copyright 2023, Gamma Remote Sensing v1.5 3-Apr-2023 clw/cm
    
    Parameters
    ----------
    DATA:
        (input) data in mosaic geometry (FLOAT or FCOMPLEX data type)
    MLI_par:
        image parameter file in mosaic geometry
    DATA_tab:
        3 column list of the output data in burst geometry, swaths are in order from near to far range
            MLI_tab line entries:  DATA   MLI_par  TOPS_par
            
            * NOTE: 1.The burst MLI_par and TOPS_par files describing the output geometry must already exist
              2.The data type (FLOAT or FCOMPLEX) specified in the MLI_par and the burst parameters (TOPS_par) must agree
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/ScanSAR_mosaic_to_burst', DATA, MLI_par, DATA_tab]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def slant_range(SLC_par, slr, logpath=None, outdir=None, shellscript=None):
    """
    | Calculate slant range for every range sample
    | Copyright 2022, Gamma Remote Sensing v1.2 8-Nov-2022 cw
    
    Parameters
    ----------
    SLC_par:
        (input) SLC or MLI image parameter file
    slr:
        (output) slant range for every sample in the image (float)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/slant_range', SLC_par, slr]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def SLC_adf(SLC, ref_SLC, ref_SLC_par, SLC_filt, mode='-', alpha='-', nfft_r='-', nfft_az='-', r_step='-', az_step='-', mwin_r='-', mwin_az='-', logpath=None, outdir=None, shellscript=None):
    """
    | Adaptive filtering of SLC data based on the local PSD of a reference SLC image
    | Copyright 2023, Gamma Remote Sensing, v1.4 18-Apr-2023 clw/cm
    
    Parameters
    ----------
    SLC:
        (input) SLC to be filtered (FCOMPLEX or SCOMPLEX)
    ref_SLC:
        (input) reference SLC
    ref_SLC_par:
        (input) reference SLC parameter file
    SLC_filt:
        (output) output filtered SLC using the power spectrum of the reference SLC
    mode:
        SLC filtering mode (enter - for default):
            * 0: 1D range PSD filter
            * 1: 1D azimuth PSD filter
            * 2: 2D range PSD \\* azimuth PSD filter 
            * 3: 2D median-filtered PSD filtering (default)
    
    alpha:
        exponent to apply to PSD value (enter - for default: 0.30)
    nfft_r:
        range filter FFT window size, 2\\*\\*N, 16->1024, (enter - for default: 128)
    nfft_az:
        azimuth filter FFT window size, 2\\*\\*N, 16->1024, (enter - for default: 128)
    r_step:
        range processing step (enter - for default: nfft_r/4)
    az_step:
        azimuth processing step (enter - for default: nfft_az/4)
    mwin_r:
        range median window size for median PSD filtering (enter - for default: 5)
    mwin_az:
        azimuth median window size for median PSD filtering (enter - for default: 5)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/SLC_adf', SLC, ref_SLC, ref_SLC_par, SLC_filt, mode, alpha, nfft_r, nfft_az, r_step, az_step, mwin_r, mwin_az]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def SLC_cat(SLC1, SLC2, SLC1_par, SLC2_par, OFF_par, SLC3, SLC3_par, dopflg='-', iflg='-', phflg='-', gainflg='-', imode='-', order='-', logpath=None, outdir=None, shellscript=None):
    """
    | Concatenate a pair of SLC images with interpolation of the second scene
    | Copyright 2024, Gamma Remote Sensing, v2.8 18-Jul-2024 clw/cm
    
    Parameters
    ----------
    SLC1:
        (input) SLC1 image (FCOMPLEX or SCOMPLEX)
    SLC2:
        (input) SLC2 image to be appended to SLC1 (same type as SLC1)
    SLC1_par:
        (input) SLC1 ISP image parameter file
    SLC2_par:
        (input) SLC2 ISP image parameter file
    OFF_par:
        (input) ISP offset parameter file containing offset polynomials between SLC1 and SLC2
    SLC3:
        (output) concatenated SLC
    SLC3_par:
        (output) ISP image parameter file for concatenated image
    dopflg:
        Doppler flag (enter - for default)
            * 0: ignore Doppler centroid information, assume 0 Hz Doppler centroid
            * 1: use Doppler centroid information for interpolation (default)
    
    iflg:
        input data type flag (enter - for default)
            * 0: input data are SLC images, use data type specified in SLC_par files (SCOMPLEX or FCOMPLEX) (default)
            * 1: input scenes are interferograms, force FCOMPLEX data type
    
    phflg:
        phase offset correction flag (enter - for default)
            * 0: no phase offset correction for SLC2 (default)
            * 1: apply constant phase offset correction to SLC2
    
    gainflg:
        gain correction flag (enter - for default)
            * 0: no gain correction for SLC2 (default)
            * 1: apply gain correction to SLC2 using calibration gain values in parameter files
            * 2: apply gain correction to SLC2 using relative intensity of overlap areas
    
    imode:
        interpolation mode for SLC2 (enter - for default)
            * 0: Lanczos interpolation (default)
            * 1: B-spline interpolation
    
    order:
        Lanczos interpolator order / B-spline degree 4 -> 9 (enter - for default: 4)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/SLC_cat', SLC1, SLC2, SLC1_par, SLC2_par, OFF_par, SLC3, SLC3_par, dopflg, iflg, phflg, gainflg, imode, order]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def SLC_cat_ScanSAR(SLC_tab1, SLC_tab2, SLC_tab3, bin_flag='-', logpath=None, outdir=None, shellscript=None):
    """
    | Concatenate sequential ScanSAR burst SLC images
    | Copyright 2024, Gamma Remote Sensing v3.5 5-Mar-2024 clw/cm
    
    Parameters
    ----------
    SLC_tab1:
        (input) 3 column list of ScanSAR SLC, swaths are listed in order from near to far range (earlier time)
            SLC_tab line entries:   SLC   SLC_par  TOPS_par
    SLC_tab2:
        (input) 3 column list of ScanSAR SLC, swaths are listed in order from near to far range (later time)
            SLC_tab line entries:   SLC   SLC_par  TOPS_par
    SLC_tab3:
        (input) 3 column list of concatenated ScanSAR SLC, swaths are listed in order from near to far range
            SLC_tab line entries:   SLC   SLC_par  TOPS_par
    bin_flag:
        binary data flag (enter - for default)
            * 0: no binary data generated (concatenate parameter files only)
            * 1: binary data generated (default)
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/SLC_cat_ScanSAR', SLC_tab1, SLC_tab2, SLC_tab3, bin_flag]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def SLC_copy(SLC_in, SLC_par_in, SLC_out, SLC_par_out, fcase='-', sc='-', roff='-', nr='-', loff='-', nl='-', swap='-', header_lines='-', logpath=None, outdir=None, shellscript=None):
    """
    | Copy SLC with options for data format conversion, segment extraction, swap real and imaginary, swap near and far range,  and azimuth spectrum shift
    | Copyright 2023, Gamma Remote Sensing, v6.1 1-May-2023 uw/clw/cm/of
    
    Parameters
    ----------
    SLC_in:
        (input) SLC (FCOMPLEX or SCOMPLEX format)
    SLC_par_in:
        (input) ISP SLC parameter file for input SLC
    SLC_out:
        (output) selected SLC section (FCOMPLEX or SCOMPLEX format)
    SLC_par_out:
        (output) ISP SLC parameter file of output SLC
    fcase:
        data format conversion (enter - for default: output format = input format)
            * 1: FCOMPLEX --> FCOMPLEX (default sc = 1.0)
            * 2: FCOMPLEX --> SCOMPLEX (default sc = 10000.0)
            * 3: SCOMPLEX --> FCOMPLEX (default sc = 0.0001)
            * 4: SCOMPLEX --> SCOMPLEX (default sc = 1.0)
    
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
            * 3: shift the SLC azimuth spectrum by 1/2 the azimuth sample rate
    
    header_lines:
        number of input file header lines (enter - for default: 0)
            * NOTE: CEOS format SLC data have 1 header line
            * NOTE: file offset pointer size (bytes): 8
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/SLC_copy', SLC_in, SLC_par_in, SLC_out, SLC_par_out, fcase, sc, roff, nr, loff, nl, swap, header_lines]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def SLC_copy_ScanSAR(SLC1_tab, SLC2_tab, BURST_tab, dtype='-', SLC2_dir='-', logpath=None, outdir=None, shellscript=None):
    """
    | Burst selection and copy from ScanSAR burst data (FCOMPLEX, SCOMPLEX)
    | Copyright 2024, Gamma Remote Sensing v3.7 29-Feb-2024 clw/cm
    
    Parameters
    ----------
    SLC1_tab:
        (input) 3 column list of ScanSAR SLC1 swaths in order from near to far range
            SLC1_tab line entries:   SLC    SLC_par   TOPS_par
    SLC2_tab:
        (input/output) 3 column list of the burst data copied from the ScanSAR swaths listed in SLC1_tab, in order from near to far range
            SLC2_tab line entries:   SLC    SLC_par   TOPS_par
            
            * NOTE: If the SLC2_tab does not yet exist, the SLC2_tab will be created with file names derived from the SLC1_tab entries and the SLC2_dir
              The new file names will have _2 appended to the root file names of the entries in SLC1_tab
    BURST_tab:
        (input) 2 column list of the first and last burst to copy from each swath, one line for each swath
            BURST_tab line entries: first_burst  last_burst
            NOTES: 1. The first burst is 1, enter - to select last physical burst
            2. If first_burst <= 0, then blank bursts are generated at the start of the output swath
            3. If last_burst exceeds the number of bursts, then blank bursts are appended to the end of the output swath
    dtype:
        output data format for complex data (enter - for default: output data has the same format as input data):
            * 0: FCOMPLEX
            * 1: SCOMPLEX
    
    SLC2_dir:
        directory for ScanSAR burst data copied from SLC1 data, ignored if the SLC2_tab already exists (enter - for default: current directory)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/SLC_copy_ScanSAR', SLC1_tab, SLC2_tab, BURST_tab, dtype, SLC2_dir]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def SLC_corners(SLC_par, terra_alt='-', kml='-', logpath=None, outdir=None, shellscript=None):
    """
    | Calculate SLC/MLI image corners in geodetic latitude and longitude (deg.)
    | Copyright 2022, Gamma Remote Sensing, v2.2 8-Nov-2022 clw/awi/cm
    
    Parameters
    ----------
    SLC_par:
        (input) ISP SLC/MLI image parameter file
    terra_alt:
        (input) average terrain altitude (enter - for default: 300.000 meters)
    kml:
        (output) kml output file (enter - for none)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/SLC_corners', SLC_par, terra_alt, kml]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def SLC_deramp(SLC1, SLC1_par, SLC2, SLC2_par, mode, dop_ph='-', logpath=None, outdir=None, shellscript=None):
    """
    | Calculate and subtract Doppler phase from an SLC image
    | Copyright 2023, Gamma Remote Sensing, v1.7 18-Apr-2023 clw
    
    Parameters
    ----------
    SLC1:
        (input) SLC data file (FCOMPLEX or SCOMPLEX format)
    SLC1_par:
        (input) SLC parameter file with Doppler information
    SLC2:
        (output) SLC with Doppler phase removed (or added)
    SLC2_par:
        (output) SLC parameter file for the output SLC
    mode:
        mode of operation:
            * 0: subtract Doppler phase ramp (deramp)
            * 1: add Doppler phase ramp (reramp)
    
    dop_ph:
        (output) Doppler phase (FLOAT) (enter - for none)
            Note: SLC1_par contains the Doppler polynomial that is used to calculate the Doppler phase ramp
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/SLC_deramp', SLC1, SLC1_par, SLC2, SLC2_par, mode, dop_ph]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def SLC_deramp_ScanSAR(SLC1_tab, SLC2_tab, mode, phflg='-', logpath=None, outdir=None, shellscript=None):
    """
    | Calculate and subtract ScanSAR or TOPS Doppler phase from burst SLC data
    | Copyright 2023, Gamma Remote Sensing, v2.1 18-Apr-2023 clw/cm
    
    Parameters
    ----------
    SLC1_tab:
        (input) 3 column list of input ScanSAR SLC, swaths are listed in order from near to far range:
            SLC_tab line entries:   SLC    SLC_par   TOPS_par
    SLC2_tab:
        (input) 3 column list of output ScanSAR SLC, swaths are listed in order from near to far range
    mode:
        mode of operation:
            * 0: subtract ScanSAR Doppler phase (deramp)
            * 1: add Doppler phase ramp (reramp)
    
    phflg:
        deramp phase flag (enter - for default)
            * 0: do not save ScanSAR Doppler phase (default)
            * 1: save ScanSAR Doppler phase, output filename is the same as the deramped SLC with extension .dph
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/SLC_deramp_ScanSAR', SLC1_tab, SLC2_tab, mode, phflg]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def SLC_deskew(SLC1, SLC1_par, SLC2, SLC2_par, mode='-', interp='-', order='-', deramp='-', ph_corr='-', sr0='-', sr2='-', logpath=None, outdir=None, shellscript=None):
    """
    | Change geometry from Doppler centroid to zero-Doppler (deskew) or vice-versa
    | Copyright 2024, Gamma Remote Sensing, v1.6 17-Oct-2024 cm/clw/uw
    
    Parameters
    ----------
    SLC1:
        (input) SLC image file (FCOMPLEX or SCOMPLEX format)
    SLC1_par:
        (input) SLC1 ISP image parameter file
    SLC2:
        (output) SLC image file in new geometry
    SLC2_par:
        (output) SLC2 ISP image parameter file
    mode:
        mode of operation (enter - for default)
            * 0: change geometry from Doppler centroid to zero-Doppler (deskew, default)
            * 1: change geometry from zero-Doppler to Doppler centroid (reskew)
    
    interp:
        interpolation method (enter - for default)
            * 0: Lanczos interpolation (default)
            * 1: B-spline interpolation
    
    order:
        Lanczos interpolator order / B-spline degree 4 -> 9 (enter - for default: 4)
    deramp:
        deramp flag (enter - for default)
            * 0: do not deramp and reramp data
            * 1: deramp data before interpolation and reramp afterwards (default)
    
    ph_corr:
        range shift phase correction flag (enter - for default)
            * 0: do not correct phase related to range shift
            * 1: correct phase related to range shift (default)
    
    sr0:
        near range distance of the resampled image in meter (enter - for default: calculated from input)
    sr2:
        far range distance of the resampled image in meter (enter - for default: calculated from input)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/SLC_deskew', SLC1, SLC1_par, SLC2, SLC2_par, mode, interp, order, deramp, ph_corr, sr0, sr2]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def SLC_freq_shift(SLC, SLC_par, SLC_shift, SLC_shift_par, freq_shift, logpath=None, outdir=None, shellscript=None):
    """
    | ISP Program GAMMA_SOFTWARE-20250625/ISP/bin/SLC_freq_shift
    | Shift the effective radar carrier frequency of an SLC image by a specified amount
    | Copyright 2022, Gamma Remote Sensing, v1.1 8-Nov-2022 clw
    
    Parameters
    ----------
    SLC:
        (input) SLC file (FCOMPLEX or SCOMPLEX)
    SLC_par:
        (input) SLC parameter file
    SLC_shift:
        (output) SLC data with shifted radar carrier frequency
    SLC_shift_par:
        (output) SLC parameter file with shifted radar carrier frequency
    freq_shift:
        radar carrier frequency shift (Hz)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/SLC_freq_shift', SLC, SLC_par, SLC_shift, SLC_shift_par, freq_shift]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def SLC_interp(SLC2, SLC1_par, SLC2_par, OFF_par, SLC2R, SLC2R_par, loff='-', nlines='-', mode='-', order='-', logpath=None, outdir=None, shellscript=None):
    """
    | SLC complex image resampling using 2-D Lanczos or B-spline interpolation
    | Copyright 2023, Gamma Remote Sensing, v4.9 18-Apr-2023 clw/cm
    
    Parameters
    ----------
    SLC2:
        (input) SLC2 image to be resampled to the geometry of the SLC1 reference image
    SLC1_par:
        (input) SLC1 ISP image parameter file
    SLC2_par:
        (input) SLC2 ISP image parameter file
    OFF_par:
        (input) ISP offset/interferogram parameter file
    SLC2R:
        (output) single-look complex image 2 coregistered to SLC1
    SLC2R_par:
        (output) SLC2R ISP image parameter file for coregistered image
    loff:
        offset to first valid output line (in SLC1 lines) (enter - for default: 0)
    nlines:
        number of valid output lines (enter - or 0 for default: to end of file)
    mode:
        interpolation mode (enter - for default)
            * 0: Lanczos (default)
            * 1: B-spline
    
    order:
        Lanczos interpolator order / B-spline degree 4 -> 9 (enter - for default: 4)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/SLC_interp', SLC2, SLC1_par, SLC2_par, OFF_par, SLC2R, SLC2R_par, loff, nlines, mode, order]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def SLC_interp_map(SLC2, SLC1_par, SLC2_par, OFF_par, SLC2R, SLC2R_par, OFF_par2, coffs_sm, loff='-', nlines='-', mode='-', order='-', logpath=None, outdir=None, shellscript=None):
    """
    | SLC image resampling using a 2-D offset map
    | Copyright 2024, Gamma Remote Sensing, v4.3 22-Aug-2024 clw/uw/cm
    
    Parameters
    ----------
    SLC2:
        (input) SLC2 image to be resampled to the reference SLC1 reference image
    SLC1_par:
        (input) SLC1 ISP image parameter file
    SLC2_par:
        (input) SLC2 ISP image parameter file
    OFF_par:
        (input) ISP offset/interferogram parameter file
    SLC2R:
        (output) single-look complex image 2 coregistered to SLC1
    SLC2R_par:
        (output) SLC2R ISP image parameter file for co-registered image
    OFF_par2:
        (input) ISP offset/interferogram parameter file used for residual offsets map (coffs_sm)
    coffs_sm:
        (input) smoothed residual range and azimuth offsets (fcomplex)
    loff:
        offset to first valid output line (in SLC1 lines) (enter - for default: 0)
    nlines:
        number of valid output lines (enter - or 0 for default: to end of file)
    mode:
        interpolation mode (enter - for default)
            * 0: Lanczos (default)
            * 1: B-spline
    
    order:
        Lanczos interpolator order / B-spline degree 4 -> 9 (enter - for default: 4)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/SLC_interp_map', SLC2, SLC1_par, SLC2_par, OFF_par, SLC2R, SLC2R_par, OFF_par2, coffs_sm, loff, nlines, mode, order]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def SLC_interp_ScanSAR(SLC2_tab, SLC2_par, SLC1_tab, SLC1_par, OFF_par, SLC2R_tab, SLC2R='-', SLC2R_par='-', mode='-', order='-', SLC2R_dir='-', burst_check='-', logpath=None, outdir=None, shellscript=None):
    """
    | Resample ScanSAR burst mode SLC using global offset polynomial
    | Copyright 2025, Gamma Remote Sensing v4.4 28-Jan-2025 clw/cm
    
    Parameters
    ----------
    SLC2_tab:
        (input) 3 column list of ScanSAR SLC2 swaths to be resampled into the geometry of SLC1 listed in order from near to far range
            SLC2_tab line entries:   SLC    SLC_par   TOPS_par
    SLC2_par:
        (input) SLC parameter file of ScanSAR SLC2 mosaic, SLC2 is generated from the ScanSAR swaths listed in SLC2_tab
    SLC1_tab:
        (input) 3 column list of the reference ScanSAR SLC swaths listed in order from near to far range
    SLC1_par:
        (input) SLC parameter file of the reference ScanSAR SLC1 mosaic, SLC1 is generated from the ScanSAR swaths listed in SLC1_tab
    OFF_par:
        (input) global ISP offset and interferogram parameter file, the offset model is determined from the ScanSAR SLC mosaics
            * NOTE: The OFF_par specifies the number of range and azimuth looks required to determine valid data bounds (burst windows)
    
    SLC2R_tab:
        (input/output) 3 column list of the resampled ScanSAR SLC2 swaths listed in order from near to far range
            * NOTE: If the SLC2R_tab does not yet exist, the entires will be created with file names derived from the filenames in SLC2_tab and the SLC2R_dir
              The file extensions of the new entries are changed from slc to rslc
    SLC2R:
        (output) mosaic generated from the resampled swaths listed in SLC2R_tab, coregistered to the reference mosaic of SLC1 (enter - for none)
    SLC2R_par:
        (output) SLC parameter file associated with the mosaic created from the resampled swaths SLC2R (enter - for none)
    mode:
        complex data interpolation mode (enter - for default)
            * 0: Lanczos (default)
            * 1: B-spline
    
    order:
        Lanczos interpolator order / B-spline degree 4 -> 9 (enter - for default: 4)
    SLC2R_dir:
        directory for resampled burst SLC2R data, ignored if the DIFF_tab already exists (enter - for default: current directory)
    burst_check:
        check and update burst parameters to match actual data (enter - for default)
            * 0: no (default)
            * 1: yes
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/SLC_interp_ScanSAR', SLC2_tab, SLC2_par, SLC1_tab, SLC1_par, OFF_par, SLC2R_tab, SLC2R, SLC2R_par, mode, order, SLC2R_dir, burst_check]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def SLC_intf(SLC1, SLC2R, SLC1_par, SLC2R_par, OFF_par, interf, rlks, azlks, loff='-', nlines='-', sps_flg='-', azf_flg='-', rp1_flg='-', rp2_flg='-', SLC1s='-', SLC2Rs='-', SLC_1s_par='-', SLC_2Rs_par='-', az_beta='-', logpath=None, outdir=None, shellscript=None):
    """
    | Calculate interferogram from co-registered SLC image data
    | Copyright 2024, Gamma Remote Sensing, v6.3 8-Mar-2024 clw/uw/cm
    
    Parameters
    ----------
    SLC1:
        (input) single-look complex image 1 (reference)
    SLC2R:
        (input) single-look complex image 2 coregistered to SLC1
    SLC1_par:
        (input) SLC1 ISP image parameter file
    SLC2R_par:
        (input) SLC2R ISP image parameter file for the co-registered image
    OFF_par:
        (input) ISP offset/interferogram parameter file
    interf:
        (output) interferogram from SLC1 and SLC2R
    rlks:
        number of range looks
    azlks:
        number of azimuth looks
    loff:
        offset to starting line relative to SLC1 for interferogram (enter - for default: 0)
    nlines:
        number of SLC lines to process (enter - for default: to end of file)
    sps_flg:
        range spectral shift flag (enter - for default)
            * 1: apply range spectral shift filter (default)
            * 0: do not apply range spectral shift filter
    
    azf_flg:
        azimuth common band filter flag (enter - for default)
            * 1: apply azimuth common-band filter (default)
            * 0: do not apply azimuth common band filter
    
    rp1_flg:
        SLC1 image range phase mode (enter - for default)
            * 0: non-zero Doppler geometry
            * 1: zero-Doppler geometry (default)
    
    rp2_flg:
        SLC2 image range phase mode (enter - for default)
            * 0: non-zero Doppler geometry
            * 1: zero-Doppler geometry (default)
    
    SLC1s:
        SLC1 after range spectral shift and azimuth common-band filtering (FCOMPLEX format) (enter - for none)
    SLC2Rs:
        SLC2R after range spectral shift and azimuth common-band filtering (FCOMPLEX format) (enter - for none)
    SLC_1s_par:
        SLC1s ISP image parameter file (enter - for none)
    SLC_2Rs_par:
        SLC2Rs ISP image parameter file (enter - for none)
    az_beta:
        azimuth common-band filter Kaiser window parameter (enter - for default: 2.120)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/SLC_intf', SLC1, SLC2R, SLC1_par, SLC2R_par, OFF_par, interf, rlks, azlks, loff, nlines, sps_flg, azf_flg, rp1_flg, rp2_flg, SLC1s, SLC2Rs, SLC_1s_par, SLC_2Rs_par, az_beta]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def SLC_intf2(SLC1, SLC2R, SLC1_par, SLC2R_par, MLI1, MLI2R, MLI1_par, MLI2R_par, interf, cc, r_dec, az_dec, rwin='-', azwin='-', wflg='-', n_ovr='-', sim_phase='-', lanczos='-', beta='-', logpath=None, outdir=None, shellscript=None):
    """
    | Calculate interferogram and MLI images from SLCs with separate averaging window dimensions and decimation factors
    | Copyright 2025, Gamma Remote Sensing, v2.3 20-May-2025 clw/cm/of
    
    Parameters
    ----------
    SLC1:
        (input) single-look complex image 1 (reference)
    SLC2R:
        (input) single-look complex image 2 coregistered to SLC1
    SLC1_par:
        (input) SLC1 image parameter file
    SLC2R_par:
        (input) SLC2R image parameter file for the co-registered image
    MLI1:
        (output) multi-look intensity image derived from SLC1 (enter - for none)
    MLI2R:
        (output) multi-look intensity image derived from SLC2R (enter - for none)
    MLI1_par:
        (output) MLI image parameter file derived from SLC1_par (enter - for none)
    MLI2R_par:
        (output) MLI image parameter file derived from SLC2R_par (enter - for none)
    interf:
        (output) complex interferogram from SLC1 and SLC2R  (enter - for none)
    cc:
        (output) interferometric correlation magnitude of SLC1 and SLC2R (enter - for none)
    r_dec:
        range decimation factor (int)
    az_dec:
        azimuth decimation factor (int)
    rwin:
        averaging window width (int) (enter - for default: r_dec)
    azwin:
        averaging window height (int) (enter - for default: az_dec)
    wflg:
        window weighting function (enter - for default):
            * 0: rectangular (default)
            * 1: Kaiser
            * 2: circular Gaussian 
    
    n_ovr:
        oversampling factor 1 -> 2 (enter - for default: 1)
    sim_phase:
        (input) simulated interferometric phase, coregistered MLI1 (FLOAT, enter - for none)
    lanczos:
        Lanczos interpolator order 5 -> 9 (enter - for default: 7)
    beta:
        Gaussian or Kaiser window parameter (enter - for default: 2.0)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/SLC_intf2', SLC1, SLC2R, SLC1_par, SLC2R_par, MLI1, MLI2R, MLI1_par, MLI2R_par, interf, cc, r_dec, az_dec, rwin, azwin, wflg, n_ovr, sim_phase, lanczos, beta]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def SLC_mosaic_range(SLC_tab, SLC, SLC_par, mode='-', order='-', logpath=None, outdir=None, shellscript=None):
    """
    | Calculate mosaic of Stripmap SLC data provided in multiple pieces in range direction (e.g. PALSAR-3)
    | Copyright 2025, Gamma Remote Sensing v1.1 5-Feb-2025 cm/clw/uw
    
    Parameters
    ----------
    SLC_tab:
        (input) 2 column list of Stripmap SLC pieces (from near to far range)
            SLC_tab line entries:   SLC   SLC_par
    SLC:
        (output) SLC mosaic image
    SLC_par:
        (output) SLC mosaic image parameter file
    mode:
        complex data interpolation mode in range (enter - for default)
            * 0: Lanczos (default)
            * 1: B-spline
            * 2: nearest neighbor
    
    order:
        Lanczos interpolator order / B-spline degree 4 -> 9 (enter - for default: 4)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/SLC_mosaic_range', SLC_tab, SLC, SLC_par, mode, order]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def SLC_mosaic_ScanSAR(SLC_tab, SLC, SLC_par, rlks, azlks, bflg='-', SLCR_tab='-', logpath=None, outdir=None, shellscript=None):
    """
    | Calculate SLC mosaic of ScanSAR SLC burst data (Sentinel-1, TerraSAR-X, RCM...)
    | Copyright 2025, Gamma Remote Sensing v5.0 14-Jan-2025 clw/awi/cm
    
    Parameters
    ----------
    SLC_tab:
        (input) 3 column list of ScanSAR SLC, swaths are listed in order from near to far range
            SLC_tab line entries:   SLC   SLC_par   TOPS_par
    SLC:
        (output) SLC mosaic image
    SLC_par:
        (output) SLC mosaic image parameter file
    rlks:
        number of range looks used to determine burst window boundaries for the mosaic
    azlks:
        number of azimuth looks used to determine burst window boundaries for the mosaic
    bflg:
        burst window calculation flag (enter - for default)
            * 0: use existing burst window parameters if they exist, otherwise calculate burst window parameters (default)
            * 1: calculate burst window parameters from burst parameters and the number of range and azimuth looks
    
    SLCR_tab:
        (input) 3 column list of the reference scene, swaths are listed in order from near to far range (enter - for none)
            SLCR_tab line entries:   SLC   SLC_par   TOPS_par
            * NOTE: When generating a mosaic of a resampled SLC, the SLC_tab of the reference scene is required
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/SLC_mosaic_ScanSAR', SLC_tab, SLC, SLC_par, rlks, azlks, bflg, SLCR_tab]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def SLC_ovr(SLC, SLC_par, SLC_ovr, SLC_ovr_par, r_ovr='-', az_ovr='-', mode='-', order='-', deramp='-', logpath=None, outdir=None, shellscript=None):
    """
    | Oversample SLC data in range and azimuth using 2-D Lanczos or B-spline interpolation
    | Copyright 2024, Gamma Remote Sensing, v1.6 1-Feb-2024 clw/cm
    
    Parameters
    ----------
    SLC:
        (input) SLC image  (FCOMPLEX or SCOMPLEX format)
    SLC_par:
        (input) SLC image parameter file
    SLC_ovr:
        (output) oversampled SLC image
    SLC_ovr_par:
        (output) oversampled SLC image parameter file
    r_ovr:
        range oversampling factor (enter - for default: 1.0)
    az_ovr:
        azimuth oversampling factor (enter - for default: 1.0)
    mode:
        interpolation mode (enter - for default)
            * 0: Lanczos interpolation (default)
            * 1: B-spline interpolation
    
    order:
        Lanczos interpolator order / B-spline degree 4 -> 9 (enter - for default: 4)
    deramp:
        deramp flag (enter - for default)
            * 0: do not deramp and reramp data
            * 1: deramp data before interpolation and reramp afterwards (default)
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/SLC_ovr', SLC, SLC_par, SLC_ovr, SLC_ovr_par, r_ovr, az_ovr, mode, order, deramp]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def SLC_phase_shift(SLC1, SLC1_par, SLC2, SLC2_par, ph_shift, logpath=None, outdir=None, shellscript=None):
    """
    | Add a constant phase from an SLC image
    | Copyright 2023, Gamma Remote Sensing, v1.3 24-Apr-2023 clw
    
    Parameters
    ----------
    SLC1:
        (input) SLC data file (fcomplex or scomplex format)
    SLC1_par:
        (input) SLC parameter file
    SLC2:
        (output) SLC with phase shift
    SLC2_par:
        (output) SLC parameter file for the output SLC
    ph_shift:
        phase shift to add to SLC phase (radians)
            * NOTE: Used to apply a constant phase shift of -1.25 radians to Sentinel-1 TOPS SLC data
              from swath IW1 acquired up to 10-Mar-2015.
              Used to apply a constant phase shift of -3.83 radians to Sentinel-1 TOPS SLC data with
              H-POL on receive (e.g. VH) acquired up to 10-Mar-2015.
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/SLC_phase_shift', SLC1, SLC1_par, SLC2, SLC2_par, ph_shift]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def SLC_RFI_filt(SLC, SLC_par, SLC_filt, rfi_thres='-', nfft_r='-', nfft_az='-', r_step='-', az_step='-', mwin_r='-', mwin_az='-', logpath=None, outdir=None, shellscript=None):
    """
    | Adaptive RFI filtering for SLC image using median spectral filtering
    | Copyright 2023, Gamma Remote Sensing, v1.6 18-Apr-2023 clw
    
    Parameters
    ----------
    SLC:
        (input) SLC to be filtered (FCOMPLEX or SCOMPLEX)
    SLC_par:
        (input) reference SLC parameter file
    SLC_filt:
        (output) output filtered SLC using the power spectrum of the reference SLC
    rfi_thres:
        RFI threshold (enter - for default: 10.00)
    nfft_r:
        range filter FFT window size, 2\\*\\*N, 16->1024, (enter - for default: 128)
    nfft_az:
        azimuth filter FFT window size, 2\\*\\*N, 16->1024, (enter - for default: 128)
    r_step:
        range processing step (enter - for default: nfft_r/4)
    az_step:
        azimuth processing step (enter - for default: nfft_az/4)
    mwin_r:
        range median window size for median PSD filtering (enter - for default: 5)
    mwin_az:
        azimuth median window size for median PSD filtering (enter - for default: 5)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/SLC_RFI_filt', SLC, SLC_par, SLC_filt, rfi_thres, nfft_r, nfft_az, r_step, az_step, mwin_r, mwin_az]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def SLC_RFI_filt2(SLC, SLC_par, SLC_filt, rfi_thres='-', method='-', f_bs='-', bs_width='-', roff='-', nr='-', azoff='-', naz='-', pltflg='-', logpath=None, outdir=None, shellscript=None):
    """
    | RFI filtering for SLC image using a band-stop filter
    | Copyright 2024, Gamma Remote Sensing, v1.5 2-Feb-2024 cm
    
    Parameters
    ----------
    SLC:
        (input) SLC to be filtered (FCOMPLEX or SCOMPLEX)
    SLC_par:
        (input) reference SLC parameter file
    SLC_filt:
        (output) output filtered SLC (same format as SLC)
    rfi_thres:
        RFI threshold in dB above reference (enter - for default: auto)
    method:
        RFI detection method (enter - for default)
            * 0: threshold above median
            * 1: threshold using spectrum symmetry (default)
    
    f_bs:
        center or seed frequency of band-stop filter in Hz (-fadc/2.0 <= f_bs < fadc/2.0, enter - for default: auto)
    bs_width:
        width of band-stop filter in Hz (enter - for default: auto)
    roff:
        offset to starting range sample to filter (enter - for default: 0)
    nr:
        number of range samples to filter (enter - for default: to end of line)
    azoff:
        offset to starting azimuth line to filter (enter - for default: 0)
    naz:
        number of azimuth lines to filter (enter - for default: to end of file)
    pltflg:
        range spectrum plotting flag (enter - for default)
            * 0: none
            * 1: output plot in PNG format (default)
            * 2: screen output plot
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/SLC_RFI_filt2', SLC, SLC_par, SLC_filt, rfi_thres, method, f_bs, bs_width, roff, nr, azoff, naz, pltflg]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def split_WB(data_in, data_par_in, data_tab, dtype, logpath=None, outdir=None, shellscript=None):
    """
    | ISP: Program GAMMA_SOFTWARE-20250625/ISP/bin/split_WB
    | Split WB mosaic image into individual beams using ISP parameter files
    | Copyright 2022, Gamma Remote Sensing, v1.4 8-Nov-2022 clw/cm
    
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
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/split_WB', data_in, data_par_in, data_tab, dtype]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def SR_to_GRD(MLI_par, OFF_par, GRD_par, in_file, out_file, rlks='-', azlks='-', interp_mode='-', grd_rsp='-', grd_azsp='-', degree='-', logpath=None, outdir=None, shellscript=None):
    """
    | Conversion to ground range for ISP MLI and INSAR data of type FLOAT
    | Copyright 2023, Gamma Remote Sensing, v2.5 18-Apr-2023 uw/clw/cm
    
    Parameters
    ----------
    MLI_par:
        (input) MLI image parameter file of the slant-range image
    OFF_par:
        (input) ISP OFF_par of the input image (enter - when the image geometry specified by the MLI_par)
    GRD_par:
        (input/output) image parameter file of output ground range image
    in_file:
        (input) slant range image (FLOAT)
    out_file:
        (output) ground range image (FLOAT)
    rlks:
        multi-looking in range (prior to resampling, enter - for default: 1)
    azlks:
        multi-looking in azimuth (prior to resampling, enter - for default: 1)
    interp_mode:
        interpolation mode (enter - for default)
            * 0: nearest-neighbor
            * 1: bicubic spline
            * 2: bicubic spline log(x)
            * 3: bicubic spline sqrt(x)
            * 4: B-spline interpolation (default B-spline degree: 3)
            * 5: B-spline interpolation sqrt(x) (default) (default B-spline degree: 3)
            * NOTE: log and sqrt interpolation modes should only be used with non-negative data!
    
    grd_rsp:
        output image ground range sample spacing (m) (enter - for default: (input image azimuth spacing) \\* azlks)
    grd_azsp:
        output image azimuth sample spacing (m) (enter - for default: (input image azimuth spacing) \\* azlks)
    degree:
        B-spline degree (2->9) (enter - for default: 3)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/SR_to_GRD', MLI_par, OFF_par, GRD_par, in_file, out_file, rlks, azlks, interp_mode, grd_rsp, grd_azsp, degree]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def subtract_phase(interf_in, phase_file, interf_out, width, factor='-', logpath=None, outdir=None, shellscript=None):
    """
    | ISP: Program GAMMA_SOFTWARE-20250625/ISP/bin/subtract_phase
    | Subtract scaled phase image from a complex interferogram
    | Copyright 2023, Gamma Remote Sensing, v3.3 19-Apr-2023 uw/clw
    
    Parameters
    ----------
    interf_in:
        (input) input interferogram (FCOMPLEX)
    phase_file:
        (input) unwrapped interferometric phase (FLOAT)
    interf_out:
        (output) output interferogram (input interferogram - scaled phase) (FCOMPLEX)
    width:
        number of samples/line
    factor:
        constant scale factor for input phase data (enter - for default: 1.0)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/subtract_phase', interf_in, phase_file, interf_out, width, factor]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def tree_cc(flag, width, mbl='-', xmin='-', xmax='-', ymin='-', ymax='-', logpath=None, outdir=None, shellscript=None):
    """
    | Phase unwrapping tree generation with low correlation search (modified ARW algorithm)
    | Copyright 2023, Gamma Remote Sensing, v3.1 18-Apr-2023 clw/uw
    
    Parameters
    ----------
    flag:
        (input) phase unwrapping flag file
    width:
        number of samples/row
    mbl:
        maximum branch length (enter - for default: 32, maximum=64) 
    xmin:
        starting range pixel offset (enter - for default: 0)
    xmax:
        last range pixel offset (enter - for default: width-1)
    ymin:
        starting azimuth row, relative to start (enter - for default: 0)
    ymax:
        last azimuth row, relative to start (enter - for default: nlines-1)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/tree_cc', flag, width, mbl, xmin, xmax, ymin, ymax]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def tree_gzw(flag, width, mbl='-', xmin='-', xmax='-', ymin='-', ymax='-', logpath=None, outdir=None, shellscript=None):
    """
    | Phase unwrapping tree generation (GZW algorithm)
    | Copyright 2023, Gamma Remote Sensing, v3.8 18-Apr-2023 clw/uw
    
    Parameters
    ----------
    flag:
        (input) phase unwrapping flag file
    width:
        number of samples/row
    mbl:
        maximum branch length (enter - for default: 32)
    xmin:
        starting range pixel offset (enter - for default: 0)
    xmax:
        last range pixel offset (enter - for default: width-1)
    ymin:
        starting azimuth row, relative to start (enter - for default: 0)
    ymax:
        last azimuth row, relative to start (enter - for default: nlines-1)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/tree_gzw', flag, width, mbl, xmin, xmax, ymin, ymax]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def unw_model(interf, unw_model, unw, width, xinit='-', yinit='-', ref_ph='-', width_model='-', logpath=None, outdir=None, shellscript=None):
    """
    | Phase unwrapping using a model of the unwrapped phase
    | Copyright 2023, Gamma Remote Sensing, v1.9 21-Sep-2023 clw/uw
    
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
        offset to phase reference location in range (col) (enter - for default: 0)
    yinit:
        offset to phase reference location in azimuth (row) (enter - for default: 0)
    ref_ph:
        reference point phase (radians) (enter - for phase at the reference point)
    width_model:
        number of samples/row of the unwrapped phase model (enter - for default: interferogram width)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/bin/unw_model', interf, unw_model, unw, width, xinit, yinit, ref_ph, width_model]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def bpf_ssi(SLC, SLC_par, SLC_flow, SLC_flow_par, SLC_fhigh, SLC_fhigh_par, rbs='-', logpath=None, outdir=None, shellscript=None):
    """
    | bpf_ssi: Apply band-pass filtering for split-spectrum interferometry
    | Copyright 2023 Gamma Remote Sensing, v1.4 19-Apr-2023 uw/cm
    
    Parameters
    ----------
    SLC:
        (input) SLC (FCOMPLEX or SCOMPLEX, SLC should not be resampled)
    SLC_par:
        (input) SLC parameter file
    SLC_flow:
        (output) low frequency band filtered SLC (FCOMPLEX or SCOMPLEX)
    SLC_flow_par:
        (output) low frequency band filtered SLC parameter file
    SLC_fhigh:
        (output) high frequency band filtered SLC (FCOMPLEX or SCOMPLEX)
    SLC_fhigh_par:
        (output) high frequency band filtered SLC parameter file (FCOMPLEX or SCOMPLEX)
    rbs:
        relative range spectrum band separation (enter - for default: 0.6666 --> lowest and highest third of processing bandwidth)
            indicate - for the output files to only calculate filtering parameters
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/bpf_ssi', SLC, SLC_par, SLC_flow, SLC_flow_par, SLC_fhigh, SLC_fhigh_par, rbs]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def bpf_ssi_S1(SLC_tab, SLC_tab_flow, SLC_tab_high, rbs='-', logpath=None, outdir=None, shellscript=None):
    """
    | bpf_ssi_S1: Apply band-pass filtering for split-spectrum interferometry for S1 TOPS data
    | Copyright 2023 Gamma Remote Sensing, v1.2 19-Apr-2023 uw/cm
    
    Parameters
    ----------
    SLC_tab:
        (input) SLC_tab
    SLC_tab_flow:
        (output) output SLC_tab filename for low frequency band filtered SLC
    SLC_tab_high:
        (output) output SLC_tab filename for high frequency band filtered SLC
    rbs:
        relative range spectrum band separation (enter - for default: 0.6666 --> lowest and highest third of processing bandwidth)
            indicate - for the output files to only calculate filtering parameters
            The filename in SLC_tab_flow and SLC_tab_high are automatically generated by adding .flow and .fhigh
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/bpf_ssi_S1', SLC_tab, SLC_tab_flow, SLC_tab_high, rbs]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def get_GAMMA_RASTER(mode, logpath=None, outdir=None, shellscript=None):
    """
    | Script to determine the default extension for raster images or the operating system type
    | Copyright 2019 Gamma Remote Sensing, v1.3 1-Apr-2019 clw/uw/cm
    
    Parameters
    ----------
    mode:
        Specify the script string output:
            * 0: raster file extension (ras, bmp, or tif)
            * 1: OS type: Linux, MINGW64_NT-10.0, CYGWIN_NT-10.0, darwin...
            * NOTE: The default raster format on Linux systems is SUN_RASTER (\\*.ras), for all other operating systems it is BMP (\\*.bmp).
              SUN_RASTER and BMP images are limited in size to 32767 x 32767. TIFF files do not have this limitation.
              To set the default image raster format for Gamma programs, set the environment variable GAMMA_RASTER as follows:
              bash:
              export GAMMA_RASTER=SUN_RASTER  #extension: ras
              export GAMMA_RASTER=BMP         #extension: bmp
              export GAMMA_RASTER=TIFF        #extension: tif
              csh,tcsh:
              setenv GAMMA_RASTER SUN_RASTER  #extension: ras
              setenv GAMMA_RASTER BMP         #extension: bmp
              setenv GAMMA_RASTER TIFF        #extension: tif
              Environment variables can be set either in processing scripts, or in the shell initialization file (e.g. .bashrc)
              Programs in the Gamma software that generate raster image files query the value of GAMMA_RASTER if it has been defined.
              This script can be called from within another script to determine the default raster image format or OS type:
              bash:        $ext=`get_GAMMA_RASTER 0`
              csh,tcsh: set ext=`get_GAMMA_RASTER 0`
              The variable $ext can then be used to specify the format of the output raster file by using it to construct
              the output file name: 
              bash:        $my_raster=$my_name"."$ext
              csh/tcsh: set my_raster=$my_name"."$ext
              OS: Linux
              GAMMA_RASTER value: TIFF
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/get_GAMMA_RASTER', mode]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def INTF_SLC(pass1, pass2, rlks, azlks, algorithm='-', cc_win='-', r_pos='-', az_pos='-', logpath=None, outdir=None, shellscript=None):
    """
    | INTF_SLC: calculate interferogram, co-registered SLC, intensity images, and correlation
    | Copyright 2023 Gamma Remote Sensing, v1.2 18-Apr-2023 clw/uw/cm
    
    Parameters
    ----------
    pass1:
        pass 1 identifier (example: pass number) reference
    pass2:
        pass 2 identifier (example: pass number)
    rlks:
        number of range looks
    azlks:
        number of azimuth looks
    algorithm:
        algorithm used to determine offsets (enter - for default)
            * 1: intensity image cross correlation (default)
            * 2: fringe visibility
    
    cc_win:
        window used for estimation of the correlation coefficient (enter - for default: 3)
    r_pos:
        range position of center of image patch for initial offset (enter - for default: image center)
    az_pos:
        azimuth position of center of image patch for initial offset (enter - for default: image center)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/INTF_SLC', pass1, pass2, rlks, azlks, algorithm, cc_win, r_pos, az_pos]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def ionosphere_check(SLC, par, rwin='-', azwin='-', thresh='-', rstep='-', azstep='-', cleaning='-', use_existing='-', logpath=None, outdir=None, shellscript=None):
    """
    | ionosphere_check: Determine azimuth spectrum sub-band range and azimuth offsets of a single SLC
    | Significant non-zero azimuth offsets are a clear indication for the presence of ionospheric effects
    | Copyright 2024 Gamma Remote Sensing, v1.8 11-Dec-2024 uw/cm
    
    Parameters
    ----------
    SLC:
        (input) SLC image (e.g. 20070214.slc)
    par:
        (input) SLC parameter file (e.g. 20070214.slc.par)
    rwin:
        range window size used in offset estimation (enter - for default: 256)
    azwin:
        azimuth window size used in offset estimation (enter - for default: 256)
    thresh:
        threshold value used in offset estimation (enter - for default: 0.1)
    rstep:
        range step used in offset estimation (enter - for default: rwin/4)
    azstep:
        azimuth step used in offset estimation (enter - for default: azwin/4)
    cleaning:
        cleaning flag (enter - for default)
            * 0: no cleaning, keep intermediate files
            * 1: delete intermediate files (default)
    
    use_existing:
        use files generated in a previous run to speed up processing (enter - for default)
            * 0: no (default)
            * 1: yes
              
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/ionosphere_check', SLC, par, rwin, azwin, thresh, rstep, azstep, cleaning, use_existing]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def make_tab(list, tab, template, logpath=None, outdir=None, shellscript=None):
    """
    | GAMMA_SOFTWARE-20250625/ISP/scripts/make_tab
    | Generate a table file from a list or multi-colum table using a text template
    | Copyright 2024, Gamma Remote Sensing, v1.1 22-Apr-2024 cm/clw/uw
    
    Parameters
    ----------
    list:
        (input) list or multi-column table (text)
    tab:
        (output) table file (text)
    template:
        template definition used to generate a line of the output table, entered between single quotes.
            Placeholders , , ... specify the columns of the input table.
            (example 1: '$1.slc $1.slc.par')
            (example 2: '$1_$2.base $1_$2.off')
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/make_tab', list, tab, template]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def mk_ptarg(RSLC_tab, cal_dir, r_samp, az_samp, osf='-', options='-', logpath=None, outdir=None, shellscript=None):
    """
    | GAMMA_SOFTWARE-20250625/ISP/scripts/mk_ptarg
    | Copyright 2023, Gamma Remote Sensing, v1.6 18-Apr-2023 clw
    | Perform point target analysis on a stack of coregistered SLCs
    
    Parameters
    ----------
    RSLC_tab:
        (input) two column list of coregistered SLC filenames and SLC parameter filenames (including paths) (ascii)
            1. SLC filename  (includes path)
            2. SLC parameter filename (includes path)
    cal_dir:
        directory for output calibration results
    r_samp:
        (input) calibration target range sample number
    az_samp:
        (input) calibration target azimuth line number
    osf:
        SLC over-sampling factor 2, 4, 8, 16, 32, 64 (enter - for default: 16)
            -s scale  (option) set image display scale factor (default: 0.3)
            -e exp    (option) set image display exponent (default: 0.5)
    options:
        not documented
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/mk_ptarg', RSLC_tab, cal_dir, r_samp, az_samp, osf, options]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def mk_ptarg_cal(CR_tab, SLC, SLC_par, cal_dir, sigma, c_rpos, c_azpos, osf='-', options='-', logpath=None, outdir=None, shellscript=None):
    """
    | GAMMA_SOFTWARE-20250625/ISP/scripts/mk_ptarg_cal
    | Copyright 2023, Gamma Remote Sensing, v2.1 18-Apr-2023 clw
    | Perform point target analysis and calibration factor evaluation for a set of point targers
    
    Parameters
    ----------
    CR_tab:
        (input) 3 column list of row and sample number of corner reflectors
            1. Corner reflector id
            2. SLC column  (includes path)
            3. SLC row    (includes path)
    SLC:
        SLC image
    SLC_par:
        SLC_parameter file
    cal_dir:
        directory for output calibration results
    sigma:
        Radar cross-section of the corner reflectors
    c_rpos:
        range sample number of the center of the region used to estimate region
    c_azpos:
        azimuth line of the center of the region used to estimate clutter
    osf:
        SLC over-sampling factor 2, 4, 8, 16, 32, 64 (enter - for default: 16)
            -s scale  (option) set image display scale factor (default: 0.2)
            -e exp    (option) set image display exponent (default: 0.5)
    options:
        not documented
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/mk_ptarg_cal', CR_tab, SLC, SLC_par, cal_dir, sigma, c_rpos, c_azpos, osf, options]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def mk_tab3(dir, ext1, ext2, ext3, tab, logpath=None, outdir=None, shellscript=None):
    """
    | Copyright 2023, Gamma Remote Sensing, v1.1 24-Apr-2023 clw
    | Generate SLC_tab, MLI_tab, or RAW_list for processing
    
    Parameters
    ----------
    dir:
        (input) directory including paths that contain the data files
    ext1:
        (input) pattern to select data files (examples: slc, raw...), (enter - for all files in the directory)
    ext2:
        (input) pattern to select parameter files that match the data (enter - for none, examples: slc.par, raw_par, raw.par)
    ext3:
        (input) pattern to select parameter files that match the data (enter - for none, examples: ppar)
    tab:
        (output) list of data filenames and associated parameter files (including paths) (text)
            * NOTE: The current directory is denoted using . 
              
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/mk_tab3', dir, ext1, ext2, ext3, tab]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def offset_plot_az(offset, r_min, r_max, r_plot, az_plot, logpath=None, outdir=None, shellscript=None):
    """
    | IPTA script: GAMMA_SOFTWARE-20250625/ISP/scripts/offset_plot_az
    | Copyright 2023, Gamma Remote Sensing, v1.4 17-Apr-2023 clw
    | extract range and azimuth offsets for a range window from an text offset file
    
    Parameters
    ----------
    offset:
        (input) list of range and azimuth offsets generated by offset_pwr (text)
    r_min:
        minimum range pixel number to extract range and azimuth offsets
    r_max:
        minimum range pixel number to extract range and azimuth offsets
    r_plot:
        range offsets xmgrace plot file
    az_plot:
        azimuth offsets xmgrace plot file
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/offset_plot_az', offset, r_min, r_max, r_plot, az_plot]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def offset_plot_r(offset, az_min, az_max, r_plot, az_plot, logpath=None, outdir=None, shellscript=None):
    """
    | IPTA script: GAMMA_SOFTWARE-20250625/ISP/scripts/offset_plot_r
    | Copyright 2004, Gamma Remote Sensing, v1.3 17-Jan-2005 clw
    | extract range and azimuth offsets for an azimuth window from an text offset file
    
    Parameters
    ----------
    offset:
        (input) list of range and azimuth offsets generated by offset_pwr (text)
    az_min:
        minimum azimuth line number to extract range and azimuth offsets
    az_max:
        minimum azimuth line number to extract range and azimuth offsets
    r_plot:
        range offsets xmgrace plot file
    az_plot:
        azimuth offsets xmgrace plot file       
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/offset_plot_r', offset, az_min, az_max, r_plot, az_plot]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def OPOD_vec(SLC_par, OPOD_dir, nstate='-', logpath=None, outdir=None, shellscript=None):
    """
    | GAMMA_SOFTWARE-20250625/ISP/scripts/OPOD_vec
    | Copyright 2025, Gamma Remote Sensing, v1.7 4-Feb-2025 clw/awi/cm
    | Extract Sentinel-1 state vectors from an OPOD file and write these state vectors to an SLC parameter file
    
    Parameters
    ----------
    SLC_par:
        (input/output) ISP SLC/MLI image parameter file
    OPOD_dir:
        (input) directory containing Sentinel-1 precise or restituted OPOD orbit data files (AUX_POEORB or AUX_RESORB)
            orbit files can be downloaded from https://s1qc.asf.alaska.edu/ or https://dataspace.copernicus.eu/
    nstate:
        number of state vectors to extract (enter - for default: include 60 sec extention at the start and end of the SLC data)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/OPOD_vec', SLC_par, OPOD_dir, nstate]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def run_all(list, command, log='-', logpath=None, outdir=None, shellscript=None):
    """
    | GAMMA_SOFTWARE-20250625/ISP/scripts/run_all
    | Run a single command iterating over arguments constructed from the elements of a list or multi-column table
    | Copyright 2025, Gamma Remote Sensing, v1.7 11-Mar-2025 clw/cm
    
    Parameters
    ----------
    list:
        (input) list or multi-column table (text)
    command:
        command template, entered between single quotes. Command arguments are constructed 
            with placeholders $1, $2, ... that specify the columns of the input table.
            (example 1: 'multi_look $1.slc $1.slc.par $1.mli $1.mli.par 5 1')
            (example 2: 'cp -r $1 $2')
    log:
        (output) log file that captures all screen output (both stdout and stderr) (enter - for none)
            Example: run_all dates 'multi_look $1.slc $1.slc.par $1.mli $1.mli.par 5 1' log
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/run_all', list, command, log]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def S1_BURST_tab(SLC1_tab, SLC2_tab, BURST_tab, logpath=None, outdir=None, shellscript=None):
    """
    | GAMMA_SOFTWARE-20250625/ISP/scripts/S1_BURST_tab
    | Copyright 2023, Gamma Remote Sensing, v1.5 18-Apr-2023 clw/cm
    | Calculate Sentinel BURST_tab based on parameters extracted from SLC parameter files listed in SLC1_tab and SLC2_tab
    | Running SLC_copy_ScanSAR using BURST_tab will generate SLC2 data with matching bursts for each swath of SLC1 and SLC2
    
    Parameters
    ----------
    SLC1_tab:
        (input) 3 column list of the reference TOPS SLC swaths in row order IW1, IW2, IW3
    SLC2_tab:
        (input) 3 column list of TOPS SLC2 swaths to be resampled to the geometry of the reference SLC1 in row order IW1, IW2, IW3.
    BURST_tab:
        (output) 2 column list of the first and last bursts to copy from each swath, one line for each swath
            BURST_tab line entries: first_burst  last_burst    Note: first burst is 1
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/S1_BURST_tab', SLC1_tab, SLC2_tab, BURST_tab]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def S1_BURST_tab_from_zipfile(zipfile_list, zipfile_ref, burst_number_table_ref='-', cleaning='-', logpath=None, outdir=None, shellscript=None):
    """
    | S1_BURST_tab_from_zipfile: Script used to generate S1_BURST_tab to support burst selection
    | Copyright 2021 Gamma Remote Sensing, v1.8 26-Jan-2021 uw/cm
    | 
    | NOTE: S1_BURST_tab_from_zipfile now calls S1_BURST_tab_from_zipfile.py
    | Using directly S1_BURST_tab_from_zipfile.py gives access to
    | additional useful options and is therefore recommended.
    | 
    
    Parameters
    ----------
    zipfile_list:
        (input) ASCII file containing S1 zip filename(s) of one data take
            indicate - to generate burst_number_table of reference TOPS SLC
    zipfile_ref:
        (input) S1 zip filename for the reference TOPS SLC
    burst_number_table_ref:
        (input) ASCII file containing first/last burst numbers selected
            indicate - to use all bursts as present in the reference TOPS SLC zipfile
    cleaning:
        flag to indicate if intermediate files are deleted (default=1: yes, 0: not deleted)
            intermediate and output filenames are generated based on the zip file names
            
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/S1_BURST_tab_from_zipfile', zipfile_list, zipfile_ref, burst_number_table_ref, cleaning]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def S1_extract_png(zipfile, logpath=None, outdir=None, shellscript=None):
    """
    | S1_extract_png: Script used to extract (and rename) quicklook (png file) from a S1 ZIP file
    | Copyright 2019 Gamma Remote Sensing, v1.1 22-Mar-2019 uw/cm
    
    Parameters
    ----------
    zipfile:
        (input) Sentinel-1 zipfile (GRD or SLC)
            
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/S1_extract_png', zipfile]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def S1_GRD_preproc(S1_list, MLI_dir, pol, log, options='-', logpath=None, outdir=None, shellscript=None):
    """
    | Preprocessing of Sentinel-1 TOPS GRD products, extract GRD data and generate MLI products
    | Copyright 2023, Gamma Remote Sensing, v1.3 18-Apr-2023 clw/cm
    
    Parameters
    ----------
    S1_list:
        (input) single column text file. Entries are directories (including path) containing Sentinel-1 TOPS GRD products
    MLI_dir:
        directory for output SLC data files and SLC parameter files
            * NOTE: output file names have the form : 20150119_hh.mli
    
    pol:
        SLC polarization to extract (hh,hv,vh,vv)
    log:
        (output) S1 GRD pre-processing log file
            -c       (option) apply radiometric calibration factor without noise subtraction
            -n       (option) apply radiometric calibration factor with noise subtraction
            -t       (option) include full timestamp YYYYMMDDtHHMMSSin SLC and SLC_par filenames, default YYYYMMDD
    options:
        not documented
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/S1_GRD_preproc', S1_list, MLI_dir, pol, log, options]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def S1_path_number(S1_zipfile, logpath=None, outdir=None, shellscript=None):
    """
    | S1_path_number: Script to determine S1 path (or track) number
    | Copyright 2025 Gamma Remote Sensing, v1.3 3-Feb-2025 uw/cm/oc
    
    Parameters
    ----------
    S1_zipfile:
        (input) S1 zip filename for the TOPS SLC
            
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/S1_path_number', S1_zipfile]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def S1_TOPS_preproc(S1_list, SLC_dir, pol, log, options='-', logpath=None, outdir=None, shellscript=None):
    """
    | Preprocessing of Sentinel-1 TOPS SLC products, extract SLC data and generate SLC_tab
    | Copyright 2023, Gamma Remote Sensing, v2.8 18-Apr-2023 clw/awi/cm
    
    Parameters
    ----------
    S1_list:
        (input) single column text file. Enteries are directories (including path) containing Sentinel-1 TOPS SLC products
    SLC_dir:
        directory for output SLC data files and SLC parameter files
            Note: output file names have the form : 20150119_iw1_hh.slc
    pol:
        SLC polarization to extract (hh,hv,vh,vv)
    log:
        (output) S1 SLC pre-processing log file
            -c          (option) apply radiometric calibration factor without noise subtraction
            -n          (option) apply radiometric calibration factor with noise subtraction
            -s          (option) output is SCOMPLEX format (default: FCOMPLEX)
            -t          (option) include full timestamp YYYYMMDDtHHMMSS in SLC and SLC_par filenames, default YYYYMMDD
            -m MLI_dir  (option) calculate MLI images and store in MLI_dir, enter . for current directory
            -r rlks     (option) number of MLI range looks (default: 10)
            -a azlks    (option) number of MLI azimuth looks (default: 2)
            -b SLC_tab  (option) SLC_tab filename, by default SLC_tab_YYMMDD or SLC_tab_YYYYMMDDtHHMMSS
    options:
        not documented
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/S1_TOPS_preproc', S1_list, SLC_dir, pol, log, options]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def SBI_INT(RSLC1, RSLC1_par, RSLC2, RSLC2_par, sbi, off, sbi_pwr, par_out, norm_sq='-', rlks='-', azlks='-', iwflg='-', cflg='-', logpath=None, outdir=None, shellscript=None):
    """
    | SBI_INT: Script to generate azimuth Split-Beam Interferogram from a coregistered interferometric SLC pair
    | Copyright 2023 Gamma Remote Sensing, v1.4 19-Apr-2023 uw/clw/cm
    
    Parameters
    ----------
    RSLC1:
        (input) master single-look complex image (FCOMPLEX or SCOMPLEX)
    RSLC1_par:
        (input) SLC ISP image parameter file of RSLC1
    RSLC2:
        (input) co-registered slave SLC image (FCOMPLEX or SCOMPLEX)
    RSLC2_par:
        (input) SLC ISP image parameter file of RSLC2
    sbi:
        (output) multi-look split-beam interferogram (FCOMPLEX)
    off:
        (output) ISP offset parameter file for multi-look split-beam interferogram (ascii)
    sbi_pwr:
        (output) multi-look reference backscatter intensity image (FLOAT)
    par_out:
        (output) SLC/MLI ISP image parameter file of sbi_pwr
    norm_sq:
        normalized squint difference parameter (enter - for default: 0.5)
    rlks:
        number of range looks in output split-beam interferogram (enter - for default: 1)
    azlks:
        number of azimuth looks in output split-beam interferogram (enter - for default: 1)
    iwflg:
        inverse weighting flag (enter - for default)
            * 0: do not remove azimuth processing spectral window (default)
            * 1: apply inverse of azimuth compression processing window
    
    cflg:
        flag to indicate if intermediate data (e.g. filtered slc) are deleted (enter - for default)
            * 0: intermediate data are deleted (default)
            * 1: intermediate data are NOT deleted
              file names for band-pass filtered SLC are generated automatically
              by adding the letter  b / f  for the backward / foward looking beam
              
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/SBI_INT', RSLC1, RSLC1_par, RSLC2, RSLC2_par, sbi, off, sbi_pwr, par_out, norm_sq, rlks, azlks, iwflg, cflg]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def ScanSAR_burst_cc_ad(DIFF_tab, MLI1_tab, MLI2R_tab, slope_tab, texture_tab, CC_tab, log, box_min='-', box_max='-', wgt_flag='-', logpath=None, outdir=None, shellscript=None):
    """
    | Estimate interferometric coherence for ScanSAR burst data using cc_ad
    | Copyright 2023, Gamma Remote Sensing, v1.2 17-Apr-2023 cm
    
    Parameters
    ----------
    DIFF_tab:
        (input) 3 column list of the DIFF swaths listed in order from near to far range 
            DIFF_tab line entries:  DIFF  MLI_par  TOPS_par
    MLI1_tab:
        (input) 3 column list of the reference ScanSAR MLI swaths listed in order from near to far range (enter - for none)
            MLI1_tab line entries:   MLI  MLI_par  TOPS_par
    MLI2R_tab:
        (input) 3 column list of ScanSAR MLI swaths listed in order from near to far range, coregistered with MLI1 (enter - for none)
            MLI2R_tab line entries:  MLI  MLI_par  TOPS_par
    slope_tab:
        (input) 1 column list of ScanSAR phase slope swaths listed in order from near to far range (enter - for none)
    texture_tab:
        (input) 1 column list of ScanSAR backscatter texture swaths listed in order from near to far range (enter - for none)
    CC_tab:
        (input/output) 3 column list of the CC swaths listed in order from near to far range
            CC_tab line entries:      CC  MLI_par  TOPS_par
            
            * NOTE: if CC_tab does not exist, it will be created in the current directory.
              The binary file will be named from the differential interferogram name, with the addition of a ".cc" extension.
              The MLI_par and TOPS_par files are copied from MLI1_tab if available, from DIFF_tab otherwise.
    log:
        (output) processing log file
    box_min:
        smallest correlation average box size (enter - for default: 3.0)
    box_max:
        largest correlation average box size  (enter - for default: 9.0)
    wgt_flag:
        weighting function (enter - for default)
            * 0: constant (default)
            * 1: gaussian
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/ScanSAR_burst_cc_ad', DIFF_tab, MLI1_tab, MLI2R_tab, slope_tab, texture_tab, CC_tab, log, box_min, box_max, wgt_flag]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def ScanSAR_burst_cc_wave(DIFF_tab, MLI1_tab, MLI2R_tab, CC_tab, log, bx='-', by='-', wflg='-', logpath=None, outdir=None, shellscript=None):
    """
    | Estimate interferometric coherence for ScanSAR burst data using cc_wave
    | Copyright 2019, Gamma Remote Sensing, v1.2 24-Apr-2019 cm
    
    Parameters
    ----------
    DIFF_tab:
        (input) 3 column list of the DIFF swaths listed in order from near to far range 
            DIFF_tab line entries:  DIFF  MLI_par  TOPS_par
    MLI1_tab:
        (input) 3 column list of the reference ScanSAR MLI swaths listed in order from near to far range (enter - for none)
            MLI1_tab line entries:   MLI  MLI_par  TOPS_par
    MLI2R_tab:
        (input) 3 column list of ScanSAR MLI swaths listed in order from near to far range, coregistered with MLI1 (enter - for none)
            MLI2R_tab line entries:  MLI  MLI_par  TOPS_par
    CC_tab:
        (input/output) 3 column list of the CC swaths listed in order from near to far range
            CC_tab line entries:      CC  MLI_par  TOPS_par
            
            * NOTE: if CC_tab does not exist, it will be created in the current directory.
              The binary file will be named from the differential interferogram name, with the addition of a ".cc" extension.
              The MLI_par and TOPS_par files are copied from MLI1_tab if available, from DIFF_tab otherwise.
    log:
        (output) processing log file
    bx:
        estimation window size in columns (enter - for default: 5.0)
    by:
        estimation window size in lines (enter - for default: 5.0)
    wflg:
        estimation window (enter - for default):
            * 0: rectangular (default)
            * 1: triangular
            * 2: Gaussian
            * 3: normalized vector sum with rectangular window
            * NOTE: This estimator does not use the MLI data, even when specified
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/ScanSAR_burst_cc_wave', DIFF_tab, MLI1_tab, MLI2R_tab, CC_tab, log, bx, by, wflg]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def SLC_copy_WB(SLC_tab, SLC2_dir, logpath=None, outdir=None, shellscript=None):
    """
    | GAMMA_SOFTWARE-20250625/ISP/scripts/SLC_copy_WB
    | Copyright 2011, Gamma Remote Sensing, v1.1 9-Apr-2011 clw
    | Create a new set of SLCs for all beams in a PALSAR WB ScanSAR image
    
    Parameters
    ----------
    SLC_tab:
        (input) two column list of input SLC files and SLC ISP image parameter files (including paths) (text)
    SLC2_dir:
        directory to contain copied segments of the input SLC data and the associated parameter files
            * NOTE: current directory is denoted using .
    
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/SLC_copy_WB', SLC_tab, SLC2_dir]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def TX_SLC_preproc(TSX_list, SLC_dir, log, logpath=None, outdir=None, shellscript=None):
    """
    | Preprocessing of TerraSAR-X TDX1 and TSX1 SLC products using par_TX_SLC
    | Copyright 2023, Gamma Remote Sensing, v1.3 17-Apr-2023 clw
    
    Parameters
    ----------
    TSX_list:
        (input) single column text file with directories (including path)
            containing path to directory containing product XML for IMAGEDATA/\\*.cos files
    SLC_dir:
        directory for output SLC data files and SLC parameter files
    log:
        (output) processing log file
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/TX_SLC_preproc', TSX_list, SLC_dir, log]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def unw_correction_filt(unw_in, unw_out, width, fsize='-', thresh1='-', thresh2='-', iterations='-', cleaning='-', logpath=None, outdir=None, shellscript=None):
    """
    | unw_correction_filt: Phase unwrapping ambiguity error correction relative to spatially filtered phase
    | Copyright 2023 Gamma Remote Sensing, v1.5 18-Apr-2023 uw/cm
    
    Parameters
    ----------
    unw_in:
        (input) unwrapped phase file to correct (float)
    unw_out:
        (output) corrected  unwrapped phase file (float)
    width:
        number of range samples per line
    fsize:
        maximum filter radius in pixels (enter - for default: 5)
    thresh1:
        upper threshold for negative phase differences (enter - for default: -3.0)
    thresh2:
        lower threshold for positive phase differences (enter - for default: 3.0)
    iterations:
        number of iterations to run (enter - for default: 1)
    cleaning:
        cleaning flag indicating if intermediary files are deleted (enter - for default)
            * 0: no
            * 1: yes (default)
              The difference between the unfiltered and spatially filtered phase (using fspf) is used
              to determine an correct phase unwrapping ambiguity errors
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/unw_correction_filt', unw_in, unw_out, width, fsize, thresh1, thresh2, iterations, cleaning]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def unw_correction_poly(unw_in, unw_out, width, poly, flag, max_iter='-', logpath=None, outdir=None, shellscript=None):
    """
    | unw_correction_poly: Phase unwrapping ambiguity error correction for polygon areas
    | Copyright 2023 Gamma Remote Sensing, v1.5 18-Apr-2023 uw/cm
    
    Parameters
    ----------
    unw_in:
        (input) unwrapped phase file to correct (FLOAT)
    unw_out:
        (output) corrected  unwrapped phase file (FLOAT)
    width:
        number of range samples per line
    poly:
        (input) polygon file (text)
    flag:
        ambiguity correction flag (1: add 2PI; -1: subtract 2PI)
    max_iter:
        maximum number of iterations done (enter - for default: 1)
            (iterations are used (a) if the ambiguity to correct is not 2PI but a
            multiple of 2PI and (b) if the ambiguity error is in an area with a
            significant phase slope)
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/unw_correction_poly', unw_in, unw_out, width, poly, flag, max_iter]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def UNWRAP(interf, cc, pwr, unwrap, flag, width, lines, corr_thr='-', pwr_thr='-', r_init='-', az_init='-', r1='-', r2='-', l1='-', l2='-', logpath=None, outdir=None, shellscript=None):
    """
    | UNWRAP: unwrap phase
    | Copyright 2023 Gamma Remote Sensing, v1.4 19-Apr-2023 clw/cm
    
    Parameters
    ----------
    interf:
        interferogram filename  (\\*.int, \\*.flt)
    cc:
        correlation filename (\\*.cc)
    pwr:
        intensity image (\\*.pwr, \\*.mli)
    unwrap:
        unwrap output file (\\*.unw)
    flag:
        unwapping flag file (\\*.flag)
    width:
        interferogram width
    lines:
        number of interferogram lines
    corr_thr:
        threshold for correlation in the unwrapping mask (enter - for default: 0.7)
    pwr_thr:
        intensity threshold for phase unwrapping neutrons, multiples of average (enter - for default: 6.0)
    r_init:
        range seed location in the interferogram (enter - for default: width/2)
    az_init:
        azimuth seed location in the interferogram (enter - for default: nlines/2)
    r1:
        starting range sample offset to unwrap (enter - for default: 0)
    r2:
        ending range sample offset to unwrap (enter - for default: width-1)
    l1:
        starting line offset to unwrap (enter - for default: 0)
    l2:
        ending line offset to unwrap (enter - for default: nlines-1)\n
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/UNWRAP', interf, cc, pwr, unwrap, flag, width, lines, corr_thr, pwr_thr, r_init, az_init, r1, r2, l1, l2]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


def UNWRAP_PAR(interf_par, interf, cc, pwr, unwrap, flag, corr_thr='-', pwr_thr='-', r_init='-', az_init='-', r1='-', r2='-', l1='-', l2='-', logpath=None, outdir=None, shellscript=None):
    """
    | UNWRAP_PAR: unwrap phase using parameters from the ISP interferogram parameter file
    | Copyright 2023 Gamma Remote Sensing, v1.3 19-Apr-2023 clw/cm
    
    Parameters
    ----------
    interf_par:
        interferogram parameter file \\*.off
    interf:
        interferogram filename  (\\*.int, \\*.flt)
    cc:
        correlation filename (\\*.cc)
    pwr:
        intensity image (\\*.pwr, \\*.mli)
    unwrap:
        unwrap output file (\\*.unw)
    flag:
        unwapping flag file (\\*.flag)
    corr_thr:
        threshold for correlation in the unwrapping mask (enter - for default: 0.7)
    pwr_thr:
        intensity threshold for phase unwrapping neutrons, multiples of average (enter - for default: 6.0)
    r_init:
        range seed location in the interferogram (enter - for default: width/2)
    az_init:
        azimuth seed location in the interferogram (enter - for default: nlines/2)
    r1:
        starting range sample offset to unwrap (enter - for default: 0)
    r2:
        ending range sample offset to unwrap (enter - for default: width-1)
    l1:
        starting line offset to unwrap (enter - for default: 0)
    l2:
        ending line offset to unwrap (enter - for default: nlines-1)\n
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    """
    cmd = ['GAMMA_SOFTWARE-20250625/ISP/scripts/UNWRAP_PAR', interf_par, interf, cc, pwr, unwrap, flag, corr_thr, pwr_thr, r_init, az_init, r1, r2, l1, l2]
    process(cmd, logpath=logpath, outdir=outdir, shellscript=shellscript)


