###############################################################################
# universal core routines for processing SAR images with GAMMA

# Copyright (c) 2014-2021, the pyroSAR Developers.

# This file is part of the pyroSAR Project. It is subject to the
# license terms in the LICENSE.txt file found in the top-level
# directory of this distribution and at
# https://github.com/johntruckenbrodt/pyroSAR/blob/master/LICENSE.txt.
# No part of the pyroSAR project, including this file, may be
# copied, modified, propagated, or distributed except according
# to the terms contained in the LICENSE.txt file.
###############################################################################

"""
This module is intended as a set of generalized processing routines for modularized GAMMA work flows.
The function parametrization is intended to be applicable to any kind of situation and input data set.
Thus, instead of choosing a specific parametrization for the data at hand,
core parameters are iterated over a set of values in order to find the one best suited for the task.
The approach of the single routines is likely to still have drawbacks and might fail in certain situations.
Testing and suggestions on improvements are very welcome.
"""
import os
import re
import shutil
import zipfile as zf
from datetime import datetime
from urllib.error import URLError

from spatialist import haversine
from spatialist.ancillary import union, finder

from ..S1 import OSV
from ..drivers import ID, identify, identify_many
from . import ISPPar, Namespace, par2hdr
from ..ancillary import multilook_factors, hasarg, groupby
from pyroSAR.examine import ExamineSnap
from .auxil import do_execute

import logging

log = logging.getLogger(__name__)

try:
    from .api import diff, disp, isp, lat
except ImportError:
    pass


def calibrate(id, directory, replace=False, logpath=None, outdir=None, shellscript=None):
    """
    radiometric calibration of SAR scenes
    
    Parameters
    ----------
    id: ~pyroSAR.drivers.ID
        an SAR scene object of type pyroSAR.ID or any subclass
    directory: str
        the directory to search for Gamma calibration candidates
    replace: bool
        replace the input images by the new files? If True, the input images will be deleted.
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format

    Returns
    -------

    """
    cname = type(id).__name__
    if cname == 'CEOS_PSR':
        for image in id.getGammaImages(directory):
            if image.endswith('_slc'):
                isp.radcal_SLC(SLC=image,
                               SLC_par=image + '.par',
                               CSLC=image + '_cal',
                               CSLC_par=image + '_cal.par',
                               K_dB=id.meta['k_dB'],
                               logpath=logpath,
                               outdir=outdir,
                               shellscript=shellscript)
                par2hdr(image + '_cal.par', image + '_cal.hdr')
    
    elif cname == 'EORC_PSR':
        for image in id.getGammaImages(directory):
            pol = re.search('[HV]{2}', os.path.basename(image)).group(0)
            if image.endswith('_mli'):
                isp.radcal_MLI(MLI=image,
                               MLI_par=image + '.par',
                               OFF_par='-',
                               CMLI=image + '_cal',
                               antenna='-',
                               rloss_flag=0,
                               ant_flag=0,
                               refarea_flag=1,
                               sc_dB=0,
                               K_dB=id.meta['k_dB'],
                               pix_area=image + '_cal_pix_ell',
                               logpath=logpath,
                               outdir=outdir,
                               shellscript=shellscript)
                par2hdr(image + '.par', image + '_cal.hdr')
                par2hdr(image + '.par', image + '_cal_pix_ell' + '.hdr')
                # rename parameter file 
                os.rename(image + '.par', image + '_cal.par')
    
    elif cname == 'ESA':
        k_db = {'ASAR': 55., 'ERS1': 58.24, 'ERS2': 59.75}[id.sensor]
        inc_ref = 90. if id.sensor == 'ASAR' else 23.
        candidates = [x for x in id.getGammaImages(directory) if re.search('_pri$', x)]
        for image in candidates:
            out = image.replace('pri', 'grd')
            isp.radcal_PRI(PRI=image,
                           PRI_par=image + '.par',
                           GRD=out,
                           GRD_par=out + '.par',
                           K_dB=k_db,
                           inc_ref=inc_ref,
                           logpath=logpath,
                           outdir=outdir,
                           shellscript=shellscript)
            par2hdr(out + '.par', out + '.hdr')
            if replace:
                for item in [image, image + '.par', image + '.hdr']:
                    if os.path.isfile(item):
                        os.remove(item)
    
    elif cname == 'SAFE':
        log.info('calibration already performed during import')
    
    else:
        raise NotImplementedError('calibration for class {} is not implemented yet'.format(cname))


def convert2gamma(id, directory, S1_tnr=True, S1_bnr=True,
                  basename_extensions=None, exist_ok=False,
                  logpath=None, outdir=None, shellscript=None):
    """
    general function for converting SAR images to GAMMA format

    Parameters
    ----------
    id: ~pyroSAR.drivers.ID
        an SAR scene object of type pyroSAR.ID or any subclass
    directory: str
        the output directory for the converted images
    S1_tnr: bool
        only Sentinel-1: should thermal noise removal be applied to the image?
    S1_bnr: bool
        only Sentinel-1 GRD: should border noise removal be applied to the image?
        This is available since version 20191203, for older versions this argument is ignored.
    basename_extensions: list of str
        names of additional parameters to append to the basename, e.g. ['orbitNumber_rel']
    exist_ok: bool
        allow existing output files and do not create new ones?
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format

    Returns
    -------

    """
    
    if not isinstance(id, ID):
        raise IOError('id must be of type pyroSAR.ID')
    
    if id.compression is not None:
        raise RuntimeError('scene is not yet unpacked')
    
    if not os.path.isdir(directory):
        os.makedirs(directory)
    
    cname = type(id).__name__
    
    if cname == 'CEOS_ERS':
        if id.sensor in ['ERS1', 'ERS2']:
            if id.product == 'SLC' \
                    and id.meta['proc_system'] in ['PGS-ERS', 'VMP-ERS', 'SPF-ERS']:
                outname_base = id.outname_base(extensions=basename_extensions)
                outname_base = '{}_{}_{}'.format(outname_base,
                                                 id.polarizations[0],
                                                 id.product.lower())
                outname = os.path.join(directory, outname_base)
                if not os.path.isfile(outname):
                    lea = id.findfiles('LEA_01.001')[0]
                    dat = id.findfiles('DAT_01.001')[0]
                    title = re.sub(r'\.PS$', '', os.path.basename(id.file))
                    
                    pars = {'CEOS_SAR_leader': lea,
                            'SLC_par': outname + '.par',
                            'CEOS_DAT': dat,
                            'SLC': outname,
                            'inlist': [title],
                            'logpath': logpath,
                            'outdir': outdir,
                            'shellscript': shellscript}
                    
                    if do_execute(pars, ['SLC', 'SLC_par'], exist_ok):
                        isp.par_ESA_ERS(**pars)
                        par2hdr(outname + '.par', outname + '.hdr')
                else:
                    log.info('scene already converted')
            else:
                raise NotImplementedError('ERS {} product of {} processor in CEOS format not implemented yet'
                                          .format(id.product, id.meta['proc_system']))
        else:
            raise NotImplementedError('sensor {} in CEOS format not implemented yet'.format(id.sensor))
    
    elif cname == 'CEOS_PSR':
        images = id.findfiles('^IMG-')
        if id.product == '1.0':
            raise RuntimeError('PALSAR level 1.0 products are not supported')
        for image in images:
            polarization = re.search('[HV]{2}', os.path.basename(image)).group(0)
            outname_base = id.outname_base(extensions=basename_extensions)
            
            pars = {'CEOS_leader': id.file,
                    'CEOS_data': image,
                    'logpath': logpath,
                    'outdir': outdir,
                    'shellscript': shellscript}
            
            if id.product == '1.1':
                outname_base = '{}_{}_slc'.format(outname_base, polarization)
                outname = os.path.join(directory, outname_base)
                
                pars['SLC'] = outname
                pars['SLC_par'] = outname + '.par'
                
                if do_execute(pars, ['SLC', 'SLC_par'], exist_ok):
                    isp.par_EORC_PALSAR(**pars)
                    par2hdr(outname + '.par', outname + '.hdr')
            else:
                outname_base = '{}_{}_mli_geo'.format(outname_base, polarization)
                outname = os.path.join(directory, outname_base)
                
                pars['MLI'] = outname
                pars['MLI_par'] = outname + '.par'
                pars['DEM_par'] = outname + '_dem.par'
                
                if do_execute(pars, ['MLI', 'MLI_par', 'DEM_par'], exist_ok):
                    diff.par_EORC_PALSAR_geo(**pars)
                    par2hdr(outname + '.par', outname + '.hdr')
    
    elif cname == 'EORC_PSR':
        images = id.findfiles('^sar.')
        facter_m = id.findfiles('facter_m.dat')
        led = id.findfiles('LED-ALOS2')
        
        for image in images:
            polarization = re.search('[HV]{2}', os.path.basename(image)).group(0)
            outname_base = id.outname_base(extensions=basename_extensions)
            outname_base = '{}_{}'.format(outname_base, polarization)
            outname = os.path.join(directory, outname_base)
            
            pars = {'facter_m': facter_m,
                    'CEOS_leader': led,
                    'SLC_par': outname + '_mli.par',
                    'pol': polarization,
                    'pls_mode': 2,
                    'KC_data': image,
                    'pwr': outname + '_mli',
                    'logpath': logpath,
                    'outdir': outdir,
                    'shellscript': shellscript}
            
            if do_execute(pars, ['pwr', 'SLC_par'], exist_ok):
                isp.par_KC_PALSAR_slr(**pars)
                par2hdr(outname + '_mli.par', outname + '_mli.hdr')
    
    elif cname == 'ESA':
        """
        the command par_ASAR also accepts a K_dB argument for calibration
        in which case the resulting image names will carry the suffix grd;
        this is not implemented here but instead in function calibrate
        """
        outname = os.path.join(directory, id.outname_base(extensions=basename_extensions))
        if not id.is_processed(directory):
            
            isp.par_ASAR(ASAR_ERS_file=os.path.basename(id.file),
                         output_name=outname,
                         outdir=os.path.dirname(id.file),
                         logpath=logpath,
                         shellscript=shellscript)
            
            os.remove(outname + '.hdr')
            for item in finder(directory, [os.path.basename(outname)], regex=True):
                ext = '.par' if item.endswith('.par') else ''
                outname_base = os.path.basename(item) \
                    .strip(ext) \
                    .replace('.', '_') \
                    .replace('PRI', 'pri') \
                    .replace('SLC', 'slc')
                outname = os.path.join(directory, outname_base + ext)
                os.rename(item, outname)
                if outname.endswith('.par'):
                    par2hdr(outname, outname.replace('.par', '.hdr'))
        else:
            if not exist_ok:
                raise IOError('scene already processed')
    
    elif cname == 'SAFE':
        if id.product == 'OCN':
            raise IOError('Sentinel-1 OCN products are not supported')
        if id.meta['category'] == 'A':
            raise IOError('Sentinel-1 annotation-only products are not supported')
        
        for xml_ann in finder(os.path.join(id.scene, 'annotation'), [id.pattern_ds], regex=True):
            base = os.path.basename(xml_ann)
            match = re.compile(id.pattern_ds).match(base)
            
            tiff = os.path.join(id.scene, 'measurement', base.replace('.xml', '.tiff'))
            xml_cal = os.path.join(id.scene, 'annotation', 'calibration', 'calibration-' + base)
            
            product = match.group('product')
            
            # specify noise calibration file
            # L1 GRD product: thermal noise already subtracted, specify xml_noise to add back thermal noise
            # SLC products: specify noise file to remove noise
            # xml_noise = '-': noise file not specified
            if (S1_tnr and product == 'slc') or (not S1_tnr and product == 'grd'):
                xml_noise = os.path.join(id.scene, 'annotation', 'calibration', 'noise-' + base)
            else:
                xml_noise = '-'
            
            fields = (id.outname_base(extensions=basename_extensions),
                      match.group('pol').upper(),
                      product)
            outname = os.path.join(directory, '_'.join(fields))
            
            pars = {'GeoTIFF': tiff,
                    'annotation_XML': xml_ann,
                    'calibration_XML': xml_cal,
                    'noise_XML': xml_noise,
                    'logpath': logpath,
                    'shellscript': shellscript,
                    'outdir': outdir}
            
            if product == 'slc':
                swath = match.group('swath').upper()
                old = '{:_<{length}}'.format(id.acquisition_mode, length=len(swath))
                outname = outname.replace(old, swath)
                pars['SLC'] = outname
                pars['SLC_par'] = outname + '.par'
                pars['TOPS_par'] = outname + '.tops_par'
                
                if do_execute(pars, ['SLC', 'SLC_par', 'TOPS_par'], exist_ok):
                    isp.par_S1_SLC(**pars)
                    par2hdr(outname + '.par', outname + '.hdr')
            else:
                if hasarg(isp.par_S1_GRD, 'edge_flag'):
                    if S1_bnr:
                        pars['edge_flag'] = 2
                    else:
                        pars['edge_flag'] = 0
                pars['MLI'] = outname
                pars['MLI_par'] = outname + '.par'
                if do_execute(pars, ['MLI', 'MLI_par'], exist_ok):
                    isp.par_S1_GRD(**pars)
                    par2hdr(outname + '.par', outname + '.hdr')
    
    elif cname == 'TSX':
        images = id.findfiles(id.pattern_ds)
        pattern = re.compile(id.pattern_ds)
        for image in images:
            pol = pattern.match(os.path.basename(image)).group('pol')
            outname_base = id.outname_base(extensions=basename_extensions)
            outname = os.path.join(directory, outname_base + '_' + pol)
            
            pars = {'annotation_XML': id.file,
                    'pol': pol,
                    'logpath': logpath,
                    'shellscript': shellscript,
                    'outdir': outdir}
            
            if id.product == 'SSC':
                outname += '_slc'
                pars['COSAR'] = image
                pars['SLC_par'] = outname + '.par'
                pars['SLC'] = outname
                if do_execute(pars, ['SLC', 'SLC_par'], exist_ok):
                    isp.par_TX_SLC(**pars)
                    par2hdr(outname + '.par', outname + '.hdr')
            
            elif id.product == 'MGD':
                outname += '_mli'
                pars['GeoTIFF'] = image
                pars['GRD_par'] = outname + '.par'
                pars['GRD'] = outname
                if do_execute(pars, ['GRD', 'GRD_par'], exist_ok):
                    isp.par_TX_GRD(**pars)
                    par2hdr(outname + '.par', outname + '.hdr')
            
            elif id.product in ['GEC', 'EEC']:
                outname += '_mli_geo'
                pars['GeoTIFF'] = image
                pars['MLI_par'] = outname + '.par'
                pars['DEM_par'] = outname + '_dem.par'
                pars['GEO'] = outname
                if do_execute(pars, ['GEO', 'MLI_par', 'DEM_par'], exist_ok):
                    diff.par_TX_geo(**pars)
                    par2hdr(outname + '.par', outname + '.hdr')
            else:
                raise RuntimeError('unknown product: {}'.format(id.product))
    
    else:
        raise NotImplementedError('conversion for class {} is not implemented yet'.format(cname))


def correctOSV(id, directory=None, osvdir=None, osvType='POE', timeout=20, logpath=None, outdir=None, shellscript=None):
    """
    correct GAMMA parameter files with orbit state vector information from dedicated OSV files;
    OSV files are downloaded automatically to either the defined `osvdir` or a sub-directory `osv` of the scene directory
    
    Parameters
    ----------
    id: ~pyroSAR.drivers.ID
        the scene to be corrected
    directory: str or None
        a directory to be scanned for files associated with the scene.
        If None, the scene is expected to be an unpacked directory in which the files are searched.
    osvdir: str
        the directory of OSV files; subdirectories POEORB and RESORB are created automatically
    osvType: {'POE', 'RES'}
        the OSV type to be used
    timeout: int or tuple or None
        the timeout in seconds for downloading OSV files as provided to :func:`requests.get`
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format
    
    Returns
    -------
    
    Examples
    --------
    
    >>> from pyroSAR import identify
    >>> from pyroSAR.gamma import correctOSV, convert2gamma
    >>> filename = 'S1A_IW_GRDH_1SDV_20150222T170750_20150222T170815_004739_005DD8_3768.zip'
    # identify the SAR scene
    >>> scene = identify(filename)
    # unpack the zipped scene to an arbitrary directory
    >>> scene.unpack('/home/test')
    >>> print(scene.scene)
    /home/test/S1A_IW_GRDH_1SDV_20150222T170750_20150222T170815_004739_005DD8_3768.SAFE
    # convert the unpacked scene to GAMMA format
    >>> convert2gamma(id=scene, directory=scene.scene)
    # correct the OSV information of the converted GAMMA images
    >>> correctOSV(id=scene, osvdir='/home/test/osv')
    
    See Also
    --------
    :meth:`pyroSAR.drivers.SAFE.getOSV`
    :class:`pyroSAR.S1.OSV`
    """
    
    if not isinstance(id, ID):
        raise IOError('id must be of type pyroSAR.ID')
    
    if id.sensor not in ['S1A', 'S1B']:
        raise IOError('this method is currently only available for Sentinel-1. Please stay tuned...')
    
    if not os.path.isdir(logpath):
        os.makedirs(logpath)
    
    if osvdir is None:
        try:
            auxdatapath = ExamineSnap().auxdatapath
        except AttributeError:
            auxdatapath = os.path.join(os.path.expanduser('~'), '.snap', 'auxdata')
        osvdir = os.path.join(auxdatapath, 'Orbits', 'Sentinel-1')
    try:
        id.getOSV(osvdir, osvType, timeout=timeout)
    except URLError:
        log.warning('..no internet access')
    
    target = directory if directory is not None else id.scene
    images = id.getGammaImages(target)
    # read parameter file entries into object
    with ISPPar(images[0] + '.par') as par:
        # extract acquisition time stamp
        timestamp = datetime.strptime(par.date, '%Y-%m-%dT%H:%M:%S.%f').strftime('%Y%m%dT%H%M%S')
    
    # find an OSV file matching the time stamp and defined OSV type(s)
    with OSV(osvdir, timeout=timeout) as osv:
        osvfile = osv.match(sensor=id.sensor, timestamp=timestamp, osvtype=osvType)
    if not osvfile:
        raise RuntimeError('no Orbit State Vector file found')
    
    if osvfile.endswith('.zip'):
        osvdir = os.path.join(id.scene, 'osv')
        with zf.ZipFile(osvfile) as zip:
            zip.extractall(path=osvdir)
        osvfile = os.path.join(osvdir, os.path.basename(osvfile).replace('.zip', ''))
    
    # update the GAMMA parameter file with the selected orbit state vectors
    log.info('correcting state vectors with file {}'.format(osvfile))
    for image in images:
        isp.S1_OPOD_vec(SLC_par=image + '.par',
                        OPOD=osvfile,
                        logpath=logpath,
                        outdir=outdir,
                        shellscript=shellscript)


def geocode(scene, dem, tmpdir, outdir, targetres, scaling='linear', func_geoback=1,
            nodata=(0, -99), osvdir=None, allow_RES_OSV=False,
            cleanup=True, export_extra=None, basename_extensions=None,
            removeS1BorderNoiseMethod='gamma', refine_lut=False):
    """
    general function for radiometric terrain correction (RTC) and geocoding of SAR backscatter images with GAMMA.
    Applies the RTC method by :cite:`Small2011` to retrieve gamma nought RTC backscatter.
    
    Parameters
    ----------
    scene: str or ~pyroSAR.drivers.ID or list
        the SAR scene(s) to be processed
    dem: str
        the reference DEM in GAMMA format
    tmpdir: str
        a temporary directory for writing intermediate files
    outdir: str
        the directory for the final GeoTIFF output files
    targetres: int
        the target resolution in meters
    scaling: {'linear', 'db'} or list
        the value scaling of the backscatter values; either 'linear', 'db' or a list of both, i.e. ['linear', 'db']
    func_geoback: {0, 1, 2, 3, 4, 5, 6, 7}
        backward geocoding interpolation mode (see GAMMA command geocode_back)
         - 0: nearest-neighbor
         - 1: bicubic spline (default)
         - 2: bicubic-spline, interpolate log(data)
         - 3: bicubic-spline, interpolate sqrt(data)
         - 4: B-spline interpolation (default B-spline degree: 5)
         - 5: B-spline interpolation sqrt(x) (default B-spline degree: 5)
         - 6: Lanczos interpolation (default Lanczos function order: 5)
         - 7: Lanczos interpolation sqrt(x) (default Lanczos function order: 5)

        NOTE: log and sqrt interpolation modes should only be used with non-negative data!
        
        NOTE: Gamma reccomendation for MLI data: "The interpolation should be performed on
        the square root of the data. A mid-order (3 to 5) B-spline interpolation is recommended."
    nodata: tuple
        the nodata values for the output files; defined as a tuple with two values, the first for linear,
        the second for logarithmic scaling
    osvdir: str
        a directory for Orbit State Vector files;
        this is currently only used by for Sentinel-1 where two subdirectories POEORB and RESORB are created;
        if set to None, a subdirectory OSV is created in the directory of the unpacked scene.
    allow_RES_OSV: bool
        also allow the less accurate RES orbit files to be used?
        Otherwise the function will raise an error if no POE file exists.
    cleanup: bool
        should all files written to the temporary directory during function execution be deleted after processing?
    export_extra: list of str or None
        a list of image file IDs to be exported to outdir
         - format is GeoTIFF if the file is geocoded and ENVI otherwise. Non-geocoded images can be converted via GAMMA
           command data2tiff yet the output was found impossible to read with GIS software
         - scaling of SAR image products is applied as defined by parameter `scaling`
         - see Notes for ID options
    basename_extensions: list of str or None
        names of additional parameters to append to the basename, e.g. ['orbitNumber_rel']
    removeS1BorderNoiseMethod: str or None
        the S1 GRD border noise removal method to be applied, See :func:`pyroSAR.S1.removeGRDBorderNoise` for details; one of the following:
         - 'ESA': the pure implementation as described by ESA
         - 'pyroSAR': the ESA method plus the custom pyroSAR refinement
         - 'gamma': the GAMMA implementation of :cite:`Ali2018`
         - None: do not remove border noise
    refine_lut: bool
        should the LUT for geocoding be refined using pixel area normalization?
    
    Returns
    -------
    
    Note
    ----
    | intermediate output files
    | DEM products are named <scene identifier>_<ID>, e.g. `S1A__IW___A_20141012T162337_inc_geo`
    | SAR products will additionally contain the polarization, e.g. `S1A__IW___A_20141012T162337_VV_grd_mli`
    | IDs in brackets are only written if selected by `export_extra`
    
    - images in range-Doppler geometry
    
      * **grd**: the ground range detected SAR intensity image
      * **grd_mli**: the multi-looked grd image with approximated target resolution
      * (**pix_ellip_sigma0**): ellipsoid-based pixel area
      * (**pix_area_sigma0**): illuminated area as obtained from integrating DEM-facets in sigma projection (command pixel_area)
      * (**pix_area_gamma0**): illuminated area as obtained from integrating DEM-facets in gamma projection (command pixel_area)
      * **pix_ratio**: pixel area normalization factor (pix_ellip_sigma0 / pix_area_gamma0)
      * **grd_mli_gamma0-rtc**: the terrain-corrected gamma0 backscatter (grd_mli * pix_ratio)
      * (**gs_ratio**): gamma-sigma ratio (pix_gamma0 / pix_sigma0)
    
    - images in map geometry
    
      * **dem_seg_geo**: dem subsetted to the extent of the intersect between input DEM and SAR image
      * (**u_geo**): zenith angle of surface normal vector n (angle between z and n)
      * (**v_geo**): orientation angle of n (between x and projection of n in xy plane)
      * **inc_geo**: local incidence angle (between surface normal and look vector)
      * (**psi_geo**): projection angle (between surface normal and image plane normal)
      * **ls_map_geo**: layover and shadow map (in map projection)
      * (**sim_sar_geo**): simulated SAR backscatter image
      * (**pix_ellip_sigma0_geo**): ellipsoid-based pixel area
      * (**pix_area_sigma0_geo**): illuminated area as obtained from integrating DEM-facets in sigma projection (command pixel_area)
      * (**pix_area_gamma0_geo**): illuminated area as obtained from integrating DEM-facets in gamma projection (command pixel_area)
      * (**pix_ratio_geo**): pixel area normalization factor (pix_ellip_sigma0 / pix_area_gamma0)
      * (**gs_ratio_geo**): gamma-sigma ratio (pix_gamma0 / pix_sigma0)
    
    - additional files
    
      * **lut_init**: initial geocoding lookup table
    
    - files specific to lookup table refinement
    
      * **lut_fine**: refined geocoding lookup table
      * **diffpar**: ISP offset/interferogram parameter file
      * **offs**: offset estimates (fcomplex)
      * **coffs**: culled range and azimuth offset estimates (fcomplex)
      * **coffsets**: culled offset estimates and cross correlation values (text format)
      * **ccp**: cross-correlation of each patch (0.0->1.0) (float)
    
    Examples
    --------
    geocode a Sentinel-1 scene and export the local incidence angle map with it
    
    >>> from pyroSAR.gamma import geocode
    >>> filename = 'S1A_IW_GRDH_1SDV_20180829T170656_20180829T170721_023464_028DE0_F7BD.zip'
    >>> geocode(scene=filename, dem='demfile', outdir='outdir', targetres=20, scaling='db',
    >>>         export_extra=['dem_seg_geo', 'inc_geo', 'ls_map_geo'])
    
    .. figure:: figures/gamma_geocode.svg
        :align: center
        
        Workflow diagram for function geocode for processing a Sentinel-1 Ground Range
        Detected (GRD) scene to radiometrically terrain corrected (RTC) gamma nought backscatter.
    
    References
    ----------
    .. bibliography:: references.bib
        :style: plain
    """
    
    scenes = scene if isinstance(scene, list) else [scene]
    if len(scenes) > 2:
        raise RuntimeError("currently only one or two scenes can be passed via argument 'scene'")
    scenes = identify_many(scenes)
    ref = scenes[0]
    
    if ref.sensor not in ['S1A', 'S1B', 'PALSAR-2']:
        raise RuntimeError(
            'this function currently only supports Sentinel-1 and PALSAR-2 Path data. Please stay tuned...')
    
    if export_extra is not None and not isinstance(export_extra, list):
        raise TypeError("parameter 'export_extra' must either be None or a list")
    
    tmpdir = os.path.join(tmpdir, ref.outname_base(extensions=basename_extensions))
    
    for dir in [tmpdir, outdir]:
        os.makedirs(dir, exist_ok=True)
    
    if ref.is_processed(outdir):
        log.info('scene {} already processed'.format(ref.outname_base(extensions=basename_extensions)))
        return
    
    shellscript = os.path.join(tmpdir, ref.outname_base(extensions=basename_extensions) + '_commands.sh')
    
    scaling = [scaling] if isinstance(scaling, str) else scaling if isinstance(scaling, list) else []
    scaling = union(scaling, ['db', 'linear'])
    if len(scaling) == 0:
        raise IOError('wrong input type for parameter scaling')
    
    for scene in scenes:
        if scene.compression is not None:
            log.info('unpacking scene')
            try:
                scene.unpack(tmpdir)
            except RuntimeError:
                log.info('scene was attempted to be processed before, exiting')
                return
        else:
            scene.scene = os.path.join(tmpdir, os.path.basename(scene.file))
            os.makedirs(scene.scene)
    
    path_log = os.path.join(tmpdir, 'logfiles')
    if not os.path.isdir(path_log):
        os.makedirs(path_log)
    
    for scene in scenes:
        if scene.sensor in ['S1A', 'S1B'] and removeS1BorderNoiseMethod in ['ESA', 'pyroSAR']:
            log.info('removing border noise')
            scene.removeGRDBorderNoise(method=removeS1BorderNoiseMethod)
    
    log.info('converting scene to GAMMA format')
    gamma_bnr = True if removeS1BorderNoiseMethod == 'gamma' else False
    for scene in scenes:
        convert2gamma(scene, directory=tmpdir, logpath=path_log, outdir=tmpdir,
                      basename_extensions=basename_extensions, shellscript=shellscript,
                      S1_bnr=gamma_bnr)
    
    for scene in scenes:
        if scene.sensor in ['S1A', 'S1B']:
            log.info('updating orbit state vectors')
            if allow_RES_OSV:
                osvtype = ['POE', 'RES']
            else:
                osvtype = 'POE'
            try:
                correctOSV(id=scene, directory=tmpdir, osvdir=osvdir, osvType=osvtype,
                           logpath=path_log, outdir=tmpdir, shellscript=shellscript)
            except RuntimeError:
                log.warning('orbit state vector correction failed for scene {}'.format(scene.scene))
                return
    
    log.info('calibrating')
    for scene in scenes:
        calibrate(id=scene, directory=tmpdir, logpath=path_log, outdir=tmpdir, shellscript=shellscript)
    
    images = []
    for scene in scenes:
        img = [x for x in scene.getGammaImages(tmpdir) if
               x.endswith('_grd') or x.endswith('_slc_cal') or x.endswith('_mli_cal')]
        images.extend(img)
    
    if len(scenes) > 1:
        images_new = []
        groups = groupby(images, 'polarization')
        for group in groups:
            out = group[0] + '_cat'
            out_par = out + '.par'
            all_exist = all([os.path.isfile(x) for x in [out, out_par]])
            if not all_exist:
                log.info('mosaicing scenes')
                isp.MLI_cat(MLI_1=group[0],
                            MLI1_par=group[0] + '.par',
                            MLI_2=group[1],
                            MLI2_par=group[1] + '.par',
                            MLI_3=out,
                            MLI3_par=out_par,
                            logpath=path_log, outdir=tmpdir, shellscript=shellscript)
                par2hdr(out_par, out + '.hdr')
            images_new.append(out)
        images = images_new
    
    products = list(images)
    
    if scene.sensor in ['S1A', 'S1B']:
        log.info('multilooking')
        for image in images:
            multilook(infile=image, outfile=image + '_mli', targetres=targetres,
                      logpath=path_log, outdir=tmpdir, shellscript=shellscript)
        
        images = [x + '_mli' for x in images]
        products.extend(images)
    
    master = images[0]
    
    # create output names for files to be written
    # appreciated files will be written
    n = Namespace(tmpdir, scene.outname_base(extensions=basename_extensions))
    n.appreciate(['dem_seg_geo', 'lut_init', 'inc_geo', 'ls_map_geo'])
    
    if export_extra is not None:
        n.appreciate(export_extra)
        pix = ['pix_area_sigma0', 'pix_area_gamma0', 'pix_ratio', 'gs_ratio', 'pix_ellip_sigma0']
        pix_geo = []
        for item in pix:
            if item + '_geo' in export_extra:
                pix_geo.append(item + '_geo')
                n.appreciate([item])
    
    if refine_lut:
        n.appreciate(['pix_area_sigma0'])
    
    ovs_lat, ovs_lon = ovs(dem + '.par', targetres)
    
    master_par = ISPPar(master + '.par')
    
    gc_map_args = {'DEM_par': dem + '.par',
                   'DEM': dem,
                   'DEM_seg_par': n.dem_seg_geo + '.par',
                   'DEM_seg': n.dem_seg_geo,
                   'lookup_table': n.lut_init,
                   'lat_ovr': ovs_lat,
                   'lon_ovr': ovs_lon,
                   'sim_sar': n.sim_sar_geo,
                   'u': n.u_geo,
                   'v': n.v_geo,
                   'inc': n.inc_geo,
                   'psi': n.psi_geo,
                   'pix': n.pix_geo,
                   'ls_map': n.ls_map_geo,
                   'frame': 8,
                   'ls_mode': 2,
                   'logpath': path_log,
                   'shellscript': shellscript,
                   'outdir': tmpdir}
    
    log.info('creating DEM products')
    if master_par.image_geometry == 'GROUND_RANGE':
        gc_map_args.update({'GRD_par': master + '.par'})
        diff.gc_map_grd(**gc_map_args)
    else:
        gc_map_args.update({'MLI_par': master + '.par',
                            'OFF_par': '-'})
        diff.gc_map(**gc_map_args)
    
    for item in ['dem_seg_geo', 'sim_sar_geo', 'u_geo', 'v_geo',
                 'psi_geo', 'pix_geo', 'inc_geo', 'ls_map_geo']:
        if n.isappreciated(item):
            mods = {'data_type': 1} if item == 'ls_map_geo' else None
            par2hdr(n.dem_seg_geo + '.par', n.get(item) + '.hdr', mods)
    
    sim_width = ISPPar(n.dem_seg_geo + '.par').width
    
    ######################################################################
    # RTC reference area computation #####################################
    ######################################################################
    
    pixel_area_wrap(master=master, namespace=n, lut=n.lut_init,
                    path_log=path_log, outdir=tmpdir, shellscript=shellscript)
    
    ######################################################################
    # lookup table Refinement ############################################
    ######################################################################
    lut_final = n.lut_init
    for image in images:
        if refine_lut:
            log.info('refining LUT of {}'.format(image))
            # Refinement of geocoding lookup table
            diff.create_diff_par(PAR_1=image + '.par',
                                 PAR_2='-',
                                 DIFF_par=image + '_diff.par',
                                 PAR_type=1,
                                 iflg=0,
                                 logpath=path_log,
                                 outdir=tmpdir,
                                 shellscript=shellscript)
            # Refinement Lookuptable
            # for "shift" data offset window size enlarged twice to 512 and 256, for data without shift 256 128
            diff.offset_pwrm(MLI_1=n.pix_area_sigma0,
                             MLI_2=image,
                             DIFF_par=image + '_diff.par',
                             offs=image + '_offs',
                             ccp=image + '_ccp',
                             rwin=512,
                             azwin=256,
                             offsets=image + '_offsets.txt',
                             n_ovr=2,
                             nr=64,
                             naz=32,
                             thres=0.2,
                             logpath=path_log,
                             outdir=tmpdir,
                             shellscript=shellscript)
            # par2hdr(master + '.par', master + '_offs' + '.hdr')
            diff.offset_fitm(offs=image + '_offs',
                             ccp=image + '_ccp',
                             DIFF_par=image + '_diff.par',
                             coffs=image + '_coffs',
                             coffsets=image + '_coffsets',
                             thres=0.2,
                             npoly=4,
                             logpath=path_log,
                             outdir=tmpdir,
                             shellscript=shellscript)
            # Updating of the look-up table
            diff.gc_map_fine(gc_in=lut_final,
                             width=sim_width,
                             DIFF_par=image + '_diff.par',
                             gc_out=lut_final + '.fine',
                             ref_flg=1,
                             logpath=path_log,
                             outdir=tmpdir,
                             shellscript=shellscript)
            # Reproduce pixel area estimate
            pixel_area_wrap(master=image, namespace=n, lut=lut_final + '.fine',
                            path_log=path_log, outdir=tmpdir, shellscript=shellscript)
            lut_final = lut_final + '.fine'
        ######################################################################
        # radiometric terrain correction and backward geocoding ##############
        ######################################################################
        lat.product(data_1=image,
                    data_2=n.pix_ratio,
                    product=image + '_gamma0-rtc',
                    width=master_par.range_samples,
                    bx=1,
                    by=1,
                    logpath=path_log,
                    outdir=tmpdir,
                    shellscript=shellscript)
        par2hdr(master + '.par', image + '_gamma0-rtc.hdr')
        diff.geocode_back(data_in=image + '_gamma0-rtc',
                          width_in=master_par.range_samples,
                          lookup_table=lut_final,
                          data_out=image + '_gamma0-rtc_geo',
                          width_out=sim_width,
                          interp_mode=func_geoback,
                          logpath=path_log,
                          outdir=tmpdir,
                          shellscript=shellscript)
        par2hdr(n.dem_seg_geo + '.par', image + '_gamma0-rtc_geo.hdr')
        products.extend([image + '_gamma0-rtc', image + '_gamma0-rtc_geo'])
    ######################################################################
    # log scaling and image export #######################################
    ######################################################################
    log.info('conversion to (dB and) GeoTIFF')
    
    def exporter(data_in, outdir, nodata, scale='linear', dtype=2):
        if scale == 'db':
            if re.search('_geo', os.path.basename(data_in)):
                width = sim_width
                refpar = n.dem_seg_geo + '.par'
            else:
                width = master_par.range_samples
                refpar = master + '.par'
            lat.linear_to_dB(data_in=data_in,
                             data_out=data_in + '_db',
                             width=width,
                             inverse_flag=0,
                             null_value=nodata,
                             logpath=path_log,
                             outdir=tmpdir,
                             shellscript=shellscript)
            par2hdr(refpar, data_in + '_db.hdr')
            data_in += '_db'
        if re.search('_geo', os.path.basename(data_in)):
            outfile = os.path.join(outdir, os.path.basename(data_in) + '.tif')
            disp.data2geotiff(DEM_par=n.dem_seg_geo + '.par',
                              data=data_in,
                              type=dtype,
                              GeoTIFF=outfile,
                              no_data=nodata,
                              logpath=path_log,
                              outdir=tmpdir,
                              shellscript=shellscript)
        
        else:
            outfile = os.path.join(outdir, os.path.basename(data_in))
            shutil.copyfile(data_in, outfile)
            shutil.copyfile(data_in + '.hdr', outfile + '.hdr')
    
    for image in images:
        for scale in scaling:
            exporter(data_in=image + '_gamma0-rtc_geo', scale=scale, dtype=2,
                     nodata=dict(zip(('linear', 'db'), nodata))[scale], outdir=outdir)
    
    if scene.sensor in ['S1A', 'S1B']:
        outname_base = scene.outname_base(extensions=basename_extensions)
        shutil.copyfile(os.path.join(scene.scene, 'manifest.safe'),
                        os.path.join(outdir, outname_base + '_manifest.safe'))
    
    if export_extra is not None:
        log.info('exporting extra products')
        for key in export_extra:
            if key in pix_geo:
                fname = n.get(key)
                diff.geocode_back(data_in=fname.replace('_geo', ''),
                                  width_in=master_par.range_samples,
                                  lookup_table=lut_final,
                                  data_out=fname,
                                  width_out=sim_width,
                                  interp_mode=func_geoback,
                                  logpath=path_log,
                                  outdir=tmpdir,
                                  shellscript=shellscript)
                par2hdr(n.dem_seg_geo + '.par', fname + '_.hdr')
            # SAR image products
            product_match = [x for x in products if x.endswith(key)]
            if len(product_match) > 0:
                for product in product_match:
                    for scale in scaling:
                        exporter(data_in=product, outdir=outdir, scale=scale, dtype=2,
                                 nodata=dict(zip(('linear', 'db'), nodata))[scale])
            # ancillary (DEM) products
            elif n.isfile(key) and key not in ['lut_init']:
                filename = n[key]
                dtype = 5 if key == 'ls_map_geo' else 2
                nodata = 0
                exporter(filename, outdir, dtype=dtype, nodata=nodata)
            else:
                log.warning('cannot export file {}'.format(key))
    
    shutil.copyfile(shellscript, os.path.join(outdir, os.path.basename(shellscript)))
    
    if cleanup:
        log.info('cleaning up temporary files')
        shutil.rmtree(tmpdir)


def ovs(parfile, targetres):
    """
    compute DEM oversampling factors for a target resolution in meters

    Parameters
    ----------
    parfile: str
        a GAMMA DEM parameter file
    targetres: int or float
        the target resolution in meters
    
    Returns
    -------
    tuple of float
        the oversampling factors for latitude and longitude
    """
    # read DEM parameter file
    dempar = ISPPar(parfile)
    
    # extract coordinates and pixel posting of the DEM
    if hasattr(dempar, 'post_north'):
        post_north, post_east = [abs(float(x)) for x in
                                 [dempar.post_north, dempar.post_east]]
    else:
        res_lat, res_lon = [abs(float(x)) for x in [dempar.post_lat, dempar.post_lon]]
        
        # compute center coordinate
        lat = float(dempar.corner_lat) - (res_lat * (dempar.nlines // 2))
        lon = float(dempar.corner_lon) + (res_lon * (dempar.width // 2))
        
        # convert DEM resolution to meters
        post_north = haversine(lat, lon, lat + res_lat, lon)
        post_east = haversine(lat, lon, lat, lon + res_lon)
    
    # compute resampling factors for the DEM
    ovs_lat = post_north / targetres
    ovs_lon = post_east / targetres
    return ovs_lat, ovs_lon


def multilook(infile, outfile, targetres, exist_ok=False, logpath=None, outdir=None, shellscript=None):
    """
    multilooking of SLC and MLI images

    if the image is in slant range the ground range resolution is computed by dividing the range pixel spacing by
    the sine of the incidence angle

    the looks in range and azimuth are chosen to approximate the target resolution by rounding the ratio between
    target resolution and ground range/azimuth pixel spacing to the nearest integer

    an ENVI HDR parameter file is automatically written for better handling in other software

    Parameters
    ----------
    infile: str
        a SAR image in GAMMA format with a parameter file of name <infile>.par
    outfile: str
        the name of the output GAMMA file
    targetres: int
        the target resolution in ground range
    exist_ok: bool
        allow existing output files and do not create new ones?
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format

    """
    # read the input parameter file
    par = ISPPar(infile + '.par')
    
    rlks, azlks = multilook_factors(sp_rg=par.range_pixel_spacing,
                                    sp_az=par.azimuth_pixel_spacing,
                                    tr_rg=targetres,
                                    tr_az=targetres,
                                    geometry=par.image_geometry,
                                    incidence=par.incidence_angle)
    
    pars = {'rlks': rlks,
            'azlks': azlks,
            'logpath': logpath,
            'shellscript': shellscript,
            'outdir': outdir}
    
    if par.image_format in ['SCOMPLEX', 'FCOMPLEX']:
        # multilooking for SLC images
        pars['SLC'] = infile
        pars['SLC_par'] = infile + '.par'
        pars['MLI'] = outfile
        pars['MLI_par'] = outfile + '.par'
        if do_execute(pars, ['MLI', 'MLI_par'], exist_ok):
            isp.multi_look(**pars)
            par2hdr(outfile + '.par', outfile + '.hdr')
    else:
        # multilooking for MLI images
        pars['MLI_in'] = infile
        pars['MLI_in_par'] = infile + '.par'
        pars['MLI_out'] = outfile
        pars['MLI_out_par'] = outfile + '.par'
        if do_execute(pars, ['MLI_out', 'MLI_out_par'], exist_ok):
            isp.multi_look_MLI(**pars)
            par2hdr(outfile + '.par', outfile + '.hdr')


def S1_deburst(burst1, burst2, burst3, name_out, rlks=5, azlks=1,
               replace=False, logpath=None, outdir=None, shellscript=None):
    """
    Debursting of Sentinel-1 SLC imagery in GAMMA
    
    The procedure consists of two steps. First antenna pattern deramping and
    then mosaicing of the single deramped bursts.
    For mosaicing, the burst boundaries are calculated from the number of looks in range (`rlks`)
    and azimuth (`azlks`), in this case 5 range looks and 1 azimuth looks.
    Alternately 10 range looks and 2 azimuth looks could be used.
    
    Parameters
    ----------
    burst1: str
        burst image 1
    burst2: str
        burst image 2
    burst3: str
        burst image 3
    name_out: str
        the name of the output file
    rlks: int
        the number of looks in range
    azlks: int
        the number of looks in azimuth
    replace: bool
        replace the burst images by the new file? If True, the three burst images will be deleted.
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format

    Returns
    -------
    
    """
    for burst in [burst1, burst2, burst3]:
        if not os.path.isfile(burst) or not os.path.isfile(burst + '.par') or not os.path.isfile(burst + '.tops_par'):
            raise IOError('input files missing; parameter files must be named e.g. {burst1}.par and {burst1}.tops_par')
    outpath = os.path.dirname(name_out)
    if not os.path.isdir(outpath):
        os.makedirs(outpath)
    tab_in = os.path.join(outpath, 'tab_deramp1')
    tab_out = os.path.join(outpath, 'tab_deramp2')
    with open(tab_in, 'w') as out1:
        with open(tab_out, 'w') as out2:
            for item in [burst1, burst2, burst3]:
                out1.write(item + '\t' + item + '.par\t' + item + '.tops_par\n')
                out2.write(item + '_drp\t' + item + '_drp.par\t' + item + '_drp.tops_par\n')
    
    isp.SLC_deramp_ScanSAR(SLC1_tab=tab_in,
                           SLC2_tab=tab_out,
                           mode=0,
                           phflg=0,
                           logpath=logpath,
                           outdir=outdir,
                           shellscript=shellscript)
    
    isp.SLC_mosaic_S1_TOPS(SLC_tab=tab_out,
                           SLC=name_out,
                           SLC_par=name_out + '.par',
                           rlks=rlks,
                           azlks=azlks,
                           logpath=logpath,
                           outdir=outdir,
                           shellscript=shellscript)
    if replace:
        for item in [burst1, burst2, burst3]:
            for subitem in [item + x for x in ['', '.par', '.tops_par']]:
                os.remove(subitem)
    for item in [burst1, burst2, burst3]:
        for subitem in [item + x for x in ['_drp', '_drp.par', '_drp.tops_par']]:
            os.remove(subitem)
    os.remove(tab_in)
    os.remove(tab_out)


def pixel_area_wrap(master, namespace, lut, path_log, outdir, shellscript):
    """
    helper function for computing pixel_area files in function geocode.
    
    Parameters
    ----------
    master: str
        the master SAR scene
    namespace: pyroSAR.gamma.auxil.Namespace
        an object collecting all output file names
    lut: str
        the name of the lookpu table
    path_log: str
        a directory to write command logfiles to
    outdir: str
        the directory to execute the command in
    shellscript: str
        a file to write the Gamma commands to in shell format

    Returns
    -------

    """
    master_par = ISPPar(master + '.par')
    
    pixel_area_args = {'MLI_par': master + '.par',
                       'DEM_par': namespace.dem_seg_geo + '.par',
                       'DEM': namespace.dem_seg_geo,
                       'lookup_table': lut,
                       'ls_map': namespace.ls_map_geo,
                       'inc_map': namespace.inc_geo,
                       'pix_sigma0': namespace.pix_area_sigma0,
                       'pix_gamma0': namespace.pix_area_gamma0,
                       'logpath': path_log,
                       'outdir': outdir,
                       'shellscript': shellscript}
    
    radcal_mli_args = {'MLI': master,
                       'MLI_par': master + '.par',
                       'OFF_par': '-',
                       'CMLI': master + '_cal',
                       'refarea_flag': 1,  # calculate sigma0, scale area by sin(inc_ang)/sin(ref_inc_ang)
                       'pix_area': namespace.pix_ellip_sigma0,
                       'logpath': path_log,
                       'outdir': outdir,
                       'shellscript': shellscript}
    
    # newer versions of GAMMA enable creating the ratio of ellipsoid based
    # pixel area and DEM-facet pixel area directly with command pixel_area
    if hasarg(diff.pixel_area, 'sig2gam_ratio'):
        namespace.appreciate(['pix_ratio'])
        pixel_area_args['sig2gam_ratio'] = namespace.pix_ratio
        diff.pixel_area(**pixel_area_args)
        
        if namespace.isappreciated('pix_ellip_sigma0'):
            isp.radcal_MLI(**radcal_mli_args)
            par2hdr(master + '.par', master + '_cal.hdr')
    else:
        # sigma0 = MLI * ellip_pix_sigma0 / pix_area_sigma0
        # gamma0 = MLI * ellip_pix_sigma0 / pix_area_gamma0
        namespace.appreciate(['pix_area_gamma0', 'pix_ellip_sigma0', 'pix_ratio'])
        # actual illuminated area as obtained from integrating DEM-facets (pix_area_sigma0 | pix_area_gamma0)
        diff.pixel_area(**pixel_area_args)
        
        # ellipsoid-based pixel area (ellip_pix_sigma0)
        isp.radcal_MLI(**radcal_mli_args)
        par2hdr(master + '.par', master + '_cal.hdr')
        
        # ratio of ellipsoid based pixel area and DEM-facet pixel area
        lat.ratio(d1=namespace.pix_ellip_sigma0,
                  d2=namespace.pix_area_gamma0,
                  ratio=namespace.pix_ratio,
                  width=master_par.range_samples,
                  bx=1,
                  by=1,
                  logpath=path_log,
                  outdir=outdir,
                  shellscript=shellscript)
    
    if namespace.isappreciated('gs_ratio'):
        lat.ratio(d1=namespace.pix_area_gamma0,
                  d2=namespace.pix_area_sigma0,
                  ratio=namespace.gs_ratio,
                  width=master_par.range_samples,
                  bx=1,
                  by=1,
                  logpath=path_log,
                  outdir=outdir,
                  shellscript=shellscript)
    
    for item in ['pix_area_sigma0', 'pix_area_gamma0',
                 'pix_ratio', 'pix_ellip_sigma0', 'gs_ratio']:
        if namespace.isappreciated(item):
            par2hdr(master + '.par', namespace[item] + '.hdr')
