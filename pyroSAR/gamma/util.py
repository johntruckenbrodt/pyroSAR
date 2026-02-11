###############################################################################
# universal core routines for processing SAR images with GAMMA

# Copyright (c) 2014-2023, the pyroSAR Developers.

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
import numpy as np
from spatialist import haversine, Raster
from spatialist.ancillary import union, finder

from ..S1 import OSV
from ..drivers import ID, identify_many
from . import ISPPar, Namespace, par2hdr
from ..ancillary import multilook_factors, hasarg, groupby, Lock
from pyroSAR.examine import ExamineSnap
from .auxil import do_execute

import logging

log = logging.getLogger(__name__)

try:
    from .api import diff, disp, isp, lat
except ImportError:
    pass


def calibrate(id, directory, return_fnames=False,
              logpath=None, outdir=None, shellscript=None):
    """
    radiometric calibration of SAR scenes
    
    Parameters
    ----------
    id: ~pyroSAR.drivers.ID
        an SAR scene object of type pyroSAR.ID or any subclass
    directory: str
        the directory to search for GAMMA calibration candidates
    return_fnames: bool
        return the names of the output image files? Default: False.
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the GAMMA commands to in shell format

    Returns
    -------
    List[str] or None
    """
    cname = type(id).__name__
    new = []
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
                new.append(image + '_cal')
    
    elif cname == 'EORC_PSR':
        for image in id.getGammaImages(directory):
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
                new.append(image + '_cal')
    
    elif cname == 'ESA':
        k_db = {'ASAR': 55., 'ERS1': 58.24, 'ERS2': 59.75}[id.sensor]
        inc_ref = 90. if id.sensor == 'ASAR' else 23.
        imgs = id.getGammaImages(directory)
        candidates = [x for x in imgs if re.search('_pri$', x)]
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
            new.append(out)
    
    elif cname == 'SAFE':
        log.info('calibration already performed during import')
    
    else:
        msg = f'calibration for class {cname} is not implemented yet'
        raise NotImplementedError(msg)
    
    if return_fnames and len(new) > 0:
        return new


def convert2gamma(id, directory, S1_tnr=True, S1_bnr=True,
                  basename_extensions=None, exist_ok=False,
                  return_fnames=False,
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
    basename_extensions: list[str] or None
        names of additional parameters to append to the basename, e.g. ['orbitNumber_rel']
    exist_ok: bool
        allow existing output files and do not create new ones?
    return_fnames: bool
        return the names of the output image files? Default: False.
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the GAMMA commands to in bash format

    Returns
    -------
    list[str] or None
        the sorted image file names if ``return_fnames=True`` and None otherwise
    """
    
    if not isinstance(id, ID):
        raise IOError('id must be of type pyroSAR.ID')
    
    if id.compression is not None:
        raise RuntimeError('scene is not yet unpacked')
    
    os.makedirs(directory, exist_ok=True)
    
    fnames = []
    
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
                    
                    with Lock(outname):
                        if do_execute(pars, ['SLC', 'SLC_par'], exist_ok):
                            isp.par_ESA_ERS(**pars)
                            par2hdr(outname + '.par', outname + '.hdr')
                    fnames.append(outname)
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
                
                with Lock(outname):
                    if do_execute(pars, ['SLC', 'SLC_par'], exist_ok):
                        isp.par_EORC_PALSAR(**pars)
                        par2hdr(outname + '.par', outname + '.hdr')
            else:
                outname_base = '{}_{}_mli_geo'.format(outname_base, polarization)
                outname = os.path.join(directory, outname_base)
                
                pars['MLI'] = outname
                pars['MLI_par'] = outname + '.par'
                pars['DEM_par'] = outname + '_dem.par'
                
                with Lock(outname):
                    if do_execute(pars, ['MLI', 'MLI_par', 'DEM_par'], exist_ok):
                        diff.par_EORC_PALSAR_geo(**pars)
                        par2hdr(outname + '.par', outname + '.hdr')
            fnames.append(outname)
    
    elif cname == 'EORC_PSR':
        images = id.findfiles('^sar.')
        facter_m = id.findfiles('facter_m.dat')
        led = id.findfiles('LED-ALOS2')
        
        for image in images:
            polarization = re.search('[HV]{2}', os.path.basename(image)).group(0)
            outname_base = id.outname_base(extensions=basename_extensions)
            outname_base = '{}_{}'.format(outname_base, polarization)
            outname = os.path.join(directory, outname_base) + '_mli'
            fnames.append(outname)
            
            pars = {'facter_m': facter_m,
                    'CEOS_leader': led,
                    'SLC_par': outname + '.par',
                    'pol': polarization,
                    'pls_mode': 2,
                    'KC_data': image,
                    'pwr': outname,
                    'logpath': logpath,
                    'outdir': outdir,
                    'shellscript': shellscript}
            
            with Lock(outname):
                if do_execute(pars, ['pwr', 'SLC_par'], exist_ok):
                    isp.par_KC_PALSAR_slr(**pars)
                    par2hdr(outname + '.par', outname + '.hdr')
    
    elif cname == 'ESA':
        """
        the command par_ASAR also accepts a K_dB argument for calibration
        in which case the resulting image names will carry the suffix grd;
        this is not implemented here but instead in function calibrate
        """
        outname = os.path.join(directory, id.outname_base(extensions=basename_extensions))
        with Lock(outname):
            
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
                fnames.append(outname)
                if outname.endswith('.par'):
                    par2hdr(outname, outname.replace('.par', '.hdr'))
    
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
            basename = '_'.join(fields)
            outname = os.path.join(directory, basename)
            
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
                base_new = basename.replace(old, swath)
                outname = os.path.join(os.path.dirname(outname), base_new)
                pars['SLC'] = outname
                pars['SLC_par'] = outname + '.par'
                pars['TOPS_par'] = outname + '.tops_par'
                with Lock(outname):
                    if do_execute(pars, ['SLC', 'SLC_par', 'TOPS_par'], exist_ok):
                        isp.par_S1_SLC(**pars)
                        par2hdr(outname + '.par', outname + '.hdr')
            else:
                if hasarg(isp.par_S1_GRD, 'edge_flag'):
                    if S1_bnr:
                        pars['edge_flag'] = 2
                    else:
                        pars['edge_flag'] = 0
                else:
                    if S1_bnr:
                        raise RuntimeError("The command par_S1_GRD of this GAMMA "
                                           "version does not support border noise "
                                           "removal. You may want to consider "
                                           "pyroSAR's own method for this task.")
                pars['MLI'] = outname
                pars['MLI_par'] = outname + '.par'
                with Lock(outname):
                    if do_execute(pars, ['MLI', 'MLI_par'], exist_ok):
                        isp.par_S1_GRD(**pars)
                        par2hdr(outname + '.par', outname + '.hdr')
            fnames.append(outname)
    
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
                with Lock(outname):
                    if do_execute(pars, ['SLC', 'SLC_par'], exist_ok):
                        isp.par_TX_SLC(**pars)
                        par2hdr(outname + '.par', outname + '.hdr')
            
            elif id.product == 'MGD':
                outname += '_mli'
                pars['GeoTIFF'] = image
                pars['GRD_par'] = outname + '.par'
                pars['GRD'] = outname
                with Lock(outname):
                    if do_execute(pars, ['GRD', 'GRD_par'], exist_ok):
                        isp.par_TX_GRD(**pars)
                        par2hdr(outname + '.par', outname + '.hdr')
            
            elif id.product in ['GEC', 'EEC']:
                outname += '_mli_geo'
                pars['GeoTIFF'] = image
                pars['MLI_par'] = outname + '.par'
                pars['DEM_par'] = outname + '_dem.par'
                pars['GEO'] = outname
                with Lock(outname):
                    if do_execute(pars, ['GEO', 'MLI_par', 'DEM_par'], exist_ok):
                        diff.par_TX_geo(**pars)
                        par2hdr(outname + '.par', outname + '.hdr')
            else:
                raise RuntimeError('unknown product: {}'.format(id.product))
            fnames.append(outname)
    
    else:
        raise NotImplementedError('conversion for class {} is not implemented yet'.format(cname))
    
    if return_fnames:
        return sorted(fnames)


def correctOSV(id, directory, osvdir=None, osvType='POE', timeout=20,
               logpath=None, outdir=None, shellscript=None, url_option=1):
    """
    correct GAMMA parameter files with orbit state vector information from dedicated OSV files;
    OSV files are downloaded automatically to either the defined `osvdir` or relative to the
    user's home directory: `~/.snap/auxdata/Orbits/Sentinel-1`.
    
    Parameters
    ----------
    id: ~pyroSAR.drivers.ID
        the scene to be corrected
    directory: str or None
        a directory to be scanned for files associated with the scene, e.g. an SLC in GAMMA format.
        If the OSV file is packed in a zip file it will be unpacked to a subdirectory `osv`.
    osvdir: str or None
        the directory of the OSV files. Default None: use the SNAP directory
        as configured via `pyroSAR.examine.ExamineSnap` or, if SNAP is not
        installed, `~/.snap/auxdata/Orbits/Sentinel-1` (SNAP default).
        Subdirectories POEORB and RESORB are created automatically.
    osvType: str or list[str]
        the OSV type (POE|RES) to be used
    timeout: int or tuple or None
        the timeout in seconds for downloading OSV files as provided to :func:`requests.get`
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the GAMMA commands to in shell format
    url_option: int
        the OSV download URL option; see :meth:`pyroSAR.S1.OSV.catch`
    
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
    
    if id.sensor not in ['S1A', 'S1B', 'S1C', 'S1D']:
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
        id.getOSV(osvdir, osvType, timeout=timeout, url_option=url_option)
    except URLError:
        log.warning('..no internet access')
    
    parfiles = finder(directory, ['*.par'])
    parfiles = [x for x in parfiles if ISPPar(x).filetype == 'isp']
    # read parameter file entries into object
    with ISPPar(parfiles[0]) as par:
        # extract acquisition time stamp
        timestamp = datetime.strptime(par.date, '%Y-%m-%dT%H:%M:%S.%f').strftime('%Y%m%dT%H%M%S')
    
    # find an OSV file matching the time stamp and defined OSV type(s)
    with OSV(osvdir, timeout=timeout) as osv:
        osvfile = osv.match(sensor=id.sensor, timestamp=timestamp, osvtype=osvType)
    if not osvfile:
        raise RuntimeError('no Orbit State Vector file found')
    
    if osvfile.endswith('.zip'):
        osvdir = os.path.join(directory, 'osv')
        with zf.ZipFile(osvfile) as zip:
            zip.extractall(path=osvdir)
        osvfile = os.path.join(osvdir, os.path.basename(osvfile).replace('.zip', ''))
    
    # update the GAMMA parameter file with the selected orbit state vectors
    log.debug('correcting state vectors with file {}'.format(osvfile))
    for par in parfiles:
        log.debug(par)
        with Lock(par.replace('.par', '')):
            isp.S1_OPOD_vec(SLC_par=par,
                            OPOD=osvfile,
                            logpath=logpath,
                            outdir=outdir,
                            shellscript=shellscript)


def gc_map_wrap(image, namespace, dem, spacing, exist_ok=False,
                logpath=None, outdir=None, shellscript=None):
    """
    helper function for computing DEM products in function geocode.

    Parameters
    ----------
    image: str
        the reference SAR image
    namespace: pyroSAR.gamma.auxil.Namespace
        an object collecting all output file names
    dem: str
        the digital elevation model
    spacing: int or float
        the target pixel spacing in meters
    exist_ok: bool
        allow existing output files and do not create new ones?
    logpath: str
        a directory to write command logfiles to
    outdir: str
        the directory to execute the command in
    shellscript: str
        a file to write the GAMMA commands to in shell format

    Returns
    -------

    """
    # compute DEM oversampling factors; will be 1 for range and
    # azimuth if the DEM spacing matches the target spacing
    ovs_lat, ovs_lon = ovs(dem + '.par', spacing)
    
    image_par = ISPPar(image + '.par')
    
    gc_map_args = {'DEM_par': dem + '.par',
                   'DEM': dem,
                   'DEM_seg_par': namespace.dem_seg_geo + '.par',
                   'DEM_seg': namespace.dem_seg_geo,
                   'lookup_table': namespace.lut_init,
                   'lat_ovr': ovs_lat,
                   'lon_ovr': ovs_lon,
                   'sim_sar': namespace.sim_sar_geo,
                   'u': namespace.u_geo,
                   'v': namespace.v_geo,
                   'inc': namespace.inc_geo,
                   'psi': namespace.psi_geo,
                   'pix': namespace.pix_geo,
                   'ls_map': namespace.ls_map_geo,
                   'frame': 8,
                   'ls_mode': 2,
                   'logpath': logpath,
                   'shellscript': shellscript,
                   'outdir': outdir}
    out_id = ['DEM_seg_par', 'DEM_seg', 'lookup_table', 'sim_sar',
              'u', 'v', 'inc', 'psi', 'pix', 'ls_map']
    
    # remove all output files to make sure they are replaced and not updated
    if not exist_ok:
        for id in out_id:
            base = gc_map_args[id]
            if base != '-':
                for suffix in ['', '.par', '.hdr']:
                    fname = base + suffix
                    if os.path.isfile(fname):
                        os.remove(fname)
    
    if image_par.image_geometry == 'GROUND_RANGE':
        gc_map_args.update({'GRD_par': image + '.par'})
        if do_execute(gc_map_args, out_id, exist_ok):
            diff.gc_map_grd(**gc_map_args)
    else:
        gc_map_args.update({'MLI_par': image + '.par'})
        if do_execute(gc_map_args, out_id, exist_ok):
            # gc_map2 is the successor of gc_map. However, earlier versions
            # did not yet come with full functionality.
            gc_map2_ok = False
            if 'gc_map2' in dir(diff):
                keys = list(gc_map_args.keys())
                keys.remove('ls_mode')
                gc_map2_ok = all([hasarg(diff.gc_map2, x) for x in keys])
            if gc_map2_ok:
                del gc_map_args['ls_mode']
                diff.gc_map2(**gc_map_args)
            else:
                # gc_map might have an argument OFF_par, which is not needed for SLC/MLI geocoding
                if hasarg(diff.gc_map, 'OFF_par'):
                    gc_map_args.update({'OFF_par': '-'})
                diff.gc_map(**gc_map_args)
    
    # create ENVI header files for all created images
    for item in ['dem_seg_geo', 'sim_sar_geo', 'u_geo', 'v_geo',
                 'psi_geo', 'pix_geo', 'inc_geo', 'ls_map_geo']:
        if namespace.isappreciated(item):
            mods = {'data_type': 1} if item == 'ls_map_geo' else None
            par2hdr(namespace.dem_seg_geo + '.par', namespace.get(item) + '.hdr', mods)


def geocode(scene, dem, tmpdir, outdir, spacing, scaling='linear', func_geoback=1,
            nodata=(0, -99), update_osv=True, osvdir=None, allow_RES_OSV=False,
            cleanup=True, export_extra=None, basename_extensions=None,
            removeS1BorderNoiseMethod='gamma', refine_lut=False, rlks=None, azlks=None,
            s1_osv_url_option=1):
    """
    general function for radiometric terrain correction (RTC) and geocoding of SAR backscatter images with GAMMA.
    Applies the RTC method by :cite:t:`Small2011` to retrieve gamma nought RTC backscatter.
    
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
    spacing: float or int
        the target pixel spacing in meters
    scaling: str or list[str]
        the value scaling of the backscatter values; either 'linear', 'db' or a list of both, i.e. ['linear', 'db']
    func_geoback: {0, 1, 2, 3, 4, 5, 6, 7}
        backward geocoding interpolation mode (see GAMMA command `geocode_back`)
        
         - 0: nearest-neighbor
         - 1: bicubic spline (default)
         - 2: bicubic-spline, interpolate log(data)
         - 3: bicubic-spline, interpolate sqrt(data)
         - 4: B-spline interpolation (default B-spline degree: 5)
         - 5: B-spline interpolation sqrt(x) (default B-spline degree: 5)
         - 6: Lanczos interpolation (default Lanczos function order: 5)
         - 7: Lanczos interpolation sqrt(x) (default Lanczos function order: 5)
        
        .. note::
        
            log and sqrt interpolation modes should only be used with non-negative data!
        
        .. note::
        
            GAMMA recommendation for MLI data: "The interpolation should be performed on
            the square root of the data. A mid-order (3 to 5) B-spline interpolation is recommended."
    nodata: tuple[float or int]
        the nodata values for the output files; defined as a tuple with two values, the first for linear,
        the second for logarithmic scaling
    update_osv: bool
        update the orbit state vectors?
    osvdir: str or None
        a directory for Orbit State Vector files;
        this is currently only used by for Sentinel-1 where two subdirectories POEORB and RESORB are created;
        if set to None, a subdirectory OSV is created in the directory of the unpacked scene.
    allow_RES_OSV: bool
        also allow the less accurate RES orbit files to be used?
        Otherwise the function will raise an error if no POE file exists.
    cleanup: bool
        should all files written to the temporary directory during function execution be deleted after processing?
    export_extra: list[str] or None
        a list of image file IDs to be exported to outdir
        
         - format is GeoTIFF if the file is geocoded and ENVI otherwise. Non-geocoded images can be converted via GAMMA
           command data2tiff yet the output was found impossible to read with GIS software
         - scaling of SAR image products is applied as defined by parameter `scaling`
         - see Notes for ID options
    basename_extensions: list[str] or None
        names of additional parameters to append to the basename, e.g. ['orbitNumber_rel']
    removeS1BorderNoiseMethod: str or None
        the S1 GRD border noise removal method to be applied, See :func:`pyroSAR.S1.removeGRDBorderNoise` for details; one of the following:
        
         - 'ESA': the pure implementation as described by ESA
         - 'pyroSAR': the ESA method plus the custom pyroSAR refinement
         - 'gamma': the GAMMA implementation of :cite:`Ali2018`
         - None: do not remove border noise
    refine_lut: bool
        should the LUT for geocoding be refined using pixel area normalization?
    rlks: int or None
        the number of range looks. If not None, overrides the computation done by function
        :func:`pyroSAR.ancillary.multilook_factors` based on the image pixel spacing and the target spacing.
    azlks: int or None
        the number of azimuth looks. Like `rlks`.
    s1_osv_url_option: int
        the OSV download URL option; see :meth:`pyroSAR.S1.OSV.catch`
    
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
    
      * **dem_seg_geo**: dem subsetted to the extent of the intersection between input DEM and SAR image
      * (**u_geo**): zenith angle of surface normal vector n (angle between z and n)
      * (**v_geo**): orientation angle of n (between x and projection of n in xy plane)
      * **inc_geo**: local incidence angle (between surface normal and look vector)
      * (**psi_geo**): projection angle (between surface normal and image plane normal)
      * **ls_map_geo**: layover and shadow map
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
    >>> geocode(scene=filename, dem='demfile', outdir='outdir', spacing=20, scaling='db',
    >>>         export_extra=['dem_seg_geo', 'inc_geo', 'ls_map_geo'])
    
    .. figure:: figures/gamma_geocode.svg
        :align: center
        
        Workflow diagram for function geocode for processing a Sentinel-1 Ground Range
        Detected (GRD) scene to radiometrically terrain corrected (RTC) gamma nought backscatter.
    
    """
    
    # experimental option to reuse intermediate products; currently affects:
    # - scene unpacking
    # - conversion to GAMMA format
    # - multilooking
    # - DEM product generation
    # - terrain flattening
    exist_ok = False
    
    scenes = scene if isinstance(scene, list) else [scene]
    if len(scenes) > 2:
        raise RuntimeError("currently only one or two scenes can be passed via argument 'scene'")
    scenes = identify_many(scenes)
    ref = scenes[0]
    
    if ref.sensor not in ['S1A', 'S1B', 'S1C', 'S1D', 'PALSAR-2']:
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
                scene.unpack(tmpdir, exist_ok=exist_ok)
            except RuntimeError:
                log.info('scene was attempted to be processed before, exiting')
                return
    
    path_log = os.path.join(tmpdir, 'logfiles')
    if not os.path.isdir(path_log):
        os.makedirs(path_log)
    
    for scene in scenes:
        if scene.sensor in ['S1A', 'S1B', 'S1C', 'S1D'] and removeS1BorderNoiseMethod in ['ESA', 'pyroSAR']:
            log.info('removing border noise')
            scene.removeGRDBorderNoise(method=removeS1BorderNoiseMethod)
    
    log.info('converting scene to GAMMA format')
    gamma_bnr = True if removeS1BorderNoiseMethod == 'gamma' else False
    images = []
    for scene in scenes:
        files = convert2gamma(scene, directory=tmpdir, logpath=path_log, outdir=tmpdir,
                              basename_extensions=basename_extensions, shellscript=shellscript,
                              S1_bnr=gamma_bnr, exist_ok=exist_ok, return_fnames=True)
        images.extend(files)
    
    if update_osv:
        for scene in scenes:
            if scene.sensor in ['S1A', 'S1B', 'S1C', 'S1D']:
                log.info('updating orbit state vectors')
                if allow_RES_OSV:
                    osvtype = ['POE', 'RES']
                else:
                    osvtype = 'POE'
                try:
                    correctOSV(id=scene, directory=tmpdir, osvdir=osvdir, osvType=osvtype,
                               url_option=s1_osv_url_option,
                               logpath=path_log, outdir=tmpdir, shellscript=shellscript)
                except RuntimeError:
                    msg = 'orbit state vector correction failed for scene {}'
                    log.warning(msg.format(scene.scene))
                    return
    
    log.info('calibrating')
    images_cal = []
    for scene in scenes:
        files = calibrate(id=scene, directory=tmpdir, return_fnames=True,
                          logpath=path_log, outdir=tmpdir, shellscript=shellscript)
        if files is not None:
            images_cal.extend(files)
    if len(images_cal) > 0:
        images = images_cal
    
    if len(scenes) > 1:
        images_new = []
        groups = groupby(images, 'polarization')
        for group in groups:
            out = group[0] + '_cat'
            out_par = out + '.par'
            all_exist = all([os.path.isfile(x) for x in [out, out_par]])
            if not all_exist:
                log.info('mosaicing scenes')
                isp.MLI_cat(MLI1=group[0],
                            MLI1_par=group[0] + '.par',
                            MLI2=group[1],
                            MLI2_par=group[1] + '.par',
                            MLI3=out,
                            MLI3_par=out_par,
                            logpath=path_log, outdir=tmpdir, shellscript=shellscript)
                par2hdr(out_par, out + '.hdr')
            images_new.append(out)
        images = images_new
    
    if scene.sensor in ['S1A', 'S1B', 'S1C', 'S1D']:
        log.info('multilooking')
        groups = groupby(images, 'polarization')
        images = []
        for group in groups:
            out = group[0].replace('IW1', 'IW_') + '_mli'
            infile = group[0] if len(group) == 1 else group
            multilook(infile=infile, outfile=out, spacing=spacing,
                      rlks=rlks, azlks=azlks, exist_ok=exist_ok,
                      logpath=path_log, outdir=tmpdir, shellscript=shellscript)
            images.append(out)
    products = list(images)
    reference = images[0]
    
    # create output names for files to be written
    # appreciated files will be written
    n = Namespace(tmpdir, scene.outname_base(extensions=basename_extensions))
    n.appreciate(['dem_seg_geo', 'lut_init', 'inc_geo', 'ls_map_geo'])
    
    pix_geo = []
    if export_extra is not None:
        n.appreciate(export_extra)
        pix = ['pix_area_sigma0', 'pix_area_gamma0', 'pix_ratio', 'gs_ratio', 'pix_ellip_sigma0']
        for item in pix:
            if item + '_geo' in export_extra:
                pix_geo.append(item + '_geo')
                n.appreciate([item])
    
    if refine_lut:
        n.appreciate(['pix_area_sigma0'])
    
    reference_par = ISPPar(reference + '.par')
    ######################################################################
    # geocoding and DEM product generation ###############################
    ######################################################################
    log.info('geocoding and creating DEM products')
    gc_map_wrap(image=reference, namespace=n, dem=dem, spacing=spacing, exist_ok=exist_ok,
                logpath=path_log, outdir=tmpdir, shellscript=shellscript)
    
    sim_width = ISPPar(n.dem_seg_geo + '.par').width
    ######################################################################
    # RTC reference area computation #####################################
    ######################################################################
    log.info('computing pixel area (for radiometric terrain correction, rtc)')
    pixel_area_wrap(image=reference, namespace=n, lut=n.lut_init, exist_ok=exist_ok,
                    logpath=path_log, outdir=tmpdir, shellscript=shellscript)
    
    ######################################################################
    # lookup table refinement ############################################
    ######################################################################
    lut_final = n.lut_init
    if refine_lut:
        log.info('refining lookup table')
        # Refinement of geocoding lookup table
        diff.create_diff_par(PAR_1=reference + '.par',
                             PAR_2='-',
                             DIFF_par=reference + '_diff.par',
                             PAR_type=1,
                             iflg=0,
                             logpath=path_log,
                             outdir=tmpdir,
                             shellscript=shellscript)
        # Refinement of lookup table
        # for "shift" data offset window size enlarged twice to 512 and 256, for data without shift 256 128
        diff.offset_pwrm(MLI_1=n.pix_area_sigma0,
                         MLI_2=reference,
                         DIFF_par=reference + '_diff.par',
                         offs=reference + '_offs',
                         ccp=reference + '_ccp',
                         rwin=512,
                         azwin=256,
                         offsets=reference + '_offsets.txt',
                         n_ovr=2,
                         nr=64,
                         naz=32,
                         thres=0.2,
                         logpath=path_log,
                         outdir=tmpdir,
                         shellscript=shellscript)
        # par2hdr(master + '.par', master + '_offs' + '.hdr')
        diff.offset_fitm(offs=reference + '_offs',
                         ccp=reference + '_ccp',
                         DIFF_par=reference + '_diff.par',
                         coffs=reference + '_coffs',
                         coffsets=reference + '_coffsets',
                         thres=0.2,
                         npoly=4,
                         logpath=path_log,
                         outdir=tmpdir,
                         shellscript=shellscript)
        # Updating of the look-up table
        diff.gc_map_fine(gc_in=lut_final,
                         width=sim_width,
                         DIFF_par=reference + '_diff.par',
                         gc_out=lut_final + '.fine',
                         ref_flg=1,
                         logpath=path_log,
                         outdir=tmpdir,
                         shellscript=shellscript)
        # Reproduce pixel area estimate
        pixel_area_wrap(image=reference, namespace=n, lut=lut_final + '.fine',
                        logpath=path_log, outdir=tmpdir, shellscript=shellscript)
        lut_final = lut_final + '.fine'
    ######################################################################
    # radiometric terrain correction and back-geocoding ##################
    ######################################################################
    log.info('applying rtc and back-geocoding')
    for image in images:
        if 'lat' in locals():
            lat.product(data_1=image,
                        data_2=n.pix_ratio,
                        product=image + '_gamma0-rtc',
                        width=reference_par.range_samples,
                        bx=1,
                        by=1,
                        logpath=path_log,
                        outdir=tmpdir,
                        shellscript=shellscript)
        else:
            lat_product(data_in1=image,
                        data_in2=n.pix_ratio,
                        data_out=image + '_gamma0-rtc')
        par2hdr(reference + '.par', image + '_gamma0-rtc.hdr')
        diff.geocode_back(data_in=image + '_gamma0-rtc',
                          width_in=reference_par.range_samples,
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
                width = reference_par.range_samples
                refpar = reference + '.par'
            if 'lat' in locals():
                lat.linear_to_dB(data_in=data_in,
                                 data_out=data_in + '_db',
                                 width=width,
                                 inverse_flag=0,
                                 null_value=nodata,
                                 logpath=path_log,
                                 outdir=tmpdir,
                                 shellscript=shellscript)
            else:
                lat_linear_to_db(data_in=data_in,
                                 data_out=data_in + '_db')
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
    
    if scene.sensor in ['S1A', 'S1B', 'S1C', 'S1D']:
        outname_base = scene.outname_base(extensions=basename_extensions)
        shutil.copyfile(os.path.join(scene.scene, 'manifest.safe'),
                        os.path.join(outdir, outname_base + '_manifest.safe'))
    
    if export_extra is not None:
        log.info('back-geocoding and exporting extra products')
        for key in export_extra:
            if key in pix_geo:
                fname = n.get(key)
                diff.geocode_back(data_in=fname.replace('_geo', ''),
                                  width_in=reference_par.range_samples,
                                  lookup_table=lut_final,
                                  data_out=fname,
                                  width_out=sim_width,
                                  interp_mode=func_geoback,
                                  logpath=path_log,
                                  outdir=tmpdir,
                                  shellscript=shellscript)
                par2hdr(n.dem_seg_geo + '.par', fname + '.hdr')
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


def lat_linear_to_db(data_in, data_out):
    """
    Alternative to LAT module command linear_to_dB.

    Parameters
    ----------
    data_in: str
        the input data file
    data_out: str
        the output data file

    Returns
    -------

    """
    with Raster(data_in) as ras:
        a1 = ras.array()
        a1[a1 <= 0] = np.nan
        out = 10 * np.log10(a1)
        tmp = data_out + '_tmp'
        ras.write(outname=tmp, array=out, format='ENVI',
                  nodata=0, dtype='float32')
    disp.swap_bytes(infile=tmp, outfile=data_out, swap_type=4)
    shutil.copy(src=data_in + '.hdr', dst=data_out + '.hdr')
    for item in [tmp, tmp + '.hdr', tmp + '.aux.xml']:
        os.remove(item)


def lat_product(data_in1, data_in2, data_out):
    """
    Alternative to LAT module command product.

    Parameters
    ----------
    data_in1: str
        input data file 1
    data_in2: str
        input data file 2
    data_out: str
        the output data file

    Returns
    -------

    """
    with Raster(data_in1) as ras:
        a1 = ras.array()
        a1[a1 == 0] = np.nan
    with Raster(data_in2) as ras:
        a2 = ras.array()
        a2[a2 == 0] = np.nan
        out = a1 * a2
        tmp = data_out + '_tmp'
        ras.write(outname=tmp, array=out, format='ENVI',
                  nodata=0, dtype='float32')
    disp.swap_bytes(infile=tmp, outfile=data_out, swap_type=4)
    shutil.copy(src=data_in1 + '.hdr', dst=data_out + '.hdr')
    for item in [tmp, tmp + '.hdr', tmp + '.aux.xml']:
        if os.path.isfile(item):
            os.remove(item)


def lat_ratio(data_in1, data_in2, data_out):
    """
    Alternative to LAT module command ratio.

    Parameters
    ----------
    data_in1: str
        input data file 1
    data_in2: str
        input data file 2
    data_out: str
        the output data file

    Returns
    -------

    """
    with Raster(data_in1) as ras:
        a1 = ras.array()
        a1[a1 == 0] = np.nan
    with Raster(data_in2) as ras:
        a2 = ras.array()
        a2[a2 == 0] = np.nan
        out = a1 / a2
        tmp = data_out + '_tmp'
        ras.write(outname=tmp, array=out, format='ENVI',
                  nodata=0, dtype='float32')
    disp.swap_bytes(infile=tmp, outfile=data_out, swap_type=4)
    shutil.copy(src=data_in1 + '.hdr', dst=data_out + '.hdr')
    for item in [tmp, tmp + '.hdr', tmp + '.aux.xml']:
        if os.path.isfile(item):
            os.remove(item)


def multilook(infile, outfile, spacing, rlks=None, azlks=None,
              exist_ok=False, logpath=None, outdir=None, shellscript=None):
    """
    Multilooking of SLC and MLI images.

    If the image is in slant range the ground range resolution is computed by dividing the range pixel spacing by
    the sine of the incidence angle.

    The looks in range and azimuth are chosen to approximate the target resolution by rounding the ratio between
    target resolution and ground range/azimuth pixel spacing to the nearest integer.

    An ENVI HDR parameter file is automatically written for better handling in other software.

    Parameters
    ----------
    infile: str or list[str]
        one of the following:

        - a SAR image in GAMMA format with a parameter file <infile>.par
        - a list of ScanSAR SLC swaths with parameter files <slc>.par and <slc>.tops_par; in this case a text file
          <outfile>_slc-tab.txt will be created, which is passed to the GAMMA command ``multi_look_ScanSAR``
    outfile: str
        the name of the output GAMMA MLI file
    spacing: int
        the target pixel spacing in ground range
    rlks: int or None
        the number of range looks. If not None, overrides the computation done by function
        :func:`pyroSAR.ancillary.multilook_factors` based on the image pixel spacing and the target spacing.
    azlks: int or None
        the number of azimuth looks. Like `rlks`.
    exist_ok: bool
        allow existing output files and do not create new ones?
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the GAMMA commands to in shell format

    See Also
    --------
    pyroSAR.ancillary.multilook_factors
    """
    # read the input parameter file
    if isinstance(infile, str):
        par = ISPPar(infile + '.par')
        range_pixel_spacing = par.range_pixel_spacing
        azimuth_pixel_spacing = par.azimuth_pixel_spacing
        incidence_angle = par.incidence_angle
        image_geometry = par.image_geometry
        image_format = par.image_format
    elif isinstance(infile, list):
        par = [ISPPar(x + '.par') for x in infile]
        range_pixel_spacings = [getattr(x, 'range_pixel_spacing') for x in par]
        range_pixel_spacing = sum(range_pixel_spacings) / len(par)
        azimuth_pixel_spacings = [getattr(x, 'azimuth_pixel_spacing') for x in par]
        azimuth_pixel_spacing = sum(azimuth_pixel_spacings) / len(par)
        incidence_angles = [getattr(x, 'incidence_angle') for x in par]
        incidence_angle = sum(incidence_angles) / len(par)
        image_geometry = par[0].image_geometry
        image_format = par[0].image_format
    else:
        raise TypeError("'infile' must be str or list")
    
    if rlks is None and azlks is None:
        rlks, azlks = multilook_factors(source_rg=range_pixel_spacing,
                                        source_az=azimuth_pixel_spacing,
                                        target=spacing,
                                        geometry=image_geometry,
                                        incidence=incidence_angle)
    if [rlks, azlks].count(None) > 0:
        raise RuntimeError("'rlks' and 'azlks' must either both be integers or None")
    
    pars = {'rlks': rlks,
            'azlks': azlks,
            'logpath': logpath,
            'shellscript': shellscript,
            'outdir': outdir}
    
    if image_format in ['SCOMPLEX', 'FCOMPLEX']:
        # multilooking of SLC images
        pars['MLI'] = outfile
        pars['MLI_par'] = outfile + '.par'
        if isinstance(infile, str):
            pars['SLC'] = infile
            pars['SLC_par'] = infile + '.par'
            if do_execute(pars, ['MLI', 'MLI_par'], exist_ok):
                isp.multi_look(**pars)
                par2hdr(outfile + '.par', outfile + '.hdr')
        else:
            slcpar = [x + '.par' for x in infile]
            topspar = [x + '.tops_par' for x in infile]
            slc_tab = outfile + '_slc-tab.txt'
            if not os.path.isfile(slc_tab) or not exist_ok:
                with open(slc_tab, 'w') as tab:
                    for item in zip(infile, slcpar, topspar):
                        tab.write(' '.join(item) + '\n')
            pars['SLC_tab'] = slc_tab
            if do_execute(pars, ['MLI', 'MLI_par'], exist_ok):
                if 'multi_look_ScanSAR' in dir(isp):
                    isp.multi_look_ScanSAR(**pars)
                else:
                    isp.multi_S1_TOPS(**pars)
                par2hdr(outfile + '.par', outfile + '.hdr')
    else:
        # multilooking of MLI images
        pars['MLI_in'] = infile
        pars['MLI_in_par'] = infile + '.par'
        pars['MLI_out'] = outfile
        pars['MLI_out_par'] = outfile + '.par'
        if do_execute(pars, ['MLI_out', 'MLI_out_par'], exist_ok):
            isp.multi_look_MLI(**pars)
            par2hdr(outfile + '.par', outfile + '.hdr')


def ovs(parfile, spacing):
    """
    compute DEM oversampling factors for a target resolution in meters

    Parameters
    ----------
    parfile: str
        a GAMMA DEM parameter file
    spacing: int or float
        the target pixel spacing in meters
    
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
    ovs_lat = post_north / spacing
    ovs_lon = post_east / spacing
    return ovs_lat, ovs_lon


def pixel_area_wrap(image, namespace, lut, exist_ok=False,
                    logpath=None, outdir=None, shellscript=None):
    """
    helper function for computing pixel_area files in function geocode.

    Parameters
    ----------
    image: str
        the reference SAR image
    namespace: pyroSAR.gamma.auxil.Namespace
        an object collecting all output file names
    lut: str
        the name of the lookup table
    exist_ok: bool
        allow existing output files and do not create new ones?
    logpath: str
        a directory to write command logfiles to
    outdir: str
        the directory to execute the command in
    shellscript: str
        a file to write the GAMMA commands to in shell format

    Returns
    -------

    """
    image_par = ISPPar(image + '.par')
    
    if namespace.isappreciated('gs_ratio'):
        namespace.appreciate(['pix_area_sigma0', 'pix_area_gamma0'])
    
    pixel_area_args = {'MLI_par': image + '.par',
                       'DEM_par': namespace.dem_seg_geo + '.par',
                       'DEM': namespace.dem_seg_geo,
                       'lookup_table': lut,
                       'ls_map': namespace.ls_map_geo,
                       'inc_map': namespace.inc_geo,
                       'pix_sigma0': namespace.pix_area_sigma0,
                       'pix_gamma0': namespace.pix_area_gamma0,
                       'logpath': logpath,
                       'outdir': outdir,
                       'shellscript': shellscript}
    
    radcal_mli_args = {'MLI': image,
                       'MLI_par': image + '.par',
                       'OFF_par': '-',
                       'CMLI': image + '_cal',
                       'refarea_flag': 1,  # calculate sigma0, scale area by sin(inc_ang)/sin(ref_inc_ang)
                       'pix_area': namespace.pix_ellip_sigma0,
                       'logpath': logpath,
                       'outdir': outdir,
                       'shellscript': shellscript}
    
    # newer versions of GAMMA enable creating the ratio of ellipsoid-based
    # pixel area and DEM-facet pixel area directly with command pixel_area
    if hasarg(diff.pixel_area, 'sig2gam_ratio'):
        namespace.appreciate(['pix_ratio'])
        pixel_area_args['sig2gam_ratio'] = namespace.pix_ratio
        if do_execute(pixel_area_args, ['pix_sigma0', 'pix_gamma0', 'sig2gam_ratio'], exist_ok):
            diff.pixel_area(**pixel_area_args)
        
        if namespace.isappreciated('pix_ellip_sigma0'):
            if do_execute(radcal_mli_args, ['pix_area'], exist_ok):
                isp.radcal_MLI(**radcal_mli_args)
                par2hdr(image + '.par', image + '_cal.hdr')
    else:
        # sigma0 = MLI * ellip_pix_sigma0 / pix_area_sigma0
        # gamma0 = MLI * ellip_pix_sigma0 / pix_area_gamma0
        namespace.appreciate(['pix_area_gamma0', 'pix_ellip_sigma0', 'pix_ratio'])
        pixel_area_args['pix_gamma0'] = namespace.pix_area_gamma0
        radcal_mli_args['pix_area'] = namespace.pix_ellip_sigma0
        
        # actual illuminated area as obtained from integrating DEM-facets (pix_area_sigma0 | pix_area_gamma0)
        if do_execute(pixel_area_args, ['pix_sigma0', 'pix_gamma0'], exist_ok):
            diff.pixel_area(**pixel_area_args)
        
        # ellipsoid-based pixel area (ellip_pix_sigma0)
        if do_execute(radcal_mli_args, ['pix_area'], exist_ok):
            isp.radcal_MLI(**radcal_mli_args)
            par2hdr(image + '.par', image + '_cal.hdr')
        
        if os.path.isfile(image + '.hdr'):
            for item in ['pix_area_sigma0', 'pix_area_gamma0', 'pix_ellip_sigma0']:
                if namespace.isappreciated(item):
                    hdr_out = namespace[item] + '.hdr'
                    c1 = not os.path.isfile(hdr_out)
                    c2 = os.path.isfile(hdr_out) and not exist_ok
                    if c1 or c2:
                        shutil.copy(src=image + '.hdr', dst=hdr_out)
        
        # ratio of ellipsoid-based pixel area and DEM-facet pixel area
        c1 = not os.path.isfile(namespace.pix_ratio)
        c2 = os.path.isfile(namespace.pix_ratio) and not exist_ok
        if c1 or c2:
            if 'lat' in locals():
                lat.ratio(d1=namespace.pix_ellip_sigma0,
                          d2=namespace.pix_area_gamma0,
                          ratio=namespace.pix_ratio,
                          width=image_par.range_samples,
                          bx=1,
                          by=1,
                          logpath=logpath,
                          outdir=outdir,
                          shellscript=shellscript)
            else:
                for item in ['pix_area_gamma0', 'pix_ellip_sigma0']:
                    par2hdr(image + '.par', namespace[item] + '.hdr')
                lat_ratio(data_in1=namespace.pix_ellip_sigma0,
                          data_in2=namespace.pix_area_gamma0,
                          data_out=namespace.pix_ratio)
    
    if namespace.isappreciated('gs_ratio'):
        c1 = not os.path.isfile(namespace.gs_ratio)
        c2 = os.path.isfile(namespace.gs_ratio) and not exist_ok
        if c1 or c2:
            if 'lat' in locals():
                lat.ratio(d1=namespace.pix_area_gamma0,
                          d2=namespace.pix_area_sigma0,
                          ratio=namespace.gs_ratio,
                          width=image_par.range_samples,
                          bx=1,
                          by=1,
                          logpath=logpath,
                          outdir=outdir,
                          shellscript=shellscript)
            else:
                for item in ['pix_area_gamma0', 'pix_area_sigma0']:
                    par2hdr(image + '.par', namespace[item] + '.hdr')
                lat_ratio(data_in1=namespace.pix_area_gamma0,
                          data_in2=namespace.pix_area_sigma0,
                          data_out=namespace.gs_ratio)
    
    for item in ['pix_area_sigma0', 'pix_area_gamma0',
                 'pix_ratio', 'pix_ellip_sigma0', 'gs_ratio']:
        if namespace.isappreciated(item):
            hdr_out = namespace[item] + '.hdr'
            c1 = not os.path.isfile(item)
            c2 = os.path.isfile(hdr_out) and not exist_ok
            if c1 or c2:
                par2hdr(image + '.par', hdr_out)


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
    
    new = 'SLC_mosaic_ScanSAR'
    old = 'SLC_mosaic_S1_TOPS'
    slc_mosaic = new if hasattr(isp, new) else old
    getattr(isp, slc_mosaic)(SLC_tab=tab_out,
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
