#!/usr/bin/env python
##############################################################
# universal core routines for processing SAR images with GAMMA
# John Truckenbrodt 2014-2018
##############################################################

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
import sys
import math
import shutil
from datetime import datetime

if sys.version_info >= (3, 0):
    from urllib.error import URLError
else:
    from urllib2 import URLError

from osgeo import ogr

from spatialist import haversine
from spatialist.ancillary import union, finder

from ..S1 import OSV
from ..drivers import ID, CEOS_ERS, CEOS_PSR, ESA, SAFE, TSX, identify
from . import ISPPar, Namespace, par2hdr

try:
    from .api import diff, disp, isp, lat
except ImportError:
    pass

ogr.UseExceptions()


def calibrate(id, directory, replace=False, logpath=None, outdir=None, shellscript=None):
    """
    
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
    if isinstance(id, CEOS_PSR):
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
    
    elif isinstance(id, ESA):
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
    
    elif isinstance(id, SAFE):
        print('calibration already performed during import')
    
    else:
        raise NotImplementedError('calibration for class {} is not implemented yet'.format(type(id).__name__))


def convert2gamma(id, directory, S1_noiseremoval=True, logpath=None, outdir=None, shellscript=None):
    """
    general function for converting SAR images to GAMMA format

    Parameters
    ----------
    id: ~pyroSAR.drivers.ID
        an SAR scene object of type pyroSAR.ID or any subclass
    directory: str
        the output directory for the converted images
    S1_noiseremoval: bool
        only Sentinel-1: should noise removal be applied to the image?
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
    
    if isinstance(id, CEOS_ERS):
        if id.sensor in ['ERS1', 'ERS2']:
            if id.product == 'SLC' and id.meta['proc_system'] in ['PGS-ERS', 'VMP-ERS', 'SPF-ERS']:
                basename = '{}_{}_{}'.format(id.outname_base(), id.polarizations[0], id.product.lower())
                outname = os.path.join(directory, basename)
                if not os.path.isfile(outname):
                    lea = id.findfiles('LEA_01.001')[0]
                    dat = id.findfiles('DAT_01.001')[0]
                    title = re.sub('\.PS$', '', os.path.basename(id.file))
                    
                    isp.par_ESA_ERS(CEOS_SAR_leader=lea,
                                    SLC_par=outname + '.par',
                                    CEOS_DAT=dat,
                                    SLC=outname,
                                    inlist=[title],
                                    logpath=logpath,
                                    outdir=outdir,
                                    shellscript=shellscript)
                else:
                    print('scene already converted')
            else:
                raise NotImplementedError('ERS {} product of {} processor in CEOS format not implemented yet'
                                          .format(id.product, id.meta['proc_system']))
        else:
            raise NotImplementedError('sensor {} in CEOS format not implemented yet'.format(id.sensor))
    
    elif isinstance(id, CEOS_PSR):
        images = id.findfiles('^IMG-')
        if id.product == '1.0':
            raise RuntimeError('PALSAR level 1.0 products are not supported')
        for image in images:
            polarization = re.search('[HV]{2}', os.path.basename(image)).group(0)
            if id.product == '1.1':
                outname_base = '{}_{}_slc'.format(id.outname_base(), polarization)
                outname = os.path.join(directory, outname_base)
                
                isp.par_EORC_PALSAR(CEOS_leader=id.file,
                                    SLC_par=outname + '.par',
                                    CEOS_data=image,
                                    SLC=outname,
                                    logpath=logpath,
                                    outdir=outdir,
                                    shellscript=shellscript)
            else:
                outname_base = '{}_{}_mli_geo'.format(id.outname_base(), polarization)
                outname = os.path.join(directory, outname_base)
                
                isp.par_EORC_PALSAR_geo(CEOS_leader=id.file,
                                        MLI_par=outname + '.par',
                                        DEM_par=outname + '_dem.par',
                                        CEOS_data=image,
                                        MLI=outname,
                                        logpath=logpath,
                                        outdir=outdir,
                                        shellscript=shellscript)
            par2hdr(outname + '.par', outname + '.hdr')
    
    elif isinstance(id, ESA):
        """
        the command par_ASAR also accepts a K_dB argument for calibration in which case the resulting image names will carry the suffix GRD;
        this is not implemented here but instead in function calibrate
        """
        outname = os.path.join(directory, id.outname_base())
        if not id.is_processed(directory):
            
            isp.par_ASAR(ASAR_ERS_file=os.path.basename(id.file),
                         output_name=outname,
                         outdir=os.path.dirname(id.file),
                         logpath=logpath,
                         shellscript=shellscript)
            
            os.remove(outname + '.hdr')
            for item in finder(directory, [os.path.basename(outname)], regex=True):
                ext = '.par' if item.endswith('.par') else ''
                base = os.path.basename(item).strip(ext)
                base = base.replace('.', '_')
                base = base.replace('PRI', 'pri')
                base = base.replace('SLC', 'slc')
                newname = os.path.join(directory, base + ext)
                os.rename(item, newname)
                if newname.endswith('.par'):
                    par2hdr(newname, newname.replace('.par', '.hdr'))
        else:
            raise IOError('scene already processed')
    
    elif isinstance(id, SAFE):
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
            if (S1_noiseremoval and product == 'slc') or (not S1_noiseremoval and product == 'grd'):
                xml_noise = os.path.join(id.scene, 'annotation', 'calibration', 'noise-' + base)
            else:
                xml_noise = '-'
            
            fields = (id.outname_base(),
                      match.group('pol').upper(),
                      product)
            name = os.path.join(directory, '_'.join(fields))
            
            pars = {'GeoTIFF': tiff,
                    'annotation_XML': xml_ann,
                    'calibration_XML': xml_cal,
                    'noise_XML': xml_noise,
                    'logpath': logpath,
                    'shellscript': shellscript,
                    'outdir': outdir}
            
            if product == 'slc':
                swath = match.group('swath').upper()
                name = name.replace('{:_<{length}}'.format(id.acquisition_mode, length=len(swath)), swath)
                pars['SLC'] = name
                pars['SLC_par'] = name + '.par'
                pars['TOPS_par'] = name + '.tops_par'
                isp.par_S1_SLC(**pars)
            else:
                pars['MLI'] = name
                pars['MLI_par'] = name + '.par'
                isp.par_S1_GRD(**pars)
            
            par2hdr(name + '.par', name + '.hdr')
    
    elif isinstance(id, TSX):
        images = id.findfiles(id.pattern_ds)
        pattern = re.compile(id.pattern_ds)
        for image in images:
            pol = pattern.match(os.path.basename(image)).group('pol')
            outname = os.path.join(directory, id.outname_base() + '_' + pol)
            
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
                isp.par_TX_SLC(**pars)
            
            elif id.product == 'MGD':
                outname += '_mli'
                pars['GeoTIFF'] = image
                pars['GRD_par'] = outname + '.par'
                pars['GRD'] = outname
                isp.par_TX_GRD(**pars)
            
            elif id.product in ['GEC', 'EEC']:
                outname += '_mli_geo'
                pars['GeoTIFF'] = image
                pars['MLI_par'] = outname + '.par'
                pars['DEM_par'] = outname + '_dem.par'
                pars['GEO'] = outname
                isp.par_TX_geo(**pars)
            else:
                raise RuntimeError('unknown product: {}'.format(id.product))
            
            par2hdr(outname + '.par', outname + '.hdr')
    else:
        raise NotImplementedError('conversion for class {} is not implemented yet'.format(type(id).__name__))


def correctOSV(id, osvdir=None, osvType='POE', logpath=None, outdir=None, shellscript=None):
    """
    correct GAMMA parameter files with orbit state vector information from dedicated OSV files
    
    Parameters
    ----------
    id: ~pyroSAR.drivers.ID
        the scene to be corrected
    osvdir: str
        the directory of OSV files; subdirectories POEORB and RESORB are created automatically
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
    
    if id.sensor not in ['S1A', 'S1B']:
        raise IOError('this method is currently only available for Sentinel-1. Please stay tuned...')
    
    if not os.path.isdir(logpath):
        os.makedirs(logpath)
    
    if osvdir is None:
        osvdir = os.path.join(id.scene, 'osv')
        if not os.path.isdir(osvdir):
            os.makedirs(osvdir)
        try:
            id.getOSV(osvdir, osvType)
        except URLError:
            print('..no internet access')
    
    images = id.getGammaImages(id.scene)
    # read parameter file entries int object
    with ISPPar(images[0] + '.par') as par:
        # extract acquisition time stamp
        timestamp = datetime(*map(int, par.date)).strftime('%Y%m%dT%H%M%S')
    
    # find an OSV file matching the time stamp and defined OSV type(s)
    with OSV(osvdir) as osv:
        osvfile = osv.match(timestamp, osvType)
    if not osvfile:
        raise RuntimeError('no Orbit State Vector file found')
    
    # update the GAMMA parameter file with the selected orbit state vectors
    print('correcting state vectors with file {}'.format(osvfile))
    for image in images:
        isp.S1_OPOD_vec(SLC_par=image + '.par',
                        OPOD=osvfile,
                        logpath=logpath,
                        outdir=outdir,
                        shellscript=shellscript)


def geocode(scene, dem, tempdir, outdir, targetres, scaling='linear', func_geoback=2,
            func_interp=0, nodata=(0, -99), sarSimCC=False, osvdir=None, allow_RES_OSV=False,
            cleanup=True, normalization_method=2):
    """
    general function for geocoding SAR images with GAMMA
    
    Parameters
    ----------
    scene: str or ~pyroSAR.drivers.ID
        the SAR scene to be processed
    dem: str
        the reference DEM in GAMMA format
    tempdir: str
        a temporary directory for writing intermediate files
    outdir: str
        the directory for the final GeoTiff output files
    targetres: int
        the target resolution in meters
    scaling: {'linear', 'db'} or list
        the value scaling of the backscatter values; either 'linear', 'db' or a list of both, i.e. ['linear', 'db']
    func_geoback: {0, 1, 2, 3}
        backward geocoding interpolation mode (see GAMMA command geocode_back)
         * 0: nearest-neighbor
         * 1: bicubic spline
         * 2: bicubic-log spline, interpolates log(data)
         * 3: bicubic-sqrt spline, interpolates sqrt(data)

        NOTE: bicubic-log spline and bicubic-sqrt spline modes should only be used with non-negative data!
    func_interp: {0, 1, 2, 3}
        output lookup table values in regions of layover, shadow, or DEM gaps (see GAMMA command gc_map)
         * 0: set to (0., 0.)
         * 1: linear interpolation across these regions
         * 2: actual value
         * 3: nn-thinned
    nodata: tuple
        the nodata values for the output files; defined as a tuple with two values, the first for linear,
        the second for logarithmic scaling
    sarSimCC: bool
        perform geocoding with SAR simulation cross correlation?
        If False, geocoding is performed with the Range-Doppler approach using orbit state vectors
    osvdir: str
        a directory for Orbit State Vector files;
        this is currently only used by for Sentinel-1 where two subdirectories POEORB and RESORB are created;
        if set to None, a subdirectory OSV is created in the directory of the unpacked scene.
    allow_RES_OSV: bool
        also allow the less accurate RES orbit files to be used?
        Otherwise the function will raise an error if no POE file exists
    cleanup: bool
        should all files written to the temporary directory during function execution be deleted after processing?
    normalization_method: {1, 2}
        the topographic normalization approach to be used
         * 1: first geocoding, then terrain flattening
         * 2: first terrain flattening, then geocoding; see `Small 2011 <https://doi.org/10.1109/Tgrs.2011.2120616>`_
    
    Returns
    -------
    
    Note
    ----
    intermediate output files (named <master_MLI>_<suffix>):
     * dem_seg: dem subsetted to the extent of the SAR image
     * lut: rough geocoding lookup table
     * lut_fine: fine geocoding lookup table
     * sim_map: simulated SAR backscatter image in DEM geometry
     * sim_sar: simulated SAR backscatter image in SAR geometry
     * u: zenith angle of surface normal vector n (angle between z and n)
     * v: orientation angle of n (between x and projection of n in xy plane)
     * inc: local incidence angle (between surface normal and look vector)
     * psi: projection angle (between surface normal and image plane normal)
     * pix: pixel area normalization factor
     * ls_map: layover and shadow map (in map projection)
     * diffpar: ISP offset/interferogram parameter file
     * offs: offset estimates (fcomplex)
     * coffs: culled range and azimuth offset estimates (fcomplex)
     * coffsets: culled offset estimates and cross correlation values (text format)
     * ccp: cross-correlation of each patch (0.0->1.0) (float)
    
    """
    
    scene = scene if isinstance(scene, ID) else identify(scene)
    
    shellscript = os.path.join(outdir, scene.outname_base() + '_commands.sh')
    
    if scene.sensor not in ['S1A', 'S1B']:
        raise IOError('this method is currently only available for Sentinel-1. Please stay tuned...')
    
    if sarSimCC:
        raise IOError('geocoding with cross correlation offset refinement is still in the making. Please stay tuned...')
    
    for dir in [tempdir, outdir]:
        if not os.path.isdir(dir):
            os.makedirs(dir)
    
    if scene.is_processed(outdir):
        print('scene {} already processed'.format(scene.outname_base()))
        return
    
    scaling = [scaling] if isinstance(scaling, str) else scaling if isinstance(scaling, list) else []
    scaling = union(scaling, ['db', 'linear'])
    if len(scaling) == 0:
        raise IOError('wrong input type for parameter scaling')
    
    if scene.compression is not None:
        print('unpacking scene..')
        try:
            scene.unpack(tempdir)
        except RuntimeError:
            print('scene was attempted to be processed before, exiting')
            return
    else:
        scene.scene = os.path.join(tempdir, os.path.basename(scene.file))
        os.makedirs(scene.scene)
    
    path_log = os.path.join(scene.scene, 'logfiles')
    if not os.path.isdir(path_log):
        os.makedirs(path_log)
    
    if scene.sensor in ['S1A', 'S1B']:
        print('removing border noise..')
        scene.removeGRDBorderNoise()
    
    print('converting scene to GAMMA format..')
    convert2gamma(scene, scene.scene, logpath=path_log, outdir=scene.scene, shellscript=shellscript)
    
    if scene.sensor in ['S1A', 'S1B']:
        print('updating orbit state vectors..')
        if allow_RES_OSV:
            osvtype = ['POE', 'RES']
        else:
            osvtype = 'POE'
        try:
            correctOSV(id=scene, osvdir=osvdir, osvType=osvtype,
                       logpath=path_log, outdir=scene.scene, shellscript=shellscript)
        except RuntimeError:
            return
    
    calibrate(scene, scene.scene, logpath=path_log, outdir=scene.scene, shellscript=shellscript)
    
    images = [x for x in scene.getGammaImages(scene.scene) if x.endswith('_grd') or x.endswith('_slc_cal')]
    
    print('multilooking..')
    for image in images:
        multilook(infile=image, outfile=image + '_mli', targetres=targetres,
                  logpath=path_log, outdir=scene.scene, shellscript=shellscript)
    
    images = [x + '_mli' for x in images]
    
    master = images[0]
    
    # create output names for files to be written
    # appreciated files will be written
    # depreciated files will be set to '-' in the GAMMA function call and are thus not written
    n = Namespace(scene.scene, scene.outname_base())
    n.appreciate(['dem_seg', 'lut_coarse', 'lut_fine', 'pix', 'ccp', 'inc', 'ls_map'])
    n.depreciate(['sim_map', 'u', 'v', 'psi'])
    
    # if sarSimCC:
    #     n.appreciate(['ls_map'])
    
    ovs_lat, ovs_lon = ovs(dem + '.par', targetres)
    
    master_par = ISPPar(master + '.par')
    
    gc_map_args = {'DEM_par': dem + '.par',
                   'DEM': dem,
                   'DEM_seg_par': n.dem_seg + '.par',
                   'DEM_seg': n.dem_seg,
                   'lookup_table': n.lut_coarse,
                   'lat_ovr': ovs_lat,
                   'lon_ovr': ovs_lon,
                   'sim_sar': n.sim_map,
                   'u': n.u,
                   'v': n.v,
                   'inc': n.inc,
                   'psi': n.psi,
                   'pix': n.pix,
                   'ls_map': n.ls_map,
                   'frame': 8,
                   'ls_mode': func_interp,
                   'logpath': path_log,
                   'shellscript': shellscript,
                   'outdir': scene.scene}
    
    print('SAR image simulation from DEM..')
    if master_par.image_geometry == 'GROUND_RANGE':
        gc_map_args.update({'GRD_par': master + '.par'})
        diff.gc_map_grd(**gc_map_args)
    else:
        gc_map_args.update({'MLI_par': master + '.par',
                            'OFF_par': '-'})
        diff.gc_map(**gc_map_args)
    
    for item in ['dem_seg', 'sim_map', 'u', 'v', 'psi', 'pix', 'inc']:
        if n.isappreciated(item):
            par2hdr(n.dem_seg + '.par', n.get(item) + '.hdr')
    
    sim_width = ISPPar(n.dem_seg + '.par').width
    
    if sarSimCC:
        raise IOError('geocoding with cross correlation offset refinement is still in the making. Please stay tuned...')
    else:
        lut_final = n.lut_coarse
    
    ######################################################################
    # normalization and backward geocoding approach 1 ####################
    ######################################################################
    print('geocoding and normalization..')
    if normalization_method == 1:
        method_suffix = 'geo_norm'
        for image in images:
            diff.geocode_back(data_in=image,
                              width_in=master_par.range_samples,
                              gc_map=lut_final,
                              data_out=image + '_geo',
                              width_out=sim_width,
                              interp_mode=func_geoback,
                              logpath=path_log,
                              outdir=scene.scene,
                              shellscript=shellscript)
            lat.product(data_1=image + '_geo',
                        data_2=n.pix,
                        product=image + '_geo_pan',
                        width=sim_width,
                        bx=1,
                        by=1,
                        logpath=path_log,
                        outdir=scene.scene,
                        shellscript=shellscript)
            lat.lin_comb(files=[image + '_geo_pan'],
                         constant=0,
                         factors=[math.cos(math.radians(master_par.incidence_angle))],
                         f_out=image + '_geo_pan_flat',
                         width=sim_width,
                         logpath=path_log,
                         outdir=scene.scene,
                         shellscript=shellscript)
            lat.sigma2gamma(pwr1=image + '_geo_pan_flat',
                            inc=n.inc,
                            gamma=image + '_{}'.format(method_suffix),
                            width=sim_width,
                            logpath=path_log,
                            outdir=scene.scene,
                            shellscript=shellscript)
            par2hdr(n.dem_seg + '.par', image + '_{}.hdr'.format(method_suffix))
    ######################################################################
    # normalization and backward geocoding approach 2 ####################
    ######################################################################
    elif normalization_method == 2:
        method_suffix = 'norm_geo'
        n.appreciate(['pixel_area_fine', 'ellipse_pixel_area', 'ratio_sigma0'])
        diff.pixel_area(MLI_par=master + '.par',
                        DEM_par=n.dem_seg + '.par',
                        DEM=n.dem_seg,
                        lookup_table=lut_final,
                        ls_map=n.ls_map,
                        inc_map=n.inc,
                        pix_sigma0=n.pixel_area_fine,
                        logpath=path_log,
                        outdir=scene.scene,
                        shellscript=shellscript)
        isp.radcal_MLI(MLI=master,
                       MLI_par=master + '.par',
                       OFF_par='-',
                       CMLI=master + '_cal',
                       refarea_flag=1,
                       pix_area=n.ellipse_pixel_area,
                       logpath=path_log,
                       outdir=scene.scene,
                       shellscript=shellscript)
        lat.ratio(d1=n.ellipse_pixel_area,
                  d2=n.pixel_area_fine,
                  ratio=n.ratio_sigma0,
                  width=master_par.range_samples,
                  bx=1,
                  by=1,
                  logpath=path_log,
                  outdir=scene.scene,
                  shellscript=shellscript)
        for image in images:
            lat.product(data_1=image,
                        data_2=n.ratio_sigma0,
                        product=image + '_pan',
                        width=master_par.range_samples,
                        bx=1,
                        by=1,
                        logpath=path_log,
                        outdir=scene.scene,
                        shellscript=shellscript)
            diff.geocode_back(data_in=image + '_pan',
                              width_in=master_par.range_samples,
                              gc_map=lut_final,
                              data_out=image + '_pan_geo',
                              width_out=sim_width,
                              interp_mode=func_geoback,
                              logpath=path_log,
                              outdir=scene.scene,
                              shellscript=shellscript)
            lat.lin_comb(files=[image + '_pan_geo'],
                         constant=0,
                         factors=[math.cos(math.radians(master_par.incidence_angle))],
                         f_out=image + '_pan_geo_flat',
                         width=sim_width,
                         logpath=path_log,
                         outdir=scene.scene,
                         shellscript=shellscript)
            lat.sigma2gamma(pwr1=image + '_pan_geo_flat',
                            inc=n.inc,
                            gamma=image + '_{}'.format(method_suffix),
                            width=sim_width,
                            logpath=path_log,
                            outdir=scene.scene,
                            shellscript=shellscript)
            par2hdr(n.dem_seg + '.par', image + '_{}.hdr'.format(method_suffix))
    else:
        raise RuntimeError('unknown option for normalization_method')
    ######################################################################
    print('conversion to (dB and) geotiff..')
    for image in images:
        for scale in scaling:
            if scale == 'db':
                nodata_out = nodata[1]
                lat.linear_to_dB(data_in=image + '_{}'.format(method_suffix),
                                 data_out=image + '_{}_db'.format(method_suffix),
                                 width=sim_width,
                                 inverse_flag=0,
                                 null_value=nodata_out,
                                 logpath=path_log,
                                 outdir=scene.scene,
                                 shellscript=shellscript)
                par2hdr(n.dem_seg + '.par', image + '_{}_db.hdr'.format(method_suffix))
            else:
                nodata_out = nodata[0]
            suffix = {'linear': '', 'db': '_db'}[scale]
            infile = image + '_{0}{1}'.format(method_suffix, suffix)
            outfile = os.path.join(outdir, os.path.basename(image) + '_{0}{1}.tif'.format(method_suffix, suffix))
            disp.data2geotiff(DEM_par=n.dem_seg + '.par',
                              data=infile,
                              type=2,
                              GeoTIFF=outfile,
                              nodata=nodata_out,
                              logpath=path_log,
                              outdir=scene.scene,
                              shellscript=shellscript)
    if scene.sensor in ['S1A', 'S1B']:
        shutil.copyfile(os.path.join(scene.scene, 'manifest.safe'),
                        os.path.join(outdir, scene.outname_base() + '_manifest.safe'))
    if cleanup:
        print('cleaning up temporary files..')
        shutil.rmtree(scene.scene)


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


def multilook(infile, outfile, targetres, logpath=None, outdir=None, shellscript=None):
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
    logpath: str or None
        a directory to write command logfiles to
    outdir: str or None
        the directory to execute the command in
    shellscript: str or None
        a file to write the Gamma commands to in shell format

    """
    # read the input parameter file
    par = ISPPar(infile + '.par')
    
    # compute the range looks
    if par.image_geometry == 'SLANT_RANGE':
        # compute the ground range resolution
        groundRangePS = par.range_pixel_spacing / (math.sin(math.radians(par.incidence_angle)))
        rlks = int(round(float(targetres) / groundRangePS))
    else:
        rlks = int(round(float(targetres) / par.range_pixel_spacing))
    # compute the azimuth looks
    azlks = int(round(float(targetres) / par.azimuth_pixel_spacing))
    
    # set the look factors to 1 if they were computed to be 0
    rlks = rlks if rlks > 0 else 1
    azlks = azlks if azlks > 0 else 1
    
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
        isp.multi_look(**pars)
    else:
        # multilooking for MLI images
        pars['MLI_in'] = infile
        pars['MLI_in_par'] = infile + '.par'
        pars['MLI_out'] = outfile
        pars['MLI_out_par'] = outfile + '.par'
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
    
    isp.S1_deramp_S1_TOPS(SLC1_tab=tab_in,
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
