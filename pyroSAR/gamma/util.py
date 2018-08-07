#!/usr/bin/env python
##############################################################
# universal core routines for processing SAR images with GAMMA
# John Truckenbrodt 2014-2018
##############################################################

"""
This module is intended as a set of generalized processing routines for modularized GAMMA work flows.
The function parametrization is intended to be applicable to any kind of situation and input data set. Thus, instead of choosing a specific parametrization for the data at hand,
core parameters are iterated over a set of values in order to find the one best suited for the task.
The approach of the single routines is likely to still have drawbacks and might fail in certain situations. Testing and suggestions on improvements are very welcome.
"""
import sys

if sys.version_info >= (3, 0):
    from urllib.error import URLError
else:
    from urllib2 import URLError

from osgeo import ogr

from ..spatial import envi
from ..drivers import *
from ..spatial import haversine

from ..ancillary import union, finder
from . import ISPPar, Namespace, process

ogr.UseExceptions()


def calibrate(id, directory, replace=False):

    if isinstance(id, CEOS_PSR):
        for image in id.getGammaImages(directory):
            if image.endswith('_slc'):
                process(
                    ['radcal_SLC', image, image + '.par', image + '_cal', image + '_cal.par',
                     '-', '-', '-', '-', '-', '-', id.meta['k_dB']])
                envi.hdr(image + '_cal.par')

    elif isinstance(id, ESA):
        k_db = {'ASAR': 55., 'ERS1': 58.24, 'ERS2': 59.75}[id.sensor]
        inc_ref = 90. if id.sensor == 'ASAR' else 23.
        candidates = [x for x in id.getGammaImages(directory) if re.search('_pri$', x)]
        for image in candidates:
            out = image.replace('pri', 'grd')
            process(['radcal_PRI', image, image + '.par', out, out + '.par', k_db, inc_ref])
            envi.hdr(out + '.par')
            if replace:
                for item in [image, image + '.par', image + '.hdr']:
                    if os.path.isfile(item):
                        os.remove(item)

    elif isinstance(id, SAFE):
        print('calibration already performed during import')

    else:
        raise NotImplementedError('calibration for class {} is not implemented yet'.format(type(id).__name__))


def convert2gamma(id, directory, S1_noiseremoval=True):
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
                    process(['par_ESA_ERS', lea, outname + '.par', dat, outname], inlist=[title])
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
                process(['par_EORC_PALSAR', id.file, outname + '.par', image, outname])
            else:
                outname_base = '{}_{}_mli_geo'.format(id.outname_base(), polarization)
                outname = os.path.join(directory, outname_base)
                process(
                    ['par_EORC_PALSAR_geo', id.file, outname + '.par', outname + '_dem.par', image, outname])
            envi.hdr(outname + '.par')

    elif isinstance(id, ESA):
        """
        the command par_ASAR also accepts a K_dB argument for calibration in which case the resulting image names will carry the suffix GRD;
        this is not implemented here but instead in function calibrate
        """
        outname = os.path.join(directory, id.outname_base())
        if not id.is_processed(directory):
            process(['par_ASAR', os.path.basename(id.file), outname], os.path.dirname(id.file))
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
                    envi.hdr(newname)
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

            if product == 'slc':
                swath = match.group('swath').upper()
                name = name.replace('{:_<{l}}'.format(id.acquisition_mode, l=len(swath)), swath)
                cmd = ['par_S1_SLC', tiff, xml_ann, xml_cal, xml_noise, name + '.par', name, name + '.tops_par']
            else:
                cmd = ['par_S1_GRD', tiff, xml_ann, xml_cal, xml_noise, name + '.par', name]

            process(cmd)
            envi.hdr(name + '.par')

    elif isinstance(id, TSX):
        images = id.findfiles(id.pattern_ds)
        pattern = re.compile(id.pattern_ds)
        for image in images:
            pol = pattern.match(os.path.basename(image)).group('pol')
            outname = os.path.join(directory, id.outname_base() + '_' + pol)
            if id.product == 'SSC':
                outname += '_slc'
                process(['par_TX_SLC', id.file, image, outname + '.par', outname, pol])
            elif id.product == 'MGD':
                outname += '_mli'
                process(['par_TX_GRD', id.file, image, outname + '.par', outname, pol])
            else:
                outname += '_mli_geo'
                process(['par_TX_geo', id.file, image, outname + '.par', outname + '_dem.par', outname, pol])
            envi.hdr(outname + '.par')

    else:
        raise NotImplementedError('conversion for class {} is not implemented yet'.format(type(id).__name__))


def correctOSV(id, osvdir=None, logpath=None, osvType='POE'):
    """
    correct GAMMA parameter files with orbit state vector information from dedicated OSV files
    :param id: a pyroSAR.ID object
    :param osvdir: the directory of OSV files; subdirectories POEORB and RESORB are created automatically
    :param logpath: a path to write logfiles to
    :param osvType: the type of orbit file either 'POE', 'RES' or a list of both
    :return: None
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
    par = ISPPar(images[0] + '.par')
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
        process(['S1_OPOD_vec', image + '.par', osvfile], logpath=logpath)
    # else:
    #     raise NotImplementedError('OSV refinement for class {} is not implemented yet'.format(type(id).__name__))


def geocode(scene, dem, tempdir, outdir, targetres, scaling='linear', func_geoback=2,
            func_interp=0, nodata=(0, -99), sarSimCC=False, osvdir=None, allow_RES_OSV=False, cleanup=True):
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
        Otherwise the function will raise an error of no POE file exists
    cleanup: bool
        should all files written to the temporary directory during function execution be deleted after processing?

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

    logdir = os.path.join(scene.scene, 'logfiles')
    if not os.path.isdir(logdir):
        os.makedirs(logdir)

    if scene.sensor in ['S1A', 'S1B']:
        print('removing border noise..')
        scene.removeGRDBorderNoise()

    print('converting scene to GAMMA format..')
    convert2gamma(scene, scene.scene)

    if scene.sensor in ['S1A', 'S1B']:
        print('updating orbit state vectors..')
        if allow_RES_OSV:
            osvtype = ['POE', 'RES']
        else:
            osvtype = 'POE'
        try:
            correctOSV(id=scene, osvdir=osvdir, logpath=logdir, osvType=osvtype)
        except RuntimeError:
            return

    calibrate(scene, scene.scene)

    images = [x for x in scene.getGammaImages(scene.scene) if x.endswith('_grd') or x.endswith('_slc_cal')]

    print('multilooking..')
    for image in images:
        multilook(image, image + '_mli', targetres)

    images = [x + '_mli' for x in images]

    master = images[0]

    # create output names for files to be written
    # appreciated files will be written
    # depreciated files will be set to '-' in the GAMMA funtion call and are thus not written
    n = Namespace(scene.scene, scene.outname_base())
    n.appreciate(['dem_seg', 'lut_coarse', 'lut_fine', 'pix', 'ccp', 'inc', 'ls_map'])
    n.depreciate(['sim_map', 'u', 'v', 'psi'])

    # if sarSimCC:
    #     n.appreciate(['ls_map'])

    ovs_lat, ovs_lon = ovs(dem + '.par', targetres)

    path_log = os.path.join(scene.scene, 'logfiles')
    if not os.path.isdir(path_log):
        os.makedirs(path_log)

    master_par = ISPPar(master + '.par')

    gc_map_args = [dem + '.par', dem, n.dem_seg + '.par', n.dem_seg, n.lut_coarse,
                   ovs_lat, ovs_lon, n.sim_map, n.u, n.v, n.inc, n.psi, n.pix, n.ls_map,
                   8, func_interp]

    print('SAR image simulation from DEM..')
    if master_par.image_geometry == 'GROUND_RANGE':
        process(['gc_map_grd', master + '.par'] + gc_map_args, logpath=path_log)
    else:
        process(['gc_map', master + '.par', '-'] + gc_map_args, logpath=path_log)

    for item in ['dem_seg', 'sim_map', 'u', 'v', 'psi', 'pix', 'inc']:
        if n.isappreciated(item):
            envi.hdr(n.dem_seg + '.par', n.get(item) + '.hdr')

    sim_width = ISPPar(n.dem_seg + '.par').width

    if sarSimCC:
        raise IOError('geocoding with cross correlation offset refinement is still in the making. Please stay tuned...')
    else:
        lut_final = n.lut_coarse

    ######################################################################
    # normalization and backward geocoding approach 1 ####################
    ######################################################################
    print('geocoding and normalization..')
    for image in images:
        process(['geocode_back', image, master_par.range_samples, lut_final, image + '_geo', sim_width, '-', func_geoback], logpath=path_log)
        process(['product', image + '_geo', n.pix, image + '_geo_pan', sim_width, 1, 1, 0], logpath=path_log)
        process(['lin_comb', 1, image + '_geo_pan', 0, math.cos(math.radians(master_par.incidence_angle)), image + '_geo_pan_flat', sim_width], logpath=path_log)
        process(['sigma2gamma', image + '_geo_pan_flat', n.inc, image + '_geo_norm', sim_width], logpath=path_log)
        envi.hdr(n.dem_seg + '.par', image + '_geo_norm.hdr')
    ######################################################################
    # normalization and backward geocoding approach 2 ####################
    ######################################################################
    # process(['pixel_area', master+'.par', dem_seg+'.par', dem_seg, lut_fine, ls_map, inc, pixel_area_fine], logpath=path_log)
    # process(['radcal_MLI', master, master+'.par', '-', master+'_cal', '-', 0, 0, 1, 0.0, '-', ellipse_pixel_area], logpath=path_log)
    # process(['ratio', ellipse_pixel_area, pixel_area_fine, ratio_sigma0, master_par.range_samples, 1, 1], logpath=path_log)
    #
    # for image in images:
    #     process(['product', image, ratio_sigma0, image+'_pan', master_par.range_samples, 1, 1], logpath=path_log)
    #     process(['geocode_back', image+'_pan', master_par.range_samples, lut_fine, image+'_pan_geo', sim_width, 0, func_geoback], logpath=path_log)
    #     process(['lin_comb', 1, image+'_pan_geo', 0, math.cos(math.radians(master_par.incidence_angle)), image+'_pan_geo_flat', sim_width], logpath=path_log)
    #     process(['sigma2gamma', image+'_pan_geo_flat', inc, image+'_geo_norm', sim_width], logpath=path_log)
    #     envi.hdr(dem_seg+'.par', image+'_geo_norm.hdr')
    ######################################################################
    print('conversion to (dB and) geotiff..')
    for image in images:
        for scale in scaling:
            if scale == 'db':
                process(['linear_to_dB', image + '_geo_norm', image + '_geo_norm_db', sim_width, 0, -99], logpath=path_log)
                envi.hdr(n.dem_seg + '.par', image + '_geo_norm_db.hdr')
                nodata_out = nodata[1]
            else:
                nodata_out = nodata[0]
            suffix = {'linear': '', 'db': '_db'}[scale]
            infile = image + '_geo_norm{}'.format(suffix)
            outfile = os.path.join(outdir, os.path.basename(image) + '_geo_norm{}.tif'.format(suffix))

            process(['data2geotiff', n.dem_seg + '.par', infile, 2, outfile, nodata_out], logpath=path_log)

    if scene.sensor in ['S1A', 'S1B']:
        shutil.copyfile(os.path.join(scene.scene, 'manifest.safe'),
                        os.path.join(outdir, scene.outname_base() + '_manifest.safe'))
    if cleanup:
        print('cleaning up temporary files..')
        shutil.rmtree(scene.scene)


def ovs(parfile, targetres):
    """
    compute DEM oversampling factors for a target resolution in meters

    :param parfile: a GAMMA DEM parameter file
    :param targetres: the target resolution in meters
    :return: two oversampling factors for latitude and longitude
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


def multilook(infile, outfile, targetres):
    """
    multilooking of SLC and MLI images

    if the image is in slant range the ground range resolution is computed by dividing the range pixel spacing by
    the sine of the incidence angle

    the looks in range and azimuth are chosen to approximate the target resolution by rounding the ratio between
    target resolution and ground range/azimuth pixel spacing to the nearest integer

    an ENVI HDR parameter file is automatically written for better handling on other software

    :param infile a SAR image in GAMMA format with a parameter file of name <infile>.par
    :param outfile the name of the output GAMMA file
    :param targetres: the target resolution in ground range


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

    if par.image_format in ['SCOMPLEX', 'FCOMPLEX']:
        # multilooking for SLC images
        process(['multi_look', infile, infile + '.par', outfile, outfile + '.par', rlks, azlks])
    else:
        # multilooking for MLI images
        process(['multi_look_MLI', infile, infile + '.par', outfile, outfile + '.par', rlks, azlks])
    envi.hdr(outfile + '.par')


def S1_deburst(burst1, burst2, burst3, name_out, rlks=5, azlks=1, replace=False, path_log=None):
    """
    debursting of S1 SLC imagery in GAMMA
    the procedure consists of two steps. First antenna pattern deramping and then mosaicing of the single deramped bursts
    for mosaicing, the burst boundaries are calculated from the number of looks in range (rlks) and azimuth (azlks), in this case 5 range looks and 1 azimuth looks.
    Alternately 10 range looks and 2 azimuth looks could be used.
    if replace is set to True, the original files will be deleted
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
    process(['SLC_deramp_S1_TOPS', tab_in, tab_out, 0, 0], logpath=path_log)
    process(['SLC_mosaic_S1_TOPS', tab_out, name_out, name_out + '.par', rlks, azlks], logpath=path_log)

    if replace:
        for item in [burst1, burst2, burst3]:
            for subitem in [item + x for x in ['', '.par', '.tops_par']]:
                os.remove(subitem)
    for item in [burst1, burst2, burst3]:
        for subitem in [item + x for x in ['_drp', '_drp.par', '_drp.tops_par']]:
            os.remove(subitem)
    os.remove(tab_in)
    os.remove(tab_out)