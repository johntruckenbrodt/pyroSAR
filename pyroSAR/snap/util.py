####################################################################
# Convenience functions for SAR image batch processing with ESA SNAP
# John Truckenbrodt, 2016-2019
####################################################################
import os
import shutil
import pyroSAR
from ..ancillary import multilook_factors
from .auxil import parse_recipe, parse_node, gpt, groupbyWorkers

from spatialist import crsConvert, Vector, Raster, bbox, intersect


def geocode(infile, outdir, t_srs=4326, tr=20, polarizations='all', shapefile=None, scaling='dB',
            geocoding_type='Range-Doppler', removeS1BorderNoise=True, removeS1ThermalNoise=True, offset=None,
            externalDEMFile=None, externalDEMNoDataValue=None, externalDEMApplyEGM=True, terrainFlattening=True,
            basename_extensions=None, test=False, export_extra=None, groupsize=2, cleanup=True,
            gpt_exceptions=None, returnWF=False,
            demResamplingMethod='BILINEAR_INTERPOLATION', imgResamplingMethod='BILINEAR_INTERPOLATION'):
    """
    wrapper function for geocoding SAR images using ESA SNAP

    Parameters
    ----------
    infile: str or ~pyroSAR.drivers.ID
        the SAR scene to be processed
    outdir: str
        The directory to write the final files to.
    t_srs: int, str or osr.SpatialReference
        A target geographic reference system in WKT, EPSG, PROJ4 or OPENGIS format.
        See function :func:`spatialist.auxil.crsConvert()` for details.
        Default: `4326 <http://spatialreference.org/ref/epsg/4326/>`_.
    tr: int or float, optional
        The target resolution in meters. Default is 20
    polarizations: list or {'VV', 'HH', 'VH', 'HV', 'all'}, optional
        The polarizations to be processed; can be a string for a single polarization e.g. 'VV' or a list of several
        polarizations e.g. ['VV', 'VH']. Default is 'all'.
    shapefile: str or :py:class:`~spatialist.vector.Vector`, optional
        A vector geometry for subsetting the SAR scene to a test site. Default is None.
    scaling: {'dB', 'db', 'linear'}, optional
        Should the output be in linear or decibel scaling? Default is 'dB'.
    geocoding_type: {'Range-Doppler', 'SAR simulation cross correlation'}, optional
        The type of geocoding applied; can be either 'Range-Doppler' (default) or 'SAR simulation cross correlation'
    removeS1BorderNoise: bool, optional
        Enables removal of S1 GRD border noise (default).
    removeS1ThermalNoise: bool, optional
        Enables removal of S1 thermal noise (default).
    offset: tuple, optional
        A tuple defining offsets for left, right, top and bottom in pixels, e.g. (100, 100, 0, 0); this variable is
        overridden if a shapefile is defined. Default is None.
    externalDEMFile: str or None, optional
        The absolute path to an external DEM file. Default is None.
    externalDEMNoDataValue: int, float or None, optional
        The no data value of the external DEM. If not specified (default) the function will try to read it from the
        specified external DEM.
    externalDEMApplyEGM: bool, optional
        Apply Earth Gravitational Model to external DEM? Default is True.
    terainFlattening: bool
        apply topographic normalization on the data?
    basename_extensions: list of str
        names of additional parameters to append to the basename, e.g. ['orbitNumber_rel']
    test: bool, optional
        If set to True the workflow xml file is only written and not executed. Default is False.
    export_extra: list or None
        a list of image file IDs to be exported to outdir. The following IDs are currently supported:
         - incidenceAngleFromEllipsoid
         - localIncidenceAngle
         - projectedLocalIncidenceAngle
         - DEM
    groupsize: int
        the number of workers executed together in one gpt call
    cleanup: bool
        should all files written to the temporary directory during function execution be deleted after processing?
    gpt_exceptions: dict
        a dictionary to override the configured GPT executable for certain operators;
        each (sub-)workflow containing this operator will be executed with the define executable;
        
         - e.g. ``{'Terrain-Flattening': '/home/user/snap/bin/gpt'}``
    returnWF: bool
        return the full name of the written workflow XML file?
    demResamplingMethod: str
        one of the following:
         - 'NEAREST_NEIGHBOUR'
         - 'BILINEAR_INTERPOLATION'
         - 'CUBIC_CONVOLUTION'
         - 'BISINC_5_POINT_INTERPOLATION'
         - 'BISINC_11_POINT_INTERPOLATION'
         - 'BISINC_21_POINT_INTERPOLATION'
         - 'BICUBIC_INTERPOLATION'
    imgResamplingMethod: str
        the resampling method for geocoding the SAR image; the options are identical to demResamplingMethod
    
    Returns
    -------
    str or None
        either the name of the workflow file if `returnWF == True` or None otherwise

    Note
    ----
    If only one polarization is selected and not extra products are defined the results are directly written to GeoTiff.
    Otherwise the results are first written to a folder containing ENVI files and then transformed to GeoTiff files
    (one for each polarization/extra product).
    If GeoTiff would directly be selected as output format for multiple polarizations then a multilayer GeoTiff
    is written by SNAP which is considered an unfavorable format
    
    
    .. figure:: snap_geocode.png
        :scale: 25%
        :align: center
        
        Workflow diagram for function geocode for processing a Sentinel-1 Ground Range
        Detected (GRD) scene to radiometrically terrain corrected (RTC) backscatter.
        An additional Subset node might be inserted in case a vector geometry is provided.

    Examples
    --------
    geocode a Sentinel-1 scene and export the local incidence angle map with it

    >>> from pyroSAR.snap import geocode
    >>> filename = 'S1A_IW_GRDH_1SDV_20180829T170656_20180829T170721_023464_028DE0_F7BD.zip'
    >>> geocode(infile=filename, outdir='outdir', tr=20, scaling='dB',
    >>>         export_extra=['DEM', 'localIncidenceAngle'], t_srs=4326)

    See Also
    --------
    :class:`pyroSAR.drivers.ID`,
    :class:`spatialist.vector.Vector`,
    :func:`spatialist.auxil.crsConvert()`
    """
    
    id = infile if isinstance(infile, pyroSAR.ID) else pyroSAR.identify(infile)
    
    if id.is_processed(outdir):
        print('scene {} already processed'.format(id.outname_base()))
        return
    # print(os.path.basename(id.scene))
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    ############################################
    # general setup
    
    if id.sensor in ['ASAR', 'ERS1', 'ERS2']:
        formatName = 'ENVISAT'
    elif id.sensor in ['S1A', 'S1B']:
        if id.product == 'SLC':
            raise RuntimeError('Sentinel-1 SLC data is not supported yet')
        formatName = 'SENTINEL-1'
    else:
        raise RuntimeError('sensor not supported (yet)')
    ######################
    # print('- assessing polarization selection')
    if isinstance(polarizations, str):
        if polarizations == 'all':
            polarizations = id.polarizations
        else:
            if polarizations in id.polarizations:
                polarizations = [polarizations]
            else:
                raise RuntimeError('polarization {} does not exists in the source product'.format(polarizations))
    elif isinstance(polarizations, list):
        polarizations = [x for x in polarizations if x in id.polarizations]
    else:
        raise RuntimeError('polarizations must be of type str or list')
    
    format = 'GeoTiff-BigTIFF' if len(polarizations) == 1 and export_extra is None else 'ENVI'
    # print(polarizations)
    # print(format)
    
    bands_int = ['Intensity_' + x for x in polarizations]
    bands_beta = ['Beta0_' + x for x in polarizations]
    bands_gamma = ['Gamma0_' + x for x in polarizations]
    ############################################
    ############################################
    # parse base workflow
    # print('- parsing base workflow')
    workflow = parse_recipe('base')
    ############################################
    # Read node configuration
    # print('-- configuring Read Node')
    read = workflow['Read']
    read.parameters['file'] = id.scene
    read.parameters['formatName'] = formatName
    ############################################
    # Remove-GRD-Border-Noise node configuration
    # print('-- configuring Remove-GRD-Border-Noise Node')
    if id.sensor in ['S1A', 'S1B'] and removeS1BorderNoise:
        bn = parse_node('Remove-GRD-Border-Noise')
        workflow.insert_node(bn, before='Read')
        bn.parameters['selectedPolarisations'] = polarizations
    ############################################
    # ThermalNoiseRemoval node configuration
    # print('-- configuring ThermalNoiseRemoval Node')
    if id.sensor in ['S1A', 'S1B'] and removeS1ThermalNoise:
        tn = parse_node('ThermalNoiseRemoval')
        workflow.insert_node(tn, before='Read')
        tn.parameters['selectedPolarisations'] = polarizations
    ############################################
    # orbit file application node configuration
    # print('-- configuring Apply-Orbit-File Node')
    orbit_lookup = {'ENVISAT': 'DELFT Precise (ENVISAT, ERS1&2) (Auto Download)',
                    'SENTINEL-1': 'Sentinel Precise (Auto Download)'}
    orbitType = orbit_lookup[formatName]
    if formatName == 'ENVISAT' and id.acquisition_mode == 'WSM':
        orbitType = 'DORIS Precise VOR (ENVISAT) (Auto Download)'
    
    orb = workflow['Apply-Orbit-File']
    orb.parameters['orbitType'] = orbitType
    ############################################
    # calibration node configuration
    # print('-- configuring Calibration Node')
    cal = workflow['Calibration']
    cal.parameters['selectedPolarisations'] = polarizations
    cal.parameters['sourceBands'] = bands_int
    ############################################
    # terrain flattening node configuration
    # print('-- configuring Terrain-Flattening Node')
    if terrainFlattening:
        tf = parse_node('Terrain-Flattening')
        workflow.insert_node(tf, before='Calibration')
        if id.sensor in ['ERS1', 'ERS2'] or (id.sensor == 'ASAR' and id.acquisition_mode != 'APP'):
            tf.parameters['sourceBands'] = 'Beta0'
        else:
            tf.parameters['sourceBands'] = bands_beta
        if externalDEMFile is None:
            tf.parameters['reGridMethod'] = True
        else:
            tf.parameters['reGridMethod'] = False
        pred_tc = 'Terrain-Flattening'
    else:
        cal.parameters['outputBetaBand'] = False
        cal.parameters['outputGammaBand'] = True
        pred_tc = 'Calibration'
    ############################################
    # configuration of node sequence for specific geocoding approaches
    # print('-- configuring geocoding approach Nodes')
    if geocoding_type == 'Range-Doppler':
        tc = parse_node('Terrain-Correction')
        workflow.insert_node(tc, before=pred_tc)
        tc.parameters['sourceBands'] = bands_gamma
    elif geocoding_type == 'SAR simulation cross correlation':
        sarsim = parse_node('SAR-Simulation')
        workflow.insert_node(sarsim, before=pred_tc)
        sarsim.parameters['sourceBands'] = bands_gamma
        
        workflow.insert_node(parse_node('Cross-Correlation'), before='SAR-Simulation')
        
        tc = parse_node('SARSim-Terrain-Correction')
        workflow.insert_node(tc, before='Cross-Correlation')
    else:
        raise RuntimeError('geocode_type not recognized')
    
    ############################################
    # Multilook node configuration
    
    try:
        image_geometry = id.meta['image_geometry']
        incidence = id.meta['incidence']
    except KeyError:
        raise RuntimeError('This function does not yet support sensor {}'.format(id.sensor))
    
    rlks, azlks = multilook_factors(sp_rg=id.spacing[0],
                                    sp_az=id.spacing[1],
                                    tr_rg=tr,
                                    tr_az=tr,
                                    geometry=image_geometry,
                                    incidence=incidence)
    
    if azlks > 1 or rlks > 1:
        workflow.insert_node(parse_node('Multilook'), before='Calibration')
        ml = workflow['Multilook']
        ml.parameters['nAzLooks'] = azlks
        ml.parameters['nRgLooks'] = rlks
        if cal.parameters['outputBetaBand']:
            ml.parameters['sourceBands'] = bands_beta
        elif cal.parameters['outputGammaBand']:
            ml.parameters['sourceBands'] = bands_gamma
    ############################################
    # specify spatial resolution and coordinate reference system of the output dataset
    # print('-- configuring CRS')
    tc.parameters['pixelSpacingInMeter'] = tr
    
    try:
        t_srs = crsConvert(t_srs, 'epsg')
    except TypeError:
        raise RuntimeError("format of parameter 't_srs' not recognized")
    
    # the EPSG code 4326 is not supported by SNAP and thus the WKT string has to be defined;
    # in all other cases defining EPSG:{code} will do
    if t_srs == 4326:
        t_srs = 'GEOGCS["WGS84(DD)",' \
                'DATUM["WGS84",' \
                'SPHEROID["WGS84", 6378137.0, 298.257223563]],' \
                'PRIMEM["Greenwich", 0.0],' \
                'UNIT["degree", 0.017453292519943295],' \
                'AXIS["Geodetic longitude", EAST],' \
                'AXIS["Geodetic latitude", NORTH]]'
    else:
        t_srs = 'EPSG:{}'.format(t_srs)
    
    tc.parameters['mapProjection'] = t_srs
    ############################################
    # (optionally) add node for conversion from linear to db scaling
    # print('-- configuring LinearToFromdB Node')
    if scaling not in ['dB', 'db', 'linear']:
        raise RuntimeError('scaling must be  a string of either "dB", "db" or "linear"')
    
    if scaling in ['dB', 'db']:
        lin2db = parse_node('lin2db')
        sourceNode = 'Terrain-Correction' if geocoding_type == 'Range-Doppler' else 'SARSim-Terrain-Correction'
        workflow.insert_node(lin2db, before=sourceNode)
        
        lin2db.parameters['sourceBands'] = bands_gamma
    
    ############################################
    # (optionally) add subset node and add bounding box coordinates of defined shapefile
    # print('-- configuring Subset Node')
    if shapefile:
        # print('--- read')
        shp = shapefile.clone() if isinstance(shapefile, Vector) else Vector(shapefile)
        # reproject the geometry to WGS 84 latlon
        # print('--- reproject')
        shp.reproject(4326)
        ext = shp.extent
        shp.close()
        # add an extra buffer of 0.01 degrees
        buffer = 0.01
        ext['xmin'] -= buffer
        ext['ymin'] -= buffer
        ext['xmax'] += buffer
        ext['ymax'] += buffer
        # print('--- create bbox')
        with bbox(ext, 4326) as bounds:
            # print('--- intersect')
            inter = intersect(id.bbox(), bounds)
            if not inter:
                raise RuntimeError('no bounding box intersection between shapefile and scene')
            # print('--- close intersect')
            inter.close()
            # print('--- get wkt')
            wkt = bounds.convert2wkt()[0]
        
        # print('--- parse XML node')
        subset = parse_node('Subset')
        # print('--- insert node')
        workflow.insert_node(subset, before='Read')
        subset.parameters['region'] = [0, 0, id.samples, id.lines]
        subset.parameters['geoRegion'] = wkt
    ############################################
    # (optionally) configure subset node for pixel offsets
    if offset and not shapefile:
        subset = parse_node('Subset')
        workflow.insert_node(subset, before='Read')
        
        # left, right, top and bottom offset in pixels
        l, r, t, b = offset
        
        subset_values = [l, t, id.samples - l - r, id.lines - t - b]
        subset.parameters['region'] = subset_values
        subset.parameters['geoRegion'] = ''
    ############################################
    # parametrize write node
    # print('-- configuring Write Node')
    # create a suffix for the output file to identify processing steps performed in the workflow
    suffix = workflow.suffix
    
    basename = os.path.join(outdir, id.outname_base(basename_extensions))
    extension = suffix if format == 'ENVI' else polarizations[0] + '_' + suffix
    outname = basename + '_' + extension
    
    write = workflow['Write']
    write.parameters['file'] = outname
    write.parameters['formatName'] = format
    ############################################
    ############################################
    if export_extra is not None:
        options = ['incidenceAngleFromEllipsoid',
                   'localIncidenceAngle',
                   'projectedLocalIncidenceAngle',
                   'DEM']
        write = parse_node('Write')
        workflow.insert_node(write, before=tc.id, resetSuccessorSource=False)
        write.parameters['file'] = outname
        write.parameters['formatName'] = format
        for item in export_extra:
            if item not in options:
                raise RuntimeError("ID '{}' not valid for argument 'export_extra'".format(item))
            key = 'save{}{}'.format(item[0].upper(), item[1:])
            tc.parameters[key] = True
    ############################################
    ############################################
    # select DEM type
    # print('-- configuring DEM')
    dempar = {'externalDEMFile': externalDEMFile,
              'externalDEMApplyEGM': externalDEMApplyEGM}
    if externalDEMFile is not None:
        if os.path.isfile(externalDEMFile):
            if externalDEMNoDataValue is None:
                with Raster(externalDEMFile) as dem:
                    dempar['externalDEMNoDataValue'] = dem.nodata
                if dempar['externalDEMNoDataValue'] is None:
                    raise RuntimeError('Cannot read NoData value from DEM file. '
                                       'Please specify externalDEMNoDataValue')
            else:
                dempar['externalDEMNoDataValue'] = externalDEMNoDataValue
            dempar['reGridMethod'] = False
        else:
            raise RuntimeError('specified externalDEMFile does not exist')
        dempar['demName'] = 'External DEM'
    else:
        # SRTM 1arcsec is only available between -58 and +60 degrees.
        # If the scene exceeds those latitudes SRTM 3arcsec is selected.
        if id.getCorners()['ymax'] > 60 or id.getCorners()['ymin'] < -58:
            dempar['demName'] = 'SRTM 3Sec'
        else:
            dempar['demName'] = 'SRTM 1Sec HGT'
        dempar['externalDEMFile'] = None
        dempar['externalDEMNoDataValue'] = 0
    
    for key, value in dempar.items():
        workflow.set_par(key, value)
    ############################################
    ############################################
    # configure the resampling methods
    
    options = ['NEAREST_NEIGHBOUR',
               'BILINEAR_INTERPOLATION',
               'CUBIC_CONVOLUTION',
               'BISINC_5_POINT_INTERPOLATION',
               'BISINC_11_POINT_INTERPOLATION',
               'BISINC_21_POINT_INTERPOLATION',
               'BICUBIC_INTERPOLATION']
    
    message = '{0} must be one of the following:\n- {1}'
    if demResamplingMethod not in options:
        raise ValueError(message.format('demResamplingMethod', '\n- '.join(options)))
    if imgResamplingMethod not in options:
        raise ValueError(message.format('imgResamplingMethod', '\n- '.join(options)))
    
    workflow.set_par('demResamplingMethod', demResamplingMethod)
    workflow.set_par('imgResamplingMethod', imgResamplingMethod)
    ############################################
    ############################################
    # write workflow to file and optionally execute it
    # print('- writing workflow to file')
    
    workflow.write(outname + '_proc')
    
    # execute the newly written workflow
    if not test:
        try:
            groups = groupbyWorkers(outname + '_proc.xml', groupsize)
            gpt(outname + '_proc.xml', groups=groups, cleanup=cleanup,
                gpt_exceptions=gpt_exceptions)
        except RuntimeError as e:
            if cleanup:
                if os.path.isdir(outname):
                    shutil.rmtree(outname)
                elif os.path.isfile(outname):
                    os.remove(outname)
                os.remove(outname + '_proc.xml')
            with open(outname + '_error.log', 'w') as log:
                log.write(str(e))
    
    if returnWF:
        return outname + '_proc.xml'
