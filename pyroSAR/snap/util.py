####################################################################
# Convenience functions for SAR image batch processing with ESA SNAP
# John Truckenbrodt, 2016-2018
####################################################################
import os

import pyroSAR
from .auxil import parse_recipe, parse_suffix, write_recipe, parse_node, insert_node, gpt

from spatialist import crsConvert, Vector, Raster, bbox, intersect


def geocode(infile, outdir, t_srs=4326, tr=20, polarizations='all', shapefile=None, scaling='dB',
            geocoding_type='Range-Doppler', removeS1BoderNoise=True, offset=None, externalDEMFile=None,
            externalDEMNoDataValue=None, externalDEMApplyEGM=True, basename_extensions=None, test=False):
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
    removeS1BoderNoise: bool, optional
        Enables removal of S1 GRD border noise (default).
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
    basename_extensions: list of str
        names of additional parameters to append to the basename, e.g. ['orbitNumber_rel']
    test: bool, optional
        If set to True the workflow xml file is only written and not executed. Default is False.

    Note
    ----
    If only one polarization is selected the results are directly written to GeoTiff.
    Otherwise the results are first written to a folder containing ENVI files and then transformed to GeoTiff files
    (one for each polarization).
    If GeoTiff would directly be selected as output format for multiple polarizations then a multilayer GeoTiff
    is written by SNAP which is considered an unfavorable format

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
    
    format = 'GeoTiff-BigTIFF' if len(polarizations) == 1 else 'ENVI'
    # print(polarizations)
    # print(format)
    
    bands_int = ','.join(['Intensity_' + x for x in polarizations])
    bands_beta = ','.join(['Beta0_' + x for x in polarizations])
    bands_gamma = ','.join(['Gamma0_' + x for x in polarizations])
    ############################################
    ############################################
    # parse base workflow
    # print('- parsing base workflow')
    workflow = parse_recipe('geocode')
    ############################################
    # Read node configuration
    # print('-- configuring Read Node')
    read = workflow.find('.//node[@id="Read"]')
    read.find('.//parameters/file').text = id.scene
    read.find('.//parameters/formatName').text = formatName
    ############################################
    # Remove-GRD-Border-Noise node configuration
    # print('-- configuring Remove-GRD-Border-Noise Node')
    if id.sensor in ['S1A', 'S1B'] and removeS1BoderNoise:
        insert_node(workflow, parse_node('Remove-GRD-Border-Noise'), before='Read')
        bn = workflow.find('.//node[@id="Remove-GRD-Border-Noise"]')
        bn.find('.//parameters/selectedPolarisations').text = ','.join(polarizations)
    ############################################
    # orbit file application node configuration
    # print('-- configuring Apply-Orbit-File Node')
    orbit_lookup = {'ENVISAT': 'DELFT Precise (ENVISAT, ERS1&2) (Auto Download)',
                    'SENTINEL-1': 'Sentinel Precise (Auto Download)'}
    orbitType = orbit_lookup[formatName]
    if formatName == 'ENVISAT' and id.acquisition_mode == 'WSM':
        orbitType = 'DORIS Precise VOR (ENVISAT) (Auto Download)'
    
    orb = workflow.find('.//node[@id="Apply-Orbit-File"]')
    orb.find('.//parameters/orbitType').text = orbitType
    ############################################
    # calibration node configuration
    # print('-- configuring Calibration Node')
    cal = workflow.find('.//node[@id="Calibration"]')
    
    cal.find('.//parameters/selectedPolarisations').text = ','.join(polarizations)
    cal.find('.//parameters/sourceBands').text = bands_int
    ############################################
    # terrain flattening node configuration
    # print('-- configuring Terrain-Flattening Node')
    tf = workflow.find('.//node[@id="Terrain-Flattening"]')
    if id.sensor in ['ERS1', 'ERS2'] or (id.sensor == 'ASAR' and id.acquisition_mode != 'APP'):
        tf.find('.//parameters/sourceBands').text = 'Beta0'
    else:
        tf.find('.//parameters/sourceBands').text = bands_beta
    ############################################
    # configuration of node sequence for specific geocoding approaches
    # print('-- configuring geocoding approach Nodes')
    if geocoding_type == 'Range-Doppler':
        insert_node(workflow, parse_node('Terrain-Correction'), before='Terrain-Flattening')
        tc = workflow.find('.//node[@id="Terrain-Correction"]')
        tc.find('.//parameters/sourceBands').text = bands_gamma
    elif geocoding_type == 'SAR simulation cross correlation':
        insert_node(workflow, parse_node('SAR-Simulation'), before='Terrain-Flattening')
        insert_node(workflow, parse_node('Cross-Correlation'), before='SAR-Simulation')
        insert_node(workflow, parse_node('SARSim-Terrain-Correction'), before='Cross-Correlation')
        tc = workflow.find('.//node[@id="SARSim-Terrain-Correction"]')
        
        sarsim = workflow.find('.//node[@id="SAR-Simulation"]')
        sarsim.find('.//parameters/sourceBands').text = bands_gamma
    else:
        raise RuntimeError('geocode_type not recognized')
    
    ############################################
    # specify spatial resolution and coordinate reference system of the output dataset
    # print('-- configuring CRS')
    tc.find('.//parameters/pixelSpacingInMeter').text = str(tr)
    
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
    
    tc.find('.//parameters/mapProjection').text = t_srs
    ############################################
    # (optionally) add node for conversion from linear to db scaling
    # print('-- configuring LinearToFromdB Node')
    if scaling not in ['dB', 'db', 'linear']:
        raise RuntimeError('scaling must be  a string of either "dB", "db" or "linear"')
    
    if scaling in ['dB', 'db']:
        lin2db = parse_node('lin2db')
        sourceNode = 'Terrain-Correction' if geocoding_type == 'Range-Doppler' else 'SARSim-Terrain-Correction'
        insert_node(workflow, lin2db, before=sourceNode)
        
        lin2db = workflow.find('.//node[@id="LinearToFromdB"]')
        lin2db.find('.//parameters/sourceBands').text = bands_gamma
    
    ############################################
    # (optionally) add subset node and add bounding box coordinates of defined shapefile
    # print('-- configuring Subset Node')
    if shapefile:
        # print('--- read')
        shp = shapefile if isinstance(shapefile, Vector) else Vector(shapefile)
        # reproject the geometry to WGS 84 latlon
        # print('--- reproject')
        shp.reproject(4326)
        ext = shp.extent
        # add an extra buffer of 0.01 degrees
        buffer = 0.01
        ext['xmin'] -= buffer
        ext['ymin'] -= buffer
        ext['xmax'] += buffer
        ext['ymax'] += buffer
        #print('--- create bbox')
        with bbox(ext, shp.srs) as bounds:
            # print('--- intersect')
            print(shapefile.srs)
            inter = intersect(id.bbox(), bounds)
            if not inter:
                raise RuntimeError('no bounding box intersection between shapefile and scene')
            # print('--- close intersect')
            inter.close()
            # print('--- get wkt')
            wkt = bounds.convert2wkt()[0]
        if isinstance(shapefile, str):
            shp.close()
        # print('--- parse XML node')
        subset = parse_node('Subset')
        # print('--- insert node')
        insert_node(workflow, subset, before='Read')
        
        subset = workflow.find('.//node[@id="Subset"]')
        subset.find('.//parameters/region').text = ','.join(map(str, [0, 0, id.samples, id.lines]))
        subset.find('.//parameters/geoRegion').text = wkt
    ############################################
    # (optionally) configure subset node for pixel offsets
    if offset and not shapefile:
        subset = parse_node('Subset')
        insert_node(workflow, subset, before='Read')
        
        # left, right, top and bottom offset in pixels
        l, r, t, b = offset
        
        subset = workflow.find('.//node[@id="Subset"]')
        subset_values = ','.join(map(str, [l, t, id.samples - l - r, id.lines - t - b]))
        subset.find('.//parameters/region').text = subset_values
        subset.find('.//parameters/geoRegion').text = ''
    ############################################
    # parametrize write node
    # print('-- configuring Write Node')
    # create a suffix for the output file to identify processing steps performed in the workflow
    suffix = parse_suffix(workflow)
    
    basename = os.path.join(outdir, id.outname_base(basename_extensions))
    extension = suffix if format == 'ENVI' else polarizations[0] + '_' + suffix
    outname = basename + '_' + extension
    
    write = workflow.find('.//node[@id="Write"]')
    write.find('.//parameters/file').text = outname
    write.find('.//parameters/formatName').text = format
    ############################################
    ############################################
    # select DEM type
    # print('-- configuring DEM')
    if externalDEMFile is not None:
        if os.path.isfile(externalDEMFile):
            if externalDEMNoDataValue is None:
                with Raster(externalDEMFile) as dem:
                    externalDEMNoDataValue = dem.nodata
                if externalDEMNoDataValue is None:
                    raise RuntimeError('Cannot read NoData value from DEM file. '
                                       'Please specify externalDEMNoDataValue')
        else:
            raise RuntimeError('specified externalDEMFile does not exist')
        demname = 'External DEM'
    else:
        # SRTM 1arcsec is only available between -58 and +60 degrees.
        # If the scene exceeds those latitudes SRTM 3arcsec is selected.
        if id.getCorners()['ymax'] > 60 or id.getCorners()['ymin'] < -58:
            demname = 'SRTM 3Sec'
        else:
            demname = 'SRTM 1Sec HGT'
        externalDEMFile = None
        externalDEMNoDataValue = 0
    
    for demName in workflow.findall('.//parameters/demName'):
        demName.text = demname
    for externalDEM in workflow.findall('.//parameters/externalDEMFile'):
        externalDEM.text = externalDEMFile
    for demNodata in workflow.findall('.//parameters/externalDEMNoDataValue'):
        demNodata.text = str(externalDEMNoDataValue)
    for egm in workflow.findall('.//parameters/externalDEMApplyEGM'):
        egm.text = str(externalDEMApplyEGM).lower()
    ############################################
    ############################################
    # write workflow to file and optionally execute it
    # print('- writing workflow to file')
    
    write_recipe(workflow, outname + '_proc')
    
    # execute the newly written workflow
    if not test:
        try:
            gpt(outname + '_proc.xml')
        except RuntimeError:
            os.remove(outname + '_proc.xml')
