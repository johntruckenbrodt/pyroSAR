import os
import re
import shutil
import subprocess as sp
import spatial
from ancillary import finder
from spatial.vector import Vector
from .auxil import suffix_lookup, parse_recipe, write_recipe, parse_node, insert_node
import pyroSAR


def geocode(infile, outdir, t_srs=None, tr=20, polarizations='all', shapefile=None, scaling='dB'):
    orbit_lookup = {'ENVISAT': 'DELFT Precise (ENVISAT, ERS1&2) (Auto Download)',
                    'SENTINEL-1': 'Sentinel Precise (Auto Download)'}

    id = infile if isinstance(infile, pyroSAR.ID) else pyroSAR.identify(infile)

    if id.sensor in ['ASAR', 'ERS1', 'ERS2']:
        formatName = 'ENVISAT'
    elif id.sensor in ['S1A', 'S1B']:
        formatName = 'SENTINEL-1'
    else:
        raise RuntimeError('sensor not supported (yet)')

    workflow = parse_recipe('geocode')

    nodes = workflow.findall('node')

    suffix = '_'.join(filter(None, [suffix_lookup[x] for x in [y.attrib['id'] for y in nodes]]))

    if polarizations == 'all':
        polarizations = id.polarizations
    else:
        polarizations = [x for x in polarizations if x in id.polarizations]

    format = 'GeoTiff-BigTIFF' if len(polarizations) == 1 else 'ENVI'

    if format == 'ENVI':
        outname = os.path.join(outdir, id.outname_base() + '_' + suffix)
    else:
        outname = os.path.join(outdir, '{}_{}_{}'.format(id.outname_base(), polarizations[0], suffix))

    ############################################
    read = workflow.find('.//node[@id="Read"]')
    read.find('.//parameters/file').text = id.scene
    read.find('.//parameters/formatName').text = formatName
    ############################################
    orb = workflow.find('.//node[@id="Apply-Orbit-File"]')
    orb.find('.//parameters/orbitType').text = orbit_lookup[formatName]
    ############################################
    cal = workflow.find('.//node[@id="Calibration"]')
    if id.sensor in ['ERS1', 'ERS2'] or (id.sensor == 'ASAR' and id.acquisition_mode != 'APP'):
        cal.find('.//parameters/selectedPolarisations').text = 'Intensity'
    else:
        cal.find('.//parameters/selectedPolarisations').text = ','.join(polarizations)
    ############################################
    tf = workflow.find('.//node[@id="Terrain-Flattening"]')
    if id.sensor in ['ERS1', 'ERS2'] or (id.sensor == 'ASAR' and id.acquisition_mode != 'APP'):
        tf.find('.//parameters/sourceBands').text = 'Beta0'
    else:
        tf.find('.//parameters/sourceBands').text = ','.join(['Beta0_' + x for x in polarizations])
    ############################################
    tc = workflow.find('.//node[@id="Terrain-Correction"]')
    tc.find('.//parameters/sourceBands').text = ','.join(['Gamma0_' + x for x in polarizations])
    tc.find('.//parameters/pixelSpacingInMeter').text = str(tr)
    if t_srs:
        t_srs = spatial.crsConvert(t_srs, 'wkt')
        tc.find('.//parameters/mapProjection').text = t_srs
    else:
        tc.find('.//parameters/mapProjection').text = \
            'GEOGCS["WGS84(DD)",' \
            'DATUM["WGS84",' \
            'SPHEROID["WGS84", 6378137.0, 298.257223563]],' \
            'PRIMEM["Greenwich", 0.0],' \
            'UNIT["degree", 0.017453292519943295],' \
            'AXIS["Geodetic longitude", EAST],' \
            'AXIS["Geodetic latitude", NORTH]]'
    ############################################
    if scaling is 'dB':
        lin2db = parse_node('lin2db')
        insert_node(workflow, 'Terrain-Correction', lin2db)

        lin2db = workflow.find('.//node[@id="LinearToFromdB"]')
        lin2db.find('.//parameters/sourceBands').text = ','.join(['Gamma0_' + x for x in polarizations])

    # dB = workflow.find('.//node[@id="LinearToFromdB"]')
    # dB.find('.//parameters/sourceBands').text = ','.join(['Gamma0_' + x for x in polarizations])
    ############################################
    write = workflow.find('.//node[@id="Write"]')
    write.find('.//parameters/file').text = outname
    write.find('.//parameters/formatName').text = format
    ############################################
    if shapefile:
        shp = shapefile if isinstance(shapefile, Vector) else Vector(shapefile)
        bounds = spatial.bbox(shp.extent, shp.wkt)
        bounds.reproject(id.projection)
        intersect = spatial.intersect(id.bbox(), bounds)
        if not intersect:
            raise RuntimeError('no bounding box intersection between shapefile and scene')
        wkt = bounds.convert2wkt()[0]
        subset = parse_node('subset')

        insert_node(workflow, 'Read', subset)

        subset = workflow.find('.//node[@id="Subset"]')
        subset.find('.//parameters/region').text = ','.join(map(str, [0, 0, id.samples, id.lines]))
        subset.find('.//parameters/geoRegion').text = wkt
        ############################################

    write_recipe(workflow, outname + '_proc')

    if format == 'GeoTiff-BigTIFF':
        cmd = ['gpt',
               # '-Dsnap.dataio.reader.tileWidth=*',
               # '-Dsnap.dataio.reader.tileHeight=1',
               '-Dsnap.dataio.bigtiff.tiling.width=256',
               '-Dsnap.dataio.bigtiff.tiling.height=256',
               # '-Dsnap.dataio.bigtiff.compression.type=LZW',
               # '-Dsnap.dataio.bigtiff.compression.quality=0.75',
               outname + '_proc.xml']
    else:
        cmd = ['gpt', outname + '_proc.xml']

    proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = proc.communicate()
    if proc.returncode != 0:
        print err
        print 'failed: ', os.path.basename(id.scene)
        os.remove(outname + '_proc.xml')
        if os.path.isfile(outname + '.tif'):
            os.remove(outname + '.tif')
        elif os.path.isdir(outname):
            shutil.rmtree(outname)
        return

    if format == 'ENVI':
        for item in finder(outname, ['*.img']):
            pol = re.search('[HV]{2}', item).group()
            name_new = os.path.join(outdir, '{}_{}_{}'.format(id.outname_base(), pol, suffix))
            sp.check_call(['gdal_translate', '-q', '-of', 'GTiff',
                           '-co', 'INTERLEAVE=BAND', '-co', 'TILED=YES',
                           item, name_new + '.tif'])
        shutil.rmtree(outname)
