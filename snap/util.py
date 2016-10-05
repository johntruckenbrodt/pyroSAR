import os
import re
import shutil
import subprocess as sp
from ancillary import finder
from auxil import suffix_lookup, parse_recipe, write_recipe
import pyroSAR


def geocode(infile, outdir, polarizations='all'):
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
    ############################################
    dB = workflow.find('.//node[@id="LinearToFromdB"]')
    dB.find('.//parameters/sourceBands').text = ','.join(['Gamma0_' + x for x in polarizations])
    ############################################
    write = workflow.find('.//node[@id="Write"]')
    write.find('.//parameters/file').text = outname
    write.find('.//parameters/formatName').text = format

    write_recipe(workflow, outname + '_proc')

    if format == 'GeoTiff-BigTIFF':
        sp.check_call(['gpt',
                       # '-Dsnap.dataio.reader.tileWidth=*',
                       # '-Dsnap.dataio.reader.tileHeight=1',
                       '-Dsnap.dataio.bigtiff.tiling.width=256',
                       '-Dsnap.dataio.bigtiff.tiling.height=256',
                       # '-Dsnap.dataio.bigtiff.compression.type=LZW',
                       # '-Dsnap.dataio.bigtiff.compression.quality=0.75',
                       outname + '_proc.xml'])
    else:
        sp.check_call(['gpt', outname + '_proc.xml'])

    if format == 'ENVI':
        for item in finder(outname, ['*.img']):
            pol = re.search('[HV]{2}', item).group()
            name_new = os.path.join(outdir, '{}_{}_{}'.format(id.outname_base(), pol, suffix))
            sp.check_call(['gdal_translate', '-q', '-of', 'GTiff',
                           '-co', 'INTERLEAVE=BAND', '-co', 'TILED=YES',
                           item, name_new + '.tif'])
        shutil.rmtree(outname)
