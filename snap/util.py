
import os
import shutil
import subprocess as sp
from auxil import suffix_lookup, parse_recipe, write_recipe
from pyroSAR import identify


def geocode(infile, outdir, polarizations='all'):

    id = identify(infile)

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
    elif format == 'GeoTiff-BigTIFF':
        outname = os.path.join(outdir, '{}_{}_{}'.format(id.outname_base(), polarizations[0], suffix))

    read = workflow.find('.//node[@id="Read"]')
    read.find('.//parameters/file').text = infile

    cal = workflow.find('.//node[@id="Calibration"]')
    cal.find('.//parameters/selectedPolarisations').text = ','.join(polarizations)

    tf = workflow.find('.//node[@id="Terrain-Flattening"]')
    tf.find('.//parameters/sourceBands').text = ','.join(['Beta0_'+x for x in polarizations])

    tc = workflow.find('.//node[@id="Terrain-Correction"]')
    tc.find('.//parameters/sourceBands').text = ','.join(['Gamma0_'+x for x in polarizations])

    write = workflow.find('.//node[@id="Write"]')
    write.find('.//parameters/file').text = outname
    write.find('.//parameters/formatName').text = format

    write_recipe(workflow, outname+'_proc')

    if format == 'GeoTiff-BigTIFF':
        sp.check_call(['gpt',
                       '-Dsnap.dataio.bigtiff.tiling.width=256',
                       '-Dsnap.dataio.bigtiff.tiling.height =256',
                       # '-Dsnap.dataio.bigtiff.compression.type=LZW',
                       # '-Dsnap.dataio.bigtiff.compression.quality=0.75',
                       outname + '_proc.xml'])
    else:
        sp.check_call(['gpt', outname+'_proc.xml'])

    if format == 'ENVI':
        for pol in polarizations:
            name_old = os.path.join(outname, 'Gamma0_{}'.format(pol))
            name_new = os.path.join(outdir, '{}_{}_{}'.format(id.outname_base(), pol, suffix))
            sp.check_call(['gdal_translate', '-q', '-of', 'GTiff',
                           '-co', 'INTERLEAVE=BAND', '-co', 'TILED=YES',
                           name_old+'.img', name_new+'.tif'])
        shutil.rmtree(outname)
