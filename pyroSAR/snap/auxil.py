###############################################################################
# pyroSAR SNAP API tools

# Copyright (c) 2017-2025, the pyroSAR Developers.

# This file is part of the pyroSAR Project. It is subject to the
# license terms in the LICENSE.txt file found in the top-level
# directory of this distribution and at
# https://github.com/johntruckenbrodt/pyroSAR/blob/master/LICENSE.txt.
# No part of the pyroSAR project, including this file, may be
# copied, modified, propagated, or distributed except according
# to the terms contained in the LICENSE.txt file.
###############################################################################
import os
import re
import copy
import shutil
import traceback
import subprocess as sp
from xml.dom import minidom
import xml.etree.ElementTree as ET

from pyroSAR import identify
from pyroSAR.examine import ExamineSnap
from pyroSAR.ancillary import windows_fileprefix, multilook_factors, Lock
from pyroSAR.auxdata import get_egm_lookup

from spatialist import Vector, Raster, vectorize, rasterize, boundary, intersect, bbox
from spatialist.auxil import gdal_translate, crsConvert
from spatialist.ancillary import finder, run

from osgeo import gdal
from osgeo.gdalconst import GA_Update

import logging

log = logging.getLogger(__name__)


def parse_recipe(name):
    """
    parse a SNAP recipe
    
    Parameters
    ----------
    name: str
        the name of the recipe; current options:
         * `blank`: a workflow without any nodes
         * `geocode`: a basic workflow containing `Read`, `Apply-Orbit-File`,
           `Calibration`, `Terrain-Flattening` and `Write` nodes

    Returns
    -------
    Workflow
        the parsed recipe
    
    Examples
    --------
    >>> from pyroSAR.snap.auxil import parse_recipe
    >>> workflow = parse_recipe('base')
    """
    name = name if name.endswith('.xml') else name + '.xml'
    absname = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'recipes', name)
    return Workflow(absname)


def parse_node(name, use_existing=True):
    """
    parse an XML node recipe. The XML representation and parameter default values are read from the docstring of an
    individual node by calling `gpt <node> -h`. The result is then written to an XML text file under
    `$HOME/.pyroSAR/snap/nodes` which is subsequently read for parsing instead of again calling `gpt`.
    
    Parameters
    ----------
    name: str
        the name of the processing node, e.g. Terrain-Correction
    use_existing: bool
        use an existing XML text file or force reparsing the gpt docstring and overwriting the XML file?

    Returns
    -------
    Node
        the parsed node
    
    Examples
    --------
    >>> tnr = parse_node('ThermalNoiseRemoval')
    >>> print(tnr.parameters)
    {'selectedPolarisations': None, 'removeThermalNoise': 'true', 'reIntroduceThermalNoise': 'false'}
    """
    snap = ExamineSnap()
    version = snap.get_version('microwavetbx')['version']
    name = name if name.endswith('.xml') else name + '.xml'
    operator = os.path.splitext(name)[0]
    nodepath = os.path.join(os.path.expanduser('~'), '.pyrosar', 'snap', 'nodes')
    abspath = os.path.join(nodepath, version)
    os.makedirs(abspath, exist_ok=True)
    absname = os.path.join(abspath, name)
    
    # remove all old XML files that were not stored in a version subdirectory
    deprecated = finder(nodepath, ['*.xml'], recursive=False)
    for item in deprecated:
        os.remove(item)
    
    with Lock(absname):
        if not os.path.isfile(absname) or not use_existing:
            gpt = snap.gpt
            
            cmd = [gpt, operator, '-h']
            
            out, err = run(cmd=cmd, void=False)
            
            if re.search('Unknown operator', out + err):
                raise RuntimeError("unknown operator '{}'".format(operator))
            
            graph = re.search('<graph id.*', out, flags=re.DOTALL).group()
            # remove placeholder values like ${value}
            graph = re.sub(r'>\${.*', '/>', graph)
            # remove <.../> placeholders
            graph = re.sub(r'<\.\.\./>.*', '', graph)
            if operator == 'BandMaths':
                graph = graph.replace('sourceProducts', 'sourceProduct')
            tree = ET.fromstring(graph)
            for elt in tree.iter():
                if elt.text in ['string', 'double', 'integer', 'float']:
                    elt.text = None
            node = tree.find('node')
            node.attrib['id'] = operator
            # add a second source product entry for multi-source nodes
            # multi-source nodes are those with an entry 'sourceProducts'
            # instead of 'sourceProduct'
            # exceptions are registered in this list:
            multisource = ['Back-Geocoding']
            if operator != 'Read' and operator != 'ProductSet-Reader':
                source = node.find('.//sources')
                child = source[0]
                if child.tag == 'sourceProducts' or operator in multisource:
                    child2 = ET.SubElement(source,
                                           'sourceProduct.1',
                                           {'refid': 'Read (2)'})
                child.tag = 'sourceProduct'
                child.attrib['refid'] = 'Read'
                child.text = None
            
            # cleanup the BandMaths node
            if operator == 'BandMaths':
                tband = tree.find('.//targetBand')
                allowed = ['name', 'type', 'expression',
                           'description', 'unit', 'noDataValue']
                invalid = [x.tag for x in tband if x.tag not in allowed]
                for tag in invalid:
                    el = tband.find(f'.//{tag}')
                    tband.remove(el)
                for item in ['targetBands', 'variables']:
                    elem = tree.find(f'.//{item}')
                    pl = elem.find('.//_.002e..')
                    elem.remove(pl)
            
            # add a class parameter and create the Node object
            value = 'com.bc.ceres.binding.dom.XppDomElement'
            tree.find('.//parameters').set('class', value)
            node = Node(node)
            
            # read the default values from the parameter documentation
            parameters = node.parameters.keys()
            out += '-P'
            for parameter in parameters:
                p1 = r'-P{}.*?-P'.format(parameter)
                p2 = r"Default\ value\ is '([a-zA-Z0-9 ._\(\)]+)'"
                r1 = re.search(p1, out, re.S)
                if r1:
                    sub = r1.group()
                    r2 = re.search(p2, sub)
                    if r2:
                        value = r2.groups()[0]
                        node.parameters[parameter] = value
                        continue
                node.parameters[parameter] = None
            
            # fill in some additional defaults
            if operator == 'BandMerge':
                node.parameters['geographicError'] = '1.0E-5'
            
            with open(absname, 'w') as xml:
                xml.write(str(node))
            return node
        else:
            with open(absname, 'r') as workflow:
                element = ET.fromstring(workflow.read())
            return Node(element)


def execute(xmlfile, cleanup=True, gpt_exceptions=None, gpt_args=None):
    """
    execute SNAP workflows via the Graph Processing Tool GPT.
    This function merely calls gpt with some additional command
    line arguments and raises a RuntimeError on fail. This
    function is used internally by function :func:`gpt`.
    
    Parameters
    ----------
    xmlfile: str
        the name of the workflow XML file
    cleanup: bool
        should all files written to the temporary directory during function execution be deleted after processing?
    gpt_exceptions: dict
        a dictionary to override the configured GPT executable for certain operators;
        each (sub-)workflow containing this operator will be executed with the define executable;
        
         - e.g. ``{'Terrain-Flattening': '/home/user/snap/bin/gpt'}``
    gpt_args: list or None
        a list of additional arguments to be passed to the GPT call
        
        - e.g. ``['-x', '-c', '2048M']`` for increased tile cache size and intermediate clearing
    
    Returns
    -------
    
    Raises
    ------
    RuntimeError
    """
    # read the file and extract some information
    workflow = Workflow(xmlfile)
    write = workflow['Write']
    outname = write.parameters['file']
    workers = [x.id for x in workflow if x.operator not in ['Read', 'Write']]
    message = ' -> '.join(workers)
    gpt_exec = None
    if gpt_exceptions is not None:
        for item, exec in gpt_exceptions.items():
            if item in workers:
                gpt_exec = exec
                message += ' (using {})'.format(exec)
                break
    log.info(message)
    # try to find the GPT executable
    if gpt_exec is None:
        try:
            gpt_exec = ExamineSnap().gpt
        except AttributeError:
            raise RuntimeError('could not find SNAP GPT executable')
    # create the list of arguments to be passed to the subprocess module calling GPT
    cmd = [gpt_exec, '-e']
    if isinstance(gpt_args, list):
        cmd.extend(gpt_args)
    if format == 'GeoTiff-BigTIFF':
        cmd.extend([
            # '-Dsnap.dataio.reader.tileWidth=*',
            # '-Dsnap.dataio.reader.tileHeight=1',
            '-Dsnap.dataio.bigtiff.tiling.width=256',
            '-Dsnap.dataio.bigtiff.tiling.height=256',
            # '-Dsnap.dataio.bigtiff.compression.type=LZW',
            # '-Dsnap.dataio.bigtiff.compression.quality=0.75'
        ])
    cmd.append(xmlfile)
    # execute the workflow
    proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = proc.communicate()
    out = out.decode('utf-8') if isinstance(out, bytes) else out
    err = err.decode('utf-8') if isinstance(err, bytes) else err
    
    # check for a message indicating an unknown parameter,
    # which can easily be removed from the workflow
    pattern = r"Error: \[NodeId: (?P<id>[a-zA-Z0-9-_]*)\] " \
              r"Operator \'[a-zA-Z0-9-_]*\': " \
              r"Unknown element \'(?P<par>[a-zA-Z]*)\'"
    match = re.search(pattern, err)
    
    if proc.returncode == 0:
        pattern = r'(?P<level>WARNING: )([a-zA-Z.]*: )(?P<message>No intersection.*)'
        match = re.search(pattern, err)
        if match is not None:
            raise RuntimeError(re.search(pattern, err).group('message'))
        return
    
    # delete unknown parameters and run the modified workflow
    elif proc.returncode == 1 and match is not None:
        replace = match.groupdict()
        with Workflow(xmlfile) as flow:
            log.info('  removing parameter {id}:{par} and executing modified workflow'.format(**replace))
            node = flow[replace['id']]
            del node.parameters[replace['par']]
            flow.write(xmlfile)
        execute(xmlfile, cleanup=cleanup, gpt_exceptions=gpt_exceptions,
                gpt_args=gpt_args)
    
    # append additional information to the error message and raise an error
    else:
        if proc.returncode == -9:
            submessage = '[{}] the process was killed by SNAP (process return code -9). ' \
                         'One possible cause is a lack of memory.'.format(os.path.basename(xmlfile))
        else:
            submessage = '{}{}\n[{}] failed with return code {}'
        if cleanup:
            if os.path.isfile(outname + '.tif'):
                os.remove(outname + '.tif')
            elif os.path.isdir(outname):
                shutil.rmtree(outname, onerror=windows_fileprefix)
            elif outname.endswith('.dim') and os.path.isfile(outname):
                os.remove(outname)
                datadir = outname.replace('.dim', '.data')
                if os.path.isdir(datadir):
                    shutil.rmtree(datadir,
                                  onerror=windows_fileprefix)
        raise RuntimeError(submessage.format(out, err, os.path.basename(xmlfile), proc.returncode))


def gpt(xmlfile, tmpdir, groups=None, cleanup=True,
        gpt_exceptions=None, gpt_args=None,
        removeS1BorderNoiseMethod='pyroSAR'):
    """
    Wrapper for ESA SNAP's Graph Processing Tool GPT.
    Input is a readily formatted workflow XML file as for example
    created by function :func:`~pyroSAR.snap.util.geocode`.
    Additional to calling GPT, this function will
    
    - (if processing Sentinel-1 GRD data with IPF version <2.9 and ``removeS1BorderNoiseMethod='pyroSAR'``)
      unpack the scene and perform the custom removal (:func:`pyroSAR.S1.removeGRDBorderNoise`).
    - if `groups` is not None:
    
      * split the workflow into sub-workflows (:func:`pyroSAR.snap.auxil.split`)
      * execute the sub-workflows (:func:`pyroSAR.snap.auxil.execute`)
    
    Note
    ----
    Depending on the parametrization this function might create two subdirectories in `tmpdir`,
    bnr for S1 GRD border noise removal and sub for sub-workflows and their intermediate outputs.
    Both are deleted if ``cleanup=True``. If `tmpdir` is empty afterward it is also deleted.
    
    Parameters
    ----------
    xmlfile: str
        the name of the workflow XML file
    tmpdir: str
        a temporary directory for storing intermediate files
    groups: list[list[str]] or None
        a list of lists each containing IDs for individual nodes. If not None, the workflow is split into
        sub-workflows executing the nodes in the respective group. These workflows and their output products
        are stored into the subdirectory sub of `tmpdir`.
    cleanup: bool
        should all temporary files be deleted after processing? First, the subdirectories bnr and sub of `tmpdir`
        are deleted. If `tmpdir` is empty afterward it is also deleted.
    gpt_exceptions: dict or None
        a dictionary to override the configured GPT executable for certain operators;
        each (sub-)workflow containing this operator will be executed with the define executable;
        
         - e.g. ``{'Terrain-Flattening': '/home/user/snap/bin/gpt'}``
    
    gpt_args: list[str] or None
        a list of additional arguments to be passed to the gpt call
        
         - e.g. ``['-x', '-c', '2048M']`` for increased tile cache size and intermediate clearing
    
    removeS1BorderNoiseMethod: str
        the border noise removal method to be applied, See :func:`pyroSAR.S1.removeGRDBorderNoise` for details;
        one of the following:
        
         - 'ESA': the pure implementation as described by ESA
         - 'pyroSAR': the ESA method plus the custom pyroSAR refinement. This is only applied if the IPF version is
           < 2.9 where additional noise removal was necessary. The output of the additional noise removal is stored
           in the subdirectory bnr of `tmpdir`.
    
    Returns
    -------
    
    Raises
    ------
    
    """
    workflow = Workflow(xmlfile)
    
    if 'ProductSet-Reader' in workflow.operators:
        read = workflow['ProductSet-Reader']
        scene = identify(read.parameters['fileList'].split(',')[0])
    else:
        read = workflow['Read']
        scene = identify(read.parameters['file'])
    
    tmpdir_bnr = os.path.join(tmpdir, 'bnr')
    tmpdir_sub = os.path.join(tmpdir, 'sub')
    
    if 'Remove-GRD-Border-Noise' in workflow.ids \
            and removeS1BorderNoiseMethod == 'pyroSAR' \
            and scene.meta['IPF_version'] < 2.9:
        if 'SliceAssembly' in workflow.operators:
            raise RuntimeError("pyroSAR's custom border noise removal is not yet implemented for multiple scene inputs")
        os.makedirs(tmpdir_bnr, exist_ok=True)
        xmlfile = os.path.join(tmpdir_bnr,
                               os.path.basename(xmlfile.replace('_bnr', '')))
        # border noise removal is done outside of SNAP and the node is thus removed from the workflow
        del workflow['Remove-GRD-Border-Noise']
        # remove the node name from the groups
        i = 0
        while i < len(groups):
            if 'Remove-GRD-Border-Noise' in groups[i]:
                del groups[i][groups[i].index('Remove-GRD-Border-Noise')]
            if len(groups[i]) == 0:
                del groups[i]
            elif len(groups[i]) == 1 and groups[i][0] == 'Read':
                # move Read into the next group if it is the only operator
                del groups[i]
                groups[i].insert(0, 'Read')
            else:
                i += 1
        # unpack the scene if necessary and perform the custom border noise removal
        log.info('unpacking scene')
        if scene.compression is not None:
            scene.unpack(tmpdir_bnr)
        log.info('removing border noise..')
        scene.removeGRDBorderNoise(method=removeS1BorderNoiseMethod)
        # change the name of the input file to that of the unpacked archive
        read.parameters['file'] = scene.scene
        # write a new workflow file
        workflow.write(xmlfile)
    
    log.info('executing node sequence{}..'.format('s' if groups is not None else ''))
    try:
        if groups is not None:
            subs = split(xmlfile=xmlfile, groups=groups, outdir=tmpdir_sub)
            for sub in subs:
                execute(sub, cleanup=cleanup, gpt_exceptions=gpt_exceptions, gpt_args=gpt_args)
        else:
            execute(xmlfile, cleanup=cleanup, gpt_exceptions=gpt_exceptions, gpt_args=gpt_args)
    except Exception:
        tb = traceback.format_exc()
        log.info(tb)
        log.info('failed: {}'.format(xmlfile))
        raise
    finally:
        if cleanup:
            for tmp in [tmpdir_bnr, tmpdir_sub]:
                if os.path.isdir(tmp):
                    shutil.rmtree(tmp, onerror=windows_fileprefix)
            if os.path.isdir(tmpdir) and not os.listdir(tmpdir):
                shutil.rmtree(tmpdir, onerror=windows_fileprefix)


def writer(xmlfile, outdir, basename_extensions=None,
           clean_edges=False, clean_edges_npixels=1):
    """
    SNAP product writing utility
    
    Parameters
    ----------
    xmlfile: str
        the name of the workflow XML file.
    outdir: str
        the directory into which to write the final files.
    basename_extensions: list of str or None
        names of additional parameters to append to the basename, e.g. ``['orbitNumber_rel']``.
    clean_edges: bool
        erode noisy image edges? See :func:`pyroSAR.snap.auxil.erode_edges`.
        Does not apply to layover-shadow mask.
    clean_edges_npixels: int
        the number of pixels to erode.

    Returns
    -------

    """
    workflow = Workflow(xmlfile)
    writers = workflow['operator=Write']
    files = list(set([x.parameters['file'] for x in writers]))
    if len(files) > 1:
        raise RuntimeError('Multiple output files are not yet supported.')
    else:
        src = files[0]
    src_format = writers[0].parameters['formatName']
    suffix = workflow.suffix()
    rtc = 'Terrain-Flattening' in workflow.operators
    dem_name = workflow.tree.find('.//demName')
    dem_nodata = None
    if dem_name is not None:
        dem_name = dem_name.text
        if dem_name == 'External DEM':
            dem_nodata = float(workflow.tree.find('.//externalDEMNoDataValue').text)
        else:
            dem_nodata_lookup = {'SRTM 1Sec HGT': -32768}
            if dem_name in dem_nodata_lookup.keys():
                dem_nodata = dem_nodata_lookup[dem_name]
    
    src_base = os.path.splitext(os.path.basename(src))[0]
    outname_base = os.path.join(outdir, src_base)
    
    if src_format in ['ENVI', 'BEAM-DIMAP']:
        message = '{}converting to GeoTIFF'
        log.info(message.format('cleaning image edges and ' if clean_edges else ''))
        translateoptions = {'options': ['-q', '-co', 'INTERLEAVE=BAND', '-co', 'TILED=YES'],
                            'format': 'GTiff'}
        if clean_edges:
            erode_edges(src=src, only_boundary=True, pixels=clean_edges_npixels)
        
        if src_format == 'BEAM-DIMAP':
            src = src.replace('.dim', '.data')
        for item in finder(src, ['*.img'], recursive=False):
            pattern = '(?P<refarea>(?:Sig|Gam)ma0)_(?P<pol>[HV]{2})'
            basename = os.path.basename(item)
            match = re.search(pattern, basename)
            if match:
                refarea, pol = match.groups()
                correction = 'elp'
                if rtc:
                    if refarea == 'Gamma0':
                        correction = 'rtc'
                    elif refarea == 'Sigma0':
                        tf = workflow['Terrain-Flattening']
                        if tf.parameters['outputSigma0']:
                            correction = 'rtc'
                suffix_new = '{0}-{1}'.format(refarea.lower(), correction)
                if 'dB' in suffix:
                    suffix_new += '_db'
                name_new = outname_base.replace(suffix, '{0}_{1}.tif'.format(pol, suffix_new))
            else:
                base = os.path.splitext(basename)[0] \
                    .replace('elevation', 'DEM')
                if re.search('scatteringArea', base):
                    base = re.sub('scatteringArea_[HV]{2}', 'scatteringArea', base)
                if re.search('gammaSigmaRatio', base):
                    base = re.sub('gammaSigmaRatio_[HV]{2}', 'gammaSigmaRatio', base)
                if re.search('NE[BGS]Z', base):
                    base = re.sub('(NE[BGS]Z)_([HV]{2})', r'\g<2>_\g<1>', base)
                if re.search('layover_shadow_mask', base):
                    base = re.sub('layover_shadow_mask_[HV]{2}', 'layoverShadowMask', base)
                name_new = outname_base.replace(suffix, '{0}.tif'.format(base))
            if re.search('elevation', basename):
                nodata = dem_nodata
            elif re.search('layoverShadowMask|layover_shadow_mask', basename):
                nodata = 255
            else:
                nodata = 0
            translateoptions['noData'] = nodata
            gdal_translate(src=item, dst=name_new, **translateoptions)
    else:
        raise RuntimeError('The output file format must be ENVI or BEAM-DIMAP.')
    ###########################################################################
    # write the Sentinel-1 manifest.safe file as addition to the actual product
    readers = workflow['operator=Read']
    for reader in readers:
        infile = reader.parameters['file']
        try:
            id = identify(infile)
            if id.sensor in ['S1A', 'S1B', 'S1C', 'S1D']:
                manifest = id.getFileObj(id.findfiles('manifest.safe')[0])
                basename = id.outname_base(basename_extensions)
                basename = '{0}_manifest.safe'.format(basename)
                outname_manifest = os.path.join(outdir, basename)
                with open(outname_manifest, 'wb') as out:
                    out.write(manifest.read())
        except RuntimeError:
            continue


def is_consistent(workflow):
    """
    check whether all nodes take either no source node or one that is in the list
    
    Parameters
    ----------
    workflow: Workflow
        the workflow to be analyzed
    Returns
    -------
    bool
        is the list of nodes consistent?
    """
    ids = workflow.ids
    check = []
    for node in workflow:
        source = node.source
        if source is None or source in ids or all([x in ids for x in source]):
            check.append(True)
        else:
            check.append(False)
    for node in workflow:
        successors = workflow.successors(node.id, recursive=True)
        operators = [workflow[x].operator for x in successors]
        if node.operator == 'Write' or 'Write' in operators:
            check.append(True)
        else:
            log.debug('node {} does not have a Write successor'.format(node.id))
            check.append(False)
    return all(check)


def split(xmlfile, groups, outdir=None):
    """
    split a workflow file into groups and write them to separate workflows including source and write target linking.
    The new workflows are written to a sub-directory `temp` of the target directory defined in the input's `Write` node.
    Each new workflow is parameterized with a `Read` and `Write` node if they don't already exist. Temporary outputs are
    written to `BEAM-DIMAP` files named after the workflow suffix sequence.
    
    Parameters
    ----------
    xmlfile: str
        the workflow to be split
    groups: list
        a list of lists each containing IDs for individual nodes
    outdir: str or None
        the directory into which to write the XML workflows and the intermediate files created by them.
        If None, the name will be created from the file name of the node with ID 'Write',
        which is treated as a directory, and a subdirectory 'tmp'.

    Returns
    -------
    list of str
        the names of the newly written temporary workflows
    
    Raises
    ------
    RuntimeError
    """
    workflow = Workflow(xmlfile)
    write = workflow['Write']
    if outdir is None:
        out = write.parameters['file']
        outdir = os.path.join(out, 'tmp')
    os.makedirs(outdir, exist_ok=True)
    
    # the temporary XML files
    outlist = []
    # the names and format of temporary products
    prod_tmp = {}
    prod_tmp_format = {}
    for position, group in enumerate(groups):
        node_lookup = {}
        log.debug('creating new workflow for group {}'.format(group))
        new = parse_recipe('blank')
        nodes = [workflow[x] for x in group]
        for node in nodes:
            id_old = node.id
            sources = node.source
            if sources is None:
                sources = []
                resetSuccessorSource = False
            elif isinstance(sources, list):
                resetSuccessorSource = False
            else:
                resetSuccessorSource = True
                sources = [sources]
            reset = []
            for source in sources:
                if source not in group:
                    read = new.insert_node(parse_node('Read'), void=False,
                                           resetSuccessorSource=resetSuccessorSource)
                    reset.append(read.id)
                    read.parameters['file'] = prod_tmp[source]
                    read.parameters['formatName'] = prod_tmp_format[source]
                    node_lookup[read.id] = source
                else:
                    reset.append(source)
            if isinstance(sources, list):
                sources_new_pos = [list(node_lookup.values()).index(x) for x in sources]
                sources_new = [list(node_lookup.keys())[x] for x in sources_new_pos]
                newnode = new.insert_node(node.copy(), before=sources_new, void=False,
                                          resetSuccessorSource=False)
            else:
                newnode = new.insert_node(node.copy(), void=False,
                                          resetSuccessorSource=False)
            node_lookup[newnode.id] = id_old
            
            if not resetSuccessorSource:
                newnode.source = reset
        
        # if possible, read the name of the SAR product for parsing names of temporary files
        # this was found necessary for SliceAssembly, which expects the names in a specific format
        products = [x.parameters['file'] for x in new['operator=Read']]
        try:
            id = identify(products[0])
            filename = os.path.basename(id.scene)
        except (RuntimeError, OSError):
            filename = os.path.basename(products[0])
        basename = os.path.splitext(filename)[0]
        basename = re.sub(r'_tmp[0-9]+', '', basename)
        
        # add a Write node to all dangling nodes
        counter = 0
        for node in new:
            dependants = [x for x in workflow.successors(node.id) if not x.startswith('Write') and not x in group]
            if node.operator != 'Read' and len(dependants) > 0:
                write = parse_node('Write')
                new.insert_node(write, before=node.id, resetSuccessorSource=False)
                id = str(position) if counter == 0 else '{}-{}'.format(position, counter)
                tmp_out = os.path.join(outdir, '{}_tmp{}.dim'.format(basename, id))
                prod_tmp[node_lookup[node.id]] = tmp_out
                prod_tmp_format[node_lookup[node.id]] = 'BEAM-DIMAP'
                write.parameters['file'] = tmp_out
                write.parameters['formatName'] = 'BEAM-DIMAP'
                counter += 1
        if not is_consistent(new):
            message = 'inconsistent group:\n {}'.format(' -> '.join(group))
            raise RuntimeError(message)
        outname = os.path.join(outdir, '{}_tmp{}.xml'.format(basename, position))
        new.write(outname)
        outlist.append(outname)
    return outlist


def groupbyWorkers(xmlfile, n=2):
    """
    split a SNAP workflow into groups containing a maximum defined number of operators.
    
    Parameters
    ----------
    xmlfile: str
        the SNAP xml workflow
    n: int
        the maximum number of worker nodes in each group; Read, Write and BandSelect are excluded.

    Returns
    -------
    list[list[str]]
        a list of lists each containing the IDs of all nodes belonging to the groups including Read and Write nodes;
        this list can e.g. be passed to function :func:`split` to split the workflow into new sub-workflow files based
        on the newly created groups or directly to function :func:`gpt`, which will call :func:`split` internally.
    """
    workflow = Workflow(xmlfile)
    workers_id = [x.id for x in workflow if x.operator not in ['Read', 'Write', 'BandSelect']]
    readers_id = [x.id for x in workflow['operator=Read']]
    writers_id = [x.id for x in workflow['operator=Write']]
    selects_id = [x.id for x in workflow['operator=BandSelect']]
    workers_groups = [workers_id[i:i + n] for i in range(0, len(workers_id), n)]
    
    # some nodes must be executed together with a preceding node. They are moved to the previous group.
    def move_group(operator):
        i = 0
        while i < len(workers_groups):
            if workers_groups[i][0].startswith(operator):
                # get the group ID of the source node
                source = workflow[workers_groups[i][0]].source
                source_group_id = [source in x for x in workers_groups].index(True)
                # move the node to the source group
                workers_groups[source_group_id].append(workers_groups[i][0])
                del workers_groups[i][0]
            # delete the group if it is empty
            if len(workers_groups[i]) == 0:
                del workers_groups[i]
            else:
                i += 1
    
    for operator in ['ThermalNoiseRemoval', 'Warp']:
        move_group(operator)
    
    # append the BandSelect nodes to the group of their source nodes
    for item in selects_id:
        source = workflow[item].source
        for group in workers_groups:
            if source in group:
                group.insert(group.index(source) + 1, item)
    nodes_groups = []
    for group in workers_groups:
        newgroup = []
        for worker in group:
            newgroup.append(worker)
            source = workflow[worker].source
            if not isinstance(source, list):
                source = [source]
            for item in source:
                if item in readers_id:
                    # append all Read nodes that are the worker's direct sources
                    newgroup.insert(newgroup.index(worker), item)
            for writer in writers_id:
                if workflow[writer].source == worker:
                    # append all Write nodes that directly have the worker as source
                    newgroup.append(writer)
        nodes_groups.append(newgroup)
    return nodes_groups


class Workflow(object):
    """
    Class for convenient handling of SNAP XML workflows
    
    Parameters
    ----------
    xmlfile: str
        the workflow XML file
    """
    
    def __init__(self, xmlfile):
        with open(xmlfile, 'r') as infile:
            self.tree = ET.fromstring(infile.read())
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        pass
    
    def __getitem__(self, item):
        pattern = '(?P<key>[a-zA-Z-_]*)=(?P<value>[a-zA-Z-_]*)'
        if isinstance(item, int):
            return self.nodes()[item]
        elif isinstance(item, str):
            if re.search(pattern, item):
                key, value = re.search(pattern, item).groups()
                return [x for x in self if getattr(x, key) == value]
            else:
                try:
                    return Node(self.tree.find('.//node[@id="{}"]'.format(item)))
                except TypeError:
                    raise KeyError('unknown key: {}'.format(item))
        else:
            raise TypeError('item must be of type int or str')
    
    def __len__(self):
        return len(self.tree.findall('node'))
    
    def __delitem__(self, key):
        if not isinstance(key, str):
            raise TypeError('key must be of type str')
        element = self.tree.find('.//node[@id="{}"]'.format(key))
        node = Node(element)
        source = node.source
        successors = [x for x in self if x.source == key]
        for node in successors:
            node.source = source
        self.tree.remove(element)
    
    def __str__(self):
        self.__optimize_appearance()
        rough_string = ET.tostring(self.tree, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        return reparsed.toprettyxml(indent='\t', newl='')
    
    def __iter__(self):
        return iter(self.nodes())
    
    def successors(self, id, recursive=False):
        """
        find the succeeding node(s) of a node
        
        Parameters
        ----------
        id: str
            the ID of the node
        recursive: bool
            find successors recursively?

        Returns
        -------
        list of str
            the ID(s) of the successors
        """
        if not isinstance(id, str):
            raise TypeError("'id' must be of type 'str', is {}".format(type(id)))
        successors = []
        for node in self:
            if node.source == id or (isinstance(node.source, list) and id in node.source):
                successors.append(node.id)
        if recursive:
            for item in successors:
                new = self.successors(item, recursive=True)
                successors.extend(new)
            successors = list(set(successors))
        return successors
    
    def __reset_successor_source(self, id):
        """
        reset the sources of nodes to that of a newly inserted one
        
        Parameters
        ----------
        id: str
            the ID of the newly inserted node

        Returns
        -------

        """
        
        def reset(id, source, excludes=None):
            if isinstance(source, list):
                for item in source:
                    successors = self.successors(item)
                    excludes = [x for x in successors if x in source]
                    reset(id, item, excludes)
            else:
                try:
                    # find the source nodes of the current node
                    if source is not None:
                        successors = self.successors(source)
                    else:
                        return  # nothing to reset
                    # delete the ID of the current node from the successors
                    if id in successors:
                        del successors[successors.index(id)]
                    if excludes is not None:
                        for item in excludes:
                            del successors[successors.index(item)]
                    for successor in successors:
                        successor_source = self[successor].source
                        if isinstance(successor_source, list):
                            successor_source[successor_source.index(source)] = id
                            self[successor].source = successor_source
                        else:
                            self[successor].source = id
                except IndexError:
                    # case where no successor exists because the new node
                    # is the new last node in the graph
                    pass
                except RuntimeError:
                    # case where the successor node is of type Read
                    pass
        
        reset(id, self[id].source)
    
    def __optimize_appearance(self):
        """
        assign grid coordinates to the nodes for display in the SNAP GraphBuilder GUI
        
        This method is applied by :meth:`__str__` for the final formatting of the XML text representation
        
        Returns
        -------

        """
        layout = self.tree.find('.//applicationData[@id="Presentation"]')
        
        counter = 0
        x = 5
        for id in self.ids:
            pres = layout.find('.//node[@id="{}"]'.format(id))
            y = 20. if counter % 2 == 0 else 160.
            if pres is None:
                pres = ET.SubElement(layout, 'node', {'id': id})
                pos = ET.SubElement(pres, 'displayPosition',
                                    {'x': "{}".format(x), 'y': "{}".format(y)})
            else:
                pres.find('displayPosition').attrib['x'] = "{}".format(x)
                pres.find('displayPosition').attrib['y'] = "{}".format(y)
            counter += 1
            x += len(id) * 8
    
    @property
    def ids(self):
        """
        
        Returns
        -------
        list
            the IDs of all nodes
        """
        return [node.id for node in self]
    
    def index(self, node):
        """
        
        Parameters
        ----------
        node: Node
            a node in the workflow

        Returns
        -------
        int
            the index position of the node in the workflow
        """
        return list(self.tree).index(node.element)
    
    def insert_node(self, node, before=None, after=None, resetSuccessorSource=True, void=True):
        """
        insert one or multiple node(s) into the workflow including setting the source to the predecessor
        and setting the ID as source of the successor.
        
        Parameters
        ----------
        node: Node or list[Node]
            the node(s) to be inserted
        before: Node, str or list
            a Node object; the ID(s) of the node(s) before the newly inserted node; a list of node IDs is intended for
            nodes that require multiple sources, e.g. sliceAssembly
        after: Node, str
            a Node object; the ID of the node after the newly inserted node
        resetSuccessorSource: bool
            reset the source of the successor node to the ID of the newly inserted node?
        void: bool
            if false, the function returns the node

        Returns
        -------
        Node or list[Node] or None
            the new node, a list of nodes, or None, depending on the `node` input and argument `void`
        """
        if isinstance(node, list):
            self.insert_node(node=node[0], before=before, after=after,
                             resetSuccessorSource=resetSuccessorSource, void=True)
            for i, item in enumerate(node[1:]):
                self.insert_node(node=item, before=node[i].id,
                                 resetSuccessorSource=resetSuccessorSource, void=True)
        else:
            ncopies = [x.operator for x in self.nodes()].count(node.operator)
            if ncopies > 0:
                node.id = '{0} ({1})'.format(node.operator, ncopies + 1)
            else:
                node.id = node.operator
            
            if isinstance(before, Node):
                before = before.id
            if isinstance(after, Node):
                after = after.id
            
            if before is None and after is None and len(self) > 0:
                before = self[len(self) - 1].id
            if before and not after:
                if isinstance(before, list):
                    indices = [self.index(self[x]) for x in before]
                    predecessor = self[before[indices.index(max(indices))]]
                else:
                    predecessor = self[before]
                log.debug('inserting node {} after {}'.format(node.id, predecessor.id))
                position = self.index(predecessor) + 1
                self.tree.insert(position, node.element)
                newnode = Node(self.tree[position])
                ####################################################
                # set the source product for the new node
                if newnode.operator != 'Read':
                    newnode.source = before
                ####################################################
                # set the source product for the node after the new node
                if resetSuccessorSource:
                    self.__reset_successor_source(newnode.id)
            ########################################################
            elif after and not before:
                successor = self[after]
                log.debug('inserting node {} before {}'.format(node.id, successor.id))
                position = self.index(successor)
                self.tree.insert(position, node.element)
                newnode = Node(self.tree[position])
                ####################################################
                # set the source product for the new node
                if newnode.operator != 'Read':
                    source = successor.source
                    newnode.source = source
                ####################################################
                # set the source product for the node after the new node
                if resetSuccessorSource:
                    self[after].source = newnode.id
            else:
                log.debug('inserting node {}'.format(node.id))
                self.tree.insert(len(self.tree) - 1, node.element)
        if not void:
            return node
    
    def nodes(self):
        """
        
        Returns
        -------
        list[Node]
            the list of :class:`Node` objects in the workflow
        """
        return [Node(x) for x in self.tree.findall('node')]
    
    @property
    def operators(self):
        """
        
        Returns
        -------
        list
            the names of the unique operators in the workflow
        """
        return sorted(list(set([node.operator for node in self])))
    
    def refresh_ids(self):
        """
        Ensure unique IDs for all nodes. If two nodes with the same ID are found one is renamed to "ID (2)".
        E.g. 2 x "Write" -> "Write", "Write (2)".
        This method is no longer used and is just kept in case there is need for it in the future.
        
        Returns
        -------

        """
        counter = {}
        for node in self:
            operator = node.operator
            if operator not in counter.keys():
                counter[operator] = 1
            else:
                counter[operator] += 1
            if counter[operator] > 1:
                new = '{} ({})'.format(operator, counter[operator])
            else:
                new = operator
            if node.id != new:
                log.debug('renaming node {} to {}'.format(node.id, new))
                node.id = new
    
    def set_par(self, key, value, exceptions=None):
        """
        set a parameter for all nodes in the workflow
        
        Parameters
        ----------
        key: str
            the parameter name
        value: bool or int or float or str
            the parameter value
        exceptions: list
            a list of node IDs whose parameters should not be changed

        Returns
        -------

        """
        for node in self:
            if exceptions is not None and node.id in exceptions:
                continue
            if key in node.parameters.keys():
                node.parameters[key] = value2str(value)
    
    def suffix(self, stop=None):
        """
        Get the SNAP operator suffix sequence
        
        Parameters
        ----------
        stop: str
            the ID of the last workflow node
        
        Returns
        -------
        str
            a file suffix created from the order of which the nodes will be executed
        """
        nodes = self.tree.findall('node')
        names = [re.sub(r'[ ]*\([0-9]+\)', '', y.attrib['id']) for y in nodes]
        names_unique = []
        for name in names:
            if name not in names_unique:
                names_unique.append(name)
            if name == stop:
                break
        config = ExamineSnap()
        suffix = '_'.join(filter(None, [config.get_suffix(x) for x in names_unique]))
        return suffix
    
    def write(self, outfile):
        """
        write the workflow to an XML file
        
        Parameters
        ----------
        outfile: str
            the name of the file to write

        Returns
        -------

        """
        outfile = outfile if outfile.endswith('.xml') else outfile + '.xml'
        log.debug('writing {}'.format(outfile))
        with open(outfile, 'w') as out:
            out.write(self.__str__())


class Node(object):
    """
    class for handling of SNAP workflow processing nodes
    
    Parameters
    ----------
    element: ~xml.etree.ElementTree.Element
        the node XML element
    """
    
    def __init__(self, element):
        if not isinstance(element, ET.Element):
            raise TypeError('element must be of type xml.etree.ElementTree.Element')
        self.element = element
    
    def __repr__(self):
        return "pyroSAR Node object '{}'".format(self.id)
    
    def __str__(self):
        rough_string = ET.tostring(self.element, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        return reparsed.toprettyxml(indent='\t', newl='')
    
    def __set_source(self, key, value):
        source = self.element.find('.//sources/{}'.format(key))
        if source is None:
            child = ET.SubElement(self.element.find('.//sources'),
                                  key, {'refid': value})
        else:
            source.attrib['refid'] = value
    
    def copy(self):
        """
        
        Returns
        -------
        Node
            a copy of the Node object
        """
        return Node(copy.deepcopy(self.element))
    
    @property
    def id(self):
        """
        
        Returns
        -------
        str
            the node ID
        """
        return self.element.attrib['id']
    
    @id.setter
    def id(self, value):
        self.element.attrib['id'] = value
    
    @property
    def operator(self):
        """
        
        Returns
        -------
        str
            the name of the node's processing operator
        """
        return self.element.find('.//operator').text
    
    @property
    def parameters(self):
        """
        
        Returns
        -------
        Par or Par_BandMath
            the processing parameters of the node
        """
        params = self.element.find('.//parameters')
        if self.operator == 'BandMaths':
            return Par_BandMath(operator=self.operator, element=params)
        else:
            return Par(operator=self.operator, element=params)
    
    @property
    def source(self):
        """
        
        Returns
        -------
        str or list
            the ID(s) of the source node(s)
        """
        sources = []
        elements = self.element.findall('.//sources/')
        for element in elements:
            if element.tag.startswith('sourceProduct'):
                sources.append(element.attrib['refid'])
        
        if len(sources) == 0:
            return None
        elif len(sources) == 1:
            return sources[0]
        else:
            return sources
    
    @source.setter
    def source(self, value):
        """
        reset the source of the node by ID
        
        Parameters
        ----------
        value: str or list
            the ID(s) of the new source node(s)

        Returns
        -------
        
        Raises
        ------
        RuntimeError
        """
        if isinstance(value, list) and len(value) == 1:
            value = value[0]
        log.debug('setting the source of node {} to {}'.format(self.id, value))
        if isinstance(value, str):
            if isinstance(self.source, list):
                raise TypeError(
                    'node {} has multiple sources, which must be reset using a list, not str'.format(self.id))
            self.__set_source('sourceProduct', value)
        elif isinstance(value, list):
            key = 'sourceProduct'
            for i, item in enumerate(value):
                self.__set_source(key, item)
                key = 'sourceProduct.{}'.format(i + 1)


class Par(object):
    """
    class for handling processing node parameters
    
    Parameters
    ----------
    operator: str
        the name of the SNAP Node operator
    element: ~xml.etree.ElementTree.Element
        the node parameter XML element
    """
    
    def __init__(self, operator, element):
        self.operator = operator
        self.__element = element
    
    def __delitem__(self, key):
        par = self.__element.find('.//{}'.format(key))
        self.__element.remove(par)
    
    def __getitem__(self, item):
        """
        
        Parameters
        ----------
        item

        Returns
        -------
        str
        """
        if item not in self.keys():
            raise KeyError('key {} does not exist'.format(item))
        return self.__element.find('.//{}'.format(item)).text
    
    def __setitem__(self, key, value):
        if key not in self.keys():
            raise KeyError("unknown key for node '{}': '{}'".format(self.operator, key))
        strval = value2str(value)
        self.__element.find('.//{}'.format(key)).text = strval
    
    def __repr__(self):
        return str(self.dict())
    
    def dict(self):
        """
        
        Returns
        -------
        dict
            the parameters as a dictionary
        """
        return dict(self.items())
    
    def items(self):
        """
        
        Returns
        -------
        list
            the parameters as (key, value) as from :meth:`dict.items()`
        """
        return list(zip(self.keys(), self.values()))
    
    def keys(self):
        """
        
        Returns
        -------
        list
            the parameter names as from :meth:`dict.keys()`
        """
        return [x.tag for x in self.__element.findall('./')]
    
    def values(self):
        """
        
        Returns
        -------
        list
            the parameter values as from :meth:`dict.values()`
        """
        return [x.text for x in self.__element.findall('./')]


class Par_BandMath(Par):
    """
    class for handling BandMaths node parameters

    Parameters
    ----------
    element: ~xml.etree.ElementTree.Element
        the node parameter XML element
    """
    
    def __init__(self, operator, element):
        self.operator = operator
        self.__element = element
        super(Par_BandMath, self).__init__(operator, element)
    
    def __getitem__(self, item):
        if item in ['variables', 'targetBands']:
            out = []
            for x in self.__element.findall('.//{}'.format(item[:-1])):
                out.append(Par(self.operator, x))
            return out
        else:
            raise ValueError("can only get items 'variables' and 'targetBands'")
    
    def clear_variables(self):
        """
        remove all `variables` elements from the node
        
        Returns
        -------

        """
        var = self.__element.find('.//variables')
        for item in var:
            var.remove(item)
    
    def add_equation(self):
        """
        add an equation element to the node
        
        Returns
        -------

        """
        eqs = self.__element.find('.//targetBands')
        eqlist = eqs.findall('.//targetBand')
        eq1 = eqlist[0]
        eq2 = copy.deepcopy(eq1)
        for item in eq2:
            item.text = None
        eqs.insert(len(eqlist), eq2)


def value2str(value):
    """
    format a parameter value to string to be inserted into a workflow
    
    Parameters
    ----------
    value: bool, int, float, list

    Returns
    -------
    str
        the string representation of the value
    """
    if isinstance(value, bool):
        strval = str(value).lower()
    elif isinstance(value, list):
        strval = ','.join(map(str, value))
    elif value is None:
        strval = value
    else:
        strval = str(value)
    return strval


def erode_edges(src, only_boundary=False, connectedness=4, pixels=1):
    """
    Erode noisy edge pixels in SNAP-processed images.
    It was discovered that images contain border pixel artifacts after `Terrain-Correction`.
    Likely this is coming from treating the value 0 as regular value instead of no data during resampling.
    This function erodes these edge pixels using :func:`scipy.ndimage.binary_erosion`.
    scipy is not a base dependency of pyroSAR and has to be installed separately.
    
    .. figure:: figures/snap_erode_edges.png
        :align: center
        
        VV gamma0 RTC backscatter image visualizing the noisy border (left) and the cleaned result (right).
        The area covers approx. 2.3 x 2.3 km². Pixel spacing is 20 m. connectedness 4, 1 pixel.
    
    Parameters
    ----------
    src: str
        a processed SAR image in BEAM-DIMAP format (.dim), a single .img file (ENVI format) or a
        directory with .img files. 0 is assumed as no data value.
    only_boundary: bool
        only erode edges at the image boundary (or also at data gaps caused by e.g. masking during Terrain-Flattening)?
    connectedness: int
        the number of pixel neighbors considered for the erosion. Either 4 or 8, translating to a
        :func:`scipy.ndimage.generate_binary_structure` `connectivity` of 1 or 2, respectively.
    pixels: int
        the number of pixels to erode from the edges. Directly translates to `iterations` of
        :func:`scipy.ndimage.iterate_structure`.
    
    Returns
    -------

    """
    images = None
    if src.endswith('.dim'):
        workdir = src.replace('.dim', '.data')
    elif src.endswith('.img'):
        images = [src]
        workdir = None
    elif os.path.isdir(src):
        workdir = src
    else:
        raise RuntimeError("'src' must be either a file in BEAM-DIMAP format (extension '.dim'), "
                           "an ENVI file with extension *.img, or a directory.")
    
    if images is None:
        images = [x for x in finder(workdir, ['*.img'], recursive=False)
                  if 'layoverShadowMask' not in x]
    if len(images) == 0:
        raise RuntimeError("could not find any files with extension '.img'")
    
    from scipy.ndimage import binary_erosion, generate_binary_structure, iterate_structure
    
    if connectedness == 4:
        connectivity = 1
    elif connectedness == 8:
        connectivity = 2
    else:
        raise ValueError('connectedness must be either 4 or 8')
    
    structure = generate_binary_structure(rank=2, connectivity=connectivity)
    if pixels > 1:
        structure = iterate_structure(structure=structure, iterations=pixels)
    
    if workdir is not None:
        fname_mask = os.path.join(workdir, 'datamask.tif')
    else:
        fname_mask = os.path.join(os.path.dirname(src), 'datamask.tif')
    write_intermediates = False  # this is intended for debugging
    
    def erosion(src, dst, structure, only_boundary, write_intermediates=False):
        with Raster(src) as ref:
            array = ref.array()
            if not os.path.isfile(dst):
                mask = array != 0
                # do not perform erosion if data only contains nodata (mask == 1)
                if len(mask[mask == 1]) == 0:
                    ref.write(outname=dst, array=mask, dtype='Byte',
                              options=['COMPRESS=DEFLATE'])
                    return array, mask
                if write_intermediates:
                    ref.write(dst.replace('.tif', '_init.tif'),
                              array=mask, dtype='Byte',
                              options=['COMPRESS=DEFLATE'])
                if only_boundary:
                    with vectorize(target=mask, reference=ref) as vec:
                        with boundary(vec, expression="value=1") as bounds:
                            with rasterize(vectorobject=bounds, reference=ref, nodata=None) as new:
                                mask = new.array()
                                if write_intermediates:
                                    vec.write(dst.replace('.tif', '_init_vectorized.gpkg'))
                                    bounds.write(dst.replace('.tif', '_boundary_vectorized.gpkg'))
                                    new.write(outname=dst.replace('.tif', '_boundary.tif'),
                                              dtype='Byte', options=['COMPRESS=DEFLATE'])
                mask = binary_erosion(input=mask, structure=structure)
                ref.write(outname=dst, array=mask, dtype='Byte',
                          options=['COMPRESS=DEFLATE'])
            else:
                with Raster(dst) as ras:
                    mask = ras.array()
        array[mask == 0] = 0
        return array, mask
    
    # make sure a backscatter image is used for creating the mask
    backscatter = [x for x in images if re.search('^(?:Sigma0_|Gamma0_|C11|C22)', os.path.basename(x))]
    images.insert(0, images.pop(images.index(backscatter[0])))
    
    mask = None
    for img in images:
        if mask is None:
            array, mask = erosion(src=img, dst=fname_mask,
                                  structure=structure, only_boundary=only_boundary,
                                  write_intermediates=write_intermediates)
        else:
            with Raster(img) as ras:
                array = ras.array()
            array[mask == 0] = 0
        # do not apply mask if it only contains 1 (valid data)
        if len(mask[mask == 0]) == 0:
            break
        ras = gdal.Open(img, GA_Update)
        band = ras.GetRasterBand(1)
        band.WriteArray(array)
        band.FlushCache()
        band = None
        ras = None


def mli_parametrize(scene, spacing=None, rlks=None, azlks=None, **kwargs):
    """
    Convenience function for parametrizing a `Multilook` node.
    
    Parameters
    ----------
    scene: pyroSAR.drivers.ID
        The SAR scene to be processed
    spacing: int or float or None
        the target pixel spacing for automatic determination of looks using function
        :func:`~pyroSAR.ancillary.multilook_factors`. Overridden by arguments `rlks` and `azlks` if they are not None.
    rlks: int or None
        the number of range looks
    azlks: int or None
        the number of azimuth looks
    bands: list[str] or None
        an optional list of bands names
    kwargs
        further keyword arguments for node parametrization. Known options:
        
         - grSquarePixel
         - outputIntensity
         - sourceBands
    
    Returns
    -------
    Node or None
        either a `Node` object if multilooking is necessary (either `rlks` or `azlks` are greater than 1) or None.
    
    See Also
    --------
    pyroSAR.ancillary.multilook_factors
    """
    try:
        image_geometry = scene.meta['image_geometry']
        incidence = scene.meta['incidence']
    except KeyError:
        msg = 'This function does not yet support {} products in {} format'
        raise RuntimeError(msg.format(scene.sensor, scene.__class__.__name__))
    
    if rlks is None and azlks is None:
        if spacing is None:
            raise RuntimeError("either 'spacing' or 'rlks' and 'azlks' must set to numeric values")
        rlks, azlks = multilook_factors(source_rg=scene.spacing[0],
                                        source_az=scene.spacing[1],
                                        target=spacing,
                                        geometry=image_geometry,
                                        incidence=incidence)
    if [rlks, azlks].count(None) > 0:
        raise RuntimeError("'rlks' and 'azlks' must either both be integers or None")
    
    if azlks > 1 or rlks > 1 or scene.sensor in ['ERS1', 'ERS2', 'ASAR']:
        ml = parse_node('Multilook')
        ml.parameters['nAzLooks'] = azlks
        ml.parameters['nRgLooks'] = rlks
        for key, val in kwargs.items():
            ml.parameters[key] = val
        return ml


def orb_parametrize(scene, formatName, allow_RES_OSV=True, url_option=1, **kwargs):
    """
    convenience function for parametrizing an `Apply-Orbit-File`.
    Required Sentinel-1 orbit files are directly downloaded.
    
    Parameters
    ----------
    scene: pyroSAR.drivers.ID
        The SAR scene to be processed
    workflow: Workflow
        the SNAP workflow object
    before: str
        the ID of the node after which the `Apply-Orbit-File` node will be inserted
    formatName: str
        the scene's data format
    allow_RES_OSV: bool
        (only applies to Sentinel-1) Also allow the less accurate RES orbit files to be used?
    url_option: int
        the OSV download URL option; see :meth:`pyroSAR.S1.OSV.catch`
    kwargs
        further keyword arguments for node parametrization. Known options:
        
         - continueOnFail
         - polyDegree
    
    Returns
    -------
    Node
        the Apply-Orbit-File node object
    """
    orbitType = None
    orbit_lookup = {'SENTINEL-1': 'Sentinel Precise (Auto Download)'}
    if formatName in orbit_lookup:
        orbitType = orbit_lookup[formatName]
    if formatName == 'ENVISAT':  # ASAR, ERS1, ERS2
        if scene.sensor == 'ASAR':
            orbitType = 'DORIS Precise VOR (ENVISAT) (Auto Download)'
        else:
            # Another option for ERS is 'DELFT Precise (ENVISAT, ERS1&2) (Auto Download)'.
            # Neither option is suitable for all products, and auto-selection can
            # only happen once a downloader (similar to S1.auxil.OSV) is written.
            orbitType = 'PRARE Precise (ERS1&2) (Auto Download)'
    if orbitType is None:
        raise RuntimeError(f'Could not determine orbit type for {formatName} format')
    
    if formatName == 'SENTINEL-1':
        osv_type = ['POE', 'RES'] if allow_RES_OSV else 'POE'
        match = scene.getOSV(osvType=osv_type, returnMatch=True, url_option=url_option)
        if match is None and allow_RES_OSV:
            scene.getOSV(osvType='RES', url_option=url_option)
            orbitType = 'Sentinel Restituted (Auto Download)'
    
    orb = parse_node('Apply-Orbit-File')
    orb.parameters['orbitType'] = orbitType
    for key, val in kwargs.items():
        orb.parameters[key] = val
    return orb


def sub_parametrize(scene, geometry=None, offset=None, buffer=0.01, copyMetadata=True, **kwargs):
    """
    convenience function for parametrizing an `Subset` node.
    
    Parameters
    ----------
    scene: pyroSAR.drivers.ID
        The SAR scene to be processed
    geometry: dict or spatialist.vector.Vector or str or None
        A vector geometry for geographic subsetting (node parameter geoRegion):
        
         - :class:`~spatialist.vector.Vector`: a vector object in arbitrary CRS
         - :class:`str`: a name of a file that can be read with :class:`~spatialist.vector.Vector` in arbitrary CRS
         - :class:`dict`: a dictionary with keys `xmin`, `xmax`, `ymin`, `ymax` in EPSG:4326 coordinates
    offset: tuple or None
        a tuple with pixel coordinates as (left, right, top, bottom)
    buffer: int or float
        an additional buffer in degrees to add around the `geometry`
    copyMetadata: bool
        copy the metadata of the source product?
    kwargs
        further keyword arguments for node parametrization. Known options:
        
         - fullSwath
         - referenceBand
         - sourceBands
         - subSamplingX
         - subSamplingY
         - tiePointGrids

    Returns
    -------
    Node
        the Subset node object
    """
    subset = parse_node('Subset')
    if geometry:
        if isinstance(geometry, dict):
            ext = geometry
        else:
            if isinstance(geometry, Vector):
                shp = geometry.clone()
            elif isinstance(geometry, str):
                shp = Vector(geometry)
            else:
                raise TypeError("argument 'geometry' must be either a dictionary, a Vector object or a filename.")
            # reproject the geometry to WGS 84 latlon
            shp.reproject(4326)
            ext = shp.extent
            shp.close()
        # add an extra buffer
        ext['xmin'] -= buffer
        ext['ymin'] -= buffer
        ext['xmax'] += buffer
        ext['ymax'] += buffer
        with bbox(ext, 4326) as bounds:
            inter = intersect(scene.bbox(), bounds)
            if not inter:
                raise RuntimeError('no bounding box intersection between shapefile and scene')
            inter.close()
            wkt = bounds.convert2wkt()[0]
        subset.parameters['region'] = [0, 0, scene.samples, scene.lines]
        subset.parameters['geoRegion'] = wkt
    #######################
    # (optionally) configure Subset node for pixel offsets
    elif offset and not geometry:
        # left, right, top and bottom offset in pixels
        l, r, t, b = offset
        
        subset_values = [l, t, scene.samples - l - r, scene.lines - t - b]
        subset.parameters['region'] = subset_values
        subset.parameters['geoRegion'] = ''
    else:
        raise RuntimeError("one of 'geometry' and 'offset' must be set")
    
    subset.parameters['copyMetadata'] = copyMetadata
    for key, val in kwargs.items():
        subset.parameters[key] = val
    return subset


def geo_parametrize(spacing, t_srs, tc_method='Range-Doppler',
                    sourceBands=None, demName='SRTM 1Sec HGT', externalDEMFile=None,
                    externalDEMNoDataValue=None, externalDEMApplyEGM=True,
                    alignToStandardGrid=False, standardGridAreaOrPoint='point',
                    standardGridOriginX=0, standardGridOriginY=0,
                    nodataValueAtSea=False, export_extra=None,
                    demResamplingMethod='BILINEAR_INTERPOLATION',
                    imgResamplingMethod='BILINEAR_INTERPOLATION',
                    **kwargs):
    """
    convenience function for parametrizing geocoding nodes.
    
    Parameters
    ----------
    workflow: Workflow
        the SNAP workflow object
    before: str
        the ID of the node after which the terrain correction node will be inserted
    tc_method: str
        the terrain correction method. Supported options:
        
         - Range-Doppler (SNAP node `Terrain-Correction`)
         - SAR simulation cross correlation
           (SNAP nodes `SAR-Simulation`->`Cross-Correlation`->`Warp`->`Terrain-Correction`)
    
    sourceBands: List[str] or None
        the image band names to geocode; default None: geocode all incoming bands.
    spacing: int or float
        The target pixel spacing in meters.
    t_srs: int or str or osgeo.osr.SpatialReference
        A target geographic reference system in WKT, EPSG, PROJ4 or OPENGIS format.
        See function :func:`spatialist.auxil.crsConvert()` for details.
    demName: str
        The name of the auto-download DEM. Default is 'SRTM 1Sec HGT'. Is ignored when `externalDEMFile` is not None.
        Supported options:
        
         - ACE2_5Min
         - ACE30
         - ASTER 1sec GDEM
         - CDEM
         - Copernicus 30m Global DEM
         - Copernicus 90m Global DEM
         - GETASSE30
         - SRTM 1Sec Grid
         - SRTM 1Sec HGT
         - SRTM 3Sec
    externalDEMFile: str or None, optional
        The absolute path to an external DEM file. Default is None. Overrides `demName`.
    externalDEMNoDataValue: int, float or None, optional
        The no data value of the external DEM. If not specified (default) the function will try to read it from the
        specified external DEM.
    externalDEMApplyEGM: bool, optional
        Apply Earth Gravitational Model to external DEM? Default is True.
    alignToStandardGrid: bool
        Align all processed images to a common grid?
    standardGridAreaOrPoint: str
        treat alignment coordinate as pixel center ('point', SNAP default) or upper left ('area').
    standardGridOriginX: int or float
        The x origin value for grid alignment
    standardGridOriginY: int or float
        The y origin value for grid alignment
    nodataValueAtSea:bool
        mask values over sea?
    export_extra: list[str] or None
        a list of ancillary layers to write. Supported options:
        
         - DEM
         - latLon
         - incidenceAngleFromEllipsoid
         - layoverShadowMask
         - localIncidenceAngle
         - projectedLocalIncidenceAngle
         - selectedSourceBand
    demResamplingMethod: str
        the DEM resampling method
    imgResamplingMethod: str
        the image resampling method
    kwargs
        further keyword arguments for node parametrization. Known options:
        
         - outputComplex
         - applyRadiometricNormalization
         - saveSigmaNought
         - saveGammaNought
         - saveBetaNought
         - incidenceAngleForSigma0
         - incidenceAngleForGamma0
         - auxFile
         - externalAuxFile
         - openShiftsFile (SAR simulation cross correlation only)
         - openResidualsFile (SAR simulation cross correlation only)
    
    Returns
    -------
    Node or list[Node]
        the Terrain-Correction node object or a list containing the objects for SAR-Simulation,
        Cross-Correlation and SARSim-Terrain-Correction.
    """
    tc = parse_node('Terrain-Correction')
    tc.parameters['nodataValueAtSea'] = nodataValueAtSea
    
    if tc_method == 'Range-Doppler':
        tc.parameters['sourceBands'] = sourceBands
        sarsim = None
        out = tc
        dem_nodes = [tc]
    elif tc_method == 'SAR simulation cross correlation':
        sarsim = parse_node('SAR-Simulation')
        sarsim.parameters['sourceBands'] = sourceBands
        cc = parse_node('Cross-Correlation')
        cc.parameters['coarseRegistrationWindowWidth'] = 64
        cc.parameters['coarseRegistrationWindowHeight'] = 64
        cc.parameters['maxIteration'] = 2
        cc.parameters['onlyGCPsOnLand'] = True
        warp = parse_node('Warp')
        dem_nodes = [sarsim, tc]
        out = [sarsim, cc, warp, tc]
    else:
        raise RuntimeError(f'tc_method not recognized: "{tc_method}"')
    
    tc.parameters['imgResamplingMethod'] = imgResamplingMethod
    
    if standardGridAreaOrPoint == 'area':
        standardGridOriginX -= spacing / 2
        standardGridOriginY += spacing / 2
    tc.parameters['alignToStandardGrid'] = alignToStandardGrid
    tc.parameters['standardGridOriginX'] = standardGridOriginX
    tc.parameters['standardGridOriginY'] = standardGridOriginY
    
    # specify spatial resolution and coordinate reference system of the output dataset
    tc.parameters['pixelSpacingInMeter'] = spacing
    
    try:
        # try to convert the CRS into EPSG code (for readability in the workflow XML)
        t_srs = crsConvert(t_srs, 'epsg')
    except TypeError:
        raise RuntimeError("format of parameter 't_srs' not recognized")
    except RuntimeError:
        # This error can occur when the CRS does not have a corresponding EPSG code.
        # In this case the original CRS representation is written to the workflow.
        pass
    
    # The EPSG code 4326 is not supported by SNAP and thus the WKT string has to be defined.
    # In all other cases defining EPSG:{code} will do.
    if t_srs == 4326:
        t_srs = 'GEOGCS["WGS84(DD)",' \
                'DATUM["WGS84",' \
                'SPHEROID["WGS84", 6378137.0, 298.257223563]],' \
                'PRIMEM["Greenwich", 0.0],' \
                'UNIT["degree", 0.017453292519943295],' \
                'AXIS["Geodetic longitude", EAST],' \
                'AXIS["Geodetic latitude", NORTH]]'
    
    if isinstance(t_srs, int):
        t_srs = 'EPSG:{}'.format(t_srs)
    
    tc.parameters['mapProjection'] = t_srs
    
    export_extra_options = \
        ['DEM', 'latLon',
         'incidenceAngleFromEllipsoid',
         'layoverShadowMask',
         'localIncidenceAngle',
         'projectedLocalIncidenceAngle',
         'selectedSourceBand']
    if export_extra is not None:
        for item in export_extra:
            if item in export_extra_options:
                key = f'save{item[0].upper()}{item[1:]}'
                tc.parameters[key] = True
    
    for dem_node in dem_nodes:
        dem_parametrize(node=dem_node, demName=demName,
                        externalDEMFile=externalDEMFile,
                        externalDEMNoDataValue=externalDEMNoDataValue,
                        externalDEMApplyEGM=externalDEMApplyEGM,
                        demResamplingMethod=demResamplingMethod)
    
    for key, val in kwargs.items():
        tc.parameters[key] = val
    return out


def dem_parametrize(workflow=None, node=None, demName='SRTM 1Sec HGT', externalDEMFile=None,
                    externalDEMNoDataValue=None, externalDEMApplyEGM=False,
                    demResamplingMethod='BILINEAR_INTERPOLATION'):
    """
    DEM parametrization for a full workflow or a single node. In the former case, all nodes with the
    DEM-relevant parameters can be modified at once, e.g. `Terrain-Flattening` and `Terrain-Correction`.
    
    Parameters
    ----------
    workflow: Workflow or None
        a SNAP workflow object
    node: Node or None
        a SNAP node object
    demName: str
        The name of the auto-download DEM. Default is 'SRTM 1Sec HGT'. Is ignored when `externalDEMFile` is not None.
        Supported options:
        
         - ACE2_5Min
         - ACE30
         - ASTER 1sec GDEM
         - CDEM
         - Copernicus 30m Global DEM
         - Copernicus 90m Global DEM
         - GETASSE30
         - SRTM 1Sec Grid
         - SRTM 1Sec HGT
         - SRTM 3Sec
    externalDEMFile: str or None, optional
        The absolute path to an external DEM file. Default is None. Overrides `demName`.
    externalDEMNoDataValue: int, float or None, optional
        The no data value of the external DEM. If not specified (default) the function will try to read it from the
        specified external DEM.
    externalDEMApplyEGM: bool, optional
        Apply Earth Gravitational Model to external DEM? Default is True.
    demResamplingMethod: str
        the DEM resampling method

    Returns
    -------

    """
    # select DEM type
    dempar = {'externalDEMFile': externalDEMFile,
              'externalDEMApplyEGM': externalDEMApplyEGM,
              'demResamplingMethod': demResamplingMethod}
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
        dempar['demName'] = demName
        dempar['externalDEMFile'] = None
        dempar['externalDEMNoDataValue'] = 0
    
    if workflow is not None:
        for key, value in dempar.items():
            workflow.set_par(key, value)
    elif node is not None:
        for key, value in dempar.items():
            if key in node.parameters.keys():
                node.parameters[key] = value
    else:
        raise RuntimeError("either 'workflow' or 'node must be defined'")
    
    # download the EGM lookup table if necessary
    if dempar['externalDEMApplyEGM']:
        get_egm_lookup(geoid='EGM96', software='SNAP')
