###############################################################################
# pyroSAR SNAP API tools

# Copyright (c) 2017-2022, the pyroSAR Developers.

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
import subprocess as sp
from xml.dom import minidom
import xml.etree.ElementTree as ET

from pyroSAR import identify
from pyroSAR.examine import ExamineSnap
from pyroSAR.ancillary import windows_fileprefix

from spatialist import Raster, vectorize, rasterize, boundary
from spatialist.auxil import gdal_translate
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
        use an existing XML text file or force re-parsing the gpt docstring and overwriting the XML file?

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
    version = snap.get_version('s1tbx')['version']
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
    
    if not os.path.isfile(absname) or not use_existing:
        gpt = snap.gpt
        
        cmd = [gpt, operator, '-h']
        
        out, err = run(cmd=cmd, void=False)
        
        if re.search('Unknown operator', out + err):
            raise RuntimeError("unknown operator '{}'".format(operator))
        
        graph = re.search('<graph id.*', out, flags=re.DOTALL).group()
        graph = re.sub(r'>\${.*', '/>', graph)  # remove placeholder values like ${value}
        graph = re.sub(r'<\.\.\./>.*', '', graph)  # remove <.../> placeholders
        if operator == 'BandMaths':
            graph = graph.replace('sourceProducts', 'sourceProduct')
        tree = ET.fromstring(graph)
        for elt in tree.iter():
            if elt.text in ['string', 'double', 'integer', 'float']:
                elt.text = None
        node = tree.find('node')
        node.attrib['id'] = operator
        # add a second source product entry for multi-source nodes
        # multi-source nodes are those with an entry 'sourceProducts' instead of 'sourceProduct'
        # exceptions are registered in this list:
        multisource = ['Back-Geocoding']
        if operator != 'Read' and operator != 'ProductSet-Reader':
            source = node.find('.//sources')
            child = source[0]
            if child.tag == 'sourceProducts' or operator in multisource:
                child2 = ET.SubElement(source, 'sourceProduct.1', {'refid': 'Read (2)'})
            child.tag = 'sourceProduct'
            child.attrib['refid'] = 'Read'
            child.text = None
        if operator == 'BandMaths':
            tband = tree.find('.//targetBand')
            for item in ['spectralWavelength', 'spectralBandwidth',
                         'scalingOffset', 'scalingFactor',
                         'validExpression', 'spectralBandIndex']:
                el = tband.find('.//{}'.format(item))
                tband.remove(el)
        tree.find('.//parameters').set('class', 'com.bc.ceres.binding.dom.XppDomElement')
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
    Depending on the parametrization this function might create two sub-directories in `tmpdir`,
    carrying a suffix  \*_bnr for S1 GRD border noise removal and \*_sub for sub-workflows and their
    intermediate outputs. Both are deleted if ``cleanup=True``. If `tmpdir` is empty afterwards, it is also deleted.
    
    Parameters
    ----------
    xmlfile: str
        the name of the workflow XML file
    tmpdir: str
        a temporary directory for storing intermediate files
    groups: list[list[str]] or None
        a list of lists each containing IDs for individual nodes
    cleanup: bool
        should all files written to the temporary directory during function execution be deleted after processing?
    gpt_exceptions: dict or None
        a dictionary to override the configured GPT executable for certain operators;
        each (sub-)workflow containing this operator will be executed with the define executable;
        
         - e.g. ``{'Terrain-Flattening': '/home/user/snap/bin/gpt'}``
    gpt_args: list or None
        a list of additional arguments to be passed to the gpt call
        
        - e.g. ``['-x', '-c', '2048M']`` for increased tile cache size and intermediate clearing
    removeS1BorderNoiseMethod: str
        the border noise removal method to be applied, See :func:`pyroSAR.S1.removeGRDBorderNoise` for details;
        one of the following:
        
         - 'ESA': the pure implementation as described by ESA
         - 'pyroSAR': the ESA method plus the custom pyroSAR refinement
    
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
    
    tmp_base = os.path.basename(tmpdir)
    tmpdir_bnr = os.path.join(tmpdir, tmp_base + '_bnr')
    tmpdir_sub = os.path.join(tmpdir, tmp_base + '_sub')
    
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
        while i < len(groups) - 1:
            if 'Remove-GRD-Border-Noise' in groups[i]:
                del groups[i][groups[i].index('Remove-GRD-Border-Noise')]
            if len(groups[i]) == 0:
                del groups[i]
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
    
    if src_format == 'BEAM-DIMAP':
        src = src.replace('.dim', '.data')
    
    src_base = os.path.splitext(os.path.basename(src))[0]
    outname_base = os.path.join(outdir, src_base)
    
    if src_format in ['ENVI', 'BEAM-DIMAP']:
        log.info('converting to GeoTIFF')
        translateoptions = {'options': ['-q', '-co', 'INTERLEAVE=BAND', '-co', 'TILED=YES'],
                            'format': 'GTiff'}
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
                name_new = outname_base.replace(suffix, '{0}.tif'.format(base))
            if re.search('elevation', basename):
                nodata = dem_nodata
            elif re.search('layoverShadowMask', basename):
                nodata = 255
            else:
                nodata = 0
            translateoptions['noData'] = nodata
            if clean_edges and not 'layover_shadow_mask' in basename:
                erode_edges(item, only_boundary=True, pixels=clean_edges_npixels)
            gdal_translate(item, name_new, translateoptions)
    else:
        raise RuntimeError('The output file format must be ENVI or BEAM-DIMAP.')
    ###########################################################################
    # write the Sentinel-1 manifest.safe file as addition to the actual product
    readers = workflow['operator=Read']
    for reader in readers:
        infile = reader.parameters['file']
        try:
            id = identify(infile)
            if id.sensor in ['S1A', 'S1B']:
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
    split SNAP workflow into groups containing a maximum defined number of operators
    
    Parameters
    ----------
    xmlfile: str
        the SNAP xml workflow
    n: int
        the maximum number of worker nodes in each group; Read and Write are excluded

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
                    newgroup.insert(newgroup.index(worker), item)
            for writer in writers_id:
                if workflow[writer].source == worker:
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
        insert a node into the workflow including setting its source to its predecessor
        and setting its ID as source of the successor.
        
        Parameters
        ----------
        node: Node
            the node to be inserted
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
        Node or None
            the new node or None, depending on arguement `void`
        """
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
            return Par_BandMath(params)
        else:
            return Par(params)
    
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
    element: ~xml.etree.ElementTree.Element
        the node parameter XML element
    """
    
    def __init__(self, element):
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
            raise KeyError('key {} does not exist'.format(key))
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
    def __init__(self, element):
        self.__element = element
        super(Par_BandMath, self).__init__(element)
    
    def __getitem__(self, item):
        if item in ['variables', 'targetBands']:
            out = []
            for x in self.__element.findall('.//{}'.format(item[:-1])):
                out.append(Par(x))
            return out
        else:
            raise ValueError("can only get items 'variables' and 'targetBands'")
    
    def clear_variables(self):
        var = self.__element.find('.//variables')
        for item in var:
            var.remove(item)
    
    def add_equation(self):
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


def erode_edges(infile, only_boundary=False, connectedness=4, pixels=1):
    """
    Erode noisy edge pixels in SNAP-processed images.
    It was discovered that images contain border pixel artifacts after `Terrain-Correction`.
    Likely this is coming from treating the value 0 as regular value instead of no data during resampling.
    This function erodes these edge pixels using :func:`scipy.ndimage.binary_erosion`.
    scipy is not a base dependency of pyroSAR and has to be installed separately.
    
    .. figure:: figures/snap_erode_edges.png
        :align: center
        
        VV gamma0 RTC backscatter image visualizing the noisy border (left) and the cleaned result (right).
        The area covers approx. 2.3 * 2.3 km². Pixel spacing is 20 m. connectedness 4, 1 pixel.
    
    Parameters
    ----------
    infile: str
        a single-layer file to modify in-place. 0 is assumed as no data value.
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
    
    with Raster(infile) as ref:
        array = ref.array()
        mask = array != 0
        if only_boundary:
            with vectorize(target=mask, reference=ref) as vec:
                with boundary(vec) as bounds:
                    with rasterize(vectorobject=bounds, reference=ref, nodata=None) as new:
                        mask = new.array()
        mask2 = binary_erosion(input=mask, structure=structure)
        array[mask2 == 0] = 0
    
    ras = gdal.Open(infile, GA_Update)
    band = ras.GetRasterBand(1)
    band.WriteArray(array)
    band.FlushCache()
    band = None
    ras = None
