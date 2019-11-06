import os
import re
import shutil
import subprocess as sp
from xml.dom import minidom
import xml.etree.ElementTree as ET

from osgeo import gdal
from osgeo.gdalconst import GA_Update

from pyroSAR import identify
from pyroSAR._dev_config import LOOKUP
from pyroSAR.examine import ExamineSnap
from pyroSAR.ancillary import parse_datasetname

from spatialist.auxil import gdal_translate
from spatialist.ancillary import finder


def parse_recipe(name):
    """
    parse a SNAP recipe
    
    Parameters
    ----------
    name: str
        the name of the recipe; current options:
         * `blank`: a workflow without any nodes
         * `geocode`: a basic workflow containing Read, Apply-Orbit-File, Calibration, Terrain-Flattening and Write nodes

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


def parse_node(name):
    """
    parse a XML node recipe
    
    Parameters
    ----------
    name: str
        the name of the processing node, e.g. Terrain-Correction

    Returns
    -------
    Node
        the parsed node
    
    Examples
    --------
    >>> ml = parse_node('ThermalNoiseRemoval')
    >>> print(ml.parameters)
    {'selectedPolarisations': None, 'removeThermalNoise': 'true', 'reIntroduceThermalNoise': 'false'}
    """
    name = name if name.endswith('.xml') else name + '.xml'
    absname = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'recipes', 'nodes', name)
    with open(absname, 'r') as workflow:
        element = ET.fromstring(workflow.read())
    return Node(element)


def execute(xmlfile, cleanup=True, gpt_exceptions=None, verbose=True):
    """
    execute SNAP workflows via the Graph Processing Tool gpt.
    This function merely calls gpt with some additional command
    line arguments and raises a RuntimeError on fail. This
    function is used internally by function :func:`gpt`, which
    should be used instead.
    
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
    verbose: bool
        print out status messages?
    
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
    infile = workflow['Read'].parameters['file']
    nodes = workflow.nodes()
    workers = [x.id for x in nodes if x.operator not in ['Read', 'Write']]
    message = ' -> '.join(workers)
    gpt_exec = None
    if gpt_exceptions is not None:
        for item, exec in gpt_exceptions.items():
            if item in workers:
                gpt_exec = exec
                message += ' (using {})'.format(exec)
                break
    if verbose:
        print(message)
    # try to find the GPT executable
    if gpt_exec is None:
        try:
            gpt_exec = ExamineSnap().gpt
        except AttributeError:
            raise RuntimeError('could not find SNAP GPT executable')
    # create the list of arguments to be passed to the subprocess module calling GPT
    if format == 'GeoTiff-BigTIFF':
        cmd = [gpt_exec,
               # '-Dsnap.dataio.reader.tileWidth=*',
               # '-Dsnap.dataio.reader.tileHeight=1',
               '-Dsnap.dataio.bigtiff.tiling.width=256',
               '-Dsnap.dataio.bigtiff.tiling.height=256',
               # '-Dsnap.dataio.bigtiff.compression.type=LZW',
               # '-Dsnap.dataio.bigtiff.compression.quality=0.75',
               xmlfile]
    else:
        cmd = [gpt_exec, xmlfile]
    # execute the workflow
    proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = proc.communicate()
    out = out.decode('utf-8') if isinstance(out, bytes) else out
    err = err.decode('utf-8') if isinstance(err, bytes) else err
    # delete intermediate files if an error occurred
    if proc.returncode == 1:
        pattern = r"Error: \[NodeId: (?P<id>[a-zA-Z0-9-_]*)\] " \
                  r"Operator \'[a-zA-Z0-9-_]*\': " \
                  r"Unknown element \'(?P<par>[a-zA-Z]*)\'"
        match = re.search(pattern, err)
        if match is not None:
            replace = match.groupdict()
            with Workflow(xmlfile) as flow:
                print('  removing parameter {id}:{par} and executing modified workflow'.format(**replace))
                node = flow[replace['id']]
                del node.parameters[replace['par']]
                flow.write(xmlfile)
            execute(xmlfile, cleanup=cleanup, gpt_exceptions=gpt_exceptions, verbose=verbose)
        else:
            if cleanup:
                if os.path.isfile(outname + '.tif'):
                    os.remove(outname + '.tif')
                elif os.path.isdir(outname):
                    shutil.rmtree(outname)
            print(out + err)
            print('failed: {}'.format(os.path.basename(infile)))
            err_match = re.search('Error: (.*)\n', out + err)
            errmessage = err_match.group(1) if err_match else err
            raise RuntimeError(errmessage)
    elif proc.returncode == -9:
        if cleanup:
            if os.path.isfile(outname + '.tif'):
                os.remove(outname + '.tif')
            elif os.path.isdir(outname):
                shutil.rmtree(outname)
        print('the process was killed by SNAP. One possible cause is a lack of memory.')
    else:
        print('process return code: {}'.format(proc.returncode))


def gpt(xmlfile, groups=None, cleanup=True, gpt_exceptions=None,
        removeS1BorderNoiseMethod='pyroSAR', basename_extensions=None):
    """
    wrapper for ESA SNAP's Graph Processing Tool GPT.
    Input is a readily formatted workflow XML file as
    created by function :func:`~pyroSAR.snap.util.geocode`.
    Additional to calling GPT, this function will
    
     * execute the workflow in groups as defined by `groups`
     * encode a nodata value into the output file if the format is GeoTiff-BigTIFF
     * convert output files to GeoTiff if the output format is ENVI
    
    Parameters
    ----------
    xmlfile: str
        the name of the workflow XML file
    groups: list
        a list of lists each containing IDs for individual nodes
    cleanup: bool
        should all files written to the temporary directory during function execution be deleted after processing?
    gpt_exceptions: dict
        a dictionary to override the configured GPT executable for certain operators;
        each (sub-)workflow containing this operator will be executed with the define executable;
        
         - e.g. ``{'Terrain-Flattening': '/home/user/snap/bin/gpt'}``
    removeS1BorderNoiseMethod: str
        the border noise removal method to be applied, See :func:`pyroSAR.S1.removeGRDBorderNoise` for details; one of the following:
         - 'ESA': the pure implementation as described by ESA
         - 'pyroSAR': the ESA method plus the custom pyroSAR refinement
    basename_extensions: list of str
        names of additional parameters to append to the basename, e.g. ['orbitNumber_rel']
    
    Returns
    -------

    """
    
    workflow = Workflow(xmlfile)
    write = workflow['Write']
    outname = write.parameters['file']
    suffix = workflow.suffix
    format = write.parameters['formatName']
    dem_name = workflow.tree.find('.//demName').text
    if dem_name == 'External DEM':
        dem_nodata = float(workflow.tree.find('.//externalDEMNoDataValue').text)
    else:
        dem_nodata = 0
    
    if 'Remove-GRD-Border-Noise' in workflow.ids and removeS1BorderNoiseMethod == 'pyroSAR':
        if 'SliceAssembly' in workflow.operators:
            raise RuntimeError("pyroSAR's custom border noise removal is not yet implemented for multiple scene inputs")
        xmlfile = os.path.join(outname,
                               os.path.basename(xmlfile.replace('_bnr', '')))
        if not os.path.isdir(outname):
            os.makedirs(outname)
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
        # identify the input scene, unpack it and perform the custom border noise removal
        read = workflow['Read']
        scene = identify(read.parameters['file'])
        print('unpacking scene')
        scene.unpack(outname)
        print('removing border noise..')
        scene.removeGRDBorderNoise(method=removeS1BorderNoiseMethod)
        # change the name of the input file to that of the unpacked archive
        read.parameters['file'] = scene.scene
        # write a new workflow file
        workflow.write(xmlfile)
    
    print('executing node sequence{}..'.format('s' if groups is not None else ''))
    if groups is not None:
        subs = split(xmlfile, groups)
        for sub in subs:
            execute(sub, cleanup=cleanup, gpt_exceptions=gpt_exceptions)
    else:
        execute(xmlfile, cleanup=cleanup, gpt_exceptions=gpt_exceptions)
    
    if format == 'ENVI':
        print('converting to GTiff')
        translateoptions = {'options': ['-q', '-co', 'INTERLEAVE=BAND', '-co', 'TILED=YES'],
                            'format': 'GTiff'}
        for item in finder(outname, ['*.img'], recursive=False):
            if re.search('[HV]{2}', item):
                pol = re.search('[HV]{2}', item).group()
                name_new = outname.replace(suffix, '{0}_{1}.tif'.format(pol, suffix))
            else:
                base = os.path.splitext(os.path.basename(item))[0] \
                    .replace('elevation', 'DEM')
                name_new = outname.replace(suffix, '{0}.tif'.format(base))
            nodata = dem_nodata if re.search('elevation', item) else 0
            translateoptions['noData'] = nodata
            gdal_translate(item, name_new, translateoptions)
    # by default the nodata value is not registered in the GTiff metadata
    elif format == 'GeoTiff-BigTIFF':
        ras = gdal.Open(outname + '.tif', GA_Update)
        for i in range(1, ras.RasterCount + 1):
            ras.GetRasterBand(i).SetNoDataValue(0)
        ras = None
    if cleanup:
        shutil.rmtree(outname)
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
                basename = '{0}_{1}_manifest.safe'.format(basename, suffix)
                outdir = os.path.dirname(outname)
                outname_manifest = os.path.join(outdir, basename)
                with open(outname_manifest, 'wb') as out:
                    out.write(manifest.read())
        except RuntimeError:
            continue
    ###########################################################################
    print('done')


def is_consistent(nodes):
    """
    check whether all nodes take either no source node or one that is in the list
    
    Parameters
    ----------
    nodes: list of Node
        a group of nodes
    Returns
    -------
    bool
        is the list of nodes consistent?
    """
    ids = [x.id for x in nodes]
    check = []
    for node in nodes:
        source = node.source
        if source is None or source in ids or all([x in ids for x in source]):
            check.append(True)
        else:
            check.append(False)
    return all(check)


def split(xmlfile, groups):
    """
    split a workflow file into groups and write them to separate workflows including source and write target linking.
    The new workflows are written to a sub-directory 'temp' of the target directory defined in the input's 'Write' node.
    Each new workflow is parametrized with a Read and Write node if they don't already exist. Temporary outputs are
    written to BEAM-DIMAP files named after the workflow suffix sequence.
    
    Parameters
    ----------
    xmlfile: str
        the workflow to be split
    groups: list
        a list of lists each containing IDs for individual nodes

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
    out = write.parameters['file']
    tmp = os.path.join(out, 'temp')
    if not os.path.isdir(tmp):
        os.makedirs(tmp)
    
    outlist = []
    prod_tmp = []
    prod_tmp_format = []
    for position, group in enumerate(groups):
        new = parse_recipe('blank')
        nodes = [workflow[x] for x in group]
        if nodes[0].operator != 'Read':
            read = new.insert_node(parse_node('Read'), void=False)
        else:
            read = new.insert_node(nodes[0], void=False)
            del nodes[0]
        
        if position != 0:
            read.parameters['file'] = prod_tmp[-1]
            read.parameters['formatName'] = prod_tmp_format[-1]
        
        for i, node in enumerate(nodes):
            if i == 0:
                new.insert_node(node, before=read.id)
            else:
                new.insert_node(node)
        
        nodes = new.nodes()
        operators = [node.operator for node in nodes]
        writers = new['operator=Write']
        formats = [write.parameters['formatName'] for write in writers]
        if (position < len(groups) - 1 and 'BEAM-DIMAP' not in formats) \
                or (operators[-1] != 'Write'):
            write = parse_node('Write')
            new.insert_node(write, before=nodes[-1].id)
            tmp_out = os.path.join(tmp, 'tmp{}.dim'.format(position))
            prod_tmp.append(tmp_out)
            prod_tmp_format.append('BEAM-DIMAP')
            write.parameters['file'] = tmp_out
            write.parameters['formatName'] = 'BEAM-DIMAP'
            operators.append('Write')
        else:
            prod_tmp.append(nodes[-1].parameters['file'])
            prod_tmp_format.append(nodes[-1].parameters['formatName'])
        nodes = new.nodes()
        if not is_consistent(nodes):
            message = 'inconsistent group:\n {}'.format(' -> '.join(group))
            raise RuntimeError(message)
        outname = os.path.join(tmp, 'tmp{}.xml'.format(position))
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
    list
        a list of lists each containing the IDs of all nodes belonging to the groups including Read and Write nodes;
        this list can e.g. be passed to function :func:`split` to split the workflow into new sub-workflow files based
        on the newly created groups
    """
    workflow = Workflow(xmlfile)
    nodes = workflow.nodes()
    workers_id = [x.id for x in nodes if x.operator not in ['Read', 'Write']]
    readers_id = [x.id for x in workflow['operator=Read']]
    writers_id = [x.id for x in workflow['operator=Write']]
    workers_groups = [workers_id[i:i + n] for i in range(0, len(workers_id), n)]
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
        elif re.search(pattern, item):
            key, value = re.search(pattern, item).groups()
            return [x for x in self.nodes() if getattr(x, key) == value]
        else:
            return Node(self.tree.find('.//node[@id="{}"]'.format(item)))
    
    def __delitem__(self, key):
        if not isinstance(key, str):
            raise TypeError('key must be of type str')
        element = self.tree.find('.//node[@id="{}"]'.format(key))
        node = Node(element)
        source = node.source
        successors = [x for x in self.nodes() if x.source == key]
        for node in successors:
            node.source = source
        self.tree.remove(element)
    
    def __str__(self):
        self.__optimize_appearance()
        rough_string = ET.tostring(self.tree, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        return reparsed.toprettyxml(indent='\t', newl='')
    
    def __find_source(self, predecessor):
        if isinstance(predecessor, list):
            return [self.__find_source(x) for x in predecessor]
        else:
            source_id = predecessor
            while True:
                # a Write node cannot be the source of another node
                if self[source_id].operator == 'Write':
                    source_id = self[source_id].source
                else:
                    break
            return source_id
    
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
        return [node.id for node in self.nodes()]
    
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
        before: str or list
            the ID(s) of the node(s) before the newly inserted node; a list of node IDs is intended for nodes that
            require multiple sources, e.g. sliceAssembly
        after: str
            the ID of the node after the newly inserted node
        resetSuccessorSource: bool
            reset the source of the successor node to the ID of the newly inserted node?
        void: bool
            if false, the function returns the node

        Returns
        -------
        Node or None
            the new node or None, depending on arguement `void`
        """
        if before and not after:
            if isinstance(before, list):
                indices = [self.index(self[x]) for x in before]
                predecessor = self[before[indices.index(max(indices))]]
            else:
                predecessor = self[before]
            # print('inserting node {} after {}'.format(node.id, predecessor.id))
            position = self.index(predecessor) + 1
            self.tree.insert(position, node.element)
            newnode = Node(self.tree[position])
            # print('new node id: {}'.format(newnode.id))
            ####################################################
            # set the source product for the new node
            if newnode.operator != 'Read':
                newnode.source = self.__find_source(before)
            ####################################################
            # set the source product for the node after the new node
            if resetSuccessorSource:
                try:
                    successor = self[position]
                    # print('resetting source of successor {} to {}'.format(successor.id, newnode.id))
                    successor.source = newnode.id
                except IndexError:
                    # case where no successor exists because the new node is the new last node in the graph
                    pass
            # else:
            #     print('no source resetting required')
        elif after and not before:
            successor = self[after]
            # print('inserting node {} before {}'.format(node.id, successor.id))
            position = self.index(successor)
            self.tree.insert(position, node.element)
            newnode = Node(self.tree[position])
            # print('new node id: {}'.format(newnode.id))
            ####################################################
            # set the source product for the new node
            try:
                if newnode.operator != 'Read':
                    newnode.source = self.__find_source(self[position - 2])
            except IndexError:
                newnode.source = None
            ####################################################
            # set the source product for the node after the new node
            if resetSuccessorSource:
                # print('resetting source of successor {} to {}'.format(successor.id, newnode.id))
                successor.source = newnode.id
            # else:
            # print('no source resetting required')
        else:
            self.tree.insert(len(self.tree) - 1, node.element)
        self.refresh_ids()
        # print([x.id for x in self.nodes()])
        if not void:
            return node
    
    def nodes(self):
        """
        
        Returns
        -------
        list
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
        return sorted(list(set([node.operator for node in self.nodes()])))
    
    def refresh_ids(self):
        counter = {}
        for node in self.nodes():
            operator = node.operator
            if operator not in counter.keys():
                counter[operator] = 1
                node.id = operator
            else:
                counter[operator] += 1
                node.id = '{} ({})'.format(operator, counter[operator])
    
    def set_par(self, key, value):
        """
        set a parameter for all nodes in the workflow
        
        Parameters
        ----------
        key: str
            the parameter name
        value

        Returns
        -------

        """
        for node in self.nodes():
            if key in node.parameters.keys():
                node.parameters[key] = value2str(value)
    
    @property
    def suffix(self):
        """
        
        Returns
        -------
        str
            a file suffix created from the order of which the nodes will be executed
        """
        nodes = self.tree.findall('node')
        names = [re.sub(r'[ ]*\([0-9]+\)', '', y.attrib['id']) for y in nodes]
        suffix = '_'.join(filter(None, [LOOKUP.snap.suffix[x] for x in names]))
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
        if source is not None:
            source.attrib['refid'] = value
        else:
            raise RuntimeError('cannot set source on {} node'.format(self.operator))
    
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
        Par
            the processing parameters of the node
        """
        params = self.element.find('.//parameters')
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
        if isinstance(value, str):
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
