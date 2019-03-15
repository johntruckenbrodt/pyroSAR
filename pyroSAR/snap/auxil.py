import ast
import json
import warnings
import pkg_resources

import os
import re
import shutil
import subprocess as sp
from xml.dom import minidom
import xml.etree.ElementTree as ET

from osgeo import gdal
from osgeo.gdalconst import GA_Update

from pyroSAR import ConfigHandler
from .._dev_config import LOOKUP

from spatialist.auxil import gdal_translate
from spatialist.ancillary import finder


def parse_recipe(name):
    name = name if name.endswith('.xml') else name + '.xml'
    absname = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'recipes', name)
    with open(absname, 'r') as workflow:
        tree = ET.fromstring(workflow.read())
    return tree


def parse_node(name):
    name = name if name.endswith('.xml') else name + '.xml'
    absname = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'recipes', 'nodes', name)
    with open(absname, 'r') as workflow:
        tree = ET.fromstring(workflow.read())
    return tree


def parse_suffix(workflow):
    nodes = workflow.findall('node')
    suffix = '_'.join(filter(None, [LOOKUP.snap.suffix[x] for x in [y.attrib['id'] for y in nodes]]))
    return suffix


def insert_node(workflow, node, before=None, after=None, resetSuccessorSource=True, void=True):
    if before and not after:
        predecessor = workflow.find('.//node[@id="{}"]'.format(before))
        position = list(workflow).index(predecessor) + 1
        workflow.insert(position, node)
        newnode = workflow[position]
        # set the source product for the new node
        newnode.find('.//sources/sourceProduct').attrib['refid'] = predecessor.attrib['id']
        # set the source product for the node after the new node
        if resetSuccessorSource:
            successor = workflow[position + 1]
            if successor.tag == 'node':
                successor.find('.//sources/sourceProduct').attrib['refid'] = newnode.attrib['id']
    elif after and not before:
        successor = workflow.find('.//node[@id="{}"]'.format(after))
        position = list(workflow).index(successor)
        workflow.insert(position, node)
        newnode = workflow[position]
        # set the source product for the new node
        predecessor = workflow[position - 1]
        source_id = predecessor.attrib['id'] if predecessor.tag == 'node' else None
        newnode.find('.//sources/sourceProduct').attrib['refid'] = source_id
        # set the source product for the node after the new node
        if resetSuccessorSource:
            successor.find('.//sources/sourceProduct').attrib['refid'] = newnode.attrib['id']
    else:
        workflow.insert(len(workflow) - 1, node)
    if not void:
        return node


def write_recipe(recipe, outfile):
    outfile = outfile if outfile.endswith('.xml') else outfile + '.xml'
    rough_string = ET.tostring(recipe, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    with open(outfile, 'w') as out:
        out.write(reparsed.toprettyxml(indent='\t', newl=''))


def execute(xmlfile):
    """
    execute SNAP workflows via the Graph Processing Tool gpt
    
    Parameters
    ----------
    xmlfile: str
        the name of the workflow XML file

    Returns
    -------
    
    Raises
    ------
    RuntimeError
    """
    # read the file and extract some information
    with open(xmlfile, 'r') as infile:
        workflow = ET.fromstring(infile.read())
    write = workflow.find('.//node[@id="Write"]')
    outname = write.find('.//parameters/file').text
    infile = workflow.find('.//node[@id="Read"]/parameters/file').text
    nodes = workflow.findall('node')
    workers = [x.attrib['id'] for x in nodes if x.find('.//operator').text not in ['Read', 'Write']]
    print('_'.join(workers))
    # try to find the GPT executable
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
    # delete intermediate files if an error occured
    if proc.returncode != 0:
        if os.path.isfile(outname + '.tif'):
            os.remove(outname + '.tif')
        elif os.path.isdir(outname):
            shutil.rmtree(outname)
        print(out + err)
        print('failed: {}'.format(os.path.basename(infile)))
        err_match = re.search('Error: (.*)\n', out + err)
        errmessage = err_match.group(1) if err_match else err
        raise RuntimeError(errmessage)


def gpt(xmlfile, groups=None):
    """
    wrapper for ESA SNAP's Graph Processing Tool GPT.
    Input is a readily formatted workflow xml file as created by function :func:`snap.util.geocode`
    
    Parameters
    ----------
    xmlfile: str
        the name of the workflow XML file
    groups: list
        a list of lists each containing IDs for individual nodes

    Returns
    -------

    """
    
    with open(xmlfile, 'r') as infile:
        workflow = ET.fromstring(infile.read())
    write = workflow.find('.//node[@id="Write"]')
    outname = write.find('.//parameters/file').text
    format = write.find('.//parameters/formatName').text
    dem_name = workflow.find('.//demName').text
    if dem_name == 'External DEM':
        dem_nodata = float(workflow.find('.//externalDEMNoDataValue').text)
    else:
        dem_nodata = 0
    print('executing node sequence{}..'.format('s' if groups is not None else ''))
    if groups is not None:
        subs = split(xmlfile, groups)
        for sub in subs:
            execute(sub)
    else:
        execute(xmlfile)
    
    if format == 'ENVI':
        print('converting to GTiff')
        suffix = parse_suffix(workflow)
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
        shutil.rmtree(outname)
    # by default the nodata value is not registered in the GTiff metadata
    elif format == 'GeoTiff-BigTIFF':
        ras = gdal.Open(outname + '.tif', GA_Update)
        for i in range(1, ras.RasterCount + 1):
            ras.GetRasterBand(i).SetNoDataValue(0)
        ras = None
    print('done')


class ExamineSnap(object):
    """
    Class to check if ESA SNAP is installed.
    Upon initialization, this class searches for relevant binaries and the accompanying
    relative directory structure, which uniquely identify an ESA SNAP installation on a system.
    First, all relevant file and folder names are read from the pyroSAR config file if it exists
    and their existence is verified.
    If this fails, a system check is performed to find relevant binaries in the system PATH variable and
    additional files and folders relative to them.
    Furthermore, a snap.auxdata.properties file is scanned for auxiliary data URLs and local storage location.
    This is used by SNAP to manage data from e.g. the SRTM mission. In case SNAP is not installed, the respective
    information is read from a default file delivered with pyroSAR. This has the advantage of using the SNAP download
    URLs and local directory structure without having SNAP installed such that it can be adapted by other SAR software.
    """
    
    def __init__(self):
        
        # define some attributes which identify SNAP
        self.identifiers = ['path', 'gpt', 'etc', 'auxdata']
        
        # try reading all necessary attributes from the config file
        # print('reading config..')
        self.__read_config()
        
        # if SNAP could not be identified from the config attributes, do a system search for it
        if not self.__is_identified():
            # print('identifying SNAP..')
            self.__identify_snap()
        
        # if the auxdatapath attribute was not yet set, create a default directory
        if not hasattr(self, 'auxdatapath'):
            self.auxdatapath = os.path.join(os.path.expanduser('~'), '.snap', 'auxdata')
            if not os.path.isdir(self.auxdatapath):
                os.makedirs(self.auxdatapath)
        
        # if the SNAP auxdata properties attribute was not yet identified,
        # point it to the default file delivered with pyroSAR
        if not hasattr(self, 'properties'):
            # print('using default properties file..')
            template = 'data/snap.auxdata.properties'
            self.properties = pkg_resources.resource_filename(__name__, template)
        
        # update the snap properties; this reads the 'properties' file and looks for any changes,
        # which are then updated for the object
        self.__update_snap_properties()
        
        # update the config file: this scans for config changes and re-writes the config file if any are found
        self.__update_config()
    
    def __is_identified(self):
        """
        Check if SNAP has been properly identified, i.e. all paths in self.identifiers
        have been detected and confirmed.
        
        Returns
        -------
        bool
        """
        return sum([hasattr(self, x) for x in self.identifiers]) == len(self.identifiers)
    
    def __identify_snap(self):
        """
        do a comprehensive search for an ESA SNAP installation
        
        Returns
        -------
        bool
            has the SNAP properties file been changed?
        """
        # create a list of possible SNAP executables
        defaults = ['snap64.exe', 'snap32.exe', 'snap.exe', 'snap']
        paths = os.environ['PATH'].split(os.path.pathsep)
        options = [os.path.join(path, option) for path in paths for option in defaults]
        options = [x for x in options if os.path.isfile(x)]
        
        if not hasattr(self, 'path') or not os.path.isfile(self.path):
            executables = options
        else:
            executables = [self.path] + options
        
        # for each possible SNAP executable, check whether additional files and directories exist relative to it
        # to confirm whether it actually is a ESA SNAP installation or something else like e.g. the Ubuntu App Manager
        for path in executables:
            if os.path.islink(path):
                path = os.path.realpath(path)
            
            # check whether a directory etc exists relative to the SNAP executable
            etc = os.path.join(os.path.dirname(os.path.dirname(path)), 'etc')
            if not os.path.isdir(etc):
                continue
            
            # check the content of the etc directory
            auxdata = os.listdir(etc)
            if 'snap.auxdata.properties' not in auxdata:
                continue
            else:
                auxdata_properties = os.path.join(etc, 'snap.auxdata.properties')
            
            # identify the gpt executable
            gpt_candidates = finder(os.path.dirname(path), ['gpt', 'gpt.exe'])
            if len(gpt_candidates) == 0:
                continue
            else:
                gpt = gpt_candidates[0]
            
            self.path = path
            self.etc = etc
            self.gpt = gpt
            self.auxdata = auxdata
            self.properties = auxdata_properties
            return
        
        warnings.warn('SNAP could not be identified. If you have installed it please add the path to the SNAP '
                      'executables (bin subdirectory) to the PATH environment. '
                      'E.g. in the Linux .bashrc file add the following line:\nexport PATH=$PATH:path/to/snap/bin"')
    
    def __read_config(self):
        """
        This method reads the config.ini to examine the snap paths.
        If the snap paths are not in the config.ini or the paths are wrong they will be automatically created.

        Returns
        -------

        """
        for attr in self.identifiers + ['auxdatapath', 'properties']:
            self.__read_config_attr(attr, 'SNAP')
        
        snap_properties = {}
        if 'OUTPUT' in ConfigHandler.sections:
            snap_properties = ConfigHandler['OUTPUT']
        if len(snap_properties.keys()) > 0:
            setattr(self, 'snap_properties', snap_properties)
    
    def __read_config_attr(self, attr, section):
        """
        read an attribute from the config file and set it as an object attribute
        
        Parameters
        ----------
        attr: str
            the attribute name
        section: str
            the config section to read the attribute from
        
        Returns
        -------
        
        """
        if section in ConfigHandler.sections:
            if attr in ConfigHandler[section].keys():
                val = ConfigHandler[section][attr]
                if attr in ['path', 'gpt', 'properties']:
                    exist = os.path.isfile(val)
                elif attr == 'auxdata':
                    val = ast.literal_eval(val)
                    exist = isinstance(val, list)
                else:
                    exist = os.path.isdir(val)
                if exist:
                    # print('setting attribute {}'.format(attr))
                    setattr(self, attr, val)
    
    def __update_config(self):
        for section in ['SNAP', 'OUTPUT', 'URL']:
            if section not in ConfigHandler.sections:
                # print('creating section {}..'.format(section))
                ConfigHandler.add_section(section)
        
        for key in self.identifiers + ['auxdatapath', 'properties']:
            if hasattr(self, key):
                self.__update_config_attr(key, getattr(self, key), 'SNAP')
        
        for key, value in self.snap_properties.items():
            self.__update_config_attr(key, value, 'OUTPUT')
    
    @staticmethod
    def __update_config_attr(attr, value, section):
        if isinstance(value, list):
            value = json.dumps(value)
        
        if attr not in ConfigHandler[section].keys() or ConfigHandler[section][attr] != value:
            # print('updating attribute {0}:{1}..'.format(section, attr))
            # print('  {0} -> {1}'.format(repr(ConfigHandler[section][attr]), repr(value)))
            ConfigHandler.set(section, key=attr, value=value, overwrite=True)
    
    def __update_snap_properties(self):
        """
        Read the snap.auxdata.properties file entries to object attributes

        Returns
        -------

        """
        pattern = r'^(?P<key>[\w\.]*)\s*=\s*(?P<value>.*)\n'
        
        if not hasattr(self, 'snap_properties'):
            self.snap_properties = {}
        
        demPath = os.path.join(self.auxdatapath, 'dem')
        landCoverPath = os.path.join(self.auxdatapath, 'LandCover')
        
        with open(self.properties, 'r') as prop:
            for line in prop:
                if re.search(pattern, line):
                    key, value = re.match(re.compile(pattern), line).groups()
                    value = value \
                        .replace('${AuxDataPath}', self.auxdatapath) \
                        .replace('${demPath}', demPath) \
                        .replace('${landCoverPath}', landCoverPath) \
                        .replace('\\', '/')
                    if not key in self.snap_properties.keys() or self.snap_properties[key] != value:
                        self.snap_properties[key] = value


def is_consistent(nodes):
    """
    check whether all nodes take either no source node or on that is in the list
    
    Parameters
    ----------
    nodes: list of ET.Element
        a group of nodes
    Returns
    -------
    bool
        is the list of nodes consistent
    """
    ids = [x.attrib['id'] for x in nodes]
    check = []
    for node in nodes:
        source = node.find('.//sources/sourceProduct')
        if source is None or source.attrib['refid'] in ids:
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

    """
    with open(xmlfile, 'r') as infile:
        workflow = ET.fromstring(infile.read())
    write = workflow.find('.//node[@id="Write"]')
    out = write.find('.//parameters/file').text
    tmp = os.path.join(out, 'temp')
    if not os.path.isdir(tmp):
        os.makedirs(tmp)
    
    outlist = []
    prod_tmp = []
    prod_tmp_format = []
    for position, group in enumerate(groups):
        new = parse_recipe('blank')
        nodes = [workflow.find('.//node[@id="{}"]'.format(x)) for x in group]
        if nodes[0].find('.//operator').text != 'Read':
            read = insert_node(new, parse_node('Read'), void=False)
        else:
            read = insert_node(new, nodes[0], void=False)
            del nodes[0]
        
        if position != 0:
            read.find('.//parameters/file').text = prod_tmp[-1]
            read.find('.//parameters/formatName').text = prod_tmp_format[-1]
        
        for i, node in enumerate(nodes):
            if i == 0:
                insert_node(new, node, before=read.attrib['id'])
            else:
                insert_node(new, node)
        
        nodes = new.findall('node')
        operators = [node.find('.//operator').text for node in nodes]
        suffix = parse_suffix(new)
        if operators[-1] != 'Write':
            write = insert_node(new, parse_node('Write'), before=nodes[-1].attrib['id'], void=False)
            tmp_out = os.path.join(tmp, 'tmp{}.dim'.format(position))
            prod_tmp.append(tmp_out)
            prod_tmp_format.append('BEAM-DIMAP')
            write.find('.//parameters/file').text = tmp_out
            write.find('.//parameters/formatName').text = 'BEAM-DIMAP'
            operators.append('Write')
        else:
            prod_tmp.append(nodes[-1].find('.//parameters/file').text)
            prod_tmp_format.append(nodes[-1].find('.//parameters/formatName').text)
        
        if not is_consistent(nodes):
            raise RuntimeError('inconsistent group:\n {}'.format('-'.format(group)))
        outname = os.path.join(tmp, 'tmp{}.xml'.format(position))
        write_recipe(new, outname)
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
    with open(xmlfile, 'r') as infile:
        workflow = ET.fromstring(infile.read())
    nodes = workflow.findall('node')
    nodes_id = [x.attrib['id'] for x in nodes]
    workers = [x for x in nodes if x.find('.//operator').text not in ['Read', 'Write']]
    workers_id = [x.attrib['id'] for x in workers]
    workers_groups = [workers_id[i:i + n] for i in range(0, len(workers), n)]
    splits = [nodes_id.index(x[0]) for x in workers_groups] + [len(nodes_id)]
    splits[0] = 0
    nodes_id_split = [nodes_id[splits[x]:splits[x + 1]] for x in range(0, len(splits) - 1)]
    return nodes_id_split
