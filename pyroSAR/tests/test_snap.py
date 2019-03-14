import os
from pyroSAR import identify
from pyroSAR.snap import geocode
from spatialist.ancillary import finder
import xml.etree.ElementTree as ET
from pyroSAR.snap.auxil import is_consistent, split, groupbyWorkers, parse_suffix, ExamineSnap


def test_installation():
    reg = ExamineSnap()
    assert os.path.isfile(reg.gpt)


def test_geocode(tmpdir, testdata):
    scene = testdata['s1']
    geocode(scene, str(tmpdir), test=True)
    xmlfile = finder(str(tmpdir), ['*.xml'])[0]
    with open(xmlfile, 'r') as infile:
        tree = ET.fromstring(infile.read())
    nodes = tree.findall('node')
    assert is_consistent(nodes) is True
    groups = groupbyWorkers(xmlfile, 2)
    split(xmlfile, groups)
    id = identify(scene)
    basename = '{}_{}'.format(id.outname_base(), parse_suffix(tree))
    procdir = os.path.join(str(tmpdir), basename)
    assert os.path.isdir(procdir)
    tempdir = os.path.join(procdir, 'temp')
    assert os.path.isdir(tempdir)
    parts = finder(tempdir, ['*.xml'])
    assert len(parts) == 4
