import os
import pytest
from pyroSAR import identify
from pyroSAR.snap import geocode
from spatialist import bbox
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


class Test_geocode_opts():
    def test_pol(self, tmpdir, testdata):
        scene = testdata['s1']
        with pytest.raises(RuntimeError):
            geocode(scene, str(tmpdir), polarizations=1, test=True)
        with pytest.raises(RuntimeError):
            geocode(scene, str(tmpdir), polarizations='foobar', test=True)
        geocode(scene, str(tmpdir), polarizations='VV', test=True)
    
    def test_geotype(self, tmpdir, testdata):
        scene = testdata['s1']
        with pytest.raises(RuntimeError):
            geocode(scene, str(tmpdir), geocoding_type='foobar', test=True)
        geocode(scene, str(tmpdir), test=True,
                geocoding_type='SAR simulation cross correlation')
    
    def test_srs(self, tmpdir, testdata):
        scene = testdata['s1']
        with pytest.raises(RuntimeError):
            geocode(scene, str(tmpdir), t_srs='foobar', test=True)
        geocode(scene, str(tmpdir), t_srs=32632, test=True)
    
    def test_scaling(self, tmpdir, testdata):
        scene = testdata['s1']
        with pytest.raises(RuntimeError):
            geocode(scene, str(tmpdir), scaling='foobar', test=True)
    
    def test_shp(self, tmpdir, testdata):
        scene = testdata['s1']
        with identify(scene).bbox() as box:
            ext = box.extent
        ext['xmax'] -= 1
        with bbox(ext, 4326) as new:
            geocode(scene, str(tmpdir), shapefile=new, test=True)
    
    def test_export_extra(self, tmpdir, testdata):
        scene = testdata['s1']
        with pytest.raises(AttributeError):
            geocode(scene, str(tmpdir), test=True,
                    export_extra=['foobar'])
        geocode(scene, str(tmpdir), test=True,
                export_extra=['localIncidenceAngle'])
    
    def test_externalDEM(self, tmpdir, testdata):
        scene = testdata['s1']
        with pytest.raises(RuntimeError):
            geocode(scene, str(tmpdir), externalDEMFile='foobar', test=True)
