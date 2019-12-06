#####################################################################
# Module for testing the functionality of the SNAP processing module
#####################################################################
import os
import pytest
from pyroSAR import identify
from pyroSAR.snap import geocode
from spatialist import bbox
from spatialist.ancillary import finder
from pyroSAR.snap.auxil import is_consistent, split, groupbyWorkers, Workflow, parse_recipe
from pyroSAR.examine import ExamineSnap


def test_installation():
    reg = ExamineSnap()
    assert os.path.isfile(reg.gpt)


def test_consistency():
    with parse_recipe('base') as wf:
        assert is_consistent(wf)


def test_geocode(tmpdir, testdata):
    scene = testdata['s1']
    geocode(scene, str(tmpdir), test=True)
    xmlfile = finder(str(tmpdir), ['*.xml'])[0]
    tree = Workflow(xmlfile)
    assert is_consistent(tree) is True
    groups = groupbyWorkers(xmlfile, 2)
    assert len(groups) == 4
    groups2 = groupbyWorkers(xmlfile, 100)
    assert len(groups2) == 1
    split(xmlfile, groups)
    id = identify(scene)
    basename = '{}_{}'.format(id.outname_base(), tree.suffix)
    procdir = os.path.join(str(tmpdir), basename)
    assert os.path.isdir(procdir)
    tempdir = os.path.join(procdir, 'temp')
    assert os.path.isdir(tempdir)
    parts = finder(tempdir, ['*.xml'])
    assert len(parts) == 4


class Test_geocode_opts():
    def test_infile_type(self, tmpdir, testdata):
        scene = testdata['s1']
        with pytest.raises(TypeError):
            geocode(infile=123, outdir=str(tmpdir), test=True)
        id = identify(scene)
        geocode(infile=id, outdir=str(tmpdir), test=True)
    
    def test_pol(self, tmpdir, testdata):
        scene = testdata['s1']
        with pytest.raises(RuntimeError):
            geocode(scene, str(tmpdir), polarizations=1, test=True)
        with pytest.raises(RuntimeError):
            geocode(scene, str(tmpdir), polarizations='foobar', test=True)
        geocode(scene, str(tmpdir), polarizations='VV', test=True)
    
    def test_pol_list(self, tmpdir, testdata):
        scene = testdata['s1']
        geocode(scene, str(tmpdir), polarizations=['VV', 'VH'], test=True)
    
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
        ext = {'xmin': 12, 'xmax': 13, 'ymin': 53, 'ymax': 54}
        with bbox(ext, 4326) as new:
            with pytest.raises(RuntimeError):
                geocode(scene, str(tmpdir), shapefile=new, test=True)
        
        with identify(scene).bbox() as box:
            ext = box.extent
        ext['xmax'] -= 1
        with bbox(ext, 4326) as new:
            geocode(scene, str(tmpdir), shapefile=new, test=True)
    
    def test_offset(self, tmpdir, testdata):
        scene = testdata['s1']
        geocode(scene, str(tmpdir), offset=(100, 100, 0, 0), test=True)
    
    def test_export_extra(self, tmpdir, testdata):
        scene = testdata['s1']
        with pytest.raises(RuntimeError):
            geocode(scene, str(tmpdir), test=True,
                    export_extra=['foobar'])
        geocode(scene, str(tmpdir), test=True,
                export_extra=['localIncidenceAngle'])
    
    def test_externalDEM(self, tmpdir, testdata):
        scene = testdata['s1']
        dem_dummy = testdata['tif']
        with pytest.raises(RuntimeError):
            geocode(scene, str(tmpdir), externalDEMFile='foobar', test=True)
        geocode(scene, str(tmpdir), externalDEMFile=dem_dummy, test=True)

    def test_speckleFilter(self, tmpdir, testdata):
        scene = testdata['s1']
        with pytest.raises(ValueError):
            geocode(scene, str(tmpdir), speckleFilter='foobar', test=True)
        geocode(scene, str(tmpdir), speckleFilter='Refined Lee', test=True)
    
    def test_refarea(self, tmpdir, testdata):
        scene = testdata['s1']
        with pytest.raises(RuntimeError):
            geocode(scene, str(tmpdir), terrainFlattening=True, refarea='beta0', test=True)
        with pytest.raises(ValueError):
            geocode(scene, str(tmpdir), terrainFlattening=False, refarea='foobar', test=True)
        geocode(scene, str(tmpdir), terrainFlattening=True, refarea='gamma0', test=True)
    
    def test_sliceassembly(self, tmpdir, testdata):
        scene1 = testdata['s1']
        scene2 = testdata['s1_2']
        wf = geocode([scene1, scene2], str(tmpdir), test=True, returnWF=True)
        for n in range(1, 4):
            groups = groupbyWorkers(wf, n=n)
            split(wf, groups)
