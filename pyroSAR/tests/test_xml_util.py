
import os
import pytest
from pyroSAR import SAFE
from pyroSAR.xml_util import XMLHandler

testdata = 'pyroSAR/tests/data/'
testfile1 = os.path.join(testdata, 'S1A_IW_GRDH_1SDV_20150222T170750_20150222T170815_004739_005DD8_3768.zip')


def test_handler(tmpdir):
    id = SAFE(testfile1)
    id.unpack(str(tmpdir))
    testfile = os.path.join(id.scene, 'manifest.safe')
    xml = XMLHandler(testfile)
    xml.restoreNamespaces()
    xml.write(os.path.join(str(tmpdir), 'test.xml'), 'w')
    with pytest.raises(RuntimeError):
        xml = XMLHandler(1)
    with pytest.raises(RuntimeError):
        xml = XMLHandler('foobar')
    with open(testfile, 'r') as infile:
        xml = XMLHandler(infile)
