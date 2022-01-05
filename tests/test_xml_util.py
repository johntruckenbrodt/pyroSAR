
import os
import pytest
from pyroSAR import SAFE
from pyroSAR.xml_util import XMLHandler


def test_handler(tmpdir, testdata):
    id = SAFE(testdata['s1'])
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
