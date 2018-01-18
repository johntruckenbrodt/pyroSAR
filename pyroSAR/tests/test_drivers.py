import pyroSAR
#import logging
#import unittest
import pytest
import os

testdir = os.getenv("TESTDATA_DIR", "pyroSAR/tests/data/")

testcases = [
    {"path": os.path.join("pyroSAR/tests/data", "S1A_IW_GRDH_1SDV_20150222T170750_20150222T170815_004739_005DD8_3768.zip"),
     "compression": "zip",
     "sensor": "S1A",
     "product": "GRD",
     "outname": 'S1A__IW___A_20150222T170750',
     "orbit": "A"},

    {"path": os.path.join(testdir, "0000022708_001001_ALOS2015976960-140909.zip"),
     "compression": "zip",
     "sensor": "PSR2",
     "product": "1.5",
     "outname": "PSR2_FBD__A_20140909T043342",
     "orbit": "A"}
]

@pytest.fixture
def scene(case):
    case['pyro'] = pyroSAR.identify(case['path'])
    return case

@pytest.mark.parametrize("case", testcases)
class Test_Metadata():
    def test_compression(self, scene):
        #scene = pyroSAR.identify(case["path"])
        assert scene['pyro'].compression == scene["compression"]

    def test_sensor(self, scene):
        assert scene['pyro'].sensor == scene["sensor"]

    def test_product(self, scene):
        assert scene['pyro'].product == scene["product"]

    def test_is_processed(self, scene):
        assert scene['pyro'].is_processed("data/") == False

    def test_outname(self, scene):
        assert scene['pyro'].outname_base() == scene["outname"]

    def test_orbit(self, scene):
        assert scene['pyro'].orbit == scene["orbit"]



def test_identify_fail():
    with pytest.raises(IOError):
        pyroSAR.identify(os.path.join(testdir, 'foobar'))

def test_export2dict():
    pass


"""
class TestMetadataS1(unittest.TestCase):
    def setUp(self):
        self.s1 = pyroSAR.identify("data/S1A_IW_GRDH_1SDV_20150222T170750_20150222T170815_004739_005DD8_3768.zip")
        #print self.s1.meta
    def tearDown(self):
        self.s1 = None

    def test_compression_zip(self):
        self.assertEqual(self.s1.compression, 'zip')

    def test_sensor_S1(self):
        self.assertEqual(self.s1.sensor, "S1A")

    def test_product(self):
        self.assertEqual(self.s1.product, "GRD")

    def test_getCorners(self):
        pass
        #self.assertEqual(self.s1.getCorners(), {'xmin':}) 

    def test_is_processed_False(self):
        self.assertFalse(self.s1.is_processed('data/'))

    def test_outname_base(self):
        self.assertEqual(self.s1.outname_base(), )

    def test_orbit(self):
        self.assertEqual(self.s1.orbit, 'A')

    #def test_

if __name__ == "__main__":
    unittest.main()
"""