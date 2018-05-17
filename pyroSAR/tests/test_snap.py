
import os
from pyroSAR.snap import geocode


def test_geocode():
    scene = 'pyroSAR/tests/data/S1A_IW_GRDH_1SDV_20150222T170750_20150222T170815_004739_005DD8_3768.zip'
    geocode(scene, 'pyroSAR/tests/data', test=True)
    os.remove('pyroSAR/tests/data/S1A__IW___A_20150222T170750_bnr_Orb_Cal_TF_TC_dB_proc.xml')
