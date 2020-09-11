import os
import pytest


@pytest.fixture
def travis():
    return 'TRAVIS' in os.environ.keys()


@pytest.fixture
def appveyor():
    return 'APPVEYOR' in os.environ.keys()


@pytest.fixture
def testdir():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')


@pytest.fixture
def testdata(testdir):
    out = {
        's1': os.path.join(testdir, 'S1A_IW_GRDH_1SDV_20150222T170750_20150222T170815_004739_005DD8_3768.zip'),
        's1_2': os.path.join(testdir, 'S1A_IW_GRDH_1SDV_20150222T170725_20150222T170750_004739_005DD8_CEAB.zip'),
        's1_3': os.path.join(testdir, 'S1A_IW_GRDH_1SDV_20150203T043109_20150203T043134_004454_00574F_6D00.zip'),
        's1_4': os.path.join(testdir, 'S1A_IW_GRDH_1SDV_20150203T043109_20150203T043134_004454_00574F_FEC3.zip'),
        # ftp://ftp.eorc.jaxa.jp/pub/ALOS-2/1501sample/310_forestbrazil/0000022708_001001_ALOS2015976960-140909.zip
        'psr2': os.path.join(testdir, '0000022708_001001_ALOS2015976960-140909.zip'),
        'tif': os.path.join(testdir, 'S1A__IW___A_20150309T173017_VV_grd_mli_geo_norm_db.tif'),
        'archive_old': os.path.join(testdir, 'archive_outdated.csv'),
        'dempar': os.path.join(testdir, 'dem.par'),
        'mlipar': os.path.join(testdir, 'mli.par')
    }
    return out


@pytest.fixture
def auxdata_dem_cases():
    cases = [('AW3D30', ['N050E010/N051E011.tar.gz']),
             ('SRTM 1Sec HGT', ['N51E011.SRTMGL1.hgt.zip']),
             ('SRTM 3Sec', ['srtm_39_02.zip']),
             ('TDX90m', ['90mdem/DEM/N51/E010/TDM1_DEM__30_N51E011.zip'])]
    return cases
