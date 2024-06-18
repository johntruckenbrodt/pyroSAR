import os
import pytest
import platform


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
        # ASAR_IMS__A_20040703T205338, product: SLC, driver: ESA
        'asar': os.path.join(testdir,
                             'ASA_IMS_1PNESA20040703_205338_000000182028_00172_12250_00001672562030318361237.N1'),
        # ERS1_IMP__A_19960808T205906, product: PRI, driver: ESA
        'ers1_esa': os.path.join(testdir, 'SAR_IMP_1PXESA19960808_205906_00000017G158_00458_26498_2615.E1'),
        # ERS1_IMS__D_19951220T024320, product: SLC, driver: CEOS_ERS
        'ers1_ceos': os.path.join(testdir, 'SAR_IMS_1PXESA19951220_024320_00000015G152_00132_23166_0252.E1.zip'),
        # PSR2_FBD__A_20140909T043342, product: 1.5, driver: CEOS_PSR
        'psr2': os.path.join(testdir, '0000022708_001001_ALOS2015976960-140909.zip'),
        # main scene for testing Sentinel-1 metadata reading and database ingestion
        's1': os.path.join(testdir, 'S1A_IW_GRDH_1SDV_20150222T170750_20150222T170815_004739_005DD8_3768.zip'),
        # for test_snap.test_slice_assembly
        's1_2': os.path.join(testdir, 'S1A_IW_GRDH_1SDV_20150222T170725_20150222T170750_004739_005DD8_CEAB.zip'),
        # for testing database duplicate handling
        's1_3': os.path.join(testdir, 'S1A_IW_GRDH_1SDV_20150203T043109_20150203T043134_004454_00574F_6D00.zip'),
        # for testing database duplicate handling
        's1_4': os.path.join(testdir, 'S1A_IW_GRDH_1SDV_20150203T043109_20150203T043134_004454_00574F_FEC3.zip'),
        # used in test_osv
        's1_orbit': os.path.join(testdir, 'S1A_IW_GRDH_1SDV_20210119T031653_20210119T031718_036201_043ED0_8255.zip'),
        'tif': os.path.join(testdir, 'S1A__IW___A_20150309T173017_VV_grd_mli_geo_norm_db.tif'),
        'archive_old_csv': os.path.join(testdir, 'archive_outdated.csv'),
        'archive_old_bbox': os.path.join(testdir, 'archive_outdated_bbox.db'),
        'dempar': os.path.join(testdir, 'dem.par'),
        'mlipar': os.path.join(testdir, 'mli.par')
    }
    return out


@pytest.fixture
def auxdata_dem_cases():
    cases = [('AW3D30', ['N050E010/N051E011.tar.gz']),
             ('SRTM 1Sec HGT', ['N51E011.SRTMGL1.hgt.zip']),
             ('SRTM 3Sec', ['srtm_39_02.zip']),
             ('TDX90m', ['DEM/N51/E010/TDM1_DEM__30_N51E011.zip'])]
    return cases


@pytest.fixture
def tmp_home(monkeypatch, tmp_path):
    home = tmp_path / 'tmp_home'
    home.mkdir()
    var = 'USERPROFILE' if platform.system() == 'Windows' else 'HOME'
    monkeypatch.setenv(var, str(home))
    assert os.path.expanduser('~') == str(home)
    yield home
