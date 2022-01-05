import os
import pytest
from pyroSAR.auxdata import dem_autoload, DEMHandler, dem_create

from spatialist import bbox


def test_handler(auxdata_dem_cases):
    with bbox({'xmin': 11.5, 'xmax': 11.9, 'ymin': 51.1, 'ymax': 51.5}, crs=4326) as box:
        with DEMHandler([box]) as handler:
            for demType, reference in auxdata_dem_cases:
                result = handler.remote_ids(demType=demType, extent=box.extent)
                assert result == reference
    
    with bbox({'xmin': -11.9, 'xmax': -11.5, 'ymin': -51.5, 'ymax': -51.1}, crs=4326) as box:
        with DEMHandler([box]) as handler:
            cases = [('AW3D30', ['S055W015/S052W012.tar.gz']),
                     ('SRTM 1Sec HGT', ['S52W012.SRTMGL1.hgt.zip']),
                     ('SRTM 3Sec', ['srtm_34_23.zip']),
                     ('TDX90m', ['90mdem/DEM/S52/W010/TDM1_DEM__30_S52W012.zip'])]
            for demType, reference in cases:
                result = handler.remote_ids(demType=demType, extent=box.extent)
                assert result == reference
    with pytest.raises(RuntimeError):
        test = DEMHandler('foobar')
    ext_utm = {'xmin': -955867, 'xmax': -915536, 'ymin': -5915518, 'ymax': -5863678}
    with bbox(ext_utm, crs=32632) as box:
        with pytest.raises(RuntimeError):
            test = DEMHandler([box])


def test_autoload(auxdata_dem_cases, travis):
    # delete all target files to test downloading them again
    home = os.path.expanduser('~')
    demdir = os.path.join(home, '.snap', 'auxdata', 'dem')
    locals = [os.path.join(demdir, x, os.path.basename(y[0])) for x, y in auxdata_dem_cases]
    for item in locals:
        if os.path.isfile(item):
            os.remove(item)
    with bbox({'xmin': 11.5, 'xmax': 11.9, 'ymin': 51, 'ymax': 51.5}, crs=4326) as box:
        # if the following is run in a loop, it is not possible to see which demType failed
        # Travis CI does not support ftp access;
        # see https://blog.travis-ci.com/2018-07-23-the-tale-of-ftp-at-travis-ci
        if not travis:
            files = dem_autoload([box], 'AW3D30')
            assert len(files) == 1
            files = dem_autoload([box], 'AW3D30', product='stk')
            assert len(files) == 1
        files = dem_autoload([box], 'SRTM 1Sec HGT')
        assert len(files) == 1
        files = dem_autoload([box], 'SRTM 3Sec')
        assert len(files) == 1
        with pytest.raises(RuntimeError):
            files = dem_autoload([box], 'TDX90m')
        with pytest.raises(RuntimeError):
            dem_autoload([box], 'AW3D30', product='foobar')


def test_dem_create(tmpdir):
    with bbox({'xmin': 11.5, 'xmax': 11.9, 'ymin': 51, 'ymax': 51.5}, crs=4326) as box:
        with pytest.raises(RuntimeError):
            files = dem_autoload([box], 'foobar')
        vrt = '/vsimem/test.vrt'
        dem_autoload([box], 'SRTM 3Sec', vrt=vrt)
    out = os.path.join(str(tmpdir), 'srtm.tif')
    dem_create(src=vrt, dst=out, t_srs=32632, tr=(90, 90))
    assert os.path.isfile(out)
