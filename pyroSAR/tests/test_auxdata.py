import pytest
from pyroSAR.auxdata import dem_autoload, DEMHandler

from spatialist import bbox


def test_handler():
    with bbox({'xmin': 11.5, 'xmax': 11.9, 'ymin': 51.1, 'ymax': 51.5}, crs=4326) as box:
        with DEMHandler([box]) as handler:
            cases = [('AW3D30', ['N050E010/N051E011.tar.gz']),
                     ('SRTM 1Sec HGT', ['N51E011.SRTMGL1.hgt.zip']),
                     ('SRTM 3Sec', ['srtm_39_02.zip']),
                     ('TDX90m', ['90mdem/DEM/N51/E010/TDM1_DEM__30_N51E011.zip'])]
            for demType, reference in cases:
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


def test_autoload(tmpdir):
    with bbox({'xmin': 11.5, 'xmax': 11.9, 'ymin': 51, 'ymax': 51.5}, crs=4326) as box:
        for type in ['AW3D30', 'SRTM 1Sec HGT', 'SRTM 3Sec']:
            files = dem_autoload([box], type)
            assert len(files) == 1
        with pytest.raises(RuntimeError):
            files = dem_autoload([box], 'TDX90m')
        with pytest.raises(RuntimeError):
            dem_autoload([box], 'AW3D30', product='foobar')
        files = dem_autoload([box], 'AW3D30', product='stk')
        assert len(files) == 1
