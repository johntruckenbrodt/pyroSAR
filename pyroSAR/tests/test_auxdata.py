import pytest
from pyroSAR.auxdata import dem_autoload, DEMHandler

from spatialist import bbox


def test_handler():
    with bbox({'xmin': 11.5, 'xmax': 11.9, 'ymin': 51.1, 'ymax': 51.5}, crs=4326) as box:
        with DEMHandler([box]) as handler:
            cases = [('AW3D30', ['N050E010_N055E015.tar.gz']),
                     ('SRTM 1Sec HGT', ['N51E011.SRTMGL1.hgt.zip']),
                     ('SRTM 3Sec', ['srtm_39_02.zip']),
                     ('TDX90m', ['90mdem/DEM/N51/E010/TDM1_DEM__30_N51E011.zip'])]
            for demType, reference in cases:
                result = handler.remote_ids(demType=demType, extent=box.extent)
                assert result == reference
    
    with bbox({'xmin': -11.9, 'xmax': -11.5, 'ymin': -51.5, 'ymax': -51.1}, crs=4326) as box:
        with DEMHandler([box]) as handler:
            cases = [('AW3D30', ['S055W015_S050W010.tar.gz']),
                     ('SRTM 1Sec HGT', ['S52W012.SRTMGL1.hgt.zip']),
                     ('SRTM 3Sec', ['srtm_34_23.zip']),
                     ('TDX90m', ['90mdem/DEM/S52/W010/TDM1_DEM__30_S52W012.zip'])]
            for demType, reference in cases:
                result = handler.remote_ids(demType=demType, extent=box.extent)
                assert result == reference


def test_autoload(tmpdir):
    with bbox({'xmin': 11.5, 'xmax': 11.9, 'ymin': 51, 'ymax': 51.5}, crs=4326) as box:
        files = dem_autoload([box], 'SRTM 1Sec HGT')
        assert len(files) == 1
        with pytest.raises(RuntimeError):
            files = dem_autoload([box], 'TDX90m')
