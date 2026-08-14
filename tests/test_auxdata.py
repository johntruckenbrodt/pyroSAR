import os
import pytest
from pyroSAR.auxdata import dem_autoload, DEMHandler, dem_create

from spatialist import bbox


@pytest.mark.parametrize(
    'cases_fixture',
    [
        'auxdata_dem_cases_northern',
        'auxdata_dem_cases_northern_antimeridian',
        'auxdata_dem_cases_southern'
    ]
)
def test_handler(cases_fixture, request):
    extent, cases = request.getfixturevalue(cases_fixture)
    with bbox(coordinates=extent, crs=4326) as box:
        with DEMHandler(box) as handler:
            assert isinstance(handler.auxdatapath, str)
            for demType, reference in cases:
                result = handler.remote_ids(dem_type=demType)
                assert result == reference


def test_handler_fail():
    with pytest.raises(RuntimeError):
        test = DEMHandler('foobar')
    ext_utm = {'xmin': -955867, 'xmax': -915536, 'ymin': -5915518, 'ymax': -5863678}
    with bbox(ext_utm, crs=32632) as box:
        with pytest.raises(RuntimeError):
            test = DEMHandler([box])


def test_autoload(travis):
    with bbox({'xmin': 11.5, 'xmax': 11.9, 'ymin': 51, 'ymax': 51.5}, crs=4326) as box:
        # if the following is run in a loop, it is not possible to see which demType failed
        # Travis CI does not support ftp access;
        # see https://blog.travis-ci.com/2018-07-23-the-tale-of-ftp-at-travis-ci
        if not travis:
            files = dem_autoload(box, 'AW3D30')
            assert len(files) == 1
            files = dem_autoload(box, 'AW3D30', product='stk')
            assert len(files) == 1
        files = dem_autoload(box, 'SRTM 1Sec HGT')
        assert len(files) == 1
        files = dem_autoload(box, 'SRTM 1Sec HGT', offline=True)
        assert len(files) == 1
        files = dem_autoload(box, 'SRTM 3Sec')
        assert len(files) == 1
        with pytest.raises(RuntimeError):
            files = dem_autoload(box, 'TDX90m')
        with pytest.raises(RuntimeError):
            dem_autoload(box, 'AW3D30', product='foobar')
    with bbox({'xmin': -30, 'xmax': -29, 'ymin': 40, 'ymax': 41}, crs=4326) as box:
        files = dem_autoload(box, 'SRTM 1Sec HGT')
        assert len(files) == 0


def test_dem_create(tmpdir):
    vrt = '/vsimem/test.vrt'
    out = os.path.join(str(tmpdir), 'srtm.tif')
    with bbox({'xmin': 11.5, 'xmax': 11.9, 'ymin': 51, 'ymax': 51.5}, crs=4326) as box:
        with pytest.raises(RuntimeError):
            files = dem_autoload(geometry=box, demType='foobar')
        dem_autoload(geometry=box, demType='SRTM 3Sec', vrt=vrt, product='dem')
        dem_create(geometry=box, src=vrt, dst=out, t_srs=32632,
                   tr=(90, 90), nodata=-32767)
    assert os.path.isfile(out)


@pytest.mark.parametrize(
    'case, extent, step, expected',
    [
        (
                'one_degree_regular_extent',
                {'xmin': 11, 'xmax': 12, 'ymin': 51, 'ymax': 51.5},
                1,
                ([51], [11]),
        ),
        (
                'five_degree_regular_extent',
                {'xmin': 11, 'xmax': 12, 'ymin': 51, 'ymax': 51.5},
                5,
                ([50], [10]),
        ),
        (
                'fifteen_degree_regular_extent',
                {'xmin': 11, 'xmax': 12, 'ymin': 51, 'ymax': 51.5},
                15,
                ([45], [0]),
        ),
        (
                'one_degree_antimeridian_extent',
                {'xmin': 179, 'xmax': -179, 'ymin': 51, 'ymax': 51.5},
                1,
                ([51], [179, -180]),
        ),
        (
                'five_degree_antimeridian_extent',
                {'xmin': 175, 'xmax': -175, 'ymin': 51, 'ymax': 51.5},
                5,
                ([50], [175, -180]),
        ),
        (
                'fifteen_degree_antimeridian_extent',
                {'xmin': 165, 'xmax': -165, 'ymin': 51, 'ymax': 51.5},
                15,
                ([45], [165, -180]),
        ),
    ],
    ids=[
        'one_degree_regular_extent',
        'five_degree_regular_extent',
        'fifteen_degree_regular_extent',
        'one_degree_antimeridian_extent',
        'five_degree_antimeridian_extent',
        'fifteen_degree_antimeridian_extent',
    ],
)
def test_intrange(case, extent, step, expected):
    with bbox(coordinates=extent, crs=4326) as box:
        with DEMHandler(geometry=box) as dem:
            assert dem.intrange(step=step) == expected
