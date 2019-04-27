import os
from pyroSAR.gamma import ISPPar, par2hdr, Namespace


def test_par(testdata, tmpdir):
    with ISPPar(testdata['dempar']) as par:
        envi = par.envidict()
        assert envi['map_info'] == ['UTM', '1.0000', '1.0000', 515353.565, 5235168.873, '20.0', '20.0',
                                    32, 'North', 'WGS-84', 'units=Meters']
        assert envi['lines'] == 6455
        assert envi['samples'] == 5927
        assert envi['interleave'] == 'bsq'
        assert envi['bands'] == 1
        assert envi['byte_order'] == 1
        assert envi['data_type'] == 4
        assert envi['file_type'] == 'ENVI Standard'
        hdrfile = os.path.join(str(tmpdir), 'dem.hdr')
    par2hdr(testdata['dempar'], hdrfile=hdrfile, modifications={'band_names': ['band1']}, nodata=0)
    assert os.path.isfile(hdrfile)
    with ISPPar(testdata['mlipar']) as par:
        assert par.date == '2014-11-15T18:18:1.309100'
        assert par.envidict()['acquisition_time'] == '2014-11-15T18:18:1.309100Z'
        print(par)


def test_namespace():
    n = Namespace(directory='/test', basename='S1A__IW___A_20180829T170656')
    n.appreciate(['inc_geo', 'ls_map'])
    assert n.isregistered('inc_geo')
    assert n.isfile('inc_geo') is False
    assert n.isappreciated('inc_geo') is True
    exp1 = os.path.join('/test', 'S1A__IW___A_20180829T170656_inc_geo')
    exp2 = os.path.join('/test', 'S1A__IW___A_20180829T170656_ls_map')
    assert n['inc_geo'] == exp1
    assert n.get('ls_map') == exp2
    n.depreciate(['inc_geo'])
    assert n.isappreciated('inc_geo') is False
    assert n['inc_geo'] == '-'
    assert n.getall() == {'inc_geo': '-', 'ls_map': exp2}
    assert n.select(['inc_geo', 'ls_map']) == ['-', exp2]
    n.depreciate(['dem_seg'])
    assert n['dem_seg'] == '-'
