import os
from pyroSAR.gamma import ISPPar, par2hdr, Namespace


def test_par(testdata, tmpdir):
    with ISPPar(testdata['dempar']) as par:
        envi = par.envidict()
        print(envi)
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


def test_namespace():
    n = Namespace(directory='/test', basename='S1A__IW___A_20180829T170656')
    n.appreciate(['inc_geo', 'ls_map'])
    assert n.isregistered('inc_geo')
    assert n['inc_geo'] == '/test/S1A__IW___A_20180829T170656_inc_geo'
    assert n.get('ls_map') == '/test/S1A__IW___A_20180829T170656_ls_map'
    n.depreciate(['inc_geo'])
    assert n.isappreciated('inc_geo') is False
    assert n['inc_geo'] == '-'
    assert n.getall() == {'inc_geo': '-',
                          'ls_map': '/test/S1A__IW___A_20180829T170656_ls_map'}
    assert n.select(['inc_geo', 'ls_map']) == ['-', '/test/S1A__IW___A_20180829T170656_ls_map']
