import pytest
from pyroSAR.ancillary import getargs
from pyroSAR.gamma import api


@pytest.mark.skipif('diff' not in dir(api), reason='requires GAMMA installation with module DIFF')
def test_args_diff():
    from pyroSAR.gamma.api import diff
    assert getargs(diff.gc_map) == ['DEM', 'DEM_par', 'DEM_seg', 'DEM_seg_par', 'MLI_par', 'OFF_par', 'frame',
                                    'inc', 'lat_ovr', 'logpath', 'lon_ovr', 'lookup_table', 'ls_map', 'ls_mode',
                                    'outdir', 'pix', 'psi', 'r_ovr', 'shellscript', 'sim_sar', 'u', 'v']
    assert getargs(diff.gc_map_grd) == ['DEM', 'DEM_par', 'DEM_seg', 'DEM_seg_par', 'GRD_par', 'frame', 'inc',
                                        'lat_ovr', 'logpath', 'lon_ovr', 'lookup_table', 'ls_map', 'ls_mode', 'outdir',
                                        'pix', 'psi', 'r_ovr', 'shellscript', 'sim_sar', 'u', 'v']
    
    args = getargs(diff.geocode_back)
    args_ref = ['data_in', 'data_out', 'dtype', 'interp_mode', 'logpath', 'lookup_table',
                             'lr_in', 'lr_out', 'nlines_out', 'order', 'outdir', 'shellscript',
                             'width_in', 'width_out']
    comp = [x in args for x in args_ref]
    assert sum(comp) == len(args_ref)
    
    assert getargs(diff.par_EORC_PALSAR_geo) == ['CEOS_data', 'CEOS_leader', 'DEM_par', 'MLI',
                                                 'MLI_par', 'cal', 'logpath', 'outdir', 'shellscript']
    assert getargs(diff.par_TX_geo) == ['DEM_par', 'GEO', 'GeoTIFF', 'MLI_par', 'annotation_XML',
                                        'logpath', 'outdir', 'pol', 'shellscript']
    
    args = getargs(diff.pixel_area)
    args_ref = ['DEM', 'DEM_par', 'MLI_par', 'area_fact', 'inc_map', 'logpath', 'lookup_table',
                           'ls_map', 'nstep', 'outdir', 'pix_gamma0', 'pix_sigma0', 'shellscript']
    comp = [x in args for x in args_ref]
    assert sum(comp) == len(args_ref)


@pytest.mark.skipif('disp' not in dir(api), reason='requires GAMMA installation with module DISP')
def test_args_disp():
    from pyroSAR.gamma.api import disp
    
    args = getargs(disp.data2geotiff)
    args_ref = ['DEM_par', 'GeoTIFF', 'data', 'logpath', 'no_data', 'outdir', 'shellscript',
                             'type']
    comp = [x in args for x in args_ref]
    assert sum(comp) == len(args_ref)


@pytest.mark.skipif('isp' not in dir(api), reason='requires GAMMA installation with module ISP')
def test_args_isp():
    from pyroSAR.gamma.api import isp
    assert getargs(isp.multi_look) == ['MLI', 'MLI_par', 'SLC', 'SLC_par', 'azlks', 'exp', 'loff', 'logpath',
                                       'nlines', 'outdir', 'rlks', 'scale', 'shellscript']
    
    args = getargs(isp.multi_look_MLI)
    args_ref = ['MLI_in', 'MLI_in_par', 'MLI_out', 'MLI_out_par', 'azlks', 'loff',
                'logpath', 'nlines', 'outdir', 'rlks', 'scale', 'shellscript']
    comp = [x in args for x in args_ref]
    assert sum(comp) == len(args_ref)
    
    assert getargs(isp.par_ASAR) == ['ASAR_ERS_file', 'K_dB', 'logpath', 'outdir', 'output_name', 'shellscript']
    assert getargs(isp.par_EORC_PALSAR) == ['CEOS_data', 'CEOS_leader', 'SLC', 'SLC_par', 'dtype',
                                            'logpath', 'outdir', 'sc_dB', 'shellscript']
    assert getargs(isp.par_ESA_ERS) == ['CEOS_DAT', 'CEOS_SAR_leader', 'SLC', 'SLC_par', 'inlist',
                                        'logpath', 'outdir', 'shellscript']
    
    args = getargs(isp.par_S1_GRD)
    args_ref = ['GRD', 'GRD_par', 'GeoTIFF', 'MLI', 'MLI_par', 'annotation_XML',
                'calibration_XML', 'eflg', 'logpath', 'noise_XML', 'noise_pwr',
                'outdir', 'rps', 'shellscript']
    comp = [x in args for x in args_ref]
    assert sum(comp) == len(args_ref)
    
    assert getargs(isp.par_S1_SLC) == ['GeoTIFF', 'SLC', 'SLC_par', 'TOPS_par', 'annotation_XML', 'calibration_XML',
                                       'dtype', 'logpath', 'noise_XML', 'noise_pwr', 'outdir', 'sc_dB', 'shellscript']
    assert getargs(isp.par_TX_GRD) == ['GRD', 'GRD_par', 'GeoTIFF', 'annotation_XML', 'logpath',
                                       'outdir', 'pol', 'shellscript']
    assert getargs(isp.par_TX_SLC) == ['COSAR', 'SLC', 'SLC_par', 'annotation_XML', 'dtype',
                                       'logpath', 'outdir', 'pol', 'shellscript']
    assert getargs(isp.radcal_MLI) == ['CMLI', 'K_dB', 'MLI', 'MLI_par', 'OFF_par', 'ant_flag', 'antenna', 'logpath',
                                       'outdir', 'pix_area', 'refarea_flag', 'rloss_flag', 'sc_dB', 'shellscript']
    assert getargs(isp.radcal_PRI) == ['GRD', 'GRD_par', 'K_dB', 'PRI', 'PRI_par',
                                       'inc_ref', 'loff', 'logpath', 'nl', 'nr',
                                       'outdir', 'roff', 'shellscript']
    assert getargs(isp.radcal_SLC) == ['CSLC', 'CSLC_par', 'K_dB', 'SLC', 'SLC_par',
                                       'ant_flag', 'antenna', 'fcase', 'logpath', 'outdir',
                                       'pix_area', 'refarea_flag', 'rloss_flag', 'sc_dB', 'shellscript']
    assert getargs(isp.S1_OPOD_vec) == ['OPOD', 'SLC_par', 'logpath', 'nstate', 'outdir', 'shellscript']
    
    args = getargs(isp.SLC_deramp_ScanSAR)
    args_ref = ['SLC1_tab', 'SLC2_tab', 'logpath', 'mode', 'outdir',
                'phflg', 'shellscript']
    comp = [x in args for x in args_ref]
    assert sum(comp) == len(args_ref)
    
    args_ref = ['SLC', 'SLCR_tab', 'SLC_par', 'SLC_tab', 'azlks', 'logpath',
                'outdir', 'rlks', 'shellscript', 'bflg']
    args = getargs(isp.SLC_mosaic_S1_TOPS)
    comp = [x in args for x in args_ref]
    assert sum(comp) == len(args_ref)


@pytest.mark.skipif('lat' not in dir(api), reason='requires GAMMA installation with module LAT')
def test_args_lat():
    from pyroSAR.gamma.api import lat
    assert getargs(lat.linear_to_dB) == ['data_in', 'data_out', 'inverse_flag', 'logpath', 'null_value', 'outdir',
                                         'shellscript', 'width']
    assert getargs(lat.product) == ['bx', 'by', 'data_1', 'data_2', 'logpath', 'outdir', 'product',
                                    'shellscript', 'wgt_flag', 'width']
    assert getargs(lat.ratio) == ['bx', 'by', 'd1', 'd2', 'logpath', 'outdir', 'ratio',
                                  'shellscript', 'wgt_flag', 'width']

    args = getargs(lat.sigma2gamma)
    args_ref = ['gamma0', 'inc', 'logpath', 'outdir', 'sigma0', 'shellscript', 'width']
    comp = [x in args for x in args_ref]
    assert sum(comp) == len(args_ref)
