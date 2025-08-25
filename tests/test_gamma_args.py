import pytest
from pyroSAR.ancillary import getargs
from pyroSAR.gamma import api


@pytest.mark.skipif('diff' not in dir(api), reason='requires GAMMA installation with module DIFF')
def test_args_diff():
    from pyroSAR.gamma.api import diff
    lookup = {
        'gc_map': ['DEM', 'DEM_par', 'DEM_seg', 'DEM_seg_par', 'MLI_par',
                   'OFF_par', 'frame', 'inc', 'lat_ovr', 'logpath', 'lon_ovr',
                   'lookup_table', 'ls_map', 'ls_mode', 'outdir', 'pix',
                   'psi', 'r_ovr', 'shellscript', 'sim_sar', 'u', 'v'],
        'gc_map_grd': ['DEM', 'DEM_par', 'DEM_seg', 'DEM_seg_par', 'GRD_par',
                       'frame', 'inc', 'lat_ovr', 'logpath', 'lon_ovr',
                       'lookup_table', 'ls_map', 'ls_mode', 'outdir', 'pix',
                       'psi', 'r_ovr', 'shellscript', 'sim_sar', 'u', 'v'],
        'geocode_back': ['data_in', 'data_out', 'dtype', 'interp_mode',
                         'logpath', 'lookup_table', 'lr_in', 'lr_out',
                         'nlines_out', 'order', 'outdir', 'shellscript',
                         'width_in', 'width_out'],
        'par_EORC_PALSAR_geo': ['CEOS_data', 'CEOS_leader', 'DEM_par', 'MLI',
                                'MLI_par', 'cal', 'logpath', 'outdir',
                                'shellscript'],
        'par_TX_geo': ['DEM_par', 'GEO', 'GeoTIFF', 'MLI_par', 'annotation_XML',
                       'logpath', 'outdir', 'pol', 'shellscript'],
        'pixel_area': ['DEM', 'DEM_par', 'MLI_par', 'area_fact', 'inc_map',
                       'logpath', 'lookup_table', 'ls_map', 'nstep', 'outdir',
                       'pix_gamma0', 'pix_sigma0', 'shellscript'],
    }
    for command, args in lookup.items():
        assert set(args).issubset(getargs(getattr(diff, command)))


@pytest.mark.skipif('disp' not in dir(api), reason='requires GAMMA installation with module DISP')
def test_args_disp():
    from pyroSAR.gamma.api import disp
    lookup = {
        'data2geotiff': ['DEM_par', 'GeoTIFF', 'data', 'logpath', 'no_data',
                         'outdir', 'shellscript', 'type']
    }
    for command, args in lookup.items():
        assert set(args).issubset(getargs(getattr(disp, command)))


@pytest.mark.skipif('isp' not in dir(api), reason='requires GAMMA installation with module ISP')
def test_args_isp():
    from pyroSAR.gamma.api import isp
    lookup = {
        'multi_look': ['MLI', 'MLI_par', 'SLC', 'SLC_par', 'azlks', 'exp',
                       'loff', 'nlines', 'rlks', 'scale'],
        'multi_look_MLI': ['MLI_in', 'MLI_in_par', 'MLI_out', 'MLI_out_par',
                           'azlks', 'loff', 'nlines', 'rlks', 'scale'],
        'par_ASAR': ['ASAR_ERS_file', 'K_dB', 'output_name'],
        'par_EORC_PALSAR': ['CEOS_data', 'CEOS_leader', 'SLC', 'SLC_par', 'dtype',
                            'sc_dB'],
        'par_ESA_ERS': ['CEOS_DAT', 'CEOS_SAR_leader', 'SLC', 'SLC_par', 'inlist'],
        'par_S1_GRD': ['GRD', 'GRD_par', 'GeoTIFF', 'MLI', 'MLI_par', 'annotation_XML',
                       'calibration_XML', 'eflg', 'noise_XML', 'noise_pwr', 'rps'],
        'par_S1_SLC': ['GeoTIFF', 'SLC', 'SLC_par', 'TOPS_par', 'annotation_XML',
                       'calibration_XML', 'dtype', 'noise_XML',
                       'noise_pwr', 'sc_dB'],
        'par_TX_GRD': ['GRD', 'GRD_par', 'GeoTIFF', 'annotation_XML', 'pol'],
        'par_TX_SLC': ['COSAR', 'SLC', 'SLC_par', 'annotation_XML', 'dtype', 'pol'],
        'radcal_MLI': ['CMLI', 'K_dB', 'MLI', 'MLI_par', 'OFF_par', 'ant_flag',
                       'antenna', 'pix_area', 'refarea_flag',
                       'rloss_flag', 'sc_dB'],
        'radcal_PRI': ['GRD', 'GRD_par', 'K_dB', 'PRI', 'PRI_par',
                       'inc_ref', 'loff', 'nl', 'nr', 'roff'],
        'radcal_SLC': ['CSLC', 'CSLC_par', 'K_dB', 'SLC', 'SLC_par',
                       'ant_flag', 'antenna', 'fcase',
                       'pix_area', 'refarea_flag', 'rloss_flag', 'sc_dB'],
        'S1_OPOD_vec': ['OPOD', 'SLC_par', 'nstate'],
        'SLC_deramp_ScanSAR': ['SLC1_tab', 'SLC2_tab', 'mode', 'phflg'],
        'SLC_mosaic_ScanSAR': ['SLC', 'SLCR_tab', 'SLC_par', 'SLC_tab',
                               'azlks', 'rlks', 'bflg']
    }
    for command, args in lookup.items():
        default = ['logpath', 'outdir', 'shellscript']
        assert set(args + default).issubset(getargs(getattr(isp, command)))


@pytest.mark.skipif('lat' not in dir(api), reason='requires GAMMA installation with module LAT')
def test_args_lat():
    from pyroSAR.gamma.api import lat
    lookup = {
        'linear_to_dB': ['data_in', 'data_out', 'inverse_flag', 'logpath', 'null_value', 'outdir',
                         'shellscript', 'width'],
        'product': ['bx', 'by', 'data_1', 'data_2', 'logpath', 'outdir', 'product',
                    'shellscript', 'wgt_flag', 'width'],
        'ratio': ['bx', 'by', 'd1', 'd2', 'logpath', 'outdir', 'ratio',
                  'shellscript', 'wgt_flag', 'width'],
        'sigma2gamma': ['gamma0', 'inc', 'logpath', 'outdir', 'sigma0', 'shellscript', 'width']
    }
    for command, args in lookup.items():
        assert set(args).issubset(getargs(getattr(lat, command)))
