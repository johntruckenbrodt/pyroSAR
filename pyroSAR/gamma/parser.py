###############################################################################
# parse Gamma command docstrings to Python functions

# Copyright (c) 2015-2021, the pyroSAR Developers.

# This file is part of the pyroSAR Project. It is subject to the
# license terms in the LICENSE.txt file found in the top-level
# directory of this distribution and at
# https://github.com/johntruckenbrodt/pyroSAR/blob/master/LICENSE.txt.
# No part of the pyroSAR project, including this file, may be
# copied, modified, propagated, or distributed except according
# to the terms contained in the LICENSE.txt file.
###############################################################################
import os
import re
import subprocess as sp
from collections import Counter
from spatialist.ancillary import finder, which, dissolve

from pyroSAR.examine import ExamineGamma


def parse_command(command, indent='    '):
    """
    Parse the help text of a Gamma command to a Python function including a docstring.
    The docstring is in rst format and can thu be parsed by e.g. sphinx.
    This function is not intended to be used by itself, but rather within function :func:`parse_module`.

    Parameters
    ----------
    command: str
        the name of the gamma command
    indent: str
        the Python function indentation string; default: four spaces

    Returns
    -------
    str
        the full Python function text

    """
    # run the command without passing arguments to just catch its usage description
    command = which(command)
    if command is None:
        raise OSError('command does not exist')
    command_base = os.path.basename(command)
    proc = sp.Popen(command, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True)
    out, err = proc.communicate()
    # sometimes the description string is split between stdout and stderr
    # for the following commands stderr contains the usage description line, which is inserted into stdout
    if command_base in ['ras_pt', 'ras_data_pt', 'rasdt_cmap_pt']:
        out = out.replace(' ***\n ', ' ***\n ' + err)
    else:
        # for all other commands stderr is just appended to stdout
        out += err
    
    pattern = r'([\w\.]+ (?:has been|was) re(?:named to|placed(?: that [ \*\n]*|) by)(?: the ISP program|) [\w\.]+)'
    match = re.search(pattern, out)
    if match:
        raise DeprecationWarning('\n' + out)
    
    if re.search(r"Can't locate FILE/Path\.pm in @INC", out):
        raise RuntimeError('unable to parse Perl script')
    ###########################################
    # fix command-specific inconsistencies in parameter naming
    # in several commands the parameter naming in the usage description line does not match that of the docstring
    parnames_lookup = {'2PASS_INT': [('OFF_PAR', 'OFF_par')],
                       'adapt_filt': [('low_snr_thr', 'low_SNR_thr')],
                       'atm_mod2': [('rpt', 'report'),
                                    ('[mode]', '[model_atm]'),
                                    ('[model]', '[model_atm]'),
                                    ('model     atm', 'model_atm atm'),
                                    ],
                       'atm_mod_2d': [('xref', 'rref'),
                                      ('yref', 'azref')],
                       'atm_mod_2d_pt': [('[sigma_min]', '[sigma_max]')],
                       'cc_monitoring': [('...', '<...>')],
                       'cct_sp_pt': [('pcct_sp_pt', 'pcct_sp')],
                       'comb_interfs': [('combi_out', 'combi_int')],
                       'coord_to_sarpix': [('north/lat', 'north_lat'),
                                           ('east/lon', 'east_lon'),
                                           ('SLC_par', '<SLC_MLI_par>'),
                                           ('SLC/MLI_par', 'SLC_MLI_par')],
                       'base_calc': [('plt_flg', 'plt_flag'),
                                     ('pltflg', 'plt_flag')],
                       'base_init': [('<base>', '<baseline>')],
                       'base_plot': [('plt_flg', 'plt_flag'),
                                     ('pltflg', 'plt_flag')],
                       'dis2hgt': [('m/cycle', 'm_cycle')],
                       'discc': [('min_corr', 'cmin'),
                                 ('max_corr', 'cmax')],
                       'disp2ras': [('<list>', '<DISP_tab>')],
                       'dis_data': [('...', '<...>')],
                       'dispwr': [('data_type', 'dtype')],
                       'DORIS_vec': [('SLC_PAR', 'SLC_par')],
                       'gc_map_fd': [('fdtab', 'fd_tab')],
                       'gc_map_grd': [('<MLI_par>', '<GRD_par>')],
                       'geocode_back': [('<gc_map>', '<lookup_table>'),
                                        ('\n  gc_map ', '\n  lookup_table ')],
                       'GRD_to_SR': [('SLC_par', 'MLI_par')],
                       'haalpha': [('<alpha> <entropy>', '<alpha2> <entropy>'),
                                   ('alpha       (output)', 'alpha2      (output)')],
                       'histogram_ras': [('mean/stdev', 'mean_stdev')],
                       'hsi_color_scale': [('[chip]', '[chip_width]')],
                       'HUYNEN_DEC': [('T11_0', 'T11'),
                                      ('<T12> <T13> <T11>', '<T11> <T12> <T13>'),
                                      ('HUYNEN_DEC:', '***')],
                       'interf_SLC': [('  SLC2_pa  ', '  SLC2_par  ')],
                       'ionosphere_mitigation': [('<SLC1> <ID1>', '<ID1>')],
                       'landsat2dem': [('<DEM>', '<image>')],
                       'line_interp': [('input file', 'data_in'),
                                       ('output file', 'data_out')],
                       'm-alpha': [('<c2 ', '<c2> ')],
                       'm-chi': [('<c2 ', '<c2> ')],
                       'm-delta': [('<c2 ', '<c2> ')],
                       'map_section': [('n1', 'north1'),
                                       ('e1', 'east1'),
                                       ('n2', 'north2'),
                                       ('e2', 'east2'),
                                       ('[coord]', '[coords]')],
                       'mask_class': [('...', '<...>')],
                       'mcf_pt': [('<azlks>', '[azlks]'),
                                  ('<rlks>', '[rlks]')],
                       'mk_2d_im_geo': [('exponent', 'exp')],
                       'mk_adf2_2d': [('[alpha_max [', '[alpha_max] ['),
                                      ('-m MLI_dir', 'mli_dir'),
                                      ('-s scale', 'scale'),
                                      ('-e exp', 'exponent'),
                                      ('-u', 'update')],
                       'mk_base_calc': [('<RSLC_tab>', '<SLC_tab>')],
                       'mk_cpd_all': [('dtab', 'data_tab')],
                       'mk_cpx_ref_2d': [('diff_tab', 'cpx_tab')],
                       'mk_dispmap2_2d': [('RMLI_image', 'MLI'),
                                          ('RMLI_par', 'MLI_par'),
                                          ('MLI_image', 'MLI'),
                                          ('DISP_tab', 'disp_tab')],
                       'mk_dispmap_2d': [('RMLI_image', 'MLI'),
                                         ('RMLI_par', 'MLI_par'),
                                         ('MLI_image', 'MLI'),
                                         ('DISP_tab', 'disp_tab')],
                       'mk_geo_data_all': [('data_geo_dir', 'geo_dir')],
                       'mk_itab': [('<offset>', '<start>')],
                       'mk_hgt_2d': [('m/cycle', 'm_cycle')],
                       'mk_pol2rec_2d': [('data_tab', 'DIFF_tab'),
                                         ('<type> <rmli>', '<dtype>'),
                                         ('<dtype> <rmli>', '<dtype>'),
                                         ('type           input', 'dtype          input'),
                                         ('\n    Options:\n', ''),
                                         ('-s scale', 'scale'),
                                         ('-e exp', 'exponent'),
                                         ('-a min', 'min'),
                                         ('-b max', 'max'),
                                         ('-R rmax', 'rmax'),
                                         ('-m mode', 'mode'),
                                         ('-u', 'update')],
                       'mk_rasdt_all': [('RMLI_image', 'MLI'),
                                        ('MLI_image', 'MLI')],
                       'mk_rasmph_all': [('RMLI_image', 'MLI'),
                                         ('MLI_image', 'MLI')],
                       'mk_unw_2d': [('unw_mask1', 'unw_mask')],
                       'mk_unw_ref_2d': [('diff_tab', 'DIFF_tab')],
                       'MLI2pt': [('MLI_TAB', 'MLI_tab'),
                                  ('pSLC_par', 'pMLI_par')],
                       'mosaic': [('<..>', '<...>'),
                                  ('DEM_parout', 'DEM_par_out')],
                       'multi_class_mapping': [('...', '<...>')],
                       'multi_look_geo': [('geo_SLC', 'SLC'),
                                          ('SLC/MLI', ('SLC_MLI'))],
                       'multi_look_MLI': [('MLI in_par', 'MLI_in_par')],
                       'offset_fit': [('interact_flag', 'interact_mode')],
                       'offset_plot_az': [('rmin', 'r_min'),
                                          ('rmax', 'r_max')],
                       'par_ASF_SLC': [('CEOS_SAR_leader', 'CEOS_leader')],
                       'par_ASAR': [('ASAR/ERS_file', 'ASAR_ERS_file')],
                       'par_EORC_JERS_SLC': [('slc', 'SLC')],
                       'par_ERSDAC_PALSAR': [('VEXCEL_SLC_par', 'ERSDAC_SLC_par')],
                       'par_ESA_JERS_SEASAT_SLC': [('[slc]', '[SLC]')],
                       'par_ICEYE_GRD': [('<GeoTIFF>', '<GeoTIFF> <XML>'),
                                         ('[mli]', '[MLI]')],
                       'par_ICEYE_SLC': [('[slc]', '[SLC]')],
                       'par_MSP': [('SLC/MLI_par', 'SLC_MLI_par')],
                       'par_SIRC': [('UTC/MET', 'UTC_MET')],
                       'par_TX_GRD': [('COSAR', 'GeoTIFF')],
                       'par_UAVSAR_SLC': [('SLC/MLC_in', 'SLC_MLC_in'),
                                          ('SLC/MLI_par', 'SLC_MLI_par'),
                                          ('SLC/MLI_out', 'SLC_MLI_out')],
                       'par_UAVSAR_geo': [('SLC/MLI_par', 'SLC_MLI_par')],
                       'phase_sim': [('sim       (', 'sim_unw   (')],
                       'product': [('wgt_flg', 'wgt_flag')],
                       'radcal_MLI': [('MLI_PAR', 'MLI_par')],
                       'radcal_PRI': [('GRD_PAR', 'GRD_par'),
                                      ('PRI_PAR', 'PRI_par')],
                       'radcal_SLC': [('SLC_PAR', 'SLC_par')],
                       'ras2jpg': [('{', '{{'),
                                   ('}', '}}')],
                       'ras_data_pt': [('pdata1', 'pdata')],
                       'ras_to_rgb': [('red channel', 'red_channel'),
                                      ('green channel', 'green_channel'),
                                      ('blue channel', 'blue_channel')],
                       'rascc_mask_thinning': [('...', '[...]')],
                       'rashgt': [('m/cycle', 'm_cycle')],
                       'rashgt_shd': [('m/cycle', 'm_cycle'),
                                      ('\n  cycle ', '\n  m_cycle ')],
                       'rasdt_cmap_pt': [('pdata1', 'pdata')],
                       'raspwr': [('hdrz', 'hdrsz')],
                       'ras_ras': [('r_lin/log', 'r_lin_log'),
                                   ('g_lin/log', 'g_lin_log'),
                                   ('b_lin/log', 'b_lin_log')],
                       'ras_ratio_dB': [('[min_cc] [max_cc] [scale] [exp]', '[min_value] [max_value] [dB_offset]')],
                       'rasSLC': [('[header]', '[hdrsz]')],
                       'ratio': [('wgt_flg', 'wgt_flag')],
                       'restore_float': [('input file', 'data_in'),
                                         ('output file', 'data_out'),
                                         ('interpolation_limit', 'interp_limit')],
                       'S1_coreg_TOPS_no_refinement': [('RLK', 'rlks'),
                                                       ('AZLK', 'azlks')],
                       'S1_OPOD_vec': [('SLC_PAR', 'SLC_par')],
                       'single_class_mapping': [('>...', '> <...>')],
                       'ScanSAR_burst_cc_ad': [('bx', 'box_min'),
                                               ('by', 'box_max')],
                       'ScanSAR_burst_to_mosaic': [('DATA_tab_ref', 'data_tab_ref'),
                                                   ('[mflg] [dtype]', '[mflg]')],
                       'ScanSAR_full_aperture_SLC': [('SLCR_dir', 'SLC2_dir')],
                       'scale_base': [('SLC-1_par-2', 'SLC1_par-2')],
                       'sigma2gamma': [('<gamma>', '<gamma0>'),
                                       ('gamma  (output)', 'gamma0  (output)'),
                                       ('pwr1', 'sigma0')],
                       'SLC_interp_lt': [('SLC-2', 'SLC2'),
                                         ('blksz', 'blk_size')],
                       'SLC_intf': [('SLC1s_par', 'SLC-1s_par'),
                                    ('SLC2Rs_par', 'SLC-2Rs_par')],
                       'SLC_intf_geo2': [('cc        (', 'CC        (')],
                       'SLC_interp_map': [('coffs2_sm', 'coffs_sm')],
                       'srtm_mosaic': [('<lon>', '<lon2>')],
                       'SSI_INT_S1': [('<SLC2> <par2>', '<SLC_tab2>')],
                       'texture': [('weights_flag', 'wgt_flag')],
                       'ts_rate': [('sim_flg', 'sim_flag')],
                       'TX_SLC_preproc': [('TX_list', 'TSX_list')],
                       'uchar2float': [('infile', 'data_in'),
                                       ('outfile', 'data_out')],
                       'validate': [('ras1', 'ras_map'),
                                    ('rasf_map', 'ras_map'),
                                    ('ras2', 'ras_inv'),
                                    ('rasf_inventory', 'ras_inv'),
                                    ('class1[1]', 'class1_1'),
                                    ('class1[2]', 'class1_2'),
                                    ('class1[n]', 'class1_n'),
                                    ('class2[1]', 'class2_1'),
                                    ('class2[2]', 'class2_2'),
                                    ('class2[n]', 'class2_n')]}
    if command_base in parnames_lookup.keys():
        for replacement in parnames_lookup[command_base]:
            out = out.replace(*replacement)
    ###########################################
    # filter header (general command description) and usage description string
    header = '\n'.join([x.strip('* ') for x in re.findall('[*]{3}.*(?:[*]{3}|)', out)])
    header = '| ' + header.replace('\n', '\n| ')
    usage = re.search('usage:.*(?=\n)', out).group()
    
    # filter required and optional arguments from usage description text
    arg_req_raw = [re.sub(r'[^\w.-]*', '', x) for x in re.findall('[^<]*<([^>]*)>', usage)]
    arg_opt_raw = [re.sub(r'[^\w.-]*', '', x) for x in re.findall(r'[^[]*\[([^]]*)\]', usage)]
    
    ###########################################
    # add parameters missing in the usage argument lists
    
    appends = {'mk_adf2_2d': ['cc_min', 'cc_max', 'mli_dir', 'scale', 'exponent', 'update'],
               'mk_pol2rec_2d': ['scale', 'exponent', 'min', 'max', 'rmax', 'mode', 'update'],
               'SLC_interp_S1_TOPS': ['mode', 'order'],
               'SLC_interp_map': ['mode', 'order']}
    
    if command_base in appends.keys():
        for var in appends[command_base]:
            arg_opt_raw.append(var)
    ###########################################
    # define parameter replacements; this is intended for parameters which are to be aggregated into a list parameter
    replacements = {'cc_monitoring': [(['nfiles', 'f1', 'f2', '...'],
                                       ['files'],
                                       ['a list of input data files (float)'])],
                    'dis_data': [(['nstack', 'pdata1', '...'],
                                  ['pdata'],
                                  ['a list of point data stack files'])],
                    'lin_comb': [(['nfiles', 'f1', 'f2', '...'],
                                  ['files'],
                                  ['a list of input data files (float)']),
                                 (['factor1', 'factor2', '...'],
                                  ['factors'],
                                  ['a list of factors to multiply the input files with'])],
                    'lin_comb_cpx': [(['nfiles', 'f1', 'f2', '...'],
                                      ['files'],
                                      ['a list of input data files (float)']),
                                     (['factor1_r', 'factor2_r', '...'],
                                      ['factors_r'],
                                      ['a list of real part factors to multiply the input files with']),
                                     (['factor1_i', 'factor2_i'],
                                      ['factors_i'],
                                      ['a list of imaginary part factors to multiply the input files with'])],
                    'mask_class': [(['n_class', 'class_1', '...', 'class_n'],
                                    ['class_values'],
                                    ['a list of class map values'])],
                    'mosaic': [(['nfiles', 'data_in1', 'DEM_par1', 'data_in2', 'DEM_par2', '...', '...'],
                                ['data_in_list', 'DEM_par_list'],
                                ['a list of input data files',
                                 'a list of DEM/MAP parameter files for each data file'])],
                    'multi_class_mapping': [(['nfiles', 'f1', 'f2', '...', 'fn'],
                                             ['files'],
                                             ['a list of input data files (float)'])],
                    'rascc_mask_thinning': [(['thresh_1', '...', 'thresh_nmax'],
                                             ['thresholds'],
                                             ['a list of thresholds sorted from smallest to '
                                              'largest scale sampling reduction'])],
                    'single_class_mapping': [(['nfiles', 'f1', '...', 'fn'],
                                              ['files'],
                                              ['a list of point data stack files']),
                                             (['lt1', 'ltn'],
                                              ['thres_lower'],
                                              ['a list of lower thresholds for the files']),
                                             (['ut1', 'utn'],
                                              ['thres_upper'],
                                              ['a list of upper thresholds for the files'])],
                    'validate': [(['nclass1', 'class1_1', 'class1_2', '...', 'class1_n'],
                                  ['classes_map'],
                                  ['a list of class values for the map data file (max. 16), 0 for all']),
                                 (['nclass2', 'class2_1', 'class2_2', '...', 'class2_n'],
                                  ['classes_inv'],
                                  ['a list of class values for the inventory data file (max. 16), 0 for all'])]}
    
    if '..' in usage and command_base not in replacements.keys():
        raise RuntimeError('the command contains multi-args which were not properly parsed')
    
    def replace(inlist, replacement):
        outlist = list(inlist)
        for old, new, description in replacement:
            if old[0] not in outlist:
                return outlist
            outlist[outlist.index(old[0])] = new
            for i in range(1, len(old)):
                if old[i] in outlist:
                    outlist.remove(old[i])
        return dissolve(outlist)
    
    arg_req = list(arg_req_raw)
    arg_opt = list(arg_opt_raw)
    
    if command_base in replacements.keys():
        arg_req = replace(arg_req, replacements[command_base])
        arg_opt = replace(arg_opt, replacements[command_base])
    
    if command_base in ['par_CS_geo', 'par_KS_geo']:
        out = re.sub('[ ]*trunk.*', '', out, flags=re.DOTALL)
    ###########################################
    # check if there are any double parameters
    
    double = [k for k, v in Counter(arg_req + arg_opt).items() if v > 1]
    if len(double) > 0:
        raise RuntimeError('double parameter{0}: {1}'.format('s' if len(double) > 1 else '', ', '.join(double)))
    ###########################################
    # add a parameter inlist for commands which take interactive input via stdin
    # the list of commands, which are interactive is hard to assess and thus likely a source of future errors
    
    inlist = ['create_dem_par', 'par_ESA_ERS']
    
    if command_base in inlist:
        arg_req.append('inlist')
    
    ######################################################################################
    # create the function argument string for the Python function
    
    # optional arguments are parametrized with '-' as default value, e.g. arg_opt='-'
    # a '-' in the parameter name is replaced with '_'
    # example: "arg1, arg2, arg3='-'"
    argstr_function = re.sub(r'([^\'])-([^\'])', r'\1_\2', ', '.join(arg_req + [x + "='-'" for x in arg_opt])) \
        .replace(', def=', ', drm=')
    
    # create the function definition string
    fun_def = 'def {name}({args_fun}, logpath=None, outdir=None, shellscript=None):' \
        .format(name=command_base.replace('-', '_'),
                args_fun=argstr_function)
    
    if command_base == '2PASS_INT':
        fun_def = fun_def.replace(command_base, 'TWO_PASS_INT')
    ######################################################################################
    # special handling of flag args
    flag_args = {'mk_adf2_2d': [('mli_dir', '-m', None),
                                ('scale', '-s', None),
                                ('exponent', '-e', None),
                                ('update', '-u', False)],
                 'mk_pol2rec_2d': [('scale', '-s', None),
                                   ('exp', '-e', None),
                                   ('min', '-a', None),
                                   ('max', '-b', None),
                                   ('rmax', '-R', None),
                                   ('mode', '-m', None),
                                   ('update', '-u', False)]}
    
    # replace arg default like arg='-' with arg=None or arg=False
    if command_base in flag_args:
        for arg in flag_args[command_base]:
            fun_def = re.sub('{}=\'-\''.format(arg[0]), '{0}={1}'.format(arg[0], arg[2]), fun_def)
    ######################################################################################
    # create the process call argument string
    
    # a '-' in the parameter name is replaced with '_'
    # e.g. 'arg1, arg2, arg3'
    # if a parameter is named 'def' (not allowed in Python) it is renamed to 'drm'
    
    # inlist is not a proc arg but a parameter passed to function process
    proc_args = arg_req + arg_opt
    if command_base in inlist:
        proc_args.remove('inlist')
    proc_args_tmp = list(proc_args)
    # insert the length of a list argument as a proc arg
    if command_base in replacements.keys() and command_base != 'rascc_mask_thinning':
        key = replacements[command_base][0][1]
        if isinstance(key, list):
            key = key[0]
        proc_args_tmp.insert(proc_args_tmp.index(key), 'len({})'.format(key))
    
    if command_base == 'validate':
        index = proc_args_tmp.index('classes_inv')
        proc_args_tmp.insert(index, 'len(classes_inv)')
    
    argstr_process = ', '.join(proc_args_tmp) \
        .replace('-', '_') \
        .replace(', def,', ', drm,')
    
    # create the process argument list string
    cmd_str = "cmd = ['{command}', {args_cmd}]".format(command=command, args_cmd=argstr_process)
    
    # special handling of optional flag args
    # the args are removed from the cmd list and flags (plus values) added if not None or True
    # e.g. '-u' if update=True or '-m /path' if mli_dir='/path'
    if command_base in flag_args:
        args = []
        for arg in flag_args[command_base]:
            cmd_str = cmd_str.replace(', {}'.format(arg[0]), '')
            args.append(arg[0])
            cmd_str += "\nif {a} is not {d}:\n{i}cmd.append('{k}')" \
                .format(i=indent, d=arg[2], k=arg[1], a=arg[0])
            if arg[2] is None:
                cmd_str += '\n{i}cmd.append({a})'.format(i=indent, a=arg[0])
    
    # create the process call string
    proc_str = "process(cmd, logpath=logpath, outdir=outdir{inlist}, shellscript=shellscript)" \
        .format(inlist=', inlist=inlist' if command_base in inlist else '')
    fun_proc = '{0}\n{1}'.format(cmd_str, proc_str)
    
    if command_base == 'lin_comb_cpx':
        fun_proc = fun_proc.replace('factors_r, factors_i', 'zip(factors_r, factors_i)')
    elif command_base == 'mosaic':
        fun_proc = fun_proc.replace('data_in_list, DEM_par_list', 'zip(data_in_list, DEM_par_list)')
    elif command_base == 'single_class_mapping':
        fun_proc = fun_proc.replace('files, thres_lower, thres_upper', 'zip(files, thres_lower, thres_upper)')
    
    ######################################################################################
    # create the function docstring
    
    # find the start of the docstring and filter the result
    doc_start = 'input parameters:[ ]*\n' if re.search('input parameters', out) else 'usage:.*(?=\n)'
    doc = '\n' + out[re.search(doc_start, out).end():]
    
    # define a pattern containing individual parameter documentations
    pattern = r'\n[ ]*[<\[]*(?P<par>{0})[>\]]*[\t ]+(?P<doc>.*)'.format(
        '|'.join(arg_req_raw + arg_opt_raw).replace('.', r'\.'))
    
    # identify the start indices of all pattern matches
    starts = [m.start(0) for m in re.finditer(pattern, doc)] + [len(out)]
    
    # filter out all individual (parameter, description) docstring tuples
    doc_items = []
    j = 0
    done = []
    for i in range(0, len(starts) - 1):
        doc_raw = doc[starts[i]:starts[i + 1]]
        doc_list = list(re.search(pattern, doc_raw, flags=re.DOTALL).groups())
        
        if doc_list[0] not in proc_args:
            if command_base in replacements.keys():
                repl = replacements[command_base][0]
                for k, item in enumerate(repl[1]):
                    if item not in done:
                        doc_items.append([item, repl[2][k]])
                        done.append(item)
                        j += 1
            continue
        
        if doc_list[0] in done:
            doc_items[-1][1] += doc_raw
            continue
        
        while doc_list[0] != proc_args[j]:
            doc_list_sub = [proc_args[j], 'not documented']
            doc_items.append(doc_list_sub)
            j += 1
        
        doc_items.append(doc_list)
        done.append(doc_items[-1][0])
        j += 1
    
    for k in range(j, len(proc_args)):
        doc_items.append([proc_args[k], 'not documented'])
    
    # add a parameter inlist to the docstring tuples
    if command_base in inlist:
        pos = [x[0] for x in doc_items].index(arg_opt[0])
        doc_items.insert(pos, ('inlist', 'a list of arguments to be passed to stdin'))
    
    # remove the replaced parameters from the argument lists
    doc_items = [x for x in doc_items if x[0] in arg_req + arg_opt]
    
    # replace parameter names which are not possible in Python syntax, i.e. containing '-' or named 'def'
    for i, item in enumerate(doc_items):
        par = item[0].replace('-', '_').replace(', def,', ', drm,')
        description = item[1]
        doc_items[i] = (par, description)
    
    if command_base in ['par_CS_geo', 'par_KS_geo']:
        doc_items.append(('MLI_par', '(output) ISP SLC/MLI parameter file (example: yyyymmdd.mli.par)'))
        doc_items.append(('DEM_par', '(output) DIFF/GEO DEM parameter file (example: yyyymmdd.dem_par)'))
        doc_items.append(('GEO', '(output) Geocoded image data file (example: yyyymmdd.geo)'))
    
    # check if all parameters are documented:
    proc_args = [x.replace('-', '_').replace(', def,', ', drm,') for x in arg_req + arg_opt]
    mismatch = [x for x in proc_args if x not in [y[0] for y in doc_items]]
    if len(mismatch) > 0:
        raise RuntimeError('parameters missing in docsring: {}'.format(', '.join(mismatch)))
    ###########################################
    # format the docstring parameter descriptions
    
    docstring_elements = ['Parameters\n----------']
    
    # do some extra formatting
    for i, item in enumerate(doc_items):
        par, description = item
        description = re.split(r'\n+\s*', description.strip('\n'))
        
        # escape * characters (which are treated as special characters for bullet lists by sphinx)
        description = [x.replace('*', r'\*') for x in description]
        
        # convert all lines starting with an integer number or 'NOTE' to bullet list items
        latest = None
        for i in range(len(description)):
            item = description[i]
            if re.search('^(?:(?:-|)[-0-9]+|NOTE):', item):
                latest = i
                # prepend '* ' and replace missing spaces after a colon: 'x:x' -> 'x: x'
                description[i] = '* ' + re.sub(r'((?:-|)[-0-9]+:)(\w+)', r'\1 \2', item)
        
        # format documentation lines coming after the last bullet list item
        # sphinx expects lines after the last bullet item to be indented by two spaces if
        # they belong to the bullet item or otherwise a blank line to mark the end of the bullet list
        if latest:
            # case if there are still lines coming after the last bullet item,
            # prepend an extra two spaces to these lines so that they are properly
            # aligned with the text of the bullet item
            if latest + 2 <= len(description):
                i = 1
                while latest + i + 1 <= len(description):
                    description[latest + i] = '  ' + description[latest + i]
                    i += 1
            # if not, then insert an extra blank line
            else:
                description[-1] = description[-1] + '\n'
        
        # parse the final documentation string for the current parameter
        description = '\n{0}{0}'.join(description).format(indent)
        doc = '{0}:\n{1}{2}'.format(par, indent, description)
        docstring_elements.append(doc)
    ###########################################
    # add docsrings of general parameters and combine the result
    
    # create docstring for parameter logpath
    doc = 'logpath: str or None\n{0}a directory to write command logfiles to'.format(indent)
    docstring_elements.append(doc)
    
    # create docstring for parameter outdir
    doc = 'outdir: str or None\n{0}the directory to execute the command in'.format(indent)
    docstring_elements.append(doc)
    
    # create docstring for parameter shellscript
    doc = 'shellscript: str or None\n{0}a file to write the Gamma commands to in shell format'.format(indent)
    docstring_elements.append(doc)
    
    # combine the complete docstring
    fun_doc = '\n{header}\n\n{doc}\n' \
        .format(header=header,
                doc='\n'.join(docstring_elements))
    ######################################################################################
    
    # combine the elements to a complete Python function string
    fun = '''{defn}\n"""{doc}"""\n{proc}'''.format(defn=fun_def, doc=fun_doc, proc=fun_proc)
    
    # indent all lines and add an extra empty line at the end
    fun = fun.replace('\n', '\n{}'.format(indent)) + '\n'
    
    return fun


def parse_module(bindir, outfile):
    """
    parse all Gamma commands of a module to functions and save them to a Python script.

    Parameters
    ----------
    bindir: str
        the `bin` directory of a module containing the commands
    outfile: str
        the name of the Python file to write

    Returns
    -------
    
    Examples
    --------
    >>> import os
    >>> from pyroSAR.gamma.parser import parse_module
    >>> outname = os.path.join(os.environ['HOME'], 'isp.py')
    >>> parse_module('/cluster/GAMMA_SOFTWARE-20161207/ISP/bin', outname)
    """
    
    if not os.path.isdir(bindir):
        raise OSError('directory does not exist: {}'.format(bindir))
    
    excludes = ['coord_trans',  # doesn't take any parameters and is interactive
                'RSAT2_SLC_preproc',  # takes option flags
                'mk_ASF_CEOS_list',  # "cannot create : Directory nonexistent"
                '2PASS_UNW',  # parameter name inconsistencies
                'mk_diff_2d',  # takes option flags
                'gamma_doc'  # opens the Gamma documentation
                ]
    failed = []
    outstring = ''
    for cmd in sorted(finder(bindir, [r'^\w+$'], regex=True), key=lambda s: s.lower()):
        basename = os.path.basename(cmd)
        if basename not in excludes:
            # print(basename)
            try:
                fun = parse_command(cmd)
            except RuntimeError as e:
                failed.append('{0}: {1}'.format(basename, str(e)))
                continue
            except DeprecationWarning:
                continue
            except:
                failed.append('{0}: {1}'.format(basename, 'error yet to be assessed'))
                continue
            outstring += fun + '\n\n'
    if len(outstring) > 0:
        if not os.path.isfile(outfile):
            with open(outfile, 'w') as out:
                out.write('from pyroSAR.gamma.auxil import process\n\n\n')
        with open(outfile, 'a') as out:
            out.write(outstring)
    if len(failed) > 0:
        print('the following functions could not be parsed:\n{0}\n({1} total)'.format('\n'.join(failed), len(failed)))


def autoparse():
    """
    automatic parsing of Gamma commands.
    This function will detect the Gamma installation via environment variable `GAMMA_HOME`, detect all available
    modules (e.g. ISP, DIFF) and parse all of the module's commands via function :func:`parse_module`.
    A new Python module will be created called `gammaparse`, which is stored under `$HOME/.pyrosar`.
    Upon importing the `pyroSAR.gamma` submodule, this function is run automatically and module `gammaparse`
    is imported as `api`.
    
    Returns
    -------

    Examples
    --------
    >>> from pyroSAR.gamma.api import diff
    >>> print('create_dem_par' in dir(diff))
    True
    """
    home = ExamineGamma().home
    target = os.path.join(os.path.expanduser('~'), '.pyrosar', 'gammaparse')
    if not os.path.isdir(target):
        os.makedirs(target)
    for module in finder(home, ['[A-Z]*'], foldermode=2):
        outfile = os.path.join(target, os.path.basename(module).lower() + '.py')
        if not os.path.isfile(outfile):
            print('parsing module {} to {}'.format(os.path.basename(module), outfile))
            for submodule in ['bin', 'scripts']:
                print('-' * 10 + '\n{}'.format(submodule))
                try:
                    parse_module(os.path.join(module, submodule), outfile)
                except OSError:
                    print('..does not exist')
            print('=' * 20)
    modules = [re.sub(r'\.py', '', os.path.basename(x)) for x in finder(target, [r'[a-z]+\.py$'], regex=True)]
    if len(modules) > 0:
        with open(os.path.join(target, '__init__.py'), 'w') as init:
            init.write('from . import {}'.format(', '.join(modules)))
