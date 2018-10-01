import os
import re
import subprocess as sp
from collections import Counter
from spatialist.ancillary import finder


def parse_command(command):
    '''
    Parse the help text of a Gamma command to a Python function including a docstring

    Parameters
    ----------
    command: str
        the name of the gamma command

    Returns
    -------
    str
        the full Python function text

    Examples
    --------
    >>> from pyroSAR.gamma.parser import parse_command
    >>> print(parse_command('offset_pwrm'))
    'def offset_pwrm(MLI_1, MLI_2, DIFF_par, offs, ccp, rwin='-', azwin='-', offsets='-', n_ovr='-', nr='-', naz='-', thres='-', c_ovr='-', pflag='-', pltflg='-', ccs='-', logpath=None):
    """
    Offset tracking between MLI images using intensity cross-correlation
    Copyright 2016, Gamma Remote Sensing, v4.8 clw 22-Octf-2016

    Args:
        MLI_1:     (input) real valued intensity image 1 (reference)
        MLI_2:     (input) real valued intensity image 2
        DIFF_par:  DIFF/GEO parameter file
        offs:      (output) offset estimates (fcomplex)
        ccp:       (output) cross-correlation of each patch (0.0->1.0) (float)
        rwin:      range patch size (range pixels, (enter - for default from offset parameter file)
        azwin:     azimuth patch size (azimuth lines, (enter - for default from offset parameter file)
        offs:ets   (output) range and azimuth offsets and cross-correlation data in text format, enter - for no output
        n_ovr:     MLI oversampling factor (integer 2**N (1,2,4,8), enter - for default: 2)
        nr:        number of offset estimates in range direction (enter - for default from offset parameter file)
        naz:       number of offset estimates in azimuth direction (enter - for default from offset parameter file)
        thres:     cross-correlation threshold (enter - for default from offset parameter file)
        c_ovr:     correlation function oversampling factor (integer 2**N (1,2,4,8) default: 4)
        pflag:     print flag
                  0:print offset summary
                  1:print all offset data
        pltflg:    plotting flag:
                  0: none (default)
                  1: screen output
                  2: screen output and PNG format plots
                  3: output plots in PDF format
        ccs:       (output) cross-correlation standard deviation of each patch (float)
    """
    process(['offset_pwrm', MLI_1, MLI_2, DIFF_par, offs, ccp, rwin, azwin, offsets, n_ovr, nr, naz, thres, c_ovr, pflag, pltflg, ccs], logpath=logpath)
    '
    '''

    proc = sp.Popen(command, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True)
    out, err = proc.communicate()
    out += err

    header = '\n'.join([x.strip('* ') for x in re.findall('[*]{3}.*[*]{3}', out)])
    print('header_raw: \n{}\n'.format(header))
    usage = re.search('usage:.*(?=\n)', out).group()
    print('usage_raw: \n{}\n'.format(usage))
    arg_req = [re.sub('[^\w.-]*', '', x) for x in re.findall('[^<]*<([^>]*)>', usage)]
    print('required args: {}\n'.format(', '.join(arg_req)))
    arg_opt = [re.sub('[^\w.-]*', '', x) for x in re.findall('[^[]*\[([^]]*)\]', usage)]
    print('optional args: {}\n'.format(', '.join(arg_opt)))
    double = [k for k, v in Counter(arg_req + arg_opt).items() if v > 1]
    print('double args: {}\n'.format(', '.join(double)))
    arg_str1 = re.sub(r'([^\'])-([^\'])', r'\1_\2', ', '.join(arg_req + [x + "='-'" for x in arg_opt])) \
        .replace(', def=', ', drm=')
    arg_str2 = ', '.join(arg_req + arg_opt) \
        .replace('-', '_') \
        .replace(', def,', ', drm,')
    doc_start = 'input parameters:[ ]*\n' if re.search('input parameters', out) else 'usage:.*(?=\n)'
    docstring = 'Args:\n' + out[re.search(doc_start, out).end():].strip('\n')
    for arg in arg_req + arg_opt:
        scan = re.finditer('\n[ ]*{}'.format(arg), docstring)
        for match in scan:
            string = match.group()
            new = re.sub(r'(\n[ ]*)({})'.format(arg), r'\1  \2:', string).replace('-', '_')
            if arg == 'def':
                new = new.replace('def', 'drm')

            # match = re.search('\n[ ]*{}'.format(arg), docstring)
            # if match:
            #     match = match.group()
            #     new = re.sub(r'(\n[ ]*)({})'.format(arg), r'\1  \2:', match).replace('-', '_')
            #     if arg == 'def':
            #         new = new.replace('def', 'drm')
            docstring = docstring.replace(string, new)
    fun = '''def {0}({1}, logpath=None):\n"""\n{2}\n\n{3}\n"""\nprocess(['{4}', {5}], logpath=logpath)''' \
        .format(os.path.basename(command).replace('-', '_'),
                arg_str1,
                header,
                docstring,
                command,
                arg_str2) \
        .replace('\n', '\n    ')
    return fun + '\n'


excludes = ['coord_trans', 'mosaic', 'lin_comb', 'lin_comb_cpx', 'validate']


def parse_module(bindir, outfile):
    """
    parse all Gamma commands of a module to functions and save them to a Python script.

    Parameters
    ----------
    bindir: str
        the `bin` directory of a module containing the commands
    outfile: str
        the name of the Pyxthon file to write

    Returns
    -------

    """
    failed = []
    outstring = 'from pyroSAR.gamma.auxil import process\n\n\n'
    for cmd in finder(bindir, ['*']):
        basename = os.path.basename(cmd)
        if basename not in excludes:
            print(basename)
            try:
                fun = parse_command(cmd)
            except (ValueError, AttributeError):
                failed.append(basename)
                continue
            outstring += fun + '\n\n'
    with open(outfile, 'w') as out:
        out.write(outstring)
    if len(failed) > 0:
        print('the following functions could not be parsed:\n{}\n'.format('\n'.join(failed)))
