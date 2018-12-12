import os
import re
import subprocess as sp
from collections import Counter
from spatialist.ancillary import finder, which

from .auxil import ExamineGamma


def parse_command(command):
    """
    Parse the help text of a Gamma command to a Python function including a docstring.
    The docstring is in rst format and can thu be parsed by e.g. sphinx.
    This function is not intended to be used by itself, but rather within function :func:`parse_module`.

    Parameters
    ----------
    command: str
        the name of the gamma command

    Returns
    -------
    str
        the full Python function text

    """
    command = which(command)
    command_base = os.path.basename(command)
    proc = sp.Popen(command, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True)
    out, err = proc.communicate()
    out += err
    
    # fix general inconsistencies in parameter naming
    parnames_lookup = {'product': ('wgt_flg', 'wgt_flag'),
                       'par_ASAR': ('ASAR/ERS_file', 'ASAR_ERS_file')}
    if command_base in parnames_lookup.keys():
        out = out.replace(*parnames_lookup[command_base])
    
    # filter header command description and usage description text
    header = '\n'.join([x.strip('* ') for x in re.findall('[*]{3}.*[*]{3}', out)])
    header = '| ' + header.replace('\n', '\n| ')
    usage = re.search('usage:.*(?=\n)', out).group()
    
    # filter required and optional arguments from usage description text
    arg_req_raw = [re.sub('[^\w.-]*', '', x) for x in re.findall('[^<]*<([^>]*)>', usage)]
    arg_opt = [re.sub('[^\w.-]*', '', x) for x in re.findall('[^[]*\[([^]]*)\]', usage)]
    
    # define parameter replacements
    replacements = {'lin_comb': [(['nfiles', 'f1', 'f2', '...'],
                                  ['files'],
                                  'a list of input data files (float)'),
                                 (['factor1', 'factor2', '...'],
                                  ['factors'],
                                  'a list of factors to multiply the input files with')]}
    
    def replace(sequence, replacement, lst):
        lst_repl = ','.join(lst).replace(','.join(sequence), ','.join(replacement))
        return lst_repl.split(',')
    
    arg_req = arg_req_raw
    if command_base in replacements.keys():
        for old, new, description in replacements['lin_comb']:
            arg_req = replace(old, new, arg_req)
    
    double = [k for k, v in Counter(arg_req + arg_opt).items() if v > 1]
    if len(double) > 0:
        raise RuntimeError('double parameter{0}: {1}'.format('s' if len(double) > 1 else '', ', '.join(double)))

    inlist = ['par_ESA_ERS']
    
    if command_base in inlist:
        arg_req.append('inlist')
    
    # print('header_raw: \n{}\n'.format(header))
    # print('usage_raw: \n{}\n'.format(usage))
    # print('required args: {}\n'.format(', '.join(arg_req)))
    # print('optional args: {}\n'.format(', '.join(arg_opt)))
    # print('double args: {}\n'.format(', '.join(double)))
    ######################################################################################
    # create the function argument string for the Python function
    
    # optional arguments are parametrized with '-' as default value, e.g. arg_opt='-'
    # a '-' in the parameter name is replaced with '_'
    # example: "arg1, arg2, arg3='-'"
    argstr_function = re.sub(r'([^\'])-([^\'])', r'\1_\2', ', '.join(arg_req + [x + "='-'" for x in arg_opt])) \
        .replace(', def=', ', drm=')
    
    # create the function definition string
    fun_def = 'def {name}({args_fun}, logpath=None, outdir=None, shellscript=None):' \
        .format(name=os.path.basename(command).replace('-', '_'),
                args_fun=argstr_function)
    ######################################################################################
    # create the process call argument string
    
    # a '-' in the parameter name is replaced with '_'
    # e.g. 'arg1, arg2, arg3'
    # if a parameter is named 'def' (not allowed in Python) it is renamed to 'drm'
    
    proc_args = arg_req + arg_opt
    
    if command_base in inlist:
        proc_args.remove('inlist')
    
    if command_base == 'lin_comb':
        proc_args.insert(0, 'len(files)')
    
    argstr_process = ', '.join(proc_args) \
        .replace('-', '_') \
        .replace(', def,', ', drm,')
    
    # create the process call string
    fun_proc = "process(['{command}', {args_cmd}], logpath=logpath, outdir=outdir{inlist}, shellscript=shellscript)" \
        .format(command=command,
                args_cmd=argstr_process,
                inlist=', inlist=inlist' if command_base in inlist else '')
    
    # print('arg_str1: \n{}\n'.format(argstr_function))
    # print('arg_str2: \n{}\n'.format(argstr_process))
    ######################################################################################
    # create the function docstring
    
    # filter out all individual (parameter, description) docstring tuples
    doc_start = 'input parameters:[ ]*\n' if re.search('input parameters', out) else 'usage:.*(?=\n)'
    doc = '\n' + out[re.search(doc_start, out).end():]
    pattern = r'\n[ ]*[<]*(?P<par>{0})[>]*[ ]+(?P<doc>.*)'.format('|'.join(arg_req_raw + arg_opt).replace('.', '\.'))
    starts = [m.start(0) for m in re.finditer(pattern, doc)] + [len(out)]
    doc_items = []
    for i in range(0, len(starts) - 1):
        doc_raw = doc[starts[i]:starts[i + 1]]
        doc_tuple = re.search(pattern, doc_raw, flags=re.DOTALL).groups()
        doc_items.append(doc_tuple)
    
    if command_base in inlist:
        pos = [x[0] for x in doc_items].index(arg_opt[0])
        doc_items.insert(pos, ('inlist', 'a list of arguments to be passed to stdin'))
    
    # insert the parameter replacements
    if command_base in replacements.keys():
        for old, new, description in replacements[command_base]:
            index = [doc_items.index(x) for x in doc_items if x[0] == old[0]][0]
            doc_items.insert(index, (new[0], description))
    
    # remove the replaced parameters from the docstring list
    doc_items = [x for x in doc_items if x[0] in arg_req + arg_opt]
    
    # replace parameter names
    for i, item in enumerate(doc_items):
        par = item[0].replace('-', '_').replace(', def,', ', drm,')
        description = item[1]
        doc_items[i] = (par, description)
    
    # check if all parameters are documented:
    proc_args = [x.replace('-', '_').replace(', def,', ', drm,') for x in arg_req + arg_opt]
    mismatch = [x for x in [y[0] for y in doc_items] if x not in proc_args]
    if len(mismatch) > 0:
        print(arg_req + arg_opt)
        print([y[0] for y in doc_items])
        raise RuntimeError('parameters missing in docsring: {}'.format(', '.join(mismatch)))
    ###########################################
    # format the docstring parameter descriptions
    
    docstring_elements = ['Parameters\n----------']
    # define the number of spaces to indent
    indent = ' ' * 4
    # do some extra formatting
    for i, item in enumerate(doc_items):
        par, description = item
        description = re.split('\n+\s*', description.strip('\n'))
        # escape * characters (which are treated as special characters for bullet lists by sphinx)
        description = [x.replace('*', '\*') for x in description]
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
    excludes = ['coord_trans', 'mosaic', 'lin_comb_cpx', 'validate']
    failed = []
    outstring = 'from pyroSAR.gamma.auxil import process\n\n\n'
    for cmd in sorted(finder(bindir, ['*']), key=lambda s: s.lower()):
        basename = os.path.basename(cmd)
        if basename not in excludes:
            # print(basename)
            try:
                fun = parse_command(cmd)
            except RuntimeError as e:
                failed.append('{0}: {1}'.format(basename, str(e)))
                continue
            outstring += fun + '\n\n'
    with open(outfile, 'w') as out:
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
            print('parsing module {}'.format(os.path.basename(module)))
            parse_module(os.path.join(module, 'bin'), outfile)
            print('=' * 20)
    modules = [re.sub('\.py', '', os.path.basename(x)) for x in finder(target, ['[a-z]+\.py$'], regex=True)]
    if len(modules) > 0:
        with open(os.path.join(target, '__init__.py'), 'w') as init:
            init.write('from . import {}'.format(', '.join(modules)))
