import os
import re
import subprocess as sp
# from collections import Counter
from spatialist.ancillary import finder, which


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
    proc = sp.Popen(command, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True)
    out, err = proc.communicate()
    out += err
    
    # filter header command description and usage description text
    header = '\n'.join([x.strip('* ') for x in re.findall('[*]{3}.*[*]{3}', out)])
    header = '| ' + header.replace('\n', '\n| ')
    usage = re.search('usage:.*(?=\n)', out).group()
    
    # filter required and optional arguements from usage description text
    arg_req = [re.sub('[^\w.-]*', '', x) for x in re.findall('[^<]*<([^>]*)>', usage)]
    arg_opt = [re.sub('[^\w.-]*', '', x) for x in re.findall('[^[]*\[([^]]*)\]', usage)]
    
    for arg in arg_req + arg_opt:
        for item in re.findall(arg, out, re.IGNORECASE):
            if item != arg:
                out = out.replace(item, arg)
    
    # double = [k for k, v in Counter(arg_req + arg_opt).items() if v > 1]
    
    # print('header_raw: \n{}\n'.format(header))
    # print('usage_raw: \n{}\n'.format(usage))
    # print('required args: {}\n'.format(', '.join(arg_req)))
    # print('optional args: {}\n'.format(', '.join(arg_opt)))
    # print('double args: {}\n'.format(', '.join(double)))
    
    # create the function argument string for the Python function
    # optional arguments are parametrized with '-' as default value, e.g. arg_opt='-'
    # a '-' in the parameter name is replaced with '_'
    # example: "arg1, arg2, arg3='-'"
    argstr_function = re.sub(r'([^\'])-([^\'])', r'\1_\2', ', '.join(arg_req + [x + "='-'" for x in arg_opt])) \
        .replace(', def=', ', drm=')
    
    argstr_process = ', '.join(arg_req + arg_opt) \
        .replace('-', '_') \
        .replace(', def,', ', drm,')
    
    # print('arg_str1: \n{}\n'.format(argstr_function))
    # print('arg_str2: \n{}\n'.format(argstr_process))
    
    # define the start of the parameter documentation string, which is either after 'input_parameters' or after
    # the usage description string
    doc_start = 'input parameters:[ ]*\n' if re.search('input parameters', out) else 'usage:.*(?=\n)'
    
    # parse the parameter documentation to a Python docstring format
    
    tabspace = ' ' * 4
    
    docstring_elements = ['Parameters\n----------']
    
    starts = []
    for x in arg_req + arg_opt:
        try:
            starts.append(re.search(r'\n[ ]*{0} .*'.format(x), out).start())
        except AttributeError:
            raise RuntimeError('cannot find parameter {}'.format(x))
    starts += [len(out)]
    # starts = [re.search(r'\n[ ]*{0}.*'.format(x), out).start() for x in arg_req + arg_opt] + [len(out)]

    pattern = r'\n[ ]*(?P<command>{0})[ ]+(?P<doc>.*)'.format('|'.join(arg_req + arg_opt))
    # print(pattern)
    
    for i in range(0, len(starts) - 1):
        doc_raw = out[starts[i]:starts[i + 1]]
        # print(repr(doc_raw))
        match = re.match(pattern, doc_raw, flags=re.DOTALL)
        if not match:
            continue
        cmd = match.group('command')
        doc_items = re.split('\n+\s*', match.group('doc').strip('\n'))
        doc_items = [x.replace('*', '\*') for x in doc_items]
        latest = None
        for i in range(len(doc_items)):
            item = doc_items[i]
            if re.search('^(?:(?:-|)[-0-9]+|NOTE):', item):
                latest = i
                doc_items[i] = '* ' + re.sub(r'((?:-|)[-0-9]+:)(\w+)', r'\1 \2', item)
        if latest:
            # case if there are still lines coming after the last bullet item, prepend an extra two spaces to these
            # lines so that they are properly aligned with the text of the bullet item
            if latest + 2 <= len(doc_items):
                i = 1
                while latest + i + 1 <= len(doc_items):
                    doc_items[latest + i] = '  ' + doc_items[latest + i]
                    i += 1
            # if not, the insert an extra blank line
            else:
                doc_items[-1] = doc_items[-1] + '\n'
        
        description = '\n{0}{0}'.join(doc_items).format(tabspace)
        doc = '{1}:\n{0}{2}'.format(tabspace, cmd, description)
        docstring_elements.append(doc)
    
    docstring = '\n'.join(docstring_elements)
    
    # combine the elements to the final Python function string
    fun = '''def {name}({args_fun}, logpath=None):\n"""\n{header}\n\n{doc}\n"""\nprocess(['{command}', {args_cmd}], logpath=logpath)''' \
        .format(name=os.path.basename(command).replace('-', '_'),
                args_fun=argstr_function,
                header=header,
                doc=docstring,
                command=command,
                args_cmd=argstr_process) \
        .replace('\n', '\n    ')
    return fun + '\n'


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
    excludes = ['coord_trans', 'mosaic', 'lin_comb', 'lin_comb_cpx', 'validate']
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
