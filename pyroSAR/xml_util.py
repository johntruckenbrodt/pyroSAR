##############################################################
# utility collection for xml files
# John Truckenbrodt 2016
##############################################################

from .ancillary import dictmerge


def getNamespaces(xml):
    infile = xml if 'readline' in dir(xml) else open(xml, 'r')
    # calling 'for line in infile:' will cause an error later on!
    namespaces = {}
    infile.seek(0)
    while True:
        line = infile.readline()
        if not line:
            break
        if 'xmlns:' in line:
            namespaces = dictmerge(namespaces,
                                   ([tuple(x.replace('xmlns:', '')
                                           .replace('"', '')
                                           .split('='))
                                     for x in line.split() if x.startswith('xmlns')]))
    # reset file pointer
    infile.seek(0)

    if isinstance(xml, str):
        infile.close()
    return namespaces
