##############################################################
# utility collection for xml files
# John Truckenbrodt 2016-2018
##############################################################

import os
import re
import ast
import xml.etree.ElementTree as ET


class XMLHandler(object):
    def __init__(self, xml):
        errormessage = 'xmlfile must be a string pointing to an existing file, ' \
                       'a string or bytes object from which an xml can be parsed or a file object'
        if 'readline' in dir(xml):
            self.infile = xml.name if hasattr(xml, 'name') else None
            xml.seek(0)
            self.text = xml.read()
            xml.seek(0)
        elif isinstance(xml, (bytes, str)):
            try:
                isfile = os.path.isfile(xml)
            except ValueError:
                isfile = False
            if isfile:
                self.infile = xml
                with open(xml, 'r') as infile:
                    self.text = infile.read()
            else:
                try:
                    tree = ET.fromstring(xml)
                    self.infile = None
                    self.text = str(xml)
                    del tree
                except ET.ParseError:
                    raise RuntimeError(errormessage)
        else:
            raise RuntimeError(errormessage)
        defs = re.findall('xmlns:[a-z0-9]+="[^"]*"', self.text)
        dictstring = '{{{}}}'.format(re.sub(r'xmlns:([a-z0-9]*)=', r'"\1":', ', '.join(defs)))
        self.namespaces = ast.literal_eval(dictstring)

    def restoreNamespaces(self):
        for key, val in self.namespaces.items():
            val_new = val.split('/')[-1]
            self.text = self.text.replace(key, val_new)

    def write(self, outname, mode):
        with open(outname, mode) as out:
            out.write(self.text)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return


def getNamespaces(xmlfile):
    with XMLHandler(xmlfile) as xml:
        return xml.namespaces
