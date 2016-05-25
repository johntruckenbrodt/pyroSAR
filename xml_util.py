##############################################################
# utility collection for xml files
# John Truckenbrodt 2016
##############################################################


def getNamespaces(xml):
    namespaces = {}
    infile = xml if "readline" in dir(xml) else open(xml, "r")
    # calling 'for line in infile:' will cause an error later on!
    while True:
        line = infile.readline()
        if not line:
            break
        # todo: find out how this can be generalized
        if line.startswith("<xfdu"):
            namespaces = dict([tuple(x.replace("xmlns:", "").replace('"', '').split("=")) for x in line.split() if x.startswith("xmlns")])
            break
    # reset file pointer
    infile.seek(0)

    if isinstance(xml, str):
        infile.close()
    return namespaces
