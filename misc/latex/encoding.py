

filename = "D:/Dokumente/Studium/Master/GEO512-Masterarbeit/Latex/Masterarbeit.bib"

with open(filename, "r") as infile:
    i = 1
    for line in infile:
        try:
            line.encode()
        except:
            print line
            print i
        i += 1









