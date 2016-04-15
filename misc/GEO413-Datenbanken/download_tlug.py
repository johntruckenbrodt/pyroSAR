import urllib2, csv, os

os.chdir("/media/user/Data/Dokumente/Studium/Master/GEO413-Geodatenbanken/GEO413_data")
with open('./TLUG_info/Pegel.csv', 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter='\t', quotechar='|')
    for row in reader:
        name_w = row[0] + "_w_28.txt"
        name_q = row[0] + "_q_28.txt"
        url = "http://www.tlug-jena.de/hw.inc/txt/" + name_w
        print url
        try:
            textfile = urllib2.urlopen(url)
            output = open("./TLUG/Wasserstand/" + name_w, 'wb')
            output.write(textfile.read())
            output.close()
        except:
            print "download error for file", name_w

        url = "http://www.tlug-jena.de/hw.inc/txt/" + name_q
        print url
        try:
            textfile = urllib2.urlopen(url)
            output = open("./TLUG/Durchfluss/" + name_q, 'wb')
            output.write(textfile.read())
            output.close()
        except:
            print "download error for file", name_q
        print "----"
