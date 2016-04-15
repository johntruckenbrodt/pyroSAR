
import os
import re
import zipfile
from urllib2 import urlopen

outdir = "/geonfs02_vol1/SRTM_1_HGT"

server = "http://e4ftl01.cr.usgs.gov/SRTM/SRTMGL1.003/2000.02.11/"

pattern = "[NS][0-9]{2}[EW][0-9]{3}"

response = urlopen(server).read()

items = sorted(set(re.findall(pattern+".SRTMGL1.hgt.zip", response)))

for item in items:

    basename = re.search(pattern, item).group()

    local_zip = os.path.join(outdir, item)
    local_unzip = os.path.join(outdir, basename + ".hgt")

    if not os.path.isfile(local_unzip):
        print item

        infile = urlopen(os.path.join(server, item))
        with open(local_zip, "wb") as outfile:
            outfile.write(infile.read())
        infile.close()
        with zipfile.ZipFile(local_zip, "r") as z:
            z.extractall(outdir)
        os.remove(local_zip)
