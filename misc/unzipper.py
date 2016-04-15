import os
import sys
import zipfile
from glob import glob


def finder(folder, matchlist):
    return list(set([f for files in [glob(os.path.join(item[0], pattern)) for item in os.walk(folder) for pattern in matchlist] for f in files]))

path = sys.argv[1]

items = finder(path, ["*.zip"])

for item in items:
    if not os.path.isdir(item[:-4]):
        print item
        with zipfile.ZipFile(item, "r") as z:
            z.extractall(os.path.dirname(item))