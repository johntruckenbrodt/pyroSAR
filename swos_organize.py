
import os
import shutil
import tarfile as tf
from ancillary import finder, multicore
from archivist import scan, tar2zip

dir_c1f = "/geonfs01_vol3/swos/data/ers_envisat"
dir_scene = "/geonfs01_vol1/ve39vem/swos"

areas = ['02_Spain_LagunaDeFuenteDePiedra',
         '04_Jordan_Azraq',
         '07_Sweden_StoreMosse-Kavsjon',
         '13_Greece_Makedonia',
         '15_Egypt_LakeBurullus',
         '17_France_Carmargue',
         '18_Kenya_LakeVictoria',
         '22_Tanzania_Kilombero']

sensors = ["ERS", "ASAR"]


def operator(package, outdir, archive):
    tar = tf.open(package)
    for item in tar.getnames():
        try:
            member = tf.open(fileobj=tar.extractfile(item))
            basename = os.path.basename(scan(member, "[NE][12]$")[0])
            print basename
            outname = os.path.join(outdir, basename)
            tar2zip(member, outname)
            member.close()
        except tf.ReadError:
            print item
            tar2zip(tar, os.path.join(outdir, item), membernames=[item])
    tar.close()
    if not os.path.isdir(archive):
        os.makedirs(archive)
    shutil.move(package, archive)

packages = []
for area in areas:
    for sensor in sensors:
        subdir = os.path.join(dir_c1f, area, sensor)
        if not os.path.isdir(subdir):
            os.makedirs(subdir)

        packages += finder(os.path.join(dir_c1f, area, sensor), ["C1F(.)*\.tar(?:\.gz|)$"], regex=True, recursive=False)

archives = [os.path.join(os.path.dirname(x), "archived") for x in packages]

print "processing {} packages".format(len(packages))

multicore(operator, 20, {"package": packages, "archive": archives}, outdir=dir_scene)
