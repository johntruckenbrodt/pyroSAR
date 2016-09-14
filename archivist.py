##############################################################
# utility collection for compressed archives
# John Truckenbrodt 2016
##############################################################

import os
import re
import calendar
import zipfile as zf
import tarfile as tf
from time import gmtime
from datetime import datetime
from ancillary import finder


# foldermodes:
# 0: no folders
# 1: folders included
def scan(archive, pattern, foldermode=1):
    files = []
    if isinstance(archive, str):
        archive = os.path.realpath(archive)
        handle = zf.ZipFile(archive, "r") if zf.is_zipfile(archive) else tf.open(archive)
    else:
        handle = archive
    if isinstance(handle, zf.ZipFile):
        files = [os.path.join(handle.filename, x) for x in handle.namelist() if re.search(pattern, x.strip("/"))]
        if foldermode == 0:
            files = [x for x in files if not x.endswith("/")]
    elif isinstance(handle, tf.TarFile):
        files = [os.path.join(handle.name, x) for x in handle.getnames() if re.search(pattern, x)]
        if foldermode == 0:
            files = [x for x in files if not handle.getmember(x).isdir()]
    if isinstance(archive, str):
        handle.close()
    return files


def compress(src, dstname, compression="zip"):
    if compression not in ["zip", "tar.gz"]:
        raise IOError("compression type not supported; use either 'zip' or 'tar.gz'")
    if isinstance(src, list) or os.path.isfile(src):
        candidates = [src] if os.path.isfile(src) else src
        targets = [os.path.basename(x) for x in candidates]
    elif os.path.isdir(src):
        candidates = [src]+finder(src, ["*"], 1)
        targets = [os.path.relpath(x, os.path.join(src, "..")) for x in candidates]
    else:
        candidates, targets = [], []
    if len(candidates) == 0:
        raise IOError("nothing found to compress")
    if not dstname.endswith("."+compression):
        dstname += "."+compression
    if compression == "zip":
        handle = zf.ZipFile(dstname, "w", zf.ZIP_DEFLATED)
        for i in range(len(candidates)):
            handle.write(candidates[i], targets[i])
        handle.close()
    elif compression == "tar.gz":
        handle = tf.open(dstname, "w:gz")
        if os.path.isdir(src):
            handle.add(candidates[0], targets[0])
        else:
            for i in range(len(candidates)):
                handle.add(candidates[i], targets[i])
        handle.close()


# adapted from 'http://stackoverflow.com/questions/18432565/python-convert-compressed-zip-to-uncompressed-tar'
def zip2tar(zipfile, tarname, membernames=None):
    zip = zipfile if isinstance(zipfile, zf.ZipFile) else zf.ZipFile(zipfile)
    if not tarname.endswith(".tar.gz"):
        tarname += ".tar.gz"
    tar = tf.open(tarname, "w:gz")
    timeshift = int((datetime.now() - datetime.utcnow()).total_seconds())
    zipmembers = zip.infolist() if membernames is None else [zipfile.getinfo(x) for x in membernames]
    for zipinfo in zipmembers:
        tarinfo = tf.TarInfo()
        tarinfo.name = zipinfo.filename
        tarinfo.size = zipinfo.file_size
        tarinfo.mtime = calendar.timegm(zipinfo.date_time) - timeshift
        tarinfo.mode = 0777 if zipinfo.filename.endswith("/") else 0666
        tarinfo.type = tf.DIRTYPE if zipinfo.filename.endswith("/") else tf.REGTYPE
        with zip.open(zipinfo.filename) as infile:
            tar.addfile(tarinfo, infile)
    if isinstance(zipfile, str):
        zip.close()
    tar.close()


def tar2zip(tarfile, zipname, membernames=None):
    if not zipname.endswith(".zip"):
        zipname += ".zip"
    zip = zf.ZipFile(zipname, "w")
    tar = tarfile if isinstance(tarfile, tf.TarFile) else tf.open(tarfile)
    timeshift = int((datetime.now() - datetime.utcnow()).total_seconds())
    tarmembers = tar.getmembers() if membernames is None else [tarfile.getmember(x) for x in membernames]
    for tarinfo in tarmembers:
        zipinfo = zf.ZipInfo()
        zipinfo.compress_type = zf.ZIP_DEFLATED
        zipinfo.filename = tarinfo.name+"/" if tarinfo.isdir() else tarinfo.name
        # zipinfo.file_size = tarinfo.size
        zipinfo.date_time = gmtime(tarinfo.mtime+timeshift)
        infile = tar.extractfile(tarinfo.name)
        zip.writestr(zipinfo, infile.read())
    if isinstance(tarfile, str):
        tar.close()
    zip.close()


# tarfile = "/geonfs01_vol1/ve39vem/ER01_SAR_IMP_1P_19911016T195412_19911016T195429_IPA_01315_0000.ESA.tar.gz"
# zipfile = "/geonfs01_vol1/ve39vem/ER01_SAR_IMP_1P_19911016T195412_19911016T195429_IPA_01315_0000.ESA.zip"


# tar2zip("E:/ER01_SAR_IMP_1P_19911016T195412_19911016T195429_IPA_01315_0000.ESA.tar.gz", "E:/test.zip")

# zip2tar("E:/ER01_SAR_IMP_1P_19911016T195412_19911016T195429_IPA_01315_0000.ESA.zip", "E:/test.tar.gz")

