##############################################################
# utility collection for compressed archives
# John Truckenbrodt 2016-2017
##############################################################

import os
import re
import sys
import fnmatch
import calendar
import zipfile as zf
import tarfile as tf
from time import gmtime
from datetime import datetime
from ancillary import finder
from StringIO import StringIO
import subprocess as sp


def addMember(archive, filename):
    """
    add file/directory to a zipped archive
    """
    handle_archive = zf.ZipFile(archive, 'a',  zf.ZIP_DEFLATED)
    base = os.path.basename(filename)
    if base not in handle_archive.namelist():
        handle_archive.write(filename, base)
    handle_archive.close()


def scan(archive, pattern, foldermode=1):
    """
    find files in a compressed archive
    output is a list with absolute directory names to files matching the defined pattern
    foldermodes:
    0: no folders
    1: folders included
    """
    files = []
    if isinstance(archive, str):
        archive = os.path.realpath(archive)
        handle = zf.ZipFile(archive, 'r') if zf.is_zipfile(archive) else tf.open(archive)
    else:
        handle = archive
    if isinstance(handle, zf.ZipFile):
        files = [os.path.join(handle.filename, x) for x in handle.namelist() if re.search(pattern, x.strip('/'))]
        if foldermode == 0:
            files = [x for x in files if not x.endswith('/')]
    elif isinstance(handle, tf.TarFile):
        files = [os.path.join(handle.name, x) for x in handle.getnames() if re.search(pattern, x)]
        if foldermode == 0:
            files = [x for x in files if not handle.getmember(x).isdir()]
    if isinstance(archive, str):
        handle.close()
    return files


def compress(src, dstname, compression='zip'):
    if compression not in ['zip', 'tar.gz']:
        raise IOError('compression type not supported; use either "zip" or "tar.gz"')
    if isinstance(src, list) or os.path.isfile(src):
        candidates = [src] if os.path.isfile(src) else src
        targets = [os.path.basename(x) for x in candidates]
    elif os.path.isdir(src):
        candidates = [src]+finder(src, ['*'], 1)
        targets = [os.path.relpath(x, os.path.join(src, '..')) for x in candidates]
    else:
        candidates, targets = [], []
    if len(candidates) == 0:
        raise IOError('nothing found to compress')
    if not dstname.endswith('.'+compression):
        dstname += '.'+compression
    if compression == 'zip':
        handle = zf.ZipFile(dstname, 'w', zf.ZIP_DEFLATED)
        for i in range(len(candidates)):
            handle.write(candidates[i], targets[i])
        handle.close()
    elif compression == 'tar.gz':
        handle = tf.open(dstname, 'w:gz')
        if os.path.isdir(src):
            handle.add(candidates[0], targets[0])
        else:
            for i in range(len(candidates)):
                handle.add(candidates[i], targets[i])
        handle.close()


def zip2tar(zipfile, tarname, membernames=None):
    """
    in-memory conversion of zip archives to tar.gz format
    adapted from 'http://stackoverflow.com/questions/18432565/python-convert-compressed-zip-to-uncompressed-tar'
    """
    zip = zipfile if isinstance(zipfile, zf.ZipFile) else zf.ZipFile(zipfile)
    if not tarname.endswith('.tar.gz'):
        tarname += '.tar.gz'
    tar = tf.open(tarname, 'w:gz')
    timeshift = int((datetime.now() - datetime.utcnow()).total_seconds())
    zipmembers = zip.infolist() if membernames is None else [zipfile.getinfo(x) for x in membernames]
    for zipinfo in zipmembers:
        tarinfo = tf.TarInfo()
        tarinfo.name = zipinfo.filename
        tarinfo.size = zipinfo.file_size
        tarinfo.mtime = calendar.timegm(zipinfo.date_time) - timeshift
        tarinfo.mode = 0777 if zipinfo.filename.endswith('/') else 0666
        tarinfo.type = tf.DIRTYPE if zipinfo.filename.endswith('/') else tf.REGTYPE
        with zip.open(zipinfo.filename) as infile:
            tar.addfile(tarinfo, infile)
    if isinstance(zipfile, str):
        zip.close()
    tar.close()


def tar2zip(tarfile, zipname, membernames=None):
    """
    in-memory conversion of tar.gz archives to zip format
    """
    if not zipname.endswith('.zip'):
        zipname += '.zip'
    zip = zf.ZipFile(zipname, 'w')
    tar = tarfile if isinstance(tarfile, tf.TarFile) else tf.open(tarfile)
    timeshift = int((datetime.now() - datetime.utcnow()).total_seconds())
    tarmembers = tar.getmembers() if membernames is None else [tarfile.getmember(x) for x in membernames]
    for tarinfo in tarmembers:
        zipinfo = zf.ZipInfo()
        zipinfo.compress_type = zf.ZIP_DEFLATED
        zipinfo.filename = tarinfo.name+'/' if tarinfo.isdir() else tarinfo.name
        # zipinfo.file_size = tarinfo.size
        zipinfo.date_time = gmtime(tarinfo.mtime+timeshift)
        infile = tar.extractfile(tarinfo.name)
        zip.writestr(zipinfo, infile.read())
    if isinstance(tarfile, str):
        tar.close()
    zip.close()


def deleteMembers(zipfile, matchlist, regex=False):
    """
    delete files/directories matching the defined patterns from a zipped archive
    """
    pattern = r'|'.join(matchlist if regex else [fnmatch.translate(x) for x in matchlist])
    members_rel = [os.path.relpath(x, zipfile) for x in scan(zipfile, pattern)]
    if len(members_rel) > 0:
        if re.search('linux', sys.platform):
            cmd = ['zip', '--quiet', '--delete', zipfile]+members_rel
            sp.check_call(cmd)
        else:
            oldzip = zf.ZipFile(zipfile, 'r')
            data = StringIO()
            newzip = zf.ZipFile(data, 'a', zf.ZIP_DEFLATED)
            if len(members_rel) > 0:
                for item in oldzip.namelist():
                    if os.path.normpath(item) not in members_rel:
                        member = oldzip.getinfo(item)
                        newzip.writestr(member, oldzip.read(item))
                oldzip.close()
                newzip.close()
                with open(zipfile, 'wb') as f:
                    f.write(data.getvalue())
