import os
import sys

from osgeo import gdal


def needsVSICurl(filename):
    """
    Check if a given URL is a http, https or a ftp server.

    Parameters
    ----------
    filename : str
        URL

    Returns
    -------
    Server : bool
        True if the URL starts with http:// or https:// or ftp://.

    """
    return filename.startswith('http://') or filename.startswith('https://') or filename.startswith('ftp://')


def iszip(filename):
    """
    Check if a filename ends with .zip or .ZIP

    Parameters
    ----------
    filename : str
        URL or filename

    Returns
    -------
    Extension : bool
        True if a file ends with .zip or .ZIP.

    """
    return filename.endswith('.zip') or filename.endswith('.ZIP')


def istgz(filename):
    """
    Check if a filename ends with .tgz or tar.gz

    Parameters
    ----------
    filename : str
        URL or filename

    Returns
    -------
    Extension : bool
        True if a file ends with with .tgz or tar.gz.

    """
    return filename.endswith('.tgz') or filename.endswith('.TGZ') or \
           filename.endswith('.tar.gz') or filename.endswith('.TAR.GZ')


def display_file(fout, dirname, prefix, filename, longformat, check_open=False):
    """
    Try to display all files in a gdal env.

    Parameters
    ----------
    fout
    dirname
    prefix
    filename
    longformat
    check_open

    Returns
    -------

    """
    statBuf = None
    filename_displayed = prefix + filename

    if dirname.endswith('/'):
        dirname_with_slash = dirname
    else:
        dirname_with_slash = dirname + '/'

    version_num = int(gdal.VersionInfo('VERSION_NUM'))
    if longformat:
        if version_num >= 1900:
            statBuf = gdal.VSIStatL(dirname_with_slash + filename,
                                    gdal.VSI_STAT_EXISTS_FLAG | gdal.VSI_STAT_NATURE_FLAG | gdal.VSI_STAT_SIZE_FLAG)
    else:
        if version_num >= 1900:
            statBuf = gdal.VSIStatL(dirname_with_slash + filename,
                                    gdal.VSI_STAT_EXISTS_FLAG | gdal.VSI_STAT_NATURE_FLAG)

    if statBuf is None and check_open:
        if version_num >= 1900:
            f = None
        else:
            f = gdal.VSIFOpenL(dirname_with_slash + filename, "rb")
        if f is None:
            sys.stderr.write('Cannot open %s\n' % (dirname_with_slash + filename))
            return
        gdal.VSIFCloseL(f)

    if statBuf is not None and statBuf.IsDirectory() and not filename_displayed.endswith('/'):
        filename_displayed = filename_displayed + "/"

    if longformat and statBuf is not None:
        import time
        bdt = time.gmtime(statBuf.mtime)
        if statBuf.IsDirectory():
            permissions = "dr-xr-xr-x"
        else:
            permissions = "-r--r--r--"
        line = "%s  1 unknown unknown %12d %04d-%02d-%02d %02d:%02d %s\n" % \
               (permissions, statBuf.size, bdt.tm_year, bdt.tm_mon, bdt.tm_mday, bdt.tm_hour, bdt.tm_min,
                filename_displayed)
    else:
        line = filename_displayed + "\n"

    try:
        fout.write(line.encode('utf-8'))
    except:
        fout.write(line)


def readDir(fout, dirname, prefix, longformat, recurse, depth, recurseInZip, recurseInTGZ, first=False):
    """
    Try to read a whole dir in gdal env.

    Parameters
    ----------
    fout
    dirname
    prefix
    longformat
    recurse
    depth
    recurseInZip
    recurseInTGZ
    first

    Returns
    -------

    """
    if depth <= 0:
        return

    if needsVSICurl(dirname):
        dirname = '/vsicurl/' + dirname
        prefix = '/vsicurl/' + prefix

    if recurseInZip and iszip(dirname) and not dirname.startswith('/vsizip'):
        dirname = '/vsizip/' + dirname
        prefix = '/vsizip/' + prefix

    if recurseInTGZ and istgz(dirname) and not dirname.startswith('/vsitar'):
        dirname = '/vsitar/' + dirname
        prefix = '/vsitar/' + prefix

    lst = gdal.ReadDir(dirname)
    if lst is None:
        if first:
            original_dirname = dirname
            (dirname, filename) = os.path.split(dirname)
            if gdal.ReadDir(dirname) is None:
                sys.stderr.write('Cannot open {0}\n'.format(original_dirname))
                return
            if dirname == '':
                dirname = '.'
                prefix = ''
            else:
                prefix = dirname + '/'
            display_file(fout, dirname, prefix, filename, longformat, True)

    else:
        for filename in lst:
            if filename == '.' or filename == '..':
                continue

            display_file(fout, dirname, prefix, filename, longformat)

            if recurse:
                new_prefix = prefix + filename
                if not new_prefix.endswith('/'):
                    new_prefix += '/'
                readDir(fout, dirname + '/' + filename, new_prefix, \
                        longformat, recurse, depth - 1, recurseInZip, recurseInTGZ)


readDir(fout=sys.stdout, dirname='D:/Kartenmaterial/test_data/data', prefix='D:/Kartenmaterial/test_data/data/', longformat=True,
        recurse=True, depth=2, recurseInZip=True, recurseInTGZ=True,
        first=False)
