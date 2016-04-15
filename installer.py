
import os
import re
import tarfile
import zipfile
from urllib2 import urlopen

root = "/homes4/geoinf/ve39vem/source"

modules = ["gdal", "geos", "hdf4", "hdf5", "netcdf", "kml", "kea", "htop", "rsgislib"]

versions = {"gdal": "2.0.1",
            "geos": "3.5.0",
            "kea": "1.4.4",
            "kml": "1.2.0",
            "htop": "1.0.2",
            "netcdf": "4.3.3.1",
            "rsgislib": "2.2.0"}

remotes = {"gdal": "http://download.osgeo.org/gdal/{0}/gdal-{1}.tar.gz".format(versions["gdal"], versions["gdal"]),
           "geos": "http://download.osgeo.org/geos/geos-{0}.tar.bz2".format(versions["geos"]),
           "hdf4": "https://www.hdfgroup.org/ftp/HDF/HDF_Current/src/",
           "hdf5": "https://www.hdfgroup.org/ftp/HDF5/current/src/",
           "htop": "http://iweb.dl.sourceforge.net/project/htop/htop/1.0.2/htop-{0}.tar.gz".format(versions["htop"]),
           "kea": "https://bitbucket.org/chchrsc/kealib/downloads/kealib-{0}.tar.gz".format(versions["kea"]),
           "kml": "https://libkml.googlecode.com/files/libkml-{0}.tar.gz".format(versions["kml"]),
           "netcdf": "ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-{0}.tar.gz".format(versions["netcdf"]),
           "rsgislib": "https://bitbucket.org/petebunting/rsgislib/downloads/rsgislib-{0}.tar.gz".format(versions["rsgislib"])}

downloaddir = os.path.join(root, "originals")

print downloaddir

if not os.path.exists(downloaddir):
    os.makedirs(downloaddir)

print "start downloading..."

locals = []

for module in modules:
    print module

    if module == "hdf4":
        remotename = os.path.join(remotes[module], re.findall("hdf-.[0-9.]*.tar.gz", urlopen(remotes[module]).read())[0])
    elif module == "hdf5":
        remotename = os.path.join(remotes[module], re.findall("hdf5-.[0-9.]*.tar.gz", urlopen(remotes[module]).read())[0])
    else:
        remotename = remotes[module]
    localname = os.path.join(downloaddir, os.path.basename(remotename))
    locals.append(localname)
    try:
        if not os.path.isfile(localname):
            infile = urlopen(remotename)
            with open(localname, "wb") as outfile:
                outfile.write(infile.read())
            infile.close()
        else:
            print "...skipped"
    except:
        print "...failed"

print "--------------------------------------------------"
print "unpacking..."

localdirs = []

for i in range(len(locals)):

    tar = tarfile.open(locals[i])
    header = os.path.commonprefix(tar.getnames())
    localdir = os.path.join(root, header)
    localdirs.append(localdir)
    if not os.path.isdir(localdir):
        print header
        tar.extractall(root)
    tar.close()

print "--------------------------------------------------"
print "configuring..."

for localdir in localdirs:
    build = os.path.join(localdir, "build")
    if not os.path.isdir(build):
        os.makedirs(build)

# ./configure --prefix=/homes4/geoinf/ve39vem/source/geos-3.5.0/build
# ./configure --prefix=/homes4/geoinf/ve39vem/source/hdf-4.2.11/build --with-zlib --with-jpeg --enable-shared --disable-netcdf --disable-fortran
# ./configure --prefix=/homes4/geoinf/ve39vem/source/hdf5-1.8.16/build --enable-shared --enable-build-all --with-zlib --with-pthread --enable-cxx
# H4DIR="/homes4/geoinf/ve39vem/source/hdf-4.2.11/build" H5DIR="/homes4/geoinf/ve39vem/source/hdf5-1.8.16/build" CPPFLAGS="-I${H5DIR}/include -I${H4DIR}/include" LDFLAGS="-L${H5DIR}/lib -L${H4DIR}/lib" ./configure --prefix=/homes4/geoinf/ve39vem/source/netcdf-4.3.3.1/build --enable-netcdf-4 --enable-shared --enable-utilities --enable-hdf4 --enable-hdf4-file-tests --enable-v2

# export LD_LIBRARY_PATH=/usr/local/lib64/:$LD_LIBRARY_PATH
# ./configure LDFLAGS=-L/usr/local/lib64 --prefix=/homes4/geoinf/ve39vem/source/gdal-2.0.1/build/ --without-python --with-geos=/homes4/geoinf/ve39vem/source/geos-3.5.0/build/bin/geos-config --with-hdf4=/homes4/geoinf/ve39vem/source/hdf-4.2.11/build/ --with-hdf5=/homes4/geoinf/ve39vem/source/hdf5-1.8.16/build/ --with-netcdf=/homes4/geoinf/ve39vem/source/netcdf-4.3.3.1/build/

# cd /homes4/geoinf/ve39vem/source/gdal-2.0.1/swig

# export LD_LIBRARY_PATH=/homes4/geoinf/ve39vem/source/gdal-2.0.1/build/lib:$LD_LIBRARY_PATH
# export PATH=/homes4/geoinf/ve39vem/source/gdal-2.0.1/build/bin/:/homes4/geoinf/ve39vem/source/gdal-2.0.1/swig/python/scripts/:$PATH

# make -j7

# cd python

# export PYTHONPATH=/homes4/geoinf/ve39vem/local/lib/:PYTHONPATH

# python setup.py install --prefix=/homes4/geoinf/ve39vem/local/
