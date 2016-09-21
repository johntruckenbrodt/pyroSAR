##############################################################
# custom installation of GDAL and optional dependencies without administrator rights
# John Truckenbrodt 2016
##############################################################

import os
import re
import sys
import tarfile
import zipfile
from urllib2 import urlopen
from install.auxil import *

# the number of CPUs used for compilation
cores = 20

# root = "/homes4/geoinf/ve39vem/gdal"
# downloaddir = os.path.join(root, "originals")
# packagedir = os.path.join(root, "packages")
# installdir = os.path.join(root, "local")

root = "/geonfs02_vol2/software/gdal"
downloaddir = os.path.join(root, "originals")
packagedir = os.path.join(root, "packages")
installdir = "/geonfs02_vol2/software/.local"

# create installation directories in case they do not exist
for dir in [root, downloaddir, packagedir, installdir]:
    if not os.path.isdir(dir):
        os.makedirs(dir)

# optional dependency packages for GDAL
modules = ["geos", "proj", "tiff", "geotiff"]

# package version to be used
versions = {"gdal": "2.1.1",
            "geos": "3.5.0",
            "kml": "1.2.0",
            "htop": "1.0.2",
            "hdf4": "4.2.12",
            "hdf5": "1.8.17",
            "netcdf": "4.3.3.1",
            "proj": "4.9.2",
            "tiff": "4.0.6",
            "geotiff": "1.4.1",
            "expat": "2.2.0"}

# package URLs
# in case of hdf4 and hdf5 the most recent versions are going to be used
remotes = {"gdal": "http://download.osgeo.org/gdal/{0}/gdal-{0}.tar.gz".format(versions["gdal"]),
           "geos": "http://download.osgeo.org/geos/geos-{0}.tar.bz2".format(versions["geos"]),
           "proj": "http://download.osgeo.org/proj/proj-{0}.tar.gz".format(versions["proj"]),
           "tiff": "http://download.osgeo.org/libtiff/tiff-{0}.tar.gz".format(versions["tiff"]),
           "geotiff": "http://download.osgeo.org/geotiff/libgeotiff/libgeotiff-{0}.tar.gz".format(versions["geotiff"]),
           "hdf4": "https://www.hdfgroup.org/ftp/HDF/releases/HDF{0}/src/hdf-{0}.tar.gz".format(versions["hdf4"]),
           "hdf5": "https://www.hdfgroup.org/ftp/HDF5/releases/hdf5-{0}/src/hdf5-{0}.tar.gz".format(versions["hdf5"]),
           "kml": "https://libkml.googlecode.com/files/libkml-{0}.tar.gz".format(versions["kml"]),
           "netcdf": "ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-{0}.tar.gz".format(versions["netcdf"]),
           "expat": "https://sourceforge.net/projects/expat/files/expat/{0}/expat-{0}.tar.bz2".format(versions["expat"])}

##################################################################
print "downloading..."

# a dictionary
locals = {}

for module in modules+["gdal"]:
    remotename = remotes[module]
    localname = os.path.join(downloaddir, os.path.basename(remotename))
    locals[module] = localname
    try:
        if not os.path.isfile(localname):
            print module
            infile = urlopen(remotename)
            with open(localname, "wb") as outfile:
                outfile.write(infile.read())
            infile.close()
    except:
        print "...failed"

##################################################################
print "--------------------------------------------------"
print "unpacking..."

localdirs = {}

for module in modules+["gdal"]:
    tar = tarfile.open(locals[module])
    header = os.path.commonprefix(tar.getnames())
    localdir = os.path.join(packagedir, header)
    localdirs[module] = localdir
    if not os.path.isdir(localdir):
        print header.strip("/")
        tar.extractall(packagedir)
    tar.close()

##################################################################
print "--------------------------------------------------"
print "installing dependencies..."

# create installation directories for the packages
build_dirs = {}
for module in modules:
    build = os.path.join(localdirs[module], "build")
    build_dirs[module] = build
    if not os.path.isdir(build):
        os.makedirs(build)

# ./configure --prefix=/homes4/geoinf/ve39vem/test_gdal/packages/hdf-4.2.12/build --with-zlib --with-jpeg --enable-shared --disable-netcdf --disable-fortran

# define specific package configuration options
build_ops = {"geos": [],
             "proj": [],
             # "hdf4": ["--with-zlib",
             #          "--with-jpeg",
             #          "--enable-shared",
             #          "--disable-netcdf",
             #          "--disable-fortran"],
             # "hdf5": ["--enable-shared",
             #          "--enable-build-all",
             #          "--with-zlib",
             #          "--with-pthread",
             #          "--enable-cxx"],
             "tiff": ["--exec-prefix={}".format(build_dirs["tiff"])],
             "geotiff": ["--with-proj={}".format(build_dirs["proj"]),
                         "--with-libtiff={}".format(build_dirs["tiff"]),
                         "--with-zlib",
                         "--with-jpeg"],
             # "netcdf": ["--enable-netcdf-4",
             #            "--enable-shared",
             #            "--enable-utilities",
             #            "--enable-hdf4",
             #            "--enable-hdf4-file-tests",
             #            "--enable-v2"],
             "gdal": ["--without-python",
                      "--with-geos={}".format(os.path.join(build_dirs["geos"], "bin/geos-config")),
                      # "--with-hdf4={}".format(build_dirs["hdf4"]),
                      # "--with-hdf5={}".format(build_dirs["hdf5"]),
                      "--with-static-proj4={}".format(build_dirs["proj"]),
                      # "--with-netcdf={}".format(build_dirs["netcdf"]),
                      "--with-libtiff={}".format(build_dirs["tiff"]),
                      "--with-geotiff={}".format(build_dirs["geotiff"])]}

# configure, compile and install the packages
for module in modules:
    localdir = localdirs[module]
    build = os.path.join(localdir, "build")
    if not os.path.isdir(os.path.join(build, "bin")):
        print module
        # define the configuration command
        cmd = ["./configure", "-q", "--prefix={}".format(build)]+build_ops[module]
        print "...configuration"
        if module == "netcdf":
            # define extra environment variables
            env = os.environ.copy()
            env["H4DIR"] = build_dirs["hdf4"]
            env["H5DIR"] = build_dirs["hdf5"]
            env["CPPFLAGS"] = "-I{0} -I{1}".format(*[os.path.join(build_dirs[x], "include") for x in ["hdf4", "hdf5"]])
            env["LDFLAGS"] = "-L{0} -L{1}".format(*[os.path.join(build_dirs[x], "lib") for x in ["hdf4", "hdf5"]])
            execute(cmd, cwd=localdir, env=env)
        else:
            execute(cmd, cwd=localdir)
        print "...compilation"
        execute(["make", "-j{}".format(cores)], cwd=localdir)
        print "...installation"
        execute(["make", "install"], cwd=localdir)
#################################
print "--------------------------------------------------"
print "installing GDAL..."

# todo: adjust to Felix' system
# define extra environment variables
env = os.environ.copy()
env["LDFLAGS"] = "-L/usr/local/lib64"

if "LD_LIBRARY_PATH" not in env.keys():
    env["LD_LIBRARY_PATH"] = "/usr/local/lib64"
else:
    env["LD_LIBRARY_PATH"] = "/usr/local/lib64" + os.pathsep + env["LD_LIBRARY_PATH"]

localdir = localdirs["gdal"]
# build_gdal = os.path.join(localdir, "build")
build_gdal = installdir

# define the configuration command
cmd = ["./configure", "-q", "--prefix={}".format(build_gdal)]+build_ops["gdal"]

print "...configuration"
execute(cmd, cwd=localdir, env=env)

print "...compilation"
execute(["make", "-j{}".format(cores)], cwd=localdir)

print "...installation"
execute(["make", "install"], cwd=localdir)
#################################
print "--------------------------------------------------"
print "installing the GDAL Python binding..."

env["LD_LIBRARY_PATH"] = os.path.join(build_gdal, "lib") + os.pathsep + env["LD_LIBRARY_PATH"]
env["PATH"] = "{0}:{1}:{2}".format(os.path.join(build_gdal, "bin"),
                                   os.path.join(localdir, "swig/python/scripts"),
                                   env["PATH"])

print "...compilation"
execute(["make", "-j{}".format(cores)], cwd=os.path.join(localdir, "swig"), env=env)

packagesubdir = os.path.join(installdir, "lib/python2.7/site-packages")

if not os.path.isdir(packagesubdir):
    os.makedirs(packagesubdir)

if "PYTHONPATH" not in env.keys():
    env["PYTHONPATH"] = packagesubdir
else:
    env["PYTHONPATH"] = packagesubdir + os.pathsep + env["PYTHONPATH"]

cmd = [sys.executable, "setup.py", "install", "--prefix=" + installdir]

print "...installation"
execute(cmd, env=env, cwd=os.path.join(localdir, "swig/python"))

# create additional symlink
execute(["ln", "-s", "libgdal.so.20.1.0", "libgdal.so.1"], cwd=os.path.join(build_gdal, "lib"))
