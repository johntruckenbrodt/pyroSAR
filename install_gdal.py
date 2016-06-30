import os
import re
import sys
import tarfile
import zipfile
from urllib2 import urlopen
import subprocess as sp

# the number of CPUs used for compilation
cores = 10

root = "/homes4/geoinf/ve39vem/test_gdal"
downloaddir = os.path.join(root, "originals")
packagedir = os.path.join(root, "packages")
installdir = os.path.join(root, "local")

for dir in [root, downloaddir, packagedir, installdir]:
    if not os.path.isdir(dir):
        os.makedirs(dir)

# optional dependency packages for GDAL
modules = ["geos", "hdf4", "hdf5", "netcdf", "proj", "tiff", "geotiff"]

# package version to be used
versions = {"gdal": "2.1.0",
            "geos": "3.5.0",
            "kea": "1.4.4",
            "kml": "1.2.0",
            "htop": "1.0.2",
            "netcdf": "4.3.3.1",
            "proj": "4.9.2",
            "tiff": "4.0.6",
            "geotiff": "1.4.1"}

# package URLs
# in case of hdf4 and hdf5 the most recent versions are going to be used
remotes = {"gdal": "http://download.osgeo.org/gdal/{0}/gdal-{1}.tar.gz".format(versions["gdal"], versions["gdal"]),
           "geos": "http://download.osgeo.org/geos/geos-{0}.tar.bz2".format(versions["geos"]),
           "hdf4": "https://www.hdfgroup.org/ftp/HDF/HDF_Current/src/",
           "hdf5": "https://www.hdfgroup.org/ftp/HDF5/current/src/",
           "kea": "https://bitbucket.org/chchrsc/kealib/downloads/kealib-{0}.tar.gz".format(versions["kea"]),
           "kml": "https://libkml.googlecode.com/files/libkml-{0}.tar.gz".format(versions["kml"]),
           "netcdf": "ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-{0}.tar.gz".format(versions["netcdf"]),
           "proj": "http://download.osgeo.org/proj/proj-{0}.tar.gz".format(versions["proj"]),
           "tiff": "http://download.osgeo.org/libtiff/tiff-{0}.tar.gz".format(versions["tiff"]),
           "geotiff": "http://download.osgeo.org/geotiff/libgeotiff/libgeotiff-{0}.tar.gz".format(versions["geotiff"])}

##################################################################
print "downloading..."

# a dictionary
locals = {}

for module in modules+["gdal"]:
    if module in ["hdf4", "hdf5"]:
        remotename = os.path.join(remotes[module], re.findall("hdf[5]*-.[0-9.]+.tar.gz", urlopen(remotes[module]).read())[0])
    else:
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
print "compiling dependencies..."

build_dirs = {}

for module in modules:
    build = os.path.join(localdirs[module], "build")
    build_dirs[module] = build
    if not os.path.isdir(build):
        os.makedirs(build)

##################################################################

# define specific package configuration options
build_ops = {"geos": [],

             "proj": [],

             "hdf4": ["--with-zlib",
                      "--with-jpeg",
                      "--enable-shared",
                      "--disable-netcdf",
                      "--disable-fortran"],

             "hdf5": ["--enable-shared",
                      "--enable-build-all",
                      "--with-zlib",
                      "--with-pthread",
                      "--enable-cxx"],

             "tiff": ["--exec-prefix={}".format(build_dirs["tiff"])],

             "geotiff": ["--with-proj={}".format(build_dirs["proj"]),
                         "--with-libtiff={}".format(build_dirs["tiff"]),
                         "--with-zlib",
                         "--with-jpeg"],

             "netcdf": ["--enable-netcdf-4",
                        "--enable-shared",
                        "--enable-utilities",
                        "--enable-hdf4",
                        "--enable-hdf4-file-tests",
                        "--enable-v2"],

             "gdal": ["--without-python",
                      "--with-geos={}".format(os.path.join(build_dirs["geos"], "bin/geos-config")),
                      "--with-hdf4={}".format(build_dirs["hdf4"]),
                      "--with-hdf5={}".format(build_dirs["hdf5"]),
                      "--with-static-proj4={}".format(build_dirs["proj"]),
                      "--with-netcdf={}".format(build_dirs["netcdf"]),
                      "--with-libtiff={}".format(build_dirs["tiff"]),
                      "--with-geotiff={}".format(build_dirs["geotiff"])]}

# configure, compile and install the packages
for module in modules:
    localdir = localdirs[module]
    config = os.path.join(localdir, "configure")
    build = os.path.join(localdir, "build")
    if not os.path.isdir(os.path.join(build, "bin")):
        print module
        # define the configuration command
        cmd = [config, "-q", "--prefix={}".format(build)]+build_ops[module]

        # configuration
        if module == "netcdf":
            # define extra environment variables
            env = os.environ.copy()
            env["H4DIR"] = build_dirs["hdf4"]
            env["H5DIR"] = build_dirs["hdf5"]
            env["CPPFLAGS"] = "-I{0} -I{1}".format(*[os.path.join(build_dirs[x], "include") for x in ["hdf4", "hdf5"]])
            env["LDFLAGS"] = "-L{0} -L{1}".format(*[os.path.join(build_dirs[x], "lib") for x in ["hdf4", "hdf5"]])
            sp.check_call(cmd, cwd=localdir, env=env)
        else:
            sp.check_call(cmd, cwd=localdir)

        # compilation
        sp.check_call(["make", "--quiet", "-C", localdir, "-j{}".format(cores)])

        # installation
        sp.check_call(["make", "--quiet", "-C", localdir, "install"])
#################################
print "--------------------------------------------------"
print "installing GDAL..."

# define extra environment variables
env = os.environ.copy()
env["LDFLAGS"] = "-L/usr/local/lib64"
env["LD_LIBRARY_PATH"] = "/usr/local/lib64" + os.pathsep + env["LD_LIBRARY_PATH"]

localdir = localdirs["gdal"]
config = os.path.join(localdir, "configure")
build_gdal = os.path.join(localdir, "build")

# define the configuration command
cmd = [config, "-q", "--prefix={}".format(build_gdal)]+build_ops["gdal"]

# configuration, compilation, installation
sp.check_call(cmd, cwd=localdir, env=env)
sp.check_call(["make", "--quiet", "-C", localdir, "-j{}".format(cores)])
sp.check_call(["make", "--quiet", "-C", localdir, "install"])
#################################
# install the python binding

env["LD_LIBRARY_PATH"] = os.path.join(build_gdal, "lib") + os.pathsep + env["LD_LIBRARY_PATH"]
env["PATH"] = "{0}:{1}:{2}".format(os.path.join(build_gdal, "bin"),
                                   os.path.join(localdir, "swig/python/scripts"),
                                   env["PATH"])

sp.check_call(["make", "-C", os.path.join(localdir, "swig"), "-j{}".format(cores)], env=env)

env["PYTHONPATH"] = os.path.join(installdir, "lib/python2.7/site-packages") + os.pathsep + env["PYTHONPATH"]

cmd = [sys.executable, os.path.join(localdir, "swig/python/setup.py"), "install", "--prefix=" + installdir]

sp.check_call(cmd, env=env, cwd=os.path.join(localdir, "swig/python"))
