import os
import re
import sys
import tarfile
import zipfile
from urllib2 import urlopen
import subprocess as sp

cores = 10

root = "/homes4/geoinf/ve39vem/test_gdal"
downloaddir = os.path.join(root, "originals")
packagedir = os.path.join(root, "packages")
installdir = os.path.join(root, "local")

for dir in [root, downloaddir, packagedir, installdir]:
    if not os.path.isdir(dir):
        os.makedirs(dir)

modules = ["gdal", "geos", "hdf4", "hdf5", "netcdf", "kml", "kea", "proj"]

versions = {"gdal": "2.1.0",
            "geos": "3.5.0",
            "kea": "1.4.4",
            "kml": "1.2.0",
            "htop": "1.0.2",
            "netcdf": "4.3.3.1",
            "rsgislib": "2.2.0",
            "proj": "4.9.2"}

remotes = {"gdal": "http://download.osgeo.org/gdal/{0}/gdal-{1}.tar.gz".format(versions["gdal"], versions["gdal"]),
           "geos": "http://download.osgeo.org/geos/geos-{0}.tar.bz2".format(versions["geos"]),
           "hdf4": "https://www.hdfgroup.org/ftp/HDF/HDF_Current/src/",
           "hdf5": "https://www.hdfgroup.org/ftp/HDF5/current/src/",
           "kea": "https://bitbucket.org/chchrsc/kealib/downloads/kealib-{0}.tar.gz".format(versions["kea"]),
           "kml": "https://libkml.googlecode.com/files/libkml-{0}.tar.gz".format(versions["kml"]),
           "netcdf": "ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-{0}.tar.gz".format(versions["netcdf"]),
           "proj": "http://download.osgeo.org/proj/proj-{0}.tar.gz".format(versions["proj"])}

print "downloading..."

locals = []

for module in modules:
    print module
    if re.search("hdf[45]", module):
        remotename = os.path.join(remotes[module], re.findall("hdf[5]*-.[0-9.]+.tar.gz", urlopen(remotes[module]).read())[0])
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
    localdir = os.path.join(packagedir, header)
    localdirs.append(localdir)
    if not os.path.isdir(localdir):
        print header.strip("/")
        tar.extractall(packagedir)
    tar.close()

print "--------------------------------------------------"
print "configuring..."

for localdir in localdirs:
    build = os.path.join(localdir, "build")
    if not os.path.isdir(build):
        os.makedirs(build)

#################################

localdir = [x for x in localdirs if re.search("geos", x)][0]
build_geos = os.path.join(localdir, "build")

if len(os.listdir(build_geos)) == 0:
    config = os.path.join(localdir, "configure")
    cmd = [config, "--prefix={}".format(build_geos)]
    sp.check_call(cmd, cwd=localdir)
    sp.check_call(["make", "-C", localdir, "-j{}".format(cores)])
    sp.check_call(["make", "-C", localdir, "install"])
#################################

localdir = [x for x in localdirs if re.search("hdf-", x)][0]
build_hdf4 = os.path.join(localdir, "build")

if len(os.listdir(build_hdf4)) == 0:
    config = os.path.join(localdir, "configure")
    cmd = [config, "--prefix={}".format(build_hdf4), "--with-zlib", "--with-jpeg", "--enable-shared", "--disable-netcdf", "--disable-fortran"]
    sp.check_call(cmd, cwd=localdir)
    sp.check_call(["make", "-C", localdir, "-j{}".format(cores)])
    sp.check_call(["make", "-C", localdir, "install"])
#################################

localdir = [x for x in localdirs if re.search("hdf5", x)][0]
build_hdf5 = os.path.join(localdir, "build")

if len(os.listdir(build_hdf5)) == 0:
    config = os.path.join(localdir, "configure")
    cmd = [config, "--prefix={}".format(build_hdf5), "--enable-shared", "--enable-build-all", "--with-zlib", "--with-pthread", "--enable-cxx"]
    sp.check_call(cmd, cwd=localdir)
    sp.check_call(["make", "-C", localdir, "-j{}".format(cores)])
    sp.check_call(["make", "-C", localdir, "install"])
#################################
# todo: check whether this is doing what it should
localdir = [x for x in localdirs if re.search("proj", x)][0]
build_proj = os.path.join(localdir, "build")

if len(os.listdir(build_proj)) == 0:
    config = os.path.join(localdir, "configure")
    cmd = [config, "--prefix={}".format(build_proj)]
    sp.check_call(cmd, cwd=localdir)
    sp.check_call(["make", "-C", localdir, "-j{}".format(cores)])
    sp.check_call(["make", "-C", localdir, "install"])
#################################

env = os.environ.copy()
env["H4DIR"] = build_hdf4
env["H5DIR"] = build_hdf5
env["CPPFLAGS"] = "-I{0} -I{1}".format(os.path.join(build_hdf4, "include"), os.path.join(build_hdf5, "include"))
env["LDFLAGS"] = "-L{0} -L{1}".format(os.path.join(build_hdf4, "lib"), os.path.join(build_hdf5, "lib"))

localdir = [x for x in localdirs if re.search("netcdf", x)][0]
build_netcdf = os.path.join(localdir, "build")

if len(os.listdir(build_netcdf)) == 0:
    config = os.path.join(localdir, "configure")
    cmd = [config, "--prefix={}".format(build_netcdf), "--enable-netcdf-4", "--enable-shared", "--enable-utilities", "--enable-hdf4", "--enable-hdf4-file-tests", "--enable-v2"]
    sp.check_call(cmd, cwd=localdir, env=env)
    sp.check_call(["make", "-C", localdir, "-j{}".format(cores)])
    sp.check_call(["make", "-C", localdir, "install"])
#################################

env["LD_LIBRARY_PATH"] = "/usr/local/lib64" + os.pathsep + env["LD_LIBRARY_PATH"]

localdir = [x for x in localdirs if re.search("gdal", x)][0]
config = os.path.join(localdir, "configure")
build_gdal = os.path.join(localdir, "build")


env["LDFLAGS"] = "-L/usr/local/lib64"

cmd = [config, "--prefix={}".format(build_gdal), "--without-python",
       "--with-geos={}".format(os.path.join(build_geos, "bin/geos-config")),
       "--with-hdf4={}".format(build_hdf4),
       "--with-hdf5={}".format(build_hdf5),
       "--with-static-proj4={}".format(build_proj),
       "--with-netcdf={}".format(build_netcdf)]

sp.check_call(cmd, cwd=localdir, env=env)
sp.check_call(["make", "-C", localdir, "-j{}".format(cores)])
sp.check_call(["make", "-C", localdir, "install"])
#################################

env["LD_LIBRARY_PATH"] = os.path.join(build_gdal, "lib") + os.pathsep + env["LD_LIBRARY_PATH"]
env["PATH"] = "{0}:{1}:{2}".format(os.path.join(build_gdal, "bin"), os.path.join(localdir, "swig/python/scripts"), env["PATH"])

sp.check_call(["make", "-C", os.path.join(localdir, "swig"), "-j{}".format(cores)], env=env)

env["PYTHONPATH"] = os.path.join(installdir, "lib/python2.7/site-packages") + os.pathsep + env["PYTHONPATH"]

sp.check_call([sys.executable, os.path.join(localdir, "swig/python/setup.py"), "install", "--prefix=" + installdir], env=env, cwd=os.path.join(localdir, "swig/python"))
