#!/usr/bin/env bash


# define a root directory for downloading packages
root= $HOME/test

# define a directory for download and unpacked packages
downloaddir=$root/originals
packagedir=$root/packages

# define the installation directory; This needs to be outside of the root directory so that the latter can be deleted in the end.
# In case installdir is set to a location outside of /usr/*, the following installation commands do not need to be run with a
# dministration rights (sudo)
#installdir=/usr/local
installdir=$HOME/local

# the version of GDAL and its dependencies
GDALVERSION=2.1.1

# these versions are not quite as important. If you use already installed them you might need to define their location
# for the configuration of GDAL
geos_version=3.5.0
proj_version=4.9.2

# define the number of threads for compilation
threads=2
########################################################################################################################
# setup environment variables and create directories

if [ -d "$root" ]; then
    if [  "$(ls -A $root)" ]; then
        echo "Error! root already exists. Please choose a fresh directory which can be deleted once finished" 1>&2
        #exit 64
    fi
fi

export PATH=$installdir/bin:$PATH
export LD_LIBRARY_PATH=$installdir/lib:$LD_LIBRARY_PATH

pythonlibdir=$installdir/lib/python2.7/site-packages
export PYTHONPATH=$pythonlibdir:$PYTHONPATH

for dir in $root $downloaddir $packagedir $pythonlibdir; do
    mkdir -p $dir
done
########################################################################################################################
# download GDAL and its dependencies

declare -a remotes=(
                "http://download.osgeo.org/gdal/$GDALVERSION/gdal-$GDALVERSION.tar.gz"
                "http://download.osgeo.org/geos/geos-$geos_version.tar.bz2"
                "http://download.osgeo.org/proj/proj-$proj_version.tar.gz"
                )

for package in "${remotes[@]}"; do
    wget $package -P $downloaddir
done
########################################################################################################################
# unpack downloaded archives

for package in $downloaddir/*tar.gz; do
    tar xfvz $package -C $packagedir
done
for package in $downloaddir/*tar.bz2; do
    tar xfvj $package -C $packagedir
done
########################################################################################################################
# install GEOS

cd $packagedir/geos*
./configure --prefix $installdir
make -j$threads
sudo make install
########################################################################################################################
# install PROJ

cd $packagedir/proj*
./configure --prefix $installdir
make -j$threads
sudo make install
########################################################################################################################
# install GDAL

# please check the output of configure to make sure that the GEOS and PROJ drivers are enabled
# otherwise you might need to define the locations of the packages

cd $packagedir/gdal*
./configure --without-python --prefix $installdir \
            --with-geos=$installdir/bin/geos-config \
            --with-static-proj4=$installdir \
            --with-libz=internal --with-pcraster=internal \
            --with-png=internal --with-pcidsk=internal \
            --with-libtiff=internal --with-geotiff=internal \
            --with-jpeg=internal --with-gif=internal \
            --with-qhull=internal --with-libjson-c=internal

make -j$threads
sudo make install
########################################################################################################################
# install GDAL Python binding
# this needs swig to be installed

cd $packagedir/gdal*/swig/python
make -j$threads
sudo python setup.py install --prefix=$installdir
########################################################################################################################
########################################################################################################################
# install pysqlite2 python package with static sqlite3 build
# this needs git to be installed

cd $packagedir
git clone https://github.com/ghaering/pysqlite.git
cd pysqlite

wget sqlite.org/2017/sqlite-amalgamation-3190300.zip

unzip sqlite-amalgamation-3190300.zip
cp sqlite-amalgamation-3190300/* .

sudo python setup.py build_static install --prefix=$installdir
########################################################################################################################
########################################################################################################################
# install spatialite

wget www.gaia-gis.it/gaia-sins/libspatialite-sources/libspatialite-amalgamation-2.4.0.tar.gz -P $downloaddir

cd $downloaddir
tar xfvz libspatialite-amalgamation-2.4.0.tar.gz -C $packagedir
cd $packagedir/libspatialite-amalgamation-2.4.0

./configure --with-geos-lib=$installdir/lib \
            --with-geos-include=$installdir/include \
            --with-proj-lib=$installdir/lib \
            --with-proj-include=$installdir/include \
            --prefix=$installdir

make -j$threads
sudo make install
########################################################################################################################
########################################################################################################################
# finishing the process

echo depending on your choice of installdir you might need to add the following lines to your .bashrc:
echo "export PATH=${installdir}/bin:$"PATH
echo "export LD_LIBRARY_PATH=${installdir}/lib:$"LD_LIBRARY_PATH
echo "export PYTHONPATH=${pythonlibdir}:$"PYTHONPATH
echo done

# deleting the root directory which is no longer needed
sudo rm -rf $root
