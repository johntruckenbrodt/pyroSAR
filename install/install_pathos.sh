#!/bin/bash

#define (and create) a main directory in which all packages will be downloaded and installed
maindir=$HOME/pathos_test

#define the python interpreter/version to be used
py=/usr/local/bin/python2.7

#the download repository
repository=dev.danse.us/packages

#define the versions of the packages to be downloaded
pathos_version=0.2a.dev-20130811
dill_version=0.2.4
pox_version=0.2.2

#create the directory into which the python packages are going to be installed and add it to the PYTHONPATH variable
mkdir -p $maindir/local/lib/python2.7/site-packages
export PYTHONPATH=$maindir/local/lib/python2.7/site-packages:$PYTHONPATH

#download the pathos package
wget $repository/pathos-$pathos_version.zip -P $maindir

#unpack the pathos package
unzip $maindir/pathos* -d $maindir

#download some extra dependencies
wget $repository/dill-$dill_version.zip -P $maindir/pathos*/external
wget $repository/pox-$pox_version.zip -P $maindir/pathos*/external

#unpack all the dependencies
for archive in $maindir/pathos*/external/*.zip; do 
	unzip $archive -d $maindir/pathos*/external; 
done

#install the dependency packages
for dir in $maindir/pathos*/external/*/; do
	(cd $dir
	$py setup.py build 
	$py setup.py install --prefix=$maindir/local)
done

#alternative installation using pip (was found to cause errors)
#$py /usr/local/bin/pip2.7 install --install-option="--prefix=${maindir}/local" dill pox
#$py /usr/local/bin/pip2.7 install --user --upgrade --install-option="--prefix${maindir}/local" dill pox

#install pathos
(cd $maindir/pathos*
$py setup.py build
$py setup.py install --prefix=$maindir/local)

#clean up
rm -rf $maindir/pathos*

echo please add the following line to your .bashrc:
echo "export PYTHONPATH=${maindir}/local/lib/python2.7/site-packages:$"PYTHONPATH
echo done
