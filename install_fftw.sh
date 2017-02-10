#!/bin/bash

#Installs a local copy of FFTW of correct versions etc into path beginning $1

VERSION='3.3.4'
cd $1
curl ftp://ftp.fftw.org/pub/fftw/fftw-$VERSION.tar.gz -o fftw-$VERSION.tar.gz
tar -xf fftw-$VERSION.tar.gz
cd ./fftw-$VERSION/

#Install it here
DIR=`pwd`
#Float version
#Run the configure and make, and delete log output if we succeed
./configure --enable-float --prefix=$DIR --with-pic --enable-type-prefix &>make_configf.log && rm make_configf.log && make &> make_makef.log && rm make_makef.log && make install &> make_makefi.log && rm make_makefi.log
#Clean up
make clean &> /dev/null
if [ -e make_configf.log -o -e make_makef.log -o -e make_makefi.log ]; then 
  echo "Something went wrong installing float version, see "$1"fftw-$VERSION/make_*f.log" 
fi

#Double version
./configure --prefix=$DIR --with-pic --enable-type-prefix &>make_configd.log && rm make_configd.log && make &> make_maked.log && rm make_maked.log && make install &> make_makedi.log && rm make_makedi.log
if [ -e make_configd.log -o -e make_maked.log -o -e make_makedi.log ]; then 
  echo "Something went wrong installing double version, see "$1"fftw-$VERSION/make_*d.log" 
fi

if [ `ls $DIR/make_*[fd].log | wc -l` -gt 0 ]; then 
  exit 1
fi

