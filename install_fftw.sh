#!/bin/bash

#Installs a local copy of FFTW of correct versions etc into path beginning $1

VERSION='3.3.4'
cd $1
curl ftp://ftp.fftw.org/pub/fftw/fftw-$VERSION.tar.gz -o fftw-$VERSION.tar.gz
tar -xf fftw-$VERSION.tar.gz
cd ./fftw-$VERSION/
echo "Got files"
#Install it here
DIR=`pwd`

DO_FLOAT=1
DO_DOUBLE=1
while [[ $# > 0 ]]; do
  if [ $1 == "--no-float" ]; then DO_FLOAT=0; fi
  if [ $1 == "--no-double" ]; then DO_DOUBLE=0; fi
  shift
done

#Float version
#Run the configure and make, and delete log output if we succeed
if [ $DO_FLOAT -eq 1 ]; then
  echo "Configuring for float precision"
  ./configure --enable-float --prefix=$DIR --with-pic --enable-type-prefix &>make_configf.log && rm make_configf.log 
  echo "Building"
  make &> make_makef.log && rm make_makef.log 
  make install &> make_makefi.log && rm make_makefi.log
  ls
  #Clean up
  make clean &> /dev/null
  if [ -e make_configf.log -o -e make_makef.log -o -e make_makefi.log ]; then 
    echo "Something went wrong installing float version, see "$1"fftw-$VERSION/make_*f.log" 
  fi
fi
if [ $DO_DOUBLE -eq 1 ]; then
  #Double version
  echo "Configuring for double precision"
  ./configure --prefix=$DIR --with-pic --enable-type-prefix &>make_configd.log && rm make_configd.log 
  echo "Building"
  make &> make_maked.log && rm make_maked.log 
  make install &> make_makedi.log && rm make_makedi.log
  if [ -e make_configd.log -o -e make_maked.log -o -e make_makedi.log ]; then 
    echo "Something went wrong installing double version, see "$1"fftw-$VERSION/make_*d.log" 
  fi
fi
if [ `ls $DIR/make_*[fd].log | wc -l` -gt 0 ]; then 
  exit 1
fi

