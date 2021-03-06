#!/bin/bash

#Do all the necessary config, i.e.
#Check for FFTW install and local install if necessary
#Update FFTW path in Makefile
#Set NO_FFT if necessary
#Identify if boost filesystem is available
#Build sdf libraries for first time
#Build the test code and run it
#Wipe and build with supplied build args

DOCSTRING=$'Script to assist in getting the code running.\nIf FFTW is installed somewhere locally, supply the path with --fftw-path <installed path>. If not, use --install-fftw <path-to-install> (must exist!) to install locally, or install with your package manager (apt, homebrew, etc). A system-install will be found automatically.\nTo build without FFTW (some functions unavailable!) use --no-fft.\nTo build ONLY a double or floating point version supply --double or --float (not both!) respectively.\nThe final version will be built using whatever is supplied in --build-args (use quotes to supply more than one). If --double or --float is given, TYPE will be set automatically. Otherwise float will be used, but an explicit TYPE in --build-args will override this. For example, --build-args "TYPE=float".\nUse --help to see this documentation.'

#Path to other scripts etc needed
scrpsdir=./files/scripts/
#Makefile to write config stuff to
configfile=./config
#Remove an existing config file content
echo -n '' > $configfile

#Process args
FFTW_PATH=""
NO_FFT=0
INSTALL_FFTW=0
INSTALL_PATH=""
BUILD_ARGS=""
NO_DOUBLE=0
NO_FLOAT=0
SELECTED_PRECISION="float"
while [[ $# > 0 ]]; do
  case $1 in
    --help)
      echo "$DOCSTRING"
      exit 0
      ;;
    --fftw-path)
      shift
      FFTW_PATH=$1
      ;;
    --install-fftw)
      shift
      INSTALL_PATH=$1
      INSTALL_FFTW=1
      ;;
    --no-fft)
      NO_FFT=1
      ;;
    --build-args)
      shift
      BUILD_ARGS=$1
      ;;
    --double)
      NO_FLOAT=1
      SELECTED_PRECISION="double"
      ;;
    --float)
      NO_DOUBLE=1
      ;;
    *)
      echo "Unknown option "$1
      ;;
  esac
  shift
done

MAKEFILE=Makefile

#Add config line for NO_FFT if required
if [ $NO_FFT -eq 1 ]; then
  echo 'NO_FFT=1' >> $configfile
fi

#If we've been told to install we don't do anything with potential set FFTW_PATH
if [ $INSTALL_FFTW -eq 0 ]; then
  #If no path given, is there a system-wide install?
  if [ "$FFTW_PATH" = "" ]; then
    SYS=`uname -s`
    HAS_LIB=0
    #Look for headers and assume good if found
    case $SYS in
      Darwin)
        if [ -e /usr/local/include/fftw3.h ]; then HAS_LIB=1;fi
      ;;
      Linux)
       if [ -e /usr/include/fftw3.h ]; then HAS_LIB=1;fi
      ;;
      *)
        echo "Unknown system. Attempting to continue"
      ;;
    esac
    if [ $HAS_LIB -eq 0 ]; then echo "FFTW not found on expected path. Attempting to continue";else echo "Found FFTW libraries";fi
  else
    #Now if an FFTW path was given & we find ./include/fftw3.h and ./lib/libfftw3(f).a there we set the path in Makefile to match. 
    #If we don't find both libs we report on it
    #If we find neither we also report and suggest using --install-fftw <desired-path>
    if [ -e $FFTW_PATH/include/fftw3.h ] && ( [ -e $FFTW_PATH/lib/libfftw3.a ] || [ -e $FFTW_PATH/lib/libfftw3f.a ] ); then
      echo "Found FFTW files at "$FFTW_PATH
      if [ ! -e $FFTW_PATH/lib/libfftw3.a ]; then echo "No double library found, use only TYPE=float and float data or reinstall FFTW without --use-float"; fi
      if [ ! -e $FFTW_PATH/lib/libfftw3f.a ]; then echo "No float library found, use only TYPE=double and double data or reinstall FFTW with --use-float"; fi
      #Now set the path in configfile
      echo 'FFTW_PATH = '$FFTW_PATH >> $configfile
    else
      echo "Cannot find FFTW files at "$FFTW_PATH" Expecting ./include/fftw3.h and ./lib/libfftw3(f).a"
      echo "Check for correct path or try installing FFTW locally using ./configure --install-fftw <path to install>"
      exit 1
    fi
  fi
else
  #If asked to install FFTW we run the install script on correct path:
  if [ -w $INSTALL_PATH ]; then
    echo "Installing FFTW in "$INSTALL_PATH
    echo "This may take some time"
    extra_args=""
    if [ $NO_DOUBLE -eq 1 ]; then extra_args="--no-double";fi
    if [ $NO_FLOAT -eq 1 ]; then extra_args="--no-float";fi

    $scrpsdir"install_fftw.sh" $INSTALL_PATH $extra_args
    if [ $? -ne 0 ]; then
      echo "Error building FFTW!!!!!"
      exit 1
    fi
  else
    echo "Cannot write to "$INSTALL_PATH" Check permissions and retry"
    exit 1
  fi
fi

#By here we should have working FFTW install.

#Build sdf libraries explicitly
echo "Building SDF file libraries"
if [ -e ./SDF/C/lib/libsdfc.a ]; then echo "SDF already present"
else
  make ./SDF/C/lib/libsdfc.a >> /dev/null && echo "SDF built successfully"
  if [ $? -ne 0 ]; then
    #If we got an error, we can try cleaning
    echo "Error building SDF. Trying black magic"
    ORIG_WD=`pwd`
    cd ./SDF/C/ && make clean && make cleanall
    cd $ORIG_WD
    make ./SDF/C/lib/libsdfc.a >> /dev/null && echo "SDF built successfully"
    if [ $? -ne 0 ]; then
      echo "Cannot build SDF"
      exit 1
    fi
  fi
fi

echo "Checking for boost libraries"
$scrpsdir/check_boost.sh &> /dev/null
if [ $? -ne 0 ]; then
  echo "Boost not found"
  #Set config to not use Boost
  echo "NO_BOOST_LIBS=1" >> $configfile
else
  echo "Boost OK"
fi

#Now build the test code and run it
make clean
echo "Building test version of code"
make MODE=test TYPE=$SELECTED_PRECISION > /dev/null 2>&1
if [ $? -ne 0 ]; then
#In case of build error, try proper clean and rebuild
  make veryclean
  make MODE=test TYPE=$SELECTED_PRECISION > /dev/null 2>&1
fi
if [ -e ./main ]; then
  echo "Running test version of code" 
  ./main > tests.out.tmp
else
  echo "Error building test code. Try make MODE=test"
  exit 1
fi
if [ $? -ne 0 ]; then
  echo "Some tests failed!!!!!"
  cat tests.out.tmp
  rm tests.out.tmp
fi

echo "Building final version of code"
make clean
if [[ ! $BUILD_ARGS == *"MODE"* ]]; then
  BUILD_ARGS+="TYPE="$SELECTED_PRECISION
fi
make $BUILD_ARGS > /dev/null 2>&1 && echo "Building utilities" && make $BUILD_ARGS utils
if [ $? -eq 0 ]; then
  echo "Build complete without errors!"
fi

