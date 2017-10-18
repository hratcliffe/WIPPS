#!/bin/bash

#Check existence of the necessary boost libs
#Specifically, we need boost_filesystem, which means also boost_system

wkdir=`dirname $0`"/"

echo "int main(){}" > $wkdir"test_boost.cpp"

#Really basic test here, one of these should work on most systems
#Don't bother identifying architecture, just try both

g++ -o $wkdir"test_boost.o" -L/usr/local/lib $wkdir"test_boost.cpp" -lboost_system -lboost_filesystem &> /dev/null
first=$?
g++ -o $wkdir"test_boost.o" -L/usr/lib $wkdir"test_boost.cpp" -lboost_system -lboost_filesystem &>/dev/null
second=$?
rm $wkdir"test_boost.o" $wkdir"test_boost.cpp" &> /dev/null

if ! (( first && second ))
 then
  exit 0
fi

exit 1
