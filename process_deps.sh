#!/bin/bash
  #dependencies .log is each entity on a new line, including continuation characters
  #We remove the lines we don't want, clean up the newlines and add blank lines after the last
  #It's ugly but it works
sed -i.bak -e '/Xcode.app/d;/boost/d;/SDF/d;/fftw3/d;/mpi/d' dependencies.log

sed -i.bak '/^ *\\/d' dependencies.log
  #remove the lines which are " \" or just "\"
sed -i.bak -e 's/$/\\/' dependencies.log
sed -i.bak -e $'s,[a-z/_]*\.o,\\\n&,' dependencies.log
  #append continuations at end

