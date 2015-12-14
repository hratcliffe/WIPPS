#!/bin/bash

sed -i.bak -e '/Xcode.app/d;/boost/d;/SDF/d;/fftw3/d;/mpi/d' dependencies.log

sed -i.bak '/^ *\\/d' dependencies.log
  #remove the lines which are " \" or just "\"
sed -i.bak -e 's/$/\\/' dependencies.log
  #append continuations at end
#sed -i.bak 's,[a-z/_]*\.o,$OBJDIR\/&,' dependencies.log
  #prepend OBJDIR string
sed -i.bak -e $'s,[a-z/_]*\.o,\\\n&,' dependencies.log

