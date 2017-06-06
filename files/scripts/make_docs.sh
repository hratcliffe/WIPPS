#!/bin/bash

#Better to strip the path part of the invokation ($0) and use that
if [ -n "$1" ]; then
  #If fed an arg, should be directory prefix for scripts
  scrpsdir=$1
else
  #Otherwise try here
  scrpsdir="./"
fi

if ! [[ `which Doxygen` ]]; then echo "Doxygen not found"; else printf "Options:\nType: ";egrep ' USE_FLOAT' ./Doxyfile >>/dev/null && echo "Float" || echo "Double";printf "Mode: ";egrep '[ \t]+RUN_TESTS_AND_EXIT' ./Doxyfile >>/dev/null && echo "Test" || echo "";printf "Docs: ";egrep -e "REFERENCED_BY_RELATION[ ]*\=[ ]*NO" ./Doxyfile >>/dev/null && echo "User" || echo "Full";
  $scrpsdir"generate_runtime_flags.sh"; doxygen Doxyfile 1> Doxy.log 2> Doxy.log.tmp; echo "Error Output:" >> Doxy.log; cat Doxy.log.tmp >> Doxy.log; rm Doxy.log.tmp; echo Processing Doxygen output...; $scrpsdir"redox.sh"; echo Running pdftex...; cd latex ; make &> ../docs.log; cd ..; echo "Docs built. See Doxy.log and docs.log for details";
fi