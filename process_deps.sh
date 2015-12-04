#!/bin/bash

sed -i.bak '/^ \\/d;/^\\/d' dependencies.log
  #remove the lines which are "  \" or just "\"
sed -i.bak -e 's/$/\\/' dependencies.log
  #append continuations at end
#sed -i.bak 's,[a-z/_]*\.o,$OBJDIR\/&,' dependencies.log
  #prepend OBJDIR string
sed -i.bak -e $'s,[a-z/_]*\.o,\\\n&,' dependencies.log

