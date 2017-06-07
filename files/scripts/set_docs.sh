#!/bin/bash
#Call with argument full to set full docs, user to set user, anything else does nothing

case "$1" in
  full)
    newval="YES"
    noop=0
  ;;
  user)
    newval="NO"
    noop=0
  ;;
  *)
    noop=1
    newval=""
  ;;
esac
echo $1
echo $newval

if [ $noop -eq 0 ]; then
  sed -i.bak "s/\(EXTRACT_PRIVATE *=\) *[A-Z]*/\1 $newval/" Doxyfile
  sed -i.bak "s/\(REFERENCED_BY_RELATION *=\) *[A-Z]*/\1 $newval/" Doxyfile
  sed -i.bak "s/\(SHOW_USED_FILES *=\) *[A-Z]*/\1 $newval/" Doxyfile
  sed -i.bak "s/\(SOURCE_BROWSER *=\) *[A-Z]*/\1 $newval/" Doxyfile
  sed -i.bak "s/\(VERBATIM_HEADERS *=\) *[A-Z]*/\1 $newval/" Doxyfile
fi