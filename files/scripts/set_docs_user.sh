#!/bin/bash

sed -i.bak 's/\(EXTRACT_PRIVATE *=\) *[A-Z]*/\1 NO/' Doxyfile
sed -i.bak 's/\(REFERENCED_BY_RELATION *=\) *[A-Z]*/\1 NO/' Doxyfile
sed -i.bak 's/\(SHOW_USED_FILES *=\) *[A-Z]*/\1 NO/' Doxyfile
sed -i.bak 's/\(SOURCE_BROWSER *=\) *[A-Z]*/\1 NO/' Doxyfile
sed -i.bak 's/\(VERBATIM_HEADERS *=\) *[A-Z]*/\1 NO/' Doxyfile