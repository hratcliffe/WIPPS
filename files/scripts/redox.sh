#!/bin/sh

#  redox.sh
#  
#
#  Created by Heather Ratcliffe on 12/01/2016.
# Script to adjust doxygen Latex output to suit us. Later may generate entirely new refman.tex

sed -i.bak 's/\\documentclass\[twoside\]{book}/\\documentclass{article}/' ./latex/refman.tex
sed -i.bak 's/\\setcounter{tocdepth}{3}/\\setcounter{tocdepth}{2}/' ./latex/refman.tex
sed -i.bak '/\\backmatter/d' ./latex/refman.tex
sed -i.bak '/\\renewcommand{\\chaptermark}\[1\]{%\n\markboth{#1}{}%/d' ./latex/refman.tex



