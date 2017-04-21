#!/bin/bash
#Extract any strings in classes called test_ containing "Set runtime_flag"
grep -r "Set runtime_flag" ./html/classtest__* | awk -F':' '{if(match($1, "classtest_")){ name = substr($1, match($1, "classtest")+5, length($1)-match($1, "class")-9); print "In class " name ":" $2}}'> ./files/tests_runtime_flags.txt
#Remove everything after the first para break
sed -i.bak 's|<\/p>.*|</p>|' ./files/tests_runtime_flags.txt
#Replace double underscores with single
sed -i.bak 's|__|_|g' ./files/tests_runtime_flags.txt
