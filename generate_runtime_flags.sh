#!/bin/bash
grep -r "Set runtime_flag" ./html/class* | awk -F':' '{if(match($1, "classtest_")){ name = substr($1, match($1, "classtest")+5, length($1)-match($1, "class")-9); print "In class " name ":" $2}}'> ./files/tests_runtime_flags.txt
