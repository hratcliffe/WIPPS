#!/bin/bash
#Create test d files that can be compared to Padie data
../main -full_d -no_level_one -d_pp
mv test_d.dat ../test_d_all_p.dat
for N in -2 -1 0 1 2
do
  echo $N
  ../main -full_d -n $N -no_level_one -d_pp
  mv test_d.dat ../test_d_${N}_p.dat
done
