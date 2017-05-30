#!/bin/bash

#-n 8 -start 0 -end 4 -rows 300 -f ./files/l1/l1 -block ay
#Will call generate_FFTs with series of parameters

#These generally constant on one run-----------------------------------
program='../generate_ffts'
rows=512
#Number of times Power of 2 means better FFT performance
fil_pref=$1
#File path
n_files=10
#Number of files to get specified number of times
n_space=8
#Number of space blocks

#These will be varied in following ranges-----------------
allblocks=('ay' 'az' 'aby' 'abz')
#Which blocks to do
start_stride=10
#First start with file 0, then file 0+start_stride etc
end_num=111
#Last file to start with (i.e. uses files up to end_num+n_files)

start=0
n_blocks=${#allblocks[@]}

current_start=$start

while [ $current_start -lt $end_num ]; do
	block_ctr=0
	while [ $block_ctr -lt $n_blocks ]; do
		block=${allblocks[$block_ctr]}
		allstring="$program -n $n_space -start $current_start -end $(($current_start+$n_files)) -rows $rows -f $fil_pref -block $block"
		echo $allstring
		$allstring
		let block_ctr=block_ctr+1
	done
	let current_start=current_start+start_stride
done


