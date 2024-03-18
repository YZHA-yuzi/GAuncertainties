#!/bin/bash

#-------------------------------------
# qsub .sh file
#-------------------------------------

for((s=1;s<=3;++s))
do
	for ((i=1;i<=5;++i))
	do 
		Rscript Simulations.R "$s" "$i"
	done
done