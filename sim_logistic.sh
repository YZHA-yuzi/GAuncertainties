#!/bin/bash

#-------------------------------------
# qsub .sh file
#-------------------------------------

for((s=1;s<=3;++s))
do
		Rscript Simulations_logistic.R "$s"
done