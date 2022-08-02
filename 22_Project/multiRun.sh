#! /bin/bash

#BSUB -L /bin/sh
#BSUB -n 23
#BSUB -M 64000
#BSUB -W 15:00
#BSUB -J multiJob500K
#BSUB -oo outMulti500K.out
#BSUB -eo errMulti500K.err

# Provide paths without trailing "/"
path=/home/khannm06/targ_projects/p053_KidneyBiobank/Data/eGFRcrea.SusztakLab.gwas
outPath=/home/khannm06

module load RHEL6-apps
module load python3/anaconda3-2.3.0

python3 multiRun.py $path $outPath