#! /bin/bash

#BSUB -L /bin/sh
#BSUB -n 1
#BSUB -M 20000
#BSUB -W 23:00
#BSUB -J chromALLJob
#BSUB -oo outChromALL.out
#BSUB -eo chromALL.err

path=/home/khannm06/targ_projects/p053_KidneyBiobank/Data/eGFRcrea.SusztakLab.gwas
outPath=/home/khannm06

module load RHEL6-apps
module load python/anaconda-4.0.0

python byChrom.py $path $outPath