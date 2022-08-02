#! /bin/bash

#BSUB -L /bin/sh
#BSUB -n 1
#BSUB -M 6000
#BSUB -W 05:00
#BSUB -oo combSortCSV.out
#BSUB -eo combSortCSV.err

module load RHEL6-apps
module load python3/anaconda3-2.3.0

python3 combSortCSV.py
