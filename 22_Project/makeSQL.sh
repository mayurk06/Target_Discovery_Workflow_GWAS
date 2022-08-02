#! /bin/bash

#BSUB -L /bin/sh
#BSUB -n 6
#BSUB -M 99000
#BSUB -W 12:00
#BSUB -J makeSQL2
#BSUB -oo makeSQL2.out
#BSUB -eo makeSQL2.err

module load RHEL6-apps
module load python3/anaconda3-2.3.0
module load ib
module load sqlite/3.34.1
# Glom QTL
python3 makeSQL.py