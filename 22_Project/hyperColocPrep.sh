#! /bin/bash

#BSUB -L /bin/sh
#BSUB -n 6
#BSUB -M 128000
#BSUB -W 45:00
#BSUB -J hyperColocPrep
#BSUB -oo hyperColocPrep.out
#BSUB -eo hyperColocPrep.err

module load RHEL6-apps
module load python3/anaconda3-2.3.0
module load ib
module load sqlite/3.34.1

dataFolderPath=/home/khannm06/users/khannam/Renal_Project/Data/HyperColoc
# mkdir $dataFolderPath
outPath=$dataFolderPath/HC_matrix

python3 hyperColocPrep.py