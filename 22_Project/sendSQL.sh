#! /bin/bash

#BSUB -L /bin/sh
#BSUB -n 6
#BSUB -M 99000
#BSUB -W 15:00
#BSUB -J sendSQLeX
#BSUB -oo sendSQLeX.out
#BSUB -eo sendSQLeX.err

module load RHEL6-apps
module load python3/anaconda3-2.3.0
module load ib
module load sqlite/3.34.1

python3 sendSQL.py