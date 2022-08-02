#! /bin/bash

# This script submits jobs to the HPC in order to create 

#BSUB -L /bin/sh
#BSUB -n 1
#BSUB -M 64000
#BSUB -W 10:00
#BSUB -J gctaPrep
#BSUB -oo gctaPrep.out
#BSUB -eo gctaPrep.err

module load RHEL6-apps
module load python3/anaconda3-2.3.0

# Path to indexSNP CSVs. Pattern match to collect only 500kbp window.
csvPaths=(/home/khannm06/renal/parallel/chroms_hits/*500K*.csv)

# Create a directory of the SNP name based on the rsID
# and create a list of the SNPs in a 500kbp flanking region.
for csv in "${csvPaths[@]}"
do
    SNPList=($(cut -d ',' -f2 ${csv}))
    chrList=($(cut -d ',' -f3 ${csv}))

    for ((i = 1; i < ${#SNPList[@]}; ++i))
    do
        lociPath=/home/khannm06/renal/GCTA/Loci/${SNPList[i]}
        mkdir $lociPath
        sleep 1
        python3 gctaPrep.py ${SNPList[i]} ${chrList[i]} $lociPath
        bsub < file_GeneID.sh
    done
done