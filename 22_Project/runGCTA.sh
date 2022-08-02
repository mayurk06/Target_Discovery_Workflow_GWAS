#! /bin/bash

#BSUB -L /bin/sh
#BSUB -n 1
#BSUB -M 16000
#BSUB -W 22:00
#BSUB -J runGCTA
#BSUB -oo runGCTA.out
#BSUB -eo runGCTA.err

# Path a directory of CSVs for each IndexSNP per chromosome
csvPaths=(/home/khannm06/renal/parallel/chroms_hits/*500K*.csv)

for csv in "${csvPaths[@]}"
do
    SNPList=($(cut -d ',' -f2 ${csv}))
    chrList=($(cut -d ',' -f3 ${csv}))

    for ((i = 1; i < ${#SNPList[@]}; ++i))
    do
        # lociPath=/home/khannm06/renal/GCTA/Loci/${SNPList[i]}
        lociPath=/home/khannm06/users/khannam/Renal_Project/Data/Loci/${SNPList[i]}
        cd $lociPath

        /hpc/grid/hgcb/workspace/projects/GeneticsTools_BinariesAndScripts/gcta_1.25.2/gcta64 \
        --bfile /hpc/grid/hgcb/workspace/projects/P002_reference_information/plink/1000G_Eur_Chr/EUR_AllChr.phase3_V5a_norm.${chrList[i]:3} \
        --chr ${chrList[i]:3} \
        --maf 0.005 \
        --cojo-collinear 0.3 \
        --cojo-p 5e-8 \
        --extract $lociPath/${SNPList[i]}.loci \
        --cojo-file /hpc/grid/wip_drm_targetsciences/projects/p053_KidneyBiobank/Data/eGFRcrea.SusztakLab.gcta \
        --cojo-slct \
        --out $lociPath/slct.eGRF.${SNPList[i]}
    done
done