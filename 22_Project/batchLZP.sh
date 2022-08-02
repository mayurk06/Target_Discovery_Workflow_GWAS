#! /bin/bash

#BSUB -L /bin/sh
#BSUB -n 1
#BSUB -M 28000
#BSUB -W 21:00
#BSUB -J lzpBatchY
#BSUB -oo lzpBatchY.out
#BSUB -eo lzpBatchY.err

ml eb/2017
ml intel
ml R
ml Python

# csvPaths=(/home/khannm06/renal/parallel/chroms_hits/*500K*.csv)
csvPaths=(/hpc/grid/wip_drm_targetsciences/users/khannam/Renal_Project/Data/chroms_hits/*500K*.csv)

for csv in "${csvPaths[@]}"
do
    SNPList=($(cut -d ',' -f2 ${csv}))
    # Skip the header line in csv
    for ((i = 1; i < ${#SNPList[@]}; ++i))
    do
        #lociPath=/home/khannm06/renal/GCTA/LociNew/${SNPList[i]}
        lociPath=/hpc/grid/wip_drm_targetsciences/users/khannam/Renal_Project/Data/LociNew/${SNPList[i]}
        mkdir $lociPath
        cd $lociPath
        echo ${SNPList[i]}
        
        /hpc/grid/hgcb/workspace/projects/GeneticsTools_BinariesAndScripts/locuszoom/bin/locuszoom \
        --metal /hpc/grid/wip_drm_targetsciences/projects/p053_KidneyBiobank/Data/eGFRcrea.SusztakLab.gwas \
        --markercol MarkerName \
        --pvalcol PVAL \
        --refsnp ${SNPList[i]} \
        --flank 500kb \
        --pop EUR \
        --build hg19 \
        --source 1000G_March2012 \
        --gwas-cat whole-cat_significant-only \
        --prefix eGFRcrea.SusztakLab.gwas.${SNPList[i]}
    done
done
