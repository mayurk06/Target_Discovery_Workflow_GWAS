#! /bin/bash

rsIDs=(/hpc/grid/wip_drm_targetsciences/users/khannam/Renal_Project/Data/Loci/rs*)
locusFile=/hpc/grid/wip_drm_targetsciences/users/khannam/Renal_Project/Data/IndexSNPs_Nearest_Gene/GWAS_500K_Loci.csv
locusList=($(cut -d ',' -f1 ${locusFile}))

# for x in ${rsIDs[@]}
for locus in ${locusList[@]}
do
    # Move gcta jma.cojo files
    rsID=$(echo ${locus} | sed 's|.*_||' )
    cp /hpc/grid/wip_drm_targetsciences/users/khannam/Renal_Project/Data/Loci/$rsID/*jma.cojo /hpc/grid/wip_drm_targetsciences/users/khannam/Renal_Project/Data/GCTA-slct/$locus.gcta.jma.cojo

    # Move locuszoom plots
    cp /hpc/grid/wip_drm_targetsciences/users/khannam/Renal_Project/Data/Loci/$rsID/*$rsID/*.pdf /hpc/grid/wip_drm_targetsciences/users/khannam/Renal_Project/Data/LZP/$locus.lzp.pdf
done
