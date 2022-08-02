#! /bin/bash

#BSUB -L /bin/sh
#BSUB -n 2
#BSUB -M 96000
#BSUB -W 20:00
#BSUB -J runHyprC
#BSUB -oo runHyprC.out
#BSUB -eo runHyprC.err

module load ib
module load R/4.1.0

hyprDataPath=(/home/khannm06/users/khannam/Renal_Project/Scripts/hyperDataFull)
hyprData=($hyprDataPath/rs*)

mkdir $hyprDataPath/HyprCOLOCResults

for csv in "${hyprData[@]}"
do
    dataBaseName=$(basename "${csv}" .csv)
    hyprDataCutName=($hyprDataPath/HyprCOLOCResults/${dataBaseName}_eQTL_TRIM.csv)
    echo $hyprDataCutName
    cat ${csv} | cut -d, -f1-9 > ${hyprDataCutName}
    Rscript HyprCOLOC.R ${hyprDataCutName}
done

