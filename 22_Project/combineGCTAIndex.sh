#! /bin/bash

# Add locus name column
# for filename in /home/khannm06/users/khannam/Renal_Project/Data/GCTA-slct/*.cojo
# do
# done

# Combine GCTA Loci variants into one file
awk '(NR == 1) || (FNR > 1)' *.cojo > Index_Loci_GWAS_GCTA.tsv

