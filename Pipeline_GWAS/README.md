# Automating the GWAS target discovery pipeline
The purpose of this pipeline is to automate the process of filtering GWAS files and to output 22 files of GWAS index loci based on the provided parameters.

## Steps
1. Copy the "Pipeline_GWAS" into your home directory or wherever you would like to run the pipeline from.
2. Modify the nf.config file with your specific values.
3. Use the following command to run the pipeline:

        nextflow run gwas.nf -config nf.config

## Nextflow configuration
A master configuration [file](nf.config) is used for the pipeline
