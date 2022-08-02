nextflow.enable.dsl = 1

gwas_sum = Channel.value(file(params.gwas_sum_path))

process FILTER_GWAS {
    tag "FILTER_GWAS ${gwas_sum_name}"
    cpus 22
    memory '64 GB'
    time 1.hour
    module 'RHEL6-apps'
    module 'python3/anaconda3-2.3.0'

    input:
        val(gwas_sum_name) from gwas_sum
        file(gwas_sum_data) from gwas_sum

    output:
        path("*Loci.csv") into gwas_loci

    script:
        """
        python3 ${params.scripts}/multiRun.py ${gwas_sum_data} temp ${params.gwas_filter_window}
        """
}

process NEAREST_GENE {
    tag "NEAREST_GENE ${chr}"
    cpus 2
    memory '8 GB'
    time 1.hour
    module 'RHEL6-apps'
    module 'python3/anaconda3-2.3.0'

    input:
        tuple val(chr), file(chr_file) from gwas_loci

    output:
        path("*Loci_Genes.csv") into gwas_loci_genes

    script:
        """
        python3 mergeGenes.py $gwas $geneMap $lociPath
        """
}