process TAXONOMY {
    label 'process_high'
    container 'quay.io/biocontainers/bioconductor-dada2:1.34.0--r44he5774e6_2'
    conda 'bioconda::bioconductor-dada2=1.34.0'
    publishDir "${params.procdir}", mode: 'copy'

    input:
    path seqs
    path db
    path species_db

    output:
    path 'taxonomy.tsv', emit: taxonomy

    script:
    """
    04-assign_taxonomy.R --seqs ${seqs} --db ${db} --species_db ${species_db} \\
        --threads ${task.cpus} --min_boot ${params.min_boot}
    """
}
