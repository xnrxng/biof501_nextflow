process TAXONOMY {
    label 'process_high'
    container 'quay.io/biocontainers/bioconductor-dada2:1.34.0--r44he5774e6_2'
    publishDir "${params.procdir}", mode: 'copy', pattern: 'taxonomy.tsv'
    publishDir "${params.outdir}/04_taxonomy", mode: 'copy', pattern: '{taxa_*.png,taxonomy_db_info.tsv}'

    input:
    path seqs
    path asv_table
    path db
    path species_db
    path metadata
    path info

    output:
    path 'taxonomy.tsv',         emit: taxonomy
    path 'taxonomy_db_info.tsv', emit: db_info
    path 'taxa_*.png',           emit: plots

    script:
    """
    04-assign_taxonomy.R --mode assign --seqs ${seqs} --asv_table ${asv_table} --db ${db} --species_db ${species_db} \\
        --level ${params.tax_level} --metadata ${metadata} --info ${info} --version_label '${params.silva_version}' --threads ${task.cpus} --min_boot ${params.min_boot} --seed ${params.seed}
    """
}


process TAXONOMY_PLOTS_FILTERED {
    label 'process_low'
    container 'quay.io/biocontainers/bioconductor-dada2:1.34.0--r44he5774e6_2'
    publishDir "${params.outdir}/04_taxonomy", mode: 'copy'

    input:
    path asv_table_filtered
    path metadata
    path info

    output:
    path 'taxa_*.png', emit: plots

    script:
    """
    04-assign_taxonomy.R --mode plot --asv_table ${asv_table_filtered} --level ${params.tax_level} --metadata ${metadata} --info ${info}
    """
}
