process FILTER_FEATURES {
    label 'process_low'
    container 'quay.io/biocontainers/bioconductor-phyloseq:1.50.0--r44hdfd78af_0'
    conda 'bioconda::bioconductor-phyloseq=1.50.0'
    publishDir "${params.procdir}", mode: 'copy', pattern: '{asv_table_filtered,taxa_counts}.tsv'
    publishDir "${params.outdir}/05_filtered", mode: 'copy', pattern: '{metadata_filtered,filter_summary}.tsv'

    input:
    path features
    path taxonomy
    path metadata

    output:
    path 'asv_table_filtered.tsv',     emit: features
    path 'taxa_counts.tsv',            emit: taxa_counts
    path 'metadata_filtered.tsv',      emit: metadata
    path 'filter_summary.tsv',         emit: summary

    script:
    """
    05-filter_features.R --features ${features} --taxonomy ${taxonomy} --metadata ${metadata} \\
        --min_depth ${params.min_depth} --level ${params.tax_level}
    """
}
