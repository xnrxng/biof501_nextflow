process PREP_METADATA {
    label 'process_low'
    container 'quay.io/biocontainers/bioconductor-phyloseq:1.50.0--r44hdfd78af_0'
    conda 'bioconda::bioconductor-phyloseq=1.50.0'
    publishDir "${params.outdir}/02_metadata", mode: 'copy', pattern: '*.tsv'
    publishDir "${params.outdir}/01_prelim_eda", mode: 'copy', pattern: '*.png'

    input:
    path metadata

    output:
    path 'metadata_clean.tsv',        emit: metadata
    path 'variable_info.tsv',         emit: info
    path 'variable_distribution.png', emit: plot

    script:
    """
    02-prep_metadata.R --metadata ${metadata} --variable '${params.variable}' --variable_type ${params.variable_type}
    """
}
