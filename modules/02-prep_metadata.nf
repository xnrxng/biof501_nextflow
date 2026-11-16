process PREP_METADATA {
    label 'process_low'
    container 'quay.io/biocontainers/bioconductor-phyloseq:1.50.0--r44hdfd78af_0'
    conda 'bioconda::bioconductor-phyloseq=1.50.0'
    publishDir "${params.outdir}/02_metadata", mode: 'copy'

    input:
    path metadata

    output:
    path 'metadata_clean.tsv', emit: metadata
    path 'variable_info.tsv',  emit: info

    script:
    """
    02-prep_metadata.R --metadata ${metadata} --variable '${params.variable}' --variable_type ${params.variable_type}
    """
}
