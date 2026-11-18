process DIFFERENTIAL {
    label 'process_low'
    container 'quay.io/biocontainers/bioconductor-deseq2:1.46.0--r44he5774e6_1'
    conda 'bioconda::bioconductor-deseq2=1.46.0'
    publishDir "${params.outdir}/07_differential_abundance", mode: 'copy'

    input:
    path counts
    path metadata
    path info
    path full_metadata

    output:
    path 'da_results.tsv', emit: results
    path '*.png',          emit: plots

    script:
    """
    07-differential.R --counts ${counts} --metadata ${metadata} --info ${info} \\
        --alpha ${params.da_alpha} --min_prev ${params.da_min_prev} --binarize ${params.da_binarize} \\
        --full_metadata ${full_metadata} --covariates '${params.da_covariates}' \\
        --reference '${params.da_reference}' --covariate_refs '${params.da_covariate_refs}'
    """
}
