process REPORT {
    label 'process_low'
    container 'quay.io/biocontainers/bioconductor-phyloseq:1.50.0--r44hdfd78af_0'
    publishDir "${params.outdir}/08_report", mode: 'copy'

    input:
    path inputs
    path 'params.json'

    output:
    path 'report.html',     emit: html
    path 'provenance.json', emit: provenance

    script:
    """
    08-make_report.R
    """
}
