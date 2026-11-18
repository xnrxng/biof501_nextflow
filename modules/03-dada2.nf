process DADA2 {
    label 'process_high'
    container 'quay.io/biocontainers/bioconductor-dada2:1.34.0--r44he5774e6_2'
    conda 'bioconda::bioconductor-dada2=1.34.0'
    publishDir "${params.procdir}", mode: 'copy', pattern: '{asv_table.tsv,asv_seqs.fasta}'
    publishDir "${params.outdir}/03_dada2_multiqc", mode: 'copy', pattern: '{*.png,read_tracking.tsv}'

    input:
    path reads

    output:
    path 'asv_table.tsv',            emit: table
    path 'asv_seqs.fasta',           emit: seqs
    path 'read_tracking.tsv',        emit: tracking
    path '*.png',                    emit: plots

    script:
    """
    03-dada2_denoise.R --threads ${task.cpus} --max_ee ${params.max_ee} \\
        --trunc_f ${params.trunc_len_f} --trunc_r ${params.trunc_len_r} \\
        --trim_left_f ${params.trim_left_f} --trim_left_r ${params.trim_left_r}
    """
}

process FASTQC {
    tag "$sra_id"
    label 'process_low'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    conda 'bioconda::fastqc=0.12.1'

    input:
    tuple val(sra_id), path(reads)

    output:
    path '*_fastqc.zip', emit: zip

    script:
    """
    fastqc --threads ${task.cpus} --quiet ${reads}
    """
}

process MULTIQC {
    label 'process_low'
    container 'quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0'
    conda 'bioconda::multiqc=1.25.1'
    publishDir "${params.outdir}/03_dada2_multiqc/multiqc", mode: 'copy'

    input:
    path fastqc_zips
    path tracking

    output:
    path 'multiqc_report.html', emit: html
    path 'multiqc_report_data',  emit: data

    script:
    """
    # DADA2 read tracking as a MultiQC custom-content table
    {
      printf '# id: dada2_tracking\n# section_name: "DADA2 read tracking"\n'
      printf '# description: "Reads per sample after each DADA2 step"\n# plot_type: "table"\n'
      cat ${tracking}
    } > dada2_tracking_mqc.tsv

    multiqc . --force --filename multiqc_report.html
    """
}
