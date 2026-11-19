process DIVERSITY {
    label 'process_low'
    container 'quay.io/biocontainers/bioconductor-phyloseq:1.50.0--r44hdfd78af_0'
    publishDir "${params.outdir}/06_diversity", mode: 'copy', pattern: '{*.png,alpha_*.tsv,beta_stats.tsv,bray_curtis_matrix.tsv,pcoa_coords.tsv}'
    publishDir "${params.procdir}", mode: 'copy', pattern: 'rarefied_filtered_counts.tsv'

    input:
    path features
    path metadata
    path info

    output:
    path 'alpha_diversity.tsv',     emit: alpha
    path 'alpha_stats.tsv',         emit: alpha_stats
    path 'beta_stats.tsv',          emit: beta_stats
    path 'pcoa_coords.tsv',         emit: pcoa
    path 'bray_curtis_matrix.tsv',  emit: bray
    path 'rarefied_filtered_counts.tsv',     emit: rarefied
    path '*.png',                   emit: plots

    script:
    """
    06-diversity.R --features ${features} --metadata ${metadata} --info ${info} --seed ${params.seed} --depth ${params.rarefaction_depth}
    """
}
