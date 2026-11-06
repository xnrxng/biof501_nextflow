process FETCH_SRA {
    tag "$sra_id"
    label 'sra_download'
    conda 'bioconda::sra-tools=3.1.1 conda-forge::pigz=2.8'
    publishDir "${params.outdir}/fastq", mode: 'copy'

    input:
    val sra_id

    output:
    tuple val(sra_id), path("${sra_id}*.fastq.gz"), emit: reads
    path "${sra_id}.download.json",                 emit: manifest

    script:
    """
    prefetch ${sra_id} --max-size ${params.max_sra_size} --output-directory .
    fasterq-dump ${sra_id} --split-files --threads ${task.cpus} --outdir .
    pigz -p ${task.cpus} ${sra_id}*.fastq

    # provenance: checksums + timestamp for every file fetched
    {
      echo '{"sra_id": "${sra_id}", "fetched_at": "'\$(date -u +%FT%TZ)'", "files": {'
      first=1
      for f in ${sra_id}*.fastq.gz; do
        [ \$first -eq 0 ] && echo ','
        echo "\\"\$f\\": \\"\$(sha256sum \$f | cut -d' ' -f1)\\""
        first=0
      done
      echo '}}'
    } > ${sra_id}.download.json
    """
}
