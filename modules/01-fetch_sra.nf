process FETCH_SRA {
    tag "$sra_id"
    label 'sra_download'
    container 'quay.io/biocontainers/sra-tools:3.1.1--h4304569_2'
    publishDir "${params.rawdir}", mode: 'copy'

    input:
    val sra_id

    output:
    tuple val(sra_id), path("${sra_id}*.fastq.gz"), emit: reads
    path "${sra_id}.download.json",                 emit: manifest

    script:
    """
    export NCBI_SETTINGS="\$PWD/user-settings.mkfg"
    printf '/LIBS/GUID = "%s"\\n/libs/cloud/report_instance_identity = "false"\\n' "${sra_id}" > "\$NCBI_SETTINGS"

    prefetch ${sra_id} --max-size ${params.max_sra_size} --output-directory .
    fasterq-dump ${sra_id} --split-files --threads ${task.cpus} --temp . --outdir .

    if [ -f ${sra_id}_1.fastq ] && [ -f ${sra_id}_2.fastq ]; then rm -f ${sra_id}.fastq; fi

    if [ ${params.max_reads} -gt 0 ]; then
      for f in ${sra_id}*.fastq; do
        head -n \$((4 * ${params.max_reads})) "\$f" > "\$f.sub" && mv "\$f.sub" "\$f"
      done
    fi
    gzip -f ${sra_id}*.fastq

    01-write_manifest.sh ${sra_id} ${params.max_reads}
    """
}
