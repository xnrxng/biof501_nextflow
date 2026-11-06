#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Stage 1: fetch raw reads from NCBI SRA for every sample ID in the metadata.
 */

include { FETCH_SRA } from './modules/fetch_sra.nf'

workflow {
    // One SRR accession per row, taken from the "#SampleID" column of the metadata
    ch_ids = Channel
        .fromPath(params.metadata, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> row['#SampleID'].trim() }
        .filter { it }
        .unique()

    FETCH_SRA(ch_ids)
}
