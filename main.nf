#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * 16S amplicon pipeline: SRA fetch -> DADA2 -> SILVA taxonomy -> filtering ->
 * alpha/beta diversity + DESeq2 + co-occurrence network -> HTML report,
 * all against one user-chosen metadata column (--variable).
 */

include { FETCH_SRA       } from './modules/01-fetch_sra.nf'
include { PREP_METADATA   } from './modules/02-prep_metadata.nf'
include { DADA2           } from './modules/03-dada2.nf'
include { FASTQC; MULTIQC } from './modules/04-qc.nf'
include { TAXONOMY        } from './modules/05-taxonomy.nf'
include { FILTER_FEATURES } from './modules/06-filter_features.nf'
include { DIVERSITY       } from './modules/07-diversity.nf'
include { DIFFERENTIAL    } from './modules/08-differential.nf'
include { NETWORK         } from './modules/09-network.nf'
include { REPORT          } from './modules/10-report.nf'

def clean(x) {
    (x ?: '').toString().replaceAll('"', '').trim()
}

workflow {
    main:
    // fail fast if the requested column is missing, before anything is downloaded
    def header = file(params.metadata).withReader { it.readLine() }
    def cols   = header.split(',').collect { it.replaceAll('"', '').trim() }
    if (!(params.variable in cols))
        error "Column '${params.variable}' not in ${params.metadata}. Columns: ${cols.join(', ')}"

    // one accession per row: the first column holding an SRR/ERR/DRR id; skip rows lacking a value.
    // splitCsv can leave the CSV quote characters on values, so strip them before matching.
    def vi    = cols.indexOf(params.variable)
    ch_ids = Channel
        .fromPath(params.metadata, checkIfExists: true)
        .splitCsv(skip: 1)
        .filter { r -> clean(r[vi]) && clean(r[vi]) != 'NA' }
        .map { r -> r.collect { x -> clean(x) }.find { x -> x ==~ /^[SED]RR[0-9]+$/ } }
        .filter { id -> id }
        .unique()
        .ifEmpty { error "No SRR/ERR/DRR accessions with a value for '${params.variable}' found in ${params.metadata}" }

    FETCH_SRA(ch_ids)
    FASTQC(FETCH_SRA.out.reads)
    PREP_METADATA(file(params.metadata, checkIfExists: true))

    DADA2(FETCH_SRA.out.reads.map { id, files -> files }.flatten().collect())
    MULTIQC(FASTQC.out.zip.collect(), DADA2.out.tracking)
    TAXONOMY(DADA2.out.seqs, file(params.silva_db, checkIfExists: true))

    FILTER_FEATURES(DADA2.out.table, TAXONOMY.out.taxonomy, PREP_METADATA.out.metadata)
    DIVERSITY(FILTER_FEATURES.out.features, FILTER_FEATURES.out.metadata, PREP_METADATA.out.info)
    DIFFERENTIAL(FILTER_FEATURES.out.taxa_counts, FILTER_FEATURES.out.metadata, PREP_METADATA.out.info)
    NETWORK(FILTER_FEATURES.out.taxa_counts, FILTER_FEATURES.out.taxa_taxonomy, FILTER_FEATURES.out.metadata, PREP_METADATA.out.info)

    def keys = ['variable', 'variable_type', 'max_reads', 'trunc_len_f', 'trunc_len_r', 'trim_left_f', 'trim_left_r', 'max_ee',
                'rawdir', 'procdir', 'silva_db', 'silva_version', 'min_boot', 'min_depth', 'tax_level', 'seed', 'da_binarize', 'da_alpha', 'da_min_prev',
                'net_min_prev', 'net_cor', 'net_p', 'net_padj_method', 'net_max_nodes']
    def pmap = params.subMap(keys) + [pipeline_name: workflow.manifest.name, pipeline_version: workflow.manifest.version]

    ch_report = Channel.empty()
        .mix(FETCH_SRA.out.manifest, PREP_METADATA.out.info,
             DADA2.out.table, DADA2.out.tracking, DADA2.out.plots, DADA2.out.notes,
             TAXONOMY.out.taxonomy, TAXONOMY.out.db_info,
             FILTER_FEATURES.out.features, FILTER_FEATURES.out.summary, FILTER_FEATURES.out.notes,
             DIVERSITY.out.alpha, DIVERSITY.out.alpha_stats, DIVERSITY.out.beta_stats, DIVERSITY.out.pcoa,
             DIVERSITY.out.summary, DIVERSITY.out.plots, DIVERSITY.out.notes,
             DIFFERENTIAL.out.results, DIFFERENTIAL.out.summary, DIFFERENTIAL.out.plots, DIFFERENTIAL.out.notes,
             NETWORK.out.nodes, NETWORK.out.edges, NETWORK.out.summary, NETWORK.out.plots, NETWORK.out.notes)
        .flatten()
        .collect()

    ch_params = Channel.value(groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(pmap)))
        .collectFile(name: 'params.json')

    REPORT(ch_report, ch_params)

    onComplete:
    log.info workflow.success ? "\nDone. Open: ${params.outdir}/10_report/report.html\n" : "\nPipeline failed; see .nextflow.log\n"
}
