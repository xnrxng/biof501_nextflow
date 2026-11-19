#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * 16S amplicon pipeline: SRA fetch -> DADA2 -> SILVA taxonomy -> filtering ->
 * alpha/beta diversity + DESeq2 -> HTML report
 */

include { FETCH_SRA       } from './modules/01-fetch_sra.nf'
include { PREP_METADATA   } from './modules/02-prep_metadata.nf'
include { DADA2; FASTQC; MULTIQC } from './modules/03-dada2.nf'
include { TAXONOMY; TAXONOMY_PLOTS_FILTERED } from './modules/04-taxonomy.nf'
include { FILTER_FEATURES } from './modules/05-filter_features.nf'
include { DIVERSITY       } from './modules/06-diversity.nf'
include { DIFFERENTIAL    } from './modules/07-differential.nf'
include { REPORT          } from './modules/08-report.nf'

def clean(x) {
    (x ?: '').toString().replaceAll('"', '').trim()
}

workflow {
    main:
    // fail fast if the requested column is missing
    def header = file(params.metadata).withReader { it.readLine() }
    def cols   = header.split(',').collect { it.replaceAll('"', '').trim() }
    if (!(params.variable in cols))
        error "Column '${params.variable}' not in ${params.metadata}. Columns: ${cols.join(', ')}"

    def levels = ['phylum', 'class', 'order', 'family', 'genus', 'species', 'asv']
    if (!(params.tax_level.toString().toLowerCase() in levels))
        error "tax_level '${params.tax_level}' not valid. Choose one of: ${levels.join(', ')}"

    // one accession per row: the first column holding an SRR/ERR/DRR id
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
    TAXONOMY(DADA2.out.seqs, DADA2.out.table, file(params.silva_db, checkIfExists: true),
             file(params.silva_species_db, checkIfExists: true),
             PREP_METADATA.out.metadata, PREP_METADATA.out.info)

    FILTER_FEATURES(DADA2.out.table, TAXONOMY.out.taxonomy, PREP_METADATA.out.metadata)
    TAXONOMY_PLOTS_FILTERED(FILTER_FEATURES.out.features, PREP_METADATA.out.metadata, PREP_METADATA.out.info)
    DIVERSITY(FILTER_FEATURES.out.features, FILTER_FEATURES.out.metadata, PREP_METADATA.out.info)
    DIFFERENTIAL(FILTER_FEATURES.out.taxa_counts, FILTER_FEATURES.out.metadata, PREP_METADATA.out.info, file(params.metadata, checkIfExists: true))

    def keys = ['variable', 'variable_type', 'max_reads', 'trunc_len_f', 'trunc_len_r', 'trim_left_f', 'trim_left_r', 'max_ee',
                'rawdir', 'procdir', 'silva_db', 'silva_species_db', 'silva_version', 'min_boot', 'min_depth', 'tax_level', 'seed', 'rarefaction_depth', 'da_binarize', 'da_reference', 'da_covariates', 'da_covariate_refs', 'da_alpha', 'da_lfc', 'da_min_prev']
    def pmap = params.subMap(keys) + [pipeline_name: workflow.manifest.name, pipeline_version: workflow.manifest.version]

    ch_report = Channel.empty()
        .mix(FETCH_SRA.out.manifest, PREP_METADATA.out.info, PREP_METADATA.out.metadata, PREP_METADATA.out.plot,
             DADA2.out.table, DADA2.out.tracking, DADA2.out.plots,
             TAXONOMY.out.taxonomy, TAXONOMY.out.db_info, TAXONOMY.out.plots, TAXONOMY_PLOTS_FILTERED.out.plots,
             FILTER_FEATURES.out.features, FILTER_FEATURES.out.summary,
             DIVERSITY.out.alpha, DIVERSITY.out.alpha_stats, DIVERSITY.out.beta_stats, DIVERSITY.out.pcoa,
             DIVERSITY.out.rarefied, DIVERSITY.out.plots,
             DIFFERENTIAL.out.results, DIFFERENTIAL.out.plots)
        .flatten()
        .collect()

    ch_params = Channel.value(groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(pmap)))
        .collectFile(name: 'params.json')

    REPORT(ch_report, ch_params)

    onComplete:
    log.info workflow.success ? "\nDone. Open: ${params.outdir}/08_report/report.html\n" : "\nPipeline failed; see .nextflow.log\n"
}
