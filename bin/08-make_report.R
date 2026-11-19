#!/usr/bin/env Rscript
# puts all key results into a single HTML report (report.html) and writes provenance.json
suppressPackageStartupMessages(library(jsonlite))

# table and image helpers
rd  <- function(f) if (file.exists(f)) read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
tbl <- function(df, n = 30) {
  if (is.null(df)) return("")
  df <- head(df, n)
  df[] <- lapply(df, function(x) ifelse(is.na(x), "", if (is.numeric(x)) ifelse(x == round(x), format(x, scientific = FALSE, trim = TRUE), format(signif(x, 4), trim = TRUE)) else x))
  paste0("<div class='scroll'><table><tr>", paste0("<th>", colnames(df), "</th>", collapse = ""), "</tr>",
         paste0("<tr>", apply(df, 1, function(r) paste0("<td>", r, "</td>", collapse = "")), "</tr>", collapse = ""), "</table></div>")
}
img <- function(f) if (file.exists(f)) sprintf("<img src='data:image/png;base64,%s'>", base64_enc(readBin(f, "raw", file.size(f)))) else ""

# inputs
info <- rd("variable_info.tsv"); vname <- info$value[info$key == "variable"]; vtype <- info$value[info$key == "type"]
nlev <- if (vtype == "categorical") length(strsplit(info$value[info$key == "levels"], ";")[[1]]) else 0
params <- fromJSON("params.json")
meta <- rd("metadata_clean.tsv"); rar <- rd("rarefied_filtered_counts.tsv")
lvl <- tolower(params$tax_level)
at_lvl <- if (lvl == "asv") "at the ASV level (no aggregation)" else paste0("aggregated to the <b>", lvl, "</b> level, as requested")

# provenance.json
outs <- c("asv_table.tsv", "asv_seqs.fasta", "taxonomy.tsv", "asv_table_filtered.tsv", "rarefied_filtered_counts.tsv", "alpha_diversity.tsv", "pcoa_coords.tsv", "da_results.tsv")
outs <- outs[file.exists(outs)]
sha <- vapply(outs, function(f) strsplit(system2("sha256sum", f, stdout = TRUE), " ")[[1]][1], "")
write(toJSON(list(pipeline = params$pipeline_name, pipeline_version = params$pipeline_version, variable = vname, variable_type = vtype, parameters = params,
                  sra_ingest = lapply(list.files(".", pattern = "\\.download\\.json$"), fromJSON), output_sha256 = as.list(sha)),
              auto_unbox = TRUE, pretty = TRUE, null = "null"), "provenance.json")

# text that depends on the variable
order_txt <- if (vtype == "numeric") paste0("Samples are ordered by ", vname, ", lowest to highest from left to right.") else paste0("Samples are grouped by the levels of ", vname, ".")
alpha_test <- if (vtype == "numeric") "a Spearman correlation between the Shannon index and the variable" else if (nlev == 2) "a Wilcoxon rank-sum test between the two groups" else "a Kruskal-Wallis test across the groups followed by pairwise Wilcoxon tests (Benjamini-Hochberg adjusted)"
da_model <- if (vtype == "numeric" && params$da_binarize == "median") "samples were split at the median of the variable and the high group was compared with the low group" else if (vtype == "numeric") "the log2 fold-change is reported per +1 standard deviation of the variable" else "each level was compared with the reference level"
covariates <- if (nzchar(params$da_covariates)) paste0("The model was adjusted for: ", params$da_covariates, ".") else "No covariates were included in the model."

css <- "body{font-family:sans-serif;max-width:1000px;margin:2rem auto;padding:0 1rem;line-height:1.5} h1,h2{border-bottom:1px solid #ddd}
img{max-width:100%} table{border-collapse:collapse;font-size:.85rem} td,th{border:1px solid #ddd;padding:.2rem .5rem} th{background:#f0f0f0} .scroll{overflow-x:auto}"

h <- c(paste0("<!doctype html><html><head><meta charset='utf-8'><title>Differential abundance analysis report</title><style>", css, "</style></head><body>"),
"<h1>Differential abundance analysis on amplicon sequencing data Nextflow pipeline results</h1>",
"<p>This is a report generated from the results of the differential abundance analysis performed on amplicon sequencing data using the Nextflow pipeline. The report includes key findings, visualizations, and statistical summaries to provide insights into the microbial communities present in the samples analyzed.</p>",

"<h2>Metadata</h2>",
"<p>The selected variable of interest's relevant information is provided in this table. It contains the name of the variable, its type, the column with the sample accessions, the number of samples used, the number dropped for having no value, and the levels if the variable is categorical.</p>",
tbl(info),
paste0("<p>In this table, you can see the first ", min(6, nrow(meta)), " rows of the values of the variable of interest, out of ", nrow(meta), " samples.</p>"),
tbl(meta, 6),
"<p>For the distribution of the variable of interest, the following plot shows how the values are spread across all samples (a histogram if the variable is numeric, a bar plot of the number of samples per level if it is categorical). It is useful to check that groups are reasonably balanced and that values are not extremely skewed before comparing them.</p>",
img("variable_distribution.png"),

"<h2>DADA2</h2>",
"<p>The DADA2 algorithm was used to infer the amplicon sequence variants (ASVs) in each sample. It filters reads by quality, learns the sequencing error rates from the data, uses them to separate real sequences from sequencing errors, merges forward and reverse reads when the data are paired-end, and removes chimeras. This table shows how many reads each sample kept at each step: raw (input), after quality filtering (filtered), after error correction (denoised), after merging read pairs (merged, paired-end data only) and after removing chimeras (nonchim). A large drop at one step points to a problem there; for example, a poor merge rate means the reads overlap too little.</p>",
tbl(rd("read_tracking.tsv")),
"<p>This plot shows the error rates DADA2 learned, one panel per type of base substitution. The points are the observed error rates, the black line is the fitted model and the red line is the rate expected from the quality scores. Error rates should decrease as quality increases and the black line should follow the points; if it does not, the denoising is less reliable.</p>",
img("error_rates.png"),
paste0("<p>This plot shows the quality of the forward reads at each position, aggregated over all samples. The grey background is the frequency of each quality score, the green line is the mean quality and the orange lines are the median (solid) and quartiles (dashed). Quality usually drops towards the end of the reads, which is why reads are truncated (forward: ", params$trunc_len_f, ", reverse: ", params$trunc_len_r, "; auto means the first position where the median quality falls below 30, 0 means no truncation).</p>"),
img("quality_profiles_forward.png"),
if (file.exists("quality_profiles_reverse.png")) "<p>If applicable, here is the same for the reverse reads. These usually degrade earlier than the forward reads.</p>",
img("quality_profiles_reverse.png"),
"<p>A MultiQC report is also available under results/03_dada2_multiqc/multiqc. Keep in mind that amplicon data usually show a high number of duplicated sequences and uneven base composition. This is expected because all reads come from the same short region of the genome, so many of them are identical by design.</p>",

"<h2>Taxonomy classification</h2>",
paste0("<p>Taxonomy classification was performed with DADA2, which assigns each ASV a taxonomy with a naive Bayesian classifier trained on a reference database, keeping a rank only if its bootstrap confidence is at least ", params$min_boot, ". Species are then assigned by exact matching to a species reference. This table shows the reference databases used (names, version, checksums and sizes), the confidence threshold and when the classification was run, so it can be reproduced.</p>"),
tbl(rd("taxonomy_db_info.tsv")),
paste0("<p>This plot shows the composition of each sample before filtering, ", at_lvl, ". Each bar is a sample and each colour a taxon; the height is the number of reads. A taxon that could not be classified at this level is named after its lowest classified rank followed by _unclassified, or Unassigned if nothing was classified. ", order_txt, " It shows the sequencing depth of each sample and which taxa dominate it.</p>"),
img("taxa_absolute_unfiltered.png"),
"<p>And this plot shows the same but for relative abundances: every sample is scaled to 100%, so the composition of samples with different depths can be compared.</p>",
img("taxa_relative_unfiltered.png"),

"<h2>Filtering</h2>",
paste0("<p>Filtering was performed on the ASV table after classification. The steps are: keeping the samples that have a value for the variable of interest, removing ASVs that are not Bacteria or Archaea or that are mitochondria or chloroplast, removing samples with fewer than ", params$min_depth, " reads, and removing ASVs left without reads. This table shows the number of ASVs, reads and samples remaining after each step.</p>"),
tbl(rd("filter_summary.tsv")),
paste0("<p>This is how the filtering affected the abundance of the taxa, shown in the same way as before (", at_lvl, "). Comparing it with the plot before filtering shows what was removed.</p>"),
img("taxa_absolute_filtered.png"),
"<p>and for relative abundances.</p>",
img("taxa_relative_filtered.png"),

"<h2>Alpha and beta diversity</h2>",
paste0("<p>Alpha diversity measures the diversity inside each sample, and beta diversity measures how different samples are from each other. They are used to find out whether the diversity or the overall composition of the communities relates to the variable of interest. These were performed on the filtered table at the ASV level (not aggregated to the ", lvl, " level), after rarefaction.</p>"),
paste0("<p>Rarefaction randomly subsamples every sample to the same number of reads so that differences in sequencing depth do not create differences in diversity. The depth used was ", format(sum(rar[[2]]), big.mark = ","), " reads, and samples with fewer reads are dropped. This plot shows how many ASVs are observed as more reads are sampled, one curve per sample, with the depth used as a dashed line. If the curves are flat by that depth, the diversity was captured; if they are still rising, part of it was missed.</p>"),
img("rarefaction_curve.png"),
paste0("<p>Alpha diversity was calculated as Shannon's diversity index, which increases with the number of ASVs and with how evenly the reads are spread between them. This plot shows the index of each sample against the variable of interest, with the result of ", alpha_test, " written on it (ns = not significant, * p &lt; 0.05, ** p &lt; 0.01, *** p &lt; 0.001).</p>"),
img("alpha_diversity.png"),
"<p>Here is a table with the alpha diversity values for each sample.</p>",
tbl(rd("alpha_diversity.tsv")),
"<p>And this table shows the statistical test: which test was used, its statistic and p-value, and an adjusted p-value when there are several comparisons.</p>",
tbl(rd("alpha_stats.tsv")),
"<p>For beta diversity, the Bray-Curtis dissimilarity was used, which compares two samples by the abundance of the ASVs they contain (0 = identical, 1 = no ASVs in common). It was calculated between the rarefied samples and summarised with a principal coordinates analysis (PCoA). This plot shows the samples on the first two axes: samples that are closer together have more similar communities, and the axis labels give the percentage of variation each axis explains.</p>",
img("beta_diversity.png"),
paste0("<p>This table shows the PERMANOVA test, which asks whether the variable explains differences in community composition: R2 is the share of the variation explained and the p-value comes from 999 permutations.", if (vtype == "categorical") " For categorical variables a test of dispersion (PERMDISP) is also given: if it is significant, a significant PERMANOVA may reflect groups with different spread rather than different centres." else "", "</p>"),
tbl(rd("beta_stats.tsv")),
"<p>The Bray-Curtis matrix and the PCoA coordinates are available in results/06_diversity as bray_curtis_matrix.tsv and pcoa_coords.tsv respectively.</p>",

"<h2>Differential abundance analysis</h2>",
paste0("<p>Differential abundance analysis was performed with DESeq2, which fits a negative binomial model to the read counts of each taxon, corrects for differences in sequencing depth, and tests whether the abundance of each taxon changes with the variable of interest; p-values are adjusted for multiple testing (Benjamini-Hochberg). It was run on the filtered counts ", at_lvl, ", so each feature is one ", if (lvl == "asv") "ASV" else lvl, ", and only features present in at least ", 100 * as.numeric(params$da_min_prev), "% of the samples were tested. Here, ", da_model, ". ", covariates, " A feature is significant when its adjusted p-value is below ", params$da_alpha, " and its absolute log2 fold-change is above ", params$da_lfc, ".</p>"),
"<p>This plot shows the effect of every tested feature (log2 fold-change on the x axis) against its significance (-log10 p-value on the y axis, higher is more significant). Significant features are coloured (one colour for higher and one for lower abundance) and named; the dashed lines mark the fold-change threshold and grey points are not significant.</p>",
img("da_volcano.png"),
"<p>and this one shows the log2 fold-change of the 20 features with the smallest p-values, with the same colours.</p>",
img("da_top_features.png"),
"<p>The full DA results can be found under results/07_differential_abundance as da_results.tsv.</p>",

"<p><b>NOTE:</b> this report only condenses the key results of the analysis. For information on the pipeline, how to run it, and inputs, please refer to the README.</p>",
"</body></html>")
writeLines(unlist(h), "report.html")
