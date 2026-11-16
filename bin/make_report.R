#!/usr/bin/env Rscript
# Self-contained HTML report. Every number is read from the TSVs produced upstream;
# nothing is hard-coded. Also writes provenance.json (sha256 + data versions + timestamps).
suppressPackageStartupMessages(library(jsonlite))
a <- commandArgs(trailingOnly = TRUE)
opt <- function(k, d = NULL) { i <- match(paste0("--", k), a); if (is.na(i)) d else a[i + 1] }
rd  <- function(f) if (file.exists(f)) read.delim(f, check.names = FALSE, stringsAsFactors = FALSE) else NULL
esc <- function(x) { x <- gsub("&", "&amp;", as.character(x), fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE); gsub(">", "&gt;", x, fixed = TRUE) }
fmt <- function(x) if (is.numeric(x)) ifelse(is.na(x), "", ifelse(abs(x) < 1e-3 & x != 0, formatC(x, format = "e", digits = 2), format(round(x, 3), nsmall = 0, trim = TRUE))) else esc(x)
tbl <- function(df, max_rows = 25) {
  if (is.null(df) || !nrow(df)) return("<p><em>(none)</em></p>")
  df <- head(df, max_rows)
  hdr <- paste0("<th>", esc(colnames(df)), "</th>", collapse = "")
  df[] <- lapply(df, fmt)
  rows <- apply(df, 1, function(r) paste0("<tr>", paste0("<td>", r, "</td>", collapse = ""), "</tr>"))
  paste0("<div class='scroll'><table><thead><tr>", hdr, "</tr></thead><tbody>", paste(rows, collapse = ""), "</tbody></table></div>")
}
kv <- function(df, k) if (is.null(df)) NA else df$value[df$key == k]
fig <- function(file, what, how, shows) {
  if (!file.exists(file)) return(sprintf("<p><em>%s not produced</em></p>", esc(file)))
  b64 <- base64_enc(readBin(file, "raw", file.size(file)))
  sprintf("<figure><img src='data:image/png;base64,%s' alt='%s'><figcaption><b>What it is:</b> %s <b>How to read it:</b> %s <b>What it shows:</b> %s</figcaption></figure>",
          b64, esc(what), what, how, shows)
}
pfmt <- function(p) ifelse(p < 1e-3, formatC(p, format = "e", digits = 1), sprintf("%.3f", p))

info <- rd("variable_info.tsv"); vname <- kv(info, "variable"); vtype <- kv(info, "type")
params <- if (file.exists("params.json")) fromJSON("params.json") else list()
manifests <- lapply(list.files(".", pattern = "\\.download\\.json$"), fromJSON)
dbinfo <- rd("taxonomy_db_info.tsv")
track <- rd("read_tracking.tsv"); fsum <- rd("filter_summary.tsv")
alpha_stats <- rd("alpha_stats.tsv"); beta_stats <- rd("beta_stats.tsv"); dsum <- rd("diversity_summary.tsv")
da <- rd("da_results.tsv"); dasum <- rd("da_summary.tsv"); nsum <- rd("network_summary.tsv")
notes <- do.call(rbind, lapply(list.files(".", pattern = "^notable_cases_.*\\.tsv$"), rd))

alpha_txt <- if (!is.null(alpha_stats)) paste(sprintf("%s: %s = %.3g, p = %s", alpha_stats$metric, ifelse(alpha_stats$test == "Spearman correlation", "rho", "H"),
                                                      alpha_stats$estimate, pfmt(alpha_stats$p_value)), collapse = "; ") else ""
beta_txt <- if (!is.null(beta_stats)) paste(sprintf("%s: R2 = %.3f, p = %s", beta_stats$distance, beta_stats$R2, pfmt(beta_stats$p_value)), collapse = "; ") else ""
da_txt <- if (!is.null(dasum)) sprintf("%s taxa tested in %s samples; %s have padj < %s and %s have raw p < 0.05. Top taxon by p-value: %s (log2FC %.2f, p = %s).",
                                       kv(dasum, "n_taxa_tested"), kv(dasum, "n_samples"), kv(dasum, "n_significant_padj"), kv(dasum, "alpha"),
                                       kv(dasum, "n_p_below_0.05"), da$taxon[1], da$log2FoldChange[1], pfmt(da$pvalue[1])) else ""
net_txt <- if (!is.null(nsum)) sprintf("%s taxa tested; %s nodes carry edges; %s edges (%s positive, %s negative) in %s modules.",
                                       kv(nsum, "n_taxa_tested"), kv(nsum, "n_nodes_with_edges"), kv(nsum, "n_edges"), kv(nsum, "n_positive"),
                                       kv(nsum, "n_negative"), kv(nsum, "n_modules")) else ""

now <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
outs <- c("asv_table.tsv", "asv_seqs.fasta", "taxonomy.tsv", "features_filtered.tsv", "alpha_diversity.tsv", "pcoa_coords.tsv",
          "da_results.tsv", "network_edges.tsv", "network_nodes.tsv")
outs <- outs[file.exists(outs)]
sha <- vapply(outs, function(f) strsplit(system2("sha256sum", f, stdout = TRUE), " ")[[1]][1], "")
prov <- list(pipeline = params$pipeline_name, pipeline_version = params$pipeline_version, generated_at = now,
             variable = vname, variable_type = vtype, parameters = params,
             taxonomy_database = if (!is.null(dbinfo)) setNames(as.list(dbinfo$value), dbinfo$key) else NULL,
             sra_ingest = manifests, output_sha256 = as.list(sha))
write(toJSON(prov, auto_unbox = TRUE, pretty = TRUE, null = "null"), "provenance.json")

css <- "body{font-family:system-ui,sans-serif;max-width:1100px;margin:2rem auto;padding:0 1rem;color:#222;line-height:1.5}
h1,h2{border-bottom:1px solid #ddd;padding-bottom:.2rem} figure{margin:1.2rem 0} img{max-width:100%}
figcaption{font-size:.9rem;color:#444;background:#f6f6f6;padding:.6rem .8rem;border-radius:4px}
table{border-collapse:collapse;font-size:.85rem} td,th{border:1px solid #ddd;padding:.25rem .5rem;text-align:left} th{background:#f0f0f0}
.scroll{overflow-x:auto} code{background:#f0f0f0;padding:0 .25rem} .warn{background:#fff7e0;border-left:4px solid #E69F00;padding:.5rem .8rem}"
h <- c(sprintf("<!doctype html><html><head><meta charset='utf-8'><title>Microbiome report: %s</title><style>%s</style></head><body>", esc(vname), css),
       sprintf("<h1>16S microbiome report &mdash; %s (%s)</h1><p>Generated %s &middot; pipeline %s %s</p>", esc(vname), vtype, now, esc(params$pipeline_name), esc(params$pipeline_version)))

h <- c(h, "<h2>1. Reads and sample filtering</h2>",
       sprintf("<p>%s samples with a value for <code>%s</code> entered the analysis. Rarefaction depth for alpha diversity / Bray-Curtis: <b>%s reads</b>.</p>",
               kv(info, "n_samples"), esc(vname), kv(dsum, "rarefaction_depth")),
       "<h3>Read tracking through DADA2</h3>", tbl(track), "<h3>Feature filtering</h3>", tbl(fsum),
       fig("quality_profiles.png", "Aggregate per-cycle quality of the forward reads.", "x = cycle, y = Phred score; green line = median, orange = quartiles.",
           sprintf("where quality falls off, which sets the truncation length (trunc_len_f = %s).", esc(params$trunc_len_f))))
h <- c(h, "<h2>2. Alpha diversity</h2>",
       fig("alpha_diversity.png", sprintf("Alpha-diversity metrics against %s.", esc(vname)),
           if (vtype == "numeric") "one point per sample; line = linear fit with 95% band." else "boxes = IQR per group, points = samples.",
           alpha_txt), tbl(alpha_stats))
h <- c(h, "<h2>3. Beta diversity</h2>",
       fig("beta_diversity.png", sprintf("PCoA of sample dissimilarities coloured by %s.", esc(vname)),
           "samples closer together have more similar communities; axis labels give variance explained.", beta_txt), tbl(beta_stats))
h <- c(h, "<h2>4. Differential abundance (DESeq2)</h2>",
       fig("da_volcano.png", "Volcano plot of per-taxon effect size vs significance.", "right/left = higher/lower abundance with the variable; higher = smaller p; colour = significance class.", da_txt),
       fig("da_top_features.png", "Effect sizes of the 20 taxa with the smallest p-values.", "bar length = log2 fold-change; colour = significance class.", da_txt),
       tbl(da, 15))
h <- c(h, "<h2>5. Co-occurrence network</h2>",
       fig("network.png", "Spearman correlation network of CLR-transformed taxa (plus the variable as a square node when numeric).",
           "nodes = taxa coloured by phylum, size = degree; blue edges = positive and orange = negative correlation, width = |rho|.", net_txt))
h <- c(h, "<h2>6. Notable cases</h2>", tbl(notes))
h <- c(h, "<h2>7. Provenance</h2>",
       sprintf("<p>Taxonomy database: <code>%s</code> (%s), sha256 <code>%s</code>.</p>", esc(kv(dbinfo, "database")), esc(kv(dbinfo, "version_label")), esc(kv(dbinfo, "sha256"))),
       "<h3>SRA ingest</h3>",
       tbl(data.frame(sra_id = vapply(manifests, function(m) m$sra_id, ""), fetched_at = vapply(manifests, function(m) m$fetched_at, ""),
                      max_reads_per_file = vapply(manifests, function(m) as.character(m$subsample_max_reads), ""),
                      sha256_first_file = vapply(manifests, function(m) substr(m$files[[1]], 1, 16), ""))),
       "<h3>Output checksums</h3>", tbl(data.frame(file = names(sha), sha256 = unname(sha))),
       "<p>Full parameters and checksums: <code>provenance.json</code>.</p>", "</body></html>")
writeLines(unlist(h), "report.html")
