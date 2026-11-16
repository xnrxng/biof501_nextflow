#!/usr/bin/env Rscript
# Differential abundance of taxa vs the variable of interest with DESeq2 (poscounts size factors).
suppressPackageStartupMessages({ library(DESeq2); library(ggplot2) })
a <- commandArgs(trailingOnly = TRUE)
opt <- function(k, d = NULL) { i <- match(paste0("--", k), a); if (is.na(i)) d else a[i + 1] }
notes <- data.frame(case = character(), observation = character(), why = character(), lesson = character())
note <- function(...) notes[nrow(notes) + 1, ] <<- list(...)
okabe <- c(blue = "#0072B2", orange = "#E69F00", grey = "grey65")

counts <- read.delim(opt("counts"), check.names = FALSE, row.names = 1)
meta   <- read.delim(opt("metadata"), stringsAsFactors = FALSE)
info   <- read.delim(opt("info"), stringsAsFactors = FALSE)
alpha_thr <- as.numeric(opt("alpha", 0.1)); min_prev <- as.numeric(opt("min_prev", 0.2))
binarize  <- opt("binarize", "none")
vname <- info$value[info$key == "variable"]; vtype <- info$value[info$key == "type"]
meta <- meta[match(colnames(counts), meta$sample_id), ]
n <- ncol(counts)

keep <- rowSums(counts > 0) >= max(2, ceiling(min_prev * n)) & rowSums(counts) >= 10
counts <- counts[keep, , drop = FALSE]
if (nrow(counts) < 2) stop("Fewer than 2 taxa pass the prevalence filter; lower --da_min_prev")

# define the model term
if (vtype == "numeric" && binarize == "median") {
  m <- as.numeric(meta$value)
  x <- factor(ifelse(m > median(m), "high", "low"), levels = c("low", "high"))
  if (nlevels(droplevels(x)) < 2) stop("median split produced a single group")
  design_desc <- sprintf("log2(high/low) where high = %s > median", vname); mode <- "two-group"
} else if (vtype == "numeric") {
  x <- as.numeric(scale(as.numeric(meta$value)))
  design_desc <- sprintf("log2 fold-change per +1 SD of %s", vname); mode <- "continuous"
} else {
  x <- factor(meta$value)
  mode <- if (nlevels(x) == 2) "two-group" else "multi-group"
  design_desc <- if (mode == "two-group") sprintf("log2(%s/%s)", levels(x)[2], levels(x)[1]) else "likelihood-ratio test across groups"
}
coldata <- data.frame(x = x, row.names = colnames(counts))
dds <- DESeqDataSetFromMatrix(round(as.matrix(counts)), coldata, ~x)
fit <- NULL
for (ft in c("parametric", "local", "mean")) {
  fit <- tryCatch(if (mode == "multi-group") DESeq(dds, sfType = "poscounts", test = "LRT", reduced = ~1, fitType = ft, quiet = TRUE)
                  else DESeq(dds, sfType = "poscounts", fitType = ft, quiet = TRUE), error = function(e) NULL)
  if (!is.null(fit)) { if (ft != "parametric") note("dispersion fit", paste("fell back to fitType =", ft),
                                                    "parametric trend did not converge on this small table", "treat dispersions with caution"); break }
}
if (is.null(fit)) stop("DESeq2 failed for every dispersion fit type")
res <- if (mode == "two-group") results(fit, contrast = c("x", levels(x)[2], levels(x)[1])) else results(fit)
res <- as.data.frame(res); res$taxon <- rownames(res)
res <- res[order(res$pvalue), c("taxon", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
write.table(res, "da_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

res$class <- ifelse(!is.na(res$padj) & res$padj < alpha_thr, sprintf("padj < %g", alpha_thr),
              ifelse(!is.na(res$pvalue) & res$pvalue < 0.05, "p < 0.05 only", "not significant"))
cls_cols <- setNames(c(okabe["orange"], okabe["blue"], okabe["grey"]), c(sprintf("padj < %g", alpha_thr), "p < 0.05 only", "not significant"))
r <- res[!is.na(res$pvalue), ]
v <- ggplot(r, aes(log2FoldChange, -log10(pvalue), colour = class)) + geom_point(size = 2.5, alpha = 0.85) +
  geom_text(data = head(r, 8), aes(label = taxon), size = 2.8, vjust = -0.8, check_overlap = TRUE, show.legend = FALSE) +
  scale_colour_manual(values = cls_cols, name = NULL) + geom_vline(xintercept = 0, linetype = 3) +
  labs(x = design_desc, y = expression(-log[10](p))) + theme_bw(base_size = 12)
ggsave("da_volcano.png", v, width = 7.5, height = 5.5, dpi = 150)

top <- head(r, 20); top$taxon <- factor(top$taxon, levels = rev(top$taxon))
b <- ggplot(top, aes(log2FoldChange, taxon, fill = class)) + geom_col() +
  scale_fill_manual(values = cls_cols, name = NULL) + labs(x = design_desc, y = NULL) + theme_bw(base_size = 12)
ggsave("da_top_features.png", b, width = 8, height = 6, dpi = 150)

n_sig <- sum(res$padj < alpha_thr, na.rm = TRUE)
write.table(data.frame(key = c("variable", "mode", "effect_size", "n_taxa_tested", "n_samples", "alpha", "n_significant_padj", "n_p_below_0.05"),
                       value = c(vname, mode, design_desc, nrow(counts), n, alpha_thr, n_sig, sum(res$pvalue < 0.05, na.rm = TRUE))),
            "da_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
if (n_sig == 0) note("no taxa pass FDR", sprintf("%d taxa tested in %d samples, 0 with padj < %g", nrow(counts), n, alpha_thr),
                     "small n gives low power; raw-p hits are exploratory only", "add samples before interpreting individual taxa")
write.table(notes, "notable_cases_da.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
