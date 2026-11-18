#!/usr/bin/env Rscript
# differential abundance of taxa vs the variable of interest with DESeq2
suppressPackageStartupMessages({ library(DESeq2); library(ggplot2) })
a <- commandArgs(trailingOnly = TRUE)
opt <- function(k, d = NULL) { i <- match(paste0("--", k), a); if (is.na(i)) d else a[i + 1] }
note <- function(case, observation, why, lesson) message("NOTE [", case, "]: ", observation, " (", why, "; ", lesson, ")")
split_arg <- function(s, sep = ",") if (is.null(s) || !nzchar(trimws(s))) character() else trimws(strsplit(s, sep, fixed = TRUE)[[1]])

counts <- read.delim(opt("counts"), check.names = FALSE, row.names = 1)
meta   <- read.delim(opt("metadata"), stringsAsFactors = FALSE)
info   <- read.delim(opt("info"), stringsAsFactors = FALSE)
alpha_thr <- as.numeric(opt("alpha", 0.1)); min_prev <- as.numeric(opt("min_prev", 0.2))
binarize  <- opt("binarize", "none")
vname <- info$value[info$key == "variable"]; vtype <- info$value[info$key == "type"]
meta <- meta[match(colnames(counts), meta$sample_id), ]

# covariates
covs     <- split_arg(opt("covariates"))
cov_refs <- split_arg(opt("covariate_refs"))                   
cov_df   <- data.frame(row.names = colnames(counts))
if (length(covs)) {
  orig <- read.csv(opt("full_metadata"), check.names = FALSE, stringsAsFactors = FALSE)
  is_acc <- function(x) length(x) > 0 && all(grepl("^[SED]RR[0-9]+$", trimws(as.character(x))))
  id_col <- Filter(function(cn) nzchar(cn) && is_acc(orig[[cn]]), colnames(orig))[1]
  if (is.na(id_col)) stop("No accession column found in the full metadata table", call. = FALSE)
  bad <- setdiff(covs, colnames(orig))
  if (length(bad)) stop("Covariate column(s) not in the metadata: ", paste(bad, collapse = ", "),
                        ". Available: ", paste(colnames(orig), collapse = ", "), call. = FALSE)
  if (vname %in% covs) stop("The variable of interest '", vname, "' cannot also be a covariate", call. = FALSE)
  refs <- setNames(sub("^[^=]*=", "", cov_refs), sub("=.*$", "", cov_refs))
  if (length(setdiff(names(refs), covs))) stop("--da_covariate_refs names a column that is not in --da_covariates: ",
                                              paste(setdiff(names(refs), covs), collapse = ", "), call. = FALSE)
  for (cv in covs) {
    raw <- orig[[cv]][match(colnames(counts), trimws(orig[[id_col]]))]
    raw[trimws(as.character(raw)) %in% c("", "NA", "N/A", "na")] <- NA
    num <- suppressWarnings(as.numeric(as.character(raw)))
    col <- paste0("cov_", make.names(cv))
    if (!anyNA(num[!is.na(raw)]) && length(unique(num[!is.na(num)])) > 2) {   
      cov_df[[col]] <- as.numeric(scale(num))
    } else {
      f <- factor(as.character(raw))
      if (cv %in% names(refs)) {
        if (!refs[[cv]] %in% levels(f)) stop("Reference '", refs[[cv]], "' is not a level of covariate '", cv, "'. Levels: ", paste(levels(f), collapse = ", "), call. = FALSE)
        f <- relevel(f, ref = refs[[cv]])
      }
      cov_df[[col]] <- f
    }
  }
  complete <- complete.cases(cov_df)
  if (any(!complete)) {
    note("samples dropped for missing covariates", paste(colnames(counts)[!complete], collapse = ", "), "DESeq2 needs every covariate for every sample",
         "fill in or drop the covariate")
    counts <- counts[, complete, drop = FALSE]; meta <- meta[complete, , drop = FALSE]; cov_df <- cov_df[complete, , drop = FALSE]
  }
  for (col in colnames(cov_df)) if (is.factor(cov_df[[col]])) cov_df[[col]] <- droplevels(cov_df[[col]])
  single <- colnames(cov_df)[vapply(cov_df, function(x) length(unique(x)) < 2, TRUE)]
  if (length(single)) stop("Covariate(s) with a single value among the analysed samples: ", paste(sub("^cov_", "", single), collapse = ", "), call. = FALSE)
}
n <- ncol(counts)

keep <- rowSums(counts > 0) >= max(2, ceiling(min_prev * n)) & rowSums(counts) >= 10
counts <- counts[keep, , drop = FALSE]
if (nrow(counts) < 2) stop("Fewer than 2 taxa pass the prevalence filter; lower --da_min_prev", call. = FALSE)

# term of interest
ref <- opt("reference", "")
if (vtype == "numeric" && binarize == "median") {
  m <- as.numeric(meta$value)
  x <- factor(ifelse(m > median(m), "high", "low"), levels = c("low", "high"))
  if (nlevels(droplevels(x)) < 2) stop("median split produced a single group", call. = FALSE)
  mode <- "two-group"; desc <- sprintf("high vs low (%s median split)", vname)
} else if (vtype == "numeric") {
  x <- as.numeric(scale(as.numeric(meta$value)))
  mode <- "continuous"; desc <- sprintf("per +1 SD of %s", vname)
} else {
  x <- droplevels(factor(meta$value))
  if (nzchar(ref)) {
    if (!ref %in% levels(x)) stop("Reference level '", ref, "' is not a level of '", vname, "'. Levels: ", paste(levels(x), collapse = ", "), call. = FALSE)
    x <- relevel(x, ref = ref)
  }
  mode <- if (nlevels(x) == 2) "two-group" else "multi-group"
}
if (nzchar(ref) && vtype == "numeric") note("reference ignored", paste("--da_reference =", ref), "the variable is numeric (or median-split)", "reference levels apply to categorical variables only")
coldata <- cbind(data.frame(x = x, row.names = colnames(counts)), cov_df)
design  <- as.formula(paste("~", paste(c(colnames(cov_df), "x"), collapse = " + ")))   # covariates first, term of interest last
if (n <= ncol(model.matrix(design, coldata))) stop("Too many model terms (", ncol(model.matrix(design, coldata)), ") for ", n, " samples; use fewer covariates", call. = FALSE)

dds <- DESeqDataSetFromMatrix(round(as.matrix(counts)), coldata, design)
fit <- NULL; err <- NULL
for (ft in c("parametric", "local", "mean")) {
  fit <- tryCatch(DESeq(dds, sfType = "poscounts", fitType = ft, quiet = TRUE), error = function(e) { err <<- conditionMessage(e); NULL })
  if (!is.null(fit)) { if (ft != "parametric") note("dispersion fit", paste("fell back to fitType =", ft), "parametric trend did not converge on this small table", "treat dispersions with caution"); break }
}
if (is.null(fit)) stop("DESeq2 failed for every dispersion fit type: ", err, call. = FALSE)

# results table
get_res <- function(contrast, label) {
  r <- as.data.frame(if (is.null(contrast)) results(fit, name = "x") else results(fit, contrast = contrast))
  r$taxon <- rownames(r); r$contrast <- label; r
}
res <- if (mode == "continuous") get_res(NULL, desc) else {
  lv <- levels(x)
  do.call(rbind, lapply(lv[-1], function(l) get_res(c("x", l, lv[1]), sprintf("%s vs %s", l, lv[1]))))
}
res <- res[order(res$contrast, res$pvalue), c("contrast", "taxon", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
write.table(res, "da_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
if (!any(res$padj < alpha_thr, na.rm = TRUE))
  note("no taxa pass FDR", sprintf("%d taxa tested in %d samples, 0 with padj < %g", nrow(counts), n, alpha_thr), "small n gives low power; raw-p hits are exploratory only", "add samples before interpreting individual taxa")

# plots
cls_levels <- c(sprintf("padj < %g", alpha_thr), "p < 0.05 only", "not significant")
res$class <- factor(ifelse(!is.na(res$padj) & res$padj < alpha_thr, cls_levels[1], ifelse(!is.na(res$pvalue) & res$pvalue < 0.05, cls_levels[2], cls_levels[3])),
                    levels = cls_levels)
r <- res[!is.na(res$pvalue), ]
multi <- length(unique(r$contrast)) > 1
adjust <- if (ncol(cov_df)) sprintf(" (adjusted for %s)", paste(covs, collapse = ", ")) else ""
v <- ggplot(r, aes(log2FoldChange, -log10(pvalue), colour = class)) + geom_point(size = 2.5, alpha = 0.85) +
  geom_text(data = do.call(rbind, lapply(split(r, r$contrast), head, 8)), aes(label = taxon), size = 2.8, vjust = -0.8, check_overlap = TRUE, show.legend = FALSE) +
  geom_vline(xintercept = 0, linetype = 3) + scale_x_continuous(expand = expansion(mult = 0.15)) + labs(x = "log2 fold-change", y = expression(-log[10](p)), colour = NULL,
                                                     title = sprintf("Differential abundance vs %s%s", vname, adjust)) +
  theme_classic(base_size = 12) + (if (multi) facet_wrap(~contrast, scales = "free") else labs(subtitle = unique(r$contrast)))
ggsave("da_volcano.png", v, width = if (multi) 11 else 7.5, height = 5.5, dpi = 150)

top <- do.call(rbind, lapply(split(r, r$contrast), head, 20))
top$key <- factor(paste(top$contrast, top$taxon, sep = "___"), levels = rev(paste(top$contrast, top$taxon, sep = "___")))
b <- ggplot(top, aes(log2FoldChange, key, fill = class)) + geom_col() +
  scale_y_discrete(labels = function(k) sub("^.*___", "", k)) +
  labs(x = "log2 fold-change", y = NULL, fill = NULL, title = sprintf("Top taxa by p-value, vs %s%s", vname, adjust)) +
  theme_classic(base_size = 12) + (if (multi) facet_wrap(~contrast, scales = "free_y") else labs(subtitle = unique(r$contrast)))
ggsave("da_top_features.png", b, width = if (multi) 12 else 8, height = 6, dpi = 150)
cat(sprintf("mode=%s samples=%d taxa_tested=%d covariates=%s\n", mode, n, nrow(counts), if (length(covs)) paste(covs, collapse = ",") else "none"))
