#!/usr/bin/env Rscript
# validate the metadata table and reduce it to (sample_id, value)

a <- commandArgs(trailingOnly = TRUE)
opt <- function(k, d = NULL) { i <- match(paste0("--", k), a); if (is.na(i)) d else a[i + 1] }

meta     <- read.csv(opt("metadata"), check.names = FALSE, stringsAsFactors = FALSE)
variable <- opt("variable")
vtype    <- opt("variable_type", "auto")

if (!variable %in% colnames(meta))
  stop("Column '", variable, "' not found in metadata. Available: ", paste(colnames(meta), collapse = ", "))

is_acc <- function(x) length(x) > 0 && all(grepl("^[SED]RR[0-9]+$", trimws(as.character(x))))
id_col <- Filter(function(cn) nzchar(cn) && is_acc(meta[[cn]]), colnames(meta))
if (!length(id_col)) stop("No column consisting of SRR/ERR/DRR accessions found in metadata")
id_col <- id_col[1]

df <- data.frame(sample_id = trimws(meta[[id_col]]), raw = meta[[variable]], stringsAsFactors = FALSE)
df$raw[trimws(as.character(df$raw)) %in% c("", "NA", "N/A", "na")] <- NA
n_in <- nrow(df)
df <- df[!is.na(df$raw) & !duplicated(df$sample_id), ]
if (nrow(df) < 2) stop("Need samples with a value for '", variable, "'; found ", nrow(df))

num <- suppressWarnings(as.numeric(as.character(df$raw)))
convertible <- !anyNA(num)
if (vtype == "auto") vtype <- if (convertible && length(unique(num)) >= 5) "numeric" else "categorical"
if (vtype == "numeric" && !convertible) stop("variable_type=numeric but '", variable, "' has non-numeric values")
df$value <- if (vtype == "numeric") num else as.character(df$raw)
if (vtype == "categorical") {
  if (length(unique(df$value)) < 2) stop("'", variable, "' has a single level; nothing to compare")
  cnt <- table(df$value)
  if (any(cnt < 2))
    stop("Each level of '", variable, "' needs at least 2 samples; too few in: ",
         paste(sprintf("%s (n=%d)", names(cnt)[cnt < 2], cnt[cnt < 2]), collapse = ", "))
} else if (nrow(df) < 4) stop("Need at least 4 samples with a value for numeric '", variable, "'; found ", nrow(df))

write.table(df[, c("sample_id", "value")], "metadata_clean.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
info <- data.frame(key = c("variable", "type", "id_column", "n_samples", "n_dropped_missing", "levels"),
                   value = c(variable, vtype, id_col, nrow(df), n_in - nrow(df),
                             if (vtype == "categorical") paste(sort(unique(df$value)), collapse = ";") else ""))
write.table(info, "variable_info.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("variable=%s type=%s id_column=%s samples=%d\n", variable, vtype, id_col, nrow(df)))
