#!/usr/bin/env Rscript
# join ASVs + taxonomy + metadata and filter to remove uninformative features and samples
# aggregate to the requested taxonomic level for the abundance-based analyses

a <- commandArgs(trailingOnly = TRUE)
opt <- function(k, d = NULL) { i <- match(paste0("--", k), a); if (is.na(i)) d else a[i + 1] }
note <- function(case, observation, why, lesson) message("NOTE [", case, "]: ", observation, " (", why, "; ", lesson, ")")

#parameters
feat <- read.delim(opt("features"), check.names = FALSE, row.names = 1)
tax  <- read.delim(opt("taxonomy"), stringsAsFactors = FALSE, row.names = 1)
meta <- read.delim(opt("metadata"), stringsAsFactors = FALSE)
min_depth <- as.numeric(opt("min_depth", 1000))
level     <- opt("level", "family")
tax <- tax[rownames(feat), , drop = FALSE]

if ("Species" %in% colnames(tax))
  tax$Species <- ifelse(!is.na(tax$Species) & !is.na(tax$Genus), paste(tax$Genus, tax$Species), NA_character_)

steps <- data.frame(step = character(), n_asv = integer(), n_reads = numeric(), n_samples = integer())
log_step <- function(s) steps[nrow(steps) + 1, ] <<- list(s, sum(rowSums(feat) > 0), sum(feat), sum(colSums(feat) > 0))
log_step("denoised ASV table")

# samples present in both the reads and the metadata
common <- intersect(colnames(feat), meta$sample_id)
miss <- setdiff(meta$sample_id, colnames(feat))
if (length(miss)) note("samples without ASVs", paste(miss, collapse = ","), "no data after denoising",
                       "check read_tracking.tsv for where the reads were lost")
feat <- feat[, common, drop = FALSE]
log_step("samples matched to metadata")

#remove non-microbial featires
bad <- is.na(tax$Kingdom) | !tax$Kingdom %in% c("Bacteria", "Archaea") |
  (!is.na(tax$Family) & tax$Family == "Mitochondria") | (!is.na(tax$Order) & tax$Order == "Chloroplast")
feat <- feat[!bad, , drop = FALSE]; tax <- tax[!bad, , drop = FALSE]
log_step("removed unassigned / eukaryote / mitochondria / chloroplast")

#remove shallow samples
shallow <- colSums(feat) < min_depth
for (s in colnames(feat)[shallow]) note(paste("sample", s, "dropped"), sprintf("%d reads < min_depth %d", colSums(feat)[s], min_depth),
                                        "too shallow for rarefaction and stable diversity estimates", "raise max_reads or lower min_depth")
feat <- feat[, !shallow, drop = FALSE]
feat <- feat[rowSums(feat) > 0, , drop = FALSE]; tax <- tax[rownames(feat), , drop = FALSE]
log_step("removed samples below min_depth")
if (ncol(feat) < 4) stop("Fewer than 4 samples remain after filtering")

# aggregate to the analysis level
vague <- c("uncultured", "Incertae Sedis", "metagenome", "unidentified")
all_ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
lvl <- tolower(level)
if (lvl %in% tolower(all_ranks)) {
  target <- all_ranks[match(lvl, tolower(all_ranks))]
  ranks  <- rev(all_ranks[seq_len(match(target, all_ranks))])          
  ranks  <- ranks[ranks %in% colnames(tax)]
  lab <- rep(NA_character_, nrow(tax))
  for (r in ranks) {                                                    
    use <- is.na(lab) & !is.na(tax[[r]]) & !tax[[r]] %in% vague
    lab[use] <- if (r == ranks[1]) tax[[r]][use] else paste0(tax[[r]][use], "_unclassified")
  }
  lab[is.na(lab)] <- "Unassigned"
} else if (lvl == "asv") lab <- rownames(tax) else
  stop("level must be one of: ", paste(c(tolower(all_ranks), "asv"), collapse = ", "))

# aggregate to the requested taxonomic level
agg <- rowsum(as.matrix(feat), group = lab)
first <- !duplicated(lab)
ttax <- data.frame(taxon = lab[first], tax[first, , drop = FALSE], row.names = NULL)

write.table(data.frame(asv_id = rownames(feat), feat, check.names = FALSE), "features_filtered.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(taxon = rownames(agg), agg, check.names = FALSE), "taxa_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ttax, "taxa_taxonomy.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(meta[meta$sample_id %in% colnames(feat), ], "metadata_filtered.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(steps, "filter_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)