#!/usr/bin/env Rscript
# taxonomy step. two modes (--mode):
#   assign: classify ASVs against a SILVA-style training set then plot per-sample composition of the unfiltered ASV table
#   plot: plot per-sample composition of the filtered ASV table 
suppressPackageStartupMessages(library(ggplot2))
a <- commandArgs(trailingOnly = TRUE)
opt <- function(k, d = NULL) { i <- match(paste0("--", k), a); if (is.na(i)) d else a[i + 1] }

# parameters and usage


usage <- paste0("\nUsage:\n",
  "  04-assign_taxonomy.R --mode assign --seqs <asv.fasta> --asv_table <asv_table.tsv> --db <train_set.fa.gz> --species_db <species_assignment.fa.gz> --level <rank|asv> --metadata <metadata_clean.tsv> --info <variable_info.tsv> [--version_label S] [--threads N] [--min_boot N]\n",
  "  04-assign_taxonomy.R --mode plot --asv_table <asv_table_filtered.tsv> --level <rank|asv> --metadata <metadata_clean.tsv> --info <variable_info.tsv>")
mode <- opt("mode", "assign")
if (!mode %in% c("assign", "plot")) stop("--mode must be 'assign' or 'plot'", usage, call. = FALSE)
need <- if (mode == "assign") c("seqs", "asv_table", "db", "species_db", "level", "metadata", "info") else c("asv_table", "level", "metadata", "info")
for (k in need) if (is.null(opt(k))) stop("Missing required argument --", k, usage, call. = FALSE)
for (k in setdiff(need, "level")) if (!file.exists(opt(k))) stop("File for --", k, " not found: ", opt(k), call. = FALSE)

level <- tolower(opt("level"))
meta  <- read.delim(opt("metadata"), stringsAsFactors = FALSE)          # sample_id, value
info  <- read.delim(opt("info"), stringsAsFactors = FALSE)              # key, value
vname <- info$value[info$key == "variable"]; vtype <- info$value[info$key == "type"]
all_ranks  <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
plot_ranks <- all_ranks[-1]                                              # Kingdom is not offered as a level
vague <- c("uncultured", "Incertae Sedis", "metagenome", "unidentified")

# composition plotting function

plot_composition <- function(tab, tax, stage) {
  sample_cols <- setdiff(colnames(tab), c("asv_id", all_ranks))
  counts <- as.matrix(tab[, sample_cols, drop = FALSE]); rownames(counts) <- tab$asv_id
  tax <- tax[match(tab$asv_id, tax$asv_id), , drop = FALSE]


  if ("Species" %in% colnames(tax))
    tax$Species <- ifelse(is.na(tax$Species) | is.na(tax$Genus), NA_character_,
                          ifelse(startsWith(tax$Species, paste0(tax$Genus, " ")), tax$Species, paste(tax$Genus, tax$Species)))


  if (level %in% tolower(plot_ranks)) {
    target <- plot_ranks[match(level, tolower(plot_ranks))]
    ranks  <- rev(plot_ranks[seq_len(match(target, plot_ranks))])
    ranks  <- ranks[ranks %in% colnames(tax)]
    lab <- rep(NA_character_, nrow(tax))
    for (r in ranks) {
      use <- is.na(lab) & !is.na(tax[[r]]) & !tax[[r]] %in% vague
      lab[use] <- if (r == ranks[1]) tax[[r]][use] else paste0(tax[[r]][use], "_unclassified")
    }
    lab[is.na(lab)] <- "Unassigned"
    level_name <- target
  } else if (level == "asv") { lab <- tab$asv_id; level_name <- "ASV"
  } else stop("level must be one of: ", paste(c(tolower(plot_ranks), "asv"), collapse = ", "), call. = FALSE)

  agg <- rowsum(counts, group = lab)
  agg <- agg[, colSums(agg) > 0, drop = FALSE]
  agg <- agg[order(rowSums(agg)), , drop = FALSE]                       
  taxa_levels <- rownames(agg)

  # sample order carries variable of interest: numeric = lowest to highest left to right; categorical = grouped by level
  m <- meta[match(colnames(agg), meta$sample_id), ]
  if (vtype == "numeric") {
    val <- suppressWarnings(as.numeric(m$value))
    ord <- order(val, na.last = TRUE)
    grp <- rep("", ncol(agg))
    slab <- ifelse(is.na(val), colnames(agg), sprintf("%s (%s)", colnames(agg), format(val, trim = TRUE, digits = 4)))
    subtitle <- sprintf("Samples ordered by %s, lowest to highest (left to right)", vname)
  } else {
    val <- as.character(m$value); val[is.na(val)] <- "no value"
    ord <- order(factor(val, levels = c(sort(setdiff(unique(val), "no value")), "no value")), colnames(agg))
    grp <- val
    slab <- colnames(agg)
    subtitle <- sprintf("Samples grouped by %s", vname)
  }
  sample_order <- colnames(agg)[ord]; names(slab) <- colnames(agg); names(grp) <- colnames(agg)

  long <- data.frame(sample = rep(colnames(agg), each = nrow(agg)), taxon = rep(rownames(agg), times = ncol(agg)),
                     count = as.vector(agg), stringsAsFactors = FALSE)
  long$taxon  <- factor(long$taxon, levels = taxa_levels)
  long$sample <- factor(long$sample, levels = sample_order)
  long$group  <- factor(grp[as.character(long$sample)], levels = unique(grp[sample_order]))
  long$rel    <- 100 * long$count / ave(long$count, long$sample, FUN = sum)

  n_tax  <- length(taxa_levels)
  n_col  <- max(1, ceiling(n_tax / 30))
  style <- function(ylab, ttl) list(
    theme_classic(base_size = 12),
    labs(x = NULL, y = ylab, fill = level_name, title = ttl, subtitle = subtitle),
    theme(plot.background = element_rect(fill = "white", colour = NA), panel.background = element_rect(fill = "white", colour = NA),
          axis.line = element_line(colour = "#0b0b0b", linewidth = 0.5), axis.ticks = element_line(colour = "#0b0b0b"),
          axis.text.x = element_text(angle = 45, hjust = 1, colour = "#52514e"), axis.text.y = element_text(colour = "#52514e"),
          legend.text = element_text(colour = "#0b0b0b", size = if (n_tax > 30) 7 else 9), legend.key.size = unit(if (n_tax > 30) 0.35 else 0.5, "cm"),
          plot.title = element_text(face = "bold", colour = "#0b0b0b"), plot.subtitle = element_text(colour = "#52514e"),
          strip.background = element_blank(), strip.text = element_text(face = "bold", colour = "#0b0b0b")),
    scale_x_discrete(labels = slab),
    guides(fill = guide_legend(ncol = n_col)),
    if (vtype != "numeric") facet_grid(. ~ group, scales = "free_x", space = "free_x"))
  bar <- function(y) geom_col(aes(y = .data[[y]]), width = 0.75, colour = "white", linewidth = 0.2)
  w <- max(9, min(24, 0.5 * ncol(agg) + 4)) + 3.2 * n_col

  p_abs <- ggplot(long, aes(sample, fill = taxon)) + bar("count") +
    scale_y_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.03))) +
    style("Reads", sprintf("Absolute abundance per sample, %s level (%s)", tolower(level_name), stage))
  p_rel <- ggplot(long, aes(sample, fill = taxon)) + bar("rel") +
    scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.01))) +
    style("Relative abundance", sprintf("Relative abundance per sample, %s level (%s)", tolower(level_name), stage))
  h <- max(5.5, 0.17 * min(n_tax, 30) + 2.5)
  ggsave(sprintf("taxa_absolute_%s.png", stage), p_abs, width = w, height = h, dpi = 150, limitsize = FALSE)
  ggsave(sprintf("taxa_relative_%s.png", stage), p_rel, width = w, height = h, dpi = 150, limitsize = FALSE)
  cat(sprintf("stage=%s level=%s taxa_shown=%d samples=%d\n", stage, level_name, n_tax, ncol(agg)))
}

tab <- read.delim(opt("asv_table"), check.names = FALSE, stringsAsFactors = FALSE)

if (mode == "plot") {
  plot_composition(tab, tab[, c("asv_id", intersect(all_ranks, colnames(tab))), drop = FALSE], "filtered")
  quit(save = "no")
}

# assign taxonomy to ASVs and plot composition of the unfiltered ASV table
suppressPackageStartupMessages({ library(dada2); library(Biostrings) })
seqs <- readDNAStringSet(opt("seqs"))
if (!length(seqs)) stop("No sequences in ", opt("seqs"), call. = FALSE)
tax  <- assignTaxonomy(as.character(seqs), opt("db"), multithread = as.integer(opt("threads", 1)),
                       minBoot = as.integer(opt("min_boot", 50)), tryRC = TRUE, verbose = FALSE)

tax  <- addSpecies(tax, opt("species_db"), tryRC = TRUE, verbose = FALSE)
out  <- data.frame(asv_id = names(seqs), tax[as.character(seqs), , drop = FALSE], row.names = NULL)
write.table(out, "taxonomy.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# provenance of the reference databases
sha256 <- function(f) strsplit(system2("sha256sum", f, stdout = TRUE), " ")[[1]][1]
write.table(data.frame(key = c("database", "species_database", "version_label", "sha256", "species_sha256", "bytes", "species_bytes", "min_boot", "classified_at"),
                       value = c(basename(opt("db")), basename(opt("species_db")), opt("version_label", "unknown"), sha256(opt("db")), sha256(opt("species_db")),
                                 file.size(opt("db")), file.size(opt("species_db")), opt("min_boot", 50), format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))),
            "taxonomy_db_info.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

plot_composition(tab, out, "unfiltered")
