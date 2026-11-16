#!/usr/bin/env Rscript
# Alpha and beta diversity against the variable of interest (continuous or categorical).
suppressPackageStartupMessages({ library(vegan); library(ggplot2) })
a <- commandArgs(trailingOnly = TRUE)
opt <- function(k, d = NULL) { i <- match(paste0("--", k), a); if (is.na(i)) d else a[i + 1] }
notes <- data.frame(case = character(), observation = character(), why = character(), lesson = character())
note <- function(...) notes[nrow(notes) + 1, ] <<- list(...)
okabe <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442")

feat <- read.delim(opt("features"), check.names = FALSE, row.names = 1)
meta <- read.delim(opt("metadata"), stringsAsFactors = FALSE)
info <- read.delim(opt("info"), stringsAsFactors = FALSE)
vname <- info$value[info$key == "variable"]; vtype <- info$value[info$key == "type"]
set.seed(as.integer(opt("seed", 42)))
meta <- meta[match(colnames(feat), meta$sample_id), ]
v <- if (vtype == "numeric") as.numeric(meta$value) else factor(meta$value)
mat <- t(as.matrix(feat))

depth <- min(rowSums(mat))
rar <- rrarefy(mat, depth)
alpha <- data.frame(sample_id = rownames(rar), value = if (vtype == "numeric") v else as.character(v),
                    observed = specnumber(rar), shannon = diversity(rar, "shannon"), simpson = diversity(rar, "simpson"))
alpha$pielou_evenness <- ifelse(alpha$observed > 1, alpha$shannon / log(alpha$observed), NA)
write.table(alpha, "alpha_diversity.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

metrics <- c("observed", "shannon", "simpson", "pielou_evenness")
stat_rows <- lapply(metrics, function(m) {
  x <- alpha[[m]]
  if (vtype == "numeric") {
    ct <- suppressWarnings(cor.test(x, v, method = "spearman"))
    data.frame(metric = m, test = "Spearman correlation", estimate = unname(ct$estimate), p_value = ct$p.value)
  } else {
    kt <- kruskal.test(x ~ v)
    data.frame(metric = m, test = "Kruskal-Wallis", estimate = unname(kt$statistic), p_value = kt$p.value)
  }
})
alpha_stats <- do.call(rbind, stat_rows)
alpha_stats$q_value_BH <- p.adjust(alpha_stats$p_value, "BH")
write.table(alpha_stats, "alpha_stats.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

long <- do.call(rbind, lapply(metrics, function(m) data.frame(sample_id = alpha$sample_id, value = v, metric = m, y = alpha[[m]])))
long$metric <- factor(long$metric, levels = metrics)
p <- ggplot(long, aes(value, y))
p <- if (vtype == "numeric") {
  p + geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = okabe[2], fill = okabe[2], alpha = 0.2) +
    geom_point(size = 3, colour = okabe[1])
} else {
  p + geom_boxplot(outlier.shape = NA, fill = NA, colour = "grey40") +
    geom_jitter(aes(colour = value), width = 0.1, size = 3, show.legend = FALSE) + scale_colour_manual(values = okabe)
}
p <- p + facet_wrap(~metric, scales = "free_y", nrow = 1) + labs(x = vname, y = "alpha-diversity value") +
  theme_bw(base_size = 12) + theme(strip.background = element_rect(fill = "grey92"))
ggsave("alpha_diversity.png", p, width = 12, height = 3.8, dpi = 150)

# beta diversity: Bray-Curtis (rarefied) and Aitchison (CLR of raw counts)
clr <- t(apply(mat + 0.5, 1, function(x) log(x) - mean(log(x))))
dists <- list("Bray-Curtis" = vegdist(rar, "bray"), "Aitchison" = dist(clr))
coords <- list(); bstats <- list()
for (nm in names(dists)) {
  d <- dists[[nm]]
  pc <- cmdscale(d, k = 2, eig = TRUE)
  ve <- 100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0])
  coords[[nm]] <- data.frame(metric = sprintf("%s (PC1 %.1f%%, PC2 %.1f%%)", nm, ve[1], ve[2]),
                             distance = nm, sample_id = rownames(pc$points), PC1 = pc$points[, 1], PC2 = pc$points[, 2],
                             value = v[match(rownames(pc$points), meta$sample_id)])
  ad <- adonis2(d ~ value, data = data.frame(value = v), permutations = 999)
  row <- data.frame(distance = nm, test = "PERMANOVA (adonis2)", R2 = ad$R2[1], statistic = ad$F[1], p_value = ad$`Pr(>F)`[1])
  if (vtype == "categorical") {
    pd <- tryCatch(permutest(betadisper(d, v), permutations = 999)$tab[1, "Pr(>F)"], error = function(e) NA)
    row$dispersion_p <- pd
  } else row$dispersion_p <- NA
  bstats[[nm]] <- row
}
coords <- do.call(rbind, coords); bstats <- do.call(rbind, bstats)
write.table(coords, "pcoa_coords.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bstats, "beta_stats.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(as.matrix(dists[["Bray-Curtis"]]), "bray_curtis_matrix.tsv", sep = "\t", quote = FALSE, col.names = NA)

b <- ggplot(coords, aes(PC1, PC2, colour = value)) + geom_point(size = 4) +
  geom_text(aes(label = sample_id), size = 2.6, vjust = -1, colour = "grey30", check_overlap = TRUE) +
  facet_wrap(~metric, scales = "free", nrow = 1) + labs(colour = vname) + theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey92"))
b <- b + if (vtype == "numeric") scale_colour_viridis_c() else scale_colour_manual(values = okabe)
ggsave("beta_diversity.png", b, width = 12, height = 4.6, dpi = 150)

if (min(rowSums(mat)) < 0.5 * median(rowSums(mat)))
  note("uneven sequencing depth", sprintf("shallowest sample %d reads vs median %d", min(rowSums(mat)), median(rowSums(mat))),
       "rarefying to the minimum discards most reads of deeper samples", "compare the Aitchison result, which uses all reads")
write.table(data.frame(key = c("rarefaction_depth", "n_samples", "n_asv"), value = c(depth, nrow(mat), ncol(mat))),
            "diversity_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(notes, "notable_cases_diversity.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
