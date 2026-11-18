#!/usr/bin/env Rscript
# rarefy, Shannon alpha diversity and Bray-Curtis beta diversity against the variable of interest
suppressPackageStartupMessages({ library(vegan); library(ggplot2) })
a <- commandArgs(trailingOnly = TRUE)
opt <- function(k, d = NULL) { i <- match(paste0("--", k), a); if (is.na(i)) d else a[i + 1] }

feat <- read.delim(opt("features"), check.names = FALSE, row.names = 1)
meta <- read.delim(opt("metadata"), stringsAsFactors = FALSE)
info <- read.delim(opt("info"), stringsAsFactors = FALSE)
feat <- feat[, colnames(feat) %in% meta$sample_id, drop = FALSE]  
vname <- info$value[info$key == "variable"]; vtype <- info$value[info$key == "type"]
set.seed(as.integer(opt("seed", 42)))
meta <- meta[match(colnames(feat), meta$sample_id), ]
mat <- t(as.matrix(feat))                                       # samples x ASVs

# rarefaction depth: user-chosen (--depth > 0) or, by default, the smallest sample; samples below a chosen depth are dropped
req   <- as.integer(opt("depth", 0))
depth <- if (req > 0) req else min(rowSums(mat))
short <- rowSums(mat) < depth
if (any(short)) {
  message("Dropping ", sum(short), " sample(s) with fewer than ", depth, " reads: ",
          paste(sprintf("%s (%d)", rownames(mat)[short], rowSums(mat)[short]), collapse = ", "))
  mat <- mat[!short, , drop = FALSE]; meta <- meta[!short, , drop = FALSE]
}
if (nrow(mat) < 4) stop("Fewer than 4 samples have at least ", depth, " reads; lower --rarefaction_depth", call. = FALSE)
v <- if (vtype == "numeric") as.numeric(meta$value) else droplevels(factor(meta$value))
if (vtype == "categorical" && nlevels(v) < 2) stop("Only one level of '", vname, "' remains after dropping shallow samples; lower --rarefaction_depth", call. = FALSE)

stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

# rarefaction
rar   <- suppressWarnings(rrarefy(mat, depth))     
write.table(data.frame(asv_id = colnames(rar), t(rar), check.names = FALSE), "rarefied_filtered_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

curve <- suppressWarnings(rarecurve(mat, step = max(1, ceiling(max(rowSums(mat)) / 60)), tidy = TRUE))  
rc <- ggplot(curve, aes(Sample, Species, colour = Site)) + geom_line(linewidth = 0.8) +
  geom_vline(xintercept = depth, linetype = "dashed") +
  annotate("text", x = depth, y = Inf, label = sprintf(" rarefaction depth = %s reads", format(depth, big.mark = ",")), hjust = 0, vjust = 1.5, size = 3.5) +
  scale_x_continuous(labels = scales::label_comma()) +
  labs(x = "Reads sampled", y = "Observed ASVs", colour = "Sample", title = "Rarefaction curves") + theme_classic(base_size = 12)
ggsave("rarefaction_curve.png", rc, width = 8, height = 5, dpi = 150)

# alpha diversity: Shannon index
alpha <- data.frame(sample_id = rownames(rar), value = if (vtype == "numeric") v else as.character(v), shannon = diversity(rar, "shannon"))
write.table(alpha, "alpha_diversity.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

x <- alpha$shannon
if (vtype == "numeric") {
  ct <- suppressWarnings(cor.test(x, v, method = "spearman", exact = FALSE))
  alpha_stats <- data.frame(metric = "shannon", test = "Spearman correlation", comparison = vname, statistic = unname(ct$statistic),
                            estimate = unname(ct$estimate), p_value = ct$p.value, p_adjusted = NA_real_)
  a_lab <- sprintf("Spearman rho = %.2f, p = %s %s", ct$estimate, fmt_p(ct$p.value), stars(ct$p.value))
} else if (nlevels(v) == 2) {
  wt <- suppressWarnings(wilcox.test(x ~ v, exact = FALSE))
  alpha_stats <- data.frame(metric = "shannon", test = "Wilcoxon rank-sum", comparison = paste(levels(v), collapse = " vs "), statistic = unname(wt$statistic),
                            estimate = NA_real_, p_value = wt$p.value, p_adjusted = NA_real_)
  a_lab <- sprintf("Wilcoxon rank-sum, W = %.1f, p = %s %s", wt$statistic, fmt_p(wt$p.value), stars(wt$p.value))
} else {
  kt <- kruskal.test(x ~ v)
  pw <- suppressWarnings(pairwise.wilcox.test(x, v, p.adjust.method = "BH", exact = FALSE))$p.value
  pr <- which(!is.na(pw), arr.ind = TRUE)
  alpha_stats <- rbind(
    data.frame(metric = "shannon", test = "Kruskal-Wallis", comparison = "all groups", statistic = unname(kt$statistic), estimate = NA_real_,
               p_value = kt$p.value, p_adjusted = NA_real_),
    data.frame(metric = "shannon", test = "Pairwise Wilcoxon (BH-adjusted)", comparison = paste(rownames(pw)[pr[, 1]], "vs", colnames(pw)[pr[, 2]]),
               statistic = NA_real_, estimate = NA_real_, p_value = NA_real_, p_adjusted = pw[pr]))
  a_lab <- sprintf("Kruskal-Wallis, chi-sq = %.2f, p = %s %s", kt$statistic, fmt_p(kt$p.value), stars(kt$p.value))
}
alpha_stats$significance <- stars(ifelse(is.na(alpha_stats$p_adjusted), alpha_stats$p_value, alpha_stats$p_adjusted))
write.table(alpha_stats, "alpha_stats.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

adf <- data.frame(value = v, shannon = x)
p <- ggplot(adf, aes(value, shannon))
p <- if (vtype == "numeric") {
  p + geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha = 0.2) + geom_point(size = 3)
} else {
  p + geom_boxplot(outlier.shape = NA, fill = NA) + geom_jitter(aes(colour = value), width = 0.1, size = 3, show.legend = FALSE)
}
p <- p + annotate("text", x = -Inf, y = Inf, label = a_lab, hjust = -0.05, vjust = 1.6, size = 3.6) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = vname, y = "Shannon diversity", title = sprintf("Shannon diversity vs %s (rarefied to %s reads)", vname, format(depth, big.mark = ","))) +
  theme_classic(base_size = 12)
ggsave("alpha_diversity.png", p, width = 6.5, height = 5, dpi = 150)

# beta diversity: Bray-Curtis dissimilarity
d  <- vegdist(rar, "bray")
write.table(as.matrix(d), "bray_curtis_matrix.tsv", sep = "\t", quote = FALSE, col.names = NA)
pc <- cmdscale(d, k = 2, eig = TRUE)
ve <- 100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0])
coords <- data.frame(sample_id = rownames(pc$points), PC1 = pc$points[, 1], PC2 = pc$points[, 2],
                     value = v[match(rownames(pc$points), meta$sample_id)])
write.table(coords, "pcoa_coords.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# PERMANOVA: numeric variable = distance-based regression on the continuous predictor; categorical = group differences in centroid
ad <- adonis2(d ~ value, data = data.frame(value = v), permutations = 999)
beta_stats <- data.frame(distance = "Bray-Curtis", test = "PERMANOVA (adonis2)", R2 = ad$R2[1], statistic = ad$F[1], p_value = ad$`Pr(>F)`[1],
                         dispersion_test = NA_character_, dispersion_p_value = NA_real_)
b_lab <- sprintf("PERMANOVA, R2 = %.3f, F = %.2f, p = %s %s", ad$R2[1], ad$F[1], fmt_p(ad$`Pr(>F)`[1]), stars(ad$`Pr(>F)`[1]))
if (vtype == "categorical") {

  pd <- tryCatch(permutest(betadisper(d, v), permutations = 999)$tab[1, "Pr(>F)"], error = function(e) NA_real_)
  beta_stats$dispersion_test <- "PERMDISP (betadisper)"; beta_stats$dispersion_p_value <- pd
  b_lab <- paste0(b_lab, sprintf("\nPERMDISP p = %s %s", fmt_p(pd), stars(pd)))
}
beta_stats$significance <- stars(beta_stats$p_value)
write.table(beta_stats, "beta_stats.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

b <- ggplot(coords, aes(PC1, PC2, colour = value)) + geom_point(size = 4) +
  geom_text(aes(label = sample_id), size = 2.6, vjust = -1, colour = "grey30", check_overlap = TRUE) +
  annotate("text", x = -Inf, y = Inf, label = b_lab, hjust = -0.05, vjust = 1.5, size = 3.6) +
  scale_x_continuous(expand = expansion(mult = 0.1)) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]), colour = vname,
       title = "PCoA of Bray-Curtis dissimilarity") + theme_classic(base_size = 12)
ggsave("beta_diversity.png", b, width = 7, height = 5.5, dpi = 150)
cat(sprintf("depth=%d samples=%d variable=%s type=%s\n", depth, nrow(mat), vname, vtype))
