#!/usr/bin/env Rscript
# Microbial co-occurrence network: Spearman correlation on CLR-transformed taxa abundances.
# A numeric variable of interest is added as an extra node so taxa linked to it stand out.
suppressPackageStartupMessages(library(igraph))
a <- commandArgs(trailingOnly = TRUE)
opt <- function(k, d = NULL) { i <- match(paste0("--", k), a); if (is.na(i)) d else a[i + 1] }
notes <- data.frame(case = character(), observation = character(), why = character(), lesson = character())
note <- function(...) notes[nrow(notes) + 1, ] <<- list(...)

counts <- read.delim(opt("counts"), check.names = FALSE, row.names = 1)
tax    <- read.delim(opt("taxonomy"), stringsAsFactors = FALSE, row.names = 1)
meta   <- read.delim(opt("metadata"), stringsAsFactors = FALSE)
info   <- read.delim(opt("info"), stringsAsFactors = FALSE)
min_prev <- as.numeric(opt("min_prev", 0.3)); cor_thr <- as.numeric(opt("cor", 0.6))
p_thr <- as.numeric(opt("p", 0.05)); padj_method <- opt("padj_method", "none"); max_nodes <- as.integer(opt("max_nodes", 60))
vname <- info$value[info$key == "variable"]; vtype <- info$value[info$key == "type"]
meta <- meta[match(colnames(counts), meta$sample_id), ]

keep <- rowSums(counts > 0) >= max(3, ceiling(min_prev * ncol(counts)))
counts <- counts[keep, , drop = FALSE]
counts <- counts[head(order(rowSums(counts), decreasing = TRUE), max_nodes), , drop = FALSE]
n_tested <- nrow(counts)
clr <- t(apply(counts + 0.5, 2, function(x) log(x) - mean(log(x))))   # samples x taxa
if (vtype == "numeric") clr <- cbind(clr, setNames(data.frame(as.numeric(meta$value)), paste0("[", vname, "]")))
clr <- as.matrix(clr)

pairs <- if (ncol(clr) >= 2) t(combn(colnames(clr), 2)) else matrix(character(), 0, 2)
edges <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
  ct <- suppressWarnings(cor.test(clr[, pairs[i, 1]], clr[, pairs[i, 2]], method = "spearman", exact = FALSE))
  data.frame(from = pairs[i, 1], to = pairs[i, 2], rho = unname(ct$estimate), p = ct$p.value)
}))
if (is.null(edges)) edges <- data.frame(from = character(), to = character(), rho = numeric(), p = numeric())
edges$p_adj <- p.adjust(edges$p, padj_method)
edges <- edges[!is.na(edges$rho) & abs(edges$rho) >= cor_thr & edges$p_adj <= p_thr, , drop = FALSE]
edges$sign <- ifelse(edges$rho > 0, "positive", "negative")
write.table(edges, "network_edges.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

nodes <- data.frame(name = colnames(clr), stringsAsFactors = FALSE)
nodes$phylum <- ifelse(nodes$name %in% rownames(tax), tax[nodes$name, "Phylum"], "variable")
nodes$phylum[is.na(nodes$phylum)] <- "Unassigned"
g <- graph_from_data_frame(edges[, c("from", "to", "rho", "sign")], directed = FALSE, vertices = nodes)
nodes$degree <- degree(g)
nodes$betweenness <- if (ecount(g) > 0) betweenness(g) else 0
nodes$module <- NA_integer_
if (ecount(g) > 0) nodes$module <- as.integer(membership(cluster_louvain(g, weights = abs(E(g)$rho))))
write.table(nodes, "network_nodes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

png("network.png", width = 1500, height = 1200, res = 150)
if (ecount(g) == 0) {
  plot.new(); text(0.5, 0.5, sprintf("No edges with |rho| >= %g and p <= %g\n(%d taxa tested)", cor_thr, p_thr, n_tested), cex = 1.2)
} else {
  g2 <- delete_vertices(g, which(degree(g) == 0))
  phyla <- sort(unique(V(g2)$phylum)); pal <- setNames(hcl.colors(max(3, length(phyla)), "Dark 3")[seq_along(phyla)], phyla)
  is_var <- V(g2)$name %in% paste0("[", vname, "]")
  pal["variable"] <- "#000000"
  set.seed(1)
  par(mar = c(1, 1, 1, 1))
  plot(g2, layout = layout_with_fr(g2, weights = abs(E(g2)$rho)),
       vertex.color = pal[V(g2)$phylum], vertex.shape = ifelse(is_var, "square", "circle"),
       vertex.size = 6 + 2.5 * degree(g2), vertex.label = V(g2)$name, vertex.label.cex = 0.65,
       vertex.label.color = "black", vertex.label.dist = 0, vertex.frame.color = "white",
       edge.color = ifelse(E(g2)$sign == "positive", "#0072B2", "#D55E00"), edge.width = 1 + 4 * abs(E(g2)$rho))
  legend("topleft", legend = names(pal)[names(pal) %in% V(g2)$phylum], fill = pal[names(pal) %in% V(g2)$phylum],
         title = "Phylum (node colour)", cex = 0.7, bty = "n")
  legend("bottomleft", legend = c("positive correlation", "negative correlation"), col = c("#0072B2", "#D55E00"),
         lwd = 3, cex = 0.7, bty = "n")
}
invisible(dev.off())

write.table(data.frame(key = c("n_taxa_tested", "n_nodes_with_edges", "n_edges", "n_positive", "n_negative", "n_modules", "min_abs_rho", "p_threshold", "p_adjust"),
                       value = c(n_tested, sum(nodes$degree > 0), nrow(edges), sum(edges$rho > 0), sum(edges$rho < 0),
                                 length(na.omit(unique(nodes$module))), cor_thr, p_thr, padj_method)),
            "network_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
if (nrow(edges) == 0) note("empty network", sprintf("0 of %d candidate pairs pass |rho| >= %g, p <= %g", nrow(pairs), cor_thr, p_thr),
                           "too few samples/taxa or thresholds too strict", "lower --net_cor / --net_min_prev or add samples")
if (padj_method == "none" && nrow(edges) > 0) note("unadjusted edge p-values", sprintf("%d edges kept at raw p <= %g", nrow(edges), p_thr),
                           "with few samples, FDR correction would remove all edges", "treat edges as hypotheses, not findings")
write.table(notes, "notable_cases_network.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
