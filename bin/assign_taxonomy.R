#!/usr/bin/env Rscript
# Naive-Bayes (RDP) classification of ASVs against a SILVA-style training set.
suppressPackageStartupMessages({ library(dada2); library(Biostrings) })
a <- commandArgs(trailingOnly = TRUE)
opt <- function(k, d = NULL) { i <- match(paste0("--", k), a); if (is.na(i)) d else a[i + 1] }

seqs <- readDNAStringSet(opt("seqs"))
db   <- opt("db")
tax  <- assignTaxonomy(as.character(seqs), db, multithread = as.integer(opt("threads", 1)),
                       minBoot = as.integer(opt("min_boot", 50)), tryRC = TRUE, verbose = FALSE)
out <- data.frame(asv_id = names(seqs), tax[as.character(seqs), , drop = FALSE], row.names = NULL)
write.table(out, "taxonomy.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

sha <- strsplit(system2("sha256sum", db, stdout = TRUE), " ")[[1]][1]
write.table(data.frame(key = c("database", "version_label", "sha256", "bytes", "classified_at"),
                       value = c(basename(db), opt("version_label", "unknown"), sha, file.size(db),
                                 format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))),
            "taxonomy_db_info.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
