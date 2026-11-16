#!/usr/bin/env Rscript
# classification of ASVs against a SILVA-style training set
suppressPackageStartupMessages({ library(dada2); library(Biostrings) })
a <- commandArgs(trailingOnly = TRUE)
opt <- function(k, d = NULL) { i <- match(paste0("--", k), a); if (is.na(i)) d else a[i + 1] }

for (k in c("seqs", "db", "species_db"))
  if (is.null(opt(k))) stop("Missing required argument --", k,
    "\nUsage: 04-assign_taxonomy.R --seqs <asv.fasta> --db <train_set.fa.gz> --species_db <species_assignment.fa.gz> [--threads N] [--min_boot N]", call. = FALSE)
for (k in c("seqs", "db", "species_db"))
  if (!file.exists(opt(k))) stop("File for --", k, " not found: ", opt(k), call. = FALSE)

seqs <- readDNAStringSet(opt("seqs"))
if (!length(seqs)) stop("No sequences in ", opt("seqs"), call. = FALSE)
db   <- opt("db")
tax  <- assignTaxonomy(as.character(seqs), db, multithread = as.integer(opt("threads", 1)),
                       minBoot = as.integer(opt("min_boot", 50)), tryRC = TRUE, verbose = FALSE)


# species classification
tax  <- addSpecies(tax, opt("species_db"), tryRC = TRUE, verbose = FALSE)
out <-data.frame(asv_id = names(seqs), tax[as.character(seqs), , drop = FALSE], row.names = NULL)
write.table(out, "taxonomy.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
