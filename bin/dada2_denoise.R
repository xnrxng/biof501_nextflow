#!/usr/bin/env Rscript
# DADA2: filter/trim -> learn errors -> denoise -> merge pairs -> remove chimeras -> ASV table.
# Reads *.fastq.gz in the working directory: <id>_1/_2.fastq.gz (paired) or <id>.fastq.gz (single).
suppressPackageStartupMessages({ library(dada2); library(ShortRead) })
a <- commandArgs(trailingOnly = TRUE)
opt <- function(k, d = NULL) { i <- match(paste0("--", k), a); if (is.na(i)) d else a[i + 1] }
notes <- data.frame(case = character(), observation = character(), why = character(), lesson = character())
note <- function(...) notes[nrow(notes) + 1, ] <<- list(...)

threads <- as.integer(opt("threads", 1))
max_ee  <- as.numeric(opt("max_ee", 2))
trim_f  <- as.integer(opt("trim_left_f", 0)); trim_r <- as.integer(opt("trim_left_r", 0))

fq      <- sort(list.files(".", pattern = "\\.fastq\\.gz$"))
samples <- unique(sub("(_1|_2)?\\.fastq\\.gz$", "", fq))
f1 <- ifelse(file.exists(paste0(samples, "_1.fastq.gz")), paste0(samples, "_1.fastq.gz"), paste0(samples, ".fastq.gz"))
f2 <- ifelse(file.exists(paste0(samples, "_2.fastq.gz")), paste0(samples, "_2.fastq.gz"), NA)
paired <- !anyNA(f2)
if (!paired && any(!is.na(f2))) {
  note("mixed layout", "some samples are paired, some single-end", "pairs cannot be merged for all samples",
       "using forward reads only for every sample")
}
names(f1) <- samples; if (paired) names(f2) <- samples

# truncation: "auto" = first cycle whose median quality drops below Q30 (0 = no truncation)
auto_trunc <- function(files, q = 30, floor_len = 50) {
  meds <- lapply(head(files, 3), function(f) {
    fs <- FastqSampler(f, 5000); on.exit(close(fs)); set.seed(1)
    apply(as(quality(yield(fs)), "matrix"), 2, median, na.rm = TRUE)
  })
  L <- max(lengths(meds))
  med <- vapply(seq_len(L), function(i) median(vapply(meds, function(x) if (i <= length(x)) x[i] else NA_real_, 0), na.rm = TRUE), 0)
  below <- which(med < q)
  len <- if (length(below)) below[1] - 1L else 0L
  if (len > 0 && len < floor_len) len <- floor_len
  len
}
tl <- function(v, files) if (identical(v, "auto")) auto_trunc(files) else as.integer(v)
trunc_f <- tl(opt("trunc_f", "auto"), f1)
trunc_r <- if (paired) tl(opt("trunc_r", "auto"), f2) else 0L
cat(sprintf("paired=%s trunc_f=%d trunc_r=%d\n", paired, trunc_f, trunc_r))

png("quality_profiles.png", width = 1400, height = 900, res = 130)
print(plotQualityProfile(f1, aggregate = TRUE) +
        ggplot2::ggtitle(if (paired) "Forward reads (all samples)" else "Reads (all samples)"))
invisible(dev.off())
if (paired) {
  png("quality_profiles_reverse.png", width = 1400, height = 900, res = 130)
  print(plotQualityProfile(f2, aggregate = TRUE) + ggplot2::ggtitle("Reverse reads (all samples)"))
  invisible(dev.off())
}

dir.create("filtered", showWarnings = FALSE)
ff <- file.path("filtered", basename(f1)); names(ff) <- samples
fr <- if (paired) setNames(file.path("filtered", basename(f2)), samples) else NULL
out <- if (paired) {
  filterAndTrim(f1, ff, f2, fr, truncLen = c(trunc_f, trunc_r), trimLeft = c(trim_f, trim_r), maxN = 0,
                maxEE = c(max_ee, max_ee), truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = threads)
} else {
  filterAndTrim(f1, ff, truncLen = trunc_f, trimLeft = trim_f, maxN = 0, maxEE = max_ee, truncQ = 2,
                rm.phix = TRUE, compress = TRUE, multithread = threads)
}
rownames(out) <- samples
keep <- file.exists(ff) & out[, "reads.out"] > 0
if (paired) keep <- keep & file.exists(fr)
for (s in samples[!keep]) note(paste("sample", s, "lost all reads"), "0 reads after filtering",
                               "sample cannot be analysed", "loosen --max_ee or the truncation lengths")
if (sum(keep) < 3) stop("Fewer than 3 samples survived filtering; adjust trunc_len_f/r and max_ee")
ff <- ff[keep]; if (paired) fr <- fr[keep]

errF <- learnErrors(ff, multithread = threads, verbose = FALSE)
ddF  <- dada(ff, err = errF, multithread = threads, verbose = FALSE)
if (paired) {
  errR <- learnErrors(fr, multithread = threads, verbose = FALSE)
  ddR  <- dada(fr, err = errR, multithread = threads, verbose = FALSE)
}
png("error_rates.png", width = 1400, height = 1000, res = 130)
print(plotErrors(errF, nominalQ = TRUE))
invisible(dev.off())

getN <- function(x) sum(getUniques(x))
if (paired) {
  mergers <- mergePairs(ddF, ff, ddR, fr, verbose = FALSE)
  merged  <- vapply(mergers, function(x) sum(x$abundance[x$accept]), 0)
  frac    <- sum(merged) / sum(vapply(ddF, getN, 0))
  if (frac < 0.2) {
    note("low merge rate", sprintf("only %.1f%% of denoised forward reads merged", 100 * frac),
         "truncation left too little forward/reverse overlap",
         "fell back to forward reads only; raise trunc_len_f/r or check the amplicon length")
    paired <- FALSE
    seqtab <- makeSequenceTable(ddF)
  } else seqtab <- makeSequenceTable(mergers)
} else seqtab <- makeSequenceTable(ddF)

seqtab_nc <- removeBimeraDenovo(seqtab, method = "consensus", multithread = threads, verbose = FALSE)
if (ncol(seqtab_nc) == 0) stop("No ASVs left after chimera removal")

ord  <- order(colSums(seqtab_nc), decreasing = TRUE)
seqs <- colnames(seqtab_nc)[ord]
ids  <- sprintf("ASV_%04d", seq_along(seqs))
tab  <- t(seqtab_nc[, seqs, drop = FALSE]); rownames(tab) <- ids
write.table(data.frame(asv_id = ids, tab, check.names = FALSE), "asv_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(rbind(paste0(">", ids), seqs)), "asv_seqs.fasta")

track <- data.frame(sample_id = samples, input = out[, "reads.in"], filtered = out[, "reads.out"])
track$denoised <- NA_integer_; track$merged <- NA_integer_
track$denoised[match(names(ddF), samples)] <- vapply(ddF, getN, 0)
if (exists("merged")) track$merged[match(names(merged), samples)] <- merged
track$nonchim <- NA_integer_; track$nonchim[match(rownames(seqtab_nc), samples)] <- rowSums(seqtab_nc)
write.table(track, "read_tracking.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(notes, "notable_cases_dada2.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
