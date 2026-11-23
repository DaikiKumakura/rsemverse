# Validation of 0.7.0 on six ENCODE K562 RSEM runs (gene and transcript files; GENCODE V29).
# Run from the package root:
#   Rscript validation/encode-k562-0.7.0.R <folder with *.genes.results and *.isoforms.results>
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, dir.exists(args[[1L]]))
suppressMessages(devtools::load_all(quiet = TRUE))
dir <- args[[1L]]

# FPKM: read as a column, and consistent with the RSEM definition FPKM = 1e3 / l_bar * TPM.
fp <- suppressMessages(rsem_read(dir, metadata = NULL, expr = "fpkm", id = "original"))
tp <- suppressMessages(rsem_read(dir, metadata = NULL, expr = "tpm", id = "original"))
ln <- suppressMessages(rsem_read(dir, metadata = NULL, expr = "length", id = "original"))
lbar <- colSums(tp * ln) / 1e6
ratio <- vapply(seq_len(ncol(fp)), function(j) {
  k <- tp[[j]] > 10
  stats::median(fp[[j]][k] / (1e3 / lbar[[j]] * tp[[j]][k]))
}, numeric(1))
cat("median FPKM / (1e3 / l_bar * TPM) over genes with TPM > 10, per sample:\n"); print(round(ratio, 4))
stopifnot(all(abs(ratio - 1) < 0.01))

# Transcript lengths are the same in every file; gene lengths are not.
chk <- suppressMessages(rsem_check(dir, level = "isoform"))
cat("\nrsem_check(level = 'isoform') findings:\n"); print(chk[, c("severity", "check", "sample_id")])
stopifnot(!"transcript_length" %in% chk$check)

# Alignment histogram is not available: the downloads have no .cnt files.

# Gene-level values for DESeq2.
a <- suppressMessages(rsem_import(dir, metadata = NULL))
zl <- attr(a, "zero_length")
cat("\ngenes with effective length 0 per sample:\n"); print(zl)
txi <- suppressMessages(rsem_as_tximport(a, zero_length = "drop"))
cat("genes kept for DESeq2:", nrow(txi$counts), "of", nrow(a$count), "\n")
stopifnot(identical(dim(txi$counts), dim(txi$length)), all(txi$length > 0))
if (requireNamespace("tximport", quietly = TRUE)) {
  files <- list.files(dir, pattern = "[.]genes[.]results$", full.names = TRUE)
  names(files) <- sub("[.]genes[.]results$", "", basename(files))
  ref <- suppressMessages(tximport::tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE))
  common <- rownames(txi$counts)
  cat("\nmax abs difference to tximport(type = 'rsem'): counts",
      max(abs(txi$counts - ref$counts[common, colnames(txi$counts)])),
      " length", max(abs(txi$length - ref$length[common, colnames(txi$length)])), "\n")
}

# RLE of the samples.
x <- suppressMessages(rsem_read(dir, metadata = NULL, expr = "tpm", id = "original"))
qc <- suppressMessages(rsem_sample_qc(x, n_top = 2000))
cat("\nrelative log expression per sample:\n"); print(qc$sample_summary[, c("sample_id", "rle_median", "rle_iqr")])
