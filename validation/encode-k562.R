# Real-data check of rsem_read() on ENCODE RSEM output.
#
# Downloads six RSEM gene quantification files (about 65 MB; K562 cells,
# GRCh38, GENCODE V29) into a work directory, verifies their MD5 sums, and
# reports how rsem_read() handles them.
#
#   Rscript validation/encode-k562.R [work-directory]
#
# Needs an installed rsemverse and org.Hs.eg.db, and internet access.

args <- commandArgs(trailingOnly = TRUE)
work <- if (length(args)) args[[1]] else file.path(tempdir(), "rsemverse-encode")
dir.create(file.path(work, "data"), recursive = TRUE, showWarnings = FALSE)

samples <- data.frame(
  sample_id = c("CTRL_1", "CTRL_2", "PTBP1_1", "PTBP1_2", "SRSF1_1", "SRSF1_2"),
  accession = c("ENCFF587DSG", "ENCFF774RTC", "ENCFF075OME", "ENCFF896CIS", "ENCFF711WSO", "ENCFF322CJL"),
  md5 = c(
    "9d64a74c8a09c01b08421a27c82664a2", "df949910f3dccaec1f6f08c021d17cd7",
    "6d6f95c21f95947fc2f2c761c06aed6d", "be3d578f82a501318ae21dc00699fcad",
    "0c0530bf1d8e1e52e97d8d0a8da3a418", "4b7891dd7da2b09a6d49bd44f64f211a"
  ),
  stringsAsFactors = FALSE
)
samples$ref <- paste0(samples$sample_id, ".genes.results")

for (i in seq_len(nrow(samples))) {
  dest <- file.path(work, "data", samples$ref[[i]])
  if (!file.exists(dest)) {
    url <- sprintf(
      "https://www.encodeproject.org/files/%s/@@download/%s.tsv",
      samples$accession[[i]], samples$accession[[i]]
    )
    utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  }
  got <- unname(tools::md5sum(dest))
  stopifnot("MD5 mismatch" = identical(got, samples$md5[[i]]))
}
metadata <- file.path(work, "metadata.csv")
writeLines(c("sample_id,ref", paste0(samples$sample_id, ",", samples$ref)), metadata)
data_dir <- file.path(work, "data")

library(rsemverse)
raw <- utils::read.delim(file.path(data_dir, samples$ref[[1]]), check.names = FALSE,
                         colClasses = c(gene_id = "character"))
n_genes <- nrow(raw)
cat("\n== Input\n")
cat("genes per file:", n_genes, "\n")
kinds <- c(
  Ensembl = sum(grepl("^ENSG[0-9]+\\.[0-9]+$", raw$gene_id)),
  PAR_Y = sum(grepl("_PAR_Y$", raw$gene_id)),
  ERCC = sum(grepl("^(gSpikein|ERCC)", raw$gene_id)),
  numeric = sum(grepl("^[0-9]+$", raw$gene_id))
)
print(kinds)

cat("\n== rsem_read(count)\n")
t_count <- system.time(count <- rsem_read(data_dir, metadata, species = "human", expr = "count"))
cat("elapsed (s):", round(t_count[["elapsed"]], 1), "\n")
cat("\n== rsem_read(tpm)\n")
tpm <- rsem_read(data_dir, metadata, species = "human", expr = "tpm")

cat("\n== Checks\n")
inputs <- lapply(samples$ref, function(f) utils::read.delim(file.path(data_dir, f), check.names = FALSE,
                                                           colClasses = c(gene_id = "character")))
count_in <- vapply(inputs, function(x) sum(x$expected_count), numeric(1))
tpm_in <- vapply(inputs, function(x) sum(x$TPM), numeric(1))
cat("column sums equal the input (count):", isTRUE(all.equal(unname(colSums(count)), count_in)), "\n")
cat("column sums equal the input (TPM):  ", isTRUE(all.equal(unname(colSums(tpm)), tpm_in)), "\n")
cat("TPM column sums (min, max):", format(range(colSums(tpm)), big.mark = ","), "\n")

# A row kept its original ID exactly when its name is a gene_id of the input
# (a symbol can never equal one: rsem_read() stops in that case).
kept <- rownames(count) %in% raw$gene_id
share <- colSums(count[kept, , drop = FALSE]) / colSums(count)
cat("rows kept under their original ID:", sum(kept), "\n")
cat("share of expected count in those rows, per sample (%):\n")
print(round(100 * share, 2))
share_tpm <- colSums(tpm[kept, , drop = FALSE]) / colSums(tpm)
cat("share of TPM in those rows, per sample (%):\n")
print(round(100 * share_tpm, 2))

cat("\nknockdown targets (mean expected count, KD / control):\n")
grp <- sub("_[0-9]+$", "", colnames(count))
for (g in c("PTBP1", "SRSF1")) {
  ratio <- mean(unlist(count[g, grp == g])) / mean(unlist(count[g, grp == "CTRL"]))
  cat(sprintf("  %s in %s knockdown: %.2f\n", g, g, ratio))
}
