# Real-data check of rsemverse on ENCODE RSEM output.
#
# Downloads six RSEM gene quantification files (about 65 MB; K562 cells,
# GRCh38, GENCODE V29) into a work directory, verifies their MD5 sums, and
# reports how rsem_read() and the functions that prepare the table handle
# them. Where edgeR, DESeq2 and limma are installed, the normalized values are
# compared with theirs.
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
samples$group <- sub("_[0-9]+$", "", samples$sample_id)
samples$replicate <- as.integer(sub("^.*_", "", samples$sample_id))
writeLines(c("sample_id,ref,group,replicate",
             paste(samples$sample_id, samples$ref, samples$group, samples$replicate, sep = ",")), metadata)
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

# ---- Functions that prepare the table (rsemverse 0.2.0) ---------------------

timed <- function(expr) {
  t <- system.time(v <- suppressMessages(force(expr)))[["elapsed"]]
  attr(v, "seconds") <- t
  v
}
secs <- function(v) round(attr(v, "seconds"), 2)

cat("
== rsem_filter
")
f_edger <- timed(rsem_filter(count, expr = "count"))
f_detect <- timed(rsem_filter(tpm, expr = "tpm"))
f_nonzero <- timed(rsem_filter(count, expr = "count", method = "nonzero"))
f_var <- timed(rsem_filter(tpm, expr = "tpm", method = "variable"))
cat("edger (count):         ", nrow(f_edger), "of", nrow(count), "genes,", secs(f_edger), "s
")
cat("detected (TPM >= 1, >= 2 of 6):", nrow(f_detect), "genes,", secs(f_detect), "s
")
cat("nonzero (count):       ", nrow(f_nonzero), "genes
")
cat("variable (top 2000):   ", nrow(f_var), "genes,", secs(f_var), "s
")
cat("share of expected count kept by edger, per sample (%):
")
print(round(100 * colSums(f_edger) / colSums(count), 2))
targets <- c("PTBP1", "SRSF1", "PTBP2")
cat("knockdown targets and PTBP2 kept by edger:", paste(targets, targets %in% rownames(f_edger), collapse = ", "), "
")

cat("
== rsem_subset (on all", nrow(count), "rows)
")
groups <- c("protein_coding", "mito", "ribo", "hb", "sex", "ercc")
member <- vapply(groups, function(g) {
  rownames(count) %in% rownames(suppressMessages(rsem_subset(count, "human", keep = g)))
}, logical(nrow(count)))
cat("genes per group:
"); print(colSums(member))
cat("share of expected count per group, per sample (%):
")
print(round(100 * t(vapply(groups, function(g) colSums(count[member[, g], , drop = FALSE]) / colSums(count),
                          numeric(ncol(count)))), 2))
cat("share of TPM per group, per sample (%):
")
print(round(100 * t(vapply(groups, function(g) colSums(tpm[member[, g], , drop = FALSE]) / colSums(tpm),
                          numeric(ncol(tpm)))), 2))
cat("chromosome Y genes with nonzero count in any sample (K562 is from a female donor):",
    sum(member[, "sex"] & rownames(count) != "XIST" & rowSums(count) > 0), "
")
cat("XIST mean TPM:", round(mean(unlist(tpm["XIST", ])), 2), "
")
s_time <- timed(rsem_subset(count, "human", keep = "protein_coding", remove = c("mito", "ribo", "hb", "sex")))
cat("protein-coding without mito/ribo/hb/sex:", nrow(s_time), "genes,", secs(s_time), "s
")

cat("
== rsem_normalize (on the edger-filtered counts)
")
methods <- c("cpm", "tmm", "upper_quartile", "deseq2", "quantile")
norm <- lapply(stats::setNames(methods, methods), function(m) timed(rsem_normalize(f_edger, "count", m)))
for (m in methods) {
  f <- attr(norm[[m]], "factors")
  cat(sprintf("%-15s %5.2f s%s
", m, secs(norm[[m]]),
              if (is.null(f)) "" else paste0("  factors: ", paste(signif(f, 4), collapse = ", "))))
}
cm <- as.matrix(f_edger)
ref <- list()
if (requireNamespace("edgeR", quietly = TRUE)) {
  ref$cpm <- edgeR::cpm(cm)
  ref$tmm <- edgeR::cpm(cm, lib.size = colSums(cm) * edgeR::calcNormFactors(cm, method = "TMM"))
  ref$upper_quartile <- edgeR::cpm(cm, lib.size = colSums(cm) * edgeR::calcNormFactors(cm, method = "upperquartile"))
  ref_keep <- suppressMessages(edgeR::filterByExpr(as.matrix(count)))
  cat("edger filter identical to edgeR::filterByExpr():", identical(rownames(f_edger), rownames(count)[ref_keep]), "
")
}
if (requireNamespace("DESeq2", quietly = TRUE)) {
  ref$deseq2 <- sweep(cm, 2L, DESeq2::estimateSizeFactorsForMatrix(cm), "/")
}
if (requireNamespace("limma", quietly = TRUE)) ref$quantile <- limma::normalizeQuantiles(cm)
for (m in names(ref)) {
  d <- max(abs(as.matrix(norm[[m]]) - ref[[m]]) / pmax(1, abs(ref[[m]])))
  cat(sprintf("largest relative difference from the reference, %-15s %.1e
", m, d))
}
cat("knockdown ratios after TMM (mean TMM-CPM, KD / control):
")
for (g in c("PTBP1", "SRSF1")) {
  v <- unlist(norm$tmm[g, ])
  cat(sprintf("  %s in %s knockdown: %.2f
", g, g, mean(v[grp == g]) / mean(v[grp == "CTRL"])))
}

cat("
== rsem_transform
")
lz <- timed(rsem_transform(norm$tmm, c("log2", "zscore")))
rk <- timed(rsem_transform(tpm, "rank"))
cat("log2 + zscore:", nrow(lz), "genes,", secs(lz), "s; row means within 1e-12 of 0:",
    all(abs(rowMeans(lz)) < 1e-12), "
")
cat("rank on all rows:", nrow(rk), "genes,", secs(rk), "s; range", paste(signif(range(as.matrix(rk)), 3), collapse = " to "), "
")

cat("
== rsem_long
")
lg <- timed(rsem_long(rsem_transform(norm$tmm, "log2"), metadata = metadata))
cat(nrow(lg), "rows,", secs(lg), "s; columns:", paste(sprintf("%s <%s>", names(lg), vapply(lg, function(v) class(v)[1], "")), collapse = ", "), "
")
