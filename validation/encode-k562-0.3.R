# Real-data check of the functions added in rsemverse 0.3.0, on the same six
# ENCODE K562 experiments as encode-k562.R.
#
# Downloads the six gene quantification files of encode-k562.R, the six
# matching transcript quantification files (same RSEM runs; about 190 MB) and
# the GENCODE V29 annotation the ENCODE reference was built from (39 MB),
# verifies their MD5 sums, and prints the numbers reported in README.md.
#
#   Rscript validation/encode-k562-0.3.R [work-directory]
#
# Needs an installed rsemverse and org.Hs.eg.db, and internet access; edgeR
# is used for the TMMwsp comparison when installed.

options(timeout = max(1200, getOption("timeout")))
args <- commandArgs(trailingOnly = TRUE)
work <- if (length(args)) args[[1]] else file.path(tempdir(), "rsemverse-encode")
dir.create(file.path(work, "data"), recursive = TRUE, showWarnings = FALSE)

files <- data.frame(
  sample_id = rep(c("CTRL_1", "CTRL_2", "PTBP1_1", "PTBP1_2", "SRSF1_1", "SRSF1_2"), 2),
  level = rep(c("genes", "isoforms"), each = 6),
  accession = c(
    "ENCFF587DSG", "ENCFF774RTC", "ENCFF075OME", "ENCFF896CIS", "ENCFF711WSO", "ENCFF322CJL",
    "ENCFF875PEA", "ENCFF180IOW", "ENCFF848LJZ", "ENCFF371ODR", "ENCFF662WOK", "ENCFF180KQR"
  ),
  md5 = c(
    "9d64a74c8a09c01b08421a27c82664a2", "df949910f3dccaec1f6f08c021d17cd7",
    "6d6f95c21f95947fc2f2c761c06aed6d", "be3d578f82a501318ae21dc00699fcad",
    "0c0530bf1d8e1e52e97d8d0a8da3a418", "4b7891dd7da2b09a6d49bd44f64f211a",
    "7f834170b04972e7e4c374bd968c0aea", "23c040be43bc28b5a9ce39db0a546f4c",
    "a0b1b43a039bdf798c7ee0b596286290", "78d1727c394a582c45f64e87d9345af0",
    "7c8ee9c79de94dd253d362b3d96e8702", "010f9f8dafd9eaa3cea99864f110e61a"
  ),
  stringsAsFactors = FALSE
)
fetch <- function(url, dest, md5) {
  if (!file.exists(dest)) utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  stopifnot("MD5 mismatch" = identical(unname(tools::md5sum(dest)), md5))
}
for (i in seq_len(nrow(files))) {
  fetch(sprintf("https://www.encodeproject.org/files/%1$s/@@download/%1$s.tsv", files$accession[[i]]),
        file.path(work, "data", sprintf("%s.%s.results", files$sample_id[[i]], files$level[[i]])), files$md5[[i]])
}
gtf <- file.path(work, "gencode.v29.annotation.gtf.gz")
fetch("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz",
      gtf, "949dca487934459e6338090040e77628")

samples <- unique(files$sample_id)
meta <- data.frame(sample_id = samples, ref = paste0(samples, ".genes.results"),
                   group = sub("_[0-9]+$", "", samples), stringsAsFactors = FALSE)
metadata <- file.path(work, "metadata.csv")
utils::write.csv(meta, metadata, row.names = FALSE, quote = FALSE)
data_dir <- file.path(work, "data")

library(rsemverse)
tm <- function(expr) { t <- system.time(v <- force(expr))[["elapsed"]]; cat(sprintf("  (%.2f s)\n", t)); v }

cat("\n== rsem_read(): files found without a sample table\n")
tpm <- rsem_read(data_dir, metadata, species = "human", expr = "tpm")
found <- tm(rsem_read(data_dir, NULL, species = "human", expr = "tpm"))
cat("same table as with the sample table (columns in sample ID order):",
    identical(as.matrix(found[, colnames(tpm)]), as.matrix(tpm)), "\n")
cat("gene_ids attribute rows:", nrow(attr(tpm, "gene_ids")), "\n")

cat("\n== Effective length\n")
len <- tm(rsem_read(data_dir, metadata, species = "human", expr = "length", id = "original"))
print(attr(len, "zero_length"))
prof <- rsem_length_profile(len)
cat("genes with length 0 in some but not all samples:", sum(prof$zero_fraction > 0 & prof$zero_fraction < 1), "\n")
cat("genes with length 0 in every sample:", sum(prof$zero_fraction == 1), "\n")
expressed <- prof$zero_fraction == 0
cat("median CV of the length over samples (genes never 0):", signif(stats::median(prof$cv[expressed]), 3), "\n")
cat("genes (never 0) with max/min length ratio > 2:", sum(prof$max_min_ratio[expressed] > 2, na.rm = TRUE), "\n")
len_sym <- rsem_read(data_dir, metadata, species = "human", expr = "length")
cat("rows with symbols:", nrow(len_sym), "(same rows as TPM:", identical(rownames(len_sym), rownames(tpm)), ")\n")

cat("\n== rsem_detection_curve()\n")
cur <- tm(rsem_detection_curve(tpm, "tpm", thresholds = c(0.1, 0.5, 1, 2, 5, 10), min_samples = c(1, 2, 6)))
print(cur, row.names = FALSE)
k <- cur$genes[cur$threshold == 1 & cur$samples_needed == 2]
f <- suppressMessages(rsem_filter(tpm, "tpm", min_value = 1, min_samples = 2))
cat("TPM >= 1 in >= 2 samples:", k, "genes; rsem_filter keeps", nrow(f), "\n")

cat("\n== Reference z-scores (control group)\n")
ref <- tm(rsem_fit_reference(f, "CTRL", metadata = metadata, group = "group"))
z <- rsem_apply_reference(f, ref, zero_sd = "drop")
zm <- as.matrix(z)
cat("control samples: mean z over genes", paste(signif(colMeans(zm[, c("CTRL_1", "CTRL_2")]), 3), collapse = ", "),
    "; each gene's two control z-scores are +-0.707:", isTRUE(all.equal(abs(zm[, "CTRL_1"]), rep(sqrt(0.5), nrow(zm)), check.attributes = FALSE)), "\n")
cat("z of the knockdown target in its knockdown: PTBP1", paste(round(zm["PTBP1", c("PTBP1_1", "PTBP1_2")], 1), collapse = ", "),
    "; SRSF1", paste(round(zm["SRSF1", c("SRSF1_1", "SRSF1_2")], 1), collapse = ", "), "\n")
z_less <- rsem_apply_reference(f[, c("CTRL_1", "CTRL_2", "PTBP1_1")], ref, zero_sd = "drop")
cat("PTBP1_1 unchanged when the other samples are left out:", identical(z_less$PTBP1_1, z$PTBP1_1), "\n")

cat("\n== TMMwsp\n")
counts <- rsem_read(data_dir, metadata, species = "human", expr = "count")
fc <- suppressMessages(rsem_filter(counts, "count"))
w <- tm(rsem_normalize(fc, "count", "tmmwsp"))
if (requireNamespace("edgeR", quietly = TRUE)) {
  y <- edgeR::calcNormFactors(edgeR::DGEList(as.matrix(fc)), method = "TMMwsp")
  ref_cpm <- edgeR::cpm(y)
  cat("largest relative difference from edgeR", as.character(utils::packageVersion("edgeR")), ":",
      signif(max(abs(as.matrix(w) - ref_cpm) / pmax(abs(ref_cpm), 1e-300)), 3), "\n")
}
w_all <- suppressMessages(rsem_normalize(counts, "count", "tmmwsp"))
cat("on all rows (zeros included), effective library sizes:", paste(signif(attr(w_all, "factors"), 4), collapse = ", "), "\n")

cat("\n== GTF gene types (GENCODE V29)\n")
ann_time <- system.time(ann <- rsemverse:::.read_annotation(gtf))[["elapsed"]]
cat(sprintf("gene rows: %d (%.1f s)\n", nrow(ann), ann_time))
p <- rsem_biotype_profile(tpm, ann, types = c("protein_coding", "lincRNA", "antisense", "lncRNA", "Mt_rRNA", "Mt_tRNA", "rRNA", "snoRNA", "misc_RNA", "unannotated", "multiple"))
share <- stats::reshape(p[, c("sample_id", "type", "share")], idvar = "type", timevar = "sample_id", direction = "wide")
share[-1] <- round(100 * share[-1], 3)
names(share) <- sub("share.", "", names(share), fixed = TRUE)
print(share, row.names = FALSE)
pc <- rsem_biotype_profile(counts, ann, types = c("rRNA", "Mt_rRNA"))
cat("share of expected count, rRNA:", paste(round(100 * pc$share[pc$type == "rRNA"], 3), collapse = ", "),
    "; Mt_rRNA:", paste(round(100 * pc$share[pc$type == "Mt_rRNA"], 2), collapse = ", "), "\n")
cleaned <- rsem_subset(tpm, remove = c("rrna", "mito_rrna"), annotation = ann)
ribo_kept <- suppressMessages(rsem_subset(tpm, "human", keep = "ribo"))
cat("ribosomal protein genes still there after removing rRNA:", all(rownames(ribo_kept) %in% rownames(cleaned)), "\n")

cat("\n== Transcripts\n")
iso <- tm(rsem_read(data_dir, NULL, species = "human", expr = "tpm", level = "isoform"))
t2g <- attr(iso, "tx2gene")
iso_ids <- rsem_read(data_dir, NULL, species = "human", expr = "tpm", level = "isoform", id = "original")
gene_orig <- rsem_read(data_dir, metadata, species = "human", expr = "tpm", id = "original")
summed <- rowsum(as.matrix(iso_ids), attr(iso_ids, "tx2gene")$gene_id)
common <- intersect(rownames(summed), rownames(gene_orig))
cat("genes in both files:", length(common), "of", nrow(gene_orig),
    "; largest |sum of transcript TPM - gene TPM|:",
    signif(max(abs(summed[common, ] - as.matrix(gene_orig)[common, colnames(summed)])), 3), "\n")
for (r in list(list(rule = "any"), list(rule = "all"), list(min_fraction = 0.1), list(highest_only = TRUE),
               list(highest_only = TRUE, rule = "mean"))) {
  s <- suppressMessages(do.call(rsem_select_isoforms, c(list(iso), r)))
  cat(sprintf("  %-40s %s transcripts of %s genes\n", paste(names(r), unlist(r), sep = "=", collapse = ", "),
              format(nrow(s), big.mark = ","), format(length(unique(attr(s, "tx2gene")$gene_id)), big.mark = ",")))
}

cat("\n== rsem_write()\n")
out <- file.path(work, "tpm-annotated.tsv.gz")
tm(rsem_write(tpm, out, annotation = ann))
back <- utils::read.delim(gzfile(out), check.names = FALSE, nrows = 3)
cat("columns:", paste(names(back), collapse = ", "), "\n")
back_all <- data.table::fread(out, data.table = FALSE)
cat("values unchanged:", isTRUE(all.equal(as.matrix(back_all[, colnames(tpm)]), as.matrix(tpm), check.attributes = FALSE)),
    "; rows with several gene IDs:", sum(grepl(";", back_all$gene_id)), "\n")
