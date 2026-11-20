# Real-data check of the changes in rsemverse 0.4.0, on the ENCODE K562 files
# that encode-k562-0.3.R downloads (run that script first, with the same
# work directory).
#
#   Rscript validation/encode-k562-0.4.0.R my-work-folder
#
# Needs an installed rsemverse, org.Hs.eg.db and edgeR.

suppressMessages(library(rsemverse))
W <- commandArgs(TRUE)[1]
stopifnot(!is.na(W), dir.exists(file.path(W, "data")))
meta <- file.path(W, "metadata.csv"); data <- file.path(W, "data")
tpm <- suppressMessages(rsem_read(data, meta, species = "human", expr = "tpm"))
f <- suppressMessages(rsem_filter(tpm, "tpm"))
g <- attr(f, "gene_ids")
cat("rows kept", nrow(f), "; gene_ids rows", nrow(g), "; all kept symbols covered:", all(rownames(f) %in% g$gene), "\n")
out <- tempfile(fileext = ".tsv.gz"); suppressMessages(rsem_write(f, out))
b <- data.table::fread(out, data.table = FALSE); cat("columns after filter:", paste(names(b)[1:3], collapse = ", "), "; rows with several IDs:", sum(grepl(";", b$gene_id)), "\n")
# gzip copies of the six gene files read to the same table
gz <- tempfile("gz"); dir.create(gz)
for (r in list.files(data, "genes.results$")) { con <- gzfile(file.path(gz, paste0(r, ".gz")), "w"); writeLines(readLines(file.path(data, r)), con); close(con) }
t0 <- system.time(z <- suppressMessages(rsem_read(gz, NULL, species = "human", expr = "tpm")))[["elapsed"]]
z <- z[, colnames(tpm)]
cat("gzip input identical to plain input:", isTRUE(all.equal(as.matrix(z), as.matrix(tpm))), sprintf("(%.1f s)\n", t0))
# group-aware edgeR filter on counts
counts <- suppressMessages(rsem_read(data, meta, species = "human", expr = "count"))
grp <- c("ctrl", "ctrl", "ptbp1", "ptbp1", "srsf1", "srsf1")
k1 <- suppressMessages(rsem_filter(counts, "count")); k2 <- suppressMessages(rsem_filter(counts, "count", group = grp))
e <- edgeR::filterByExpr(as.matrix(counts), group = grp)
cat("edger rule: all-in-one", nrow(k1), "genes; with group", nrow(k2), "; equals edgeR filterByExpr(group):", identical(rownames(k2), rownames(counts)[e]), "\n")
# biotype share against the unfiltered total
ann <- rsemverse:::.read_annotation(file.path(W, "gencode.v29.annotation.gtf.gz"))
p_all <- suppressMessages(rsem_biotype_profile(tpm, ann, types = "Mt_rRNA"))
p_flt <- suppressMessages(rsem_biotype_profile(f, ann, types = "Mt_rRNA"))
p_tot <- suppressMessages(rsem_biotype_profile(f, ann, types = "Mt_rRNA", total = colSums(tpm)))
cat("Mt_rRNA share of TPM (CTRL_1): unfiltered", round(100*p_all$share[1], 3), "| filtered matrix as denominator", round(100*p_flt$share[1], 3), "| filtered with total of unfiltered", round(100*p_tot$share[1], 3), "\n")
# isoform tables keep tx2gene through a filter
iso <- suppressMessages(rsem_read(data, NULL, species = "human", expr = "tpm", level = "isoform"))
fi <- suppressMessages(rsem_filter(iso, "tpm", min_value = 1, min_samples = 2))
cat("isoforms kept", nrow(fi), "of", nrow(iso), "; tx2gene rows", nrow(attr(fi, "tx2gene")), "; all kept transcripts covered:", all(rownames(fi) %in% attr(fi, "tx2gene")$transcript_id), "\n")

# IDs in the style of --append-names, made from the real files (the symbol is
# appended to the genes that have one); the file is read back to the same table
ids <- attr(tpm, "gene_ids")
ap <- tempfile("an"); dir.create(ap)
for (r in list.files(data, "genes.results$")) {
  d <- utils::read.delim(file.path(data, r), check.names = FALSE, colClasses = c(gene_id = "character"))
  i <- match(d$gene_id, ids$gene_id)
  d$gene_id <- ifelse(ids$gene[i] != ids$gene_id[i], paste0(d$gene_id, "_", ids$gene[i]), d$gene_id)
  utils::write.table(d, file.path(ap, r), sep = "	", quote = FALSE, row.names = FALSE)
}
a <- suppressMessages(rsem_read(ap, NULL, species = "human", expr = "tpm"))[, colnames(tpm)]
cat("append-names style IDs: same rows:", setequal(rownames(a), rownames(tpm)),
    "; same values:", isTRUE(all.equal(as.matrix(a[rownames(tpm), ]), as.matrix(tpm))), "
")

# sample QC on the filtered table: equal to prcomp() and cor(), and the
# replicates of each condition are each other's nearest samples
qc <- suppressMessages(rsem_sample_qc(f, metadata = meta, n_top = 2000))
l <- log2(as.matrix(f) + 1)
v <- apply(l, 1, stats::var)
top <- sort(order(-v)[1:2000])
p <- stats::prcomp(t(l[top, ]), center = TRUE, scale. = FALSE)
cat("sample QC: PC1-PC5 variance share", paste(round(100 * qc$pca_variance$proportion[1:5], 1), collapse = ", "), "%;",
    "scores equal prcomp:", isTRUE(all.equal(unname(as.matrix(qc$pca_scores[, paste0("PC", 1:5)])), unname(p$x[, 1:5]))),
    "; correlation equal cor():", isTRUE(all.equal(qc$correlation_matrix, stats::cor(l[rowSums(as.matrix(f) >= 1) >= 1, ], method = "spearman"))), "\n")
print(qc$sample_summary[, c("sample_id", "detected_genes", "median_correlation", "most_similar")], row.names = FALSE)
cat("most similar sample is the replicate:", all(sub("_[0-9]+$", "", qc$sample_summary$most_similar) == sub("_[0-9]+$", "", qc$sample_summary$sample_id)), "\n")
