# Generates the synthetic example in inst/extdata/example/.
#
# The gene IDs are real human Ensembl gene IDs (taken from org.Hs.eg.db); the
# expression values are synthetic. The example also contains
#   * a pair of Ensembl IDs that share one gene symbol,
#   * a spike-in and a numeric ID that have no symbol.
#
#   Rscript data-raw/make_example.R

suppressMessages(library(org.Hs.eg.db))
set.seed(20261110)

symbols <- c(
  "TP53", "GAPDH", "ACTB", "MYC", "EGFR", "BRCA1", "BRCA2", "KRAS", "PTEN",
  "CDKN2A", "RB1", "CD3E", "CD4", "CD8A", "IL2", "IFNG", "TNF", "VEGFA",
  "STAT3", "JAK2", "MTOR", "AKT1", "PIK3CA", "BRAF", "NRAS", "ERBB2", "ALB",
  "INS", "HBB", "CD19", "MS4A1", "PDCD1", "CD274", "CTLA4", "LAG3", "TIGIT",
  "FOXP3", "GZMB", "PRF1", "IL6"
)
single <- suppressMessages(AnnotationDbi::mapIds(
  org.Hs.eg.db, keys = symbols, keytype = "SYMBOL", column = "ENSEMBL", multiVals = "first"
))
single <- single[!is.na(single)]

# One symbol that is used by exactly two Ensembl IDs (first in alphabetical order).
all_map <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys = keys(org.Hs.eg.db, "ENSEMBL"), keytype = "ENSEMBL", columns = "SYMBOL"
))
all_map <- unique(all_map[!is.na(all_map$SYMBOL), ])
n_ids <- table(all_map$SYMBOL)
shared_symbol <- sort(names(n_ids)[n_ids == 2L])[1]
shared <- all_map$ENSEMBL[all_map$SYMBOL == shared_symbol]
cat("shared symbol:", shared_symbol, "->", paste(shared, collapse = ", "), "\n")

stable <- c(unname(single), shared)
gene_id <- c(paste0(stable, ".", sample(1:19, length(stable), replace = TRUE)), "ERCC-00002", "10904")

n_samples <- 3L
base <- c(stats::rgamma(length(gene_id), shape = 1.5, rate = 0.01) + 5)
counts <- vapply(seq_len(n_samples), function(i) {
  stats::rnbinom(length(gene_id), mu = base * runif(1, 0.8, 1.2), size = 20)
}, numeric(length(gene_id)))
counts <- counts + 0.5 * (counts > 0)   # RSEM expected counts are not integers
effective_length <- round(stats::runif(length(gene_id), 500, 5000)) + 0.5

out_dir <- file.path("inst", "extdata", "example")
dir.create(file.path(out_dir, "data"), recursive = TRUE, showWarnings = FALSE)
sample_id <- paste0("S", seq_len(n_samples))
for (i in seq_len(n_samples)) {
  rate <- counts[, i] / effective_length
  tpm <- rate / sum(rate) * 1e6
  tab <- data.frame(
    gene_id = gene_id,
    `transcript_id(s)` = sub("G", "T", gene_id),
    length = effective_length + 150,
    effective_length = effective_length,
    expected_count = counts[, i],
    TPM = round(tpm, 4),
    FPKM = round(tpm * 0.9, 4),
    check.names = FALSE
  )
  utils::write.table(
    tab, file.path(out_dir, "data", paste0(sample_id[i], ".genes.results")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}
writeLines(
  c("sample_id,ref", paste0(sample_id, ",", sample_id, ".genes.results")),
  file.path(out_dir, "metadata.csv")
)
cat("wrote", length(gene_id), "genes x", n_samples, "samples to", out_dir, "\n")
