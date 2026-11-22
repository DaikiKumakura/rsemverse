# Validation of rsem_reconcile() and rsem_id_report() on six ENCODE K562 RSEM
# runs (gene and transcript files; GENCODE V29). Run from the package root:
#   Rscript validation/encode-k562-0.6.0.R <folder with *.genes.results and *.isoforms.results>
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, dir.exists(args[[1L]]))
suppressMessages(devtools::load_all(quiet = TRUE))
dir <- args[[1L]]

rec <- rsem_reconcile(dir)
print(rec[, c("sample_id", "column", "n_genes", "n_beyond_rounding", "max_abs_diff", "genes_without_isoforms", "transcript_mapping_mismatch")])
stopifnot(all(rec$n_beyond_rounding == 0), all(rec$genes_without_isoforms == 0), all(rec$isoform_genes_missing == 0),
          all(rec$transcript_mapping_mismatch == 0))

# A different sample's isoform file must be caught.
tmp <- tempfile(); dir.create(tmp)
files <- list.files(dir, pattern = "results$", full.names = TRUE)
ok <- file.copy(files, tmp)
ids <- sub("\\.genes\\.results$", "", basename(list.files(dir, pattern = "genes.results$")))
swap <- file.path(tmp, paste0(ids[[1L]], ".isoforms.results"))
file.copy(file.path(dir, paste0(ids[[2L]], ".isoforms.results")), swap, overwrite = TRUE)
bad <- suppressMessages(rsem_reconcile(tmp, columns = "TPM"))
cat("\nIsoform file of sample 2 under the name of sample 1:\n")
print(bad[, c("sample_id", "n_genes", "n_beyond_rounding", "max_abs_diff")])
stopifnot(bad$n_beyond_rounding[bad$sample_id == ids[[1L]]] > 0, all(bad$n_beyond_rounding[bad$sample_id != ids[[1L]]] == 0))

chk <- rsem_check(dir)
cat("\nrsem_check findings:", nrow(chk), "\n")

x <- suppressMessages(rsem_read(dir, metadata = NULL, species = "human", expr = "tpm"))
rep <- rsem_id_report(x)
cat("\nrows that combine several gene IDs:", nrow(rep), "\n")
print(utils::head(rep, 5L))
