# Real-data check of the changes in rsemverse 0.5.0, on the ENCODE K562 files
# that encode-k562-0.3.R downloaded (run that script first, with the same
# work directory).
#
#   Rscript validation/encode-k562-0.5.0.R my-work-folder
#
# Needs an installed rsemverse and org.Hs.eg.db.

suppressMessages(library(rsemverse)); W <- commandArgs(TRUE)[1]
data <- file.path(W, "data"); meta <- file.path(W, "metadata.csv")
t0 <- system.time(chk <- suppressMessages(rsem_check(data)))[["elapsed"]]
cat(sprintf("rsem_check on the six ENCODE gene files (%.1f s): %d findings\n", t0, nrow(chk)))
print(chk[, c("severity", "check", "sample_id", "n")], row.names = FALSE)
tp <- suppressMessages(rsem_check(data, level = "isoform"))
cat("isoform files:", nrow(tp), "findings;", paste(unique(tp$check), collapse = ", "), "\n")
# validate_reference on the real files: the transcript lists are equal in all six
x <- suppressMessages(rsem_read(data, meta, species = "human", expr = "tpm", validate_reference = TRUE))
cat("validate_reference = TRUE accepted the six files:", ncol(x) == 6, "\n")
# a file from another reference: change the transcript list of one gene in one sample
tmp <- tempfile("mix"); dir.create(tmp)
fs <- list.files(data, "genes.results$")
for (f in fs) file.copy(file.path(data, f), file.path(tmp, f))
d <- data.table::fread(file.path(tmp, fs[2]), data.table = FALSE, colClasses = "character", sep = "\t")
d[["transcript_id(s)"]][5] <- paste0(d[["transcript_id(s)"]][5], ",ENSTFAKE")
data.table::fwrite(d, file.path(tmp, fs[2]), sep = "\t", quote = FALSE)
cat("one gene with another transcript list: ", tryCatch({suppressMessages(rsem_read(tmp, NULL, species = "human", expr = "tpm")); "read"}, error = function(e) conditionMessage(e)), "\n", sep = "")
print(suppressMessages(rsem_check(tmp))[, c("severity", "check", "sample_id", "n")], row.names = FALSE)
# shuffled row order in one file
sh <- tempfile("shuf"); dir.create(sh)
for (f in fs) file.copy(file.path(data, f), file.path(sh, f))
d <- data.table::fread(file.path(sh, fs[3]), data.table = FALSE, colClasses = list(character = c("gene_id", "transcript_id(s)")), sep = "\t")
set.seed(1); d <- d[sample(nrow(d)), ]
data.table::fwrite(d, file.path(sh, fs[3]), sep = "\t", quote = FALSE)
cat("shuffled file, default: ", tryCatch({suppressMessages(rsem_read(sh, NULL, species = "human", expr = "tpm")); "read"}, error = function(e) "stopped"), "; allow_reorder = TRUE: ",
    isTRUE(all.equal(as.matrix(suppressMessages(rsem_read(sh, NULL, species = "human", expr = "tpm", allow_reorder = TRUE))),
                     as.matrix(suppressMessages(rsem_read(data, NULL, species = "human", expr = "tpm"))))), "\n", sep = "")
# name from the GENCODE annotation instead of the OrgDb
ann <- rsemverse:::.read_annotation(file.path(W, "gencode.v29.annotation.gtf.gz"), need_type = FALSE)
a <- suppressMessages(rsem_read(data, NULL, expr = "tpm", id = "annotation", annotation = ann))
o <- suppressMessages(rsem_read(data, NULL, species = "human", expr = "tpm"))
ids <- attr(a, "gene_ids")
cat("names from GENCODE: rows", nrow(a), "(OrgDb:", nrow(o), "); genes named:", sum(ids$gene != ids$gene_id), "of", nrow(ids),
    "; column sums equal:", isTRUE(all.equal(colSums(a), colSums(o))), "\n")
common <- intersect(rownames(a), rownames(o))
cat("rows with the same name in both:", length(common), "; same values:", isTRUE(all.equal(as.matrix(a[common, ]), as.matrix(o[common, ]))), "\n")

# names from GENCODE against the OrgDb: where the values differ
ann <- rsemverse:::.read_annotation(file.path(W, "gencode.v29.annotation.gtf.gz"), need_type = FALSE)
a <- suppressMessages(rsem_read(data, NULL, expr = "tpm", id = "annotation", annotation = ann))
o <- suppressMessages(rsem_read(data, NULL, species = "human", expr = "tpm"))
ga <- attr(a, "gene_ids"); go <- attr(o, "gene_ids")
common <- intersect(rownames(a), rownames(o))
d <- common[abs(a[common, "CTRL_1"] - o[common, "CTRL_1"]) > 1e-9]
na <- table(ga$gene)[d]; no <- table(go$gene)[d]
cat("common names", length(common), "; differing", length(d), "\n")
cat("differing rows where the annotation merged >1 ID:", sum(na > 1), "; where the OrgDb merged >1 ID:", sum(no > 1), "; neither:", sum(na == 1 & no == 1), "\n")
cat("annotation rows merging >1 ID in total:", sum(table(ga$gene) > 1), " (of which PAR_Y copies:", sum(grepl("_PAR_Y$", ga$gene_id)), "IDs)\n")
cat("example differing, single in both:", head(d[na == 1 & no == 1]), "\n")
