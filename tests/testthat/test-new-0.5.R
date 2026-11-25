# rsem_check(), reference validation, annotation names, posterior columns,
# reordered files and the log2 guard.

write_rsem5 <- function(dir, name, ids, len = 500, tpm = 1, count = tpm, tx = "T", extra = NULL) {
  tab <- data.frame(gene_id = ids, `transcript_id(s)` = tx, length = len + 100, effective_length = len,
                    expected_count = count, TPM = tpm, FPKM = tpm, check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(extra)) for (n in names(extra)) tab[[n]] <- extra[[n]]
  dir.create(dirname(file.path(dir, name)), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(tab, file.path(dir, name), sep = "\t", quote = FALSE, row.names = FALSE)
}

test_that("files from different references are refused, as sets of transcripts", {
  d <- tempfile("ref"); dir.create(d)
  ids <- c("G1", "G2", "G3")
  write_rsem5(d, "S1.genes.results", ids, tx = c("T1,T2", "T3", "T4"))
  write_rsem5(d, "S2.genes.results", ids, tx = c("T2,T1", "T3", "T4"))   # same sets, other order inside
  x <- suppressMessages(rsem_read(d, NULL, expr = "tpm", id = "original"))
  expect_equal(colnames(x), c("S1", "S2"))
  write_rsem5(d, "S3.genes.results", ids, tx = c("T1", "T3", "T5"))
  expect_error(rsem_read(d, NULL, expr = "tpm", id = "original"), "different transcript_id\\(s\\) for 2 genes.*G1")
  y <- suppressMessages(rsem_read(d, NULL, expr = "tpm", id = "original", validate_reference = FALSE))
  expect_equal(ncol(y), 3L)
  expect_error(rsem_read(d, NULL, expr = "tpm", id = "original", validate_reference = NA), "TRUE or FALSE")
})

test_that("allow_reorder accepts the same IDs in another order, and only then", {
  d <- tempfile("order"); dir.create(d)
  write_rsem5(d, "A.genes.results", c("G1", "G2", "G3"), tpm = c(1, 2, 3), tx = c("T1", "T2", "T3"))
  write_rsem5(d, "B.genes.results", c("G3", "G1", "G2"), tpm = c(30, 10, 20), tx = c("T3", "T1", "T2"))
  expect_error(rsem_read(d, NULL, expr = "tpm", id = "original"), "different order; allow_reorder = TRUE")
  x <- suppressMessages(rsem_read(d, NULL, expr = "tpm", id = "original", allow_reorder = TRUE))
  expect_equal(rownames(x), c("G1", "G2", "G3"))
  expect_equal(x$B, c(10, 20, 30))
  write_rsem5(d, "C.genes.results", c("G3", "G1", "GX"), tpm = 1, tx = "T")
  expect_error(rsem_read(d, NULL, expr = "tpm", id = "original", allow_reorder = TRUE), "same genes in the same order")
})

test_that("names come from an annotation, also for IDs with appended names and underscores", {
  ann <- data.frame(
    gene_id = c("gene1000", "gene_with_us", "ENSG00000141510", "gene2000"),
    gene_name = c("TP53", "ODD_NAME", "TP53", "GAPDH"), stringsAsFactors = FALSE
  )
  d <- tempfile("ann"); dir.create(d)
  ids <- c("gene1000_TP53", "gene_with_us_ODD_NAME", "ENSG00000141510.18_TP53", "gene2000", "other_gene", "gene3000",
           "gene2000_2", "gene1000_WRONG", "gene2000_GAPDH_x")
  write_rsem5(d, "S1.genes.results", ids, tpm = 1:9)
  write_rsem5(d, "S2.genes.results", ids, tpm = 10 * (1:9))
  x <- suppressMessages(rsem_read(d, NULL, expr = "tpm", id = "annotation", annotation = ann))
  # gene2000_2, gene1000_WRONG and gene2000_GAPDH_x are not <gene_id>_<gene_name> of a gene: they keep their IDs
  expect_equal(rownames(x), c("TP53", "ODD_NAME", "GAPDH", "other_gene", "gene3000", "gene2000_2", "gene1000_WRONG", "gene2000_GAPDH_x"))
  expect_equal(x$S1, c(1 + 3, 2, 4, 5, 6, 7, 8, 9))         # rows with the same name are summed
  g <- attr(x, "gene_ids")
  expect_equal(g$gene_id, ids)
  expect_equal(g$gene[g$gene_id == "gene_with_us_ODD_NAME"], "ODD_NAME")
  expect_error(rsem_read(d, NULL, expr = "tpm", id = "annotation"), "needs `annotation`")
  expect_error(rsem_read(d, NULL, species = "human", expr = "tpm", annotation = ann), "id = \"annotation\" only")
  # an ID is never cut where the annotation does not have the shortened ID
  expect_equal(rsemverse:::.names_from_annotation(c("gene_with_us", "gene_with", "gene1000_TP53", "gene1000_tp53"), ann),
               c("ODD_NAME", NA, "TP53", NA))
  # a duplicated gene_id in the annotation and a gene without a name do not break the lookup
  ann2 <- rbind(ann, data.frame(gene_id = c("gene1000", "gene4000"), gene_name = c("OTHER", NA)))
  expect_equal(rsemverse:::.names_from_annotation(c("gene1000", "gene4000", "gene4000_X"), ann2), c("TP53", NA, NA))
  # transcript level
  e <- tempfile("anniso"); dir.create(e)
  tab <- data.frame(transcript_id = c("T1", "T2"), gene_id = c("gene1000_TP53", "gene2000"), length = 100, effective_length = 90,
                    expected_count = 1, TPM = c(1, 2), FPKM = 1, IsoPct = 100)
  utils::write.table(tab, file.path(e, "S1.isoforms.results"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(tab, file.path(e, "S2.isoforms.results"), sep = "\t", quote = FALSE, row.names = FALSE)
  iso <- suppressMessages(rsem_read(e, NULL, expr = "tpm", level = "isoform", id = "annotation", annotation = ann))
  expect_equal(attr(iso, "tx2gene")$gene, c("TP53", "GAPDH"))
})

test_that("posterior columns are read, and the non-additive ones are never merged", {
  d <- tempfile("pme"); dir.create(d)
  ids <- c("G1", "G2")
  for (s in c("S1", "S2")) {
    write_rsem5(d, paste0(s, ".genes.results"), ids, tpm = c(1, 2), tx = c("T1", "T2"), extra = list(
      posterior_mean_count = c(5, 6), posterior_standard_deviation_of_count = c(0.5, 0.6), pme_TPM = c(1.1, 2.1),
      pme_FPKM = c(1, 2), TPM_ci_lower_bound = c(0.8, 1.7), TPM_ci_upper_bound = c(1.4, 2.5)
    ))
  }
  rd <- function(expr, ...) suppressMessages(rsem_read(d, NULL, expr = expr, id = "original", ...))
  expect_equal(rd("pme_count")$S1, c(5, 6))
  expect_equal(rd("pme_tpm")$S2, c(1.1, 2.1))
  expect_equal(rd("count_sd")$S1, c(0.5, 0.6))
  expect_equal(rd("tpm_ci_lower")$S1, c(0.8, 1.7))
  expect_equal(rd("tpm_ci_upper")$S1, c(1.4, 2.5))
  expect_error(rsem_read(d, NULL, species = "human", expr = "tpm_ci_lower"), "cannot be summed")
  expect_error(rsem_read(d, NULL, expr = "isopct", id = "original"), "must be one of")
  skip_if_not_installed("org.Hs.eg.db")
  # posterior means add up like counts when symbols are merged
  e <- tempfile("pme2"); dir.create(e)
  ids2 <- c("ENSG00000094914.1", "ENSG00000291836.1")
  for (s in c("S1", "S2")) {
    write_rsem5(e, paste0(s, ".genes.results"), ids2, tpm = c(1, 2), tx = c("T1", "T2"), extra = list(posterior_mean_count = c(5, 6)))
  }
  sym <- rsemverse:::.lookup_symbols(rsemverse:::.strip_ensembl_version(ids2), "org.Hs.eg.db")
  x <- suppressMessages(rsem_read(e, NULL, species = "human", expr = "pme_count"))
  expect_equal(sum(x$S1), 11)
  # a file without posterior columns says which option to use
  plain <- tempfile("plain"); dir.create(plain)
  write_rsem5(plain, "S1.genes.results", ids, tx = c("T1", "T2"))
  write_rsem5(plain, "S2.genes.results", ids, tx = c("T1", "T2"))
  expect_error(rsem_read(plain, NULL, expr = "pme_tpm", id = "original"), "--calc-pme or --calc-ci")
  expect_error(rsem_read(plain, NULL, expr = "tpm_ci_upper", id = "original"), "--calc-ci")
})

test_that("rsem_sample_qc refuses a second log2 and the history survives the filters", {
  set.seed(2)
  m <- as.data.frame(matrix(round(stats::rexp(80 * 5, 0.1), 1), 80, 5, dimnames = list(paste0("g", 1:80), paste0("S", 1:5))))
  lg <- suppressMessages(rsem_transform(m, "log2"))
  expect_error(rsem_sample_qc(lg, transform = "log2"), "already log2-transformed")
  expect_error(rsem_sample_qc(suppressMessages(rsem_filter(lg, "tpm", method = "nonzero")), transform = "log2"), "already log2-transformed")
  qc <- suppressMessages(rsem_sample_qc(lg, transform = "none", n_top = 30))
  ref <- suppressMessages(rsem_sample_qc(m, n_top = 30))
  expect_equal(qc$pca_variance$variance, ref$pca_variance$variance)
  expect_equal(qc$correlation_matrix, ref$correlation_matrix)
  # z-scores against a reference keep the history, so they are not logged again either
  z <- suppressMessages(rsem_apply_reference(m, suppressMessages(rsem_fit_reference(m, c("S1", "S2", "S3")))))
  expect_equal(attr(z, "transform"), c("log2", "zscore"))
  expect_error(rsem_sample_qc(z, transform = "log2"), "already log2-transformed")
  expect_equal(suppressMessages(rsem_sample_qc(z, n_top = 30))$parameters$transform, "none")  # the default sees it
})

test_that("rsem_check reports every problem at once", {
  d <- tempfile("chk"); dir.create(d)
  ids <- c("G1", "G2", "G3", "G4")
  tpm_ok <- c(250000, 250000, 250000, 250000)
  write_rsem5(d, "S1.genes.results", ids, len = c(500, 0, 300, 200), tpm = tpm_ok, tx = c("T1", "T2", "T3", "T4"))
  write_rsem5(d, "S2.genes.results", ids[c(2, 1, 3, 4)], len = 500, tpm = tpm_ok, tx = c("T2", "T1", "T3", "T4"))
  write_rsem5(d, "S3.genes.results", ids, len = 500, tpm = c(1, 1, 1, 1), tx = c("T1", "T2", "T9", "T4"))
  write_rsem5(d, "S4.genes.results", c("G1", "G2", "G3", "GX"), len = 500, tpm = tpm_ok, tx = c("T1", "T2", "T3", "T4"))
  write_rsem5(d, "S5.genes.results", ids, len = 500, tpm = c(-1, 250000, 500000, 250001), tx = c("T1", "T2", "T3", "T4"))
  writeLines("gene_id\tx", file.path(d, "S6.genes.results"))
  chk <- suppressMessages(rsem_check(d))
  has <- function(sev, check, sid) any(chk$severity == sev & chk$check == check & chk$sample_id %in% sid)
  expect_true(has("WARNING", "zero_effective_length", "S1"))
  expect_equal(chk$n[chk$check == "zero_effective_length"], 1L)
  expect_true(has("INFO", "id_order", "S2"))
  expect_true(has("WARNING", "tpm_total", "S3"))
  expect_true(has("ERROR", "reference_mismatch", "S3"))
  expect_true(has("ERROR", "id_set", "S4"))
  expect_true(has("ERROR", "negative_values", "S5"))
  expect_true(has("WARNING", "required_columns", "S6"))
  expect_equal(chk$severity, chk$severity[order(match(chk$severity, c("ERROR", "WARNING", "INFO")))])
  expect_message(rsem_check(d), "6 files; 5 errors, 3 warnings, 1 notes")
  # a clean directory and the bundled example
  ex <- suppressMessages(rsem_check(file.path(example_dir(), "data")))
  expect_false(any(ex$severity == "ERROR"))
  clean <- tempfile("clean"); dir.create(clean)
  write_rsem5(clean, "A.genes.results", ids, len = 500, tpm = tpm_ok, count = 5, tx = c("T1", "T2", "T3", "T4"))
  write_rsem5(clean, "B.genes.results", ids, len = 500, tpm = tpm_ok, count = 7, tx = c("T1", "T2", "T3", "T4"))
  none <- suppressMessages(rsem_check(clean))
  expect_equal(nrow(none), 0L)
  expect_named(none, c("severity", "check", "sample_id", "n", "detail"))
  # problems with the files themselves
  expect_equal(suppressMessages(rsem_check(tempfile("none")))$check, "files")
  dup <- tempfile("dup"); dir.create(dup)
  write_rsem5(dup, "A.genes.results", c("G1", "G1"), tpm = 1)
  write_rsem5(dup, "B.genes.results", c("G1", "G1"), tpm = 1)
  expect_true(any(suppressMessages(rsem_check(dup))$check == "duplicated_ids"))
  # isoform level: a different gene for a transcript
  iso <- tempfile("iso"); dir.create(iso)
  mk <- function(s, genes) {
    tab <- data.frame(transcript_id = c("T1", "T2"), gene_id = genes, length = 100, effective_length = 90,
                      expected_count = 5, TPM = c(5e5, 5e5), FPKM = 1, IsoPct = 50)
    utils::write.table(tab, file.path(iso, paste0(s, ".isoforms.results")), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  mk("S1", c("G1", "G1")); mk("S2", c("G1", "G2"))
  ci <- suppressMessages(rsem_check(iso, level = "isoform"))
  expect_equal(ci$check, "reference_mismatch")
  expect_equal(ci$sample_id, "S2")
})

test_that("rsem_check boundaries: transcript order inside a gene, TPM total tolerance, zero length with counts", {
  ids <- c("G1", "G2", "G3", "G4")
  tpm <- c(250000, 250000, 250000, 250000)
  d <- tempfile("edge"); dir.create(d)
  write_rsem5(d, "A.genes.results", ids, tpm = tpm, tx = c("T1,T2", "T3", "T4", "T5,T6,T7"))
  # the same sets listed in another order inside each gene are not a difference
  write_rsem5(d, "B.genes.results", ids, tpm = tpm, tx = c("T2,T1", "T3", "T4", "T7,T5,T6"))
  expect_equal(nrow(suppressMessages(rsem_check(d))), 0L)
  # a different set is one, also when only one transcript is added or dropped
  write_rsem5(d, "C.genes.results", ids, tpm = tpm, tx = c("T1,T2", "T3", "T4", "T5,T6"))
  chk <- suppressMessages(rsem_check(d))
  expect_equal(chk$check, "reference_mismatch")
  expect_equal(chk$n, 1L)
  expect_equal(chk$sample_id, "C")
  # TPM total: rounding in the files (a few units) is accepted, 0.2% off is a warning
  t <- tempfile("tpm"); dir.create(t)
  write_rsem5(t, "A.genes.results", ids, tpm = tpm + c(0.4, -0.3, 0.2, 0.1), tx = c("T1", "T2", "T3", "T4"))
  write_rsem5(t, "B.genes.results", ids, tpm = tpm + c(900, 0, 0, 0), tx = c("T1", "T2", "T3", "T4"))
  write_rsem5(t, "C.genes.results", ids, tpm = tpm + c(2000, 0, 0, 0), tx = c("T1", "T2", "T3", "T4"))
  chk <- suppressMessages(rsem_check(t))
  expect_equal(chk$sample_id[chk$check == "tpm_total"], "C")
  # a gene with effective length 0 and a positive count is reported as such and never changed
  z <- tempfile("zero"); dir.create(z)
  write_rsem5(z, "A.genes.results", ids, len = c(500, 0, 300, 200), tpm = tpm, count = c(10, 5, 7, 3), tx = c("T1", "T2", "T3", "T4"))
  write_rsem5(z, "B.genes.results", ids, len = c(500, 0, 300, 200), tpm = tpm, count = c(11, 6, 8, 4), tx = c("T1", "T2", "T3", "T4"))
  chk <- suppressMessages(rsem_check(z))
  expect_equal(chk$check, c("zero_effective_length", "zero_effective_length"))
  expect_true(all(chk$severity == "WARNING"))
  len <- suppressMessages(rsem_read(z, NULL, expr = "length", id = "original"))
  cnt <- suppressMessages(rsem_read(z, NULL, expr = "count", id = "original"))
  expect_equal(len$A, c(500, 0, 300, 200))
  expect_equal(cnt$A, c(10, 5, 7, 3))
  expect_equal(attr(len, "zero_length"), c(A = 1, B = 1))
})
