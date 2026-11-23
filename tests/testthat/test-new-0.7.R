tab <- "\t"
write_genes <- function(dir, sid, df) {
  utils::write.table(df, file.path(dir, paste0(sid, ".genes.results")), sep = tab, quote = FALSE, row.names = FALSE)
}
write_iso <- function(dir, sid, df) {
  utils::write.table(df, file.path(dir, paste0(sid, ".isoforms.results")), sep = tab, quote = FALSE, row.names = FALSE)
}
genes_df <- function(count = c(10, 5), tpm = c(600000, 400000), fpkm = c(6, 4), len = c(900, 800), extra = list()) {
  g <- data.frame(
    gene_id = c("G1", "G2"), "transcript_id(s)" = c("T1,T2", "T3"), length = 1000, effective_length = len,
    expected_count = count, TPM = tpm, FPKM = fpkm, check.names = FALSE, stringsAsFactors = FALSE
  )
  for (nm in names(extra)) g[[nm]] <- extra[[nm]]
  g
}
iso_df <- function(len = c(1000, 1500, 800)) {
  data.frame(
    transcript_id = c("T1", "T2", "T3"), gene_id = c("G1", "G1", "G2"), length = len, effective_length = 700,
    expected_count = c(4, 6, 5), TPM = c(250000, 350000, 400000), FPKM = c(2, 3, 4), IsoPct = c(40, 60, 100),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

test_that("rsem_read reads FPKM and the posterior columns RSEM writes with --calc-ci", {
  d <- withr::local_tempdir()
  ex <- list(
    posterior_mean_count = c(9.5, 5.2), pme_TPM = c(590000, 410000), pme_FPKM = c(5.9, 4.1),
    TPM_ci_lower_bound = c(580000, 400000), TPM_ci_upper_bound = c(610000, 420000),
    FPKM_ci_lower_bound = c(5.8, 4), FPKM_ci_upper_bound = c(6.1, 4.2),
    TPM_coefficient_of_quartile_variation = c(0.01, NaN), FPKM_coefficient_of_quartile_variation = c(0.02, 0.03)
  )
  write_genes(d, "s1", genes_df(extra = ex))
  write_genes(d, "s2", genes_df(extra = ex))
  rd <- function(e) suppressMessages(rsem_read(d, metadata = NULL, expr = e, id = "original"))
  expect_equal(rd("fpkm")$s1, c(6, 4))
  expect_equal(rd("pme_fpkm")$s1, c(5.9, 4.1))
  expect_equal(rd("fpkm_ci_upper")$s1, c(6.1, 4.2))
  cqv <- rd("tpm_cqv")
  expect_equal(cqv$s1[[1]], 0.01)
  expect_true(is.na(cqv$s1[[2]]))
  expect_equal(rd("fpkm_cqv")$s2, c(0.02, 0.03))
  expect_error(suppressMessages(rsem_read(d, metadata = NULL, expr = "tpm_cqv", id = "symbol", species = "human")), "cannot be summed")
  d2 <- withr::local_tempdir()
  write_genes(d2, "s1", genes_df(extra = list(posterior_mean_count = c(9.5, 5.2), pme_TPM = c(1, 2))))
  write_genes(d2, "s2", genes_df(extra = list(posterior_mean_count = c(9.5, 5.2), pme_TPM = c(1, 2))))
  expect_error(suppressMessages(rsem_read(d2, metadata = NULL, expr = "fpkm_cqv", id = "original")), "--calc-ci")
  expect_error(suppressMessages(rsem_read(d2, metadata = NULL, expr = "pme_fpkm", id = "original")), "--calc-pme or --calc-ci")
})

test_that("rsem_read reads IsoPct_from_pme_TPM for transcripts only", {
  d <- withr::local_tempdir()
  iso <- iso_df(); iso$IsoPct_from_pme_TPM <- c(41, 59, 100)
  write_iso(d, "s1", iso); write_iso(d, "s2", iso)
  v <- suppressMessages(rsem_read(d, metadata = NULL, level = "isoform", expr = "isopct_pme", id = "original"))
  expect_equal(v$s1, c(41, 59, 100))
  write_genes(d, "s1", genes_df()); write_genes(d, "s2", genes_df())
  expect_error(suppressMessages(rsem_read(d, metadata = NULL, expr = "isopct_pme", id = "original")), "must be one of")
})

test_that("rsem_check reports damaged posterior columns and transcripts with another length", {
  d <- withr::local_tempdir()
  ex <- list(TPM_ci_lower_bound = c(700000, 400000), TPM_ci_upper_bound = c(610000, 420000), pme_TPM = c(5, -1))
  write_genes(d, "s1", genes_df(extra = ex))
  res <- suppressMessages(rsem_check(d))
  expect_true("posterior_interval" %in% res$check)
  expect_true("posterior_values" %in% res$check)
  d2 <- withr::local_tempdir()
  write_iso(d2, "s1", iso_df()); write_iso(d2, "s2", iso_df(len = c(1000, 1500, 790)))
  res2 <- suppressMessages(rsem_check(d2, level = "isoform"))
  w <- res2[res2$check == "transcript_length", ]
  expect_equal(nrow(w), 1L)
  expect_equal(w$sample_id, "s2")
  expect_match(w$detail, "T3")
  g <- withr::local_tempdir()
  a <- genes_df(); b <- genes_df(); b$length <- c(1100, 900)  # gene length is a weighted mean and may vary
  write_genes(g, "s1", a); write_genes(g, "s2", b)
  expect_false("transcript_length" %in% suppressMessages(rsem_check(g))$check)
})

test_that("rsem_mapping_qc refuses impossible .cnt values and handles no alignable reads", {
  d <- withr::local_tempdir()
  w <- function(name, lines) writeLines(lines, file.path(d, paste0(name, ".cnt")))
  w("ok", c("100 850 50 1000", "800 50 120", "1500 2", "0\t100", "1\t700", "2\t150", "Inf\t50"))
  w("zero", c("100 0 0 100", "0 0 0", "0 2", "0\t100", "Inf\t0"))
  res <- suppressMessages(rsem_mapping_qc(d))
  expect_true(is.na(res$unique_gene_fraction[res$sample_id == "zero"]))
  d2 <- withr::local_tempdir()
  w2 <- function(lines) { writeLines(lines, file.path(d2, "x.cnt")); d2 }
  expect_error(suppressMessages(rsem_mapping_qc(w2(c("Inf 850 50 Inf", "800 50 120", "1500 2")))), "not an RSEM .cnt")
  expect_error(suppressMessages(rsem_mapping_qc(w2(c("100 850 50 1000", "800 50 120", "1500 7")))), "read type")
  expect_error(suppressMessages(rsem_mapping_qc(w2(c("100 850 50 1000", "800 50 900", "1500 2")))), "nUncertain")
  prof <- suppressMessages(rsem_alignment_profile(d))
  p <- prof[prof$sample_id == "ok", ]
  expect_equal(p$n_alignments, c(0, 1, 2, Inf))
  expect_equal(p$n_reads, c(100, 700, 150, 50))
  expect_equal(sum(p$fraction), 1)
  w("bad", c("100 850 50 1000", "800 50 120", "1500 2", "0\t10", "Inf\t5"))
  expect_warning(suppressMessages(rsem_alignment_profile(d)), "do not add up")
})

test_that("rsem_write writes GCT 1.2 for GSEA", {
  x <- data.frame(S1 = c(1.5, 2), S2 = c(3, NA), row.names = c("TP53", "GAPDH"))
  f <- tempfile(fileext = ".gct")
  suppressMessages(rsem_write(x, f, format = "gct"))
  lines <- readLines(f)
  expect_equal(lines[1:3], c("#1.2", "2\t2", "Name\tDescription\tS1\tS2"))
  expect_equal(lines[4], "TP53\tNA\t1.5\t3")
  expect_equal(lines[5], "GAPDH\tNA\t2\tNA")
  expect_error(suppressMessages(rsem_write(x, "x.gct.gz", format = "gct")), "plain text")
})

test_that("rsem_sample_qc gives relative log expression for log2 values", {
  set.seed(3)
  m <- matrix(rexp(300 * 6, 0.05), 300, 6, dimnames = list(paste0("g", 1:300), paste0("S", 1:6)))
  m[, "S6"] <- m[, "S6"] * 4
  qc <- suppressMessages(rsem_sample_qc(as.data.frame(m), n_top = 100))
  s <- qc$sample_summary
  expect_true(all(c("rle_median", "rle_iqr") %in% names(s)))
  expect_gt(s$rle_median[s$sample_id == "S6"], 1)
  lg <- log2(m + 1)
  dev <- lg - apply(lg, 1, median)
  expect_equal(s$rle_median[s$sample_id == "S1"], median(dev[, "S1"]), tolerance = 1e-12)
  qn <- suppressMessages(rsem_sample_qc(as.data.frame(m), n_top = 100, transform = "none"))
  expect_true(all(is.na(qn$sample_summary$rle_median)))
})

test_that("rsem_import and rsem_as_tximport give the list that DESeq2 takes", {
  d <- withr::local_tempdir()
  write_genes(d, "s1", genes_df()); write_genes(d, "s2", genes_df(count = c(11, 6)))
  a <- suppressMessages(rsem_import(d, metadata = NULL))
  expect_s3_class(a, "rsem_assays")
  expect_equal(names(a), c("count", "tpm", "length"))
  txi <- rsem_as_tximport(a)
  expect_equal(names(txi), c("abundance", "counts", "length", "countsFromAbundance"))
  expect_equal(txi$counts[, "s2"], c(G1 = 11, G2 = 6))
  expect_equal(txi$length[, "s1"], c(G1 = 900, G2 = 800))
  expect_equal(txi$countsFromAbundance, "no")
  d2 <- withr::local_tempdir()
  write_genes(d2, "s1", genes_df(len = c(900, 0))); write_genes(d2, "s2", genes_df())
  a2 <- suppressMessages(rsem_import(d2, metadata = NULL))
  expect_error(rsem_as_tximport(a2), "effective length of 0")
  dropped <- suppressMessages(rsem_as_tximport(a2, zero_length = "drop"))
  expect_equal(rownames(dropped$counts), "G1")
  expect_equal(nrow(rsem_as_tximport(a2, zero_length = "keep")$counts), 2L)
  expect_error(rsem_as_tximport(list()), "rsem_import")
  b <- suppressMessages(rsem_import(d, metadata = NULL, assays = "count"))
  expect_error(rsem_as_tximport(b), "needs the assays")
})

test_that("rsem_reconcile says when the transcripts could not be compared", {
  d <- withr::local_tempdir()
  g <- genes_df(count = c(10, 5)); g[["transcript_id(s)"]] <- NULL
  write_genes(d, "s1", g); write_iso(d, "s1", iso_df())
  res <- suppressMessages(rsem_reconcile(d, columns = "expected_count"))
  expect_false(res$transcript_mapping_checked)
  expect_message(rsem_reconcile(d, columns = "expected_count"), "were not compared")
})

test_that("a GTF without gene rows gives names from its exon rows", {
  gtf <- tempfile(fileext = ".gtf")
  writeLines(c(
    'chr1\tt\texon\t1\t10\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "ALPHA";',
    'chr1\tt\texon\t20\t30\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "ALPHA";',
    'chr1\tt\texon\t50\t80\t.\t+\t.\tgene_id "G2"; transcript_id "T2"; gene_name "BETA_X";'
  ), gtf)
  a <- rsemverse:::.read_annotation(gtf, need_type = FALSE)
  expect_equal(a$gene_id, c("G1", "G2"))
  expect_equal(a$gene_name, c("ALPHA", "BETA_X"))
})

test_that("files written by RSEM 1.3.3 itself (--calc-pme --calc-ci --append-names) are read", {
  d <- testthat::test_path("fixtures", "rsem-real")
  rd <- function(e, ...) suppressMessages(rsem_read(d, metadata = NULL, expr = e, id = "original", ...))
  expect_equal(rd("count")$sampA, c(900, 800, 50, 400))
  expect_equal(rownames(rd("count")), c("GA_ALPHA", "GB_BETA", "GC_GAMMA", "GD_DELTA_X"))
  expect_equal(rd("tpm_ci_lower")$sampA, c(410511, 290644, 0, 236636))
  expect_equal(rd("tpm_cqv")$sampA[[1]], 0.0177118)
  expect_equal(rd("fpkm_cqv")$sampA[[3]], 0)
  iso <- suppressMessages(rsem_read(d, metadata = NULL, level = "isoform", expr = "isopct_pme", id = "original"))
  expect_equal(iso$sampA, c(51.67, 48.33, 100, 0, 100))
  chk <- suppressMessages(rsem_check(d))
  expect_equal(unique(chk$severity), "WARNING")
  expect_equal(unique(chk$check), "zero_effective_length")
  expect_equal(unique(suppressMessages(rsem_reconcile(d))$n_beyond_rounding), 0L)
  m <- suppressMessages(rsem_mapping_qc(d, sample_from = "parent_dir"))
  expect_equal(m$n_total, c(2150, 2150)); expect_equal(m$n_uncertain, c(323, 323)); expect_equal(m$read_type, c(0, 0))
  p <- suppressMessages(rsem_alignment_profile(d, sample_from = "parent_dir"))
  expect_equal(p$n_reads[p$sample_id == "sampA"], c(0, 1827, 323, 0))
  nm <- suppressMessages(rsem_read(d, metadata = NULL, expr = "count", id = "annotation", annotation = file.path(d, "genes.gtf")))
  expect_equal(rownames(nm), c("ALPHA", "BETA", "GAMMA", "DELTA_X"))
  expect_equal(nm$sampA, c(900, 800, 50, 400))
})

test_that("rsem_model_qc reads the .model files RSEM 1.3.3 wrote", {
  d <- testthat::test_path("fixtures", "rsem-real")
  q <- suppressMessages(rsem_model_qc(d, sample_from = "parent_dir"))
  s <- q$sample_summary
  expect_equal(s$sample_id, c("sampA", "sampB", "sampC"))
  expect_equal(s$rspd_estimated, c(FALSE, FALSE, TRUE))
  expect_equal(s$rspd_bins, c(NA, NA, 20L))
  expect_equal(s$model_type, c(0, 0, 0)); expect_equal(s$forward_prob, c(0.5, 0.5, 0.5))
  expect_equal(c(s$length_lower[1], s$length_upper[1]), c(49, 250))
  a <- q$length_distribution[q$length_distribution$sample_id == "sampA", ]
  expect_equal(nrow(a), 201L); expect_equal(sum(a$probability), 1, tolerance = 1e-6)
  expect_equal(s$length_mean[1], sum(a$length * a$probability))
  r <- q$rspd
  expect_equal(unique(r$sample_id), "sampC"); expect_equal(nrow(r), 20L); expect_equal(sum(r$probability), 1, tolerance = 1e-6)
})

test_that("rsem_model_qc refuses files that are not .model files", {
  d <- withr::local_tempdir()
  f <- file.path(d, "x.model")
  writeLines(c("0", "", "0.5", "", "9 13 3", "0.1 0.4 0.4 0.1"), f)
  expect_error(suppressMessages(rsem_model_qc(d)), "span")
  writeLines(c("0", "", "0.5", "", "9 13 4", "0.1 0.4 0.4 0.5", "0", "0"), f)
  expect_error(suppressMessages(rsem_model_qc(d)), "add up to 1")
  writeLines(c("9", "0.5"), f)
  expect_error(suppressMessages(rsem_model_qc(d)), "model type")
  writeLines(c("0", "", "0.5", "", "9 13 4", "0.1 0.4 0.4 0.1", "0", "0", "1", "2", "0.25 0.75"), f)
  ex <- suppressMessages(rsem_model_qc(d))
  expect_equal(ex$sample_summary$length_mean, sum(10:13 * c(0.1, 0.4, 0.4, 0.1)))
  expect_equal(ex$sample_summary$length_mode, 11)
})
