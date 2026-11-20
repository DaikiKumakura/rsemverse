# Fixes and additions after 0.3.0: row annotation is carried, --append-names,
# gzip input, group/design in rsem_filter(), .cnt statistics.

write_rsem2 <- function(dir, name, ids, len, tpm, count = tpm) {
  tab <- data.frame(gene_id = ids, `transcript_id(s)` = "T", length = len + 100, effective_length = len,
                    expected_count = count, TPM = tpm, FPKM = tpm, check.names = FALSE)
  dir.create(dirname(file.path(dir, name)), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(tab, file.path(dir, name), sep = "\t", quote = FALSE, row.names = FALSE)
}

test_that("gene IDs and transcript tables survive the chained functions", {
  skip_if_not_installed("org.Hs.eg.db")
  x <- suppressMessages(read_example("tpm"))
  ids <- attr(x, "gene_ids")
  expect_false(is.null(ids))
  f <- suppressMessages(rsem_filter(x, "tpm", method = "detected", min_value = 5000, min_samples = 1))
  expect_lt(nrow(f), nrow(x))
  g <- attr(f, "gene_ids")
  expect_setequal(g$gene, rownames(f))
  expect_equal(g$gene_id, ids$gene_id[ids$gene %in% rownames(f)])
  results <- list(
    suppressMessages(rsem_normalize(x, "tpm", "quantile")),
    suppressMessages(rsem_transform(x, "log2")),
    suppressMessages(rsem_subset(x, "human", keep = "protein_coding")),
    suppressMessages(rsem_apply_reference(x, suppressMessages(rsem_fit_reference(x, c("S1", "S2")))))
  )
  for (res in results) {
    gi <- attr(res, "gene_ids")
    expect_false(is.null(gi))
    expect_setequal(gi$gene, rownames(res))
  }
  out <- tempfile(fileext = ".tsv")
  suppressMessages(rsem_write(f, out))
  back <- utils::read.delim(out, check.names = FALSE)
  expect_true("gene_id" %in% names(back))

  iso <- data.frame(S1 = c(5, 1, 9), S2 = c(4, 2, 8), row.names = c("T1", "T2", "T3"))
  attr(iso, "tx2gene") <- data.frame(transcript_id = rownames(iso), gene_id = c("G1", "G1", "G2"), gene = c("a", "a", "b"))
  k <- suppressMessages(rsem_filter(iso, "tpm", method = "detected", min_value = 3, min_samples = 1))
  expect_equal(attr(k, "tx2gene")$transcript_id, c("T1", "T3"))
})

test_that("--append-names IDs are read and any organism works with the original IDs", {
  expect_equal(
    rsemverse:::.strip_ensembl_version(c(
      "ENSG00000141510.18_TP53", "ENSMUSG00000000001.5_Gnai3", "ENSG00000141510.18",
      "ENSG00000182378.14_PAR_Y", "ERCC-00002", "10904"
    )),
    c("ENSG00000141510", "ENSMUSG00000000001", "ENSG00000141510", "ENSG00000182378.14_PAR_Y", "ERCC-00002", "10904")
  )
  e <- tempfile("plant"); dir.create(e)
  write_rsem2(e, "P1.genes.results", c("AT1G01010", "AT1G01020"), c(500, 600), c(1, 2))
  write_rsem2(e, "P2.genes.results", c("AT1G01010", "AT1G01020"), c(500, 600), c(3, 4))
  y <- suppressMessages(rsem_read(e, NULL, expr = "tpm", id = "original"))
  expect_equal(rownames(y), c("AT1G01010", "AT1G01020"))
  expect_error(rsem_read(e, NULL, expr = "tpm"), "species")
  skip_if_not_installed("org.Hs.eg.db")
  expect_error(rsem_read(e, NULL, species = "human", expr = "tpm"), "do not look like Ensembl")
  d <- tempfile("an"); dir.create(d)
  gid <- c("ENSG00000141510.18_TP53", "ENSG00000111640.18_GAPDH", "ERCC-00002")
  write_rsem2(d, "S1.genes.results", gid, c(1000, 900, 500), c(10, 20, 3))
  write_rsem2(d, "S2.genes.results", gid, c(1000, 900, 500), c(11, 19, 4))
  x <- suppressMessages(rsem_read(d, NULL, species = "human", expr = "tpm"))
  expect_equal(x["TP53", "S1"], 10)
  expect_equal(x["GAPDH", "S2"], 19)
  expect_true("ERCC-00002" %in% rownames(x))
  expect_equal(attr(x, "gene_ids")$gene_id, gid)
})

test_that("gzip-compressed RSEM files are found and read, also without a sample table", {
  d <- tempfile("gz"); dir.create(d)
  ids <- c("ENSG00000141510.5", "ENSG00000111640.18")
  write_rsem2(d, "A.genes.results", ids, c(10, 20), c(1, 2))
  write_rsem2(d, "B.genes.results", ids, c(10, 20), c(3, 4))
  for (f in c("A", "B")) {
    p <- file.path(d, paste0(f, ".genes.results"))
    con <- gzfile(paste0(p, ".gz"), "w"); writeLines(readLines(p), con); close(con); unlink(p)
  }
  x <- suppressMessages(rsem_read(d, NULL, expr = "tpm", id = "original"))
  expect_equal(colnames(x), c("A", "B"))
  expect_equal(x$B, c(3, 4))
  writeLines(c("sample_id,ref", "A,A.genes.results.gz", "B,B.genes.results.gz"), file.path(d, "m.csv"))
  y <- suppressMessages(rsem_read(d, file.path(d, "m.csv"), expr = "count", id = "original"))
  expect_equal(y$A, c(1, 2))
  expect_equal(rsemverse:::.rsem_header(file.path(d, "A.genes.results.gz"))[1], "gene_id")
})

test_that("isoform files from different references are refused", {
  d <- tempfile("isoref"); dir.create(d)
  mk <- function(s, genes) {
    tab <- data.frame(transcript_id = c("T1", "T2"), gene_id = genes, length = 100, effective_length = 90,
                      expected_count = 5, TPM = c(1, 2), FPKM = 1, IsoPct = 50)
    utils::write.table(tab, file.path(d, paste0(s, ".isoforms.results")), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  mk("S1", c("G1", "G1")); mk("S2", c("G1", "G2"))
  expect_error(rsem_read(d, NULL, expr = "tpm", level = "isoform", id = "original"), "different gene")
})

test_that("rsem_filter group and design follow edgeR filterByExpr", {
  skip_if_not_installed("edgeR")
  set.seed(5)
  m <- matrix(stats::rnbinom(300 * 13, mu = 8, size = 1), 300, 13,
              dimnames = list(paste0("g", 1:300), paste0("s", 1:13)))
  m[1:40, 1:3] <- stats::rnbinom(120, mu = 60, size = 2)
  x <- as.data.frame(m)
  grp <- rep(c("ctrl", "trt"), c(3, 10))
  kept <- function(...) rownames(suppressMessages(rsem_filter(x, "count", ...)))
  edge <- function(...) rownames(m)[suppressMessages(edgeR::filterByExpr(m, ...))]
  expect_equal(kept(), edge())
  expect_equal(kept(group = grp), edge(group = grp))
  expect_gt(length(kept(group = grp)), length(kept()))
  expect_equal(kept(group = stats::setNames(grp, colnames(m))), edge(group = grp))
  dm <- stats::model.matrix(~ factor(grp))
  expect_equal(kept(design = dm), edge(design = dm))
  g3 <- rep(c("a", "b", "c"), c(2, 5, 6))
  expect_equal(kept(group = g3), edge(group = g3))
  expect_error(rsem_filter(x, "count", group = grp[-1]), "one value per sample")
  expect_error(rsem_filter(x, "count", group = grp, design = dm), "not both")
  expect_error(rsem_filter(x, "count", method = "nonzero", group = grp), "edger")
  expect_message(rsem_filter(x, "count", group = grp), "smallest of 2 groups, 3 samples")
  big <- matrix(stats::rnbinom(200 * 26, mu = 10, size = 1), 200, 26,
                dimnames = list(paste0("g", 1:200), paste0("s", 1:26)))
  g26 <- rep(c("a", "b"), 13)
  expect_equal(
    rownames(suppressMessages(rsem_filter(as.data.frame(big), "count", group = g26))),
    rownames(big)[suppressMessages(edgeR::filterByExpr(big, group = g26))]
  )
})

test_that("rsem_biotype_profile uses given totals and gives NA for empty samples", {
  ann <- data.frame(gene_id = c("ENSG01", "ENSG02"), gene_name = c("RPL3", "RNA5-8SN1"),
                    gene_type = c("protein_coding", "rRNA"))
  x <- data.frame(S1 = c(10, 5), S2 = c(0, 0), row.names = c("RPL3", "RNA5-8SN1"))
  p <- suppressMessages(rsem_biotype_profile(x, ann, types = "rRNA"))
  expect_equal(p$share, c(1 / 3, NA))
  expect_message(rsem_biotype_profile(x, ann, types = "rRNA"), "total 0, share NA: S2")
  p2 <- suppressMessages(rsem_biotype_profile(x[2, ], ann, types = "rRNA", total = c(S1 = 100, S2 = 50)))
  expect_equal(p2$share, c(0.05, 0))
  expect_message(rsem_biotype_profile(x[2, ], ann, types = "rRNA", total = c(100, 50)), "totals given in `total`")
  expect_error(rsem_biotype_profile(x, ann, total = 1), "one number per sample")
})

cnt_file <- function(dir, name, n0 = 100, n1 = 850, n2 = 50, u = 800, mu = 50, unc = 120, hits = 1500, type = 2) {
  dir.create(dirname(file.path(dir, name)), recursive = TRUE, showWarnings = FALSE)
  writeLines(c(sprintf("%d %d %d %d   # counts", n0, n1, n2, n0 + n1 + n2), sprintf("%d %d %d", u, mu, unc),
               sprintf("%d %d", hits, type), "0\t100", "1\t700", "2\t100", "Inf\t50"), file.path(dir, name))
}

test_that("rsem_mapping_qc reads .cnt files and checks their identities", {
  d <- tempfile("cnt"); dir.create(d)
  cnt_file(d, "B.stat/B.cnt", n0 = 200, n1 = 700, n2 = 100, u = 600, mu = 100, unc = 50, hits = 900)
  cnt_file(d, "A.stat/A.cnt")
  q <- suppressMessages(rsem_mapping_qc(d))
  expect_equal(q$sample_id, c("A", "B"))
  a <- q[q$sample_id == "A", ]
  expect_equal(a$n_total, 1000)
  expect_equal(a$alignable_fraction, 0.85)
  expect_equal(a$unalignable_fraction, 0.1)
  expect_equal(a$filtered_fraction, 0.05)
  expect_equal(a$unique_gene_fraction, 800 / 850)
  expect_equal(a$multi_gene_fraction, 50 / 850)
  expect_equal(a$uncertain_fraction, 120 / 850)
  expect_equal(a$alignments_per_read, 1500 / 850)
  expect_equal(a$read_type, 2)
  expect_equal(q$alignable_fraction[q$sample_id == "B"], 0.7)
  p <- suppressMessages(rsem_mapping_qc(d, sample_from = "parent_dir"))
  expect_equal(p$sample_id, c("A", "B"))
  writeLines(c("sample_id,ref", "X,A.stat/A.cnt"), file.path(d, "m.csv"))
  expect_equal(suppressMessages(rsem_mapping_qc(d, file.path(d, "m.csv")))$sample_id, "X")
  expect_message(rsem_mapping_qc(d), "alignable fraction 0.7 to 0.85")
  bad <- tempfile("bad"); dir.create(bad)
  writeLines(c("100 850 50 999", "800 50 120", "1500 2"), file.path(bad, "S.cnt"))
  expect_error(rsem_mapping_qc(bad), "N_tot")
  writeLines(c("100 850 50 1000", "700 50 120", "1500 2"), file.path(bad, "S.cnt"))
  expect_error(rsem_mapping_qc(bad), "nUnique \\+ nMulti")
  writeLines(c("not a cnt file"), file.path(bad, "S.cnt"))
  expect_error(rsem_mapping_qc(bad), "not an RSEM .cnt file")
  expect_error(rsem_mapping_qc(tempfile("none")), "existing directory")
})
