# Functions added in 0.3.0.

toy <- function() {
  set.seed(3)
  m <- matrix(round(stats::rexp(60 * 6, 1 / 20), 2), 60, 6,
              dimnames = list(sprintf("G%02d", 1:60), sprintf("S%d", 1:6)))
  m[1:5, ] <- 0
  m[6, 1:3] <- 7
  as.data.frame(m)
}

# ---- detection curve ----
test_that("rsem_detection_curve counts match a hand computation and rsem_filter", {
  x <- toy()
  m <- as.matrix(x)
  cur <- suppressMessages(rsem_detection_curve(x, "tpm", thresholds = c(0, 1, 10, 30), min_samples = c(1, 0.5, 6)))
  expect_equal(nrow(cur), 12L)
  for (i in seq_len(nrow(cur))) {
    k <- cur$samples_needed[i]
    expect_equal(cur$genes[i], sum(rowSums(m >= cur$threshold[i]) >= k))
    kept <- suppressMessages(rsem_filter(x, "tpm", method = "detected", min_value = cur$threshold[i],
                                         min_samples = cur$min_samples[i]))
    expect_equal(cur$genes[i], nrow(kept))
  }
  expect_equal(cur$samples_needed[cur$min_samples == 0.5][1], 3)
  expect_equal(cur$genes_mean[cur$threshold == 10][1], sum(rowMeans(m) >= 10))
  # more demanding thresholds never keep more genes
  for (k in unique(cur$min_samples)) expect_false(is.unsorted(rev(cur$genes[cur$min_samples == k])))
  expect_equal(cur$share, cur$genes / 60)
  expect_error(rsem_detection_curve(x, "tpm", min_samples = 7), "larger than the number of samples")
  expect_error(rsem_detection_curve(x, "tpm", thresholds = -1), "not negative")
})

# ---- reference z-scores ----
test_that("reference z-scores use only the reference samples", {
  x <- toy()
  ref <- suppressMessages(rsem_fit_reference(x, c("S1", "S2", "S3")))
  l <- log2(as.matrix(x) + 1)
  mu <- rowMeans(l[, 1:3])
  s <- apply(l[, 1:3], 1, stats::sd)
  z <- suppressMessages(rsem_apply_reference(x, ref))
  ok <- s > 0
  expect_equal(as.matrix(z)[ok, ], ((l - mu) / s)[ok, ], tolerance = 1e-12)
  expect_true(all(is.na(as.matrix(z)[!ok, ])))
  expect_equal(sum(!ok), 6L)  # five all-zero genes and G06 (7 in S1-S3)
  # adding or removing other samples does not change a sample's value
  z2 <- suppressMessages(rsem_apply_reference(x[, c("S4", "S5")], ref))
  expect_equal(as.matrix(z2), as.matrix(z)[, c("S4", "S5")])
  # drop
  zd <- suppressMessages(rsem_apply_reference(x, ref, zero_sd = "drop"))
  expect_equal(nrow(zd), sum(ok))
  # wrapper and transform = "none"
  expect_equal(suppressMessages(rsem_reference_scale(x, c("S1", "S2", "S3"))), z, ignore_attr = TRUE)
  r0 <- suppressMessages(rsem_fit_reference(x, c("S1", "S2"), transform = "none"))
  expect_equal(unname(r0$mean), unname(rowMeans(as.matrix(x)[, 1:2])))
})

test_that("reference groups come from the sample table and are checked", {
  x <- toy()
  meta <- data.frame(sample_id = sprintf("S%d", 1:6), condition = rep(c("ctrl", "drug"), each = 3))
  a <- suppressMessages(rsem_fit_reference(x, "ctrl", metadata = meta, group = "condition"))
  b <- suppressMessages(rsem_fit_reference(x, c("S1", "S2", "S3")))
  expect_equal(a, b)
  expect_error(rsem_fit_reference(x, "S1"), "at least 2")
  expect_error(rsem_fit_reference(x, "none", metadata = meta, group = "condition"), "not found")
  expect_error(rsem_fit_reference(x, "ctrl", metadata = meta), "together")
  expect_error(rsem_fit_reference(x, "S9"), "not in `x`")
  y <- x[-(1:10), ]
  z <- suppressMessages(rsem_apply_reference(y, a))
  expect_equal(nrow(z), 50L)
  expect_message(rsem_apply_reference(y, a), "reference genes missing from x: 10")
  expect_error(rsem_apply_reference(x, list()), "rsem_fit_reference")
  expect_output(print(a), "3 reference samples")
})

# ---- TMMwsp ----
test_that("TMMwsp equals edgeR", {
  skip_if_not_installed("edgeR")
  set.seed(11)
  m <- matrix(stats::rnbinom(400 * 5, mu = 30, size = 0.5), 400, 5,
              dimnames = list(paste0("g", 1:400), paste0("s", 1:5)))
  m[sample(length(m), 600)] <- 0
  x <- as.data.frame(m)
  out <- suppressMessages(rsem_normalize(x, "count", "tmmwsp"))
  y <- edgeR::DGEList(m)
  y <- edgeR::calcNormFactors(y, method = "TMMwsp")
  expect_equal(unname(attr(out, "factors")), unname(y$samples$lib.size * y$samples$norm.factors), tolerance = 1e-12)
  expect_equal(as.matrix(out), edgeR::cpm(y), tolerance = 1e-12, ignore_attr = TRUE)
  expect_error(rsem_normalize(x, "tpm", "tmmwsp"), "needs counts")
})

# ---- GTF annotation ----
gtf_file <- function() {
  f <- tempfile(fileext = ".gtf")
  row <- function(id, name, type, feature = "gene") {
    paste("chr1", "SRC", feature, 1, 100, ".", "+", ".",
          sprintf('gene_id "%s"; gene_type "%s"; gene_name "%s";', id, type, name), sep = "\t")
  }
  writeLines(c(
    "##description: test",
    row("ENSG00000000001.3", "RPL3", "protein_coding"),
    row("ENSG00000000002.1", "RNA5-8SN1", "rRNA"),
    row("ENSG00000000002.1", "RNA5-8SN1", "rRNA", feature = "exon"),
    row("ENSG00000000003.1", "MT-RNR1", "Mt_rRNA"),
    row("ENSG00000000004.1", "SNORD3A", "snoRNA"),
    row("ENSG00000000005.2", "MALAT1", "lncRNA")
  ), f)
  f
}

test_that("the GTF reader takes gene rows and strips versions", {
  a <- rsemverse:::.read_annotation(gtf_file())
  expect_equal(nrow(a), 5L)
  expect_equal(a$gene_id[2], "ENSG00000000002")
  expect_equal(a$gene_type, c("protein_coding", "rRNA", "Mt_rRNA", "snoRNA", "lncRNA"))
  expect_equal(a$gene_name[5], "MALAT1")
  gz <- tempfile(fileext = ".gtf.gz")
  con <- gzfile(gz, "w"); writeLines(readLines(gtf_file()), con); close(con)
  expect_equal(rsemverse:::.read_annotation(gz), a)
})

test_that("rsem_subset removes rRNA, not ribosomal proteins, with a GTF", {
  x <- data.frame(S1 = 1:6, S2 = 2:7,
                  row.names = c("RPL3", "RNA5-8SN1", "ENSG00000000003", "SNORD3A", "MALAT1", "OTHER"))
  out <- suppressMessages(rsem_subset(x, remove = c("rrna", "mito_rrna"), annotation = gtf_file()))
  expect_equal(rownames(out), c("RPL3", "SNORD3A", "MALAT1", "OTHER"))
  kept <- suppressMessages(rsem_subset(x, keep = "lncrna", annotation = gtf_file()))
  expect_equal(rownames(kept), "MALAT1")
  expect_error(rsem_subset(x, remove = "mirna", annotation = gtf_file()), "cannot be evaluated")
  expect_error(rsem_subset(x, remove = "rrna"), "give `annotation`")
  expect_error(rsem_subset(x, remove = "mito"), "species")
})

test_that("rsem_biotype_profile gives shares and NA for types not in the annotation", {
  x <- data.frame(S1 = c(10, 5, 1, 4), S2 = c(20, 0, 2, 8), row.names = c("RPL3", "RNA5-8SN1", "MALAT1", "XYZ"))
  p <- suppressMessages(rsem_biotype_profile(x, gtf_file()))
  s1 <- p[p$sample_id == "S1", ]
  expect_equal(sum(s1$share), 1)
  expect_equal(s1$share[s1$type == "rRNA"], 0.25)
  expect_equal(s1$share[s1$type == "unannotated"], 0.2)
  p2 <- suppressMessages(rsem_biotype_profile(x, gtf_file(), types = c("rRNA", "miRNA")))
  expect_true(all(is.na(p2$share[p2$type == "miRNA"])))
  expect_message(rsem_biotype_profile(x, gtf_file(), types = "miRNA"), "not measurable")
})

# ---- rsem_write ----
test_that("rsem_write writes gene IDs and annotation, also compressed", {
  x <- data.frame(S1 = c(1.5, 2), S2 = c(3, 4), row.names = c("RPL3", "MALAT1"))
  ids <- data.frame(gene = c("RPL3", "MALAT1", "MALAT1"),
                    gene_id = c("ENSG00000000001.3", "ENSG00000000005.2", "ENSG00000000009.1"))
  f <- tempfile(fileext = ".tsv.gz")
  suppressMessages(rsem_write(x, f, ids = ids, annotation = gtf_file()))
  back <- utils::read.delim(gzfile(f), check.names = FALSE)
  expect_equal(names(back), c("gene", "gene_id", "gene_type", "S1", "S2"))
  expect_equal(back$gene_id[2], "ENSG00000000005.2;ENSG00000000009.1")
  expect_equal(back$gene_type, c("protein_coding", "lncRNA"))
  expect_equal(back$S2, c(3, 4))
  expect_error(rsem_write(x, tempfile(), ids = ids[1, ]), "does not cover")
})

# ---- isoforms ----
test_that("rsem_select_isoforms follows the Trinity rules and combines samples as asked", {
  x <- data.frame(S1 = c(10, 1, 0.5, 4, 0), S2 = c(8, 3, 0, 0, 6),
                  row.names = c("T1", "T2", "T3", "T4", "T5"))
  t2g <- data.frame(transcript_id = rownames(x), gene_id = c("G1", "G1", "G1", "G2", "G2"))
  sel <- function(...) rownames(suppressMessages(rsem_select_isoforms(x, t2g, ...)))
  expect_equal(sel(min_tpm = 1), c("T1", "T2", "T4", "T5"))
  expect_equal(sel(min_tpm = 1, rule = "all"), c("T1", "T2"))
  expect_equal(sel(min_tpm = 1, min_fraction = 0.2), c("T1", "T2", "T4", "T5"))  # T2: 3/8 in S2
  expect_equal(sel(min_tpm = 1, min_fraction = 0.2, rule = "all"), "T1")
  expect_equal(sel(highest_only = TRUE), c("T1", "T4", "T5"))       # T4 tops S1, T5 tops S2
  expect_equal(sel(highest_only = TRUE, rule = "mean"), c("T1", "T5"))  # means 2 and 3
  expect_equal(sel(min_tpm = 1, rule = "prevalence", prevalence = 0.5), c("T1", "T2", "T4", "T5"))
  expect_error(rsem_select_isoforms(x, t2g, min_fraction = 0.1, highest_only = TRUE), "cannot be combined")
  expect_error(rsem_select_isoforms(x, t2g[-1, ]), "not in `tx2gene`")
  out <- suppressMessages(rsem_select_isoforms(x, t2g, highest_only = TRUE, rule = "mean"))
  expect_equal(attr(out, "tx2gene")$gene_id, c("G1", "G2"))
  expect_message(rsem_select_isoforms(x, t2g), "not re-quantified")
})

# ---- length profile ----
test_that("rsem_length_profile summarizes lengths and finds zeros", {
  len <- data.frame(S1 = c(1000, 0, 500, 0), S2 = c(1100, 200, 500, 0), S3 = c(900, 400, 500, 0),
                    row.names = c("A", "B", "C", "D"))
  p <- suppressMessages(rsem_length_profile(len))
  expect_equal(p$zero_fraction, c(0, 1 / 3, 0, 1))
  expect_equal(p$max_min_ratio, c(1100 / 900, 2, 1, NA))
  expect_equal(p$cv[1], stats::sd(c(1000, 1100, 900)) / 1000)
  expect_true(is.na(p$cv[4]))
  expect_equal(p$median, c(1000, 200, 500, 0))
})

# ---- rsem_read: discovery, lengths, original IDs, isoforms ----
write_rsem <- function(dir, name, ids, len, tpm, count = tpm) {
  tab <- data.frame(gene_id = ids, `transcript_id(s)` = "T", length = len + 100, effective_length = len,
                    expected_count = count, TPM = tpm, FPKM = tpm, check.names = FALSE)
  dir.create(dirname(file.path(dir, name)), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(tab, file.path(dir, name), sep = "\t", quote = FALSE, row.names = FALSE)
}
ids <- c("ENSG00000141510.5", "ENSG00000111640.18", "ENSG00000000003.15", "ERCC-00002")

test_that("rsem_read finds files by file name or folder and refuses duplicates", {
  d <- tempfile("disc"); dir.create(d)
  write_rsem(d, "B.genes.results", ids, c(1, 2, 3, 4), c(1, 2, 3, 4))
  write_rsem(d, "sub/A.genes.results", ids, c(1, 2, 3, 4), c(5, 6, 7, 8))
  x <- suppressMessages(rsem_read(d, NULL, species = "human", expr = "tpm", id = "original"))
  expect_equal(colnames(x), c("A", "B"))
  expect_equal(x$A, c(5, 6, 7, 8))
  expect_message(rsem_read(d, NULL, species = "human", expr = "tpm", id = "original"), "A <- sub/A.genes.results")
  expect_error(rsem_read(d, NULL, species = "human", expr = "tpm", sample_from = "parent_dir"), "own folder")
  write_rsem(d, "sub2/B.genes.results", ids, c(1, 2, 3, 4), c(1, 2, 3, 4))
  expect_error(rsem_read(d, NULL, species = "human", expr = "tpm"), "same sample ID")
  e <- tempfile("dirs"); dir.create(e)
  write_rsem(e, "P1/rsem.genes.results", ids, 1:4, 1:4)
  write_rsem(e, "P2/rsem.genes.results", ids, 1:4, 4:1)
  y <- suppressMessages(rsem_read(e, NULL, species = "human", expr = "tpm", id = "original", sample_from = "parent_dir"))
  expect_equal(colnames(y), c("P1", "P2"))
  empty <- tempfile("empty"); dir.create(empty)
  expect_error(rsem_read(empty, NULL, species = "human", expr = "tpm", level = "isoform"), "No")
})

test_that("rsem_read reads effective lengths, reports zeros and keeps original IDs", {
  d <- tempfile("len"); dir.create(d)
  write_rsem(d, "S1.genes.results", ids, c(1000, 0, 300, 50), c(10, 0, 1, 2))
  write_rsem(d, "S2.genes.results", ids, c(1200, 0, 0, 50), c(10, 0, 0, 2))
  msg <- testthat::capture_messages(x <- rsem_read(d, NULL, species = "human", expr = "length", id = "original"))
  expect_true(any(grepl("2 in at least one sample, 1 in every sample", msg)))
  expect_equal(attr(x, "zero_length"), c(S1 = 1, S2 = 2))
  expect_equal(rownames(x), ids)
  skip_if_not_installed("org.Hs.eg.db")
  s <- suppressMessages(rsem_read(d, NULL, species = "human", expr = "length"))
  expect_equal(s["TP53", "S2"], 1200)
  g <- attr(s, "gene_ids")
  expect_equal(g$gene_id, ids)
  expect_equal(g$gene[1], "TP53")
})

test_that("effective lengths of rows merged into one symbol are TPM-weighted", {
  len <- matrix(c(100, 300, 50, 70), 2, 2, dimnames = list(c("a", "b"), c("S1", "S2")))
  tpm <- matrix(c(1, 3, 0, 0), 2, 2, dimnames = dimnames(len))
  out <- rsemverse:::.weighted_length_by_symbol(len, tpm, c("X", "X"))
  expect_equal(out["X", "S1"], (100 * 1 + 300 * 3) / 4)
  expect_equal(out["X", "S2"], 60)
})

test_that("rsem_read reads isoforms with a tx2gene table", {
  d <- tempfile("iso"); dir.create(d)
  for (s in c("S1", "S2")) {
    tab <- data.frame(transcript_id = c("ENST01.1", "ENST02.1", "ENST03.2"),
                      gene_id = c("ENSG00000141510.5", "ENSG00000141510.5", "ENSG00000111640.18"),
                      length = 1000, effective_length = 900, expected_count = 5,
                      TPM = if (s == "S1") c(3, 1, 6) else c(0, 4, 2), FPKM = 1, IsoPct = c(75, 25, 100))
    utils::write.table(tab, file.path(d, paste0(s, ".isoforms.results")), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  x <- suppressMessages(rsem_read(d, NULL, species = "human", expr = "tpm", level = "isoform", id = "original"))
  expect_equal(rownames(x), c("ENST01.1", "ENST02.1", "ENST03.2"))
  expect_equal(x$S2, c(0, 4, 2))
  t2g <- attr(x, "tx2gene")
  expect_equal(t2g$gene_id[3], "ENSG00000111640.18")
  p <- suppressMessages(rsem_read(d, NULL, species = "human", expr = "isopct", level = "isoform", id = "original"))
  expect_equal(p$S1, c(75, 25, 100))
  sel <- suppressMessages(rsem_select_isoforms(x, highest_only = TRUE, rule = "mean"))
  expect_equal(rownames(sel), c("ENST02.1", "ENST03.2"))
  expect_error(rsem_read(d, NULL, species = "human", expr = "isopct"), "must be one of")
})

test_that("GTF types are found through the original IDs when names differ", {
  x <- data.frame(S1 = c(4, 6), S2 = c(1, 1), row.names = c("NEWNAME", "RPL3"))
  ids <- data.frame(gene = c("NEWNAME", "RPL3"), gene_id = c("ENSG00000000002.1", "ENSG00000000001.3"))
  p0 <- suppressMessages(rsem_biotype_profile(x, gtf_file(), types = "rRNA", ids = NULL))
  expect_equal(p0$share, c(0, 0))
  attr(x, "gene_ids") <- ids
  p <- suppressMessages(rsem_biotype_profile(x, gtf_file(), types = "rRNA"))
  expect_equal(p$share, c(0.4, 0.5))
  out <- suppressMessages(rsem_subset(x, remove = "rrna", annotation = gtf_file()))
  expect_equal(rownames(out), "RPL3")
})
