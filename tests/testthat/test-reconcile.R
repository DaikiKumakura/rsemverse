write_pair <- function(dir, sid, genes, iso) {
  g <- data.frame(
    gene_id = genes$id, "transcript_id(s)" = genes$tx, length = 1000, effective_length = 900,
    expected_count = genes$count, TPM = genes$tpm, FPKM = genes$fpkm, check.names = FALSE
  )
  i <- data.frame(
    transcript_id = iso$tx, gene_id = iso$gene, length = 1000, effective_length = 900,
    expected_count = iso$count, TPM = iso$tpm, FPKM = iso$fpkm, IsoPct = 0
  )
  utils::write.table(g, file.path(dir, paste0(sid, ".genes.results")), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(i, file.path(dir, paste0(sid, ".isoforms.results")), sep = "\t", quote = FALSE, row.names = FALSE)
}

base_genes <- function() list(
  id = c("G1", "G2"), tx = c("T1,T2", "T3"), count = c(10.01, 5), tpm = c(600000, 400000), fpkm = c(1.5, 1)
)
base_iso <- function() list(
  tx = c("T1", "T2", "T3"), gene = c("G1", "G1", "G2"), count = c(4.00, 6.01, 5),
  tpm = c(250000.00, 350000.00, 400000), fpkm = c(0.5, 1, 1)
)

test_that("rsem_reconcile accepts differences explained by two-decimal rounding", {
  d <- withr::local_tempdir()
  iso <- base_iso(); iso$count <- c(4.004, 6.006, 5)  # sums to 10.010; gene value 10.01
  write_pair(d, "s1", base_genes(), iso)
  res <- suppressMessages(rsem_reconcile(d))
  expect_equal(unique(res$column), c("expected_count", "TPM", "FPKM"))
  expect_true(all(res$n_beyond_rounding == 0))
  expect_true(all(res$genes_without_isoforms == 0))
})

test_that("rsem_reconcile finds a gene whose transcripts do not add up", {
  d <- withr::local_tempdir()
  iso <- base_iso(); iso$tpm <- c(250000, 350000, 300000)
  write_pair(d, "s1", base_genes(), iso)
  res <- suppressMessages(rsem_reconcile(d, columns = "TPM"))
  expect_equal(res$n_beyond_rounding, 1L)
  expect_equal(res$worst_gene, "G2")
  expect_equal(res$max_abs_diff, 1e5)
})

test_that("rsem_reconcile reports genes without isoforms and stops on different samples", {
  d <- withr::local_tempdir()
  gg <- base_genes(); gg$id <- c("G1", "G9")
  write_pair(d, "s1", gg, base_iso())
  res <- suppressMessages(rsem_reconcile(d, columns = "TPM"))
  expect_equal(res$genes_without_isoforms, 1L)
  expect_equal(res$isoform_genes_missing, 1L)
  file.copy(file.path(d, "s1.genes.results"), file.path(d, "s2.genes.results"))
  expect_error(suppressMessages(rsem_reconcile(d)), "same samples")
  expect_error(suppressMessages(rsem_reconcile(d, columns = "IsoPct")), "columns")
})

test_that("rsem_id_report lists rows that combine several IDs", {
  x <- data.frame(s1 = c(3, 1), row.names = c("A", "B"))
  attr(x, "gene_ids") <- data.frame(gene = c("A", "A", "B"), gene_id = c("g1", "g2", "g3"))
  res <- suppressMessages(rsem_id_report(x))
  expect_equal(res$gene, "A")
  expect_equal(res$n_ids, 2L)
  expect_equal(res$gene_ids, "g1,g2")
  attr(x, "gene_ids") <- NULL
  expect_error(rsem_id_report(x), "gene_ids")
})

test_that("rsem_check reports the rounding allowance when the TPM total is off", {
  d <- withr::local_tempdir()
  g <- base_genes(); g$tpm <- c(600000, 300000)
  write_pair(d, "s1", g, base_iso())
  res <- suppressMessages(rsem_check(d))
  w <- res[res$check == "tpm_total", ]
  expect_equal(nrow(w), 1L)
  expect_match(w$detail, "rounding of 2 rows explains at most 0.01")
})
