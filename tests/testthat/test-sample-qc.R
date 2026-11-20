qc_toy <- function() {
  set.seed(21)
  m <- matrix(round(stats::rexp(120 * 6, 1 / 15), 2), 120, 6,
              dimnames = list(sprintf("g%03d", 1:120), sprintf("S%d", 1:6)))
  m[1:10, ] <- 0
  m[11:14, 4:6] <- m[11:14, 4:6] * 8
  as.data.frame(m)
}

test_that("rsem_sample_qc agrees with prcomp and cor", {
  x <- qc_toy()
  m <- as.matrix(x)
  qc <- suppressMessages(rsem_sample_qc(x, n_top = 50, n_pcs = 4))
  l <- log2(m + 1)
  v <- apply(l, 1, stats::var)
  v[rowSums(m != 0) == 0] <- -Inf
  top <- sort(order(-v)[1:50])
  p <- stats::prcomp(t(l[top, ]), center = TRUE, scale. = FALSE)
  expect_equal(unname(as.matrix(qc$pca_scores[, c("PC1", "PC2", "PC3", "PC4")])), unname(p$x[, 1:4]), tolerance = 1e-12)
  expect_equal(unname(qc$pca_variance$variance), p$sdev[1:4]^2)
  expect_equal(qc$pca_variance$cumulative, cumsum(p$sdev^2)[1:4] / sum(p$sdev^2))
  expect_equal(unname(qc$pca_loadings), unname(p$rotation[, 1:4]), tolerance = 1e-12)
  expect_equal(rownames(qc$pca_loadings), rownames(m)[top])
  det <- rowSums(m >= 1) >= 1
  expect_equal(qc$correlation_matrix, stats::cor(l[det, ], method = "spearman"))
  expect_equal(qc$sample_distance, 1 - qc$correlation_matrix)
  expect_equal(qc$parameters$pca_genes, 50L)
  expect_equal(qc$parameters$correlation_gene_count, sum(det))
  s <- qc$sample_summary
  expect_equal(s$total, unname(colSums(m)))
  expect_equal(s$detected_genes, unname(colSums(m >= 1)))
  expect_equal(s$detected_share, unname(colSums(m >= 1)) / 120)
  off <- qc$correlation_matrix; diag(off) <- NA
  expect_equal(s$median_correlation, unname(apply(off, 2, stats::median, na.rm = TRUE)))
  expect_equal(s$most_similar, colnames(m)[apply(off, 2, which.max)])
  expect_s3_class(qc, "rsem_sample_qc")
  expect_output(print(qc), "rsem_sample_qc: 120 genes x 6 samples")
})

test_that("options change the genes and the correlation as documented", {
  x <- qc_toy()
  m <- as.matrix(x)
  a <- suppressMessages(rsem_sample_qc(x, n_top = 30, correlation_genes = "variable", correlation = "pearson"))
  expect_equal(a$parameters$correlation_gene_count, 30L)
  top <- rownames(a$pca_loadings)
  expect_equal(a$correlation_matrix, stats::cor(log2(m[top, ] + 1)), tolerance = 1e-12)
  b <- suppressMessages(rsem_sample_qc(x, correlation_genes = "all"))
  expect_equal(b$parameters$correlation_gene_count, 110L)
  s <- suppressMessages(rsem_sample_qc(x, n_top = 40, pca_scale = TRUE, n_pcs = 2))
  l <- log2(m + 1); v <- apply(l, 1, stats::var); v[rowSums(m != 0) == 0] <- -Inf
  p <- stats::prcomp(t(l[sort(order(-v)[1:40]), ]), center = TRUE, scale. = TRUE)
  expect_equal(unname(as.matrix(s$pca_scores[, c("PC1", "PC2")])), unname(p$x[, 1:2]), tolerance = 1e-12)
  n <- suppressMessages(rsem_sample_qc(log2(x + 1), transform = "none", n_top = 40))
  expect_equal(n$pca_variance$variance[1:2], suppressMessages(rsem_sample_qc(x, n_top = 40))$pca_variance$variance[1:2])
  expect_equal(nrow(suppressMessages(rsem_sample_qc(x, n_pcs = 50))$pca_variance), 6L)
})

test_that("the sample table is joined and checked; bad input stops", {
  x <- qc_toy()
  meta <- data.frame(sample_id = rev(colnames(x)), group = rev(rep(c("a", "b"), 3)), age = 1:6)
  qc <- suppressMessages(rsem_sample_qc(x, metadata = meta, n_top = 30))
  expect_equal(qc$pca_scores$sample_id, colnames(x))
  expect_equal(qc$sample_summary$group, rep(c("a", "b"), 3))
  expect_equal(qc$sample_summary$age, 6:1)
  expect_error(rsem_sample_qc(x, metadata = meta[-1, ]), "not in the sample table")
  expect_error(rsem_sample_qc(x[, 1:2]), "at least three samples")
  expect_error(rsem_sample_qc(log2(x + 1) - 3, n_top = 10), "negative")
  const <- as.data.frame(matrix(5, 20, 4, dimnames = list(paste0("g", 1:20), paste0("S", 1:4))))
  expect_error(rsem_sample_qc(const), "Fewer than two genes vary")
  expect_error(rsem_sample_qc(x, correlation = "none"), "must be one of")
  expect_error(rsem_sample_qc(x, n_top = 1), "at least 2")
})

test_that("a sample with all zeros gives NA correlations instead of an error", {
  x <- qc_toy()
  x$S6 <- 0
  qc <- suppressMessages(rsem_sample_qc(x, n_top = 30))
  expect_true(all(is.na(qc$correlation_matrix[, "S6"])))
  expect_true(is.na(qc$sample_summary$median_correlation[6]))
  expect_equal(qc$sample_summary$detected_genes[6], 0)
})
