golden <- function() readRDS(test_path("fixtures", "golden.rds"))

# ---- input handling -------------------------------------------------------

test_that("the table can have genes as row names, as a first column, or be a matrix", {
  m <- matrix(c(0, 5, 20, 1, 8, 30), 3, 2, dimnames = list(c("A", "B", "C"), c("S1", "S2")))
  by_row <- suppressMessages(rsem_filter(as.data.frame(m), "count", method = "nonzero"))
  first_col <- data.frame(gene = rownames(m), S1 = m[, 1], S2 = m[, 2])
  by_col <- suppressMessages(rsem_filter(first_col, "count", method = "nonzero"))
  by_mat <- suppressMessages(rsem_filter(m, "count", method = "nonzero"))
  expect_identical(by_row, by_col)
  expect_identical(by_row, by_mat)
  expect_equal(rownames(by_row), c("A", "B", "C"))
})

test_that("bad tables are refused with a reason", {
  m <- matrix(1:4, 2, 2, dimnames = list(c("A", "B"), c("S1", "S2")))
  expect_error(rsem_filter(list(1), "count"), "data frame or a numeric matrix")
  expect_no_error(suppressMessages(rsem_filter(m, "count", method = "nonzero")))
  bad <- m; bad[1, 1] <- NA
  expect_error(rsem_filter(bad, "count"), "missing values")
  bad <- m; bad[1, 1] <- -1
  expect_error(rsem_filter(bad, "count"), "negative values")
  bad <- m; rownames(bad) <- c("A", "A")
  expect_error(rsem_filter(bad, "count"), "duplicated gene names: A")
  bad <- m; colnames(bad) <- NULL
  expect_error(rsem_filter(bad, "count"), "sample names")
  expect_error(rsem_filter(data.frame(gene = "A", S1 = "x"), "count"), "numeric sample columns; not numeric: S1")
  expect_error(rsem_filter(m), "`expr` must be")
  expect_error(rsem_filter(m, "fpkm"), "`expr` must be one of")
})

# ---- rsem_filter ------------------------------------------------------------

test_that("\"edger\" keeps the same genes as edgeR::filterByExpr()", {
  g <- golden()
  for (set in c("six", "twelve")) {
    x <- as.data.frame(g[[set]]$counts)
    expect_message(kept <- rsem_filter(x, "count"), "rsem_filter \\(edger\\)")
    expect_equal(rownames(kept), rownames(x)[g[[set]]$filter_edger], label = set)
    kept5 <- suppressMessages(rsem_filter(x, "count", min_value = 5))
    expect_equal(rownames(kept5), rownames(x)[g[[set]]$filter_edger_min5], label = set)
  }
})

test_that("\"edger\" reports its rule and refuses TPM", {
  x <- as.data.frame(golden()$twelve$counts)
  msg <- testthat::capture_messages(rsem_filter(x, "count"))
  expect_match(msg, "in >= 11.4 of 12 samples and total count >= 15", fixed = TRUE)
  expect_error(rsem_filter(x, "tpm", method = "edger"), "needs counts")
})

test_that("\"detected\" counts the samples above a threshold", {
  m <- matrix(
    c(0, 0, 0, 0, 0,
      2, 0, 0, 0, 0,
      2, 2, 0, 0, 0,
      2, 2, 2, 2, 2,
      0.5, 0.5, 0.5, 0.5, 0.5),
    5, 5, byrow = TRUE, dimnames = list(paste0("g", 1:5), paste0("S", 1:5))
  )
  # Default for TPM: TPM >= 1 in at least 20% of samples (1 of 5).
  expect_equal(rownames(suppressMessages(rsem_filter(m, "tpm"))), c("g2", "g3", "g4"))
  expect_equal(rownames(suppressMessages(rsem_filter(m, "tpm", min_samples = 2))), c("g3", "g4"))
  expect_equal(rownames(suppressMessages(rsem_filter(m, "tpm", min_samples = 0.5))), "g4")
  expect_equal(rownames(suppressMessages(rsem_filter(m, "tpm", min_value = 0.5, min_samples = 5))), c("g4", "g5"))
  # Counts default to a threshold of 10.
  expect_equal(nrow(suppressMessages(rsem_filter(m, "count", method = "detected"))), 0L)
  expect_error(rsem_filter(m, "tpm", min_samples = 6), "larger than the number of samples")
  expect_match(testthat::capture_messages(rsem_filter(m, "tpm")), "TPM >= 1 in >= 1 of 5 samples", fixed = TRUE)
})

test_that("\"nonzero\" removes only genes that are zero everywhere", {
  m <- matrix(c(0, 0, 0, 1, 0, 0), 3, 2, dimnames = list(c("a", "b", "c"), c("S1", "S2")))
  expect_equal(rownames(suppressMessages(rsem_filter(m, "count", method = "nonzero"))), "a")
})

test_that("\"variable\" keeps the most variable genes in the original order", {
  m <- rbind(
    flat = c(10, 10, 10, 10),
    small = c(10, 12, 9, 11),
    big = c(1, 100, 1, 100),
    mid = c(5, 20, 5, 20),
    zero = c(0, 0, 0, 0)
  )
  colnames(m) <- paste0("S", 1:4)
  expect_equal(rownames(suppressMessages(rsem_filter(m, "tpm", method = "variable", n_top = 2))), c("big", "mid"))
  # More than the available genes: all genes except those that are zero everywhere.
  expect_equal(rownames(suppressMessages(rsem_filter(m, "tpm", method = "variable", n_top = 100))), c("flat", "small", "big", "mid"))
  expect_error(rsem_filter(m[, 1, drop = FALSE], "tpm", method = "variable"), "at least two samples")
  expect_error(rsem_filter(m, "tpm", method = "variable", n_top = 1.5), "whole number")
})

test_that("filtering does not change values or sample order", {
  x <- as.data.frame(golden()$six$counts)
  kept <- suppressMessages(rsem_filter(x, "count"))
  expect_equal(colnames(kept), colnames(x))
  expect_identical(as.matrix(kept), as.matrix(x[rownames(kept), ]))
})

# ---- rsem_normalize ---------------------------------------------------------

test_that("each method matches edgeR, DESeq2 and limma", {
  g <- golden()
  for (set in c("six", "twelve")) {
    ref <- g[[set]]
    x <- as.data.frame(ref$counts)
    for (method in c("cpm", "tmm", "upper_quartile", "deseq2", "quantile")) {
      out <- suppressMessages(rsem_normalize(x, "count", method))
      expect_equal(as.matrix(out), ref[[method]], tolerance = 1e-10, label = paste(set, method))
      expect_equal(dimnames(out), dimnames(x))
    }
    tmm <- suppressMessages(rsem_normalize(x, "count", "tmm"))
    expect_equal(unname(attr(tmm, "factors")), unname(colSums(ref$counts) * ref$tmm_factors), tolerance = 1e-10)
    msg <- testthat::capture_messages(dds <- rsem_normalize(x, "count", "deseq2"))
    expect_match(msg, "size factors: S01 [0-9.]+, S02")
    expect_equal(names(attr(dds, "factors")), colnames(x))
    expect_equal(unname(attr(dds, "factors")), unname(ref$deseq2_size_factors), tolerance = 1e-10)
  }
})

test_that("CPM columns sum to one million and quantile columns share one distribution", {
  x <- as.data.frame(golden()$six$counts)
  cpm <- suppressMessages(rsem_normalize(x, "count", "cpm"))
  expect_equal(unname(colSums(cpm)), rep(1e6, ncol(x)))
  m <- matrix(c(5, 2, 3, 4, 1, 6, 9, 8, 7), 3, 3, dimnames = list(c("a", "b", "c"), c("S1", "S2", "S3")))
  q <- as.matrix(suppressMessages(rsem_normalize(m, "tpm", "quantile")))
  expect_equal(apply(q, 2, sort), matrix(rep(rowMeans(apply(m, 2, sort)), 3), 3), ignore_attr = TRUE)
})

test_that("methods that need counts are refused for TPM, with the reason", {
  m <- matrix(1:6, 3, 2, dimnames = list(c("a", "b", "c"), c("S1", "S2")))
  for (method in c("cpm", "tmm", "upper_quartile", "deseq2")) {
    expect_error(rsem_normalize(m, "tpm", method), "needs counts")
  }
  expect_error(rsem_normalize(m, "count"), "`method` is required")
  expect_error(rsem_normalize(m, "count", "vst"), "`method` must be one of")
})

test_that("degenerate inputs give a clear error", {
  m <- matrix(c(0, 3, 4, 0), 2, 2, dimnames = list(c("a", "b"), c("S1", "S2")))
  expect_error(rsem_normalize(m, "count", "deseq2"), "nonzero in every sample")
  z <- cbind(m, S3 = 0)
  expect_error(rsem_normalize(z, "count", "cpm"), "nonzero total count; zero: S3")
  expect_error(rsem_filter(z, "count"), "nonzero total count; zero: S3")
})

test_that("a single sample is left unchanged by the scaling methods", {
  m <- matrix(c(10, 20, 70), 3, 1, dimnames = list(c("a", "b", "c"), "S1"))
  expect_equal(unname(suppressMessages(rsem_normalize(m, "count", "tmm"))[, 1]), c(1e5, 2e5, 7e5))
  expect_equal(unname(suppressMessages(rsem_normalize(m, "count", "quantile"))[, 1]), c(10, 20, 70))
})
