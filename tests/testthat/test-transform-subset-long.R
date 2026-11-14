small <- function() {
  matrix(
    c(0, 1, 3,
      4, 4, 4,
      10, 0, 30),
    3, 3, byrow = TRUE,
    dimnames = list(c("g1", "g2", "g3"), c("S1", "S2", "S3"))
  )
}

# ---- rsem_transform ---------------------------------------------------------

test_that("\"log2\" adds the pseudocount and records the step", {
  m <- small()
  out <- suppressMessages(rsem_transform(m, "log2"))
  expect_equal(as.matrix(out), log2(m + 1))
  expect_equal(attr(out, "transform"), "log2")
  expect_equal(as.matrix(suppressMessages(rsem_transform(m, "log2", pseudocount = 0.5))), log2(m + 0.5))
  expect_error(rsem_transform(m, "log2", pseudocount = 0), "-Inf")
  expect_error(rsem_transform(out, "log2"), "already log2-transformed")
})

test_that("\"zscore\" gives mean 0 and standard deviation 1 per gene, and 0 for constant genes", {
  m <- small()
  z <- as.matrix(suppressMessages(rsem_transform(m, c("log2", "zscore"))))
  ref <- t(scale(t(log2(m + 1))))
  ref["g2", ] <- 0
  expect_equal(z, ref, ignore_attr = TRUE)
  expect_equal(unname(rowMeans(z)), c(0, 0, 0))
  expect_equal(unname(apply(z[c("g1", "g3"), ], 1, stats::sd)), c(1, 1))
  msg <- testthat::capture_messages(rsem_transform(m, c("log2", "zscore")))
  expect_match(msg, "log2(value + 1), then z-score per gene (1 constant genes set to 0)", fixed = TRUE)
  expect_error(rsem_transform(m[, 1, drop = FALSE], "zscore"), "at least two samples")
  expect_error(rsem_transform(m, c("zscore", "log2")), "cannot follow")
})

test_that("\"center\" subtracts the gene mean", {
  m <- small()
  expect_equal(as.matrix(suppressMessages(rsem_transform(m, "center"))), m - rowMeans(m))
})

test_that("\"rank\" ranks genes within each sample, ties sharing the average", {
  m <- small()
  r <- as.matrix(suppressMessages(rsem_transform(m, "rank")))
  expect_equal(r[, "S1"], c(g1 = 1, g2 = 2, g3 = 3) / 3)
  expect_equal(r[, "S2"], c(g1 = 2, g2 = 3, g3 = 1) / 3)
  m2 <- m; m2[, "S3"] <- c(5, 5, 9)
  expect_equal(as.matrix(suppressMessages(rsem_transform(m2, "rank")))[, "S3"], c(g1 = 1.5, g2 = 1.5, g3 = 3) / 3)
  # Scale-free: the same ranks after a log transform.
  expect_equal(r, as.matrix(suppressMessages(rsem_transform(m, c("log2", "rank")))), ignore_attr = TRUE)
})

test_that("unknown methods are refused", {
  expect_error(rsem_transform(small()), "`method` is required")
  expect_error(rsem_transform(small(), "vst"), "one or more of")
})

# ---- rsem_subset ------------------------------------------------------------

human_table <- function() {
  genes <- c("TP53", "GAPDH", "MT-CO1", "MT-ND1", "RPL13A", "RPS6", "RPS6KA1", "RPL10P1",
             "HBB", "HBA1", "HBEGF", "XIST", "DDX3Y", "CD99", "MALAT1", "ERCC-00002",
             "ENSG00000141510", "NOT_A_GENE")
  m <- matrix(seq_along(genes) * 1.0, length(genes), 2, dimnames = list(genes, c("S1", "S2")))
  as.data.frame(m)
}

subset_names <- function(...) rownames(suppressMessages(rsem_subset(...)))

test_that("each built-in group matches the expected human genes", {
  skip_if_not_installed("org.Hs.eg.db")
  x <- human_table()
  expect_setequal(subset_names(x, "human", keep = "mito"), c("MT-CO1", "MT-ND1"))
  expect_setequal(subset_names(x, "human", keep = "ribo"), c("RPL13A", "RPS6"))
  expect_setequal(subset_names(x, "human", keep = "hb"), c("HBB", "HBA1"))
  # CD99 is pseudoautosomal (on X and Y) and is not a sex gene.
  expect_setequal(subset_names(x, "human", keep = "sex"), c("XIST", "DDX3Y"))
  expect_setequal(subset_names(x, "human", keep = "ercc"), "ERCC-00002")
  pc <- subset_names(x, "human", keep = "protein_coding")
  expect_true(all(c("TP53", "GAPDH", "MT-CO1", "RPS6KA1", "HBEGF", "DDX3Y", "CD99", "ENSG00000141510") %in% pc))
  expect_false(any(c("MALAT1", "XIST", "RPL10P1", "ERCC-00002", "NOT_A_GENE") %in% pc))
})

test_that("keep is applied before remove, and the order is kept", {
  skip_if_not_installed("org.Hs.eg.db")
  x <- human_table()
  out <- suppressMessages(rsem_subset(x, "human", keep = "protein_coding", remove = c("mito", "ribo", "hb", "sex")))
  expect_equal(rownames(out), c("TP53", "GAPDH", "RPS6KA1", "HBEGF", "CD99", "ENSG00000141510"))
  expect_identical(as.matrix(out), as.matrix(x[rownames(out), ]))
  removed <- subset_names(x, "human", remove = c("mito", "ercc"))
  expect_equal(removed, setdiff(rownames(x), c("MT-CO1", "MT-ND1", "ERCC-00002")))
  msg <- testthat::capture_messages(rsem_subset(x, "human", keep = "protein_coding", remove = "mito"))
  expect_match(msg, "kept 11 of 18 genes; keep protein_coding; remove mito (genes per group: protein_coding 13, mito 2)", fixed = TRUE)
})

test_that("mouse groups use mouse symbols", {
  skip_if_not_installed("org.Mm.eg.db")
  genes <- c("Trp53", "mt-Co1", "Rpl13a", "Rps6ka1", "Hbb-bs", "Hba-a1", "Xist", "Ddx3y", "Malat1")
  x <- as.data.frame(matrix(1, length(genes), 2, dimnames = list(genes, c("A", "B"))))
  expect_setequal(subset_names(x, "mouse", keep = "mito"), "mt-Co1")
  expect_setequal(subset_names(x, "mouse", keep = "ribo"), "Rpl13a")
  expect_setequal(subset_names(x, "mouse", keep = "hb"), c("Hbb-bs", "Hba-a1"))
  expect_setequal(subset_names(x, "mouse", keep = "sex"), c("Xist", "Ddx3y"))
  expect_false("Malat1" %in% subset_names(x, "mouse", keep = "protein_coding"))
})

test_that("rsem_subset checks its arguments", {
  x <- human_table()
  expect_error(rsem_subset(x, keep = "mito"), "`species` must be")
  expect_error(rsem_subset(x, "human"), "at least one gene group")
  expect_error(rsem_subset(x, "human", keep = "lncRNA"), "`keep` must be one or more of")
  expect_error(rsem_subset(x, "human", remove = NA_character_), "`remove` must be one or more of")
})

test_that("groups that need no annotation work without the OrgDb lookup", {
  x <- human_table()
  expect_setequal(subset_names(x, "human", keep = "ercc"), "ERCC-00002")
})

# ---- rsem_long --------------------------------------------------------------

test_that("rsem_long returns one row per gene and sample", {
  m <- small()
  long <- suppressMessages(rsem_long(m))
  expect_equal(names(long), c("gene", "sample_id", "value"))
  expect_equal(nrow(long), 9L)
  expect_equal(levels(long$sample_id), c("S1", "S2", "S3"))
  expect_equal(long$gene[1:3], rep("g1", 3))
  for (i in seq_len(nrow(long))) {
    expect_equal(long$value[i], m[long$gene[i], as.character(long$sample_id[i])])
  }
})

test_that("every column of the sample table is added with its own type", {
  dir <- tempfile("rsemverse-long-")
  dir.create(dir)
  meta <- file.path(dir, "metadata.csv")
  writeLines(c(
    "sample_id,ref,group,age,dose,responder",
    "S3,S3.genes.results,treated,61,0.5,TRUE",
    "S1,S1.genes.results,control,45,0,FALSE",
    "S2,S2.genes.results,treated,52,1.5,TRUE",
    "S9,S9.genes.results,control,70,0,FALSE"
  ), meta)
  msg <- testthat::capture_messages(long <- rsem_long(small(), metadata = meta))
  expect_equal(names(long), c("gene", "sample_id", "value", "ref", "group", "age", "dose", "responder"))
  expect_type(long$group, "character")
  expect_true(is.numeric(long$age) && is.numeric(long$dose))
  expect_type(long$responder, "logical")
  s2 <- long[long$sample_id == "S2", ]
  expect_true(all(s2$group == "treated" & s2$age == 52 & s2$dose == 1.5))
  expect_match(msg, "1 sample(s) in the sample table but not in the data were left out", fixed = TRUE)
  expect_match(msg, "group (text), age (number)", fixed = TRUE)
})

test_that("sample IDs stay text and a data frame works as the sample table", {
  m <- small(); colnames(m) <- c("001", "002", "010")
  meta <- data.frame(sample_id = c("001", "002", "010"), batch = c("a", "a", "b"))
  long <- suppressMessages(rsem_long(m, metadata = meta))
  expect_equal(levels(long$sample_id), c("001", "002", "010"))
  expect_equal(unique(long$batch[long$sample_id == "010"]), "b")
})

test_that("missing samples and clashing columns are refused", {
  m <- small()
  expect_error(rsem_long(m, metadata = data.frame(sample_id = c("S1", "S2"))), "not in the sample table: S3")
  expect_error(rsem_long(m, metadata = data.frame(id = "S1")), "'sample_id' column")
  expect_error(rsem_long(m, metadata = data.frame(sample_id = c("S1", "S1", "S2", "S3"))), "duplicated sample_id")
  expect_error(rsem_long(m, metadata = data.frame(sample_id = c("S1", "S2", "S3"), value = 1)), "rename them")
  expect_error(rsem_long(m, metadata = "no-such.csv"), "existing CSV")
})

test_that("negative values from rsem_transform() are accepted", {
  z <- suppressMessages(rsem_transform(small(), c("log2", "zscore")))
  long <- suppressMessages(rsem_long(z))
  expect_true(any(long$value < 0))
})
