test_that("arguments are required and checked", {
  d <- file.path(example_dir(), "data")
  m <- file.path(example_dir(), "metadata.csv")
  expect_error(rsem_read(d, m, expr = "count"), "`species` must be")
  expect_error(rsem_read(d, m, species = "human"), "`expr` must be")
  expect_error(rsem_read(d, m, species = "rat", expr = "count"), "`species` must be one of")
  expect_error(rsem_read(d, m, species = "human", expr = "fpkm"), "`expr` must be one of")
  expect_error(rsem_read(d, m, species = c("human", "mouse"), expr = "count"), "`species` must be one of")
  expect_error(rsem_read("no-such-dir", m, species = "human", expr = "count"), "existing directory")
  expect_error(rsem_read(d, "no-such.csv", species = "human", expr = "count"), "existing CSV")
})

test_that("the example returns a symbol-by-sample data frame", {
  skip_if_not_installed("org.Hs.eg.db")
  msg <- testthat::capture_messages(df <- read_example("count"))
  expect_length(msg, 1L)
  expect_match(msg, "rsem_read: 3 samples, 44 genes -> 43 rows (species = human, expr = count, org.Hs.eg.db", fixed = TRUE)
  expect_match(msg, "with symbol: 42 genes (1 symbols merged from several IDs); kept as original IDs: 2 genes", fixed = TRUE)
  expect_s3_class(df, "data.frame")
  expect_equal(colnames(df), c("S1", "S2", "S3"))
  expect_false(anyDuplicated(rownames(df)) > 0)
  expect_true(all(vapply(df, is.numeric, logical(1))))
  expect_true(all(c("TP53", "GAPDH", "AAAS") %in% rownames(df)))
})

test_that("versioned Ensembl IDs are looked up, and the rest keep their ID", {
  skip_if_not_installed("org.Hs.eg.db")
  df <- suppressMessages(read_example("count"))
  raw <- example_raw("S1")
  expect_true(any(grepl("^ENSG[0-9]+\\.[0-9]+$", raw$gene_id)))
  expect_false(any(grepl("^ENSG[0-9]+", rownames(df))))
  expect_true(all(c("ERCC-00002", "10904") %in% rownames(df)))
  expect_equal(df["TP53", "S1"], raw$expected_count[grepl("^ENSG00000141510\\.", raw$gene_id)])
})

test_that("rows sharing a symbol are summed", {
  skip_if_not_installed("org.Hs.eg.db")
  raw <- example_raw("S2")
  pair <- raw$expected_count[sub("\\.[0-9]+$", "", raw$gene_id) %in% c("ENSG00000094914", "ENSG00000291836")]
  expect_length(pair, 2L)
  df <- suppressMessages(read_example("count"))
  expect_equal(df["AAAS", "S2"], sum(pair))

  tpm <- suppressMessages(read_example("tpm"))
  raw_tpm <- raw$TPM[sub("\\.[0-9]+$", "", raw$gene_id) %in% c("ENSG00000094914", "ENSG00000291836")]
  expect_equal(tpm["AAAS", "S2"], sum(raw_tpm))
})

test_that("nothing is lost: column sums equal the input", {
  skip_if_not_installed("org.Hs.eg.db")
  for (s in c("S1", "S2", "S3")) {
    raw <- example_raw(s)
    expect_equal(sum(suppressMessages(read_example("count"))[[s]]), sum(raw$expected_count))
    expect_equal(sum(suppressMessages(read_example("tpm"))[[s]]), sum(raw$TPM))
  }
})

test_that("count and tpm read different columns", {
  skip_if_not_installed("org.Hs.eg.db")
  expect_false(isTRUE(all.equal(suppressMessages(read_example("count")), suppressMessages(read_example("tpm")))))
})

test_that("values are not rounded", {
  skip_if_not_installed("org.Hs.eg.db")
  ids <- c("ENSG00000141510.18", "ENSG00000111640.5")
  v <- matrix(c(1.25, 2.5, 3.125, 4.75), 2, dimnames = list(ids, c("A", "B")))
  p <- make_dataset(v)
  df <- suppressMessages(rsem_read(p$data, p$metadata, species = "human", expr = "count"))
  expect_equal(df["TP53", "A"], 1.25)
  expect_equal(df["GAPDH", "B"], 4.75)
})

test_that("columns follow the metadata order and sample IDs stay text", {
  skip_if_not_installed("org.Hs.eg.db")
  ids <- c("ENSG00000141510.18", "ENSG00000111640.5")
  v <- values_matrix(ids, c("001", "002", "010"))
  p <- make_dataset(v)
  writeLines(
    c("sample_id,ref,group", "010,010.genes.results,x", "001,001.genes.results,y", "002,002.genes.results,x"),
    p$metadata
  )
  df <- suppressMessages(rsem_read(p$data, p$metadata, species = "human", expr = "count"))
  expect_equal(colnames(df), c("010", "001", "002"))
  expect_equal(df["TP53", "001"], unname(v["ENSG00000141510.18", "001"]))
})

test_that("mouse data use org.Mm.eg.db", {
  skip_if_not_installed("org.Mm.eg.db")
  ids <- c("ENSMUSG00000059552.7", "ERCC-00002")
  v <- values_matrix(ids)
  p <- make_dataset(v)
  df <- suppressMessages(rsem_read(p$data, p$metadata, species = "mouse", expr = "count"))
  expect_equal(rownames(df), c("Trp53", "ERCC-00002"))
  expect_error(
    rsem_read(p$data, p$metadata, species = "human", expr = "count"),
    "species = \"human\", but the gene IDs look like mouse"
  )
})

test_that("the species must match the gene IDs", {
  ids <- c("ENSG00000141510.18", "ENSG00000111640.5")
  p <- make_dataset(values_matrix(ids))
  expect_error(
    rsem_read(p$data, p$metadata, species = "mouse", expr = "count"),
    "looks? like human|look like human"
  )
  q <- make_dataset(values_matrix(c("geneA", "geneB")))
  expect_error(
    rsem_read(q$data, q$metadata, species = "human", expr = "count"),
    "do not look like Ensembl"
  )
})

test_that("Ensembl-like IDs unknown to the OrgDb are kept, not an error", {
  skip_if_not_installed("org.Hs.eg.db")
  ids <- c("ENSG00099999991.1", "ENSG00099999992.2")
  p <- make_dataset(values_matrix(ids))
  msg <- testthat::capture_messages(
    df <- rsem_read(p$data, p$metadata, species = "human", expr = "count")
  )
  expect_equal(rownames(df), ids)
  expect_match(msg, "with symbol: 0 genes", fixed = TRUE)
})

test_that("a gene ID equal to a gene symbol is not merged silently", {
  skip_if_not_installed("org.Hs.eg.db")
  p <- make_dataset(values_matrix(c("ENSG00000141510.18", "TP53")))
  expect_error(
    rsem_read(p$data, p$metadata, species = "human", expr = "count"),
    "would be merged"
  )
})

test_that("metadata problems are reported", {
  ids <- c("ENSG00000141510.18", "ENSG00000111640.5")
  p <- make_dataset(values_matrix(ids))
  call <- function() rsem_read(p$data, p$metadata, species = "human", expr = "count")

  writeLines(c("sample_id,path", "A,A.genes.results"), p$metadata)
  expect_error(call(), "missing: ref")
  writeLines(c("id,ref", "A,A.genes.results"), p$metadata)
  expect_error(call(), "missing: sample_id")
  writeLines(c("sample_id,ref", "A,A.genes.results", "A,B.genes.results"), p$metadata)
  expect_error(call(), "duplicated sample_id: A|duplicated sample_id values: A")
  writeLines(c("sample_id,ref", "A,A.genes.results", "B,A.genes.results"), p$metadata)
  expect_error(call(), "same file is listed")
  writeLines(c("sample_id,ref", ",A.genes.results"), p$metadata)
  expect_error(call(), "missing or empty sample_id")
  writeLines(c("sample_id,ref", "A,"), p$metadata)
  expect_error(call(), "missing or empty ref")
  writeLines(c("sample_id,ref", "A,nope.genes.results"), p$metadata)
  expect_error(call(), "1 file\\(s\\) not found.*A -> nope.genes.results")
  writeLines("sample_id,ref", p$metadata)
  expect_error(call(), "no rows")
})

test_that("RSEM file problems are reported with the file name", {
  ids <- c("ENSG00000141510.18", "ENSG00000111640.5", "ENSG00000075624.14")
  rewrite <- function(p, sample, fun) {
    f <- file.path(p$data, paste0(sample, ".genes.results"))
    tab <- utils::read.delim(f, check.names = FALSE, colClasses = c(gene_id = "character"))
    utils::write.table(fun(tab), f, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  call <- function(p) rsem_read(p$data, p$metadata, species = "human", expr = "count")

  p <- make_dataset(values_matrix(ids))
  rewrite(p, "B", function(t) t[, setdiff(names(t), "expected_count")])
  expect_error(call(p), "B.genes.results \\(sample B\\) is missing the column\\(s\\): expected_count")

  p <- make_dataset(values_matrix(ids))
  rewrite(p, "B", function(t) { t$expected_count <- as.character(t$expected_count); t$expected_count[2] <- "abc"; t })
  expect_error(call(p), "non-numeric")

  p <- make_dataset(values_matrix(ids))
  rewrite(p, "B", function(t) { t$expected_count[2] <- NA; t })
  expect_error(call(p), "missing values")

  p <- make_dataset(values_matrix(ids))
  rewrite(p, "A", function(t) { t$gene_id[2] <- t$gene_id[1]; t })
  expect_error(call(p), "duplicated gene_id")

  p <- make_dataset(values_matrix(ids))
  rewrite(p, "B", function(t) t[c(2, 1, 3), ])
  expect_error(call(p), "same genes in the same order as A.genes.results")

  p <- make_dataset(values_matrix(ids))
  rewrite(p, "B", function(t) { t$gene_id[1] <- "ENSG00000141510.17"; t })
  expect_error(call(p), "same genes in the same order")
})

test_that("extra columns in the RSEM files and the metadata are ignored", {
  skip_if_not_installed("org.Hs.eg.db")
  ids <- c("ENSG00000141510.18", "ENSG00000111640.5")
  p <- make_dataset(values_matrix(ids))
  f <- file.path(p$data, "A.genes.results")
  tab <- utils::read.delim(f, check.names = FALSE, colClasses = c(gene_id = "character"))
  tab$posterior_mean_count <- 1
  utils::write.table(tab, f, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(c("sample_id,ref,group", "A,A.genes.results,x", "B,B.genes.results,y"), p$metadata)
  df <- suppressMessages(rsem_read(p$data, p$metadata, species = "human", expr = "tpm"))
  expect_equal(dim(df), c(2L, 2L))
})

test_that("a missing OrgDb is explained", {
  # Only meaningful where the mouse package is absent.
  skip_if(requireNamespace("org.Mm.eg.db", quietly = TRUE))
  p <- make_dataset(values_matrix(c("ENSMUSG00000059552.7", "ENSMUSG00000064370.1")))
  expect_error(
    rsem_read(p$data, p$metadata, species = "mouse", expr = "count"),
    "BiocManager::install\\(\"org.Mm.eg.db\"\\)"
  )
})
