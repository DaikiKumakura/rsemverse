example_dir <- function() {
  system.file("extdata", "example", package = "rsemverse")
}

# Read the example directly, without rsem_read(), for independent expectations.
example_raw <- function(sample = "S1") {
  utils::read.delim(
    file.path(example_dir(), "data", paste0(sample, ".genes.results")),
    check.names = FALSE, colClasses = c(gene_id = "character")
  )
}

# Write RSEM files for `values` (genes x samples; rows named by gene_id) and a
# metadata CSV that lists them. Returns the paths.
make_dataset <- function(values, column = "expected_count", dir = tempfile("rsemverse-")) {
  dir.create(file.path(dir, "data"), recursive = TRUE)
  for (s in colnames(values)) {
    tab <- data.frame(
      gene_id = rownames(values),
      `transcript_id(s)` = "T",
      length = 1000,
      effective_length = 800,
      expected_count = 0,
      TPM = 0,
      FPKM = 0,
      check.names = FALSE
    )
    tab[[column]] <- values[, s]
    utils::write.table(
      tab, file.path(dir, "data", paste0(s, ".genes.results")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
  writeLines(
    c("sample_id,ref", paste0(colnames(values), ",", colnames(values), ".genes.results")),
    file.path(dir, "metadata.csv")
  )
  list(data = file.path(dir, "data"), metadata = file.path(dir, "metadata.csv"), dir = dir)
}

values_matrix <- function(ids, samples = c("A", "B"), seed = 1L) {
  set.seed(seed)
  m <- matrix(
    round(stats::runif(length(ids) * length(samples), 0, 500), 2),
    length(ids), length(samples), dimnames = list(ids, samples)
  )
  m
}

read_example <- function(expr = "count", species = "human") {
  rsem_read(
    data = file.path(example_dir(), "data"),
    metadata = file.path(example_dir(), "metadata.csv"),
    species = species, expr = expr
  )
}
