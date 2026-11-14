# Records reference results from edgeR, DESeq2 and limma for the tests of
# rsem_filter() and rsem_normalize(). The package itself does not use these
# packages. Run from the package root: Rscript data-raw/make_golden.R
SEED <- 20261114
N_GENES <- 1000

simulate_counts <- function(n_samples, seed) {
  set.seed(seed)
  mu <- exp(stats::rnorm(N_GENES, mean = 3, sd = 2.2))
  lib <- stats::runif(n_samples, 0.5, 2)
  m <- vapply(seq_len(n_samples), function(j) {
    stats::rnbinom(N_GENES, mu = mu * lib[j], size = 5)
  }, numeric(N_GENES))
  # Composition difference: a few genes take up much more of two samples.
  m[1:20, 1:2] <- m[1:20, 1:2] * 30
  m[sample(N_GENES, 150), 1] <- 0
  dimnames(m) <- list(sprintf("G%04d", seq_len(N_GENES)), sprintf("S%02d", seq_len(n_samples)))
  m
}

reference <- function(m) {
  dds_sf <- DESeq2::estimateSizeFactorsForMatrix(m)
  tmm_nf <- edgeR::calcNormFactors(m, method = "TMM")
  uq_nf <- edgeR::calcNormFactors(m, method = "upperquartile")
  list(
    counts = m,
    filter_edger = suppressMessages(edgeR::filterByExpr(m)),
    filter_edger_min5 = suppressMessages(edgeR::filterByExpr(m, min.count = 5)),
    cpm = edgeR::cpm(m),
    tmm = edgeR::cpm(m, lib.size = colSums(m) * tmm_nf),
    tmm_factors = tmm_nf,
    upper_quartile = edgeR::cpm(m, lib.size = colSums(m) * uq_nf),
    deseq2 = sweep(m, 2L, dds_sf, "/"),
    deseq2_size_factors = dds_sf,
    quantile = limma::normalizeQuantiles(m)
  )
}

golden <- list(
  six = reference(simulate_counts(6, SEED)),
  twelve = reference(simulate_counts(12, SEED + 1)),
  versions = c(
    edgeR = as.character(utils::packageVersion("edgeR")),
    DESeq2 = as.character(utils::packageVersion("DESeq2")),
    limma = as.character(utils::packageVersion("limma"))
  )
)
saveRDS(golden, file.path("tests", "testthat", "fixtures", "golden.rds"), version = 3)
print(golden$versions)
cat("kept by filterByExpr:", sum(golden$six$filter_edger), "of", N_GENES, "(6 samples);",
    sum(golden$twelve$filter_edger), "(12 samples)\n")
