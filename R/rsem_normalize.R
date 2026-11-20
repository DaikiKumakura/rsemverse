#' Normalize between samples
#'
#' Makes samples comparable by correcting for sequencing depth and, for some
#' methods, RNA composition. The genes and samples stay the same; only the
#' values change. A message reports the method and the per-sample factors.
#'
#' @section Methods:
#' For counts (`expr = "count"`):
#' * `"cpm"`: counts per million, `count / library size * 1e6`. Corrects for
#'   depth only.
#' * `"tmm"`: CPM with TMM normalization factors (trimmed mean of M values;
#'   Robinson and Oshlack 2010), as `edgeR::calcNormFactors()` then
#'   `edgeR::cpm()`. Also corrects for differences in RNA composition, for
#'   example when a few genes take up a large part of one sample.
#' * `"tmmwsp"`: CPM with TMMwsp factors ("TMM with singleton pairing"), as
#'   `edgeR::calcNormFactors(method = "TMMwsp")`. A variant of TMM for data
#'   with many zeros: genes that are zero in one of the two samples compared
#'   are paired by size instead of discarded, and the reference sample is the
#'   one with the largest sum of square-root counts.
#' * `"upper_quartile"`: CPM with upper-quartile factors (Bullard et al.
#'   2010), as `edgeR::calcNormFactors(method = "upperquartile")`.
#' * `"deseq2"`: counts divided by the DESeq2 median-of-ratios size factors
#'   (Anders and Huber 2010), the values of `DESeq2::counts(dds, normalized =
#'   TRUE)`. Needs genes that are nonzero in every sample.
#' * `"quantile"`: quantile normalization, as `limma::normalizeQuantiles()`.
#'
#' For TPM (`expr = "tpm"`), TPM already corrects for depth and gene length,
#' so only `"quantile"` is offered. It makes every sample have the same
#' distribution of values, which is a strong assumption: use it only when the
#' samples are expected to differ in a few genes.
#'
#' The results are the same as those of the packages named above (checked in
#' the tests), computed here without them. Low-count genes affect TMM and the
#' size factors, so filter first ([rsem_filter()]). Methods for differential
#' expression (edgeR, DESeq2) need the raw counts, not these values. For a log
#' scale, pass the result to [rsem_transform()].
#'
#' @inheritParams rsem_filter
#' @param method For counts: `"cpm"`, `"tmm"`, `"tmmwsp"`, `"upper_quartile"`, `"deseq2"`
#'   or `"quantile"`. For TPM: `"quantile"`. Required.
#'
#' @return A data frame of normalized values in the same layout as `x`.
#'   The per-sample factors are in the attribute `"factors"`: library size
#'   for `"cpm"`, effective library size (library size times the
#'   normalization factor) for `"tmm"`, `"tmmwsp"` and `"upper_quartile"`, and size
#'   factors for `"deseq2"`.
#'
#' @references
#' Robinson MD, Oshlack A (2010). A scaling normalization method for
#' differential expression analysis of RNA-seq data. *Genome Biology* 11:R25.
#'
#' Anders S, Huber W (2010). Differential expression analysis for sequence
#' count data. *Genome Biology* 11:R106.
#'
#' Bullard JH, Purdom E, Hansen KD, Dudoit S (2010). Evaluation of
#' statistical methods for normalization and differential expression in
#' mRNA-Seq experiments. *BMC Bioinformatics* 11:94.
#'
#' Bolstad BM, Irizarry RA, Astrand M, Speed TP (2003). A comparison of
#' normalization methods for high density oligonucleotide array data based
#' on variance and bias. *Bioinformatics* 19:185-193.
#'
#' @examples
#' if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   ex <- system.file("extdata", "example", package = "rsemverse")
#'   counts <- rsem_read(file.path(ex, "data"), file.path(ex, "metadata.csv"),
#'                       species = "human", expr = "count")
#'   tmm <- rsem_normalize(counts, expr = "count", method = "tmm")
#'   attr(tmm, "factors")
#'   rsem_transform(tmm, "log2")    # log2(TMM-CPM + 1)
#' }
#' @export
rsem_normalize <- function(x, expr, method) {
  if (missing(expr)) .stopf("`expr` must be \"count\" or \"tpm\".")
  expr <- .choice(expr, c("count", "tpm"), "expr")
  if (missing(method)) .stopf("`method` is required; see ?rsem_normalize.")
  allowed <- if (expr == "count") c("cpm", "tmm", "tmmwsp", "upper_quartile", "deseq2", "quantile") else "quantile"
  if (is.character(method) && length(method) == 1L && !is.na(method) && expr == "tpm" &&
      method %in% c("cpm", "tmm", "tmmwsp", "upper_quartile", "deseq2")) {
    .stopf("method = \"%s\" needs counts (expr = \"count\"); TPM already corrects for depth. For TPM, use \"quantile\".", method)
  }
  method <- .choice(method, allowed, "method")
  m <- .expr_matrix(x)
  lib <- colSums(m)
  if (method != "quantile" && any(lib <= 0)) {
    .stopf("Every sample needs a nonzero total count; zero: %s", paste(colnames(m)[lib <= 0], collapse = ", "))
  }

  factors <- NULL
  out <- switch(method,
    cpm = {
      factors <- lib
      sweep(m, 2L, lib, "/") * 1e6
    },
    tmm = ,
    tmmwsp = ,
    upper_quartile = {
      nf <- switch(method, tmm = .tmm_factors(m, lib), tmmwsp = .tmmwsp_factors(m, lib),
                   upper_quartile = .upper_quartile_factors(m, lib))
      factors <- lib * nf
      sweep(m, 2L, factors, "/") * 1e6
    },
    deseq2 = {
      factors <- .median_ratio_factors(m)
      sweep(m, 2L, factors, "/")
    },
    quantile = .quantile_normalize(m)
  )
  dimnames(out) <- dimnames(m)
  if (!is.null(factors)) factors <- stats::setNames(as.numeric(factors), colnames(m))

  message(sprintf(
    "rsem_normalize (%s, %s): %s genes x %d samples%s",
    method, expr, .fmt(nrow(m)), ncol(m),
    if (is.null(factors)) "" else paste0("; ", c(
      cpm = "library sizes", tmm = "effective library sizes", tmmwsp = "effective library sizes", upper_quartile = "effective library sizes",
      deseq2 = "size factors"
    )[[method]], ": ", paste(sprintf("%s %s", names(factors), signif(factors, 4)), collapse = ", "))
  ))
  res <- .as_table(out, from = x)
  if (!is.null(factors)) attr(res, "factors") <- factors
  res
}

# TMM normalization factors (Robinson and Oshlack 2010), with the defaults of
# edgeR::calcNormFactors(): 30% trim of log ratios, 5% trim of mean
# expression, precision weights, and the reference sample whose upper
# quartile is closest to the mean upper quartile.
.tmm_factors <- function(m, lib) {
  n <- ncol(m)
  x <- m[rowSums(m > 0) > 0, , drop = FALSE]
  if (!nrow(x) || n == 1L) return(rep(1, n))
  uq <- apply(x, 2L, stats::quantile, probs = 0.75) / lib
  ref <- if (stats::median(uq) < 1e-20) which.max(colSums(sqrt(x))) else which.min(abs(uq - mean(uq)))
  f <- vapply(seq_len(n), function(i) {
    .tmm_one(x[, i], x[, ref], lib[[i]], lib[[ref]])
  }, numeric(1))
  f / exp(mean(log(f)))
}

.tmm_one <- function(obs, ref, n_obs, n_ref, logratio_trim = 0.3, sum_trim = 0.05) {
  p_obs <- obs / n_obs
  p_ref <- ref / n_ref
  m_val <- log2(p_obs / p_ref)
  a_val <- (log2(p_obs) + log2(p_ref)) / 2
  w_inv <- (n_obs - obs) / n_obs / obs + (n_ref - ref) / n_ref / ref
  ok <- is.finite(m_val) & is.finite(a_val)
  m_val <- m_val[ok]
  a_val <- a_val[ok]
  w_inv <- w_inv[ok]
  if (!length(m_val) || max(abs(m_val)) < 1e-6) return(1)
  g <- length(m_val)
  lo_m <- floor(g * logratio_trim) + 1
  lo_a <- floor(g * sum_trim) + 1
  rm <- rank(m_val)
  ra <- rank(a_val)
  keep <- rm >= lo_m & rm <= g + 1 - lo_m & ra >= lo_a & ra <= g + 1 - lo_a
  f <- sum(m_val[keep] / w_inv[keep], na.rm = TRUE) / sum(1 / w_inv[keep], na.rm = TRUE)
  if (is.na(f)) f <- 0
  2^f
}

# TMMwsp factors (edgeR::calcNormFactors(method = "TMMwsp"), edgeR 4.x):
# the reference is the sample with the largest sum of square-root counts;
# genes zero in both samples are dropped, and genes zero in only one sample
# are paired by decreasing size instead of being dropped.
.tmmwsp_factors <- function(m, lib) {
  n <- ncol(m)
  x <- m[rowSums(m > 0) > 0, , drop = FALSE]
  if (!nrow(x) || n == 1L) return(rep(1, n))
  ref <- which.max(colSums(sqrt(x)))
  f <- vapply(seq_len(n), function(i) .tmmwsp_one(x[, i], x[, ref], lib[[i]], lib[[ref]]), numeric(1))
  f / exp(mean(log(f)))
}

.tmmwsp_one <- function(obs, ref, n_obs, n_ref, logratio_trim = 0.3, sum_trim = 0.05) {
  eps <- 1e-14
  npos <- 2L * (obs > eps) + (ref > eps)
  drop <- npos == 0L
  obs <- obs[!drop]; ref <- ref[!drop]; npos <- npos[!drop]
  zero_obs <- npos == 1L
  zero_ref <- npos == 2L
  k <- zero_obs | zero_ref
  singles <- min(sum(zero_obs), sum(zero_ref))
  if (singles > 0L) {
    obs <- c(obs[!k], sort(obs[k], decreasing = TRUE)[seq_len(singles)])
    ref <- c(ref[!k], sort(ref[k], decreasing = TRUE)[seq_len(singles)])
  } else {
    obs <- obs[!k]; ref <- ref[!k]
  }
  g <- length(obs)
  if (!g) return(1)
  p_obs <- obs / n_obs
  p_ref <- ref / n_ref
  m_val <- log2(p_obs / p_ref)
  a_val <- 0.5 * log2(p_obs * p_ref)
  if (max(abs(m_val)) < 1e-6) return(1)
  m_shrunk <- log2(((obs + 0.5) / (n_obs + 0.5)) / ((ref + 0.5) / (n_ref + 0.5)))
  o_m <- order(m_val, m_shrunk)
  o_a <- order(a_val)
  lo_m <- as.integer(g * logratio_trim) + 1L
  lo_a <- as.integer(g * sum_trim) + 1L
  keep_m <- keep_a <- rep_len(FALSE, g)
  keep_m[o_m[lo_m:(g + 1L - lo_m)]] <- TRUE
  keep_a[o_a[lo_a:(g + 1L - lo_a)]] <- TRUE
  keep <- keep_m & keep_a
  p_obs <- p_obs[keep]
  p_ref <- p_ref[keep]
  v <- (1 - p_obs) / p_obs / n_obs + (1 - p_ref) / p_ref / n_ref
  w <- (1 + 1e-6) / (v + 1e-6)
  2^(sum(w * m_val[keep]) / sum(w))
}

# Upper-quartile factors: the 75th percentile of the counts of genes that are
# nonzero in at least one sample, divided by the library size, scaled to a
# geometric mean of 1.
.upper_quartile_factors <- function(m, lib) {
  x <- m[rowSums(m > 0) > 0, , drop = FALSE]
  if (!nrow(x) || ncol(m) == 1L) return(rep(1, ncol(m)))
  f <- apply(x, 2L, stats::quantile, probs = 0.75) / lib
  if (any(f <= 0)) {
    .stopf("The upper quartile of the counts is zero in: %s. Filter out low-count genes first, or use \"tmm\".", paste(colnames(m)[f <= 0], collapse = ", "))
  }
  f / exp(mean(log(f)))
}

# Median-of-ratios size factors (Anders and Huber 2010): for each sample, the
# median ratio of its counts to the geometric mean of each gene over samples,
# using genes that are nonzero in every sample.
.median_ratio_factors <- function(m) {
  log_geo <- rowMeans(log(m))
  usable <- is.finite(log_geo)
  if (!any(usable)) {
    .stopf("method = \"deseq2\" needs genes that are nonzero in every sample; every gene has a zero somewhere.")
  }
  vapply(seq_len(ncol(m)), function(j) {
    exp(stats::median(log(m[usable, j]) - log_geo[usable]))
  }, numeric(1))
}

# Quantile normalization with tied values sharing the average of their
# positions (as limma::normalizeQuantiles(ties = TRUE)).
.quantile_normalize <- function(m) {
  if (ncol(m) == 1L) return(m)
  g <- nrow(m)
  target <- rowMeans(apply(m, 2L, sort))
  if (g == 1L) return(matrix(target, 1L, ncol(m)))
  pos <- (seq_len(g) - 1) / (g - 1)
  apply(m, 2L, function(v) {
    stats::approx(pos, target, (rank(v) - 1) / (g - 1), ties = list("ordered", mean))$y
  })
}
