#' Sample-level quality control: PCA, sample correlation and a summary
#'
#' Returns the numbers behind the usual first look at a gene by sample table:
#' principal components, correlation between samples and, for each sample,
#' its total, the number of detected genes and how well it agrees with the
#' other samples. Nothing is plotted or removed; the results are tables that
#' can be passed to any plotting code.
#'
#' @section Genes used:
#' The principal components use the `n_top` genes with the largest variance of
#' the transformed values (genes that are zero everywhere are never used), so
#' that the components reflect the genes that vary. The correlation uses the
#' genes chosen by `correlation_genes`, which is a separate choice because
#' restricting the correlation to the most variable genes changes the
#' correlations a great deal:
#' * `"detected"` (the default): genes with a value of at least `min_value`
#'   in at least one sample;
#' * `"all"`: every gene that is not zero everywhere;
#' * `"variable"`: the same genes as the principal components.
#'
#' @section Choices that are not made silently:
#' * `transform = "log2"` (the default) uses `log2(value + 1)`; `"none"` uses
#'   the values as they are (for values that are already on a log scale). A
#'   table made with [rsem_transform()] remembers that it is log2-transformed,
#'   and `transform = "log2"` then stops instead of transforming twice.
#'   Rank correlation does not depend on the transformation.
#' * The values are centred per gene before the PCA, but not scaled
#'   (`pca_scale = FALSE`): scaling gives every gene the same weight, so lowly
#'   expressed, noisy genes count as much as highly expressed ones.
#' * Samples are not flagged or removed. `median_correlation` and the other
#'   columns are numbers to look at, and what counts as an outlier depends on
#'   the study.
#'
#' @param x The gene by sample table, as returned by [rsem_read()],
#'   [rsem_filter()] or [rsem_normalize()]. Needs at least three samples.
#' @param metadata Optional sample table (path to a CSV or a data frame with a
#'   `sample_id` column), as in [rsem_long()]; its columns are added to
#'   `pca_scores` and `sample_summary`.
#' @param transform `"log2"` (the default) or `"none"`.
#' @param pseudocount For `"log2"`. Default 1.
#' @param min_value Value a gene must reach to count as detected. Default 1
#'   (TPM 1 or count 1).
#' @param n_top Number of genes for the PCA. Default 2000.
#' @param n_pcs Number of components to return. Default 10 (fewer when there
#'   are fewer).
#' @param pca_scale Scale genes to unit variance before the PCA. Default
#'   `FALSE`.
#' @param correlation `"spearman"` (the default), `"pearson"` or `"kendall"`
#'   (see [stats::cor()]).
#' @param correlation_genes `"detected"`, `"all"` or `"variable"`; see Genes
#'   used.
#'
#' @return A list of class `rsem_sample_qc`:
#'   * `pca_scores`: one row per sample, with `sample_id`, the components
#'     `PC1`, `PC2`, ... and the metadata columns;
#'   * `pca_variance`: `pc`, `sdev`, `variance`, `proportion` and
#'     `cumulative` (of the variance of the genes used);
#'   * `pca_loadings`: genes by components (the rotation of
#'     [stats::prcomp()]);
#'   * `correlation_matrix`: samples by samples;
#'   * `sample_distance`: `1 - correlation_matrix`;
#'   * `sample_summary`: per sample `total` (sum of the values), `detected_genes`
#'     and `detected_share` (of all genes in `x`), `median_correlation`,
#'     `min_correlation` and `most_similar` (the sample with the highest
#'     correlation), `rle_median` and `rle_iqr` (relative log expression: the
#'     median and the interquartile range of the deviations of a sample's log2
#'     values from the median of each gene over the samples, for genes detected
#'     in at least one sample; around 0 and similar between samples for
#'     comparable libraries; `NA` unless the values are log2) and the metadata
#'     columns;
#'   * `parameters`: the settings and the numbers of genes used.
#'
#' @examples
#' set.seed(1)
#' m <- matrix(rexp(200 * 6, 0.05), 200, 6,
#'             dimnames = list(paste0("g", 1:200), paste0("S", 1:6)))
#' qc <- rsem_sample_qc(as.data.frame(m), n_top = 100)
#' qc$pca_variance[1:3, ]
#' qc$sample_summary
#' @export
rsem_sample_qc <- function(x, metadata = NULL, transform = "log2", pseudocount = 1, min_value = 1,
                           n_top = 2000, n_pcs = 10, pca_scale = FALSE,
                           correlation = "spearman", correlation_genes = "detected") {
  transform <- .choice(transform, c("log2", "none"), "transform")
  correlation <- .choice(correlation, c("spearman", "pearson", "kendall"), "correlation")
  correlation_genes <- .choice(correlation_genes, c("detected", "all", "variable"), "correlation_genes")
  min_value <- .number(min_value, "min_value", lower = 0)
  n_top <- .number(n_top, "n_top", lower = 2, integer = TRUE)
  n_pcs <- .number(n_pcs, "n_pcs", lower = 1, integer = TRUE)
  if (!is.logical(pca_scale) || length(pca_scale) != 1L || is.na(pca_scale)) .stopf("`pca_scale` must be TRUE or FALSE.")
  if (transform == "log2" && "log2" %in% attr(x, "transform", exact = TRUE)) {
    .stopf("`x` is already log2-transformed (rsem_transform()); use transform = \"none\" so that it is not transformed twice.")
  }
  m <- .expr_matrix(x, nonnegative = transform == "log2")
  n <- ncol(m)
  if (n < 3L) .stopf("rsem_sample_qc() needs at least three samples; `x` has %d.", n)
  if (transform == "log2") {
    pseudocount <- .number(pseudocount, "pseudocount", lower = 0)
    if (pseudocount == 0 && any(m == 0)) .stopf("log2 with pseudocount = 0 would turn zeros into -Inf; use a positive pseudocount.")
  }
  meta <- NULL
  if (!is.null(metadata)) {
    meta <- .read_sample_table(metadata)
    absent <- setdiff(colnames(m), meta$sample_id)
    if (length(absent)) .stopf("Sample(s) of `x` not in the sample table: %s", paste(utils::head(absent, 5L), collapse = ", "))
    meta <- meta[match(colnames(m), meta$sample_id), , drop = FALSE]
    rownames(meta) <- NULL
  }

  v <- if (transform == "log2") log2(m + pseudocount) else m
  nonzero <- rowSums(m != 0) > 0
  vars <- apply(v, 1L, stats::var)
  vars[!nonzero | !is.finite(vars)] <- -Inf
  candidates <- which(vars > 0)
  if (length(candidates) < 2L) .stopf("Fewer than two genes vary between the samples; there is nothing to analyse.")
  pca_idx <- candidates[order(-vars[candidates], seq_along(candidates))][seq_len(min(n_top, length(candidates)))]
  pca_idx <- sort(pca_idx)

  pc <- stats::prcomp(t(v[pca_idx, , drop = FALSE]), center = TRUE, scale. = pca_scale)
  k <- min(n_pcs, ncol(pc$x))
  pc_names <- paste0("PC", seq_len(k))
  scores <- data.frame(sample_id = colnames(m), pc$x[, seq_len(k), drop = FALSE], stringsAsFactors = FALSE)
  names(scores)[-1L] <- pc_names
  variance <- pc$sdev^2
  pca_variance <- data.frame(
    pc = paste0("PC", seq_along(variance)), sdev = pc$sdev, variance = variance,
    proportion = variance / sum(variance), cumulative = cumsum(variance) / sum(variance),
    stringsAsFactors = FALSE
  )[seq_len(k), , drop = FALSE]
  loadings <- pc$rotation[, seq_len(k), drop = FALSE]
  colnames(loadings) <- pc_names

  cor_idx <- switch(correlation_genes,
    detected = which(rowSums(m >= min_value) >= 1L),
    all = which(nonzero),
    variable = pca_idx
  )
  if (length(cor_idx) < 3L) .stopf("Fewer than three genes are available for the correlation (correlation_genes = \"%s\").", correlation_genes)
  cmat <- suppressWarnings(stats::cor(v[cor_idx, , drop = FALSE], method = correlation))
  dimnames(cmat) <- list(colnames(m), colnames(m))
  # a sample with the same value in every gene has no correlation, not even with itself
  flat <- apply(v[cor_idx, , drop = FALSE], 2L, function(z) length(unique(z)) < 2L)
  cmat[flat, ] <- NA_real_
  cmat[, flat] <- NA_real_
  off <- cmat
  diag(off) <- NA_real_
  # Relative log expression (Gandolfo & Speed 2018): the deviation of each gene in a sample from
  # the median of that gene over the samples, on the log scale; the boxes of clean data sit at 0 and
  # are equally wide. Only on log2 values, and only for genes detected in at least one sample.
  rle_median <- rle_iqr <- rep(NA_real_, n)
  if (transform == "log2" || "log2" %in% attr(x, "transform", exact = TRUE)) {
    det_idx <- which(rowSums(m >= min_value) >= 1L)
    if (length(det_idx) >= 3L) {
      vd <- v[det_idx, , drop = FALSE]
      dev <- vd - apply(vd, 1L, stats::median)
      rle_median <- apply(dev, 2L, stats::median)
      rle_iqr <- apply(dev, 2L, stats::IQR)
    }
  }
  summary <- data.frame(
    sample_id = colnames(m), total = colSums(m), detected_genes = colSums(m >= min_value),
    detected_share = colSums(m >= min_value) / nrow(m),
    median_correlation = apply(off, 2L, stats::median, na.rm = TRUE),
    min_correlation = apply(off, 2L, function(z) if (all(is.na(z))) NA_real_ else min(z, na.rm = TRUE)),
    most_similar = apply(off, 2L, function(z) if (all(is.na(z))) NA_character_ else colnames(m)[which.max(z)]),
    rle_median = unname(rle_median), rle_iqr = unname(rle_iqr),
    stringsAsFactors = FALSE
  )
  rownames(summary) <- NULL
  if (!is.null(meta)) {
    extra <- meta[, setdiff(names(meta), "sample_id"), drop = FALSE]
    if (ncol(extra)) {
      scores <- cbind(scores, extra)
      summary <- cbind(summary, extra)
    }
  }
  parameters <- list(
    transform = transform, pseudocount = if (transform == "log2") pseudocount else NA_real_, min_value = min_value,
    n_top = n_top, pca_genes = length(pca_idx), pca_scale = pca_scale, correlation = correlation,
    correlation_genes = correlation_genes, correlation_gene_count = length(cor_idx),
    n_genes = nrow(m), n_samples = n
  )
  message(sprintf(
    "rsem_sample_qc: %s genes x %d samples; PCA on %s genes (PC1 %s%%, PC2 %s%% of their variance); %s correlation on %s genes (%s), median between samples %s",
    .fmt(nrow(m)), n, .fmt(length(pca_idx)), format(round(100 * pca_variance$proportion[1L], 1)),
    format(round(100 * pca_variance$proportion[min(2L, k)], 1)), correlation, .fmt(length(cor_idx)), correlation_genes,
    format(signif(stats::median(off, na.rm = TRUE), 3))
  ))
  structure(
    list(pca_scores = scores, pca_variance = pca_variance, pca_loadings = loadings, correlation_matrix = cmat,
         sample_distance = 1 - cmat, sample_summary = summary, parameters = parameters),
    class = "rsem_sample_qc"
  )
}

#' @export
print.rsem_sample_qc <- function(x, ...) {
  p <- x$parameters
  cat(sprintf(
    "rsem_sample_qc: %s genes x %d samples\n  PCA: %s genes, %s; first components explain %s%%\n  %s correlation on %s genes (%s)\n",
    .fmt(p$n_genes), p$n_samples, .fmt(p$pca_genes), if (p$pca_scale) "scaled" else "centred, not scaled",
    paste(round(100 * x$pca_variance$proportion[seq_len(min(3L, nrow(x$pca_variance)))], 1), collapse = ", "),
    p$correlation, .fmt(p$correlation_gene_count), p$correlation_genes
  ))
  cat("  elements: pca_scores, pca_variance, pca_loadings, correlation_matrix, sample_distance, sample_summary, parameters\n")
  invisible(x)
}
