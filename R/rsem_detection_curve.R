#' How many genes each detection threshold keeps
#'
#' Counts the genes that reach each threshold in enough samples, for a grid of
#' thresholds and sample requirements, so that the effect of a threshold can
#' be seen before [rsem_filter()] removes genes. Nothing is removed.
#'
#' For a threshold \eqn{\tau} and a required number of samples \eqn{k}, the
#' count is the number of genes with a value of at least \eqn{\tau} in at
#' least \eqn{k} samples. This is exactly the rule of
#' `rsem_filter(method = "detected", min_value = tau, min_samples = k)`, so
#' each row tells how many genes that call would keep. The column
#' `genes_mean` gives the number of genes whose mean over samples reaches the
#' threshold, the other rule often used in published scripts.
#'
#' @inheritParams rsem_filter
#' @param thresholds Values to try. Default `c(0.1, 0.5, 1, 2, 5, 10)` for
#'   TPM and `c(1, 5, 10, 20, 50, 100)` for counts.
#' @param min_samples Sample requirements to try; a number below 1 is a
#'   proportion of the samples, as in [rsem_filter()]. Default `c(1, 0.2,
#'   0.5)`: at least one sample, 20% and half of the samples.
#'
#' @return A data frame with one row per threshold and sample requirement:
#'   `threshold`, `min_samples` (as given), `samples_needed`, `genes` (genes
#'   kept), `share` (of all genes in `x`) and `genes_mean` (genes whose mean
#'   reaches the threshold; the same in every row of a threshold).
#'
#' @examples
#' if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   ex <- system.file("extdata", "example", package = "rsemverse")
#'   tpm <- rsem_read(file.path(ex, "data"), file.path(ex, "metadata.csv"),
#'                    species = "human", expr = "tpm")
#'   rsem_detection_curve(tpm, expr = "tpm")
#' }
#' @export
rsem_detection_curve <- function(x, expr, thresholds = NULL, min_samples = c(1, 0.2, 0.5)) {
  if (missing(expr)) .stopf("`expr` must be \"count\" or \"tpm\".")
  expr <- .choice(expr, c("count", "tpm"), "expr")
  m <- .expr_matrix(x)
  n <- ncol(m)
  if (is.null(thresholds)) thresholds <- if (expr == "tpm") c(0.1, 0.5, 1, 2, 5, 10) else c(1, 5, 10, 20, 50, 100)
  if (!is.numeric(thresholds) || !length(thresholds) || anyNA(thresholds) || any(thresholds < 0)) {
    .stopf("`thresholds` must be numbers that are not negative.")
  }
  if (!is.numeric(min_samples) || !length(min_samples)) .stopf("`min_samples` must be one or more numbers.")
  needed <- vapply(min_samples, .needed_samples, numeric(1), n = n)
  thresholds <- sort(unique(thresholds))
  mu <- rowMeans(m)

  rows <- lapply(thresholds, function(t) {
    passing <- rowSums(m >= t)
    data.frame(
      threshold = t, min_samples = min_samples, samples_needed = needed,
      genes = vapply(needed, function(k) sum(passing >= k), integer(1)),
      genes_mean = sum(mu >= t), stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$share <- out$genes / nrow(m)
  out <- out[, c("threshold", "min_samples", "samples_needed", "genes", "share", "genes_mean")]
  rownames(out) <- NULL
  message(sprintf(
    "rsem_detection_curve: %s genes x %d samples; %d thresholds x %d sample requirements (%s)",
    .fmt(nrow(m)), n, length(thresholds), length(min_samples), if (expr == "tpm") "TPM" else "count"
  ))
  out
}
