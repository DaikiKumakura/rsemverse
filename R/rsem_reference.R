#' Z-scores against a fixed reference group
#'
#' [rsem_fit_reference()] computes each gene's mean and standard deviation in
#' a reference group of samples (for example controls or vehicle), and
#' [rsem_apply_reference()] expresses any sample as the number of reference
#' standard deviations it lies from the reference mean:
#' \deqn{z_{gs} = (x_{gs} - \mu_{g,R}) / \sigma_{g,R}}
#' where \eqn{\mu_{g,R}} and \eqn{\sigma_{g,R}} (with the \eqn{n - 1}
#' denominator) are taken over the reference samples \eqn{R} only.
#' [rsem_reference_scale()] does both on one table.
#'
#' Unlike `rsem_transform(x, "zscore")`, which centres on the mean of all
#' samples in the table, the reference stays fixed: adding or removing other
#' samples does not change anyone's value, and the fitted reference can be
#' applied to new data (for example a validation cohort, without fitting
#' anything on it).
#'
#' @section Choices that are not made silently:
#' * A gene with no variation in the reference (standard deviation 0, for
#'   example zero in every reference sample) has no z-score. With
#'   `zero_sd = "na"` (the default) its row is `NA`; with `"drop"` it is
#'   removed. Nothing is smoothed or borrowed from other genes.
#' * The reference needs at least two samples; one sample has no standard
#'   deviation.
#' * Genes of the reference missing from the new table are reported; genes of
#'   the new table missing from the reference are removed and reported.
#' * The new table must be on the same scale as the table the reference was
#'   fitted on (the same `expr` and normalization); the transformation stored
#'   in the reference (`log2` with its pseudocount, or none) is applied again.
#'
#' @param x The gene by sample table, as returned by [rsem_read()],
#'   [rsem_filter()] or [rsem_normalize()], on the original scale (not
#'   log-transformed) when `transform = "log2"`.
#' @param reference For [rsem_fit_reference()]: the reference samples, as
#'   sample IDs, or, with `metadata` and `group`, as one or more values of the
#'   column `group`. For [rsem_apply_reference()]: the object returned by
#'   [rsem_fit_reference()].
#' @param metadata Optional sample table (path to a CSV or a data frame with a
#'   `sample_id` column), as in [rsem_long()].
#' @param group Column of `metadata` that holds the groups.
#' @param transform `"log2"` (the default; `log2(value + pseudocount)` before
#'   the mean and standard deviation) or `"none"`.
#' @param pseudocount For `"log2"`. Default 1.
#' @param zero_sd `"na"` (the default) or `"drop"`, for genes whose reference
#'   standard deviation is 0.
#' @param ... Passed from [rsem_reference_scale()] to [rsem_fit_reference()].
#'
#' @return [rsem_fit_reference()]: an object of class `rsem_reference`, a list
#'   with `mean`, `sd` (named by gene), `samples` (the reference sample IDs),
#'   `transform` and `pseudocount`. [rsem_apply_reference()] and
#'   [rsem_reference_scale()]: a data frame of z-scores in the layout of `x`,
#'   with the reference in the attribute `"reference"`.
#'
#' @examples
#' if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   ex <- system.file("extdata", "example", package = "rsemverse")
#'   tpm <- rsem_read(file.path(ex, "data"), file.path(ex, "metadata.csv"),
#'                    species = "human", expr = "tpm")
#'   ref <- rsem_fit_reference(tpm, reference = c("S1", "S2"))
#'   z <- rsem_apply_reference(tpm, ref)
#'   head(z)
#' }
#' @export
rsem_fit_reference <- function(x, reference, metadata = NULL, group = NULL, transform = "log2", pseudocount = 1) {
  if (missing(reference)) .stopf("`reference` is required: sample IDs, or group values with `metadata` and `group`.")
  transform <- .choice(transform, c("log2", "none"), "transform")
  m <- .expr_matrix(x, nonnegative = transform == "log2")
  if (!is.character(reference) || !length(reference) || anyNA(reference)) {
    .stopf("`reference` must be one or more sample IDs or group values.")
  }
  if (!is.null(metadata) || !is.null(group)) {
    if (is.null(metadata) || is.null(group)) .stopf("`metadata` and `group` must be given together.")
    meta <- .read_sample_table(metadata)
    if (!is.character(group) || length(group) != 1L || !group %in% names(meta)) {
      .stopf("`group` must be a column of the sample table.")
    }
    absent <- setdiff(reference, as.character(meta[[group]]))
    if (length(absent)) .stopf("Value(s) not found in column '%s': %s", group, paste(absent, collapse = ", "))
    samples <- meta$sample_id[as.character(meta[[group]]) %in% reference]
    samples <- samples[samples %in% colnames(m)]
  } else {
    absent <- setdiff(reference, colnames(m))
    if (length(absent)) .stopf("Reference sample(s) not in `x`: %s", paste(absent, collapse = ", "))
    samples <- unique(reference)
  }
  if (length(samples) < 2L) {
    .stopf("The reference has %d sample(s) in `x`; at least 2 are needed for a standard deviation.", length(samples))
  }
  pseudocount <- .number(pseudocount, "pseudocount", lower = 0)
  r <- .ref_transform(m[, samples, drop = FALSE], transform, pseudocount)
  mu <- rowMeans(r)
  s <- sqrt(rowSums((r - mu)^2) / (ncol(r) - 1))
  s[s <= sqrt(.Machine$double.eps) * pmax(1, abs(mu))] <- 0
  message(sprintf(
    "rsem_fit_reference: %s genes, %d reference samples (%s); %s; %s genes with standard deviation 0",
    .fmt(length(mu)), length(samples), paste(samples, collapse = ", "),
    if (transform == "log2") sprintf("log2(value + %s)", format(pseudocount)) else "no transformation",
    .fmt(sum(s == 0))
  ))
  structure(
    list(mean = mu, sd = s, samples = samples, transform = transform, pseudocount = pseudocount),
    class = "rsem_reference"
  )
}

#' @rdname rsem_fit_reference
#' @export
rsem_apply_reference <- function(x, reference, zero_sd = "na") {
  if (!inherits(reference, "rsem_reference")) .stopf("`reference` must be the result of rsem_fit_reference().")
  zero_sd <- .choice(zero_sd, c("na", "drop"), "zero_sd")
  m <- .expr_matrix(x, nonnegative = reference$transform == "log2")
  genes <- names(reference$mean)
  common <- rownames(m)[rownames(m) %in% genes]
  if (!length(common)) .stopf("No gene of `x` is in the reference.")
  dropped <- nrow(m) - length(common)
  missing_ref <- length(setdiff(genes, rownames(m)))
  r <- .ref_transform(m[common, , drop = FALSE], reference$transform, reference$pseudocount)
  s <- reference$sd[common]
  z <- (r - reference$mean[common]) / s
  z[s == 0, ] <- NA_real_
  if (zero_sd == "drop") z <- z[s != 0, , drop = FALSE]
  message(sprintf(
    paste0(
      "rsem_apply_reference: %s genes x %d samples against %d reference samples; ",
      "standard deviation 0: %s genes (%s); not in the reference: %s genes removed; ",
      "reference genes missing from x: %s"
    ),
    .fmt(nrow(z)), ncol(z), length(reference$samples), .fmt(sum(s == 0)),
    if (zero_sd == "na") "set to NA" else "removed", .fmt(dropped), .fmt(missing_ref)
  ))
  res <- .as_table(z, from = x)
  attr(res, "reference") <- reference
  res
}

#' @rdname rsem_fit_reference
#' @export
rsem_reference_scale <- function(x, reference, ..., zero_sd = "na") {
  rsem_apply_reference(x, rsem_fit_reference(x, reference, ...), zero_sd = zero_sd)
}

#' @export
print.rsem_reference <- function(x, ...) {
  cat(sprintf(
    "rsem_reference: %s genes, %d reference samples (%s), %s, %s genes with standard deviation 0\n",
    .fmt(length(x$mean)), length(x$samples), paste(x$samples, collapse = ", "),
    if (x$transform == "log2") sprintf("log2(value + %s)", format(x$pseudocount)) else "no transformation",
    .fmt(sum(x$sd == 0))
  ))
  invisible(x)
}

.ref_transform <- function(m, transform, pseudocount) {
  if (transform == "none") return(m)
  if (pseudocount == 0 && any(m == 0)) .stopf("log2 with pseudocount = 0 would turn zeros into -Inf; use a positive pseudocount.")
  log2(m + pseudocount)
}
