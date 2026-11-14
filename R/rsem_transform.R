#' Transform values for plots, clustering and scores
#'
#' Applies one or more transformations in the order given. The genes and
#' samples stay the same. A message reports what was done.
#'
#' @section Methods:
#' * `"log2"`: `log2(value + pseudocount)`, with `pseudocount = 1` by
#'   default. Compresses the range of TPM, CPM or normalized counts so that
#'   highly expressed genes do not dominate plots, distances and PCA. Needs
#'   values that are not negative.
#' * `"zscore"`: for each gene, subtracts its mean over samples and divides by
#'   its standard deviation, so every gene has mean 0 and standard deviation
#'   1. Used for heatmaps and for comparing genes on one scale. Apply it after
#'   `"log2"`. A gene with the same value in every sample gets 0. Needs at
#'   least two samples.
#' * `"center"`: subtracts each gene's mean over samples, without dividing by
#'   the standard deviation.
#' * `"rank"`: within each sample, the rank of each gene divided by the number
#'   of genes, from near 0 (lowest) to 1 (highest); tied values share their
#'   average rank. It does not depend on the scale of the values, so samples
#'   measured or normalized differently can be compared. Rank-based scores
#'   such as singscore start from this.
#'
#' For example, `rsem_transform(tpm, c("log2", "zscore"))` gives the usual
#' heatmap scale.
#'
#' @param x The gene by sample table, as returned by [rsem_read()],
#'   [rsem_filter()] or [rsem_normalize()].
#' @param method One or more of `"log2"`, `"zscore"`, `"center"` and `"rank"`,
#'   applied in this order. Required.
#' @param pseudocount For `"log2"`: the value added before the logarithm.
#'   Default 1.
#'
#' @return A data frame of transformed values in the same layout as `x`.
#'   The attribute `"transform"` lists the methods applied.
#'
#' @examples
#' if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   ex <- system.file("extdata", "example", package = "rsemverse")
#'   tpm <- rsem_read(file.path(ex, "data"), file.path(ex, "metadata.csv"),
#'                    species = "human", expr = "tpm")
#'   log_tpm <- rsem_transform(tpm, "log2")               # log2(TPM + 1)
#'   heat <- rsem_transform(tpm, c("log2", "zscore"))     # heatmap scale
#'   ranks <- rsem_transform(tpm, "rank")
#' }
#' @export
rsem_transform <- function(x, method, pseudocount = 1) {
  if (missing(method)) .stopf("`method` is required; see ?rsem_transform.")
  choices <- c("log2", "zscore", "center", "rank")
  if (!is.character(method) || !length(method) || anyNA(method) || !all(method %in% choices)) {
    .stopf("`method` must be one or more of: %s.", paste(sprintf("\"%s\"", choices), collapse = ", "))
  }
  m <- .expr_matrix(x, nonnegative = FALSE)
  dn <- dimnames(m)
  done <- c(attr(x, "transform", exact = TRUE))
  notes <- character(0)

  for (step in method) {
    if (step == "log2") {
      pseudocount <- .number(pseudocount, "pseudocount", lower = 0)
      if (any(m < 0)) .stopf("\"log2\" needs values that are not negative; it cannot follow \"zscore\" or \"center\".")
      if (pseudocount == 0 && any(m == 0)) .stopf("\"log2\" with pseudocount = 0 would turn zeros into -Inf; use a positive pseudocount.")
      if ("log2" %in% done) .stopf("The values are already log2-transformed.")
      m <- log2(m + pseudocount)
      notes <- c(notes, sprintf("log2(value + %s)", format(pseudocount)))
    } else if (step %in% c("zscore", "center")) {
      if (ncol(m) < 2L) .stopf("\"%s\" needs at least two samples.", step)
      mu <- rowMeans(m)
      m <- m - mu
      if (step == "zscore") {
        s <- sqrt(rowSums(m^2) / (ncol(m) - 1))
        constant <- s <= sqrt(.Machine$double.eps) * pmax(1, abs(mu))
        m[!constant, ] <- m[!constant, , drop = FALSE] / s[!constant]
        m[constant, ] <- 0
        notes <- c(notes, sprintf("z-score per gene (%s constant genes set to 0)", .fmt(sum(constant))))
      } else {
        notes <- c(notes, "centered per gene")
      }
    } else {
      m <- apply(m, 2L, rank) / nrow(m)
      if (!is.matrix(m)) m <- matrix(m, nrow = 1L)
      notes <- c(notes, "rank within sample / number of genes")
    }
    done <- c(done, step)
  }
  dimnames(m) <- dn

  message(sprintf("rsem_transform: %s genes x %d samples; %s", .fmt(nrow(m)), ncol(m), paste(notes, collapse = ", then ")))
  res <- .as_table(m)
  attr(res, "transform") <- done
  res
}
