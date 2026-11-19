#' How much each gene's effective length varies between samples
#'
#' Summarizes the gene-level effective lengths from
#' `rsem_read(expr = "length")`. RSEM's gene effective length is an average of
#' its transcripts' effective lengths weighted by their abundance, so it
#' changes between samples when transcript usage or the fragment length
#' distribution changes, and it is 0 where every transcript is shorter than
#' the mean fragment length. Large variation is therefore not by itself a sign
#' of an error or of poor quality.
#'
#' Genes with length 0 in some samples are the ones that stop
#' `tximport`/`DESeq2::DESeqDataSetFromTximport()`; `zero_fraction > 0` finds
#' them.
#'
#' @param x Effective lengths, as returned by `rsem_read(expr = "length")`.
#'
#' @return A data frame with one row per gene: `gene`, `median`, `mean`,
#'   `cv` (standard deviation / mean, over all samples), `zero_fraction`
#'   (share of samples with length 0) and `max_min_ratio` (largest / smallest
#'   nonzero length; `NA` with fewer than two nonzero samples).
#'
#' @examples
#' len <- data.frame(S1 = c(1000, 0, 500), S2 = c(1100, 200, 500),
#'                   row.names = c("A", "B", "C"))
#' rsem_length_profile(len)
#' @export
rsem_length_profile <- function(x) {
  m <- .expr_matrix(x)
  mu <- rowMeans(m)
  s <- if (ncol(m) > 1L) sqrt(rowSums((m - mu)^2) / (ncol(m) - 1)) else rep(NA_real_, nrow(m))
  nz <- m
  nz[nz == 0] <- NA_real_
  n_nz <- rowSums(!is.na(nz))
  ratio <- suppressWarnings(apply(nz, 1L, max, na.rm = TRUE) / apply(nz, 1L, min, na.rm = TRUE))
  ratio[n_nz < 2L] <- NA_real_
  out <- data.frame(
    gene = rownames(m), median = apply(m, 1L, stats::median), mean = mu,
    cv = ifelse(mu > 0, s / mu, NA_real_), zero_fraction = rowMeans(m == 0),
    max_min_ratio = unname(ratio), stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  message(sprintf(
    "rsem_length_profile: %s genes x %d samples; length 0 in some samples: %s genes, in every sample: %s genes",
    .fmt(nrow(m)), ncol(m), .fmt(sum(out$zero_fraction > 0 & out$zero_fraction < 1)), .fmt(sum(out$zero_fraction == 1))
  ))
  out
}
