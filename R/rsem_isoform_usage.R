#' How the isoforms of each gene share its expression
#'
#' For every gene with at least two transcripts, and every sample, describes the split of the gene's
#' TPM among its transcripts: the share of the dominant transcript, the number of transcripts above a
#' share, and the Shannon entropy of the shares, with `exp(entropy)` as the effective number of
#' isoforms. These are descriptions of one sample; they are not a test of a change of isoform usage
#' between conditions (see `edgeR::catchRSEM()`, DRIMSeq or DEXSeq for that).
#'
#' The share of a transcript is its TPM divided by the TPM of its gene (the sum over its
#' transcripts). When the gene has no TPM in a sample the shares do not exist and every number is `NA`, not 0.
#' `IsoPct` of RSEM is the same share in percent; a gene with one transcript has `IsoPct` 100 whether it
#' is expressed or not, and is left out here.
#'
#' @param x A table of transcript TPM made by `rsem_read(level = "isoform", expr = "tpm")`.
#' @param tx2gene The transcript-to-gene table (`transcript_id`, `gene_id`); by default the attribute
#'   that [rsem_read()] sets.
#' @param min_share A transcript counts as expressed when its share is at least this (default 0.05).
#' @return A data frame with one row per gene and sample: `gene_id`, `sample_id`, `n_isoforms`
#'   (in the reference), `gene_tpm`, `dominant`, `dominant_share`, `n_expressed` (shares of at least
#'   `min_share`), `entropy` (natural logarithm) and `effective_isoforms`.
#' @examples
#' x <- data.frame(S1 = c(30, 10, 0, 0), S2 = c(0, 0, 0, 5), row.names = c("t1", "t2", "t3", "t4"))
#' tx <- data.frame(transcript_id = rownames(x), gene_id = c("g1", "g1", "g2", "g2"))
#' rsem_isoform_usage(x, tx2gene = tx)
#' @export
rsem_isoform_usage <- function(x, tx2gene = attr(x, "tx2gene", exact = TRUE), min_share = 0.05) {
  if (!is.data.frame(tx2gene) || !all(c("transcript_id", "gene_id") %in% names(tx2gene))) {
    .stopf("`tx2gene` must have the columns transcript_id and gene_id; read the transcripts with rsem_read(level = \"isoform\").")
  }
  min_share <- .number(min_share, "min_share", lower = 0, upper = 1)
  m <- .expr_matrix(x)
  gene_of <- tx2gene$gene_id[match(rownames(m), tx2gene$transcript_id)]
  if (anyNA(gene_of)) .stopf("%d transcripts of `x` are not in `tx2gene` (e.g. %s).", sum(is.na(gene_of)), rownames(m)[which(is.na(gene_of))[[1L]]])
  size <- table(gene_of)
  multi <- names(size)[size >= 2L]
  rows <- list()
  for (g in multi) {
    idx <- which(gene_of == g)
    sub <- m[idx, , drop = FALSE]
    tot <- colSums(sub)
    for (j in seq_len(ncol(m))) {
      if (tot[[j]] > 0) {
        p <- sub[, j] / tot[[j]]
        h <- -sum(ifelse(p > 0, p * log(p), 0))
        rows[[length(rows) + 1L]] <- data.frame(
          gene_id = g, sample_id = colnames(m)[[j]], n_isoforms = length(idx), gene_tpm = tot[[j]],
          dominant = rownames(sub)[[which.max(p)]], dominant_share = max(p), n_expressed = sum(p >= min_share),
          entropy = h, effective_isoforms = exp(h), stringsAsFactors = FALSE
        )
      } else {
        rows[[length(rows) + 1L]] <- data.frame(
          gene_id = g, sample_id = colnames(m)[[j]], n_isoforms = length(idx), gene_tpm = 0, dominant = NA_character_,
          dominant_share = NA_real_, n_expressed = NA_integer_, entropy = NA_real_, effective_isoforms = NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame(
    gene_id = character(0), sample_id = character(0), n_isoforms = integer(0), gene_tpm = numeric(0), dominant = character(0),
    dominant_share = numeric(0), n_expressed = integer(0), entropy = numeric(0), effective_isoforms = numeric(0)
  )
  rownames(out) <- NULL
  message(sprintf(
    "rsem_isoform_usage: %d genes with several transcripts, %d samples; %d gene-sample pairs without TPM",
    length(multi), ncol(m), sum(is.na(out$entropy))
  ))
  out
}
