#' Select transcripts by expression
#'
#' Keeps the transcripts that reach a TPM threshold and, optionally, a share
#' of the most expressed transcript of the same gene, or only the most
#' expressed transcript of each gene. The rules follow Trinity's
#' `filter_low_expr_transcripts.pl` (`--min_expr_any`, `--min_pct_dom_iso`,
#' `--highest_iso_only`); as there, the share rule and the "most expressed
#' only" rule cannot be combined.
#'
#' @section How samples are combined:
#' In each sample, a transcript passes when its TPM is at least `min_tpm`
#' and, with `min_fraction`, its TPM is at least `min_fraction` times the
#' highest TPM among the transcripts of its gene in that sample (with
#' `highest_only`, when it is that highest transcript). The per-sample results
#' are then combined by `rule`:
#' * `"any"` (the default): passes in at least one sample;
#' * `"all"`: passes in every sample;
#' * `"prevalence"`: passes in at least the share `prevalence` of samples;
#' * `"mean"`: the rule is applied once to the mean TPM over samples.
#'
#' The most expressed transcript can differ between samples, so `"any"` with
#' `highest_only` may keep several transcripts of a gene, while `"mean"`
#' keeps exactly one per expressed gene.
#'
#' @section Caution:
#' Removing transcripts by expression can remove biologically important
#' ones, and keeping one transcript per gene is not a general preprocessing
#' step for differential expression. The values of the kept transcripts are
#' the original estimates: they were not re-quantified without the removed
#' transcripts. To use the selection as a new reference, build the reference
#' again and rerun RSEM.
#'
#' @param x Transcript TPM, as returned by `rsem_read(level = "isoform",
#'   expr = "tpm")`.
#' @param tx2gene Transcript-to-gene table with the columns `transcript_id`
#'   and `gene_id`. Default: the attribute `"tx2gene"` of `x`.
#' @param min_tpm TPM a transcript must reach. Default 1.
#' @param min_fraction Optional share (0 to 1) of the gene's highest TPM.
#' @param highest_only Keep only the most expressed transcript of each gene.
#' @param rule `"any"`, `"all"`, `"prevalence"` or `"mean"`.
#' @param prevalence For `rule = "prevalence"`: share of samples. Default 0.5.
#'
#' @return The kept transcripts, in the layout and order of `x`, with the
#'   matching rows of `tx2gene` in the attribute `"tx2gene"`.
#'
#' @examples
#' x <- data.frame(S1 = c(10, 1, 0.5), S2 = c(8, 3, 0), row.names = c("T1", "T2", "T3"))
#' t2g <- data.frame(transcript_id = c("T1", "T2", "T3"), gene_id = "G1")
#' rsem_select_isoforms(x, t2g, min_tpm = 1, min_fraction = 0.2)
#' rsem_select_isoforms(x, t2g, highest_only = TRUE, rule = "mean")
#' @export
rsem_select_isoforms <- function(x, tx2gene = attr(x, "tx2gene", exact = TRUE), min_tpm = 1, min_fraction = NULL,
                                 highest_only = FALSE, rule = "any", prevalence = 0.5) {
  m <- .expr_matrix(x)
  if (is.null(tx2gene) || !is.data.frame(tx2gene) || !all(c("transcript_id", "gene_id") %in% names(tx2gene))) {
    .stopf("`tx2gene` must be a data frame with the columns transcript_id and gene_id (as attr(x, \"tx2gene\") from rsem_read(level = \"isoform\")).")
  }
  rule <- .choice(rule, c("any", "all", "prevalence", "mean"), "rule")
  min_tpm <- .number(min_tpm, "min_tpm", lower = 0)
  if (!is.logical(highest_only) || length(highest_only) != 1L || is.na(highest_only)) .stopf("`highest_only` must be TRUE or FALSE.")
  if (!is.null(min_fraction)) {
    min_fraction <- .number(min_fraction, "min_fraction", lower = 0, upper = 1)
    if (highest_only) .stopf("`min_fraction` and `highest_only = TRUE` cannot be combined; choose one.")
  }
  if (rule == "prevalence") prevalence <- .number(prevalence, "prevalence", lower = 0, upper = 1)
  gene <- as.character(tx2gene$gene_id[match(rownames(m), tx2gene$transcript_id)])
  if (anyNA(gene)) .stopf("%s transcript(s) of `x` are not in `tx2gene`.", .fmt(sum(is.na(gene))))

  v <- if (rule == "mean") matrix(rowMeans(m), ncol = 1L, dimnames = list(rownames(m), "mean")) else m
  gmax <- apply(v, 2L, function(col) stats::ave(col, gene, FUN = max))
  if (!is.matrix(gmax)) gmax <- matrix(gmax, ncol = ncol(v))
  pass <- v >= min_tpm
  if (!is.null(min_fraction)) pass <- pass & v >= min_fraction * gmax
  if (highest_only) {
    if (rule == "mean") {
      # one transcript per gene: the first of tied top transcripts
      top <- which(v[, 1] == gmax[, 1])
      first <- matrix(FALSE, nrow(v), 1L)
      first[top[!duplicated(gene[top])], 1] <- TRUE
      pass <- pass & first
    } else {
      pass <- pass & v == gmax
    }
  }
  keep <- switch(rule,
    any = rowSums(pass) >= 1,
    all = rowSums(pass) == ncol(pass),
    prevalence = rowMeans(pass) >= prevalence - 1e-12,
    mean = pass[, 1]
  )
  message(sprintf(
    "rsem_select_isoforms: kept %s of %s transcripts (%s of %s genes); TPM >= %s%s%s, rule = %s. Values are the original estimates, not re-quantified.",
    .fmt(sum(keep)), .fmt(nrow(m)), .fmt(length(unique(gene[keep]))), .fmt(length(unique(gene))), format(min_tpm),
    if (!is.null(min_fraction)) sprintf(", >= %s of the gene's top transcript", format(min_fraction)) else "",
    if (highest_only) ", top transcript only" else "", rule
  ))
  res <- .as_table(m[keep, , drop = FALSE])
  attr(res, "tx2gene") <- tx2gene[match(rownames(res), tx2gene$transcript_id), , drop = FALSE]
  rownames(attr(res, "tx2gene")) <- NULL
  res
}
