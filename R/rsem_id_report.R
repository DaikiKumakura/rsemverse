#' Which gene IDs were merged into the same row
#'
#' `rsem_read()` sums the rows of gene IDs that have the same symbol (or the
#' same name from an annotation). `rsem_id_report()` lists these rows, so that
#' the merging can be inspected, for example before using summed counts or
#' averaged effective lengths downstream.
#'
#' The attribute `"gene_ids"` that [rsem_read()] sets, and the functions that
#' keep it, are the source. For `expr = "length"` the merged value is a
#' TPM-weighted mean of the lengths of different genes; it is not a length that
#' RSEM estimated for one gene, so it is not suitable as a length offset.
#'
#' @param x A table from [rsem_read()] or a function that keeps its `"gene_ids"`
#'   attribute.
#' @return A data frame with one row per merged row: `gene` (the row name),
#'   `n_ids`, and `gene_ids` (the IDs, separated by commas). It has no rows
#'   when nothing was merged.
#' @examples
#' \dontrun{
#' x <- rsem_read("rsem", species = "human")
#' rsem_id_report(x)
#' }
#' @export
rsem_id_report <- function(x) {
  ids <- attr(x, "gene_ids")
  if (!is.data.frame(ids) || !all(c("gene", "gene_id") %in% names(ids))) {
    .stopf("`x` has no \"gene_ids\" attribute; use a table from rsem_read().")
  }
  keep <- ids$gene %in% rownames(x)
  ids <- ids[keep, , drop = FALSE]
  n <- table(ids$gene)
  merged <- names(n)[n > 1L]
  out <- data.frame(
    gene = merged,
    n_ids = as.integer(n[merged]),
    gene_ids = vapply(merged, function(g) paste(ids$gene_id[ids$gene == g], collapse = ","), character(1)),
    stringsAsFactors = FALSE
  )
  out <- out[order(-out$n_ids, out$gene), , drop = FALSE]
  rownames(out) <- NULL
  message(sprintf("rsem_id_report: %d of %d rows combine several gene IDs", nrow(out), nrow(x)))
  out
}
