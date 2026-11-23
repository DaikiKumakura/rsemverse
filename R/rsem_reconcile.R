#' Compare gene-level and isoform-level results of the same samples
#'
#' RSEM defines a gene's abundance as the sum of the abundances of its
#' transcripts. `rsem_reconcile()` checks this for every sample, which
#' detects files that belong to different samples or runs, truncated files, and
#' gene files that were quantified against another reference than the isoform
#' files.
#'
#' RSEM writes `expected_count`, `TPM` and `FPKM` with two decimals, and it
#' rounds the gene and the transcript values separately. The sum of the
#' printed transcript values can therefore differ from the printed gene value
#' by at most `0.005 * (n + 1)` for a gene with `n` transcripts; a difference
#' above this is not explained by rounding. `tolerance_extra` is added to that
#' bound.
#'
#' @param data Folder with the `*.genes.results` files (searched recursively).
#' @param isoform_data Folder with the `*.isoforms.results` files; the default
#'   is `data`.
#' @param columns Value columns to compare. Any of `"expected_count"`, `"TPM"`
#'   and `"FPKM"`.
#' @param sample_from `"filename"` or `"parent_dir"`, as in [rsem_read()].
#' @param tolerance_extra Added to the rounding bound, to allow for files that
#'   were written with more digits (default `1e-6`).
#'
#' @return A data frame with one row per sample and column: `n_genes` compared,
#'   `n_beyond_rounding` (genes whose difference is above the bound),
#'   `max_abs_diff`, the `worst_gene`, and the numbers of genes only in the gene
#'   file (`genes_without_isoforms`) or only in the isoform file
#'   (`isoform_genes_missing`), and `transcript_mapping_mismatch`, the number of
#'   genes whose `transcript_id(s)` in the gene file are not the transcripts of
#'   that gene in the isoform file (a change of reference that the sums cannot
#'   show), and `transcript_mapping_checked`, `FALSE` when the gene file has no
#'   `transcript_id(s)` column so that this comparison could not be made. Samples found in only one of the folders are
#'   reported as a stop.
#' @examples
#' \dontrun{
#' rsem_reconcile("rsem")
#' }
#' @export
rsem_reconcile <- function(data = "data", isoform_data = data,
                           columns = c("expected_count", "TPM", "FPKM"),
                           sample_from = "filename", tolerance_extra = 1e-6) {
  if (!is.character(columns) || !length(columns) || !all(columns %in% c("expected_count", "TPM", "FPKM"))) {
    .stopf("`columns` must be one or more of \"expected_count\", \"TPM\", \"FPKM\".")
  }
  columns <- unique(columns)
  sample_from <- .choice(sample_from, c("filename", "parent_dir"), "sample_from")
  if (!is.numeric(tolerance_extra) || length(tolerance_extra) != 1L || !is.finite(tolerance_extra) || tolerance_extra < 0) {
    .stopf("`tolerance_extra` must be a single finite number that is not negative.")
  }
  g <- suppressMessages(.discover_samples(data, "gene", sample_from))
  i <- suppressMessages(.discover_samples(isoform_data, "isoform", sample_from))
  if (!setequal(g$sample_id, i$sample_id)) {
    .stopf(
      "Gene files and isoform files do not have the same samples. Only genes: %s. Only isoforms: %s.",
      toString(setdiff(g$sample_id, i$sample_id)), toString(setdiff(i$sample_id, g$sample_id))
    )
  }

  rows <- list()
  for (sid in g$sample_id) {
    gp <- g$path[g$sample_id == sid]
    ip <- i$path[i$sample_id == sid]
    for (pth in c(gp, ip)) {
      hdr <- .rsem_header(pth)
      need_id <- if (pth == gp) "gene_id" else c("gene_id", "transcript_id")
      if (!all(need_id %in% hdr)) .stopf("%s has no '%s' column.", pth, paste(setdiff(need_id, hdr), collapse = "', '"))
    }
    G <- .fread_rsem(gp, colClasses = list(character = "gene_id"))
    I <- .fread_rsem(ip, colClasses = list(character = c("transcript_id", "gene_id")))
    need <- setdiff(columns, intersect(names(G), names(I)))
    if (length(need)) .stopf("Column(s) %s are missing in the files of sample '%s'.", toString(need), sid)
    if (anyNA(G$gene_id) || anyNA(I$gene_id) || anyNA(I$transcript_id) ||
        any(!nzchar(G$gene_id)) || any(!nzchar(I$gene_id)) || any(!nzchar(I$transcript_id))) {
      .stopf("Sample '%s': empty or missing IDs; run rsem_check() first.", sid)
    }
    if (anyDuplicated(G$gene_id)) .stopf("Sample '%s': the gene file has duplicated gene_id values.", sid)
    if (anyDuplicated(I$transcript_id)) .stopf("Sample '%s': the isoform file has duplicated transcript_id values.", sid)
    for (col in columns) {
      if (!is.numeric(G[[col]]) || !is.numeric(I[[col]]) || any(!is.finite(G[[col]])) || any(!is.finite(I[[col]]))) {
        .stopf("Sample '%s': column '%s' has missing, infinite or non-numeric values; run rsem_check() first.", sid, col)
      }
      if (any(G[[col]] < 0) || any(I[[col]] < 0)) {
        .stopf("Sample '%s': column '%s' has negative values (RSEM writes none); run rsem_check() first.", sid, col)
      }
    }
    n_iso <- table(I$gene_id)
    only_gene <- setdiff(G$gene_id, I$gene_id)
    only_iso <- setdiff(I$gene_id, G$gene_id)
    common <- intersect(G$gene_id, I$gene_id)
    # The transcripts that the gene file lists for a gene must be the transcripts
    # of that gene in the isoform file, whatever the values.
    tx_mismatch <- 0L
    tx_checked <- "transcript_id(s)" %in% names(G)
    if (tx_checked) {
      listed <- strsplit(G[["transcript_id(s)"]][match(common, G$gene_id)], ",", fixed = TRUE)
      in_iso <- split(I$transcript_id, factor(I$gene_id, levels = common))
      tx_mismatch <- sum(!mapply(function(a, b) setequal(a, b), listed, in_iso[common]))
    }
    for (col in columns) {
      s <- rowsum(I[[col]], I$gene_id, reorder = FALSE)
      gv <- G[[col]][match(common, G$gene_id)]
      iv <- s[common, 1L]
      d <- abs(gv - iv)
      bound <- 0.005 * (as.vector(n_iso[common]) + 1) + tolerance_extra
      beyond <- d > bound
      rows[[length(rows) + 1L]] <- data.frame(
        sample_id = sid, column = col, n_genes = length(common),
        n_beyond_rounding = sum(beyond),
        max_abs_diff = if (length(d)) max(d) else NA_real_,
        worst_gene = if (length(d)) common[[which.max(d - bound)]] else NA_character_,
        genes_without_isoforms = length(only_gene), isoform_genes_missing = length(only_iso),
        transcript_mapping_mismatch = tx_mismatch, transcript_mapping_checked = tx_checked,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  bad <- out$n_beyond_rounding > 0 | out$genes_without_isoforms > 0 | out$isoform_genes_missing > 0 |
    out$transcript_mapping_mismatch > 0
  message(sprintf(
    "rsem_reconcile: %d samples, %s; %d of %d sample-column checks are not explained by rounding",
    nrow(g), paste(columns, collapse = ", "), sum(bad), nrow(out)
  ))
  if (!all(out$transcript_mapping_checked)) {
    message("  the transcripts of each gene were not compared for ", sum(!out$transcript_mapping_checked) / length(columns),
            " sample(s): the gene file has no 'transcript_id(s)' column")
  }
  out
}
