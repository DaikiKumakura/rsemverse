#' Read several values of the same RSEM files at once
#'
#' Reads more than one column of the RSEM files (for example the expected
#' count, the TPM and the effective length) into tables with the same genes and
#' samples. The values that differential expression tools need together are in
#' one object, and the files are checked once for the same genes in the same
#' order, as in [rsem_read()], which does the reading for each value.
#'
#' @inheritParams rsem_read
#' @param assays One or more of the `expr` values of [rsem_read()]; the default
#'   is the expected count, the TPM and the effective length.
#' @param id `"original"` (the default, the gene IDs of the files),
#'   `"symbol"` or `"annotation"` as in [rsem_read()]. Unlike there, the
#'   default is the original IDs, because the next step is usually a tool that
#'   needs one row per RSEM gene.
#' @return A list of class `rsem_assays` with one table per assay, named as in
#'   `assays`. Every table has the same row and column names, and the
#'   attributes of [rsem_read()] (`"gene_ids"`, `"zero_length"`, `"tx2gene"`) are
#'   carried on the list.
#' @examples
#' \dontrun{
#' a <- rsem_import("rsem", assays = c("count", "tpm", "length"))
#' txi <- rsem_as_tximport(a)
#' # dds <- DESeq2::DESeqDataSetFromTximport(txi, colData = samples, design = ~ group)
#' }
#' @export
rsem_import <- function(data = "data", metadata = "metadata.csv", assays = c("count", "tpm", "length"),
                        species = NULL, level = "gene", id = "original", sample_from = "filename",
                        validate_reference = TRUE, allow_reorder = FALSE, annotation = NULL) {
  if (!is.character(assays) || !length(assays) || anyNA(assays) || anyDuplicated(assays)) {
    .stopf("`assays` must be one or more different values of `expr` in rsem_read().")
  }
  tabs <- lapply(assays, function(a) {
    suppressMessages(rsem_read(
      data = data, metadata = metadata, species = species, expr = a, level = level, id = id,
      sample_from = sample_from, validate_reference = validate_reference, allow_reorder = allow_reorder,
      annotation = annotation
    ))
  })
  names(tabs) <- assays
  ref_rows <- rownames(tabs[[1L]])
  ref_cols <- colnames(tabs[[1L]])
  same <- vapply(tabs, function(t) identical(rownames(t), ref_rows) && identical(colnames(t), ref_cols), logical(1))
  if (!all(same)) .stopf("Internal error: the assays %s do not have the same rows and samples.", toString(names(tabs)[!same]))
  out <- structure(tabs, class = c("rsem_assays", "list"))
  for (a in c("gene_ids", "zero_length", "tx2gene")) {
    for (t in tabs) {
      v <- attr(t, a, exact = TRUE)
      if (!is.null(v)) { attr(out, a) <- v; break }
    }
  }
  message(sprintf(
    "rsem_import: %d samples, %d rows; assays: %s",
    length(ref_cols), length(ref_rows), paste(assays, collapse = ", ")
  ))
  out
}

#' @export
print.rsem_assays <- function(x, ...) {
  cat(sprintf("<rsem_assays> %d assays (%s), %d rows x %d samples\n", length(x), paste(names(x), collapse = ", "), nrow(x[[1L]]), ncol(x[[1L]])))
  invisible(x)
}

#' The list that tximport returns, from RSEM values read by rsem_import()
#'
#' Builds the list `list(abundance, counts, length, countsFromAbundance)` that
#' `DESeq2::DESeqDataSetFromTximport()` and `edgeR`'s offset recipe take, from
#' the expected counts, TPM and effective lengths of RSEM files, without
#' changing a value. It is the same list that `tximport::tximport(type =
#' "rsem", txIn = FALSE)` returns for the same files (`countsFromAbundance =
#' "no"`); the differential expression itself, and the length offset, stay
#' with those packages.
#'
#' DESeq2 stops when a gene has an effective length of 0, which RSEM writes
#' for a transcript shorter than the mean fragment length (see
#' [rsem_read()]). `zero_length` decides what happens: `"stop"` (the default)
#' names the genes, `"drop"` removes every gene that has a zero length in any
#' sample (the counts of these genes are lost; the number is printed), and
#' `"keep"` leaves them. Replacing 0 by another length is not offered: it
#' changes the offset of those genes only.
#'
#' Rows that combine several gene IDs (from `id = "symbol"`) are refused:
#' their length is an average over different genes, not a length that RSEM
#' estimated.
#'
#' @param x The result of [rsem_import()] with the assays `"count"`, `"tpm"`
#'   and `"length"`.
#' @param zero_length `"stop"`, `"drop"` or `"keep"`, see Details.
#' @return A list with the matrices `abundance` (TPM), `counts` and `length`,
#'   and `countsFromAbundance = "no"`.
#' @export
rsem_as_tximport <- function(x, zero_length = "stop") {
  if (!inherits(x, "rsem_assays")) .stopf("`x` must be the result of rsem_import().")
  zero_length <- .choice(zero_length, c("stop", "drop", "keep"), "zero_length")
  need <- c("count", "tpm", "length")
  if (!all(need %in% names(x))) .stopf("`x` needs the assays %s; it has %s.", toString(need), toString(names(x)))
  ids <- attr(x, "gene_ids", exact = TRUE)
  if (is.data.frame(ids) && anyDuplicated(ids$gene[ids$gene %in% rownames(x$count)])) {
    .stopf("Some rows combine several gene IDs (rsem_id_report() lists them); their length is not an RSEM length. Use rsem_import(id = \"original\").")
  }
  counts <- as.matrix(x$count)
  abundance <- as.matrix(x$tpm)
  len <- as.matrix(x$length)
  zero_rows <- which(rowSums(len == 0) > 0L)
  if (length(zero_rows)) {
    if (zero_length == "stop") {
      .stopf(
        "%d genes have an effective length of 0 in at least one sample (e.g. %s), and DESeqDataSetFromTximport() stops on them. Use zero_length = \"drop\" to remove them or \"keep\" to leave them.",
        length(zero_rows), paste(utils::head(rownames(len)[zero_rows], 3L), collapse = ", ")
      )
    }
    if (zero_length == "drop") {
      keep <- setdiff(seq_len(nrow(len)), zero_rows)
      message(sprintf(
        "rsem_as_tximport: %d genes with an effective length of 0 removed (their expected counts add up to %s of %s)",
        length(zero_rows), format(round(sum(counts[zero_rows, , drop = FALSE]), 1)), format(round(sum(counts), 1))
      ))
      counts <- counts[keep, , drop = FALSE]; abundance <- abundance[keep, , drop = FALSE]; len <- len[keep, , drop = FALSE]
    }
  }
  list(abundance = abundance, counts = counts, length = len, countsFromAbundance = "no")
}
