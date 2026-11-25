#' A SummarizedExperiment from RSEM values
#'
#' Wraps the values of [rsem_import()] (or one table) with the sample table and the gene IDs in a
#' `SummarizedExperiment`, the object that most Bioconductor packages take. Nothing is calculated or
#' changed: each assay is one of the tables, `colData` is the sample table, `rowData` has the original
#' gene IDs behind each row (several joined by `";"` where rows were merged), and the steps taken so far
#' (see [rsem_steps()]) are kept in `metadata(se)$steps`. For DESeq2 on RSEM gene results, the list of
#' [rsem_as_tximport()] is the more direct way.
#'
#' The `SummarizedExperiment` package is needed (it is not installed with rsemverse).
#'
#' @param x The result of [rsem_import()], or one table from [rsem_read()] (the assay is then named after
#'   its `expr`, or `"values"`).
#' @param metadata The sample table, a path to a CSV file or a data frame with a `sample_id` column; the
#'   rows are put in the order of the samples. Optional.
#' @return A `SummarizedExperiment`.
#' @examples
#' \dontrun{
#' se <- rsem_import("rsem", metadata = NULL) |> rsem_as_summarizedexperiment("metadata.csv")
#' }
#' @export
rsem_as_summarizedexperiment <- function(x, metadata = NULL) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    .stopf("rsem_as_summarizedexperiment() needs the SummarizedExperiment package: BiocManager::install(\"SummarizedExperiment\").")
  }
  tabs <- if (inherits(x, "rsem_assays")) {
    unclass(x)
  } else if (is.data.frame(x)) {
    stats::setNames(list(x), if (is.null(attr(x, "expr", exact = TRUE))) "values" else attr(x, "expr", exact = TRUE))
  } else {
    .stopf("`x` must be the result of rsem_import() or a table from rsem_read().")
  }
  first <- tabs[[1L]]
  assays <- lapply(tabs, function(t) as.matrix(t))
  col_data <- S4Vectors::DataFrame(row.names = colnames(first))
  if (!is.null(metadata)) {
    meta <- .read_sample_table(metadata)
    absent <- setdiff(colnames(first), meta$sample_id)
    if (length(absent)) .stopf("Sample(s) not in the sample table: %s", paste(utils::head(absent, 5L), collapse = ", "))
    meta <- meta[match(colnames(first), meta$sample_id), , drop = FALSE]
    rownames(meta) <- meta$sample_id
    col_data <- S4Vectors::DataFrame(meta, check.names = FALSE)
  }
  ids <- attr(x, "gene_ids", exact = TRUE)
  if (is.null(ids) && is.data.frame(x)) ids <- attr(x, "gene_ids", exact = TRUE)
  row_data <- S4Vectors::DataFrame(row.names = rownames(first))
  if (is.data.frame(ids) && all(c("gene", "gene_id") %in% names(ids))) {
    joined <- tapply(ids$gene_id, ids$gene, paste, collapse = ";")
    row_data$gene_id <- as.vector(joined[rownames(first)])
  }
  steps <- attr(x, "steps", exact = TRUE)
  SummarizedExperiment::SummarizedExperiment(
    assays = assays, colData = col_data, rowData = row_data,
    metadata = list(steps = if (is.null(steps)) character(0) else steps)
  )
}
