#' Turn the table into long format with the sample information
#'
#' Returns one row per gene and sample, with the columns `gene`, `sample_id`
#' and `value`, followed by every column of the sample table. This is the
#' layout that ggplot2 and dplyr work with: a group or a continuous variable
#' in the sample table can go straight to colour, facets or the x axis.
#'
#' @section Sample table:
#' The sample table is the same `metadata.csv` that [rsem_read()] uses (or a
#' data frame). It needs a `sample_id` column; every other column is added as
#' it is, with no setting per column: columns of numbers stay numeric (for
#' example age or dose) and other columns stay text (for example group or
#' batch). `sample_id` is always read as text, so `001` stays `001`.
#'
#' Every sample in `x` must be in the sample table. Rows of the sample table
#' whose sample is not in `x` (for example samples removed earlier) are left
#' out, and a message says so.
#'
#' `sample_id` is a factor in the column order of `x`, so plots show the
#' samples in that order. Each gene's rows are next to each other.
#'
#' @param x The gene by sample table, as returned by [rsem_read()] or the
#'   other functions of this package.
#' @param metadata The sample table: a path to a CSV file, or a data frame,
#'   with a `sample_id` column. Optional; without it only `gene`, `sample_id`
#'   and `value` are returned.
#'
#' @return A data frame with `nrow(x) * ncol(x)` rows.
#'
#' @examples
#' if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   ex <- system.file("extdata", "example", package = "rsemverse")
#'   meta <- file.path(ex, "metadata.csv")
#'   tpm <- rsem_read(file.path(ex, "data"), meta, species = "human", expr = "tpm")
#'   long <- rsem_long(rsem_transform(tpm, "log2"), metadata = meta)
#'   head(long)
#'   # ggplot2::ggplot(subset(long, gene == "TP53"),
#'   #                 ggplot2::aes(sample_id, value)) + ggplot2::geom_col()
#' }
#' @export
rsem_long <- function(x, metadata = NULL) {
  m <- .expr_matrix(x, nonnegative = FALSE)
  samples <- colnames(m)
  out <- data.frame(
    gene = rep(rownames(m), each = ncol(m)),
    sample_id = factor(rep(samples, times = nrow(m)), levels = samples),
    value = as.vector(t(m)),
    stringsAsFactors = FALSE
  )
  if (is.null(metadata)) {
    message(sprintf("rsem_long: %s rows (%s genes x %d samples)", .fmt(nrow(out)), .fmt(nrow(m)), ncol(m)))
    return(out)
  }

  meta <- .read_sample_table(metadata)
  absent <- setdiff(samples, meta$sample_id)
  if (length(absent)) {
    .stopf("These samples are not in the sample table: %s", paste(utils::head(absent, 10L), collapse = ", "))
  }
  extra <- setdiff(meta$sample_id, samples)
  meta <- meta[match(samples, meta$sample_id), , drop = FALSE]
  clash <- intersect(setdiff(names(meta), "sample_id"), c("gene", "value"))
  if (length(clash)) {
    .stopf("The sample table has column(s) named %s, which are used for the output; rename them.", paste(clash, collapse = ", "))
  }
  info <- meta[as.integer(out$sample_id), setdiff(names(meta), "sample_id"), drop = FALSE]
  rownames(info) <- NULL
  out <- cbind(out, info)

  message(sprintf(
    "rsem_long: %s rows (%s genes x %d samples) with %d sample column(s): %s%s",
    .fmt(nrow(out)), .fmt(nrow(m)), ncol(m), ncol(info),
    if (ncol(info)) paste(sprintf("%s (%s)", names(info), vapply(info, .kind, character(1))), collapse = ", ") else "none",
    if (length(extra)) sprintf("; %d sample(s) in the sample table but not in the data were left out", length(extra)) else ""
  ))
  out
}

.kind <- function(v) {
  if (is.numeric(v)) "number" else if (is.logical(v)) "TRUE/FALSE" else "text"
}

.read_sample_table <- function(metadata) {
  if (is.character(metadata) && length(metadata) == 1L && !is.na(metadata)) {
    if (!file.exists(metadata) || dir.exists(metadata)) .stopf("`metadata` must be an existing CSV file: %s", metadata)
    meta <- data.table::fread(
      metadata, sep = ",", header = TRUE, encoding = "UTF-8", data.table = FALSE,
      check.names = FALSE, showProgress = FALSE, colClasses = list(character = "sample_id"),
      strip.white = TRUE
    )
  } else if (is.data.frame(metadata)) {
    meta <- as.data.frame(metadata, stringsAsFactors = FALSE)
  } else {
    .stopf("`metadata` must be a path to a CSV file or a data frame.")
  }
  if (!"sample_id" %in% names(meta)) .stopf("The sample table must have a 'sample_id' column.")
  meta$sample_id <- trimws(as.character(meta$sample_id))
  if (anyNA(meta$sample_id) || any(!nzchar(meta$sample_id))) .stopf("The sample table contains missing or empty sample_id values.")
  if (anyDuplicated(meta$sample_id)) {
    .stopf("The sample table contains duplicated sample_id values: %s", paste(unique(meta$sample_id[duplicated(meta$sample_id)]), collapse = ", "))
  }
  meta
}
