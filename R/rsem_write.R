#' Write the table with gene IDs and annotation
#'
#' Writes the gene by sample table as a tab-separated file, with the gene in
#' the first column, optionally followed by the original gene IDs behind each
#' row and the gene type from a GTF file, then one column per sample. A file
#' name ending in `.gz` is compressed. The file is for sharing and for other
#' tools; to continue in R, keep using the table itself.
#'
#' @param x The gene by sample table. Missing values, such as the `NA` that
#'   [rsem_apply_reference()] gives a gene with no spread in the reference, are
#'   written as `NA`.
#' @param file Output path (`.tsv` or `.tsv.gz`, or `.gct` with `format = "gct"`).
#' @param format `"tsv"` (the default) or `"gct"`: the GCT 1.2 format that GSEA,
#'   ssGSEA and GenePattern read (the lines `#1.2` and the numbers of genes and
#'   samples, then `Name` (the gene), `Description` (`NA`) and one column per
#'   sample). Values are written as they are: GSEA expects normalized values
#'   (not raw counts), and ssGSEA gene-level TPM. The gene IDs and the gene
#'   type are not written in the GCT file.
#' @param ids The attribute `"gene_ids"` of the table [rsem_read()] returned,
#'   `attr(tpm, "gene_ids")`: adds the column `gene_id` with the original IDs
#'   of each row, joined by `";"` where several IDs were summed into one
#'   symbol. The other functions return new tables without this attribute, so
#'   take it from the table [rsem_read()] gave. By default the attribute of
#'   `x` is used when it is there.
#' @param annotation Optional GTF file or data frame, as in
#'   [rsem_biotype_profile()]: adds the column `gene_type`.
#'
#' @return `file`, invisibly.
#'
#' @examples
#' x <- data.frame(S1 = c(1.5, 2), S2 = c(3, 4), row.names = c("TP53", "GAPDH"))
#' f <- tempfile(fileext = ".tsv")
#' rsem_write(x, f)
#' readLines(f)
#' @export
rsem_write <- function(x, file, ids = attr(x, "gene_ids", exact = TRUE), annotation = NULL, format = "tsv") {
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) .stopf("`file` must be a file path.")
  format <- .choice(format, c("tsv", "gct"), "format")
  m <- .expr_matrix(x, nonnegative = FALSE, allow_na = TRUE)
  if (format == "gct") {
    if (grepl("\\.gz$", file)) .stopf("A GCT file is written as plain text; use a file name ending in .gct.")
    con <- file(file, "w")
    on.exit(close(con))
    writeLines(c("#1.2", paste(nrow(m), ncol(m), sep = "	"), paste(c("Name", "Description", colnames(m)), collapse = "	")), con)
    body <- data.frame(Name = rownames(m), Description = "NA", as.data.frame(m, check.names = FALSE), check.names = FALSE, stringsAsFactors = FALSE)
    utils::write.table(body, con, sep = "	", quote = FALSE, row.names = FALSE, col.names = FALSE, na = "NA")
    message(sprintf("rsem_write: %s genes x %d samples to %s (GCT 1.2)", .fmt(nrow(m)), ncol(m), file))
    return(invisible(file))
  }
  out <- data.frame(gene = rownames(m), stringsAsFactors = FALSE)
  if (!is.null(ids)) {
    if (!is.data.frame(ids) || !all(c("gene", "gene_id") %in% names(ids))) {
      .stopf("`ids` must be the attribute \"gene_ids\" of a table from rsem_read() (columns gene and gene_id).")
    }
    joined <- tapply(ids$gene_id, ids$gene, paste, collapse = ";")
    out$gene_id <- unname(joined[out$gene])
    if (anyNA(out$gene_id)) .stopf("`ids` does not cover %s row(s) of `x`; it belongs to another table.", .fmt(sum(is.na(out$gene_id))))
  }
  if (!is.null(annotation)) out$gene_type <- .annotate_rows(out$gene, .read_annotation(annotation), ids)
  out <- cbind(out, as.data.frame(m, check.names = FALSE), stringsAsFactors = FALSE)
  data.table::fwrite(out, file, sep = "\t", quote = FALSE, na = "NA", compress = if (grepl("\\.gz$", file)) "gzip" else "none")
  message(sprintf(
    "rsem_write: %s genes x %d samples to %s (columns: %s)",
    .fmt(nrow(m)), ncol(m), file, paste(c(setdiff(names(out), colnames(m)), "samples"), collapse = ", ")
  ))
  invisible(file)
}
