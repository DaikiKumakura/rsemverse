#' Read RSEM gene-level results as a gene symbol by sample table
#'
#' Reads the `*.genes.results` files listed in a metadata CSV and returns one
#' expression value per gene symbol and sample, either the RSEM expected count
#' or the TPM.
#'
#' @section Input:
#' `metadata` is a CSV with the columns `sample_id` and `ref`; other columns
#' are ignored. `ref` is the name of the sample's RSEM file inside `data`
#' (that is, the file read is `file.path(data, ref)`).
#'
#' ```text
#' sample_id,ref
#' S1,S1.genes.results
#' S2,S2.genes.results
#' ```
#'
#' Every file needs the columns `gene_id` and `expected_count` or `TPM`
#' (whichever `expr` selects); other columns are ignored. All files must
#' contain the same genes in the same order, without duplicates.
#'
#' @section Gene symbols:
#' Symbols come from `org.Hs.eg.db` (human) or `org.Mm.eg.db` (mouse). The
#' version suffix of an Ensembl gene ID (`ENSG00000141510.18`) is removed
#' before the lookup, because the OrgDb keys have no version. The rules are:
#'
#' * An Ensembl ID with several symbols gets the symbol of its lowest Entrez
#'   ID.
#' * Rows that end up with the same symbol are **summed** (counts and TPM are
#'   both additive).
#' * A gene without a symbol (not in the OrgDb, or not an Ensembl ID, such as
#'   a spike-in or an ID ending in `_PAR_Y`) keeps its `gene_id` from the
#'   file, so no expression is lost. If such an ID equals a gene symbol, the
#'   function stops instead of merging the two.
#'
#' Rows appear in the order in which they first occur in the files. Values are
#' not rounded and TPM is not renormalized.
#'
#' @param data Directory that contains the RSEM files.
#' @param metadata Path to the metadata CSV.
#' @param species `"human"` or `"mouse"`. Required. The function stops if the
#'   gene IDs are of the other species.
#' @param expr `"count"` (`expected_count`) or `"tpm"` (`TPM`). Required.
#'
#' @return A data frame with one row per gene symbol (as row names) and one
#'   column per sample, in the order of the metadata. A message reports how
#'   many genes received a symbol.
#'
#' @examples
#' if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   ex <- system.file("extdata", "example", package = "rsemverse")
#'   tpm <- rsem_read(
#'     data = file.path(ex, "data"),
#'     metadata = file.path(ex, "metadata.csv"),
#'     species = "human",
#'     expr = "tpm"
#'   )
#'   head(tpm)
#' }
#' @export
rsem_read <- function(data = "data", metadata = "metadata.csv", species, expr) {
  if (missing(species)) .stopf("`species` must be \"human\" or \"mouse\".")
  if (missing(expr)) .stopf("`expr` must be \"count\" or \"tpm\".")
  species <- .choice(species, c("human", "mouse"), "species")
  expr <- .choice(expr, c("count", "tpm"), "expr")

  samples <- .read_samples(data, metadata)
  mat <- .read_matrix(samples, c(count = "expected_count", tpm = "TPM")[[expr]])

  stable <- .strip_ensembl_version(rownames(mat))
  .check_species(stable, species)
  orgdb <- .orgdb_package(species)
  symbol <- .lookup_symbols(stable, orgdb)

  out <- .sum_by_symbol(mat, symbol)

  has_symbol <- !is.na(symbol)
  message(sprintf(
    paste0(
      "rsem_read: %d samples, %d genes -> %d rows (species = %s, expr = %s, %s %s)\n",
      "  with symbol: %d genes (%d symbols merged from several IDs); ",
      "kept as original IDs: %d genes"
    ),
    ncol(mat), nrow(mat), nrow(out), species, expr, orgdb,
    as.character(utils::packageVersion(orgdb)),
    sum(has_symbol), sum(table(symbol[has_symbol]) > 1L), sum(!has_symbol)
  ))

  as.data.frame(out, check.names = FALSE)
}

.stopf <- function(fmt, ...) {
  stop(sprintf(fmt, ...), call. = FALSE)
}

# Bioconductor functions print progress messages; warnings and errors still
# propagate.
.quietly <- function(expr) {
  suppressMessages(expr)
}

.choice <- function(x, choices, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !x %in% choices) {
    .stopf("`%s` must be one of: %s.", name, paste(sprintf("\"%s\"", choices), collapse = ", "))
  }
  x
}

# Only Ensembl gene IDs lose their version suffix; other IDs (spike-ins,
# numeric IDs, `..._PAR_Y`) are returned unchanged.
.strip_ensembl_version <- function(x) {
  sub("^(ENS[A-Z]*G[0-9]+)\\.[0-9]+$", "\\1", as.character(x))
}

.read_samples <- function(data, metadata) {
  if (!is.character(data) || length(data) != 1L || is.na(data) || !dir.exists(data)) {
    .stopf("`data` must be an existing directory: %s", toString(data))
  }
  if (!is.character(metadata) || length(metadata) != 1L || is.na(metadata) ||
      !file.exists(metadata) || dir.exists(metadata)) {
    .stopf("`metadata` must be an existing CSV file: %s", toString(metadata))
  }

  tab <- data.table::fread(
    metadata,
    sep = ",",
    header = TRUE,
    colClasses = "character",
    encoding = "UTF-8",
    data.table = FALSE,
    check.names = FALSE,
    showProgress = FALSE
  )
  missing_cols <- setdiff(c("sample_id", "ref"), names(tab))
  if (length(missing_cols)) {
    .stopf(
      "The metadata must have the columns 'sample_id' and 'ref'; missing: %s",
      paste(missing_cols, collapse = ", ")
    )
  }
  sample_id <- trimws(tab$sample_id)
  ref <- trimws(tab$ref)
  if (!length(sample_id)) .stopf("The metadata has no rows.")
  if (anyNA(sample_id) || any(!nzchar(sample_id))) {
    .stopf("The metadata contains missing or empty sample_id values.")
  }
  if (anyNA(ref) || any(!nzchar(ref))) {
    .stopf("The metadata contains missing or empty ref values.")
  }
  if (anyDuplicated(sample_id)) {
    .stopf(
      "The metadata contains duplicated sample_id values: %s",
      paste(unique(sample_id[duplicated(sample_id)]), collapse = ", ")
    )
  }
  if (anyDuplicated(ref)) {
    .stopf(
      "The same file is listed for more than one sample: %s",
      paste(unique(ref[duplicated(ref)]), collapse = ", ")
    )
  }

  path <- file.path(data, ref)
  missing_file <- !file.exists(path) | dir.exists(path)
  if (any(missing_file)) {
    shown <- utils::head(which(missing_file), 5L)
    .stopf(
      "%d file(s) not found in '%s': %s%s",
      sum(missing_file), data,
      paste(sprintf("%s -> %s", sample_id[shown], ref[shown]), collapse = "; "),
      if (sum(missing_file) > length(shown)) "; ..." else ""
    )
  }
  data.frame(sample_id = sample_id, ref = ref, path = path, stringsAsFactors = FALSE)
}

.read_matrix <- function(samples, column) {
  required <- c("gene_id", column)
  gene_id <- NULL
  mat <- NULL
  for (i in seq_len(nrow(samples))) {
    label <- sprintf("%s (sample %s)", samples$ref[[i]], samples$sample_id[[i]])
    header <- names(data.table::fread(
      samples$path[[i]], sep = "\t", nrows = 0L, check.names = FALSE, showProgress = FALSE
    ))
    missing_cols <- setdiff(required, header)
    if (length(missing_cols)) {
      .stopf("%s is missing the column(s): %s", label, paste(missing_cols, collapse = ", "))
    }

    x <- suppressWarnings(data.table::fread(
      samples$path[[i]],
      sep = "\t",
      select = required,
      colClasses = list(character = "gene_id", numeric = column),
      data.table = FALSE,
      showProgress = FALSE
    ))
    if (!is.numeric(x[[column]])) .stopf("%s has non-numeric values in '%s'.", label, column)
    if (anyNA(x[[column]]) || anyNA(x$gene_id)) .stopf("%s contains missing values.", label)
    if (anyDuplicated(x$gene_id)) .stopf("%s contains duplicated gene_id values.", label)

    if (is.null(gene_id)) {
      gene_id <- x$gene_id
      mat <- matrix(
        NA_real_, length(gene_id), nrow(samples),
        dimnames = list(gene_id, samples$sample_id)
      )
    } else if (!identical(x$gene_id, gene_id)) {
      .stopf(
        "%s does not have the same genes in the same order as %s.",
        label, samples$ref[[1]]
      )
    }
    mat[, i] <- x[[column]]
  }
  mat
}

.check_species <- function(stable, species) {
  share <- c(
    human = mean(grepl("^ENSG[0-9]+$", stable)),
    mouse = mean(grepl("^ENSMUSG[0-9]+$", stable))
  )
  if (max(share) < 0.5) {
    .stopf("The gene_id values do not look like Ensembl gene IDs (ENSG... or ENSMUSG...).")
  }
  detected <- names(which.max(share))
  if (!identical(detected, species)) {
    .stopf("species = \"%s\", but the gene IDs look like %s IDs.", species, detected)
  }
  invisible(TRUE)
}

.orgdb_package <- function(species) {
  c(human = "org.Hs.eg.db", mouse = "org.Mm.eg.db")[[species]]
}

# Gene symbol for each stable Ensembl ID, NA where the OrgDb has none.
.lookup_symbols <- function(stable, orgdb) {
  if (!.quietly(requireNamespace(orgdb, quietly = TRUE))) {
    .stopf("The package '%s' is required. Install it with BiocManager::install(\"%s\").", orgdb, orgdb)
  }
  db <- .quietly(getExportedValue(orgdb, orgdb))
  # select() fails when none of the keys exist, so keep only keys the OrgDb has.
  keys <- intersect(
    unique(stable[grepl("^ENS[A-Z]*G[0-9]+$", stable)]),
    .quietly(AnnotationDbi::keys(db, keytype = "ENSEMBL"))
  )
  if (!length(keys)) return(rep(NA_character_, length(stable)))
  ann <- .quietly(AnnotationDbi::select(
    db,
    keys = keys,
    keytype = "ENSEMBL",
    columns = c("SYMBOL", "ENTREZID")
  ))
  ann <- ann[!is.na(ann$SYMBOL) & nzchar(ann$SYMBOL), , drop = FALSE]
  ann <- ann[order(ann$ENSEMBL, suppressWarnings(as.numeric(ann$ENTREZID))), , drop = FALSE]
  ann <- ann[!duplicated(ann$ENSEMBL), , drop = FALSE]
  ann$SYMBOL[match(stable, ann$ENSEMBL)]
}

.sum_by_symbol <- function(mat, symbol) {
  unmapped_id <- rownames(mat)[is.na(symbol)]
  clash <- intersect(unmapped_id, symbol[!is.na(symbol)])
  if (length(clash)) {
    .stopf(
      "A gene_id without a symbol equals a gene symbol and would be merged with it: %s",
      paste(utils::head(clash, 5L), collapse = ", ")
    )
  }
  key <- ifelse(is.na(symbol), rownames(mat), symbol)
  out <- rowsum(mat, group = key, reorder = FALSE)
  if (!isTRUE(all.equal(colSums(out), colSums(mat), tolerance = 1e-9))) {
    .stopf("Internal error: the column sums changed while merging rows.")
  }
  out
}

#' rsemverse: read RSEM results as gene symbol tables
#'
#' [rsem_read()] reads RSEM gene results into a gene symbol by sample table.
#' The other functions take that table and return a table of the same layout
#' (or a long table), so they can be chained:
#'
#' * [rsem_filter()]: remove lowly expressed or uninformative genes.
#' * [rsem_normalize()]: normalize between samples (CPM, TMM, upper
#'   quartile, DESeq2 median of ratios, quantile).
#' * [rsem_transform()]: log2, z-score, centering, ranks.
#' * [rsem_subset()]: keep or remove built-in gene groups (protein-coding,
#'   mitochondrial, ribosomal, hemoglobin, sex, ERCC).
#' * [rsem_long()]: one row per gene and sample, with the sample table.
#'
#' @keywords internal
"_PACKAGE"
