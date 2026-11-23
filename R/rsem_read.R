#' Read RSEM results as a gene (or transcript) by sample table
#'
#' Reads RSEM `*.genes.results` (or `*.isoforms.results`) files and returns
#' one value per gene symbol and sample: the expected count, the TPM, the
#' effective length or, for transcripts, the isoform percentage.
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
#' With `metadata = NULL` the files are found instead: every
#' `*.genes.results` (or `*.isoforms.results`, with or without `.gz`) file under `data`, searched
#' recursively. The sample ID is the file name without the suffix
#' (`sample_from = "filename"`) or the name of the folder that holds the file
#' (`sample_from = "parent_dir"`, for layouts such as `S1/rsem.genes.results`).
#' The files and the sample IDs found are printed, sorted by sample ID, and
#' duplicated sample IDs stop the function rather than overwrite each other.
#' Groups are never guessed from file names; give them in the sample table of
#' [rsem_long()] or [rsem_fit_reference()].
#'
#' Transcript files must give every transcript the same gene in all samples; a
#' file from a different reference stops the function.
#'
#' Every file needs the column `gene_id` (`transcript_id` for transcripts) and
#' the column `expr` selects; other columns are ignored. All files must
#' contain the same rows in the same order, without duplicates. A listed file
#' that is missing stops the function; no sample is dropped silently.
#'
#' @section Gene symbols:
#' With `id = "symbol"` (the default), symbols come from `org.Hs.eg.db`
#' (human) or `org.Mm.eg.db` (mouse). The version suffix of an Ensembl gene ID
#' (`ENSG00000141510.18`) is removed before the lookup, because the OrgDb
#' keys have no version. The rules are:
#'
#' * An Ensembl ID with several symbols gets the symbol of its lowest Entrez
#'   ID.
#' * Rows that end up with the same symbol are **summed** (counts and TPM are
#'   both additive), never reduced to their maximum. Effective lengths are not
#'   additive: they are averaged with the TPM of each row as weight, the way
#'   RSEM averages the lengths of a gene's transcripts (a plain mean where the
#'   rows have no TPM in that sample).
#' * A gene without a symbol (not in the OrgDb, or not an Ensembl ID, such as
#'   a spike-in or an ID ending in `_PAR_Y`) keeps its `gene_id` from the
#'   file, so no expression is lost. If such an ID equals a gene symbol, the
#'   function stops instead of merging the two.
#'
#' The original IDs behind every row are kept in the attribute `"gene_ids"`
#' (a data frame with the columns `gene` and `gene_id`); [rsem_write()] can
#' write them next to the values. Summing happens on the original scale,
#' before any log transformation: `log2(1 + TPM_a + TPM_b)` is not the sum of
#' `log2(1 + TPM_a)` and `log2(1 + TPM_b)`.
#'
#' With `id = "original"`, rows keep the `gene_id` of the file and no OrgDb or
#' `species` is needed, so results for any organism or a de novo transcriptome
#' can be read.
#'
#' Files made with `rsem-calculate-expression --append-names` have IDs such
#' as `ENSG00000141510.18_TP53`; the version and the appended name are removed
#' before the lookup, and rows without a symbol keep the full ID of the file.
#' Files may be gzip-compressed (`*.genes.results.gz`).
#'
#' Rows appear in the order in which they first occur in the files. Values are
#' not rounded and TPM is not renormalized.
#'
#' @section Names from an annotation file:
#' RefSeq-based references have IDs such as `gene1000`, which say nothing
#' about the gene; `rsem-calculate-expression --append-names` appends the name
#' (`gene1000_TP53`). With `id = "annotation"` and the GTF of the reference
#' (or a data frame), each row is named by the `gene_name` of its `gene_id`.
#' A file ID is matched to the annotation as it is, or without its Ensembl
#' version, or as `<gene_id>_<gene_name>` (the form `--append-names` writes):
#' the part before an underscore must be a `gene_id` of the annotation and the
#' part after it exactly that gene's `gene_name`. An ID such as `gene1_2` is
#' therefore not taken for `gene1`, and underscores inside IDs or names are
#' not cut at the wrong place. Rows that match nothing, or
#' have no name, keep their ID. Rows with the same name are summed, as for
#' symbols (this includes the `_PAR_Y` copies of pseudoautosomal genes in
#' GENCODE, which the OrgDb route keeps under their IDs), and
#' `attr(x, "gene_ids")` keeps the original IDs. Names from an annotation and
#' from the OrgDb can differ for some genes.
#'
#' @section Checks on the files:
#' Every file must list the same genes (transcripts) as the first, in the same
#' order, and, with `validate_reference`, the same transcripts for each gene
#' (gene-level files) or the same gene for each transcript (isoform files):
#' the same gene ID with a different set of transcripts means the files were
#' quantified against different references, and their values are not
#' comparable. [rsem_check()] reports such problems for all files at once
#' instead of stopping at the first.
#'
#' @section Posterior estimates:
#' With `--calc-pme` or `--calc-ci`, RSEM adds posterior means, their standard
#' deviation and credibility intervals. `"pme_count"`, `"pme_tpm"` and
#' `"pme_fpkm"` are posterior means and add up over rows like counts. The
#' standard deviation, the interval bounds, the coefficient of quartile
#' variation and `"isopct_pme"` do not add up (the interval of a sum is not the sum
#' of the intervals), so they need `id = "original"` or `level = "isoform"`
#' and are never merged. A file without the column stops with the option to
#' use. A coefficient of quartile variation that RSEM could not compute (both
#' quartiles are 0, written as `nan`) is `NA`.
#'
#' @section Effective length:
#' `expr = "length"` reads `effective_length`. RSEM sets it to 0 for a
#' transcript shorter than the mean fragment length, and the gene-level value
#' is an average over transcripts weighted by their abundance, so a gene can
#' have length 0 in some samples and a positive length in others.
#' `tximport` and `DESeq2::DESeqDataSetFromTximport()` stop on such zeros. The
#' function reports how many genes have a zero length and keeps them as 0
#' (attribute `"zero_length"`: the number of zero genes per sample). Common
#' choices are to remove those genes or to import the transcript-level
#' results; replacing 0 by 1 changes the length offset of those genes. See
#' [rsem_length_profile()].
#'
#' @section Transcripts:
#' With `level = "isoform"`, the `*.isoforms.results` files are read; rows are
#' transcript IDs as in the file and are never merged. `expr` can also be
#' `"isopct"` (the transcript's share of its gene's abundance, in percent) and
#' `"isopct_pme"` (the same from the posterior mean TPM).
#' The attribute `"tx2gene"` holds the transcript-to-gene table
#' (`transcript_id`, `gene_id`, and `gene`: the symbol with `id = "symbol"`,
#' else the `gene_id`), which [rsem_select_isoforms()] uses.
#'
#' @param data Directory that contains the RSEM files.
#' @param metadata Path to the metadata CSV, or `NULL` to find the files in
#'   `data` (see Input).
#' @param species `"human"` or `"mouse"`. Required with `id = "symbol"`; the
#'   function stops if the gene IDs are of the other species. Not used with
#'   `id = "original"`, which works for any organism.
#' @param expr `"count"` (`expected_count`), `"tpm"` (`TPM`), `"fpkm"`
#'   (`FPKM`), `"length"` (`effective_length`) or, for transcripts, `"isopct"`
#'   (`IsoPct`); with `--calc-pme` or `--calc-ci` also the posterior values
#'   `"pme_count"` (`posterior_mean_count`), `"pme_tpm"` (`pme_TPM`),
#'   `"pme_fpkm"` (`pme_FPKM`), `"count_sd"`
#'   (`posterior_standard_deviation_of_count`) and, for transcripts,
#'   `"isopct_pme"` (`IsoPct_from_pme_TPM`); with `--calc-ci` the credibility
#'   interval `"tpm_ci_lower"`, `"tpm_ci_upper"`, `"fpkm_ci_lower"`,
#'   `"fpkm_ci_upper"` and the coefficient of quartile variation `"tpm_cqv"`,
#'   `"fpkm_cqv"` (`(Q3 - Q1) / (Q3 + Q1)` of the posterior samples; the
#'   level of the interval is set by `--ci-credibility-level`, 95% by default,
#'   and is not recorded in the file). See Posterior estimates. Required.
#' @param level `"gene"` (`*.genes.results`, the default) or `"isoform"`
#'   (`*.isoforms.results`).
#' @param id `"symbol"` (the default), `"original"` (the IDs of the file) or
#'   `"annotation"` (names from `annotation`).
#' @param validate_reference Check that every file lists the same transcripts
#'   for each gene (the column `transcript_id(s)`), as set; files made with
#'   different references stop the function. Default `TRUE`; skipped when the
#'   column is missing.
#' @param allow_reorder Accept files that have the same IDs in a different
#'   order, and put them in the order of the first file. Default `FALSE`: a
#'   different order stops the function, because it usually means the files do
#'   not belong together.
#' @param annotation With `id = "annotation"`: a GTF file or a data frame with
#'   the columns `gene_id` and `gene_name`, as in [rsem_biotype_profile()]; for
#'   a RefSeq GFF3 reference, the GTF that `rsem-prepare-reference` writes.
#' @param sample_from With `metadata = NULL`: `"filename"` (the default) or
#'   `"parent_dir"`.
#'
#' @return A data frame with one row per gene symbol (or gene ID, or
#'   transcript; as row names) and one column per sample, in the order of the
#'   metadata. A message reports how many genes received a symbol.
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
#'   # find the files instead of listing them
#'   tpm2 <- rsem_read(file.path(ex, "data"), metadata = NULL,
#'                     species = "human", expr = "tpm")
#' }
#' @export
rsem_read <- function(data = "data", metadata = "metadata.csv", species, expr,
                      level = "gene", id = "symbol", sample_from = "filename",
                      validate_reference = TRUE, allow_reorder = FALSE, annotation = NULL) {
  id <- .choice(id, c("symbol", "original", "annotation"), "id")
  if (id == "symbol" && missing(species)) .stopf("`species` must be \"human\" or \"mouse\".")
  if (missing(expr)) .stopf("`expr` must be one of: %s.", paste(sprintf("\"%s\"", .read_exprs), collapse = ", "))
  # with the IDs of the files or an annotation file, no annotation package is
  # used and any organism works
  species <- if (id != "symbol") NULL else .choice(species, c("human", "mouse"), "species")
  level <- .choice(level, c("gene", "isoform"), "level")
  expr <- .choice(expr, c(setdiff(.read_exprs, .isoform_only), if (level == "isoform") .isoform_only), "expr")
  column <- .read_columns[[expr]]
  if (!is.logical(validate_reference) || length(validate_reference) != 1L || is.na(validate_reference)) .stopf("`validate_reference` must be TRUE or FALSE.")
  if (!is.logical(allow_reorder) || length(allow_reorder) != 1L || is.na(allow_reorder)) .stopf("`allow_reorder` must be TRUE or FALSE.")
  if (id == "annotation" && is.null(annotation)) .stopf("id = \"annotation\" needs `annotation`: a GTF file or a data frame with gene_id and gene_name.")
  if (id != "annotation" && !is.null(annotation)) .stopf("`annotation` is used with id = \"annotation\" only.")
  if (expr %in% .not_additive && id != "original" && level == "gene") {
    .stopf("expr = \"%s\" cannot be summed over rows that share a name; use id = \"original\" (or level = \"isoform\").", expr)
  }

  samples <- if (is.null(metadata)) {
    .discover_samples(data, level, .choice(sample_from, c("filename", "parent_dir"), "sample_from"))
  } else {
    .read_samples(data, metadata)
  }
  if (expr %in% .posterior_exprs && !column %in% .rsem_header(samples$path[[1L]])) {
    .stopf(
      "%s has no column '%s'; run rsem-calculate-expression with %s.", samples$ref[[1L]], column,
      if (expr %in% .needs_calc_ci) "--calc-ci" else "--calc-pme or --calc-ci"
    )
  }
  ann <- if (id == "annotation") .read_annotation(annotation, need_type = FALSE) else NULL
  if (level == "isoform") return(.read_isoforms(samples, column, species, id, expr, ann, allow_reorder))

  tx_col <- "transcript_id(s)"
  check_tx <- validate_reference && tx_col %in% .rsem_header(samples$path[[1L]])
  mat <- .read_matrix(samples, column, allow_reorder = allow_reorder, also = if (check_tx) tx_col, also_set = TRUE)
  stable <- .strip_ensembl_version(rownames(mat))
  if (id == "symbol") .check_species(stable, species)
  zero_note <- ""
  zl <- NULL
  if (expr == "length") {
    zl <- colSums(mat == 0)
    zero_note <- sprintf(
      "\n  genes with effective length 0: %s in at least one sample, %s in every sample",
      .fmt(sum(rowSums(mat == 0) > 0)), .fmt(sum(rowSums(mat == 0) == ncol(mat)))
    )
  }

  if (id == "original") {
    message(sprintf(
      "rsem_read: %d samples, %d genes (expr = %s, original IDs)%s",
      ncol(mat), nrow(mat), expr, zero_note
    ))
    res <- as.data.frame(mat, check.names = FALSE)
    attr(res, "gene_ids") <- data.frame(gene = rownames(mat), gene_id = rownames(mat), stringsAsFactors = FALSE)
    if (!is.null(zl)) attr(res, "zero_length") <- zl
    return(res)
  }

  if (id == "annotation") {
    symbol <- .names_from_annotation(rownames(mat), ann)
    source <- "annotation file"
  } else {
    orgdb <- .orgdb_package(species)
    symbol <- .lookup_symbols(stable, orgdb)
    source <- sprintf("%s %s", orgdb, as.character(utils::packageVersion(orgdb)))
  }
  out <- if (expr == "length") {
    .weighted_length_by_symbol(mat, .read_matrix(samples, "TPM", allow_reorder = allow_reorder), symbol)
  } else {
    .sum_by_symbol(mat, symbol)
  }

  has_symbol <- !is.na(symbol)
  message(sprintf(
    paste0(
      "rsem_read: %d samples, %d genes -> %d rows (%sexpr = %s, %s)\n",
      "  with symbol: %d genes (%d symbols merged from several IDs); ",
      "kept as original IDs: %d genes%s"
    ),
    ncol(mat), nrow(mat), nrow(out), if (is.null(species)) "" else paste0("species = ", species, ", "), expr, source,
    sum(has_symbol), sum(table(symbol[has_symbol]) > 1L), sum(!has_symbol), zero_note
  ))

  res <- as.data.frame(out, check.names = FALSE)
  attr(res, "gene_ids") <- data.frame(
    gene = ifelse(is.na(symbol), rownames(mat), symbol),
    gene_id = rownames(mat), stringsAsFactors = FALSE
  )
  if (!is.null(zl)) attr(res, "zero_length") <- zl
  res
}

# What rsem_read() can return, and the RSEM column each value comes from. The
# posterior columns exist only with --calc-pme or --calc-ci; posterior means
# (like counts and TPM) add up over rows, standard deviations and interval
# bounds do not.
.read_columns <- list(
  count = "expected_count", tpm = "TPM", fpkm = "FPKM", length = "effective_length", isopct = "IsoPct",
  pme_count = "posterior_mean_count", pme_tpm = "pme_TPM", pme_fpkm = "pme_FPKM",
  count_sd = "posterior_standard_deviation_of_count",
  tpm_ci_lower = "TPM_ci_lower_bound", tpm_ci_upper = "TPM_ci_upper_bound",
  fpkm_ci_lower = "FPKM_ci_lower_bound", fpkm_ci_upper = "FPKM_ci_upper_bound",
  tpm_cqv = "TPM_coefficient_of_quartile_variation", fpkm_cqv = "FPKM_coefficient_of_quartile_variation",
  isopct_pme = "IsoPct_from_pme_TPM"
)
.read_exprs <- names(.read_columns)
.posterior_exprs <- c(
  "pme_count", "pme_tpm", "pme_fpkm", "count_sd", "tpm_ci_lower", "tpm_ci_upper", "fpkm_ci_lower", "fpkm_ci_upper",
  "tpm_cqv", "fpkm_cqv", "isopct_pme"
)
# Cannot be summed over rows that share a name (a sum of intervals, deviations,
# quartile variations or percentages is not the value of the sum).
.not_additive <- c(
  "count_sd", "tpm_ci_lower", "tpm_ci_upper", "fpkm_ci_lower", "fpkm_ci_upper", "tpm_cqv", "fpkm_cqv", "isopct_pme"
)
# Written only with --calc-ci (the others with --calc-pme or --calc-ci).
.needs_calc_ci <- c("tpm_ci_lower", "tpm_ci_upper", "fpkm_ci_lower", "fpkm_ci_upper", "tpm_cqv", "fpkm_cqv")
.isoform_only <- c("isopct", "isopct_pme")

# Name of each row from an annotation (gene_id -> gene_name). A file ID matches
# the annotation exactly or after its Ensembl version is dropped. An ID with
# the name that `rsem-calculate-expression --append-names` adds (`<gene_id>_<
# gene_name>`) is accepted only when the part after some `_` is exactly the
# gene_name of the gene whose gene_id is the part before it; so `gene1_2` is
# not read as the gene `gene1` just because that gene exists, and underscores
# inside IDs or names are never cut at the wrong place.
.names_from_annotation <- function(ids, ann) {
  ann <- ann[!is.na(ann$gene_name) & nzchar(ann$gene_name), , drop = FALSE]
  ann <- ann[!duplicated(ann$gene_id), , drop = FALSE]
  drop_version <- function(x) sub("^(ENS[A-Z]*G[0-9]+)\\.[0-9]+$", "\\1", x)
  key <- rep(NA_character_, length(ids))
  exact <- ids %in% ann$gene_id
  key[exact] <- ids[exact]
  plain <- is.na(key) & drop_version(ids) %in% ann$gene_id
  key[plain] <- drop_version(ids[plain])
  name_of <- stats::setNames(ann$gene_name, ann$gene_id)
  for (i in which(is.na(key) & grepl("_", ids, fixed = TRUE))) {
    for (p in gregexpr("_", ids[[i]], fixed = TRUE)[[1L]]) {
      pre <- substr(ids[[i]], 1L, p - 1L)
      rest <- substr(ids[[i]], p + 1L, nchar(ids[[i]]))
      g <- c(pre, drop_version(pre))
      g <- g[g %in% ann$gene_id][1L]
      if (!is.na(g) && identical(unname(name_of[[g]]), rest)) {
        key[[i]] <- g
        break
      }
    }
  }
  ann$gene_name[match(key, ann$gene_id)]
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

# Only Ensembl gene IDs lose their version suffix and the gene name that
# `rsem-calculate-expression --append-names` adds (`ENSG00000141510.18_TP53`);
# other IDs (spike-ins, numeric IDs, `..._PAR_Y`) are returned unchanged.
.strip_ensembl_version <- function(x) {
  sub("^(ENS[A-Z]*G[0-9]+)(\\.[0-9]+)?(_(?!PAR_Y).+)?$", "\\1", as.character(x), perl = TRUE)
}

# RSEM result files may be gzip-compressed. Read them through a connection so
# that no extra package is needed for .gz files.
.rsem_header <- function(path) {
  line <- if (grepl("\\.gz$", path)) {
    con <- gzfile(path, "r"); on.exit(close(con)); readLines(con, n = 1L, warn = FALSE)
  } else {
    con <- file(path, "r"); on.exit(close(con)); readLines(con, n = 1L, warn = FALSE)
  }
  strsplit(line, "\t", fixed = TRUE)[[1L]]
}
.fread_rsem <- function(path, ...) {
  if (grepl("\\.gz$", path)) {
    con <- gzfile(path, "r"); on.exit(close(con))
    data.table::fread(text = readLines(con, warn = FALSE), sep = "\t", data.table = FALSE, showProgress = FALSE, ...)
  } else {
    data.table::fread(path, sep = "\t", data.table = FALSE, showProgress = FALSE, ...)
  }
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

# Find RSEM files under `data` and name the samples by file or folder name.
.discover_samples <- function(data, level, sample_from) {
  if (!is.character(data) || length(data) != 1L || is.na(data) || !dir.exists(data)) {
    .stopf("`data` must be an existing directory: %s", toString(data))
  }
  suffix <- c(gene = ".genes.results", isoform = ".isoforms.results", cnt = ".cnt", model = ".model", theta = ".theta")[[level]]
  pattern <- paste0(gsub(".", "\\.", suffix, fixed = TRUE), if (level %in% c("cnt", "model", "theta")) "$" else "(\\.gz)?$")
  path <- sort(list.files(data, pattern = pattern, recursive = TRUE), method = "radix")
  if (!length(path)) .stopf("No '*%s' files found under '%s'.", suffix, data)
  if (sample_from == "filename") {
    sample_id <- sub(pattern, "", basename(path))
  } else {
    if (any(dirname(path) == ".")) {
      .stopf("sample_from = \"parent_dir\" needs every file in its own folder under '%s'; some are directly in it.", data)
    }
    sample_id <- basename(dirname(path))
    if (level %in% c("cnt", "model", "theta")) sample_id <- sub("\\.stat$", "", sample_id)  # RSEM writes <name>.stat/<name>.cnt
  }
  if (anyDuplicated(sample_id)) {
    d <- unique(sample_id[duplicated(sample_id)])
    .stopf(
      "Several files give the same sample ID (%s): %s", paste(d, collapse = ", "),
      paste(path[sample_id %in% d], collapse = "; ")
    )
  }
  o <- order(sample_id, method = "radix")
  out <- data.frame(
    sample_id = sample_id[o], ref = path[o], path = file.path(data, path[o]),
    stringsAsFactors = FALSE
  )
  message(sprintf(
    "%s: found %d files under '%s' (sample ID from the %s):\n%s",
    if (level == "cnt") "rsem_mapping_qc" else if (level == "model") "rsem_model_qc" else if (level == "theta") "rsem_noise_qc" else "rsem_read", nrow(out), data, if (sample_from == "filename") "file name" else "folder name",
    paste(sprintf("  %s <- %s", out$sample_id, out$ref), collapse = "\n")
  ))
  out
}

# `also`: a second ID column that must be the same in every file (the gene of
# each transcript), checked against the first file.
.read_matrix <- function(samples, column, key = "gene_id", also = NULL, allow_reorder = FALSE, also_set = FALSE) {
  required <- c(key, column, also)
  row_id <- NULL
  also_id <- NULL
  mat <- NULL
  classes <- list(c(key, also), column)
  names(classes) <- c("character", "numeric")
  for (i in seq_len(nrow(samples))) {
    label <- sprintf("%s (sample %s)", samples$ref[[i]], samples$sample_id[[i]])
    header <- .rsem_header(samples$path[[i]])
    missing_cols <- setdiff(required, header)
    if (length(missing_cols)) {
      .stopf("%s is missing the column(s): %s", label, paste(missing_cols, collapse = ", "))
    }

    # (Q3 - Q1) / (Q3 + Q1) is not defined when both quartiles are 0; such a value
    # (written as nan) is kept as NA for the coefficient of quartile variation only.
    nan_ok <- grepl("coefficient_of_quartile_variation$", column)
    x <- suppressWarnings(.fread_rsem(
      samples$path[[i]], select = required, colClasses = classes,
      na.strings = if (nan_ok) c("NA", "nan", "-nan", "NaN", "-NaN") else "NA"
    ))
    if (!is.numeric(x[[column]])) .stopf("%s has non-numeric values in '%s'.", label, column)
    if ((!nan_ok && anyNA(x[[column]])) || anyNA(x[[key]])) .stopf("%s contains missing values.", label)
    if (any(is.infinite(x[[column]])) || any(x[[column]] < 0, na.rm = TRUE)) {
      .stopf("%s has infinite or negative values in '%s' (RSEM writes none); the file may be damaged. rsem_check() lists the rows.", label, column)
    }
    if (anyDuplicated(x[[key]])) .stopf("%s contains duplicated %s values.", label, key)

    if (is.null(row_id)) {
      row_id <- x[[key]]
      mat <- matrix(
        NA_real_, length(row_id), nrow(samples),
        dimnames = list(row_id, samples$sample_id)
      )
    } else if (!identical(x[[key]], row_id)) {
      if (allow_reorder && setequal(x[[key]], row_id)) {
        x <- x[match(row_id, x[[key]]), , drop = FALSE]
      } else {
        .stopf(
          "%s does not have the same %s in the same order as %s%s.",
          label, if (key == "gene_id") "genes" else "transcripts", samples$ref[[1]],
          if (setequal(x[[key]], row_id)) " (the same IDs in a different order; allow_reorder = TRUE accepts this)" else ""
        )
      }
    }
    if (!is.null(also)) {
      if (is.null(also_id)) {
        also_id <- x[[also]]
      } else if (!identical(x[[also]], also_id)) {
        same <- if (also_set) {
          canon <- function(v) vapply(strsplit(v, ",", fixed = TRUE), function(z) paste(sort(z), collapse = ","), character(1))
          diff <- which(x[[also]] != also_id)
          canon(x[[also]][diff]) == canon(also_id[diff])
        } else {
          FALSE
        }
        if (!all(same)) {
          bad <- if (also_set) which(x[[also]] != also_id)[!same] else which(x[[also]] != also_id)
          .stopf(
            "%s has a different %s for %s %s than %s (e.g. %s): the files come from different references; `validate_reference = FALSE` skips this check.",
            label, also, .fmt(length(bad)), if (key == "gene_id") "genes" else "transcripts", samples$ref[[1]], x[[key]][bad[[1L]]]
          )
        }
      }
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

.symbol_keys <- function(mat, symbol) {
  unmapped_id <- rownames(mat)[is.na(symbol)]
  clash <- intersect(unmapped_id, symbol[!is.na(symbol)])
  if (length(clash)) {
    .stopf(
      "A gene_id without a symbol equals a gene symbol and would be merged with it: %s",
      paste(utils::head(clash, 5L), collapse = ", ")
    )
  }
  ifelse(is.na(symbol), rownames(mat), symbol)
}

.sum_by_symbol <- function(mat, symbol) {
  key <- .symbol_keys(mat, symbol)
  out <- rowsum(mat, group = key, reorder = FALSE)
  if (!isTRUE(all.equal(colSums(out), colSums(mat), tolerance = 1e-9))) {
    .stopf("Internal error: the column sums changed while merging rows.")
  }
  out
}

# Effective lengths of rows that share a symbol: TPM-weighted mean, or the
# plain mean where none of the rows has TPM in that sample.
.weighted_length_by_symbol <- function(len, tpm, symbol) {
  key <- .symbol_keys(len, symbol)
  w <- rowsum(tpm, key, reorder = FALSE)
  wsum <- rowsum(len * tpm, key, reorder = FALSE)
  n <- as.vector(table(key)[rownames(w)])
  plain <- rowsum(len, key, reorder = FALSE) / n
  out <- plain
  pos <- w > 0
  out[pos] <- wsum[pos] / w[pos]
  out
}

.read_isoforms <- function(samples, column, species, id, expr, ann = NULL, allow_reorder = FALSE) {
  mat <- .read_matrix(samples, column, key = "transcript_id", also = "gene_id", allow_reorder = allow_reorder)
  g <- .fread_rsem(samples$path[[1]], select = c("transcript_id", "gene_id"), colClasses = "character")
  if (!identical(g$transcript_id, rownames(mat))) .stopf("Internal error: transcript order changed.")
  stable <- .strip_ensembl_version(g$gene_id)
  if (id == "symbol") .check_species(stable, species)
  gene <- g$gene_id
  if (id == "symbol") {
    s <- .lookup_symbols(stable, .orgdb_package(species))
    gene <- ifelse(is.na(s), g$gene_id, s)
  } else if (id == "annotation") {
    s <- .names_from_annotation(g$gene_id, ann)
    gene <- ifelse(is.na(s), g$gene_id, s)
  }
  message(sprintf(
    "rsem_read: %d samples, %s transcripts of %s genes (%sexpr = %s, level = isoform)",
    ncol(mat), .fmt(nrow(mat)), .fmt(length(unique(g$gene_id))),
    if (is.null(species)) paste0(if (id == "annotation") "names from the annotation file" else "original IDs", ", ") else paste0("species = ", species, ", "), expr
  ))
  res <- as.data.frame(mat, check.names = FALSE)
  attr(res, "tx2gene") <- data.frame(
    transcript_id = g$transcript_id, gene_id = g$gene_id, gene = gene,
    stringsAsFactors = FALSE
  )
  res
}
