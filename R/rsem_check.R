#' Check RSEM result files before reading them
#'
#' Looks at every RSEM result file and reports what would make the data
#' unreliable or stop a later step, all at once instead of stopping at the
#' first problem. Nothing is read into a table and nothing is changed. Use it
#' before [rsem_read()] on a new data set, or when [rsem_read()] stops.
#'
#' @section Checks:
#' Per file: the files can be read; the ID column (`gene_id` or
#' `transcript_id`) is present; `expected_count`, `TPM` and
#' `effective_length` are present; the values are numeric, not missing, finite
#' and not negative; the IDs are not duplicated; the effective lengths that
#' are 0 (the genes that stop `DESeqDataSetFromTximport()`); a TPM total that
#' is not 1e6 (RSEM's TPM adds up to one million per sample, within 0.1%
#' here); a sample without any expected count.
#'
#' Across files, compared with the first readable file: the same IDs (an
#' error where IDs are missing or extra; information where only the order
#' differs, which `rsem_read(allow_reorder = TRUE)` accepts); the same
#' transcripts for each gene (`transcript_id(s)`, compared as sets) in
#' gene-level files, or the same gene for each transcript in isoform files,
#' which differ when the files were quantified against different references.
#'
#' @section Severity:
#' `"ERROR"`: [rsem_read()] stops or the files cannot be compared.
#' `"WARNING"`: the files can be read but the values need care.
#' `"INFO"`: noted, nothing to fix.
#'
#' @param data,metadata,sample_from As in [rsem_read()]; `metadata = NULL`
#'   finds the files.
#' @param level `"gene"` (`*.genes.results`) or `"isoform"`
#'   (`*.isoforms.results`).
#'
#' @return A data frame with one row per finding: `severity`, `check`,
#'   `sample_id` (`NA` for a finding about all files), `n` (how many rows or
#'   genes) and `detail`. Zero rows when nothing was found.
#'
#' @examples
#' ex <- system.file("extdata", "example", package = "rsemverse")
#' rsem_check(file.path(ex, "data"))
#' @export
rsem_check <- function(data = "data", metadata = NULL, level = "gene", sample_from = "filename") {
  level <- .choice(level, c("gene", "isoform"), "level")
  key <- if (level == "gene") "gene_id" else "transcript_id"
  second <- if (level == "gene") "transcript_id(s)" else "gene_id"
  found <- list()
  add <- function(severity, check, sample_id, detail, n = NA_integer_) {
    found[[length(found) + 1L]] <<- data.frame(
      severity = severity, check = check, sample_id = sample_id, n = as.integer(n), detail = detail,
      stringsAsFactors = FALSE
    )
  }
  finish <- function(n_files) {
    out <- if (length(found)) do.call(rbind, found) else data.frame(
      severity = character(0), check = character(0), sample_id = character(0), n = integer(0), detail = character(0),
      stringsAsFactors = FALSE
    )
    out <- out[order(match(out$severity, c("ERROR", "WARNING", "INFO")), out$sample_id, out$check, na.last = TRUE), , drop = FALSE]
    rownames(out) <- NULL
    message(sprintf(
      "rsem_check: %d files; %d errors, %d warnings, %d notes", n_files,
      sum(out$severity == "ERROR"), sum(out$severity == "WARNING"), sum(out$severity == "INFO")
    ))
    out
  }

  samples <- tryCatch(
    suppressMessages(
      if (is.null(metadata)) .discover_samples(data, level, .choice(sample_from, c("filename", "parent_dir"), "sample_from"))
      else .read_samples(data, metadata)
    ),
    error = function(e) e
  )
  if (inherits(samples, "error")) {
    add("ERROR", "files", NA_character_, conditionMessage(samples))
    return(finish(0L))
  }

  ref <- NULL
  value_cols <- c("expected_count", "TPM", "effective_length")
  for (i in seq_len(nrow(samples))) {
    sid <- samples$sample_id[[i]]
    path <- samples$path[[i]]
    header <- tryCatch(.rsem_header(path), error = function(e) NULL)
    if (is.null(header) || !length(header) || anyNA(header)) {
      add("ERROR", "unreadable", sid, sprintf("%s cannot be read as a tab-separated RSEM file.", samples$ref[[i]]))
      next
    }
    if (!key %in% header) {
      add("ERROR", "required_columns", sid, sprintf("%s has no '%s' column.", samples$ref[[i]], key))
      next
    }
    missing_cols <- setdiff(value_cols, header)
    if (length(missing_cols)) {
      add("WARNING", "required_columns", sid, sprintf("missing column(s): %s", paste(missing_cols, collapse = ", ")), length(missing_cols))
    }
    id_cols <- intersect(c(key, second), header)
    x <- tryCatch(
      suppressWarnings(.fread_rsem(path, colClasses = list(character = id_cols), na.strings = c("NA", "nan", "-nan", "NaN", "-NaN"))),
      error = function(e) NULL
    )
    if (is.null(x)) {
      add("ERROR", "unreadable", sid, sprintf("%s cannot be read as a tab-separated RSEM file.", samples$ref[[i]]))
      next
    }
    ids <- x[[key]]
    n_dup <- sum(duplicated(ids))
    if (anyNA(ids) || any(!nzchar(ids))) add("ERROR", "missing_ids", sid, "empty or missing IDs", sum(is.na(ids) | !nzchar(ids)))
    if (n_dup) add("ERROR", "duplicated_ids", sid, sprintf("%d duplicated %s values", n_dup, key), n_dup)
    if (second %in% names(x)) {
      sec <- x[[second]]
      n_sec <- sum(is.na(sec) | !nzchar(sec))
      if (n_sec) add("ERROR", "missing_ids", sid, sprintf("empty or missing '%s' for %d rows", second, n_sec), n_sec)
    }

    for (col in intersect(value_cols, names(x))) {
      v <- x[[col]]
      if (!is.numeric(v)) {
        add("ERROR", "non_numeric", sid, sprintf("column '%s' is not numeric", col))
        next
      }
      if (anyNA(v)) add("ERROR", "missing_values", sid, sprintf("column '%s' has missing values", col), sum(is.na(v)))
      if (any(is.infinite(v))) add("ERROR", "nonfinite_values", sid, sprintf("column '%s' has infinite values", col), sum(is.infinite(v)))
      if (any(v < 0, na.rm = TRUE)) add("ERROR", "negative_values", sid, sprintf("column '%s' has negative values", col), sum(v < 0, na.rm = TRUE))
      if (col == "effective_length" && any(v == 0, na.rm = TRUE)) {
        add("WARNING", "zero_effective_length", sid,
            "RSEM sets the effective length of a transcript shorter than the mean fragment length to 0; tximport and DESeqDataSetFromTximport() stop on zero lengths",
            sum(v == 0, na.rm = TRUE))
      }
      if (col == "TPM" && !anyNA(v)) {
        total <- sum(v)
        # RSEM prints TPM with two decimals: the total can be off by up to 0.005 per row
        # from rounding alone (the typical error is about 0.003 * sqrt(rows)), far below 1e3.
        if (abs(total - 1e6) > 1e3) add("WARNING", "tpm_total", sid, sprintf(
          "TPM adds up to %s, not 1e6 (difference %s; rounding of %d rows explains at most %s)",
          format(signif(total, 7), scientific = FALSE), format(signif(total - 1e6, 4), scientific = FALSE),
          length(v), format(signif(0.005 * length(v), 3), scientific = FALSE)
        ), length(v))
      }
      if (col == "expected_count" && !anyNA(v) && sum(v) == 0) add("WARNING", "no_counts", sid, "the expected counts are 0 for every row")
    }

    # Posterior columns written by --calc-pme / --calc-ci: finite, not negative; a lower bound
    # must not be above its upper bound. (A coefficient of quartile variation can be nan.)
    post_cols <- intersect(unlist(.read_columns[.posterior_exprs], use.names = FALSE), names(x))
    for (col in post_cols) {
      v <- x[[col]]
      if (!is.numeric(v)) {
        add("ERROR", "posterior_values", sid, sprintf("column '%s' is not numeric", col))
        next
      }
      nan_ok <- grepl("coefficient_of_quartile_variation$", col)
      if (!nan_ok && anyNA(v)) add("ERROR", "posterior_values", sid, sprintf("column '%s' has missing values", col), sum(is.na(v)))
      if (any(is.infinite(v)) || any(v < 0, na.rm = TRUE)) {
        add("ERROR", "posterior_values", sid, sprintf("column '%s' has infinite or negative values", col), sum(is.infinite(v) | v < 0, na.rm = TRUE))
      }
    }
    for (pair in list(c("TPM_ci_lower_bound", "TPM_ci_upper_bound"), c("FPKM_ci_lower_bound", "FPKM_ci_upper_bound"))) {
      if (all(pair %in% names(x)) && is.numeric(x[[pair[[1L]]]]) && is.numeric(x[[pair[[2L]]]])) {
        bad <- sum(x[[pair[[1L]]]] > x[[pair[[2L]]]], na.rm = TRUE)
        if (bad) add("ERROR", "posterior_interval", sid, sprintf("%s is above %s for %d rows", pair[[1L]], pair[[2L]], bad), bad)
      }
    }

    if (is.null(ref)) {
      ref <- list(
        ids = ids, second = if (second %in% names(x)) x[[second]] else NULL, ref = samples$ref[[i]], sid = sid,
        length = if (level == "isoform" && "length" %in% names(x) && is.numeric(x[["length"]])) x[["length"]] else NULL
      )
      next
    }
    if (!identical(ids, ref$ids)) {
      only_here <- setdiff(ids, ref$ids)
      only_ref <- setdiff(ref$ids, ids)
      if (length(only_here) || length(only_ref)) {
        add("ERROR", "id_set", sid, sprintf(
          "%d %s not in %s, %d of %s missing (e.g. %s)", length(only_here), key, ref$ref, length(only_ref), ref$ref,
          paste(utils::head(c(only_here, only_ref), 2L), collapse = ", ")
        ), length(only_here) + length(only_ref))
        next
      }
      add("INFO", "id_order", sid, sprintf("the same IDs as %s in a different order; rsem_read(allow_reorder = TRUE) accepts this", ref$ref))
    }
    # The sequence length of a transcript is a property of the reference, not of the sample
    # (the gene-level length is an abundance-weighted mean and does vary), so a transcript
    # whose length differs from the first file comes from another sequence.
    if (level == "isoform" && !is.null(ref$length) && "length" %in% names(x) && is.numeric(x[["length"]])) {
      len_here <- x[["length"]][match(ref$ids, ids)]
      diff_len <- which(abs(len_here - ref$length) > 0.5)
      if (length(diff_len)) {
        add("ERROR", "transcript_length", sid, sprintf(
          "the length of %d transcripts differs from %s (e.g. %s: %s vs %s): the files seem to come from different references",
          length(diff_len), ref$ref, ref$ids[[diff_len[[1L]]]], format(len_here[[diff_len[[1L]]]]), format(ref$length[[diff_len[[1L]]]])
        ), length(diff_len))
      }
    }
    if (second %in% names(x) && !is.null(ref$second)) {
      here <- x[[second]][match(ref$ids, ids)]
      differ <- which(here != ref$second)
      if (length(differ)) {
        if (level == "gene") {
          canon <- function(v) vapply(strsplit(v, ",", fixed = TRUE), function(z) paste(sort(z), collapse = ","), character(1))
          differ <- differ[canon(here[differ]) != canon(ref$second[differ])]
        }
        if (length(differ)) {
          add("ERROR", "reference_mismatch", sid, sprintf(
            "%s differs from %s for %d %s (e.g. %s): the files seem to come from different references",
            second, ref$ref, length(differ), if (level == "gene") "genes" else "transcripts", ref$ids[differ[[1L]]]
          ), length(differ))
        }
      }
    } else if (!is.null(ref$second) && !second %in% names(x)) {
      add("WARNING", "reference_check", sid, sprintf("no '%s' column, so the reference cannot be compared with %s", second, ref$ref))
    }
  }
  finish(nrow(samples))
}
