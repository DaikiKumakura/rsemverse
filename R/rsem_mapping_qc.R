#' Read assignment statistics from RSEM `.cnt` files
#'
#' Reads the `*.cnt` file that `rsem-calculate-expression` writes (in
#' `<name>.stat/<name>.cnt`) for every sample and returns one row per sample
#' with the read counts and fractions that describe how the reads were
#' assigned. The numbers come from the alignments RSEM was given, not from the
#' abundance estimates.
#'
#' @section What the file contains:
#' Following RSEM's `cnt_file_description.txt`: `N0` reads that could not be
#' aligned, `N1` alignable reads, `N2` reads filtered because they have too
#' many alignments, and `N_tot = N0 + N1 + N2`; `nUnique` reads aligned
#' uniquely to one gene and `nMulti` reads aligned to several genes
#' (`nUnique + nMulti = N1`); `nUncertain` reads aligned to several places in
#' the reference, which includes ambiguity between isoforms of the same gene,
#' so it is not a subset of `nMulti`; `nHits`, the total number of alignments;
#' and the read type (0 single-end without quality, 1 single-end with quality,
#' 2 paired-end without quality, 3 paired-end with quality). The histogram of
#' alignments per read that follows in the file is not used.
#'
#' A file that breaks these identities (`N_tot` is not `N0 + N1 + N2`, or
#' `nUnique + nMulti` is not `N1`) stops the function instead of returning
#' numbers that cannot be right.
#'
#' @section Files:
#' As in [rsem_read()]: `metadata = NULL` finds every `*.cnt` file under
#' `data` (the sample is named after the file, or with `sample_from =
#' "parent_dir"` after its folder without the `.stat` suffix), or a metadata
#' CSV with `sample_id` and `ref` names the files.
#'
#' @inheritParams rsem_read
#'
#' @return A data frame with one row per sample: `sample_id`, the counts
#'   `n_total`, `n_unalignable`, `n_alignable`, `n_filtered`, `n_unique`,
#'   `n_multi`, `n_uncertain`, `n_hits`, `read_type`, and the fractions
#'   `alignable_fraction` (`N1 / N_tot`), `unalignable_fraction`
#'   (`N0 / N_tot`), `filtered_fraction` (`N2 / N_tot`), `unique_gene_fraction`
#'   (`nUnique / N1`), `multi_gene_fraction` (`nMulti / N1`),
#'   `uncertain_fraction` (`nUncertain / N1`) and `alignments_per_read`
#'   (`nHits / N1`). Fractions are `NA` where the denominator is 0.
#'
#' @examples
#' f <- tempfile(fileext = ".cnt")
#' writeLines(c("100 850 50 1000", "800 50 120", "1500 2", "0\t100", "1\t700", "Inf\t50"), f)
#' rsem_mapping_qc(dirname(f), metadata = NULL)[, c("sample_id", "n_total", "alignable_fraction")]
#' @export
rsem_mapping_qc <- function(data = "data", metadata = NULL, sample_from = "filename") {
  samples <- if (is.null(metadata)) {
    .discover_samples(data, "cnt", .choice(sample_from, c("filename", "parent_dir"), "sample_from"))
  } else {
    .read_samples(data, metadata)
  }
  rows <- lapply(seq_len(nrow(samples)), function(i) .read_cnt(samples$path[[i]], samples$sample_id[[i]], samples$ref[[i]]))
  out <- do.call(rbind, rows)
  ratio <- function(a, b) ifelse(b > 0, a / b, NA_real_)
  out$alignable_fraction <- ratio(out$n_alignable, out$n_total)
  out$unalignable_fraction <- ratio(out$n_unalignable, out$n_total)
  out$filtered_fraction <- ratio(out$n_filtered, out$n_total)
  out$unique_gene_fraction <- ratio(out$n_unique, out$n_alignable)
  out$multi_gene_fraction <- ratio(out$n_multi, out$n_alignable)
  out$uncertain_fraction <- ratio(out$n_uncertain, out$n_alignable)
  out$alignments_per_read <- ratio(out$n_hits, out$n_alignable)
  rownames(out) <- NULL
  message(sprintf(
    "rsem_mapping_qc: %d samples; alignable fraction %s to %s; unique-gene fraction of alignable reads %s to %s",
    nrow(out), format(signif(min(out$alignable_fraction, na.rm = TRUE), 3)), format(signif(max(out$alignable_fraction, na.rm = TRUE), 3)),
    format(signif(min(out$unique_gene_fraction, na.rm = TRUE), 3)), format(signif(max(out$unique_gene_fraction, na.rm = TRUE), 3))
  ))
  out
}

# First three lines of a .cnt file (comments after '#' are ignored).
.read_cnt <- function(path, sample_id, ref) {
  label <- sprintf("%s (sample %s)", ref, sample_id)
  lines <- readLines(path, n = 3L, warn = FALSE)
  nums <- lapply(sub("#.*$", "", lines), function(l) suppressWarnings(as.numeric(strsplit(trimws(l), "[[:space:]]+")[[1L]])))
  ok <- length(nums) == 3L && length(nums[[1L]]) == 4L && length(nums[[2L]]) >= 3L && length(nums[[3L]]) >= 2L &&
    !anyNA(unlist(nums[1:3])) && all(unlist(nums[1:3]) >= 0)
  if (!ok) .stopf("%s is not an RSEM .cnt file: expected the lines 'N0 N1 N2 N_tot', 'nUnique nMulti nUncertain' and 'nHits read_type'.", label)
  n <- nums[[1L]]; u <- nums[[2L]]; h <- nums[[3L]]
  if (abs(n[1L] + n[2L] + n[3L] - n[4L]) > 0.5) .stopf("%s: N_tot (%s) is not N0 + N1 + N2 (%s).", label, format(n[4L]), format(n[1L] + n[2L] + n[3L]))
  if (abs(u[1L] + u[2L] - n[2L]) > 0.5) .stopf("%s: nUnique + nMulti (%s) is not N1 (%s).", label, format(u[1L] + u[2L]), format(n[2L]))
  data.frame(
    sample_id = sample_id, n_total = n[4L], n_unalignable = n[1L], n_alignable = n[2L], n_filtered = n[3L],
    n_unique = u[1L], n_multi = u[2L], n_uncertain = u[3L], n_hits = h[1L], read_type = h[2L],
    stringsAsFactors = FALSE
  )
}
