#' The fragment length distribution and read start positions learned by RSEM
#'
#' Reads the `*.model` file that `rsem-calculate-expression` writes (in
#' `<name>.stat/<name>.model`) for every sample, and returns what RSEM learned
#' about the library: the fragment length distribution (or the read length
#' distribution of single-end reads) and, when it was estimated, the read start
#' position distribution (RSPD), as numbers that can be compared between
#' samples. The plots of one sample are made by RSEM's own `rsem-plot-model`;
#' this function does not draw.
#'
#' @section What is read:
#' Following RSEM's `model_file_description.txt`: the model type (0 single-end
#' without quality scores, 1 single-end with, 2 paired-end without, 3
#' paired-end with), the probability of a read from the forward strand, the
#' fragment length distribution (`(lower, upper]`, one probability per length),
#' the optional read length distribution of single-end models, and the read
#' start position distribution (bins of a transcript). The sections after
#' these (quality scores, sequencing error) are not read. In the files RSEM
#' 1.3.3 wrote for the test data of this package, the number of bins and their
#' probabilities follow the flag only when the RSPD was estimated
#' (`--estimate-rspd`); without it the flag is followed by the next section.
#'
#' `forward_prob` is the strand probability that the run was given (`--strandedness`,
#' or 0.5 for none), written into the model; it is not learned from the reads.
#'
#' Without `--estimate-rspd` RSEM assumes that reads start uniformly along a
#' transcript, so a flat RSPD is not evidence that the library has no 3'
#' bias; `rspd_estimated` says which case it is. With single-end reads and no
#' `--fragment-length-mean`, the distribution is that of the reads, not of the
#' fragments.
#'
#' @inheritParams rsem_mapping_qc
#' @return A list of data frames:
#'   * `sample_summary`: per sample `model_type`, `forward_prob`,
#'     `length_lower` and `length_upper` (bounds of the distribution),
#'     `length_mean`, `length_sd` and `length_mode` of the fragment (or read)
#'     lengths, `rspd_estimated` and `rspd_bins`;
#'   * `length_distribution`: `sample_id`, `length`, `probability`;
#'   * `rspd`: `sample_id`, `bin` and `probability`, only for the samples in
#'     which it was estimated.
#' @examples
#' f <- tempfile(fileext = ".model")
#' writeLines(c("0", "", "0.5", "", "9 13 4", "0.1 0.4 0.4 0.1", "", "0", "", "1", "2", "0.25 0.75"), f)
#' rsem_model_qc(dirname(f), metadata = NULL)$sample_summary
#' @export
rsem_model_qc <- function(data = "data", metadata = NULL, sample_from = "filename") {
  samples <- if (is.null(metadata)) {
    .discover_samples(data, "model", .choice(sample_from, c("filename", "parent_dir"), "sample_from"))
  } else {
    .read_samples(data, metadata)
  }
  parsed <- lapply(seq_len(nrow(samples)), function(i) .read_model(samples$path[[i]], samples$sample_id[[i]], samples$ref[[i]]))
  summary <- do.call(rbind, lapply(parsed, `[[`, "summary"))
  length_distribution <- do.call(rbind, lapply(parsed, `[[`, "length"))
  rspd_tabs <- Filter(Negate(is.null), lapply(parsed, `[[`, "rspd"))
  rspd <- if (length(rspd_tabs)) do.call(rbind, rspd_tabs) else data.frame(sample_id = character(0), bin = integer(0), probability = numeric(0))
  rownames(summary) <- rownames(length_distribution) <- rownames(rspd) <- NULL
  message(sprintf(
    "rsem_model_qc: %d samples; mean fragment length %s to %s; read start positions estimated for %d",
    nrow(summary), format(signif(min(summary$length_mean), 4)), format(signif(max(summary$length_mean), 4)), sum(summary$rspd_estimated)
  ))
  list(sample_summary = summary, length_distribution = length_distribution, rspd = rspd)
}

# The first sections of a .model file. Blank lines separate sections, the numbers are read in order.
.read_model <- function(path, sample_id, ref) {
  label <- sprintf("%s (sample %s)", ref, sample_id)
  tok <- scan(path, what = "", quiet = TRUE)
  pos <- 0L
  take <- function(n, what) {
    if (pos + n > length(tok)) .stopf("%s is not an RSEM .model file: it ends inside %s.", label, what)
    out <- suppressWarnings(as.numeric(tok[pos + seq_len(n)]))
    if (anyNA(out) || any(!is.finite(out))) .stopf("%s is not an RSEM .model file: %s has values that are not numbers.", label, what)
    pos <<- pos + n
    out
  }
  dist <- function(what) {
    b <- take(3L, paste(what, "bounds"))
    if (b[3L] != b[2L] - b[1L] || b[3L] < 1) .stopf("%s: the span of the %s (%s) is not upper - lower.", label, what, format(b[3L]))
    p <- take(b[3L], paste("the", what))
    if (any(p < 0) || any(p > 1)) .stopf("%s: the %s has a probability outside 0 to 1.", label, what)
    if (abs(sum(p) - 1) > 1e-3) .stopf("%s: the probabilities of the %s do not add up to 1.", label, what)
    list(lower = b[1L], upper = b[2L], prob = p, len = b[1L] + seq_len(b[3L]))
  }
  model_type <- take(1L, "the model type")
  if (!model_type %in% 0:3) .stopf("%s: the model type (%s) is not 0, 1, 2 or 3.", label, format(model_type))
  forward <- take(1L, "the forward probability")
  if (forward < 0 || forward > 1) .stopf("%s: the forward probability (%s) is outside 0 to 1.", label, format(forward))
  fl <- dist("fragment length distribution")
  if (model_type < 2) {
    has_optional <- take(1L, "the optional read length flag")
    if (has_optional == 1) dist("read length distribution")
  }
  est <- take(1L, "the RSPD flag")
  if (!est %in% 0:1) .stopf("%s: the RSPD flag (%s) is not 0 or 1.", label, format(est))
  rspd <- NULL
  bins <- NA_integer_
  if (est == 1) {
    bins <- take(1L, "the number of RSPD bins")
    if (bins < 1 || bins != round(bins)) .stopf("%s: the number of RSPD bins (%s) is not a positive integer.", label, format(bins))
    p <- take(bins, "the RSPD")
    if (any(p < 0) || any(p > 1)) .stopf("%s: the RSPD has a probability outside 0 to 1.", label)
    if (abs(sum(p) - 1) > 1e-3) .stopf("%s: the RSPD probabilities do not add up to 1.", label)
    rspd <- data.frame(sample_id = sample_id, bin = seq_len(bins), probability = p, stringsAsFactors = FALSE)
  }
  mu <- sum(fl$len * fl$prob)
  list(
    summary = data.frame(
      sample_id = sample_id, model_type = model_type, forward_prob = forward, length_lower = fl$lower, length_upper = fl$upper,
      length_mean = mu, length_sd = sqrt(sum(fl$prob * (fl$len - mu)^2)), length_mode = fl$len[[which.max(fl$prob)]],
      rspd_estimated = est == 1, rspd_bins = as.integer(bins), stringsAsFactors = FALSE
    ),
    length = data.frame(sample_id = sample_id, length = fl$len, probability = fl$prob, stringsAsFactors = FALSE),
    rspd = rspd
  )
}

#' The share of reads that RSEM explains as background noise
#'
#' Reads the `*.theta` file that `rsem-calculate-expression` writes (in
#' `<name>.stat/<name>.theta`) and returns the first value of the estimated
#' read-generating probabilities: the share of reads from the noise isoform that
#' RSEM adds to its model. The file has the number of values on its first
#' line, the probabilities before a correction on the second line, and after
#' it on the third line (the correction takes the reads of transcripts with an
#' effective length of 0 out); the third line is the one RSEM's simulator
#' uses as the noise share.
#'
#' The noise share is not the share of unalignable reads of
#' [rsem_mapping_qc()]: in the toy data of this package, where the only
#' reads that align nowhere are the added ones, the second-line value equals
#' `N0 / N_tot` and the third-line value `N0 / (N_tot - reads of the transcript
#' with effective length 0)`, but in real data noise also takes alignable reads
#' that no transcript explains. It depends on the model and the reference, and
#' is for comparing samples of one data set, not a quality threshold.
#'
#' @inheritParams rsem_mapping_qc
#' @return A data frame with one row per sample: `sample_id`,
#'   `n_components` (the number of values per line, the noise and the
#'   transcripts), `noise_fraction` (third line) and `noise_fraction_uncorrected`
#'   (second line).
#' @examples
#' f <- tempfile(fileext = ".theta")
#' writeLines(c("3", "0.1 0.5 0.4", "0.125 0.5 0.375"), f)
#' rsem_noise_qc(dirname(f), metadata = NULL)
#' @export
rsem_noise_qc <- function(data = "data", metadata = NULL, sample_from = "filename") {
  samples <- if (is.null(metadata)) {
    .discover_samples(data, "theta", .choice(sample_from, c("filename", "parent_dir"), "sample_from"))
  } else {
    .read_samples(data, metadata)
  }
  rows <- lapply(seq_len(nrow(samples)), function(i) {
    label <- sprintf("%s (sample %s)", samples$ref[[i]], samples$sample_id[[i]])
    lines <- readLines(samples$path[[i]], n = 3L, warn = FALSE)
    num <- lapply(lines, function(l) suppressWarnings(as.numeric(strsplit(trimws(l), "[[:space:]]+")[[1L]])))
    n <- if (length(num) >= 1L && length(num[[1L]]) == 1L) num[[1L]] else NA_real_
    ok <- length(num) == 3L && !is.na(n) && n >= 1 && n == round(n) && length(num[[2L]]) == n && length(num[[3L]]) == n &&
      !anyNA(unlist(num)) && all(unlist(num[2:3]) >= 0 & unlist(num[2:3]) <= 1)
    if (!ok) .stopf("%s is not an RSEM .theta file: expected the number of values and two lines of that many probabilities.", label)
    if (abs(sum(num[[2L]]) - 1) > 1e-3 || abs(sum(num[[3L]]) - 1) > 1e-3) .stopf("%s: the probabilities do not add up to 1.", label)
    data.frame(
      sample_id = samples$sample_id[[i]], n_components = as.integer(n), noise_fraction = num[[3L]][[1L]],
      noise_fraction_uncorrected = num[[2L]][[1L]], stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  message(sprintf("rsem_noise_qc: %d samples; noise fraction %s to %s", nrow(out), format(signif(min(out$noise_fraction), 3)), format(signif(max(out$noise_fraction), 3))))
  out
}
