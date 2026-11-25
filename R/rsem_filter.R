#' Remove lowly expressed or uninformative genes
#'
#' Keeps the genes that pass one of four common filters that need no sample
#' groups. The table keeps its values and row order; only rows are removed.
#' A message reports the rule and how many genes were kept.
#'
#' @section Methods:
#' * `"edger"` (counts only; the default for counts): the rule of
#'   `edgeR::filterByExpr()` with all samples treated as one group. A gene is
#'   kept when its count per million (CPM) reaches `min_value / median library
#'   size * 1e6` in at least *k* samples and its total count is at least 15.
#'   *k* is the number of samples, or `10 + (n - 10) * 0.7` when there are more
#'   than 10. `min_value` is the count a gene needs in a sample of median size
#'   (default 10, the edgeR default). With `group`, *k* is the size of the
#'   smallest group (then reduced in the same way when it is above 10), and
#'   with `design` it is `1 / max(hat(design))`, as in `filterByExpr()`. The
#'   groups only set how many samples a gene must be expressed in; genes are
#'   not chosen by comparing groups, so the filter stays independent of the
#'   later test. Without them, a study with 3 controls and 10 treated samples
#'   asks for expression in all 13, which can remove genes that only the
#'   smaller group expresses.
#' * `"detected"` (the default for TPM): a gene is kept when its value is at
#'   least `min_value` in at least `min_samples` samples. The defaults are
#'   TPM >= 1 or count >= 10, in at least 20% of the samples (at least one
#'   sample). This is the form used for bulk RNA-seq (for example, GTEx keeps
#'   genes above a TPM and a count threshold in at least 20% of samples) and
#'   for single-cell data ("expressed in at least n cells").
#' * `"nonzero"`: removes genes that are zero in every sample.
#' * `"variable"`: keeps the `n_top` genes with the largest variance of
#'   `log2(value + 1)` across samples (default 2000), as used before PCA,
#'   clustering or heatmaps. Genes that are zero everywhere are never kept.
#'   Rows stay in their original order, not in the order of variance.
#'
#' Filters can be chained, for example `"edger"` and then `"variable"`.
#'
#' @param x The gene by sample table, as returned by [rsem_read()]: gene
#'   names as row names (or in the first column) and one numeric column per
#'   sample.
#' @param expr `"count"` or `"tpm"`: what the values are. Required.
#' @param method `"edger"`, `"detected"`, `"nonzero"` or `"variable"`. The
#'   default is `"edger"` for counts and `"detected"` for TPM.
#' @param min_value For `"edger"`, the count needed in a sample of median size
#'   (default 10). For `"detected"`, the value a gene must reach (default 1
#'   for TPM, 10 for counts).
#' @param min_samples For `"detected"`: in how many samples the gene must
#'   reach `min_value`. A number below 1 is read as a proportion of the
#'   samples (`0.5` = half). Default 0.2.
#' @param n_top For `"variable"`: how many genes to keep. Default 2000.
#' @param group For `"edger"`: the group of each sample (a vector in the
#'   column order of `x`, or named by sample), as in `edgeR::filterByExpr()`.
#' @param design For `"edger"`: a design matrix with one row per sample (in
#'   the column order of `x`) instead of `group`.
#'
#' @return A data frame with the kept genes, in the same layout as `x`
#'   (genes as row names, samples as columns).
#'
#' @references
#' Chen Y, Lun ATL, Smyth GK (2016). From reads to genes to pathways:
#' differential expression analysis of RNA-Seq experiments using Rsubread
#' and the edgeR quasi-likelihood pipeline. *F1000Research* 5:1438.
#'
#' GTEx Consortium (2020). The GTEx Consortium atlas of genetic regulatory
#' effects across human tissues. *Science* 369:1318-1330.
#'
#' @examples
#' if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   ex <- system.file("extdata", "example", package = "rsemverse")
#'   counts <- rsem_read(file.path(ex, "data"), file.path(ex, "metadata.csv"),
#'                       species = "human", expr = "count")
#'   kept <- rsem_filter(counts, expr = "count")                 # edgeR rule
#'   tpm <- rsem_read(file.path(ex, "data"), file.path(ex, "metadata.csv"),
#'                    species = "human", expr = "tpm")
#'   rsem_filter(tpm, expr = "tpm", min_value = 1, min_samples = 2)
#'   rsem_filter(tpm, expr = "tpm", method = "variable", n_top = 10)
#' }
#' @export
rsem_filter <- function(x, expr = attr(x, "expr", exact = TRUE), method = NULL, min_value = NULL, min_samples = 0.2, n_top = 2000,
                        group = NULL, design = NULL) {
  # a normalized or log-transformed table has no expr; "variable" and "nonzero" do not need one
  if (is.null(expr) && is.character(method) && length(method) == 1L && method %in% c("variable", "nonzero")) expr <- "tpm"
  expr <- .need_expr(expr)
  expr <- .choice(expr, c("count", "tpm"), "expr")
  if (is.null(method)) method <- if (expr == "count") "edger" else "detected"
  method <- .choice(method, c("edger", "detected", "nonzero", "variable"), "method")
  m <- .expr_matrix(x)
  n <- ncol(m)
  if ((!is.null(group) || !is.null(design)) && method != "edger") {
    .stopf("`group` and `design` are used by method = \"edger\" only.")
  }
  if (!is.null(group) && !is.null(design)) .stopf("Give `group` or `design`, not both.")

  if (method == "edger") {
    if (expr != "count") {
      .stopf("method = \"edger\" needs counts (expr = \"count\"); for TPM use \"detected\".")
    }
    min_count <- if (is.null(min_value)) 10 else .number(min_value, "min_value", lower = 0)
    lib <- colSums(m)
    if (any(lib <= 0)) .stopf("Every sample needs a nonzero total count; zero: %s", paste(colnames(m)[lib <= 0], collapse = ", "))
    k <- n
    how <- sprintf("%d samples", n)
    if (!is.null(group)) {
      if (!is.null(names(group))) {
        if (!setequal(names(group), colnames(m)) || anyDuplicated(names(group))) .stopf("The names of `group` must be the sample names of `x`.")
        group <- group[colnames(m)]
      }
      if (length(group) != n || anyNA(group)) .stopf("`group` needs one value per sample (%d) and no missing values.", n)
      sizes <- tabulate(as.factor(group))
      k <- min(sizes[sizes > 0L])
      how <- sprintf("the smallest of %d groups, %d samples", sum(sizes > 0L), k)
    } else if (!is.null(design)) {
      design <- as.matrix(design)
      if (nrow(design) != n || anyNA(design)) .stopf("`design` needs one row per sample (%d) and no missing values.", n)
      k <- 1 / max(stats::hat(design))
      how <- sprintf("1 / max leverage of the design = %s", format(signif(k, 4)))
    }
    if (k > 10) k <- 10 + (k - 10) * 0.7
    cutoff <- min_count / stats::median(lib) * 1e6
    cpm <- sweep(m, 2L, lib, "/") * 1e6
    tol <- 1e-14
    keep <- rowSums(cpm >= cutoff) >= k - tol & rowSums(m) >= 15 - tol
    rule <- sprintf("CPM >= %s in >= %s of %d samples%s and total count >= 15", signif(cutoff, 3), format(signif(k, 6)), n,
                    if (is.null(group) && is.null(design)) "" else sprintf(" (%s)", how))
  } else if (method == "detected") {
    threshold <- if (is.null(min_value)) c(count = 10, tpm = 1)[[expr]] else .number(min_value, "min_value", lower = 0)
    needed <- .needed_samples(min_samples, n)
    keep <- rowSums(m >= threshold) >= needed
    rule <- sprintf("%s >= %s in >= %s of %d samples", if (expr == "tpm") "TPM" else "count", format(threshold), format(needed), n)
  } else if (method == "nonzero") {
    keep <- rowSums(m > 0) > 0
    rule <- "nonzero in at least one sample"
  } else {
    n_top <- .number(n_top, "n_top", lower = 1, integer = TRUE)
    if (n < 2L) .stopf("method = \"variable\" needs at least two samples.")
    lm <- log2(m + 1)
    v <- apply(lm, 1L, stats::var)
    v[rowSums(m > 0) == 0] <- -Inf
    candidates <- which(is.finite(v))
    chosen <- candidates[order(-v[candidates])][seq_len(min(n_top, length(candidates)))]
    keep <- seq_len(nrow(m)) %in% chosen
    rule <- sprintf("the %s genes with the largest variance of log2(%s + 1)", .fmt(sum(keep)), expr)
  }

  message(sprintf(
    "rsem_filter (%s): kept %s of %s genes (%s); removed %s",
    method, .fmt(sum(keep)), .fmt(nrow(m)), rule, .fmt(sum(!keep))
  ))
  .as_table(m[keep, , drop = FALSE], from = x, step = sprintf("filter: %s, kept %d of %d genes", method, sum(keep), nrow(m)))
}
