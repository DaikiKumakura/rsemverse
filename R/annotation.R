# Gene annotation from a GTF file (or a data frame), used for gene types that
# the OrgDb does not separate: rRNA versus ribosomal protein genes, snoRNA,
# snRNA, miRNA, lncRNA, mitochondrial rRNA and tRNA.

# Gene groups that need a GTF annotation, and the gene types (GENCODE
# `gene_type`, Ensembl and NCBI `gene_biotype`) that make them up.
.gtf_groups <- list(
  rrna = c("rRNA", "rRNA_pseudogene"),
  mito_rrna = "Mt_rRNA",
  mito_trna = "Mt_tRNA",
  snorna = "snoRNA",
  snrna = "snRNA",
  mirna = "miRNA",
  lncrna = c("lncRNA", "lincRNA")
)

# Gene rows of a GTF as a data frame (gene_id without version, gene_name,
# gene_type). `annotation` may also be a data frame with these columns.
.read_annotation <- function(annotation, need_type = TRUE) {
  if (is.data.frame(annotation)) {
    a <- as.data.frame(annotation, stringsAsFactors = FALSE)
    need <- c("gene_id", if (need_type) "gene_type")
    if (!all(need %in% names(a))) .stopf("An annotation data frame needs the columns gene_id and %s (gene_name is optional).", if (need_type) "gene_type" else "gene_name")
    if (!"gene_name" %in% names(a)) a$gene_name <- NA_character_
    if (!"gene_type" %in% names(a)) a$gene_type <- NA_character_
    a <- a[, c("gene_id", "gene_name", "gene_type")]
  } else if (is.character(annotation) && length(annotation) == 1L && !is.na(annotation)) {
    if (!file.exists(annotation) || dir.exists(annotation)) .stopf("`annotation` must be an existing GTF file: %s", annotation)
    con <- if (grepl("\\.gz$", annotation)) gzfile(annotation, "r") else file(annotation, "r")
    on.exit(close(con))
    # The gene rows give the names and types; a gene without a gene row, or without a name or type
    # there, gets them from its transcript rows and then its exon rows (RSEM needs only exon rows,
    # so a GTF made for it may have no gene rows, and a GTF may have them for some genes only).
    get <- function(attrs, keys) {
      out <- rep(NA_character_, length(attrs))
      for (k in keys) {
        hit <- is.na(out) & grepl(paste0("(^|; ?)", k, " \""), attrs)
        out[hit] <- sub(paste0('^.*?(^|; ?)', k, ' "([^"]*)".*$'), "\\2", attrs[hit], perl = TRUE)
      }
      out
    }
    first_known <- function(df) {  # df ordered by priority; the first non-missing value per gene_id
      ids <- unique(df$gene_id)
      out <- data.frame(gene_id = ids, stringsAsFactors = FALSE)
      for (f in c("gene_name", "gene_type")) {
        d <- df[!is.na(df[[f]]), , drop = FALSE]
        d <- d[!duplicated(d$gene_id), , drop = FALSE]
        out[[f]] <- d[[f]][match(ids, d$gene_id)]
      }
      out
    }
    parts <- list()
    done <- character(0)
    repeat {
      chunk <- readLines(con, n = 500000L, warn = FALSE)
      if (!length(chunk)) break
      # Fixed-string searches first: they are much faster than a regular expression on every row.
      is_gene <- grepl("\tgene\t", chunk, fixed = TRUE)
      is_tx <- !is_gene & grepl("\ttranscript\t", chunk, fixed = TRUE)
      # every other feature row (exon, CDS, UTR, ...) is used like an exon row
      use <- which(is_gene | is_tx | (nzchar(chunk) & substr(chunk, 1L, 1L) != "#"))
      if (!length(use)) next
      rows <- chunk[use]
      feat <- ifelse(is_gene[use], "gene", ifelse(is_tx[use], "transcript", "exon"))
      parse_rows <- function(r, f) {
        attrs <- sub("^([^\t]*\t){8}", "", r)
        df <- data.frame(
          gene_id = get(attrs, "gene_id"), gene_name = get(attrs, c("gene_name", "gene")),
          gene_type = get(attrs, c("gene_type", "gene_biotype")), prio = match(f, c("gene", "transcript", "exon")),
          stringsAsFactors = FALSE
        )
        df <- df[!is.na(df$gene_id), , drop = FALSE]
        first_known(df[order(df$prio), , drop = FALSE])
      }
      # The gene rows first; the transcript and exon rows of a gene that is already complete are not parsed.
      g <- feat == "gene"
      if (any(g)) {
        kn <- parse_rows(rows[g], feat[g])
        parts[[length(parts) + 1L]] <- kn
        done <- union(done, kn$gene_id[!is.na(kn$gene_name) & !is.na(kn$gene_type)])
      }
      rows <- rows[!g]; feat <- feat[!g]
      if (length(rows) && length(done)) {
        p <- regexpr("gene_id \"[^\"]*\"", rows, perl = TRUE)
        id <- substr(rows, p + 9L, p + attr(p, "match.length") - 2L)
        rows <- rows[!(id %in% done)]; feat <- feat[!(id %in% done)]
      }
      if (length(rows)) {
        kn <- parse_rows(rows, feat)
        parts[[length(parts) + 1L]] <- kn
        done <- union(done, kn$gene_id[!is.na(kn$gene_name) & !is.na(kn$gene_type)])
      }
    }
    if (!length(parts)) .stopf("No 'gene', 'transcript' or 'exon' rows found in the GTF file %s.", annotation)
    a <- first_known(do.call(rbind, parts))
    if (need_type && all(is.na(a$gene_type))) .stopf("The GTF file has no gene_type or gene_biotype attribute.")
  } else {
    .stopf("`annotation` must be a path to a GTF file or a data frame.")
  }
  a$gene_id <- .strip_ensembl_version(a$gene_id)
  a[!is.na(a$gene_id), , drop = FALSE]
}

# Gene type for each row: through the original gene IDs behind the row when
# `ids` (attr "gene_ids" of rsem_read()) is given, else Ensembl IDs by ID; the
# remaining names by gene_name. A row with several types gets all of them
# joined by ";".
.annotate_rows <- function(genes, ann, ids = NULL) {
  stable <- .strip_ensembl_version(genes)
  join <- function(v) paste(sort(unique(v[!is.na(v)])), collapse = ";")
  by_id <- tapply(ann$gene_type, ann$gene_id, join)
  ok_name <- !is.na(ann$gene_name)
  by_name <- tapply(ann$gene_type[ok_name], ann$gene_name[ok_name], join)
  type <- unname(by_id[stable])
  if (!is.null(ids)) {
    if (!is.data.frame(ids) || !all(c("gene", "gene_id") %in% names(ids))) {
      .stopf("`ids` must be the attribute \"gene_ids\" of a table from rsem_read() (columns gene and gene_id).")
    }
    t_id <- by_id[.strip_ensembl_version(ids$gene_id)]
    hit <- !is.na(t_id)
    via <- tapply(unname(t_id[hit]), ids$gene[hit], function(v) join(unlist(strsplit(v, ";", fixed = TRUE))))
    type <- ifelse(is.na(type), unname(via[genes]), type)
  }
  miss <- is.na(type)
  type[miss] <- unname(by_name[genes[miss]])
  type[!is.na(type) & !nzchar(type)] <- NA_character_
  type
}

#' Share of each gene type in every sample
#'
#' Sums the values of each gene type (from a GTF file) per sample and divides
#' by the sample total, for example to see how much of the TPM goes to rRNA,
#' mitochondrial genes or lncRNA. Genes that the annotation does not know are
#' counted as `"unannotated"`, and genes whose name has several types as
#' `"multiple"`.
#'
#' A type that the annotation does not contain cannot be measured: if a
#' requested type has no gene in the annotation, its share is `NA`, not 0, and
#' a message says so. For example, an annotation without rRNA genes says
#' nothing about rRNA contamination. Reads from genes outside the annotation
#' were never quantified by RSEM and are not in any share.
#'
#' @param x The gene by sample table (TPM or counts), as returned by
#'   [rsem_read()].
#' @param annotation Path to a GTF file (plain or `.gz`; GENCODE, Ensembl or
#'   NCBI, with `gene_type` or `gene_biotype`), or a data frame with the
#'   columns `gene_id`, `gene_type` and optionally `gene_name`. It should be
#'   the annotation the RSEM reference was built from.
#' @param ids The attribute `"gene_ids"` of the table [rsem_read()] returned
#'   (default: that attribute of `x`, when present). Rows are then matched to
#'   the annotation through their original gene IDs, which avoids losing
#'   genes whose OrgDb symbol differs from the GTF gene name. Without it,
#'   rows are matched by Ensembl ID or by gene name.
#' @param total Optional total per sample used as the denominator of `share`:
#'   a numeric vector in the column order of `x`, or named by sample. The
#'   default is the sum of `x`, so shares are relative to the genes that are in
#'   `x`: after [rsem_filter()] or [rsem_subset()] the denominator is the
#'   total of the remaining genes. For QC shares (rRNA, mitochondrial) pass the
#'   totals of the unfiltered table, `colSums(unfiltered)`.
#' @param types Gene types to report (values of `gene_type`, such as
#'   `"rRNA"` or `"protein_coding"`). Default: every type found in `x`.
#'
#' @return A data frame with one row per sample and type: `sample_id`,
#'   `type`, `genes` (genes of that type in `x`), `value` (sum) and `share`
#'   (value / sample total). A sample whose total is 0 has `NA` shares.
#'
#' @examples
#' ann <- data.frame(gene_id = c("ENSG01", "ENSG02", "ENSG03"),
#'                   gene_name = c("A", "B", "C"),
#'                   gene_type = c("protein_coding", "rRNA", "lncRNA"))
#' x <- data.frame(S1 = c(10, 5, 1), S2 = c(20, 0, 2), row.names = c("A", "B", "C"))
#' rsem_biotype_profile(x, ann)
#' @export
rsem_biotype_profile <- function(x, annotation, types = NULL, ids = attr(x, "gene_ids", exact = TRUE),
                                 total = NULL) {
  m <- .expr_matrix(x)
  total_arg <- total
  ann <- .read_annotation(annotation)
  type <- .annotate_rows(rownames(m), ann, ids)
  type[is.na(type)] <- "unannotated"
  type[grepl(";", type)] <- "multiple"
  if (is.null(types)) {
    types <- sort(unique(type), method = "radix")
  } else if (!is.character(types) || !length(types) || anyNA(types)) {
    .stopf("`types` must be one or more gene types.")
  }
  absent <- setdiff(types, c(ann$gene_type, "unannotated", "multiple"))
  if (is.null(total)) {
    total <- colSums(m)
  } else {
    if (!is.numeric(total) || length(total) != ncol(m) || anyNA(total) || any(total < 0)) {
      .stopf("`total` needs one number per sample (%d), not negative.", ncol(m))
    }
    if (!is.null(names(total))) {
      if (!setequal(names(total), colnames(m))) .stopf("The names of `total` must be the sample names of `x`.")
      total <- total[colnames(m)]
    }
  }
  total <- as.numeric(total)
  zero <- total == 0
  rows <- lapply(types, function(t) {
    sel <- type == t
    v <- colSums(m[sel, , drop = FALSE])
    data.frame(
      sample_id = colnames(m), type = t, genes = sum(sel), value = if (t %in% absent) NA_real_ else v,
      share = if (t %in% absent) NA_real_ else ifelse(zero, NA_real_, v / total), stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  message(sprintf(
    "rsem_biotype_profile: %s genes x %d samples; shares of the %s; %s of the genes not in the annotation%s%s",
    .fmt(nrow(m)), ncol(m), if (is.null(total_arg)) "total of the genes in `x`" else "totals given in `total`",
    .fmt(sum(type == "unannotated")),
    if (any(zero)) paste0("; total 0, share NA: ", paste(colnames(m)[zero], collapse = ", ")) else "",
    if (length(absent)) paste0("; not in the annotation, so not measurable (NA): ", paste(absent, collapse = ", ")) else ""
  ))
  out
}
