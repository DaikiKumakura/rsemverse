#' Keep or remove groups of genes
#'
#' Selects genes by built-in gene groups, so that no gene list has to be
#' prepared. Gene types and chromosomes come from the same annotation package
#' as [rsem_read()] (`org.Hs.eg.db` or `org.Mm.eg.db`); the other groups are
#' fixed symbol patterns. Values are not changed. A message reports how many
#' genes each group matched.
#'
#' @section Gene groups:
#' * `"protein_coding"`: genes whose type is `protein-coding` in the
#'   annotation package (NCBI Gene). Many analyses use only these genes, for
#'   example to compare with proteomics, to reduce noise from lowly expressed
#'   non-coding RNAs and pseudogenes, or because a signature was built on them.
#' * `"mito"`: the 37 genes of the mitochondrial genome (chromosome MT;
#'   symbols `MT-` in human and `mt-` in mouse). The rRNA and tRNA genes have
#'   no Ensembl link in the OrgDb, so [rsem_read()] keeps them under their
#'   Ensembl IDs; a built-in list of the 37 Ensembl IDs (Ensembl release 116)
#'   finds them as well, together with the symbols the OrgDb gives to these
#'   IDs (the ID of `MT-TS1` maps to `MIR12136`, a microRNA inside it). A high share of these reads is a sign
#'   of damaged cells or degraded RNA.
#' * `"ribo"`: cytosolic ribosomal protein genes (`RPL*`, `RPS*`, `RPLP0-2`,
#'   `RPSA`; for mouse `Rpl*`, `Rps*`), protein-coding only. Mitochondrial
#'   ribosomal genes (`MRPL*`, `MRPS*`), kinases such as `RPS6KA1` and
#'   pseudogenes such as `RPL10P1` are not included. They are highly expressed and often removed before clustering.
#' * `"hb"`: hemoglobin genes (`HBA1`, `HBA2`, `HBB`, `HBD`, `HBE1`, `HBG1`,
#'   `HBG2`, `HBM`, `HBQ1`, `HBZ`; for mouse `Hba-*`, `Hbb-*`, `Hbq1*`).
#'   protein-coding only. They signal blood contamination. `HBEGF`,
#'   hemoglobin pseudogenes and the locus control regions are not included.
#' * `"sex"`: genes on chromosome Y and `XIST`. They separate samples by sex
#'   and are often removed when sex is not of interest. Pseudoautosomal genes
#'   (on both X and Y, expressed in both sexes) are not included.
#' * `"ercc"`: ERCC spike-in controls (IDs starting with `ERCC-`).
#'
#' `keep` is applied first (a gene is kept when it is in any of the listed
#' groups), then `remove`. Genes are matched by symbol; rows that kept their
#' Ensembl ID in [rsem_read()] are matched by that ID. Genes that the
#' annotation package does not know are not in `"protein_coding"`,
#' `"ribo"`, `"hb"` or `"sex"`; for example, chromosome Y genes that
#' [rsem_read()] kept under their Ensembl ID are not found by `"sex"`.
#'
#' @inheritParams rsem_filter
#' @param species `"human"` or `"mouse"`. Required.
#' @param keep,remove Gene groups: one or more of `"protein_coding"`,
#'   `"mito"`, `"ribo"`, `"hb"`, `"sex"` and `"ercc"`. At least one of the
#'   two is required.
#'
#' @return A data frame with the selected genes, in the same layout and order
#'   as `x`.
#'
#' @examples
#' if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   ex <- system.file("extdata", "example", package = "rsemverse")
#'   tpm <- rsem_read(file.path(ex, "data"), file.path(ex, "metadata.csv"),
#'                    species = "human", expr = "tpm")
#'   pc <- rsem_subset(tpm, "human", keep = "protein_coding")
#'   rsem_subset(tpm, "human", keep = "protein_coding", remove = c("mito", "ribo"))
#' }
#' @export
rsem_subset <- function(x, species, keep = NULL, remove = NULL) {
  if (missing(species)) .stopf("`species` must be \"human\" or \"mouse\".")
  species <- .choice(species, c("human", "mouse"), "species")
  groups <- c("protein_coding", "mito", "ribo", "hb", "sex", "ercc")
  check <- function(g, name) {
    if (is.null(g)) return(character(0))
    if (!is.character(g) || !length(g) || anyNA(g) || !all(g %in% groups)) {
      .stopf("`%s` must be one or more of: %s.", name, paste(sprintf("\"%s\"", groups), collapse = ", "))
    }
    unique(g)
  }
  keep <- check(keep, "keep")
  remove <- check(remove, "remove")
  if (!length(keep) && !length(remove)) .stopf("Give at least one gene group in `keep` or `remove`.")
  m <- .expr_matrix(x, nonnegative = FALSE)

  member <- .gene_groups(rownames(m), species, union(keep, remove))
  report <- vapply(union(keep, remove), function(g) sprintf("%s %s", g, .fmt(sum(member[, g]))), character(1))

  selected <- if (length(keep)) rowSums(member[, keep, drop = FALSE]) > 0 else rep(TRUE, nrow(m))
  if (length(remove)) selected <- selected & rowSums(member[, remove, drop = FALSE]) == 0

  message(sprintf(
    "rsem_subset: kept %s of %s genes%s%s (genes per group: %s)",
    .fmt(sum(selected)), .fmt(nrow(m)),
    if (length(keep)) paste0("; keep ", paste(keep, collapse = " or ")) else "",
    if (length(remove)) paste0("; remove ", paste(remove, collapse = ", ")) else "",
    paste(report, collapse = ", ")
  ))
  .as_table(m[selected, , drop = FALSE])
}

# Logical matrix: genes x requested groups.
.gene_groups <- function(genes, species, groups) {
  out <- matrix(FALSE, length(genes), length(groups), dimnames = list(genes, groups))
  need_db <- intersect(groups, c("protein_coding", "mito", "ribo", "hb", "sex"))
  ann <- if (length(need_db)) .gene_annotation(genes, species, need_db) else NULL

  coding <- !is.null(ann$type) & !is.na(ann$type) & ann$type == "protein-coding"
  for (g in groups) {
    out[, g] <- switch(g,
      protein_coding = coding,
      mito = grepl("^MT-", genes, ignore.case = TRUE) | (!is.na(ann$chr) & ann$chr == "MT") |
        .strip_ensembl_version(genes) %in% .mito_ensembl[[species]] | genes %in% .mito_symbols(species),
      # Pseudogenes (RPL10P1, Hba-ps3) and regulatory regions (HBB-LCR) share
      # these prefixes, so both groups are limited to protein-coding genes.
      ribo = coding & grepl("^RP[SL][0-9]+[A-Z]?[0-9]?$|^RPLP[0-2]$|^RPSA$", genes, ignore.case = TRUE),
      hb = coding & grepl("^HB(A[12]|B|D|E1|G[12]|M|Q1|Z)$|^HB[AB]-|^HBQ1[AB]?$", genes, ignore.case = TRUE),
      sex = toupper(genes) == "XIST" | (!is.na(ann$chr) & ann$chr == "Y"),
      ercc = grepl("^ERCC-", genes)
    )
  }
  out
}

# Symbols that the OrgDb gives to the mitochondrial Ensembl IDs. rsem_read()
# names rows by these, and one of them is not an MT- symbol: the ID of MT-TS1
# maps to MIR12136, a microRNA inside it.
.mito_symbols <- function(species) {
  orgdb <- .orgdb_package(species)
  db <- .quietly(getExportedValue(orgdb, orgdb))
  keys <- intersect(.mito_ensembl[[species]], .quietly(AnnotationDbi::keys(db, keytype = "ENSEMBL")))
  if (!length(keys)) return(character(0))
  a <- .quietly(AnnotationDbi::select(db, keys = keys, keytype = "ENSEMBL", columns = "SYMBOL"))
  unique(a$SYMBOL[!is.na(a$SYMBOL)])
}

# Gene type and chromosome for each row: symbols are looked up as SYMBOL and
# Ensembl IDs (rows without a symbol in rsem_read()) as ENSEMBL. A symbol with
# several entries is protein-coding when any entry is.
.gene_annotation <- function(genes, species, groups) {
  orgdb <- .orgdb_package(species)
  if (!.quietly(requireNamespace(orgdb, quietly = TRUE))) {
    .stopf("The package '%s' is required. Install it with BiocManager::install(\"%s\").", orgdb, orgdb)
  }
  db <- .quietly(getExportedValue(orgdb, orgdb))
  columns <- c(if (any(c("protein_coding", "ribo", "hb") %in% groups)) "GENETYPE", if (any(c("mito", "sex") %in% groups)) "CHR")
  # CHR is no longer listed by columns() but is still served by select().
  if ("GENETYPE" %in% columns && !"GENETYPE" %in% .quietly(AnnotationDbi::columns(db))) {
    .stopf("'%s' does not provide GENETYPE (gene type), which these gene groups need.", orgdb)
  }

  stable <- .strip_ensembl_version(genes)
  is_ens <- grepl("^ENS[A-Z]*G[0-9]+$", stable)
  type <- chr <- rep(NA_character_, length(genes))
  for (keytype in c("SYMBOL", "ENSEMBL")) {
    idx <- if (keytype == "SYMBOL") which(!is_ens) else which(is_ens)
    keys <- intersect(unique(stable[idx]), .quietly(AnnotationDbi::keys(db, keytype = keytype)))
    if (!length(keys)) next
    # Chromosome via CHR is marked as deprecated in AnnotationDbi but still
    # served by the OrgDb packages; the warning is not useful to the user.
    a <- tryCatch(
      suppressWarnings(.quietly(AnnotationDbi::select(db, keys = keys, keytype = keytype, columns = columns))),
      error = function(e) {
        .stopf("Could not read %s from '%s': %s", paste(columns, collapse = " and "), orgdb, conditionMessage(e))
      }
    )
    key <- a[[keytype]]
    if ("GENETYPE" %in% columns) {
      pc <- tapply(a$GENETYPE == "protein-coding" & !is.na(a$GENETYPE), key, any)
      first <- tapply(a$GENETYPE, key, function(v) v[!is.na(v)][1])
      t <- ifelse(pc, "protein-coding", first)
      type[idx] <- t[stable[idx]]
    }
    if ("CHR" %in% columns) {
      ch <- tapply(a$CHR, key, function(v) {
        v <- unique(v[!is.na(v)])
        # Pseudoautosomal genes are on X and Y and are expressed in both
        # sexes, so only genes found on Y alone count as chromosome Y.
        if (length(v) && all(v == "Y")) "Y" else if ("MT" %in% v) "MT" else setdiff(v, "Y")[1]
      })
      chr[idx] <- ch[stable[idx]]
    }
  }
  data.frame(type = unname(type), chr = unname(chr), stringsAsFactors = FALSE)
}
