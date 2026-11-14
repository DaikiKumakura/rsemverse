# rsemverse 0.2.0

New functions that take the table from `rsem_read()` and return a table of
the same layout, so they can be chained. None of them needs a gene list or a
sample grouping from the user.

* `rsem_filter()` removes genes: the edgeR `filterByExpr()` rule with all
  samples as one group (default for counts), a value threshold in a number or
  share of samples (default for TPM: TPM >= 1 in 20% of samples), genes that
  are zero everywhere, or all but the most variable genes.
* `rsem_normalize()` normalizes between samples: CPM, TMM, upper quartile and
  DESeq2 median of ratios for counts, and quantile normalization. The tests
  check that the results equal those of edgeR 4.8.2, DESeq2 1.50.2 and limma
  3.66.0; the package does not depend on them.
* `rsem_transform()` applies log2 (with a pseudocount), per-gene z-scores or
  centering, and within-sample ranks, in the order given.
* `rsem_subset()` keeps or removes built-in gene groups: protein-coding genes
  (gene type from the OrgDb), mitochondrial-genome genes, cytosolic ribosomal
  protein genes, hemoglobin genes, chromosome Y genes with XIST, and ERCC
  spike-ins. Mitochondrial rRNA and tRNA genes, which have no Ensembl link in
  the OrgDb and keep their Ensembl ID in `rsem_read()`, are found through a
  built-in list of the 37 mitochondrial Ensembl IDs.
* `rsem_long()` returns one row per gene and sample with every column of the
  sample table added as it is (numbers stay numeric, other columns stay text).

On the six ENCODE K562 files in `validation/`, every function ran on all
59,427 rows, and the normalized values were within 3e-16 (relative) of edgeR,
DESeq2 and limma. See `validation/README.md`.

# rsemverse 0.1.0

First version.

* `rsem_read()` reads RSEM `*.genes.results` files listed in a metadata CSV
  (`sample_id`, `ref`) and returns expected counts or TPM per gene symbol and
  sample, using `org.Hs.eg.db` (human) or `org.Mm.eg.db` (mouse). Ensembl
  versions are removed before the lookup, and rows that share a symbol are
  summed.
