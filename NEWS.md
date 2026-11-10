# rsemverse 0.1.0

First version.

* `rsem_read()` reads RSEM `*.genes.results` files listed in a metadata CSV
  (`sample_id`, `ref`) and returns expected counts or TPM per gene symbol and
  sample, using `org.Hs.eg.db` (human) or `org.Mm.eg.db` (mouse). Ensembl
  versions are removed before the lookup, and rows that share a symbol are
  summed.
