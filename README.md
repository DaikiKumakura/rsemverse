# rsemverse

<!-- badges: start -->
[![R-universe](https://daikikumakura.r-universe.dev/badges/rsemverse)](https://daikikumakura.r-universe.dev/rsemverse)
[![R-CMD-check](https://github.com/DaikiKumakura/rsemverse/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/DaikiKumakura/rsemverse/actions/workflows/R-CMD-check.yaml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24210459.svg)](https://doi.org/10.5281/zenodo.24210459)
<!-- badges: end -->

**RSEM results in R, one step at a time, with the pipe.** rsemverse reads
`*.genes.results` and `*.isoforms.results` files, checks them, and turns them into the
tables that plots, scores and downstream tools need. Every step takes a table and returns a
table, so steps chain with `|>` the way tidyverse functions do, and each step prints one line
about what it did.

<a id="install"></a>

## Installation

```r
install.packages("BiocManager")
BiocManager::install("DaikiKumakura/rsemverse")
BiocManager::install("org.Hs.eg.db")   # human; org.Mm.eg.db for mouse (only for gene symbols)
```

## Quick start

Three small RSEM files come with the package, so this runs as it is:

```r
library(rsemverse)
ex   <- system.file("extdata", "example", package = "rsemverse")
meta <- file.path(ex, "metadata.csv")

logcpm <- file.path(ex, "data") |>
  rsem_read(meta, species = "human", expr = "count") |>
  rsem_filter() |>
  rsem_subset("human", keep = "protein_coding", remove = c("mito", "ribo")) |>
  rsem_normalize(method = "tmm") |>
  rsem_transform("log2")
```

## Documentation

[Package website](https://daikikumakura.github.io/rsemverse/) · [Function reference](https://daikikumakura.github.io/rsemverse/reference/) · [Changelog](https://daikikumakura.github.io/rsemverse/news/) · [Releases](https://github.com/DaikiKumakura/rsemverse/releases) · [R-universe](https://daikikumakura.r-universe.dev/rsemverse) · [Usage guide](https://github.com/DaikiKumakura/rsemverse/blob/main/README-GUIDE.md)

## Validation and limitations

- **Validation:** Compared selected filtering, normalization and import results with established R packages and checked real RSEM output.
- **Results:** Selected outputs agreed with the corresponding reference implementations; the checks also include ENCODE K562 data and a small RSEM run.
- **Limitations:** These comparisons cover the recorded inputs, settings and package versions.
- **Reproduction:** [validation](https://github.com/DaikiKumakura/rsemverse/tree/main/validation).

## Citation

Kumakura, D. (2026). *rsemverse: From RSEM Results to Analysis-Ready Gene Tables* [Software]. [DOI: 10.5281/zenodo.24210459](https://doi.org/10.5281/zenodo.24210459).

Citation metadata: [CITATION.cff](CITATION.cff).

## License

This project is licensed under MIT. See [LICENSE.md](LICENSE.md) for the license terms.
