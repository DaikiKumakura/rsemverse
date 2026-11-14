# Real-data check on ENCODE RSEM files

`encode-k562.R` runs `rsem_read()` and the functions that prepare the table
on six ENCODE gene quantification files (human K562 cells, GRCh38, GENCODE
V29) and prints the numbers below. It downloads about 65 MB and checks the MD5
sums. Where edgeR, DESeq2 and limma are installed, it compares the normalized
values with theirs.

```bash
Rscript validation/encode-k562.R my-work-folder
```

Files: control (ENCSR129RWD: ENCFF587DSG, ENCFF774RTC), PTBP1 knockdown
(ENCSR527IVX: ENCFF075OME, ENCFF896CIS), SRSF1 knockdown (ENCSR066VOO:
ENCFF711WSO, ENCFF322CJL). Both knockdown experiments use the control above.

## Result (November 10, 2026, org.Hs.eg.db 3.22.0)

| | |
| --- | --- |
| Genes per file | 59,526: 58,735 versioned Ensembl IDs, 97 ERCC spike-ins, 649 numeric IDs, 45 ending in `_PAR_Y` |
| Rows returned | 59,427 |
| With a symbol | 36,385 genes (89 symbols merged from several IDs) |
| Kept under the original ID | 23,141 genes (38.9%): 22,350 Ensembl IDs without a symbol in the OrgDb, plus 791 non-Ensembl IDs |
| Column sums | equal to the input, for counts and for TPM |
| Share of the total in the kept rows | 3.76–4.90% of the expected count, 5.89–7.99% of the TPM, per sample |
| Run time | 3.6 s for six samples |
| Knockdown targets (mean count, knockdown / control) | PTBP1 0.35 in the PTBP1 knockdown; SRSF1 0.19 in the SRSF1 knockdown |

The kept genes were not examined further. The knockdown ratios only show that
sample columns and symbol rows line up; the experiment has two replicates per
group and is not analysed here.

## Functions that prepare the table (November 14, 2026, rsemverse 0.2.0)

Run on Windows with R 4.5.3, org.Hs.eg.db 3.22.0, edgeR 4.8.2, DESeq2 1.50.2
and limma 3.66.0. Times are for six samples.

| Step | Result |
| --- | --- |
| `rsem_filter()` `"edger"` (counts) | 12,808 of 59,427 genes kept, the same genes as `edgeR::filterByExpr()`; they hold 99.75–99.78% of the expected count per sample; PTBP1, SRSF1 and PTBP2 kept |
| `rsem_filter()` `"detected"` (TPM >= 1 in >= 2 of 6) | 13,993 genes |
| `rsem_filter()` `"nonzero"` / `"variable"` | 26,384 genes / 2,000 genes (0.48 s) |
| `rsem_subset()` genes per group | protein-coding 19,328; mito 37; ribo 86; hb 10; sex 256; ERCC 0 (these files have no ERCC rows under that name) |
| Share of expected count per sample | mito 6.56–7.60%, ribo 6.41–8.29%, hb 0.37–0.69%, sex 0.05–0.10% |
| Share of TPM per sample | mito 13.32–17.40%, ribo 16.49–21.73%, hb 1.02–1.83% |
| `rsem_subset()` protein-coding without mito, ribo, hb and sex | 19,182 genes (2.7 s) |
| `rsem_normalize()` on the filtered counts | CPM, TMM, upper quartile, DESeq2 and quantile within 3e-16 (relative) of edgeR, DESeq2 and limma; each under 0.1 s |
| Knockdown targets after TMM (mean TMM-CPM, knockdown / control) | PTBP1 0.35, SRSF1 0.19 (the same as on raw counts) |
| `rsem_transform()` | log2 then z-score on 12,808 genes; rank on all 59,427 rows (0.08 s) |
| `rsem_long()` | 76,848 rows with `group` (text) and `replicate` (integer) from the sample table |

The first run found only 13 of the 37 mitochondrial genes: the OrgDb links
only the protein-coding ones to Ensembl IDs, so `rsem_read()` kept the rRNA
and tRNA genes under their IDs, and the ID of MT-TS1 maps to MIR12136, a
microRNA inside it. `rsem_subset()` now also uses a built-in list of the 37
Ensembl IDs and the symbols the OrgDb gives them. K562 cells come from a
female donor and express hemoglobin; 3 chromosome Y genes have a nonzero
count in any sample and XIST has a mean TPM of 114.
