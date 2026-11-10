# Real-data check on ENCODE RSEM files

`encode-k562.R` runs `rsem_read()` on six ENCODE gene quantification files
(human K562 cells, GRCh38, GENCODE V29) and prints the numbers below. It
downloads about 65 MB and checks the MD5 sums.

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
