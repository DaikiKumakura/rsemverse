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

## Functions added in 0.3.0 (November 20, 2026, rsemverse 0.3.0)

`encode-k562-0.3.R` uses the six gene files above, the six transcript files
of the same RSEM runs (ENCFF875PEA, ENCFF180IOW, ENCFF848LJZ, ENCFF371ODR,
ENCFF662WOK, ENCFF180KQR) and the GENCODE V29 annotation the ENCODE reference
was built from (`gencode.v29.annotation.gtf.gz`). All downloads are checked
by MD5. Run on Windows with R 4.5.3, org.Hs.eg.db 3.22.0 and edgeR 4.8.2.

```bash
Rscript validation/encode-k562-0.3.R my-work-folder
```

| Check | Result |
| --- | --- |
| `rsem_read(metadata = NULL)` | found the six files, named them by file name, and returned the same table as with the sample table (0.72 s) |
| Effective length 0 | 5,908 genes in every sample; 48 more in three of the six samples (5,956 or 5,908 per sample). Median coefficient of variation of the length across samples, for genes never 0: 0.0165; 3,179 such genes vary more than twofold |
| `rsem_detection_curve()` | TPM >= 1 in >= 2 samples: 13,993 genes, the number `rsem_filter()` keeps with the same settings; 15,002 genes in >= 1 sample, 11,741 in all 6 |
| `rsem_fit_reference()` on the controls (log2 TPM + 1, the 13,993 genes) | 136 genes have no variation between the two controls and were removed; each control has z = ±0.707 for every gene, as two reference samples imply; a sample's values do not change when other samples are left out |
| Knockdown targets, z against controls | PTBP1 −1,815 and −1,997 in its knockdown; SRSF1 −27.5 and −27.6. The large PTBP1 values come from two controls with almost the same PTBP1 level (a tiny standard deviation); with two reference samples, z-scores are a ranking aid, not a test |
| `rsem_normalize("tmmwsp")` | on the 12,808 genes kept by the edgeR rule, within 2.6e-16 (relative) of `edgeR::calcNormFactors(method = "TMMwsp")` then `edgeR::cpm()` |
| GTF reading | 58,721 gene rows (19 s for the compressed file) |
| `rsem_biotype_profile()` (share of TPM) | protein-coding 90.6–92.8%; lincRNA 0.67–0.89%; antisense 0.32–0.39%; Mt_rRNA 1.01–1.80%; rRNA 0.08–0.45%; snoRNA 0.05–0.10%; Mt_tRNA 0 (tRNAs are shorter than the fragments, so RSEM gives them length 0 and TPM 0); not in the annotation 2.1–3.5%, mostly the 97 ERCC spike-ins (2.4% of CTRL_1) |
| `lncRNA` in GENCODE V29 | reported as `NA`, not 0: V29 still splits long non-coding genes into `lincRNA`, `antisense` and others. `rsem_subset("lncrna")` uses `lncRNA` and `lincRNA` |
| Matching through the original IDs | before matching through `attr(x, "gene_ids")`, 4.4–5.9% of the TPM was "not in the annotation" because OrgDb symbols and GENCODE gene names differ; afterwards 2.1–3.5% |
| `rsem_subset(remove = c("rrna", "mito_rrna"))` | removed 552 rRNA and 2 mitochondrial rRNA genes; every ribosomal protein gene stayed |
| `rsem_read(level = "isoform")` | 207,507 transcripts of 59,526 genes (0.96 s). Summed per gene, transcript TPM equals the gene file's TPM within 0.05 TPM (rounding in the files) |
| `rsem_select_isoforms()` (TPM >= 1) | any sample: 59,603 transcripts of 14,449 genes; all samples: 25,948 of 10,810; with >= 10% of the gene's top transcript: 47,049; top transcript in any sample: 19,106; top transcript by mean: 12,610, one per gene |
| `rsem_write(annotation = gtf)` | `.tsv.gz` with `gene`, `gene_id`, `gene_type` and the samples; values read back unchanged; 89 rows list several gene IDs |
