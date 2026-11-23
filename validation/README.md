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

## Changes in 0.4.0 (November 21, 2026, rsemverse 0.4.0)

`encode-k562-0.4.0.R` runs on the files that `encode-k562-0.3.R` downloaded.
Windows, R 4.5.3, org.Hs.eg.db 3.22.0, edgeR 4.8.2.

| Check | Result |
| --- | --- |
| Gene IDs after `rsem_filter()` (TPM >= 1 in >= 20% of samples) | 13,993 rows kept; `attr(x, "gene_ids")` has 14,034 rows (one per original ID, 37 rows of the table sum several IDs) and covers every kept symbol; `rsem_write()` of the filtered table writes `gene`, `gene_id` and the samples. Before the fix the attribute was gone |
| Gzip input | the six gene files compressed with gzip and found with `metadata = NULL` gave the same table as the plain files (2.6–3.0 s) |
| `--append-names` style IDs | files rewritten with the symbol appended to every gene that has one (`ENSG...._TP53`; made from the real files, not written by RSEM) gave the same 59,427 rows and values |
| `rsem_filter(group = )` on counts | one group: 12,808 genes (as before); groups of two samples each (control, PTBP1, SRSF1): 14,755 genes, equal to `edgeR::filterByExpr(group = )` |
| `rsem_biotype_profile(total = )` | Mt_rRNA share of TPM in CTRL_1: 1.575% on the unfiltered table; 1.578% when the filtered table is its own denominator; 1.575% with `total = colSums(unfiltered)` |
| Isoforms through `rsem_filter()` | 49,759 of 207,507 transcripts kept (TPM >= 1 in >= 2 samples), with the same 49,759 rows in `attr(x, "tx2gene")` |
| `rsem_mapping_qc()` | not run on real data: the ENCODE downloads have no RSEM `.cnt` files. The tests use files written to the format in RSEM's `cnt_file_description.txt`, and a file whose totals do not add up stops the function |
| `rsem_sample_qc()` on the filtered table (13,993 genes, 2,000 for the PCA) | PC1–PC5 explain 41.7, 32.8, 9.7, 8.1 and 7.7% of the variance of those genes; scores and loadings equal `stats::prcomp()` and the Spearman matrix equals `stats::cor()` (12,000+ detected genes); each sample's most similar sample is its own replicate (median correlation with the others 0.944–0.958) |

## Changes in 0.5.0 (November 21, 2026, rsemverse 0.5.0)

`encode-k562-0.5.0.R`, on the files that `encode-k562-0.3.R` downloaded.
Windows, R 4.5.3, org.Hs.eg.db 3.22.0.

| Check | Result |
| --- | --- |
| `rsem_check()` on the six gene files | no errors; one warning per file for genes with effective length 0 (5,956 or 5,908); no TPM total off 1e6 |
| `rsem_check(level = "isoform")` on the six transcript files | only the zero-length warnings |
| `rsem_read(validate_reference = TRUE)` | accepted the six files (the transcript lists of all 59,526 genes are equal) |
| One gene given an extra transcript in one file | `rsem_read()` stopped, naming the file and the gene; `rsem_check()` reported one `reference_mismatch` error for that sample |
| One file with its rows shuffled | stopped by default; with `allow_reorder = TRUE` the table equals the one from the original files |
| `rsem_read(id = "annotation")` with the GENCODE V29 annotation | 57,938 rows from 59,526 genes (58,721 genes named); column sums equal to the OrgDb route. Of the 31,544 names present in both routes, 87 differ in value: the two sources name some genes differently, and GENCODE gives the `_PAR_Y` copies the name of the X copy, so they are summed (216 names merge several IDs, 45 of them `_PAR_Y` IDs) |
| Posterior columns and `.cnt` files | not run on real data: the ENCODE downloads were made without `--calc-pme` and have no `.cnt` files. The tests use files written to RSEM's documented format |

With the stricter name matching of 0.5.1 (an appended name must be exactly the `gene_name` of the gene before the underscore), the GENCODE V29 numbers above did not change: 57,938 rows, 58,721 genes named, 87 of 31,544 shared names differing in value.

## 0.6.0 (`validation/encode-k562-0.6.0.R`, November 22, 2026)

On the same six ENCODE K562 runs with their transcript files:

| Check | Result |
|---|---|
| `rsem_reconcile()`, 3 columns x 6 samples | no gene beyond the rounding bound; the largest difference was 0.05 (TPM) and 0.04 (count) |
| Isoform file of another sample under one sample's name | 20,405 of 59,526 genes beyond the bound for that sample, none for the others |
| `rsem_id_report()` with `org.Hs.eg.db` | 89 rows combine several gene IDs (largest: 7) |
| TPM total of a run | within 0.1 of 1e6 |

## 0.7.0 (`validation/encode-k562-0.7.0.R`, November 23, 2026)

On the same six ENCODE K562 runs:

| Check | Result |
|---|---|
| `rsem_read(expr = "fpkm")` against `1e3 / l_bar * TPM` (documentation) | median ratio 1.0000 in all six samples (genes with TPM > 10) |
| `rsem_check(level = "isoform")`, transcript lengths | no transcript has another length in any file |
| `rsem_as_tximport(zero_length = "drop")` against `tximport(type = "rsem")` | counts and lengths identical (difference 0); 53,570 of 59,526 genes remain (5,908 to 5,956 per sample have effective length 0) |
| `rsem_sample_qc()` relative log expression | median -0.025 to 0.020, IQR 0.28 to 0.39 |

Then RSEM 1.3.3 itself was built and run on a toy data set (`--calc-pme --calc-ci --append-names`; `--estimate-rspd`); its output is in `tests/testthat/fixtures/rsem-real/`. `rsem_read()` read every posterior column, the interval and quartile variation, `IsoPct_from_pme_TPM`, the names after `--append-names` (including an underscore inside the name), the `.cnt` file and its histogram (the reads in it add up to `N_tot`), and the `.model` files; `rsem_reconcile()` found no gene beyond rounding. A gene with a 60 nt transcript had effective length 0, TPM 0 and 50 expected counts. Not read: the `.theta` file and the sequencing error sections of the `.model` file.

## 0.7.2 (November 23, 2026)

RSEM 1.3.3 was run once more with 600 unalignable reads added (`tests/testthat/fixtures/rsem-real/sampD.stat`). The `.theta` file then had 0.2182 (= 600 / 2750) as the first value of its second line and 0.2222 (= 600 / 2700) as the first value of its third line; `rsem_noise_qc()` reads both. GENCODE V29 read through `.read_annotation()`: 58,721 genes with name and type, 23.8 s (16.4 s of it reading the gzip file).
