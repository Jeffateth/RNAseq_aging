# Data for your own analysis

Everything the tutorial's analysis part starts from: 27 cerebellum samples, quantified with **RSEM on STAR alignments** (GENCODE M39 / GRCm39; [why this method](../docs/pipeline.md#choosing-the-counting-method)). Read the files; do not edit them.

| File | Content |
|---|---|
| `rsem_expected_counts.csv.gz` | expected counts, genes × 27 samples (not whole numbers; round for DESeq2) |
| `rsem_tpm.csv.gz` | TPM (for PCA, clustering, ANOVA, deconvolution) |
| `rsem_effective_length.csv.gz` | effective gene length; with the counts it rebuilds DESeq2's length correction |
| `samples.csv` | one row per sample, in the column order of the tables: `SRR` (the table's column names), `group`, `age_months`, `sex`, `cohort`, `mouse_id` |
| `genes_M39.csv` | `gene_id` → `symbol` |

Load with `read.csv(..., row.names = 1)` (R) or `pd.read_csv(..., index_col=0)` (Python). Gene ids keep their version (`ENSMUSG00000000001.6`).

**Filter genes first.** The tables hold all 78,289 annotated genes, most with no or few reads. We kept genes with ≥ 10 counts in at least as many samples as the smallest group.

## Groups
Compare only within a cohort; never pool cohorts (they differ for technical reasons, and AL is female, PBS male).

| `cohort` | `group` (mice) | Compare |
|---|---|---|
| Aging | `aging 3mo`, `aging 18mo`, `aging 21mo` (4 each, 2 male + 2 female) | 18 or 21 vs 3 months, with sex as covariate |
| Dietary Intervention | `AL` (4), `aDR` (5); female, 19 months | aDR vs AL |
| Young Plasma Injection | `PBS` (3), `YMP` (3); male, 19 months | YMP vs PBS |

**Outlier.** SRR21355389 (aging 21mo) is potentially blood-contaminated: red-blood-cell and neutrophil genes are 50–150× the other aging samples (TPM: Hba-a1 4,499 vs median 57; Hbb-bs 4,242 vs 49; S100a8 3,959 vs 25). We exclude it from the aging analysis ([evidence](../docs/pipeline.md)).

## More
- `../results/euler/rsem_results.tar.gz` (72 MB): RSEM's per-run files, for `tximport` as in the tutorial. Check with `shasum -a 256 -c ../metadata/data_checksums.txt`; unpack with `tar xzf` into `../results/`.
- Our own results (DESeq2 tables, aging vs rejuvenation): [../results/README.md](../results/README.md).
