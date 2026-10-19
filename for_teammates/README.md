# Data for your own analysis

Everything you need for the tutorial's analysis part (PCA, clustering, differential expression, enrichment, deconvolution). Quantification is **RSEM on STAR alignments**, GENCODE M39 (GRCm39), 27 cerebellum samples ([why this method](../docs/pipeline.md#choosing-the-counting-method)). Read these files; do not edit them.

| File | Content |
|---|---|
| `rsem_expected_counts.csv.gz` | RSEM expected counts, genes × 27 runs (not whole numbers; round them for DESeq2) |
| `rsem_tpm.csv.gz` | TPM, genes × 27 runs (for PCA, clustering, ANOVA, deconvolution) |
| `rsem_effective_length.csv.gz` | effective gene length; with the counts it rebuilds DESeq2's length correction |
| `samples.csv` | one row per run, in the column order of the tables above: `SRR` (the column names of the tables), `group`, `age_months`, `sex`, `cohort`, `mouse_id` |
| `genes_M39.csv` | `gene_id` → `symbol` |

Load the tables with `read.csv(..., row.names = 1)` in R or `pd.read_csv(..., index_col=0)` in Python. Gene ids keep their version (`ENSMUSG00000000001.6`).

**Filter genes first.** The tables hold all 78,289 annotated genes, most with no or very few reads; the tutorial's expression filter applies (we kept genes with ≥ 10 counts in at least as many samples as the smallest group).

**Groups** (`group`; compare only within a cohort):

| `cohort` | `group` values (mice) | Compare |
|---|---|---|
| Aging | `aging 3mo` (4), `aging 18mo` (4), `aging 21mo` (4) | 18 or 21 vs 3 months; model sex as a covariate (2 male + 2 female per age) |
| Dietary Intervention | `AL` (4), `aDR` (5); all female, 19 months | aDR vs AL |
| Young Plasma Injection | `PBS` (3), `YMP` (3); all male, 19 months | YMP vs PBS |

**Caution.** The cohorts differ for technical reasons (and AL is female, PBS male), so never pool them. SRR21355389 (aging 21mo) is blood-contaminated; we exclude it from the aging analysis.

**Optional, larger file (not in git):** `rsem_results.tar.gz` (72 MB, RSEM's per-run files for `tximport`) from the GitHub release `data-v1`; check it with `shasum -a 256 -c ../metadata/data_checksums.txt`, unpack with `tar xzf rsem_results.tar.gz -C ../results/`.

Our own results (DESeq2 tables, aging vs rejuvenation) are in `../results/`, described in [../results/README.md](../results/README.md).
