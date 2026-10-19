# Data for your own analysis

Starting point for your own analysis, with the approach of your choice: 27 mouse cerebellum samples, quantified with **RSEM on STAR alignments** (GENCODE M39 / GRCm39; [why this method](../docs/pipeline.md#choosing-the-counting-method)). Read the files; do not edit them.

| File | Content |
|---|---|
| `rsem_expected_counts.csv.gz` | estimated number of reads per gene and sample (rows = genes, columns = 27 samples). Not whole numbers, so round them if a tool needs integers (e.g. DESeq2) |
| `rsem_tpm.csv.gz` | TPM = transcripts per million: expression normalised for gene length and sequencing depth, so genes and samples are comparable. Same layout |
| `rsem_effective_length.csv.gz` | length of each gene that RSEM used for each sample, in bases (it differs per sample). Only needed to correct counts for gene length yourself, e.g. with `tximport` / DESeq2. Same layout |
| `samples.csv` | one row per sample, in the same order as the table columns. `SRR`: sample id (= column name in the tables); `group`: experimental group; `age_months`: age of the mouse; `sex`; `cohort`: which experiment it belongs to; `mouse_id`: the paper's animal number |
| `genes_M39.csv` | translates each row name of the tables (`gene_id`, an Ensembl id such as `ENSMUSG00000000001.6`) into the usual gene name (`symbol`, e.g. `Gnai3`) |

Load with `read.csv(..., row.names = 1)` (R) or `pd.read_csv(..., index_col=0)` (Python). The row names are Ensembl gene ids with a version suffix (`.6`); match them to symbols via `genes_M39.csv`.

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
- `../results/euler/rsem_results.tar.gz` (72 MB): RSEM's per-run output files, the input for `tximport` (as in the tutorial) if you prefer to build the gene table yourself. Check with `shasum -a 256 -c ../metadata/data_checksums.txt`; unpack with `tar xzf` into `../results/`.
- Our own results (DESeq2 tables, aging vs rejuvenation): [../results/README.md](../results/README.md).
