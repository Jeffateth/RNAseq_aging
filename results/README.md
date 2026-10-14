# Results: files and how to read them

Small outputs kept in git; large or regenerable files are git-ignored (`.gitignore`). Terms: [../docs/glossary.md](../docs/glossary.md).

| Folder | Content | In git |
|---|---|---|
| `qc/` | fastp report per run (`<SRR>.fastp.json` in git; the `.html` version is regenerated); `polyG_sensitivity_SRR21355398.csv`: genes whose kallisto counts changed > 5 % when fastp's poly-G trimming was added (one-run test, notebook 02; the pipeline now always uses it) | yes |
| `qc/raw/`, `qc/multiqc_raw.html` | FastQC reports on the untrimmed reads, combined MultiQC report (notebook 02 or 01) | no |
| `kallisto/` | per-run kallisto estimates (notebook 02 or 01) | no |
| `rsem/` | per-run RSEM gene results from Euler (notebook 01) | no |
| `star/` | STAR counts and mapping statistics from Euler (notebook 01) | counts + QC table |
| `de/<method>/` | DESeq2 table per comparison for the method set in notebook 03 | no (regenerated) |
| `compare/` | kallisto vs STAR vs RSEM (notebook 04) | summary tables |
| `aging_vs_rejuvenation/` | the main analysis (notebook 05) | yes |

Gene identifiers everywhere: GENCODE M39 Ensembl ids with version (e.g. `ENSMUSG00000102693.1`) and the gene symbol.

## Reading the comparisons (direction)
Every comparison is *first group vs second group*, and a **positive** log2 fold change (log2FC) or enrichment score (NES) means **higher in the first group**.

| Name in the files | First vs second group | Positive means higher in |
|---|---|---|
| `aging_21v3` | 21-month vs 3-month mice | 21-month mice |
| `aging_18v3` | 18-month vs 3-month mice | 18-month mice |
| `aDR` | aDR vs AL (dietary restriction vs free food) | aDR mice |
| `YMP` | YMP vs PBS (young plasma vs saline) | YMP mice |

Rule of thumb: **significant = `padj` < 0.05** (adjusted for the number of tests). A **reversal** is a gene or pathway that moves in the opposite direction in a treatment than in aging.

## Column dictionary

### Differential expression tables (`de/<comparison>.csv`, `compare/star_<comparison>.csv`)
One row per gene tested.

| Column | Meaning |
|---|---|
| `gene_id`, `symbol` | GENCODE gene id and gene symbol |
| `baseMean` | average normalised read count across all samples (low = weakly expressed) |
| `log2FoldChange` | log2 of the ratio first group / second group; +1 = twice as high, −1 = half |
| `lfcSE` | standard error of that log2 fold change |
| `stat` | Wald statistic = `log2FoldChange` / `lfcSE`; its sign is the direction, its size the evidence; used to rank genes for GSEA |
| `pvalue` | raw p-value of the test that the log2 fold change is 0 |
| `padj` | p-value adjusted for ~15,000 genes tested at once (Benjamini–Hochberg). **`NA`** = not tested: too few counts (DESeq2's automatic filter) or flagged as a count outlier |

### `aging_vs_rejuvenation/gsea_hallmark.csv` (pathway enrichment)
One row per method × comparison × Hallmark gene set (50 sets × 3 comparisons × 2 methods = 300 rows). **GSEA** ranks all genes by the Wald statistic and asks whether a gene set sits mostly at the top (up) or bottom (down) of the ranking.

| Column | Meaning |
|---|---|
| `method` | quantification method used upstream: `kallisto` or `STAR` |
| `contrast` | `aging_21v3`, `aDR` or `YMP` (see direction table) |
| `pathway` | MSigDB Hallmark gene set, e.g. `HALLMARK_INTERFERON_ALPHA_RESPONSE` |
| `pval`, `padj` | permutation p-value and its adjustment over the 50 sets; **significant = `padj` < 0.05** |
| `ES` | enrichment score: how far the gene set is concentrated at one end of the ranking |
| `NES` | **normalised** enrichment score: `ES` corrected for set size so sets are comparable. **Positive = pathway up in the first group, negative = down**; the larger |NES|, the stronger the shift |
| `size` | number of the set's genes present in the ranked list (sets of 15–500 genes are tested) |
| `log2err` | accuracy of `pval` (expected error of its estimate); a technical value, ignore for interpretation |

Example: `aging_21v3`, `HALLMARK_INTERFERON_ALPHA_RESPONSE`, NES +2.5, padj < 0.05 = interferon-response genes are higher in 21-month than in 3-month mice.

```python
import pandas as pd
g = pd.read_csv("results/aging_vs_rejuvenation/gsea_hallmark.csv")
# pathways significant in aging (kallisto), strongest first
g[(g.method == "kallisto") & (g.contrast == "aging_21v3") & (g.padj < 0.05)].sort_values("NES")[["pathway", "NES", "padj"]]
```

### `aging_vs_rejuvenation/gsea_concordance.csv`
For the pathways significant in aging: does a treatment reverse them? One row per gene-set collection × method × treatment.

| Column | Meaning |
|---|---|
| `collection` | `Hallmark` (50 sets) or `GO_BP` (~7,500 Gene Ontology biological-process sets) |
| `intervention` | `aDR` or `YMP` |
| `aging_sets` | number of sets significant in aging |
| `significant_in_both` | of these, significant in the treatment too |
| `reversed` / `same_direction` | of those, treatment NES has the **opposite** / **same** sign as the aging NES |
| `NES_spearman_all_sets` | rank correlation of NES (aging vs treatment) over all sets; **negative = reversal** |

### `aging_vs_rejuvenation/overlap_with_aging.csv`
Do a treatment's DEGs overlap the aging DEGs (21 vs 3 months)? Only genes tested in both comparisons count.

| Column | Meaning |
|---|---|
| `method`, `intervention` | quantification method; `aDR` or `YMP` |
| `genes` | genes tested in both comparisons |
| `aging_DEG`, `intervention_DEG` | DEGs (padj < 0.05) in aging / in the treatment |
| `shared` | DEGs in both; `expected_by_chance` = how many a random overlap would give |
| `shared_reversal` / `shared_same` | shared DEGs moving opposite to / like aging |
| `fisher_p` | p-value that the overlap exceeds chance (Fisher's exact test) |
| `lfc_spearman_all` | rank correlation of log2FC (aging vs treatment) over all genes |
| `lfc_spearman_agingDEG` | the same over aging DEGs only; **negative = reversal** |

### `aging_vs_rejuvenation/signature_score_tests.csv`
Aging-signature score (mean z-score of aging-up genes minus aging-down genes) in treated vs control mice.

| Column | Meaning |
|---|---|
| `treated`, `control` | the two groups (`aDR` vs `AL`, `YMP` vs `PBS`); `n_t`, `n_c` = mice per group |
| `diff` | mean score treated − control; **negative = rejuvenating direction** |
| `t_test_p` | Welch t-test p-value (rough, with 3–5 mice per group) |

### `aging_vs_rejuvenation/paper_CAS_named_genes.csv`
The paper's aging-signature genes named in its Fig. 2B: `lfc_*` and `padj_*` for aging (21 vs 3 months), aDR and YMP (method set in `SHOW`, default kallisto). Many treatment values are not significant (`padj` near 1).

### `compare/` (kallisto, STAR and RSEM)
| File / column | Meaning |
|---|---|
| `counts_per_sample.csv`: `spearman_<a>_<b>`, `kallisto_M`, `STAR_M`, `RSEM_M` | per-sample rank correlation of gene counts for each pair of methods; total counts in millions; STAR counts only uniquely mapped reads, so its total is lower |
| `de_agreement.csv`: `method_a`, `method_b`, `contrast` | one row per pair of methods and comparison |
| `total_DEG_a`, `total_DEG_b` | all DEGs of each method |
| `genes_tested_both`, `DEG_a`, `DEG_b`, `shared` | restricted to genes with a valid p-value in both methods |
| `jaccard` | shared ÷ combined DEGs (1 = identical lists) |
| `same_direction` | share of shared DEGs with the same sign (1 = all) |
| `lfc_spearman_all` / `_sig` | log2FC rank correlation over all genes / over genes significant in either method |

### `compare/method_evidence.csv` (choosing a method)
One row per method (`kallisto`, `STAR`, `RSEM`); criteria and caveats are in [../docs/pipeline.md](../docs/pipeline.md#choosing-the-counting-method).

| Column | Meaning |
|---|---|
| `counted_pairs_pct_median`, `min`, `max` | share of trimmed read pairs counted in a gene, over the 27 runs |
| `genes_detected` | genes with ≥ 10 counts in ≥ 3 runs |
| `replicate_r_median` | median over groups of the median pairwise correlation of log2 CPM between mice of a group |
| `dispersion_aging`, `_aDR`, `_YMP` | median DESeq2 dispersion over genes tested by all methods |
| `DEG_without`, `DEG_with`, `jaccard_outlier` | aging DEGs without / with outlier SRR21355389 and their overlap |
| `aging_genes_up_of_5`, `aDR_genes_up_of_23` | positive controls: paper genes that go up |
| `mean_jaccard_with_others` | mean DEG overlap with the other two methods |

### `star/`
| File | Content |
|---|---|
| `counts_unstranded.csv.gz` | STAR read counts per gene (rows, `gene_id`) × run (columns, SRR); integers |
| `star_qc.csv` | per run (`SRR`): `input_pairs_M` read pairs (millions); `unique_%` uniquely mapped; `multi_%` mapped to several places; `assigned_to_gene_%` counted in a gene; `frac_read1_same_strand` ≈ 0.5 means unstranded library; `group` |

### Pipeline outputs (not in git)
- `kallisto/<SRR>/abundance.tsv`: `target_id` (transcript, with gene id and symbol), `length`, `eff_length` (length adjusted for fragment size), `est_counts` (estimated reads), `tpm` (transcripts per million). `run_info.json` holds the pseudoalignment rate (`p_pseudoaligned`).
- `qc/<SRR>.fastp.json`: reads and bases before/after trimming, adapter statistics, duplication rate.
