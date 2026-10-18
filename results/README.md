# Results: files and how to read them

Small outputs are kept in git; large or regenerable files are git-ignored. Terms: [../docs/glossary.md](../docs/glossary.md). Gene ids are GENCODE M39 Ensembl ids with version (e.g. `ENSMUSG00000102693.1`) plus the symbol.

| Folder | Content | In git |
|---|---|---|
| `qc/` | fastp report per run (`<SRR>.fastp.json`; the `.html` is regenerated), `polyG_sensitivity_SRR21355398.csv` (genes whose kallisto counts changed > 5 % with poly-G trimming, one-run test); `qc/raw/`, `multiqc_raw.html`: FastQC of the raw reads | json, csv |
| `kallisto/`, `rsem/` | per-run kallisto and RSEM output from Euler (notebook 01) | no |
| `star/` | `counts_unstranded.csv.gz` (STAR counts, genes × runs) and `star_qc.csv` (per run: `input_pairs_M`, `unique_%`, `multi_%`, `STAR_assigned_to_gene_%`, `frac_read1_same_strand` ≈ 0.5 = unstranded) | yes |
| `de/<method>/` | DESeq2 table per comparison for the method set in notebook 03 | no |
| `compare/` | kallisto vs STAR vs RSEM (notebook 04) | summary tables |
| `aging_vs_rejuvenation/` | the main analysis (notebook 05) | yes |

## Direction and significance
Every comparison is *first group vs second*; **positive** log2FC or NES = **higher in the first group**. Significant = `padj` < 0.05. A **reversal** moves opposite to aging.

| Name | Comparison |
|---|---|
| `aging_21v3`, `aging_18v3` | 21 (18) months vs 3 months |
| `aDR` | dietary restriction vs free food (AL) |
| `YMP` | young plasma vs saline (PBS) |

## Tables
**DE tables** (`de/<method>/*.csv`, `compare/{star,rsem}_*.csv`), one row per gene tested: `gene_id`, `symbol`, `baseMean` (average normalised count), `log2FoldChange`, `lfcSE`, `stat` (Wald = log2FC / SE, ranks genes for GSEA), `pvalue`, `padj` (Benjamini–Hochberg; `NA` = too few counts or a count outlier, not tested).

**`aging_vs_rejuvenation/`** (one row per method × treatment unless noted):
- `overlap_with_aging.csv`: `genes` tested in both, `aging_DEG`, `intervention_DEG`, `shared` (and `expected_by_chance`), `shared_reversal` / `shared_same`, `fisher_p`, `lfc_spearman_all`, `lfc_spearman_agingDEG` (**negative = reversal**).
- `signature_score_tests.csv`: `treated`, `control`, `n_t`, `n_c`, `diff` (treated − control; **negative = rejuvenating direction**), `t_test_p` (Welch, rough with 3–5 mice).
- `gsea_hallmark.csv`: one row per method × comparison × Hallmark set (50 × 3 × 3 = 450): `pathway`, `pval`, `padj` (over the 50 sets), `ES`, `NES` (**positive = pathway up in the first group**), `size`, `log2err` (ignore). Example: `aging_21v3`, `HALLMARK_INTERFERON_ALPHA_RESPONSE`, NES +2.5, padj < 0.05 = interferon genes are higher in 21- than in 3-month mice.
- `gsea_concordance.csv`: of the sets significant in aging, per collection (`Hallmark`, `GO_BP`) × method × treatment: `aging_sets`, `significant_in_both`, `reversed`, `same_direction`, `NES_spearman_all_sets` (negative = reversal).
- `paper_CAS_named_genes.csv`: the paper's Fig. 2B aging genes with `lfc_*` / `padj_*` for aging, aDR, YMP (method set in `SHOW`, default kallisto).

**`compare/`:**
- `counts_per_sample.csv`: `spearman_<a>_<b>` per pair of methods and total counts in millions per method (STAR counts only unique reads, so its total is lower).
- `de_agreement.csv`: per pair (`method_a`, `method_b`) and `contrast`: `total_DEG_a/b` (all DEGs), then on genes with a valid p-value in both (`genes_tested_both`): `DEG_a/b`, `shared`, `jaccard` (shared ÷ combined), `same_direction` (share of shared DEGs with the same sign), `lfc_spearman_all` / `_sig` (log2FC rank correlation over all genes / genes significant in either).
- `method_evidence.csv`: one row per method, the criteria of [../docs/pipeline.md](../docs/pipeline.md#choosing-the-counting-method): `counted_pairs_pct_median/min/max`, `genes_detected`, `replicate_r_median`, `dispersion_aging/aDR/YMP`, `DEG_without/with` and `jaccard_outlier` (aging DEGs without / with SRR21355389), `aging_genes_up_of_5`, `aDR_genes_up_of_23` (positive controls), `mean_jaccard_with_others`.

**Pipeline outputs (not in git):** `kallisto/<SRR>/abundance.tsv` (`target_id`, `length`, `eff_length`, `est_counts`, `tpm`; `run_info.json` has `p_pseudoaligned`); `rsem/<SRR>/<SRR>.genes.results` (expected counts, TPM, effective length per gene).
