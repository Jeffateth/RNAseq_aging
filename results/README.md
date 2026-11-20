# Results: files and how to read them

Small outputs are kept in git; large or regenerable files are git-ignored. Terms: [../docs/glossary.md](../docs/glossary.md). Gene ids are GENCODE M39 Ensembl ids with version (e.g. `ENSMUSG00000102693.1`) plus the symbol.

| Folder | Content | In git |
|---|---|---|
| `qc/` | fastp report per run (`<SRR>.fastp.json`; the `.html` is regenerated); `adapter_content_before_after.png` / `.svg` (the slide figure, notebook 02); `qc/raw/`, `qc/trimmed/`, `multiqc_raw.html`: FastQC of the raw and of the trimmed reads | json, csv, figure |
| `kallisto/`, `rsem/` | per-run kallisto and RSEM output unpacked from the Euler tarballs (notebook 01) | no |
| `euler/` | the result tarballs from Euler (RSEM, STAR, QC; the 626 MB kallisto one is not in git) | yes |
| `star/` | `counts_unstranded.csv.gz` (STAR counts, genes × runs) and `star_qc.csv` (per run: `input_pairs_M`, `unique_%`, `multi_%`, `STAR_assigned_to_gene_%`, `frac_read1_same_strand` ≈ 0.5 = unstranded) | yes |
| `de/<method>/` | DESeq2 table per comparison for the method set in notebook 03 (RSEM); `deg_annotation.csv` = every DEG (padj < 0.05) with `symbol`, `name` and `type` (org.Mm.eg.db; NA = no entry), `log2FC`, `padj`, `direction` | `deg_annotation.csv` only |
| `compare/` | sensitivity checks: kallisto vs RSEM, poly-G trimming (notebook 04) | summary tables |
| `aging_vs_rejuvenation/` | the main analysis (notebook 05) | yes |
| `common/` | reference outputs of the team's common setup (RSEM, notebook 05): `de_<comparison>.csv` (all tested genes: `log2FoldChange`, `lfc_shrunk` (apeglm), `stat`, `pvalue`, `padj`), `gsea_hallmark.csv`, `gsea_gobp.csv`, `reversal.csv` (Spearman correlation and % reversed over the aging DEGs) | yes |
| `diagnostics/` | why DESeq2 instead of ANCOVA: `deseq2_vs_ancova.png` (aging DEGs found by ANCOVA (~ sex + age on log2 counts) vs DESeq2, and their overlap); made by `scripts/deseq2_vs_ancova.R` | yes |

## Direction and significance
Every comparison is *first group vs second*; **positive** log2FC or NES = **higher in the first group**. Significant = `padj` < 0.05. A **reversal** moves opposite to aging.

| Name | Comparison |
|---|---|
| `aging` | old (18 + 21 months) vs young (3 months) |
| `aDR` | dietary restriction vs free food (AL) |
| `YMP` | young plasma vs saline (PBS) |

## Tables
**DE tables** (`de/<method>/*.csv`, `compare/rsem_*.csv`), one row per gene tested: `gene_id`, `symbol`, `baseMean` (average normalised count), `log2FoldChange`, `lfcSE`, `stat` (Wald = log2FC / SE, ranks genes for GSEA), `pvalue`, `padj` (Benjamini–Hochberg; `NA` = too few counts or a count outlier, not tested).

**`aging_vs_rejuvenation/`** (one row per method × treatment unless noted):
- `overlap_with_aging.csv`: `genes` tested in both, `aging_DEG`, `intervention_DEG`, `shared` (and `expected_by_chance`), `shared_reversal` / `shared_same`, `fisher_p`, `lfc_spearman_all`, `lfc_spearman_agingDEG` (**negative = reversal**).
- `signature_score_tests.csv`: `treated`, `control`, `n_t`, `n_c`, `diff` (treated − control; **negative = rejuvenating direction**), `t_test_p` (Welch, rough with 3–5 mice).
- `gsea_hallmark.csv`: one row per method × comparison × Hallmark set (50 × 3 × 2 = 300): `pathway`, `pval`, `padj` (over the 50 sets), `ES`, `NES` (**positive = pathway up in the first group**), `size`, `log2err` (ignore). Example: `aging`, `HALLMARK_INTERFERON_ALPHA_RESPONSE`, NES +2.6, padj < 0.05 = interferon genes are higher in old than in young mice.
- `gsea_concordance.csv`: of the sets significant in aging, per collection (`Hallmark`, `GO_BP`) × method × treatment: `aging_sets`, `significant_in_both`, `reversed`, `same_direction`, `NES_spearman_all_sets` (negative = reversal).
- `zoom_21v3.csv`: the same tests with a 21 vs 3 month (late aging) signature, per method × treatment: `aging_DEG`, `shared_reversed` / `shared_like_aging`, `lfc_cor_with_aging`, `score_diff`, `score_p`, `hallmark_reversed` / `_same`, `gobp_reversed` / `_same`.
- `paper_CAS_named_genes.csv`: the paper's Fig. 2B aging genes with `lfc_*` / `padj_*` for aging, aDR, YMP (method set in `SHOW`, default RSEM).

**`compare/`:**
- `counts_per_sample.csv`: per sample, `spearman_kallisto_RSEM` (rank correlation of the gene counts) and total counts in millions per method.
- `de_agreement.csv`: kallisto (`method_a`) vs RSEM (`method_b`) per `contrast`: `total_DEG_a/b` (all DEGs), then on genes with a valid p-value in both (`genes_tested_both`): `DEG_a/b`, `shared`, `jaccard` (shared ÷ combined), `same_direction` (share of shared DEGs with the same sign), `lfc_spearman_all` / `_sig` (log2FC rank correlation over all genes / genes significant in either).
- `tpm_variation.csv`: technical (kallisto vs RSEM, same mouse) vs biological (two mice, same group and sex) variation of per-gene log2 TPM ratios, one row per kind of comparison: `n` comparisons, `max_abs_shift` (global shift removed by centring), medians over comparisons of `sd`, `median_abs` (median |log2FC|, also per expression level `TPM …`) and `over_2fold_pct` (% of genes differing > 2-fold). `offset removed` = each gene's mean method difference subtracted.
- `polyg_counts.csv`: per method and run, counts with vs without poly-G trimming: `spearman`, `counts_with_M` / `counts_without_M`, `genes_changed_5pct` / `genes_changed_2fold` (% of genes with mean ≥ 10 counts). `polyg_de_agreement.csv`: kallisto DE results with (a) vs without (b) poly-G trimming, same columns as `de_agreement.csv`.
- `method_evidence.csv`: one row per method, the criteria of [../docs/pipeline.md](../docs/pipeline.md#choosing-the-counting-method): `counted_pairs_pct_median/min/max`, `genes_detected`, `replicate_r_median`, `dispersion_aging/aDR/YMP`, `DEG_without/with` and `jaccard_outlier` (aging DEGs without / with SRR21355389), `aging_genes_up_of_5`, `aDR_genes_up_of_23` (positive controls).

**Pipeline outputs (not in git):** `kallisto/<SRR>/abundance.tsv` (`target_id`, `length`, `eff_length`, `est_counts`, `tpm`; `run_info.json` has `p_pseudoaligned`); `rsem/<SRR>/<SRR>.genes.results` (expected counts, TPM, effective length per gene).
