# Changes from Lea's original notebook

Original: Lea's `Script.ipynb` at `8a133ad` (copy in git history, commit `9e240d7`). Much of it is kept as it was.

| | Original | Refactored |
|---|---|---|
| DEGs (aging / aDR / YMP) | 219 / 171 / 37 | 175 / 118 / 50 |
| Hallmark sets significant | 17 / 20 / 23 | 15 / 16 / 24 |
| ρ aging vs aDR / YMP | −0.03 / −0.20 (raw log2FC) | 0.13 / −0.09 (shrunken) |

## Kept
- Exclusion of SRR21355389 with the haemoglobin table; within-cohort comparisons; pooled 18 + 21 months
- Wald test, padj < 0.05, `alpha = 0.05`
- GSEA: mouse Hallmark and GO BP, fgsea on the Wald statistic, sizes 15–500, seed 42
- From 5b: shrunken fold changes, reversed set = significant in both with opposite NES, a gene counts only if it is also an intervention DEG
- From 5a: the aging vs intervention scatter, ρ, % opposite
- The aging signature category table; the 21 vs 18 result as reason for pooling

## Changed (agreed setup)
| Original | Now | Why |
|---|---|---|
| `DESeqDataSetFromMatrix` | `DESeqDataSetFromTximport` with gene lengths | corrects per-sample gene length |
| aging `~ age_bin` | `~ sex + age_bin` | both sexes in the aging cohort |
| ≥ 10 reads in ≥ 3 samples | ≥ smallest-group size (4 / 4 / 3) | agreed filter |
| volcano: padj and \|log2FC\| > 1 | padj only | one DEG definition |
| heatmaps on log TPM | one heatmap on VST | same normalisation as DESeq2 |
| GO BP bar plots | GO BP table only | Hallmark is the main figure |

## Removed
| Removed | Why |
|---|---|
| Histograms, correlation trees, extra heatmaps, 21 vs 18 volcano | PCA and one heatmap cover the question |
| 5a binomial and correlation p-values, GO of reversed genes | genes are not independent; wrong background (both noted by Lea) |
| 5b camera / GSEA of the aging gene sets | valid but an extra method; the scatter answers the same question descriptively |
| 5b controls (random sets, other signatures, cell-type scores) | each needs its own justification; listed as limitations instead |
| Install cell, `setwd()`, cache | pixi environment and relative paths; runs in a few minutes |

Everything removed can be restored from git history.
