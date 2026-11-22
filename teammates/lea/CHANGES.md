# What changed from Lea's original notebook

Original: Lea's `Script.ipynb` at commit `8a133ad` of her repository (copy in this repo's git history, commit `9e240d7`). Refactored: `Script.ipynb` in this folder. The refactor is a proposal for a common baseline, not a judgement that the original was wrong. Several of Lea's choices were already the agreed ones and are kept unchanged.

## Results side by side
| | Original | Refactored (agreed setup) | Main reason for the difference |
|---|---|---|---|
| Aging DEGs | 219 | 175 (138 up / 37 down) | sex term added; gene-length correction; filter by smallest group |
| aDR DEGs | 171 | 118 (46 / 72) | gene-length correction; filter ≥ 10 reads in ≥ 4 samples |
| YMP DEGs | 37 | 50 (22 / 28) | gene-length correction |
| Hallmark sets (padj < 0.05) | 17 / 20 / 23 | 15 / 16 / 24 | follows from the DE changes (ranking and seed unchanged) |
| Aging vs intervention, ρ (aDR / YMP) | −0.03 / −0.20 (raw log2FC, 216 genes) | 0.13 / −0.09 (shrunken log2FC, 175 genes) | shrunken fold changes, new aging gene list |

## Retained
| Element | Notes |
|---|---|
| Exclusion of SRR21355389 with the haemoglobin table | same genes; table now also shows fold vs the other aging samples |
| Within-cohort comparisons, pooled 18 + 21 months | unchanged |
| Wald test, `alpha = 0.05`, DEG = padj < 0.05 | unchanged |
| GSEA: mouse-native MSigDB Hallmark and GO BP, fgsea on the Wald statistic, sets 15–500, seed 42 | identical to Lea's setup; GO BP now a supplementary table, not figures |
| apeglm-shrunken fold changes (from 5b step 1) | now also used in the volcano plots |
| 5b's definition of a reversed pathway (significant in both, opposite NES sign) | used for the pathway count in Section 4 |
| 5b's rule that a gene counts as reversed only if it is also an intervention DEG | used to highlight genes in the Section 5 scatter |
| 5a's scatter of aging vs intervention log2FC, Spearman ρ, % opposite sign | kept as descriptive numbers |
| Aging signature category table (microglia, complement, MHC, cathepsins, *Gfap*) | kept in Section 3, checked against the new DEG list |
| 21 vs 18 month result (16 DEGs) | kept as a one-line justification for pooling |

## Changed to match the agreed setup
| Original | Refactored | Why |
|---|---|---|
| `DESeqDataSetFromMatrix` on rounded counts | `DESeqDataSetFromTximport` with RSEM effective lengths | corrects for per-sample gene-length differences |
| Aging design `~ age_bin` | `~ sex + age_bin` | the aging cohort has both sexes; the term keeps sex differences out of the age test |
| Gene filter ≥ 10 reads in ≥ 3 samples | ≥ 10 reads in ≥ smallest-group size (4 / 4 / 3) | a gene can then be expressed in a whole group |
| Volcano colour: padj < 0.05 and \|log2FC\| > 1 | padj < 0.05 only | same DEG definition everywhere (no fold-change cutoff) |
| Heatmaps on log TPM | heatmap on VST counts | same normalisation as DESeq2 |

## Simplified or removed
| Removed | Why |
|---|---|
| Expression histograms, Pearson / Spearman clustering trees | the PCA shows the same structure; not needed for the research question |
| Per-cohort PCA | the sex term is explained in one sentence instead |
| Heatmaps for aDR and YMP, 21 vs 18 volcano | one aging heatmap is enough; 21 vs 18 reduced to the pooling sentence |
| GO BP bar plots | Hallmark's 50 sets are the main figure; GO BP GSEA itself is kept as a supplementary table (`output/gsea_gobp.csv`) |
| 5a binomial test and correlation p-values | aging genes are not independent, so the p-values are too small (Lea's own note in 5a) |
| 5a GO enrichment of reversed genes | wrong background (Lea's own note); not needed |
| 5b step 2: camera and GSEA of the aging_up / aging_down sets | valid, but an extra method to justify; the scatter, ρ and % opposite answer the same question descriptively |
| 5b controls: random gene sets, alternative aging signatures, cell-type marker scores | each needs its own justification; sex and cell-type composition are listed as limitations instead |
| Package installation cell, local `setwd()`, cache helper | the repo's pixi environment and relative paths make it run anywhere; the whole notebook takes a few minutes |

Everything removed is still in git history (commits `9e240d7` and `9d815b9`) and can be restored if the group decides it is needed.
