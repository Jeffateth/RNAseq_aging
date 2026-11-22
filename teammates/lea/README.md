# Lea's analysis on the team setup (proposed baseline)

Lea's notebook ([lhasch/System-Genomics-Project](https://github.com/lhasch/System-Genomics-Project) @ `8a133ad`) reduced to the DESeq2 setup we agreed on 20 November, so we can all use the same DEG lists and gene set results. A proposal for review: nothing is final until we agree. Method details are in the notebook.

**Run:** `pixi run jupyter nbconvert --to notebook --execute --inplace Script.ipynb` (a few minutes; outputs in `output/`). `Script.R` is the same as a plain script.

## Setup
RSEM counts + gene lengths · exclude SRR21355389 (blood), keep SRR23914572 · aging old (18 + 21 mo) vs young, `~ sex + age_bin`; aDR vs AL and YMP vs PBS, `~ condition` · filter ≥ 10 reads in ≥ smallest-group size · Wald, padj < 0.05, no fold-change cutoff · apeglm-shrunken fold changes for plots · GSEA (fgsea, Wald ranking, seed 42): Hallmark main, GO BP supplementary · aging vs intervention: descriptive (ρ, % opposite, median |log2FC|), no test.

## Results
| | Aging | aDR | YMP | Lea's original (aging / aDR / YMP) |
|---|---|---|---|---|
| DEGs (up / down) | 175 (138 / 37) | 118 (46 / 72) | 50 (22 / 28) | 219 / 171 / 37 |
| Hallmark sets significant (of 50) | 15 | 16 | 24 | 17 / 20 / 23 |
| GO BP terms significant (of ~3,700) | 382 | 444 | 309 | |
| Of the 15 aging Hallmark sets: significant in intervention | – | 6, all opposite | 11: 9 opposite, 2 same | |
| Aging DEGs: median \|log2FC\| | 0.70 | 0.01 | 0.01 | |
| Aging DEGs: ρ / % opposite | – | 0.13 / 51 % | −0.09 / 62 % | ρ −0.03 / −0.20 (raw log2FC) |
| Aging DEGs also intervention DEGs | – | 5, all same direction | 2, both opposite | |

**In one line:** both interventions oppose some aging pathways, but neither moves the aging genes themselves.

**Figures (`output/`):** `pca_all.png` (cohorts separate) · `volcano.png` · `heatmap_aging.png` (top 30 aging DEGs) · `gsea_hallmark.png` (15 aging sets in aging / aDR / YMP) · `intervention_vs_aging.png`. **Tables:** `de_*.csv` (all tested genes, for importing DEG lists), `gsea_hallmark.csv`, `gsea_gobp.csv`, `intervention_vs_aging.csv`, `deg_summary.csv`.

**Limitations:** 3–5 mice per group · separate experiments with different sexes (aDR female, YMP male) and ages, so aDR and YMP are never compared with each other · bulk tissue: dissection (meninges, choroid plexus, white matter) can look like expression change · near-zero fold change ≠ no change.

## Changes from Lea's original
- **Changed to the agreed setup:** gene lengths (`DESeqDataSetFromTximport` instead of `FromMatrix`), sex term for aging, filter by smallest-group size, DEG = padj only (no |log2FC| > 1 in volcanoes), one heatmap on VST instead of log TPM.
- **Kept:** sample exclusion and haemoglobin table; pooled 18 + 21 months (her 21 vs 18 result as the reason); GSEA settings (identical); from 5b: shrunken fold changes, reversed set = significant in both with opposite NES, a gene counts only if also an intervention DEG; from 5a: the scatter, ρ, % opposite; the aging category table.
- **Removed:** histograms, correlation trees, extra heatmaps, GO bar plots (GO BP kept as a table) · 5a binomial / correlation p-values and GO of reversed genes (genes not independent, wrong background, as Lea noted) · 5b camera / GSEA of the aging gene sets (valid, but an extra method; the scatter answers the same question) · 5b controls (listed as limitations instead) · install cell, `setwd()`, cache. All restorable from git history (commit `9e240d7`).

## Open decisions
| Decision | Suggestion |
|---|---|
| Use these outputs as the shared DEG lists and gene sets? | yes, after review |
| DEG counts on slides | agreed counts everywhere (no 1.5-fold cutoff) |
| Myelin results | "myelin-related genes higher / lower", not "more myelin"; not changed with age |
| 16 h fasting before sacrifice (aDR slide) | cite the paper's methods or remove |
| Literature for interpretations | add before finalising |
| Who presents what | to agree |

## Slide check (Nick's slides 4–8 and the GSEA bubble figure)
| Slide | Change needed |
|---|---|
| 4 Aging overview | 134 → 175 DEGs (138 / 37), remove fold-change lines; recompute per-age counts with this setup; define the "approx. increase" reference; "Abbreviations" |
| GSEA bubble figure | remove arrows for aDR inflammatory response (padj 0.07) and oligodendrocyte differentiation in aging (0.6); label GO vs Hallmark |
| 5 Immune | holds; astrocyte activation not significant (0.09); exact term names; "consistent with" instead of "driven by" |
| 6 Metabolism | holds; sterol biosynthesis not significant (0.2) |
| 7 Diet | 53 → 118 DEGs (46 / 72); response to starvation, angiogenesis, sprouting angiogenesis not significant; fasting needs a source |
| 8 Young plasma | 42 → 50 DEGs (22 / 28); type I interferon and EMT not significant; "higher myelin-related genes" |

GO values for the slides: table in notebook Section 4; all terms in `output/gsea_gobp.csv`.
