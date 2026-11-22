# Handover: common baseline for review

Lea's notebook (`Script.ipynb`) reduced to the setup we agreed on 20 November, so we can all use the same DEG lists and gene set results. A proposal: nothing is final until we agree. Method details are in the notebook; changes from Lea's original are in `CHANGES.md`.

## Setup
RSEM counts + gene lengths · exclude SRR21355389 (blood), keep SRR23914572 · aging old (18 + 21 mo) vs young, `~ sex + age_bin`; aDR vs AL and YMP vs PBS, `~ condition` · filter ≥ 10 reads in ≥ smallest-group size · Wald, padj < 0.05, no fold-change cutoff · apeglm-shrunken fold changes for plots · GSEA (fgsea, Wald ranking, seed 42): Hallmark main, GO BP supplementary · aging vs intervention: descriptive (ρ, % opposite, median |log2FC|), no test.

## Results
| | Aging | aDR | YMP |
|---|---|---|---|
| DEGs (up / down) | 175 (138 / 37) | 118 (46 / 72) | 50 (22 / 28) |
| Hallmark sets significant (of 50) | 15 | 16 | 24 |
| GO BP terms significant (of ~3,700) | 382 | 444 | 309 |
| Of the 15 aging Hallmark sets: significant in intervention | – | 6, all opposite | 11: 9 opposite, 2 same |
| Aging DEGs: median \|log2FC\| | 0.70 | 0.01 | 0.01 |
| Aging DEGs: ρ / % opposite | – | 0.13 / 51 % | −0.09 / 62 % |
| Aging DEGs also intervention DEGs | – | 5, all same direction | 2, both opposite |

**In one line:** both interventions oppose some aging pathways, but neither moves the aging genes themselves.

## Figures (`output/`)
`pca_all.png` (cohorts separate) · `volcano.png` (DEGs) · `heatmap_aging.png` (top 30 aging DEGs per mouse) · `gsea_hallmark.png` (15 aging sets, NES in aging / aDR / YMP) · `intervention_vs_aging.png` (aging vs intervention log2FC). Tables: `de_*.csv` (all tested genes, for importing DEG lists), `gsea_hallmark.csv`, `gsea_gobp.csv`, `intervention_vs_aging.csv`, `deg_summary.csv`.

## Limitations
3–5 mice per group · separate experiments with different sexes (aDR female, YMP male) and ages, so aDR and YMP are never compared with each other · bulk tissue: dissection (meninges, choroid plexus, white matter) can look like expression change · near-zero fold change ≠ no change.

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

GO values for the slides: table in notebook Section 4, all terms in `output/gsea_gobp.csv`.
