# Handover: common baseline for review

**What this is:** Lea's notebook (`Script.ipynb`), reduced to the DESeq2 setup we agreed on 20 November 2026, so everyone can import the same DEG lists and gene set results (as Lea suggested). It is a proposal: please review it, and nothing is final until we agree. Changes from Lea's original: `CHANGES.md`.

**Run:** from this folder, `pixi run jupyter nbconvert --to notebook --execute --inplace Script.ipynb` (a few minutes; the first run downloads MSigDB 2026.1 once). Input: `../../for_teammates/`. All outputs go to `output/`.

## Agreed setup (as implemented)
| Choice | Setting |
|---|---|
| Input | RSEM expected counts + effective lengths (`DESeqDataSetFromTximport`) |
| Samples | exclude SRR21355389 (blood contamination: haemoglobin 11–15× above any other aging sample); keep SRR23914572 |
| Comparisons | aging: old (18 + 21 mo, n = 7) vs young (3 mo, n = 4), `~ sex + age_bin`; aDR (5) vs AL (4), `~ condition`; YMP (3) vs PBS (3), `~ condition` |
| Gene filter | ≥ 10 reads in at least the smallest-group size (4 / 4 / 3) |
| Test and DEG | Wald test, `results(alpha = 0.05)`, DEG = padj < 0.05, no fold-change cutoff |
| Fold changes in plots and comparisons | apeglm-shrunken |
| Gene sets | mouse MSigDB Hallmark, fgsea on the Wald statistic (all tested genes), sizes 15–500, seed 42 |
| Aging vs intervention | the 175 aging DEGs: scatter of shrunken log2FC, Spearman ρ, % opposite sign, median \|log2FC\|; descriptive, no test |

**On Lea's question ("are these 3 comparisons relevant and the tests right?"):** yes for the three comparisons (aging, aDR vs AL, YMP vs PBS, each within its cohort). Tests: the Wald test for genes and fgsea for gene sets are the only tests. The aging-vs-intervention comparison uses descriptive numbers only, because aging genes are correlated, so a binomial or correlation p-value would overstate the evidence.

## Main results
| | Aging | aDR | YMP |
|---|---|---|---|
| DEGs (up / down) | 175 (138 / 37) | 118 (46 / 72) | 50 (22 / 28) |
| Hallmark sets significant (of 50) | 15 | 16 | 24 (22 up) |
| The 15 aging sets, also significant in the intervention | – | 6, all opposite to aging | 11: 9 opposite, 2 same |
| Aging DEGs: median \|shrunken log2FC\| | 0.70 (with age) | 0.01 | 0.01 |
| Aging DEGs: Spearman ρ / % opposite sign | – | 0.13 / 51 % | −0.09 / 62 % |
| Aging DEGs also DEG in the intervention | – | 5, all same direction as aging | 2, both opposite |

- **Aging:** neuroinflammation (complement, activated microglia, MHC class I, cathepsins, *Gfap*); interferon and inflammatory pathways up; oxidative phosphorylation, mTORC1 and other metabolic pathways down.
- **aDR:** glucocorticoid, circadian and immediate early genes; at pathway level opposite to aging (complement, coagulation, KRAS down; OXPHOS, MYC, TNFα up).
- **YMP:** few DEGs, partly meningeal genes (*Ptgds, Ogn, Slc47a1*); interferon down and metabolism up, within a general upward shift of pathways.
- **Conclusion:** both interventions oppose some aging pathways, but neither reverses the individual aging genes (they barely move: median change 0.01 vs 0.70 with age).

## Five figures (`output/`)
| File | Shows |
|---|---|
| `pca_all.png` | samples separate by cohort (PC1 41 %), so each comparison stays within its cohort |
| `volcano.png` | DEGs per comparison (shrunken log2FC vs −log10 padj) |
| `heatmap_aging.png` | top 30 aging DEGs per mouse: young low, 18 and 21 months alike |
| `gsea_hallmark.png` | the 15 Hallmark sets significant in aging, NES in aging / aDR / YMP; non-significant tiles faded |
| `intervention_vs_aging.png` | aging vs intervention log2FC for the 175 aging DEGs; genes significant in both highlighted |

Tables: `de_aging.csv`, `de_aDR.csv`, `de_YMP.csv` (all tested genes; for importing DEG lists), `deg_summary.csv`, `gsea_hallmark.csv` (all 50 sets × 3), `intervention_vs_aging.csv`.

## Limitations
- 3–5 mice per group; single-gene fold changes from 3 vs 3 mice are uncertain.
- Separate experiments: different sequencing runs, sexes (aDR female, YMP male) and ages (19 months). aDR and YMP cannot be compared with each other; sex-specific aging and batch can distort the comparison with aging.
- Bulk tissue: dissection differences (meninges, choroid plexus, white matter content) can appear as expression changes. *Cdr1* in aDR is on/off.
- A shrunken log2FC near 0 means little evidence of change, not proof of none.

## Decisions still open
| Decision | Options | Suggestion |
|---|---|---|
| Use this notebook's outputs as the common source for DEG lists and gene sets? | yes / adapt / another notebook | yes, after Lea's and Nick's review |
| GO BP | drop, or add as a supplementary table with the same ranking and seed | supplementary table only; Hallmark stays the main figure |
| DEG definition on the slides | agreed (no fold-change cutoff) or label a 1.5-fold subset explicitly | use the agreed counts everywhere |
| Myelin / oligodendrocyte gene sets (GO BP) | include or not | if included, say "myelin-related genes higher / lower", not "more myelin"; not changed with aging here |
| 16 h fasting before sacrifice (aDR slide) | cite the paper's methods, or remove | needs a source |
| Literature support for interpretations | add citations per claim | add before finalising |
| Who presents what | | to agree |

## Presentation check (Nick's slides 4–8 and the combined GSEA figure, as shared in chat)
| Slide | Current | Needs updating to match this baseline |
|---|---|---|
| 4 Aging overview | 134 DEGs (117 / 17), 1.5-fold lines in volcano; 18 vs 3: 89, 21 vs 3: 74, 21 vs 18: 0 | 175 (138 / 37), no fold-change lines; recompute the per-age counts with the same setup; state the reference for "approx. increase" (we report 11–15× vs the highest other aging sample); "Abreviations" typo |
| Combined GSEA figure (aging / diet / plasma) | Hallmark + GO mixed; reversal arrows | aDR "inflammatory response reversed" not significant (padj 0.07); oligodendrocyte differentiation not significant in aging (padj 0.6), so the YMP arrow goes; label GO vs Hallmark |
| 5 Aging activates immune genes | renamed sets; "driven by complement and microglia" | holds: interferon, complement, allograft, inflammatory, antigen presentation, synapse pruning, macrophage activation. Astrocyte activation not significant (0.09). Give exact term names; "consistent with" instead of "driven by" |
| 6 Aging shifts metabolism | sterol biosynthesis −1.53 | OXPHOS, mTORC1, fatty acid hold; sterol biosynthesis not significant (0.2); typo "Fatty acid cholesterol" |
| 7 Diet restriction | 53 DEGs (16 / 37); fasting; starvation, angiogenesis | 118 (46 / 72); circadian, OXPHOS, myelin terms hold; response to starvation (0.2), angiogenesis (0.05), sprouting angiogenesis (0.2) not significant; fasting needs a source |
| 8 Young plasma | 42 DEGs (15 / 27); "more myelin" | 50 (22 / 28); interferon down and metabolism up hold; response to type I interferon (0.09) and EMT (0.3) not significant; "higher myelin-related genes" |

Values for the check come from `../../results/common/` (same setup; GO BP computed there with the same ranking and seed).
