# Aging and rejuvenating interventions: summary

Refactored version of Lea's notebook (`Script.ipynb`), run with the team's common DESeq2 setup (agreed 20 November 2026). All numbers below come from `output/`.

## Main findings
| Comparison | DEGs (up / down) | Hallmark sets (padj < 0.05) | vs aging: genes (175 aging DEGs; median |log2FC| with age 0.70) | vs aging: pathways (15 aging sets) |
|---|---|---|---|---|
| Aging, old (18 + 21 mo) vs young (3 mo) | 175 (138 / 37) | 15 | – | – |
| aDR vs AL | 118 (46 / 72) | 16 | median |log2FC| 0.01; ρ = 0.13, 51 % opposite; 5 DEGs in both, all same direction as aging | 6 significant in both, all opposite |
| YMP vs PBS | 50 (22 / 28) | 24 (22 up) | median |log2FC| 0.01; ρ = −0.09, 62 % opposite; 2 DEGs in both, both opposite | 11 significant in both: 9 opposite, 2 same |

- **Aging:** neuroinflammation: complement (*C4b, C1q*), activated microglia (*Clec7a, Cst7, Trem2, Tyrobp, Cd68*), MHC class I (*H2-K1, H2-Q6, B2m*), cathepsins, *Gfap*. Interferon and inflammatory pathways go up; oxidative phosphorylation, mTORC1 and other metabolic pathways go down.
- **aDR:** small effects (glucocorticoid, circadian, immediate early genes). At the pathway level it opposes aging (complement, coagulation and KRAS down; oxidative phosphorylation, MYC and TNFα/NF-κB up). At the gene level it does not oppose aging.
- **YMP:** few DEGs, several of them meningeal genes (*Ptgds, Ogn, Slc47a1*), which may reflect dissection differences. Most of its pathways go up, so its "opposite to aging" pathway calls partly reflect a general upward shift.
- **Neither intervention** measurably changes the aging genes: their median shrunken log2FC is 0.01 after aDR or YMP, against 0.70 with age. With changes this small, the % opposite sign is largely noise.
- **Gene level vs pathway level:** both interventions shift some aging-related Hallmark pathways in the opposite direction. Pathways summarise small coordinated shifts across many genes, mostly not aging DEGs, so pathway-level opposition does **not** mean the aging genes themselves are reversed.

## Limitations
- 3–5 mice per group. Single-gene fold changes from 3 vs 3 mice are uncertain.
- The cohorts are separate experiments with different sequencing runs, sexes (aging both, aDR female, YMP male) and ages (intervention mice are 19 months). aDR and YMP are therefore each compared with aging only, never with each other. Batch differences and sex-specific aging changes can weaken or distort that comparison.
- Bulk tissue: dissection differences (meninges, choroid plexus in SRR23914572) can appear as expression changes. *Cdr1* in aDR is on/off (0 reads in all AL mice).
- The aging-vs-intervention comparison is descriptive. A % opposite above 50 % is not tested for significance, because aging genes are correlated, not independent.

## Methodological decisions
| Decision | Reason |
|---|---|
| RSEM expected counts + effective lengths, `DESeqDataSetFromTximport` | corrects for gene-length differences between samples |
| Exclude SRR21355389 | haemoglobin genes 11–15× above every other aging sample (blood contamination) |
| Keep SRR23914572 | choroid plexus genes, not blood (team decision) |
| Aging: 18 + 21 months pooled, `~ sex + age_bin`; aDR, YMP: `~ condition` | team setup; sex is the main axis of variation within the aging cohort |
| ≥ 10 reads in the smallest-group size (4 / 4 / 3); Wald; `alpha = 0.05`; DEG = padj < 0.05, no fold-change cutoff | team setup |
| apeglm-shrunken log2FC for plots and comparisons | stable fold changes for low-count genes |
| Mouse MSigDB Hallmark (release 2026.1), fgsea on the Wald statistic, seed 42 | all tested genes ranked, not only DEGs |
| No significance test of the % opposite | aging genes are correlated, not independent |

Removed from the original notebook: correlation clustering, expression histograms, the 21 vs 18 month analysis, heatmaps for aDR and YMP, GO BP, camera, permutation and binomial tests, alternative aging signatures, cell-type scores, and GO of reversed genes. They remain available in git history (commit `9d815b9`).

## Outputs (`output/`)
| File | Content |
|---|---|
| `de_aging.csv`, `de_aDR.csv`, `de_YMP.csv` | all tested genes: `log2FoldChange`, `lfc_shrunk`, `stat`, `pvalue`, `padj` |
| `deg_summary.csv` | tested genes, DEGs (up / down) per comparison, with the reference counts |
| `gsea_hallmark.csv` | Hallmark GSEA per comparison: `NES`, `pval`, `padj`, `size` |
| `intervention_vs_aging.csv` | per intervention: Spearman ρ, % opposite, median |log2FC| (aging, intervention), DEGs in both (opposite / same) |
| `pca_all.png`, `pca_within.png`, `volcano.png`, `heatmap_aging.png`, `gsea_hallmark.png`, `intervention_vs_aging.png` | figures |
