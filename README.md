# Mouse cerebellum: aging and rejuvenation (RNA-seq)

**Question.** Dietary restriction (aDR) and young blood plasma (YMP) are reported to rejuvenate the aged mouse brain. In the cerebellum, do they push gene expression back towards the young state, and how does that relate to the gene changes of normal aging?

**Data.** 27 bulk RNA-seq samples (one mouse each) from Hahn *et al.*, *Cell* 2023: 12 aging mice (3, 18, 21 months) and 15 mice of a 19-month rejuvenation experiment, re-processed from the raw reads ([docs/data.md](docs/data.md)).

**Approach.** Reads → QC and trimming → gene counts from STAR alignments with RSEM (GENCODE M39 / GRCm39; kallisto as sensitivity check) → DESeq2 → treatment effects compared with the aging signal ([docs/pipeline.md](docs/pipeline.md)).

## Main results (cerebellum, RSEM)
Aging = old (18 + 21 months, 7 mice) vs young (3 months, 4 mice); 176 aging DEGs. A follow-up with 21 vs 3 months isolates late aging.

| | aDR vs AL (9 females) | YMP vs PBS (6 males) |
|---|---|---|
| DEGs | 118 (46 up, 72 down) | 50 (22 up, 28 down) |
| DEGs shared with aging | 5 (1.4 expected), all **like aging** | 2 (0.5 expected), both **opposite** to aging |
| Aging-signature score, treated − control (Welch p) | +0.28 (p = 0.28) | −0.44 (p = 0.13) |
| Hallmark pathways significant in aging and in the treatment | 7: **reversed 7** | 12: **reversed 9**, same 3 |
| Late aging (21 vs 3): signature score | **+0.81** (p = 0.013), aging-like | **−0.78** (p = 0.044), reversal |

kallisto, run as a sensitivity check, gives different counts (e.g. 139 DEGs for aDR, 41 for YMP) but the same conclusions (`results/aging_vs_rejuvenation/`).

- **Aging:** microglial / interferon / complement genes rise by 18 months and stay up; metabolism (oxidative phosphorylation, mTORC1, fatty acids) falls; myelin genes fall late (at 21 months only).
- **Young plasma is partly rejuvenating:** it reverses most aging pathways it changes (interferon down, metabolism up); its genes and mouse scores lean young, significantly only against late aging.
- **Dietary restriction is mixed:** it reverses the immune and metabolic aging pathways and induces the paper's circadian genes (19 of 23 up), but deepens the late loss of myelin genes, so against late aging its mice look older.
- **Limits:** 3–5 mice per group; treatments at 19 months; aDR female-only, YMP male-only; p-values are indicative.

## Notebooks
Run in numerical order. The heavy steps run on the Euler cluster (01).

| # | Notebook | Kernel | Does |
|---|---|---|---|
| 01 | [euler_pipeline](notebooks/01_euler_pipeline.ipynb) | Python | runs the heavy steps on Euler and copies the results back; describes the cluster scripts (`euler/`) |
| 02 | [qc_trimming](notebooks/02_qc_trimming.ipynb) | Python | QC and trimming decisions with their evidence |
| 03 | [differential_expression](notebooks/03_differential_expression.ipynb) | R | what changes in each comparison: sample overview, outlier check, volcano plots for aging, aDR, YMP, annotated DEGs, 18 vs 21 months |
| 04 | [compare_methods](notebooks/04_compare_methods.ipynb) | R | sensitivity checks: does kallisto vs RSEM, or skipping poly-G trimming, change the results? |
| 05 | [aging_vs_rejuvenation](notebooks/05_aging_vs_rejuvenation.ipynb) | R | the main analysis, treatments vs aging: overlap, signature score, pathways, biological themes |

## Reproduce
1. **Software:** laptop `pixi install` (exact versions in `pixi.lock`); Euler `bash euler/00_setup.sh` (tool versions pinned there).
2. **Data and reference:** `euler/01_download.sh` fetches the 27 runs from NCBI SRA (`metadata/srr_list.txt`) and GENCODE M39.
3. **Processing:** `euler/submit.sh` → `euler/05_check.sh` → `euler/04_collect.sh`, then copy the tarballs to `results/euler/` (commands in notebook 01).
4. **Analysis:** notebooks 01–05 in order. Pathway tests are seeded; MSigDB mouse-native gene sets (msigdbr, db_species = MM) are fetched on first use, so a later MSigDB release can change pathway results.

The scripts are plain bash and SLURM and should run on any SLURM cluster (not tested elsewhere). To skip steps 2–3, use the tarballs in `results/euler/` (checksums: `metadata/data_checksums.txt`). A few one-run tests quoted in `docs/pipeline.md` (marked there) are not scripted.

```bash
pixi install                              # R/Bioconductor, fastp, Jupyter, ... (pixi.toml)
pixi run R -e 'IRkernel::installspec()'   # once: makes the R kernel visible to VS Code / Jupyter
pixi run jupyter lab                      # or open notebooks/ in VS Code
brew install sratoolkit                   # only for the adapter check in notebook 02
```
Notebook 01 needs an ETH account and the VPN. Finished work is skipped, so interrupted runs resume.

## Data for your own analysis
Start from [for_teammates/](for_teammates/README.md): RSEM expression tables, a sample table (age, sex, group) and a gene-name table (about 10 MB). RSEM's per-run files (input for `tximport`, which sums them to gene level) are `results/euler/rsem_results.tar.gz` (72 MB). The kallisto tarballs (626 MB, and 312 MB for the poly-G check; over GitHub's file limit) are only on the laptop that ran the pipeline. Not available: the reads (54 GB) and the alignments (BAM files were not kept).

## Files
- **Docs:** [docs/data.md](docs/data.md) (study, design), [docs/pipeline.md](docs/pipeline.md) (each step and why, method choice), [docs/glossary.md](docs/glossary.md), [results/README.md](results/README.md) (how to read the result tables), [docs/Tutorial.md](docs/Tutorial.md) (the course tutorial this follows, Z. He, ETH Zurich).
- **Code:** `notebooks/`, `scripts/`, `R/helpers.R`, `euler/`, `pixi.toml` / `pixi.lock`.
- **Data:** `for_teammates/`, `metadata/`, `results/`. Not in git: `raw_data/`, `results/kallisto*/`, `results/rsem*/`, the kallisto tarballs in `results/euler/`, FastQC reports.

Paper: Hahn *et al.* 2023, *Cell* 186:4117 ([doi](https://doi.org/10.1016/j.cell.2023.07.027)); authors' code: github.com/OliInTheValley/SpatioTemporal_Analysis.
