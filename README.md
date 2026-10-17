# Mouse cerebellum: aging and rejuvenation (RNA-seq)

**Question.** Dietary restriction (aDR) and young blood plasma (YMP) are reported to rejuvenate the aged mouse brain. In the cerebellum, do they push gene expression back towards the young state, and how does that relate to the gene changes of normal aging?

**Data.** 27 bulk RNA-seq samples (one mouse each) from Hahn *et al.*, *Cell* 2023: 12 aging mice (3, 18, 21 months) and 15 mice of a 19-month rejuvenation experiment, re-processed from the raw reads ([docs/data.md](docs/data.md)).

**Approach.** reads → QC and trimming → gene counts by three counting methods (kallisto; STAR; RSEM on STAR's alignments), GENCODE M39 / GRCm39 → DESeq2 → treatment effects compared with the aging signal ([docs/pipeline.md](docs/pipeline.md)).

## Main results (cerebellum)
The three methods are shown as kallisto / STAR / RSEM; they agree in direction everywhere. Which one is primary is decided on the evidence in [docs/pipeline.md](docs/pipeline.md#choosing-the-counting-method) (pending).

| | aDR vs AL (9 females) | YMP vs PBS (6 males) |
|---|---|---|
| Differentially expressed genes (DEGs) | 139 (64 up) / 185 (70 up) / 118 (46 up) | 41 (17 up) / 26 (9 up) / 50 (22 up) |
| DEGs shared with aging (21 vs 3 months: 134 / 168 / 116 DEGs) | 6 / 13 / 5 (1.0 / 2.0 / 1.0 expected by chance); **like aging** 6 / 12 / 5, opposite 0 / 1 / 0 | 2 / 3 / 2 (0.3 expected); all **opposite** to aging |
| log2 fold-change correlation with aging, among aging DEGs | **+0.39 / +0.48 / +0.44** (aging-like) | **−0.48 / −0.41 / −0.33** (reversal) |
| Aging-signature score, treated minus control (Welch p) | **+0.73 / +0.77 / +0.82** (p = 0.016 / 0.018 / 0.011) | **−0.79 / −0.73 / −0.82** (p = 0.004 / 0.048 / 0.043) |
| Hallmark pathways significant in aging and in the treatment | 8 / 7 / 7: reversed 5 / 4 / 5, same 3 / 3 / 2 | 14 / 15 / 11: **reversed 11 / 13 / 10**, same 3 / 2 / 1 |

- **Young plasma reverses part of the cerebellar aging signal:** interferon-response pathways (up in aging) go down; mTORC1, fatty-acid and cholesterol metabolism, p53, hypoxia and oxidative phosphorylation (down in aging) go up.
- **Dietary restriction does not:** it induces the circadian-clock genes the paper reports (19 of 23 up in all three methods; Cry1, Dbp, Ciart, Cirbp) and reverses some pathways, but not the interferon response, and its aging-signature score rises.
- **Limits:** 3–5 mice per group; treatments at 19 months vs aging at 21; aDR female-only, YMP male-only; p-values are indicative.

## Notebooks
Run in numerical order. The heavy steps (QC, trimming, kallisto, STAR, RSEM) run on the Euler cluster (01).

| # | Notebook | Kernel | Does |
|---|---|---|---|
| 00 | [acquire_data](notebooks/00_acquire_data.ipynb) | Python | *optional:* downloads the 27 runs (~54 GB) to the laptop, for local runs and the adapter check |
| 01 | [euler_pipeline](notebooks/01_euler_pipeline.ipynb) | Python | **runs the heavy steps on Euler** and copies the results back; the cluster commands and scripts (`euler/`) are described here |
| 02 | [quantify_kallisto](notebooks/02_quantify_kallisto.ipynb) | Python | QC and trimming decisions, kallisto summary (`RUN_LOCALLY = True` runs the steps on the laptop: hours) |
| 03 | [differential_expression](notebooks/03_differential_expression.ipynb) | R | outlier check, DESeq2 for aging, aDR, YMP (method set in `METHOD`) |
| 04 | [compare_methods](notebooks/04_compare_methods.ipynb) | R | agreement of the three methods and the evidence for choosing one |
| 05 | [aging_vs_rejuvenation](notebooks/05_aging_vs_rejuvenation.ipynb) | R | the main analysis: overlap, signature score, pathways |

## Run it
```bash
pixi install                              # kallisto, fastp, FastQC, R/Bioconductor, Jupyter (pixi.toml)
pixi run R -e 'IRkernel::installspec()'   # once: makes the R kernel visible to VS Code / Jupyter
pixi run jupyter lab                      # or open notebooks/ in VS Code and pick the kernel named above
brew install sratoolkit                   # only for notebook 00 / local runs
```
Notebook 01 needs an ETH account and the VPN. Finished work is always skipped, so interrupted runs resume.

## Files
- **Docs:** [docs/data.md](docs/data.md) (study, design, identifiers), [docs/pipeline.md](docs/pipeline.md) (each step and why, method choice), [docs/glossary.md](docs/glossary.md), [results/README.md](results/README.md) (how to read the result tables), [docs/Tutorial.md](docs/Tutorial.md) (the course tutorial this follows, Z. He, ETH Zurich).
- **Code:** `notebooks/`, `scripts/` (`helpers.py`, `build_metadata.py`), `R/helpers.R`, `euler/` (SLURM scripts), `pixi.toml` / `pixi.lock` (software).
- **Data:** `metadata/` (sample table, NCBI source tables), `results/` (small outputs in git).
- **Not in git (regenerated):** `raw_data/`, `reference/`, `results/kallisto/`, `results/rsem/`, FastQC reports, STAR alignments. Local and unused: `processed_data/` (the authors' GEO count tables).

Paper: Hahn *et al.* 2023, *Cell* 186:4117 ([doi](https://doi.org/10.1016/j.cell.2023.07.027)); authors' code: github.com/OliInTheValley/SpatioTemporal_Analysis.
