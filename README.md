# Mouse cerebellum: aging and rejuvenation (RNA-seq)

**Question.** Dietary restriction (aDR) and young blood plasma (YMP) are reported to rejuvenate the aged mouse brain. In the cerebellum, do they push gene expression back towards the young state, and how does that relate to the gene changes of normal aging?

**Data.** 27 bulk RNA-seq samples (one mouse each) from Hahn *et al.*, *Cell* 2023: 12 aging mice (3, 18, 21 months) and 15 mice of a 19-month rejuvenation experiment. Raw reads are re-processed from scratch ([docs/data.md](docs/data.md)).

**Approach.** reads → quality control and trimming → gene counts by two independent methods (kallisto, STAR) → differential expression with DESeq2 → comparison of treatment effects with the aging signal ([docs/pipeline.md](docs/pipeline.md)).

## Main results (cerebellum)
The direction of every result is the same with both counting methods (kallisto numbers shown).

| | aDR vs AL (9 females) | YMP vs PBS (6 males) |
|---|---|---|
| Differentially expressed genes (DEGs) | 142 (65 up, 77 down) | 38 (19 up, 19 down) |
| DEGs shared with aging (21 vs 3 months: 131 DEGs, 123 testable here) | 7 (1.1 expected by chance); 6 move **like aging**, 1 opposite | 3 (0.2 expected); all **opposite** to aging |
| log2 fold-change correlation with aging, among aging DEGs | **+0.42** (aging-like) | **−0.38** (reversal) |
| Aging-signature score, treated minus control | **+0.68** (p = 0.005) | **−0.80** (p = 0.008) |
| Pathways significant in aging and in the treatment (Hallmark) | 8: 5 reversed, 3 same | 12: **11 reversed**, 1 same |

- **Young plasma reverses part of the cerebellar aging signal:** interferon-response pathways (up in aging) go down; mTORC1, fatty-acid and cholesterol metabolism, p53, hypoxia and oxidative phosphorylation (down in aging) go up.
- **Dietary restriction does not:** it induces the circadian-clock genes the paper reports (19 of 23 up; Cry1, Dbp, Ciart, Cirbp) and reverses some pathways, but not the interferon response, and its aging-signature score rises.
- **Limits:** 3–5 mice per group; treatments at 19 months vs aging at 21; aDR female-only, YMP male-only; p-values are indicative.

## Workflow
Run the notebooks in numerical order. The heavy steps (QC, trimming, kallisto, STAR) run on the Euler cluster, so notebook 01 comes first.
```
NCBI SRA --> .sra --> FASTQ --> FastQC, fastp --> kallisto --+   QC, trimming, kallisto and STAR
 (27 runs)       \---------> fastp --> STAR ----------------+   run on the Euler cluster (01)
                                                            v
        02 QC/trimming decisions --> 03 DESeq2 per cohort --> 04 kallisto vs STAR --> 05 aging vs rejuvenation
```

| # | Notebook | Kernel | Does | Time |
|---|---|---|---|---|
| 00 | [acquire_data](notebooks/00_acquire_data.ipynb) | Python | *optional:* downloads the 27 runs (~54 GB) to the laptop, for local runs and the adapter check (Euler downloads its own copy) | 30–60 min |
| 01 | [euler_pipeline](notebooks/01_euler_pipeline.ipynb) | Python | **runs QC, trimming (adapter and poly-G), kallisto and STAR on the ETH cluster (Euler)**; copies the results back | STAR array ~2 h (kallisto/QC array not yet timed); the laptop stays idle |
| 02 | [quantify_kallisto](notebooks/02_quantify_kallisto.ipynb) | Python | QC and trimming decisions, kallisto summary; **reads** the Euler results (set `RUN_LOCALLY = True` to run the steps on the laptop: hours) | minutes |
| 03 | [differential_expression](notebooks/03_differential_expression.ipynb) | R | outlier check, DESeq2 for aging, aDR, YMP | 1 min |
| 04 | [compare_kallisto_star](notebooks/04_compare_kallisto_star.ipynb) | R | agreement of the two methods | 3 min |
| 05 | [aging_vs_rejuvenation](notebooks/05_aging_vs_rejuvenation.ipynb) | R | the main analysis: overlap, signature score, pathways | 4 min |

STAR counts are in the repo (`results/star/counts_unstranded.csv.gz`), but kallisto output is not (650 MB): copy it back from Euler (01).

## Run it
```bash
pixi install                              # kallisto, fastp, FastQC, R/Bioconductor, Jupyter (pixi.toml)
pixi run R -e 'IRkernel::installspec()'   # once: makes the R kernel visible to VS Code / Jupyter
pixi run jupyter lab                      # or open notebooks/ in VS Code; pick the kernel named in the table
brew install sratoolkit                   # only for notebook 00 / local runs (prefetch, fasterq-dump)
```
Notebook 01 needs an ETH account and the VPN ([euler/README.md](euler/README.md)). Finished work is always skipped, so interrupted runs resume.

## Documentation
| File | Content |
|---|---|
| [docs/data.md](docs/data.md) | what the data are: study, design, samples, identifiers |
| [docs/pipeline.md](docs/pipeline.md) | the pipeline step by step, with the decision and alternatives at each step |
| [docs/glossary.md](docs/glossary.md) | terms and abbreviations |
| [results/README.md](results/README.md) | how to read the result tables (columns, direction, significance) |
| [euler/README.md](euler/README.md) | the cluster scripts |
| [docs/Tutorial.md](docs/Tutorial.md) | the course tutorial the workflow follows (Z. He, ETH Zurich) |

## Layout
```
notebooks/   00–05, the analysis          scripts/  helpers.py (Python), build_metadata.py (sample table)
R/helpers.R  R functions for 03, 04, 05   euler/    SLURM scripts for the cluster
metadata/    sample table and NCBI source tables   results/  small outputs kept in git
pixi.toml, pixi.lock   software environment
```
Not in git (regenerated): `raw_data/` (reads), `reference/` (mouse genome and annotation), `results/kallisto/`, FastQC reports, STAR alignments. Local only and unused: `processed_data/` (the authors' GEO count tables).

Paper: Hahn *et al.* 2023, *Cell* 186:4117 ([doi](https://doi.org/10.1016/j.cell.2023.07.027)); authors' code: github.com/OliInTheValley/SpatioTemporal_Analysis.
