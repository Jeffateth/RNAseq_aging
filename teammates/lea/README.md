# Lea's analysis, refactored to the team setup

`Script.ipynb` is Lea's notebook (from [lhasch/System-Genomics-Project](https://github.com/lhasch/System-Genomics-Project), commit `8a133ad`), refactored to the team's common DESeq2 setup and cut down to five sections: import and QC, exploratory analysis, differential expression, Hallmark GSEA, and aging vs interventions. The original version is in git history (commit `9d815b9`).

- **Run:** from this folder with the repo's pixi environment (R kernel `ir`). Input: `../../for_teammates/` (RSEM tables, in git). Outputs: `output/`.
- `SUMMARY.md`: main findings, limitations and methodological decisions.
- `Script.R`: the same notebook as a plain R script (text cells as comments).
- Results match the team reference in `../../results/common/`: DEGs (175 / 118 / 50), Hallmark sets (15 / 16 / 24) and NES values.
