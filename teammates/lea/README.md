# Lea's analysis, refactored to the team setup (proposed common baseline)

`Script.ipynb` is Lea's notebook (from [lhasch/System-Genomics-Project](https://github.com/lhasch/System-Genomics-Project), commit `8a133ad`), reduced to the team's agreed DESeq2 setup in five sections: import and QC, exploratory analysis, differential expression, Hallmark GSEA, aging vs interventions. It is offered for review, not as a final analysis.

| File | Content |
|---|---|
| `HANDOVER.md` | start here: agreed setup, results, figures, limitations, open decisions, presentation check |
| `CHANGES.md` | what was retained, changed or removed from Lea's original, and why |
| `Script.ipynb` / `Script.R` | the notebook / the same as a plain R script |
| `output/` | DE tables, Hallmark and GO BP results, comparison table, five figures |

Run from this folder with the repo's pixi environment: `pixi run jupyter nbconvert --to notebook --execute --inplace Script.ipynb`. The original notebook is in git history (commits `9e240d7`, `9d815b9`).
