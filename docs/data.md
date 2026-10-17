# The data

**Study.** Hahn *et al.*, "Atlas of the aging mouse brain reveals white matter as vulnerable foci", *Cell* 186:4117 (2023), [doi](https://doi.org/10.1016/j.cell.2023.07.027): RNA-seq of 15 mouse brain regions at 7 ages and after two anti-aging ("rejuvenation") treatments. **We use one region, the cerebellum**: 27 samples from 27 inbred C57BL/6JN mice (one punch of a frozen brain slice per mouse).

| Experiment | Accession | Mice | Groups | Treatment |
|---|---|---|---|---|
| **Aging** | [PRJNA875066](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA875066) (GEO GSE212336) | 12 (2 male + 2 female per age) | 3, 18, 21 months, 4 each | none |
| **aDR** (acute dietary restriction) | [PRJNA946653](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA946653) (GEO GSE227689) | 9 females, 19 months | AL ×4, aDR ×5 | from 18 months, food cut to 90 %, then 75 % of the controls' amount, for ~4–5 weeks |
| **YMP** (young plasma) | same | 6 males, 19 months | PBS ×3, YMP ×3 | from 18 months, plasma of 2-month-old mice or saline every 3–4 days for ~4 weeks |

Each treatment is compared only with its own control: the controls differ in sex (AL female, PBS male) and the rejuvenation mice come from a different experiment than the aging mice.

**One sample** = one mouse = one sequencing run (accession `SRR…`): bulk RNA-seq (cell types averaged), Smart-seq2 poly-A library with Tn5 fragmentation, Illumina NovaSeq 6000, **paired-end 2 × 150 bp**, **unstranded** (NCBI; STAR's strand check gives 0.5). 6.8–30.2 million read pairs per run, ~54 GB as `.sra`.

**Identifiers.** BioProject (PRJNA…) = a study; GEO series (GSE…) / sample (GSM…) = the same in the Gene Expression Omnibus; BioSample (SAMN…) = the biological sample; **SRR** = one sequencing run, the unit we download; mouse ID = the paper's animal number (in sample titles, "…37ID").

**Metadata (`metadata/`).** `samples.csv` has one row per run (SRR, GSM, BioSample, sex, age, mouse ID, group, depth) and is the file all code reads. It is built by `scripts/build_metadata.py` from `runinfo_raw.csv` and `biosample_raw.csv` (NCBI tables as fetched); `srr_list.txt` lists the 27 runs.

**Not used.** The other 14 regions, the paper's single-nucleus and spatial data, and the authors' GEO count tables. The paper's full 82-gene "common aging score" list (its Table S1) was not available, so we use the genes named in its figures.
