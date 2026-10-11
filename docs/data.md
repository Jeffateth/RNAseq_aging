# The data

## Study
Hahn *et al.*, "Atlas of the aging mouse brain reveals white matter as vulnerable foci", *Cell* 186:4117 (2023), [doi](https://doi.org/10.1016/j.cell.2023.07.027). Gene expression (RNA-seq) in 15 mouse brain regions across 7 ages and after two anti-aging ("rejuvenation") treatments. **We use one region, the cerebellum** (hindbrain, movement coordination): 27 samples from 27 mice.

| Dataset | Accession | Content |
|---|---|---|
| Aging | BioProject [PRJNA875066](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA875066) (GEO GSE212336) | untreated mice of increasing age |
| Rejuvenation | BioProject [PRJNA946653](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA946653) (GEO GSE227689) | 19-month-old mice given a treatment or its control |

## Design
All mice are **C57BL/6JN** (inbred, so near-identical genetically). A biopsy punch of cerebellum from a frozen brain slice is one sample; one mouse is sampled once, so the groups are different animals.

| Experiment | Mice | Groups | Treatment |
|---|---|---|---|
| **Aging** | 12 (2 male + 2 female per age) | 3, 18, 21 months, 4 mice each | none |
| **aDR** (acute dietary restriction) | 9 females, 19 months | AL ×4, aDR ×5 | from 18 months: food cut to 90 % for a week, then to 75 % of the controls' amount; ~4–5 weeks |
| **YMP** (young mouse plasma) | 6 males, 19 months | PBS ×3, YMP ×3 | from 18 months: injections every 3–4 days for ~4 weeks of plasma from 2-month-old mice (YMP) or saline (PBS) |

AL = *ad libitum* (free food, the control); PBS = phosphate-buffered saline (the injection control). Each treatment is compared only with its own control: the controls differ in sex (AL female, PBS male) and come from a different experiment than the aging mice.

## What one sample is
- **Bulk RNA-seq:** the RNA of all cells in the punch is sequenced together; each gene's value averages over neurons, glia and other cells.
- 1 mouse = 1 sample = 1 sequencing **run** (accession **SRR**).
- Poly-A mRNA library by **Smart-seq2** (cDNA from mRNA, fragmented with the Tn5 enzyme), sequenced on an **Illumina NovaSeq 6000**.
- **Paired-end, 2 × 150 bp:** each fragment is read from both ends (a read *pair*). All 27 runs are `PAIRED` in NCBI's table.
- **Unstranded:** reads do not record the RNA's strand (NCBI: "Unstranded mRNA-seq"; STAR's strand check gives 0.5).
- **Depth:** 6.8–30.2 million read pairs per run, about 54 GB on disk in NCBI's SRA format.

## Identifiers
| ID | Meaning | Example |
|---|---|---|
| BioProject (PRJNA…) | one study at NCBI | PRJNA875066 |
| GEO series (GSE…) / sample (GSM…) | the same study / one sample, in the Gene Expression Omnibus | GSE212336 / GSM6527152 |
| BioSample (SAMN…) | NCBI record of the biological sample | SAMN30602300 |
| **SRR** | one sequencing run: the unit we download | SRR21355392 |
| Mouse ID | the paper's animal number (also in sample titles, "…37ID") | 37 |

## The 27 samples
| Group | Run (SRR) | GEO sample | Mouse ID | Sex | Age (mo) | Read pairs (M) |
|---|---|---|---|---|---|---|
| Aging, 3 months | SRR21355401 | GSM6527141 | 29 | M | 3 | 15.6 |
| Aging, 3 months | SRR21355396 | GSM6527147 | 33 | M | 3 | 27.1 |
| Aging, 3 months | SRR21355392 | GSM6527152 | 37 | F | 3 | 17.0 |
| Aging, 3 months | SRR21355388 | GSM6527157 | 41 | F | 3 | 17.7 |
| Aging, 18 months | SRR21355398 | GSM6527145 | 31 | M | 18 | 6.8 |
| Aging, 18 months | SRR21355394 | GSM6527149 | 35 | M | 18 | 8.2 |
| Aging, 18 months | SRR21355390 | GSM6527154 | 39 | F | 18 | 29.5 |
| Aging, 18 months | SRR21355386 | GSM6527159 | 43 | F | 18 | 9.8 |
| Aging, 21 months | SRR21355397 | GSM6527146 | 32 | M | 21 | 14.9 |
| Aging, 21 months | SRR21355393 | GSM6527151 | 36 | M | 21 | 10.3 |
| Aging, 21 months | SRR21355389 | GSM6527155 | 40 | F | 21 | 17.7 |
| Aging, 21 months | SRR21355385 | GSM6527160 | 44 | F | 21 | 16.3 |
| Diet control (AL) | SRR23914577 | GSM7105321 | 63 | F | 19 | 29.1 |
| Diet control (AL) | SRR23914574 | GSM7105324 | 66 | F | 19 | 24.8 |
| Diet control (AL) | SRR23914572 | GSM7105326 | 68 | F | 19 | 15.0 |
| Diet control (AL) | SRR23914571 | GSM7105327 | 69 | F | 19 | 17.0 |
| Dietary restriction (aDR) | SRR23914579 | GSM7105319 | 61 | F | 19 | 16.0 |
| Dietary restriction (aDR) | SRR23914578 | GSM7105320 | 62 | F | 19 | 12.0 |
| Dietary restriction (aDR) | SRR23914576 | GSM7105322 | 64 | F | 19 | 14.5 |
| Dietary restriction (aDR) | SRR23914575 | GSM7105323 | 65 | F | 19 | 28.5 |
| Dietary restriction (aDR) | SRR23914573 | GSM7105325 | 67 | F | 19 | 13.9 |
| Plasma control (PBS) | SRR23914455 | GSM7105448 | 72 | M | 19 | 30.2 |
| Plasma control (PBS) | SRR23914452 | GSM7105451 | 80 | M | 19 | 16.3 |
| Plasma control (PBS) | SRR23914450 | GSM7105453 | 84 | M | 19 | 23.6 |
| Young plasma (YMP) | SRR23914454 | GSM7105449 | 73 | M | 19 | 24.5 |
| Young plasma (YMP) | SRR23914453 | GSM7105450 | 77 | M | 19 | 15.4 |
| Young plasma (YMP) | SRR23914451 | GSM7105452 | 81 | M | 19 | 19.3 |

Read pairs (M) = millions, before filtering. Sex: M = male, F = female.

## Not used
The other 14 brain regions; the paper's single-nucleus and spatial data; the authors' GEO count tables. The paper's full 82-gene "common aging score" list (supplementary Table S1) was not available, so we use the genes named in its figures.

## Metadata files (`metadata/`)
| File | Content |
|---|---|
| `runinfo_raw.csv`, `biosample_raw.csv` | NCBI tables as fetched (run statistics; sex, strain, treatment, mouse ID) |
| `srr_list.txt` | the 27 assigned run accessions |
| `samples.csv` | **one row per run**, built from the two tables above by `scripts/build_metadata.py`; the file all code reads (`sra_size_MB` is in MiB) |
