# Pipeline and decisions

Each step with the choice made and why. † = tested or checked on the data. Terms: [glossary.md](glossary.md); data: [data.md](data.md).

```
27 runs (SRA) → QC → trim → count (RSEM; kallisto as check) → DESeq2 per cohort → compare methods → aging vs treatment
 Euler (01)       01/02    01/02      01                               03                 04                05
```

| # | Step | What we did | Why |
|---|---|---|---|
| 1 | **Data** | 27 cerebellum runs from two BioProjects, `prefetch` from NCBI SRA (`euler/01_download.sh`) | follows the course tutorial |
| 2 | **QC** | FastQC on **all reads** of every run, before and after trimming (`qc` and `qctrim` stages), MultiQC. Q30 in 83–87 % of bases, Nextera **adapter in 17–34 % of reads**, duplication 3–8 % † | all reads cost only minutes; a 2 M-pair subset gave the same verdicts † |
| 3 | **Trim** | fastp, **same settings for all runs**: adapters auto-detected, `--trim_poly_g`, minimum length 20, default quality filter (drops 4–6 % of pairs); duplicates kept | Same settings keep groups technically comparable †. Duplicates are mostly genuine and cannot be told from PCR copies without UMIs. Adapter cutting alone leaves poly-G tails in 0.8–2.2 % of reads, and a few G-rich genes (Gm19220, Gm26822, Grip2) had inflated kallisto counts †. Minimum length 20 vs the tutorial's 25 differs by 0.004 % of pairs †. Poly-A is not trimmed: poly-A tails are rare (0.2–0.4 % of reads) and do not match the transcripts †. No quality end-trimming (the paper's Trim Galore did): quality is already high. After trimming, FastQC on all 27 runs: Nextera adapter 0.01–0.05 % of reads (PASS in all 54 files), poly-G tails gone, quality equal or slightly better † |
| 4 | **Reference** | GENCODE **M39** (GRCm39): primary-assembly genome, comprehensive annotation (78,289 genes), transcripts | newest release (latest on 14 Oct 2026); one annotation for all methods so gene ids match. The paper used GRCm38, so coordinates and some genes differ. An earlier vM25 run was replaced |
| 5 | **Count** | Three methods on identical reads, run on Euler (notebook 01): **kallisto** (`--threads 8`), **STAR** (unique reads, unstranded column), **RSEM** (EM on STAR's transcriptome alignments) | RSEM counts on STAR's alignments, so STAR and RSEM differ only in counting; kallisto differs in mapping too. Share of pairs counted: kallisto 65–72.5 %, STAR 47–59 %, RSEM 47–60 %. The rules differ (e.g. multi-mappers), so the rates are not directly comparable. Threads are fixed because they change kallisto slightly †. Strand check 0.5 = unstranded † |
| 6 | **Gene counts** | `tximport` sums kallisto / RSEM to genes; STAR counts are used as they are | DESeq2 needs gene counts |
| 7 | **Differential expression** | DESeq2 (Wald). Cohorts **separately**, each treatment against its **own control**. Aging `~ sex + age group` (18 vs 3, 21 vs 3 months); treatments `~ group`. Genes need ≥ 10 counts in ≥ (smallest group) samples; significant = padj < 0.05, no fold-change cutoff | PC1 and the sample clustering split the cohorts completely, and the controls differ in sex †. With 3–5 mice a fold-change cutoff would remove most calls. The aging DEGs are also clustered by their course over 3 → 18 → 21 months (notebook 03, 2c) |
| 8 | **Outlier** | **Exclude** SRR21355389 from the aging analysis (switch in notebook 03) | potentially blood-contaminated: red-blood-cell and neutrophil genes are 50–150× the median of the other aging samples in TPM (Hba-a1 4,499 vs 57; Hbb-bs 4,242 vs 49; Alas2 49 vs 3.5; S100a8 3,959 vs 25; S100a9 1,756 vs 15), Col1a1 119 vs 1.6 (meninges/vessels). This is an inference from expression only; some other aging samples also have more haemoglobin, so it is a gradient and the exclusion is a judgement call. Keeping it roughly doubles the aging DEGs (RSEM / kallisto: 239 / 279 vs 114 / 132); the lists overlap at Jaccard 0.43 / 0.42 † |
| 9 | **Compare methods** | kallisto vs RSEM through the same DESeq2 code (notebook 04, `results/compare/`); technical vs biological variation on TPM: kallisto vs RSEM on the same mouse against two mice of the same group and sex | Both estimate transcript abundance with EM and give counts and TPM, so they compare on the same level; kallisto is carried through step 10 as well. Once a fixed per-gene offset between the methods is removed (it cancels within a method), the method difference is smaller than the mouse difference for every mouse (SD 0.30–0.46 vs 0.47–0.70 log2) † |
| 10 | **Aging vs treatment** | signature = 21 vs 3 month DEGs; (i) DEG overlap and direction (Fisher); (ii) per-mouse signature score, z-scored **within each treated/control pair**, Welch t-test; (iii) fgsea on the Wald statistic (Hallmark, GO:BP); (iv) the paper's named genes | the paper's full 82-gene list was not available. A rank test cannot reach p < 0.1 with 3 vs 3 mice, hence Welch; p-values are indicative |

### Choosing the counting method
No method is assumed best. Notebook 04 (section 4) measures the criteria below for kallisto and RSEM (`results/compare/method_evidence.csv`). No single one decides.

| Criterion | Measures | Better | Caveat |
|---|---|---|---|
| Counted read pairs | share of trimmed pairs counted in a gene | higher | data use, not correctness |
| Genes detected | genes with ≥ 10 counts in ≥ 3 runs | not decisive | more genes, more low-count noise |
| Replicate agreement | median pairwise correlation of log2 CPM within a group | higher | precision, not bias |
| Dispersion | median DESeq2 dispersion | lower | count scales differ |
| Stability to the outlier | DEG overlap with and without SRR21355389 | higher | one mouse |
| Positive controls | known aging / aDR genes go up | all up | sanity check only, circular as a selector |

Cost: kallisto index 11 min, 4–13 min per run, peak 3–23.7 GB; STAR index ~30 GB and 22 min, 5–13 min per run; RSEM 13–97 min per run, peak 6–14.8 GB. References (check before citing): Bray *et al.* 2016 *Nat Biotechnol* 34:525 (kallisto); Li & Dewey 2011 *BMC Bioinformatics* 12:323 (RSEM); Soneson, Love & Robinson 2015 *F1000Research* 4:1521 (tximport).

**Decision (18 October 2026): STAR alignments with RSEM.** RSEM combines two inherent advantages: STAR's base-to-base, splice-aware alignment of every read to the genome, and RSEM's EM estimate of how many reads come from each isoform (which also corrects for gene length). kallisto skips the alignment (pseudoalignment), and STAR `GeneCounts` counts a read pair only if it is uniquely mapped and overlaps the exons of exactly one gene, so multi-mapping and ambiguous pairs are discarded. Time and computation are not constraints here, so the more thorough route is affordable, and no benchmark with a known truth is needed to choose it. For transparency, kallisto is kept as a sensitivity check (notebooks 04, 05): both agree in direction for every shared DEG. STAR `GeneCounts` are still produced (notebook 01) but no longer used in the analysis.

Notebook 03 runs on RSEM (`METHOD`) and notebook 05 shows RSEM (`SHOW`); change either to repeat the analysis with another method.

### fastp vs cutadapt † (one run, SRR21355398, early exploratory test on vM25; not scripted in the repo)
| Trimming | Pairs left | Pseudoaligned |
|---|---|---|
| **fastp** (ours) | 6.39 M | **65.2 %** |
| cutadapt, adapter only (keeps all reads) | 6.77 M | 62.8 % |
| cutadapt, Trim Galore-like (paper's tool; also cuts a 1-base adapter match and quality < 30 ends) | 6.58 M | 63.9 % |

Gene counts agree (Spearman 0.999), but 1,500–1,800 of 14,800 genes move by > 5 % (median 1–2 %), and *Dock10* went from 262 to 18,920 counts under the Trim Galore-like settings (not investigated). One run cannot show which tool is better; trimming at all raised the share by 2–3 points (five runs).
