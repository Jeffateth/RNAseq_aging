# Pipeline and decisions

Each step with the choice made and why. † = tested or checked on the data. Terms: [glossary.md](glossary.md); data: [data.md](data.md).

```
27 runs (SRA) → QC → trim → count (kallisto, STAR, RSEM) → DESeq2 per cohort → compare methods → aging vs treatment
 Euler (01)       01/02    01/02      01/02                       03                 04                05
```

| # | Step | What we did | Why |
|---|---|---|---|
| 1 | **Data** | 27 cerebellum runs from two BioProjects, `prefetch` from NCBI SRA (`euler/01_download.sh`) | follows the course tutorial |
| 2 | **QC** | FastQC on **all reads** of every run, before and after trimming (`qc` and `qctrim` stages), MultiQC. Q30 in 83–87 % of bases, Nextera **adapter in 17–34 % of reads**, duplication 3–8 % † | all reads cost only minutes; a 2 M-pair subset gave the same verdicts † |
| 3 | **Trim** | fastp, **same settings for all runs**: adapters auto-detected, `--trim_poly_g`, minimum length 20, default quality filter (drops 4–6 % of pairs); duplicates kept | same for all runs, so no technical difference between groups †. Duplicates are mostly genuine and cannot be told from PCR copies without UMIs. Adapter cutting alone leaves the G tail in 0.8–2.2 % of poly-G reads †, and a few G-rich genes (Gm19220, Gm26822, Grip2) had inflated kallisto counts †. No Q30 end-trimming (the paper's Trim Galore): quality is already high |
| 4 | **Reference** | GENCODE **M39** (GRCm39): primary-assembly genome, comprehensive annotation (78,289 genes), transcripts | newest release (latest on 14 Oct 2026); one annotation for all methods so gene ids match. The paper used GRCm38, so coordinates and some genes differ. An earlier vM25 run was replaced |
| 5 | **Count** | Three methods on identical reads, run on Euler (notebook 01): **kallisto** (`--threads 8`), **STAR** (unique reads, unstranded column), **RSEM** (EM on STAR's transcriptome alignments) | RSEM is a counting step on STAR's alignments: STAR vs RSEM differ only in counting, kallisto also in mapping. kallisto pseudoaligns 65–72.5 % of pairs; STAR counts 47–59 % in genes, RSEM 47–60 %. The rates follow different rules (e.g. multi-mappers counted or not), so they are not directly comparable; reads outside annotated transcripts are counted by no method. Threads are fixed because they change kallisto slightly †; strand check 0.5 = unstranded † |
| 6 | **Gene counts** | `tximport` sums kallisto / RSEM to genes; STAR counts are used as they are | DESeq2 needs gene counts |
| 7 | **Differential expression** | DESeq2 (Wald). Cohorts **separately**, each treatment against its **own control**. Aging `~ sex + age group` (18 vs 3, 21 vs 3 months); treatments `~ group`. Genes need ≥ 10 counts in ≥ (smallest group) samples; significant = padj < 0.05, no fold-change cutoff | PC1 splits the cohorts completely and the controls differ in sex †. With 3–5 mice a fold-change cutoff would remove most calls |
| 8 | **Outlier** | **Exclude** SRR21355389 from the aging analysis (switch in notebook 03) | top genes are haemoglobins, neutrophil genes, collagens: blood/meningeal contamination. Keeping it roughly doubles the aging DEGs (kallisto / STAR / RSEM: 279 / 307 / 239 vs 132 / 168 / 114); the lists overlap at Jaccard 0.42 / 0.48 / 0.43 † |
| 9 | **Compare methods** | every pair of methods through the same DESeq2 code (notebook 04, `results/compare/`) | all three are carried through step 10 |
| 10 | **Aging vs treatment** | signature = 21 vs 3 month DEGs; (i) DEG overlap and direction (Fisher); (ii) per-mouse signature score, z-scored **within each treated/control pair**, Welch t-test; (iii) fgsea on the Wald statistic (Hallmark, GO:BP); (iv) the paper's named genes | the paper's full 82-gene list was not available. A rank test cannot reach p < 0.1 with 3 vs 3 mice, hence Welch; p-values are indicative |

### Choosing the counting method
No method is assumed best. Notebook 04 (section 5) measures the criteria below, saved to `results/compare/method_evidence.csv`. No single one decides; a method is a poor choice only if clearly worse on several.

| Criterion | Measures | Better | Caveat |
|---|---|---|---|
| Counted read pairs | share of trimmed pairs counted in a gene | higher | data use, not correctness |
| Genes detected | genes with ≥ 10 counts in ≥ 3 runs | not decisive | more genes, more low-count noise |
| Replicate agreement | median pairwise correlation of log2 CPM within a group | higher | precision, not bias |
| Dispersion | median DESeq2 dispersion | lower | count scales differ |
| Stability to the outlier | DEG overlap with and without SRR21355389 | higher | one mouse |
| Agreement with the others | mean Jaccard of DEG calls | higher | not accuracy; STAR and RSEM share the alignment |
| Positive controls | known aging / aDR genes go up | all up | sanity check only, circular as a selector |

Practical: kallisto needs only the transcript sequences (index 11 min, 4–13 min per run, peak 3–23.7 GB); STAR needs ~30 GB for its index (22 min) and takes 5–13 min per run; RSEM needs STAR's alignments plus a reference (13–97 min per run, peak 6–14.8 GB). References (check before citing): Bray *et al.* 2016 *Nat Biotechnol* 34:525 (kallisto); Li & Dewey 2011 *BMC Bioinformatics* 12:323 (RSEM); Soneson, Love & Robinson 2015 *F1000Research* 4:1521 (tximport).

**Decision:** pending; record the chosen method and its evidence here. Until then notebook 03 uses `METHOD`, notebook 05 draws figures for `SHOW`, and its tables cover all three.

### How kallisto gets to gene counts
kallisto works on **transcripts**; which transcript belongs to which gene comes from the annotation (M39: 481,956 transcripts of 78,289 genes; each transcript name carries its gene id). It (1) indexes the transcript sequences once, (2) pseudoaligns each read pair to the transcripts it is compatible with without placing it on the genome, (3) estimates reads per transcript, sharing reads that fit several isoforms by EM. `tximport` then sums a gene's transcripts (e.g. isoforms of 800 and 200 reads give a gene of 1,000).

### fastp vs cutadapt † (one run, SRR21355398, early exploratory test on vM25; not scripted in the repo)
| Trimming | Pairs left | Pseudoaligned |
|---|---|---|
| **fastp** (ours) | 6.39 M | **65.2 %** |
| cutadapt, adapter only (keeps all reads) | 6.77 M | 62.8 % |
| cutadapt, Trim Galore-like (paper's tool; also cuts a 1-base adapter match and quality < 30 ends) | 6.58 M | 63.9 % |

Gene counts agree (Spearman 0.999), but 1,500–1,800 of 14,800 genes move by > 5 % (median 1–2 %), and *Dock10* went from 262 to 18,920 counts under the Trim Galore-like settings (not investigated). One run cannot show which tool is better; trimming at all raised the share by 2–3 points (five runs).
