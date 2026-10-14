# Pipeline and decisions

From raw reads to the aging-vs-treatment result, with the choice made at each step and why. † = tested against an alternative or checked on the data; unmarked = reasoned choice. Terms: [glossary.md](glossary.md). Data: [data.md](data.md).

```
27 runs (SRA) → QC → trim → count (kallisto, STAR) → DESeq2 per cohort → compare methods → aging vs treatment
 notebook 00      01/02    01/02      01/02                 03                 04                05
```

## Steps
| # | Step | What we did | Why, and alternatives |
|---|---|---|---|
| 1 | **Data** | 27 cerebellum runs from two BioProjects, `prefetch` from NCBI SRA, checksum-validated; `.sra` kept (54 GB), FASTQ temporary | Follows the course tutorial; compact. ENA FASTQ would be ~25 % larger |
| 2 | **QC** | FastQC on **all reads** of every run, combined with MultiQC. Result: base quality good (83–87 % of bases ≥ Q30), Nextera **adapter in 17–34 % of reads**, duplication 3–8 % † | All reads, not a subset: quality differs by sample and it costs minutes per run. A 2 M-pair subset gave the same verdicts † |
| 3 | **Trim** | fastp, **identical settings for all 27 runs**: adapters found automatically from the read pair, poly-G trimming (`--trim_poly_g`), minimum length 20, default quality filter (drops 4–6 % of pairs). Duplicates are kept | *All runs alike:* adapters fail in every sample, and trimming only some would add a technical difference between groups †. *Duplicates:* in bulk RNA-seq they are mostly genuine, and without UMIs PCR copies cannot be told apart. *Poly-G:* most G runs sit behind the adapter and go with it, but 1.5–3 % of poly-G reads keep theirs after adapter cutting alone †; `--trim_poly_g` removes these too. Without it a few G-rich genes (Gm19220, Gm26822, Grip2) had inflated kallisto counts (none differentially expressed); with it they drop, e.g. Gm19220 2576 → 94 †. *No Q30 end-trimming* as in the paper's Trim Galore: quality is already high |
| 4 | **Reference** | GENCODE **M39** (genome **GRCm39**): primary-assembly genome, comprehensive annotation (78,289 genes), transcripts | The current mouse genome and annotation. One annotation for all three methods so gene ids match. The paper used the older genome (GRCm38) and does not state its annotation release, so coordinates and some genes differ from the paper's. An earlier run of this pipeline used GENCODE vM25 (GRCm38); its results were replaced |
| 5 | **Count** | Three methods on identical trimmed reads: **kallisto** (transcript estimates, `--threads 8`), **STAR** (per-gene counts, unique reads, unstranded column) and **RSEM** (expected counts from STAR's transcriptome alignments, EM shares multi-mapping reads). **FastQC, fastp, kallisto, STAR and RSEM are submitted to Euler**, not run on the laptop (notebook 01; notebook 02 only runs the first steps locally if `RUN_LOCALLY = True`) | kallisto fits on a laptop, STAR needs ~30 GB RAM and is the paper's method; RSEM is the tutorial's way to handle reads that fit several isoforms. It is a counting step on STAR's alignments, not an independent pipeline: STAR vs RSEM differs only in counting, kallisto vs either also in mapping. Strand check 0.5 = unstranded †. Thread count is fixed because it changes kallisto's estimates slightly †. kallisto assigns only 62–69 % of pairs: the rest map to introns or unannotated sequence (STAR: 78–89 % map to the genome, 51–61 % fall in genes) and are simply not counted |
| 6 | **Gene counts** | `tximport` sums kallisto transcripts to genes | Needed by DESeq2 |
| 7 | **Differential expression** | DESeq2 (Wald test). Cohorts **separately**, each treatment against its **own control**. Aging: `~ sex + age group` (18 vs 3, 21 vs 3 months); treatments: `~ group`. Genes need ≥ 10 counts in ≥ (smallest group size) samples; significant = padj < 0.05, no fold-change cutoff | *Separate cohorts:* PC1 splits them completely, and controls differ in sex (AL female, PBS male) †. DESeq2 is the paper's and tutorial's choice; edgeR/limma were not compared. Only three ages exist, so pairwise contrasts, not a slope. A fold-change cutoff would remove most calls with 3–5 mice |
| 8 | **Outlier** | **Exclude** SRR21355389 from the aging analysis (one switch in notebook 03) | Top genes are haemoglobins, neutrophil genes, collagens: blood/meningeal contamination. Keeping it gives 269 instead of 130 aging DEGs; 118 are shared (log2FC rank correlation 0.98) † |
| 9 | **Compare methods** | kallisto, STAR and RSEM, every pair, through the same DESeq2 code | Agreement of DEG calls and log2FC per pair (notebook 04, `results/compare/`). All three methods are carried through step 10; which one is primary is decided by the evidence in section "Choosing the counting method" |
| 10 | **Aging vs treatment** | Aging signature = 21 vs 3 month DEGs. (i) overlap of treatment DEGs and direction (Fisher test); (ii) per-mouse signature score, z-scored **within each treated/control pair**, Welch t-test; (iii) fgsea on the Wald statistic with Hallmark and GO biological-process sets; (iv) the paper's named genes | The paper's full 82-gene list was not available. Pair-wise z-scoring keeps cohort differences out. A rank test cannot reach p < 0.1 with 3 vs 3 mice, hence Welch; read p-values as indicative |

### Choosing the counting method
We do not assume a best method. Three are run on identical trimmed reads (step 5) and compared with measurements fixed in advance (notebook 04, section 6, saved to `results/compare/method_evidence.csv`). No single measurement decides; a method is a poor choice only if it is clearly worse on several.

| Criterion | What it measures | Better | Caveat |
|---|---|---|---|
| Counted read pairs (%) | share of trimmed pairs that end up in a gene count | higher | says how much data is used, not whether it is assigned correctly; counting reads that fit several places raises it |
| Genes detected | genes with ≥ 10 counts in ≥ 3 runs | not decisive | more genes also means more low-count noise |
| Replicate agreement | median pairwise correlation of log2 CPM between the mice of a group | higher | measures precision, not bias |
| Dispersion | median DESeq2 gene-wise dispersion (unexplained variation within groups) | lower | count scales differ between methods |
| Stability to the outlier | overlap of the aging DEGs with and without SRR21355389 | higher | one mouse only |
| Agreement with the other methods | mean Jaccard of DEG calls over all pairs | higher | agreement is not accuracy; STAR and RSEM share the alignment, so they are not independent |
| Positive controls | known aging and aDR genes go the expected way | all up | sanity check only; using it to choose would be circular |

**Practical points** (not computed in the notebook): kallisto needs only the transcript sequences, a few GB RAM and minutes for the index, no genome alignment. STAR needs ~30 GB RAM for the genome index (22 min on Euler) and 3.5–13 min per run (earlier run). RSEM needs STAR's transcriptome alignments plus its own reference; its run time is not measured yet. kallisto and RSEM share reads between isoforms by EM; STAR GeneCounts counts uniquely mapped reads only. Pointers (check before citing): Bray *et al.* 2016, *Nat Biotechnol* 34:525 (kallisto); Li & Dewey 2011, *BMC Bioinformatics* 12:323 (RSEM, also linked in the tutorial); Soneson, Love & Robinson 2015, *F1000Research* 4:1521 (transcript-level estimates summarised with tximport for gene-level analysis).

**Decision:** pending. Record the chosen method and the evidence for it here after notebook 04 has run on the M39 results. Until then notebook 03 runs on the method set in `METHOD`, notebook 05 draws figures for the method set in `SHOW`, and its tables cover all three.

### How kallisto gets from reads to gene counts
kallisto works on **transcripts** (the mRNA isoforms a gene can make), not on genes. Which transcripts belong to which gene is **known from the annotation** (GENCODE M39: 481,956 transcripts of 78,289 genes; each transcript name in our reference carries its gene id), and kallisto does not discover it. Its job is to estimate how much of each transcript is present:
1. **Index:** builds a searchable table of the known transcript sequences (once).
2. **Pseudoalignment:** for each read pair, finds the transcripts it is *compatible* with, without placing it exactly on the genome.
3. **Quantification:** estimates reads per transcript. A read compatible with several transcripts (isoforms sharing an exon) is shared out between them statistically (expectation–maximization).
4. **Gene level:** `tximport` (notebook 03) sums each gene's transcript estimates; DESeq2 works on these gene counts. STAR, by contrast, counts reads per gene directly.

Illustration (made-up values): gene A with isoforms A1 = 800 and A2 = 200 estimated reads gives gene A = 1,000 reads.

### fastp vs cutadapt † (one run, SRR21355398)
Question: does the choice of trimming tool change the counts? We trimmed the same run three ways and ran kallisto on each result.

| Trimming | What it does | Read pairs left for kallisto | Share kallisto could assign to a transcript |
|---|---|---|---|
| **fastp** (ours) | cuts adapters, drops low-quality reads | 6.39 M | **65.2 %** |
| cutadapt, adapter only | cuts adapters, keeps all reads | 6.77 M | 62.8 % |
| cutadapt, Trim Galore-like (the paper's tool) | cuts adapters, also when only the first adapter base is visible at the read end (so it sometimes cuts one genuine base); cuts low-quality read ends (quality < 30, error chance > 1 in 1,000) | 6.58 M | 63.9 % |

- **Gene counts are almost the same:** the ranking of genes agrees between the three (Spearman 0.999).
- **Individual genes still move:** 1,500–1,800 of 14,800 genes change by more than 5 % (typically 1–2 %). One gene, *Dock10*, went from 262 to 18,920 counts under the Trim Galore-like settings; we did not investigate why.
- **Where the tools differ:** mainly in what they discard. fastp drops about 5–6 % of the read pairs for low quality; adapter-only cutadapt drops none.
- **Takeaway:** the choice matters a little, but one run cannot show which tool is better. Trimming at all (any tool) raised the assigned share by 2–3 points compared with untrimmed reads (five runs).
