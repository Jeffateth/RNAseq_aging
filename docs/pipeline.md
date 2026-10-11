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
| 4 | **Reference** | GENCODE vM25 (genome GRCm38/mm10): genome, annotation, transcripts | Last release on the paper-era genome version; one annotation for both methods so gene ids match. The paper does not state its version, so this is a choice |
| 5 | **Count** | **kallisto** (transcript estimates, `--threads 8`) and **STAR** (per-gene counts, unique reads, unstranded column), on identical trimmed reads. **FastQC, fastp, kallisto and STAR are submitted to Euler**, not run on the laptop (notebook 01; notebook 02 only runs them locally if `RUN_LOCALLY = True`) | kallisto fits on a laptop, STAR needs ~30 GB RAM and is the paper's method. Strand check 0.5 = unstranded †. Thread count is fixed because it changes kallisto's estimates slightly †. kallisto assigns only 62–69 % of pairs: the rest map to introns or unannotated sequence (STAR: 78–89 % map to the genome, 51–61 % fall in genes) and are simply not counted |
| 6 | **Gene counts** | `tximport` sums kallisto transcripts to genes | Needed by DESeq2 |
| 7 | **Differential expression** | DESeq2 (Wald test). Cohorts **separately**, each treatment against its **own control**. Aging: `~ sex + age group` (18 vs 3, 21 vs 3 months); treatments: `~ group`. Genes need ≥ 10 counts in ≥ (smallest group size) samples; significant = padj < 0.05, no fold-change cutoff | *Separate cohorts:* PC1 splits them completely, and controls differ in sex (AL female, PBS male) †. DESeq2 is the paper's and tutorial's choice; edgeR/limma were not compared. Only three ages exist, so pairwise contrasts, not a slope. A fold-change cutoff would remove most calls with 3–5 mice |
| 8 | **Outlier** | **Exclude** SRR21355389 from the aging analysis (one switch in notebook 03) | Top genes are haemoglobins, neutrophil genes, collagens: blood/meningeal contamination. Keeping it gives 269 instead of 130 aging DEGs; 118 are shared (log2FC rank correlation 0.98) † |
| 9 | **Compare methods** | kallisto vs STAR through the same DESeq2 code | Same direction for every shared DEG; DEG overlap (Jaccard) 0.37–0.58; log2FC rank correlation 0.86–0.96 among significant genes †. Both methods are carried through step 10 |
| 10 | **Aging vs treatment** | Aging signature = 21 vs 3 month DEGs. (i) overlap of treatment DEGs and direction (Fisher test); (ii) per-mouse signature score, z-scored **within each treated/control pair**, Welch t-test; (iii) fgsea on the Wald statistic with Hallmark and GO biological-process sets; (iv) the paper's named genes | The paper's full 82-gene list was not available. Pair-wise z-scoring keeps cohort differences out. A rank test cannot reach p < 0.1 with 3 vs 3 mice, hence Welch; read p-values as indicative |

### How kallisto gets from reads to gene counts
kallisto works on **transcripts** (the mRNA isoforms a gene can make), not on genes. Which transcripts belong to which gene is **known from the annotation** (GENCODE vM25: 142,604 transcripts, 55,401 genes; each transcript name in our reference carries its gene id), and kallisto does not discover it. Its job is to estimate how much of each transcript is present:
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
