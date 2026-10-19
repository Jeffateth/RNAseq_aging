# Glossary

## Study
| Term | Meaning |
|---|---|
| **AL / aDR** | *ad libitum* (free food, the control) / acute dietary restriction |
| **YMP / PBS** | young mouse plasma injections / saline control |
| **Rejuvenation** | a treatment that shifts aged expression back towards the young pattern |
| **Cohort / batch effect** | animals from different experiments differ for technical reasons, so the aging and rejuvenation cohorts are never pooled |
| **CAS** | the paper's "common aging score", 82 genes that change with age in most brain regions |

## Sequencing data
| Term | Meaning |
|---|---|
| **Bulk RNA-seq** | RNA of many cells pooled; one value per gene per sample |
| **Read / paired-end / mate** | one sequenced stretch (150 bases) / each fragment read from both ends / the two reads of a fragment (`<SRR>_1.fastq`, `<SRR>_2.fastq`; the same position in both is one fragment; at ~225 bases per fragment the mates overlap) |
| **Spot, read pair** | SRA's name for one sequenced cluster = one read pair here. NCBI's `Spots` = number of pairs; `AvgLength` 300 = both mates (150 + 150) |
| **Unstranded** | reads do not record which RNA strand they came from |
| **SRA / SRR, GEO, BioProject / BioSample** | NCBI's read archive / one run in it; Gene Expression Omnibus (GSE series, GSM sample); NCBI records of a study / a sample |
| **`.sra`, FASTQ** | NCBI's compressed raw format / plain-text reads with quality scores |
| **Phred quality (Q30)** | per-base confidence; Q30 = 1 error in 1,000 bases |

## Preprocessing
| Term | Meaning |
|---|---|
| **FastQC / MultiQC** | per-sample quality report / one report for all samples. PASS / WARN / FAIL come from FastQC's fixed thresholds (e.g. adapter content WARN > 5 %, FAIL > 10 % of reads); they are generic, so some WARNs (duplication, base content) are normal for RNA-seq |
| **Insert** | the part of a fragment that comes from the sample (cDNA), between the two adapters; a read longer than the insert runs on into the adapter |
| **Adapter** | synthetic sequence attached to fragments (here Nextera, `CTGTCTCTTATACACATCT…`); it must be cut off, and fastp cuts from the adapter to the read end |
| **Trimming (fastp)** | removing adapters, poly-G tails and low-quality reads |
| **Poly-G** | a run of G's at a read's end: the NovaSeq has no dye for G, so "no signal" past the end of a molecule is called G, with high quality. Mostly behind the adapter; `--trim_poly_g` removes the rest |
| **Duplicates** | identical read pairs; in bulk RNA-seq mostly genuine (highly expressed genes), so kept |

## Quantification
| Term | Meaning |
|---|---|
| **Reference genome, transcriptome, annotation** | the DNA sequence / the sequences of all known mRNAs / the GTF file saying where genes and exons are. We use GENCODE M39 (GRCm39), the latest on 14 October 2026 |
| **Transcript vs gene** | a gene can make several transcripts (isoforms); the link comes from the annotation, not from the reads |
| **kallisto** | builds an index of the transcript sequences, finds for each read pair the transcripts it is compatible with (**pseudoalignment**; the rate is the share of pairs compatible with at least one) and estimates reads per transcript, sharing reads that fit several isoforms |
| **STAR** | aligner placing each read at its exact genome position, across exon–exon junctions; ~30 GB RAM for the index |
| **GeneCounts** | STAR's per-gene counts of uniquely mapped reads (we use the unstranded column); multi-mapped reads (fitting several places) are dropped |
| **RSEM** | takes reads aligned to transcripts (here STAR's) and shares multi-mapping reads between isoforms by expectation-maximisation (EM); gives expected counts (not whole numbers), TPM, effective length |
| **tximport** | R package summing transcript estimates into gene counts for DESeq2 |
| **Counts, TPM** | reads per gene / transcripts per million (normalised for length and depth); DESeq2 needs counts |
| **Index** | pre-built searchable form of a reference |

## Statistics and pathways
| Term | Meaning |
|---|---|
| **DESeq2, DEG** | R package testing per gene whether counts differ between groups / a differentially expressed gene (padj < 0.05) |
| **log2FC** | log2 fold change: +1 = twice as high, −1 = half |
| **p / padj** | p-value / adjusted for testing ~15,000 genes at once (Benjamini–Hochberg) |
| **Wald statistic** | log2FC ÷ its standard error; the sign is the direction |
| **Dispersion** | variation of a gene's counts between mice of a group beyond counting noise; lower = more power |
| **Covariate** | a known factor (here sex) in the model so it does not blur the comparison |
| **vst, PCA** | variance-stabilising transformation for plotting / 2-D view of the main differences between samples |
| **Spearman ρ, Jaccard, Fisher** | rank correlation / shared ÷ combined size of two gene lists / test that two lists overlap more than by chance |
| **Outlier** | sample far from its group; SRR21355389 (potentially blood-contaminated, [pipeline.md](pipeline.md) step 8) is excluded from the aging analysis |
| **Gene set, GSEA (fgsea), NES** | genes sharing a function / test whether a set sits toward the top or bottom of genes ranked by the Wald statistic / its score (positive = up) |
| **MSigDB, Hallmark, GO:BP** | gene-set database / its 50 curated sets / Gene Ontology biological-process terms |
| **Reversal** | a treatment changes a gene or pathway opposite to how aging did |
| **Aging-signature score** | per sample: mean z-score of aging-up genes minus that of aging-down genes |

## Computing
**Euler / SLURM:** ETH's cluster / its scheduler (`sbatch`, `squeue`); jobs run on compute nodes, you log in on a login node. **Scratch:** temporary storage, files older than 15 days are deleted. **pixi / conda:** isolated tool environments (`pixi.toml`). **Kernel:** the program running a notebook's code (Python or R). **tmux:** keeps a terminal alive after you disconnect. **VPN:** ETH network access from outside.
