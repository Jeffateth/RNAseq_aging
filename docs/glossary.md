# Glossary

Terms are grouped by where they appear in the workflow. Abbreviations are spelled out at first use in each notebook and listed here.

## Study and design
| Term | Meaning |
|---|---|
| **Cerebellum** | hindbrain region; the only region analysed here |
| **Aging cohort** | the 12 untreated mice aged 3, 18 and 21 months |
| **AL** | *ad libitum*: free access to food (diet control group) |
| **aDR** | acute dietary restriction: food reduced to 75% of AL for several weeks |
| **YMP / PBS** | young mouse plasma injections / saline injections (control) |
| **Rejuvenation** | here: a treatment that shifts aged gene expression back towards the young pattern |
| **Cohort / batch effect** | animals processed in different experiments differ for technical reasons; hence the aging and rejuvenation cohorts are never pooled |
| **C57BL/6JN** | inbred mouse strain (substrain "N" from the NIA colony) |
| **CAS** | the paper's "common aging score": a signature of 82 genes that change with age in most brain regions |

## Sequencing data
| Term | Meaning |
|---|---|
| **RNA-seq** | measuring gene expression by sequencing RNA converted to cDNA |
| **Bulk RNA-seq** | RNA of many cells pooled; one value per gene per sample |
| **Read** | one sequenced stretch of cDNA (here 150 bases) |
| **Paired-end** | each fragment is read from both ends: two reads (the *mates*) per fragment |
| **Mate** | one of the two reads of a fragment: mate 1 (`<SRR>_1.fastq`) is read from one end, mate 2 (`<SRR>_2.fastq`) from the other end on the opposite strand; the same position in both files is the same fragment. Fragments here are ~225 bases and reads 150, so the mates overlap |
| **Read pair** | the two mates of one fragment; "depth" counts pairs |
| **Spot** | SRA's term for one sequenced cluster on the flow cell; for paired-end data one spot = one read pair. `Spots` in the NCBI run table = number of read pairs |
| **Spots_with_mates** | spots that have both mates; equal to `Spots` when every fragment was read from both ends (true for all 27 runs) |
| **AvgLength (300)** | SRA reports the length of a whole *spot*: mate 1 + mate 2 = 150 + 150 bases, not the length of one read |
| **Bases** | `Spots` × `AvgLength` (total sequenced bases) |
| **Unstranded** | no record of which RNA strand a read came from |
| **Smart-seq2, Tn5** | the library preparation method; Tn5 is the enzyme that cuts and tags cDNA |
| **Illumina NovaSeq** | the sequencing instrument |
| **SRA / SRR** | NCBI's Sequence Read Archive / one run in it (accession `SRR…`) |
| **GEO (GSE / GSM)** | Gene Expression Omnibus: series / sample identifiers |
| **BioProject / BioSample** | NCBI records for a study / a biological sample |
| **`.sra`, FASTQ** | NCBI's compressed raw format / the plain-text format of reads plus quality scores (`.fastq`, `.fastq.gz`) |
| **Phred quality (Q30)** | per-base confidence; Q30 = 1 error in 1,000 bases |

## Preprocessing
| Term | Meaning |
|---|---|
| **QC** | quality control |
| **FastQC / MultiQC** | per-sample quality report / one report combining all samples |
| **Insert** | the part of a sequenced fragment that comes from the sample (a piece of a transcript, as cDNA), between the two adapters. Reads longer than the insert run on into the adapter |
| **Adapter** | synthetic sequence attached to fragments; read into when the fragment is shorter than the read length, so it must be removed (here: Nextera transposase sequence, `CTGTCTCTTATACACATCT…`); fastp cuts from the adapter to the read end |
| **Trimming (fastp)** | removing adapters and low-quality ends from reads; fastp also filters low-quality reads |
| **Poly-G** | run of G's at a read's end: the NovaSeq has no dye for G, so "no signal" (past the end of a molecule) is called G, with high quality. Mostly behind the adapter (cut with it); fastp's `--trim_poly_g` removes the rest |
| **Duplicates** | identical read pairs; in bulk RNA-seq mostly genuine (highly expressed genes), so they are kept |
| **PASS / WARN / FAIL** | FastQC's own verdict per check (written to each report's `summary.txt`), from fixed thresholds in FastQC's `limits.txt` (e.g. adapter content: WARN > 5 % of reads, FAIL > 10 %; thresholds in notebook 02). They are generic, not RNA-seq specific, so some WARNs (duplication, base content) are normal for RNA-seq |
| **Sequence Duplication, Adapter Content, …** | the FastQC checks used to decide on trimming |

## Quantification
| Term | Meaning |
|---|---|
| **Reference genome / transcriptome / annotation** | the mouse DNA sequence / the sequences of all known mRNAs / the file (GTF) saying where genes and exons are. We use GENCODE release vM25 (mouse genome version GRCm38 = mm10) |
| **Transcript vs gene** | one gene can make several transcripts (isoforms). The transcript → gene link comes from the annotation, not from the sequencing data; kallisto estimates per transcript, `tximport` sums to genes |
| **kallisto** | estimates how much of each known *transcript* is present: builds an index of the transcript sequences, finds for each read pair the transcripts it is compatible with (pseudoalignment), then estimates reads per transcript, sharing out reads that fit several isoforms; fast, small memory |
| **Pseudoalignment (rate)** | finding which transcripts a read is compatible with, without placing it on the genome; the rate is the share of read pairs compatible with at least one transcript |
| **STAR** | aligner: places each read at its exact position in the genome, including across exon–exon junctions; needs ~30 GB RAM for the index |
| **Index** | pre-built searchable form of the reference |
| **GeneCounts** | STAR's per-gene read counts (unique reads only); we use the unstranded column |
| **tximport** | R package that sums kallisto's transcript estimates into gene counts for DESeq2, using the annotation's transcript → gene link |
| **Counts, TPM** | reads per gene / transcripts per million (counts normalised for gene length and depth); DESeq2 needs counts |
| **Multi-mapped read** | a read that fits several places; STAR's GeneCounts drops these |

## Statistics
| Term | Meaning |
|---|---|
| **DESeq2** | R package testing, per gene, whether counts differ between groups |
| **DEG** | differentially expressed gene (here: adjusted p < 0.05) |
| **log2FC** | log2 fold change: +1 = twice as high, −1 = half |
| **p / padj** | p-value / p-value adjusted for testing ~15,000 genes at once (Benjamini–Hochberg); padj < 0.05 means ≤ 5% of such calls are expected to be false |
| **Wald statistic** | DESeq2's per-gene test statistic (log2FC ÷ its standard error); sign shows direction |
| **Covariate** | a known factor (here sex) included in the model so it does not blur the comparison |
| **vst** | variance-stabilising transformation: puts counts on a scale where PCA and distances behave well |
| **PCA** | principal component analysis: a 2-D view of the main differences between samples |
| **Spearman correlation (ρ)** | rank-based correlation; robust to outliers |
| **Fisher's exact test** | tests whether two gene lists overlap more than by chance |
| **Jaccard index** | shared ÷ combined size of two gene lists |
| **Outlier** | sample far from its group; SRR21355389 (blood/meninges contamination) is excluded from the aging analysis |

## Pathway analysis
| Term | Meaning |
|---|---|
| **Gene set / pathway** | a list of genes sharing a function |
| **MSigDB, Hallmark, GO:BP** | database of gene sets / its 50 curated "hallmark" sets / Gene Ontology biological-process terms |
| **GSEA (fgsea)** | gene set enrichment analysis: tests whether a set's genes sit toward the top or bottom of the genes ranked by the Wald statistic |
| **NES** | normalised enrichment score: positive = set up-regulated, negative = down-regulated |
| **Reversal** | a treatment changes a gene or pathway opposite to how aging changed it |
| **Aging-signature score** | per sample: mean z-score of aging-up genes minus mean z-score of aging-down genes |

## Computing
| Term | Meaning |
|---|---|
| **Euler / SLURM** | ETH's computing cluster / its job scheduler (`sbatch`, `squeue`) |
| **Login node, compute node** | where you log in and prepare jobs / where jobs run |
| **Scratch** | temporary cluster storage; files older than 15 days are deleted |
| **conda / pixi** | package managers that install the tools in an isolated environment (`pixi.toml`) |
| **Jupyter notebook, kernel** | document mixing text, code and output / the program that runs its code (Python or R here) |
| **tmux** | keeps a terminal session alive after you disconnect |
| **VPN** | ETH network access from outside the campus |
