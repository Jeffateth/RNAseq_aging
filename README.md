# RNAseq_aging – mouse cerebellum bulk RNA-seq (Hahn et al., Cell 2023)

Local, reproducible exploration of 27 SRA runs from two BioProjects. All 27 runs are **cerebellum**, paired-end, Illumina NovaSeq 6000, bulk mRNA-seq (BioSample/SRA metadata, verified per run).

| BioProject | GEO | Design | Runs | Raw size |
|---|---|---|---|---|
| PRJNA875066 | GSE212336 | Aging: 3 / 18 / 21 mo, 2M + 2F per age | 12 | 19.8 GB |
| PRJNA946653 | GSE227689 | 19 mo rejuvenation: aDR (5F) vs AL (4F); young plasma YMP (3M) vs PBS (3M) | 15 | 30.6 GB |

BioProject membership was read from NCBI runinfo, not inferred from SRR prefixes (SRR21355* → PRJNA875066, SRR239144*/SRR239145* → PRJNA946653).
No metadata is missing for the fields below. Age is not a structured field in BioSample; it is parsed from the sample title (e.g. "Cerebellum Male 18mo 31ID").

## Layout
```
notebooks/     01_acquire_data (prefetch)  02_quantify_kallisto (FastQC, fastp, kallisto)  03_differential_expression (R, DESeq2)
scripts/       helpers.py  (idempotent Python pipeline functions)      R/helpers.R  (R functions)
metadata/      srr_list.txt, runinfo_raw.csv, biosample_raw.csv, sra_metadata*.csv, GSE212336_series_matrix.txt.gz
processed_data/ GEO count tables (unmodified, see SOURCES.txt) + cerebellum_counts_27runs.csv.gz
raw_data/      .sra (kept), temporary FASTQ (deleted after quantification)   [git-ignored]
reference/     GENCODE vM25 transcriptome + annotation, kallisto index       [git-ignored]
results/       qc/ (FastQC, MultiQC, fastp), kallisto/ (git-ignored), de/ (DESeq2 tables per contrast), figures
explore_rnaseq.ipynb   first look at the GEO count tables
pixi.toml      environment (kallisto, fastp, R/Bioconductor, Jupyter)
```

## Storage estimate (raw data – NOT downloaded)
- Combined `.sra` download: **50.4 GB** (491 M read pairs). Estimated gz FASTQ: ~65 GB; `fasterq-dump` temp space ~7× the `.sra` of one run (max ≈ 22 GB, one run at a time).
- Free disk: 128 GB. Downloading all raw data (~50 GB + ~65 GB FASTQ) would nearly fill it unless `.sra` files are deleted after conversion.
- 16 GB RAM is fine for download/conversion; STAR mapping of a mouse genome needs ~30 GB for index generation, which is too much locally (use Euler, or kallisto).
- **Raw reads are not needed for initial exploration**: GEO provides count tables (25 MB total).

## Processed data
- `GSE212336_BulkSeq_Counttable.txt.gz` (21,076 genes × 847 samples, 15 regions) 
- `GSE227689_BulkSeqRejuvenation_Counttable.txt.gz` (21,076 × 238)
- Columns are `BulkSeq_<region>_<mouseID>`; cerebellum = `cer`. Mouse IDs match `mouse_id` in `sra_metadata.csv`; all 27 matched.
- Raw integer counts, no missing values, no duplicated genes.

## Observations from the notebook
- PC1 (31%) separates the two projects entirely: strong cohort/batch effect. Do not compare across projects.
- Aging cohort: PC2 separates sex; one 21-mo female (SRR21355389, mouse 40) is a strong PC1 outlier; inspect before DE.
- Group sizes are small (3–5 per group), so power is limited.

## Reproduce
```bash
pixi install                 # environment (kallisto, fastp, DESeq2, ...); needs Homebrew sra-tools for prefetch
pixi run jupyter lab         # then run notebooks/01, 02, 03 in order (03 uses the R kernel; start Jupyter/VS Code kernels from the pixi env)
```
Raw reads are processed locally with kallisto; only STAR (genome index ≈ 30 GB RAM) would need Euler.
