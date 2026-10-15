# The heavy steps on the ETH cluster (Euler)

Scripts for full-read QC, kallisto, STAR and RSEM (notebook 01). Run in this order on Euler; see the notebook for the full commands and for why each step is needed.

| Script | Runs on | Purpose |
|---|---|---|
| `00_setup.sh` | login node | create or update the conda environment `rnaseq_star` (sra-tools, STAR, fastp, samtools, RSEM, kallisto, FastQC; versions pinned) |
| `01_download.sh` | login node, in `tmux` | download the mouse genome, annotation and transcripts (GENCODE M39, GRCm39) and the 27 `.sra` runs; resumable |
| `02_star_index.sh` | SLURM job | build the STAR genome index (~30 GB RAM) |
| `02b_kallisto_index.sh` | SLURM job | build the kallisto index (a few GB RAM) |
| `02c_rsem_reference.sh` | SLURM job | build the RSEM reference from the same genome and annotation |
| `03_process_array.sh` | SLURM array, 1 task per run | per run: convert to FASTQ, then (by `STAGES`) FastQC on the raw reads, fastp trimming, kallisto, STAR with gene counting, RSEM on STAR's transcriptome alignments; stages with existing output are skipped; delete FASTQ |
| `04_collect.sh` | login node | pack results into `star_results`, `kallisto_results` (~0.7 GB), `rsem_results` (~0.4 GB) and `qc_results` tarballs for copying to the laptop |
| `submit.sh` | login node | submits the three indices/references, then a light array (FastQC + kallisto, little memory) and a STAR array and an RSEM array (task i waits for STAR task i) with the right dependencies |
| `common.sh` | sourced by the others | paths on scratch storage, loads conda and the environment |
| `tunnel.sh` | laptop | optional: keeps one SSH connection to Euler open so later `ssh`/`scp` need no new login |

Scripts begin with `#SBATCH` lines: SLURM reads them as job requests (cores, memory, time limit, log name). They are comments to bash.
Euler needs the ETH network or VPN. Results go to `$SCRATCH/rnaseq_aging/` (files there are deleted after 15 days).
