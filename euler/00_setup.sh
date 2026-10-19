#!/bin/bash
# Once, on the Euler login node:  git clone https://github.com/Jeffateth/RNAseq_aging && bash RNAseq_aging/euler/00_setup.sh
# Safe to re-run: creates the environment `rnaseq_star`, or adds/updates these packages in it.
# --override-channels: only conda-forge + bioconda, so conda's Terms-of-Service check on its default channels is not triggered.
# Versions are pinned to what the results were produced with (sra-tools 3.1 is also what the tutorial needs).
set -euo pipefail
PKGS="sra-tools=3.1.1 star=2.7.11b fastp=1.4.0 samtools=1.24 rsem=1.3.3 kallisto=0.52.0 fastqc=0.13.0"
if conda env list | grep -qE '^rnaseq_star[[:space:]]'; then
  conda install -n rnaseq_star --override-channels -c conda-forge -c bioconda -y $PKGS
else
  conda create -n rnaseq_star --override-channels -c conda-forge -c bioconda -y $PKGS
fi
