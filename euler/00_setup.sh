#!/bin/bash
# Once, on the Euler login node:  git clone https://github.com/Jeffateth/RNAseq_aging && bash RNAseq_aging/euler/00_setup.sh
# --override-channels: use only conda-forge + bioconda, so conda's Anaconda Terms-of-Service check on its default channels is not triggered.
# sra-tools=3.1 is pinned (Tutorial: avoids prefetch/fasterq-dump version conflicts); rsem is optional (isoform-level quantification).
set -euo pipefail
conda create -n rnaseq_star --override-channels -c conda-forge -c bioconda sra-tools=3.1 star fastp samtools rsem -y
