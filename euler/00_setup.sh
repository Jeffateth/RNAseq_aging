#!/bin/bash
# Once, on the Euler login node:  git clone https://github.com/Jeffateth/RNAseq_aging && bash RNAseq_aging/euler/00_setup.sh
set -euo pipefail
conda env create -f "$(dirname "$0")/environment.yml" --solver=libmamba
