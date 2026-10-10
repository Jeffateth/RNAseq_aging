# Sourced by the other scripts. Override REPO / WORK in the environment if needed.
REPO="${REPO:-$HOME/RNAseq_aging}"
WORK="${WORK:-$SCRATCH/rnaseq_aging}"          # scratch: 2.5 TB, files deleted after 15 days -> copy results back
REF="$WORK/ref"; SRA="$WORK/sra"; TMP="$WORK/tmp"; MAP="$WORK/mapping"; QC="$WORK/qc"; INDEX="$WORK/star_index"
GENOME="$REF/GRCm38.primary_assembly.genome.fa"; GTF="$REF/gencode.vM25.annotation.gtf"   # GENCODE vM25 = what kallisto used
mkdir -p "$REF" "$SRA" "$TMP" "$MAP" "$QC" "$INDEX" "$WORK/logs"
set +u; eval "$(conda shell.bash hook)"; conda activate rnaseq_star; set -u   # conda scripts are not set -u safe
