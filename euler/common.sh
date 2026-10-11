# Sourced by the other scripts. Override REPO / WORK in the environment if needed.
REPO="${REPO:-$HOME/RNAseq_aging}"
WORK="${WORK:-$SCRATCH/rnaseq_aging}"          # scratch: 2.5 TB, files deleted after 15 days -> copy results back
REF="$WORK/ref"; SRA="$WORK/sra"; TMP="$WORK/tmp"; QC="$WORK/qc"
MAP="$WORK/mapping"; INDEX="$WORK/star_index"; KAL="$WORK/kallisto"
GENOME="$REF/GRCm38.primary_assembly.genome.fa"; GTF="$REF/gencode.vM25.annotation.gtf"   # GENCODE vM25: same annotation for STAR and kallisto
TRANSCRIPTS="$REF/gencode.vM25.transcripts.fa.gz"; KAL_INDEX="$REF/gencode.vM25.kallisto.idx"
mkdir -p "$REF" "$SRA" "$TMP" "$QC/raw" "$MAP" "$INDEX" "$KAL" "$WORK/logs"
# load conda directly: a batch job does not inherit it from the login shell (SKIP_CONDA=1 only for local testing)
if [[ -z "${SKIP_CONDA:-}" ]]; then set +u; source "${CONDA_ROOT:-$HOME/miniconda3}/etc/profile.d/conda.sh"; conda activate rnaseq_star; set -u; fi
