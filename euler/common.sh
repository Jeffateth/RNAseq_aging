# Sourced by the other scripts. Override REPO / WORK in the environment if needed.
REPO="${REPO:-$HOME/RNAseq_aging}"
WORK="${WORK:-$SCRATCH/rnaseq_aging}"          # scratch: 2.5 TB, files deleted after 15 days -> copy results back
# Reference: GENCODE M39 on GRCm39, one annotation for kallisto, STAR and RSEM. Outputs that depend on it carry the suffix _M39,
# so results of the earlier GENCODE vM25 (GRCm38) run on the same scratch space are never mixed with them. The .sra files and the
# FastQC reports of the raw reads do not depend on the reference and are shared.
SRA="$WORK/sra"; TMP="$WORK/tmp"; QC="$WORK/qc"
# POLYG=off: the sensitivity run without poly-G trimming (everything else identical); its counting outputs get the suffix _nopolyg.
POLYG="${POLYG:-on}"; [[ "$POLYG" == off ]] && SFX=_nopolyg || SFX=""
REF="$WORK/ref_M39"; INDEX="$WORK/star_index_M39"; MAP="$WORK/mapping_M39$SFX"; KAL="$WORK/kallisto_M39$SFX"; RSEM="$WORK/rsem_M39$SFX"
GENOME="$REF/GRCm39.primary_assembly.genome.fa"; GTF="$REF/gencode.vM39.annotation.gtf"
TRANSCRIPTS="$REF/gencode.vM39.transcripts.fa.gz"; KAL_INDEX="$REF/gencode.vM39.kallisto.idx"
RSEM_REF="$REF/rsem_M39/rsem_M39"                                   # prefix of the RSEM reference files
mkdir -p "$REF/rsem_M39" "$SRA" "$TMP" "$QC/raw" "$MAP" "$INDEX" "$KAL" "$RSEM" "$WORK/logs"
# load conda directly: a batch job does not inherit it from the login shell (SKIP_CONDA=1 only for local testing)
if [[ -z "${SKIP_CONDA:-}" ]]; then set +u; source "${CONDA_ROOT:-$HOME/miniconda3}/etc/profile.d/conda.sh"; conda activate rnaseq_star; set -u; fi
