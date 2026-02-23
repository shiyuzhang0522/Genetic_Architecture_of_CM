# Annotate WGS variant lists using VEP with LoFTEE and dbNSFP plugins under Singularity
#!/usr/bin/env bash
#SBATCH --job-name=vep_anno
#SBATCH --partition=cpuQ
#SBATCH --qos=cpuq
#SBATCH --account=pi_dengguangtong
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=190464
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# ==============================================================================
# VEP annotation (GRCh38) with LoFTEE + dbNSFP under Singularity
#
# Usage:
#   sbatch vep_annotate.sh 1
#   sbatch vep_annotate.sh 22
#
# Notes:
#   - This script expects an input site-only VCF per chromosome
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# 0) Environment
# ------------------------------------------------------------------------------
module purge
module load singularity

# ------------------------------------------------------------------------------
# 1) CLI
# ------------------------------------------------------------------------------
if [[ $# -ne 1 ]]; then
  echo "ERROR: Expected 1 argument: <chromosome_number>" >&2
  echo "Usage: $0 <chromosome_number>" >&2
  exit 1
fi

CHR="$1"
if ! [[ "$CHR" =~ ^([1-9]|1[0-9]|2[0-2])$ ]]; then
  echo "ERROR: chromosome_number must be an integer in [1..22]. Got: $CHR" >&2
  exit 1
fi

# ------------------------------------------------------------------------------
# 2) CONFIG (edit here)
# ------------------------------------------------------------------------------
SCRATCH_DIR="/public/home/hpc8301200407/Shelley"

INPUT_DIR="${SCRATCH_DIR}/plink2"
OUTPUT_DIR="${SCRATCH_DIR}/variant.annotation/VEP_new"
LOG_DIR="${OUTPUT_DIR}/log"

CACHE_DIR="/public/home/hpc8301200407/tool/ensembl-vep/cache.annotations"
SIF_IMAGE="/public/home/hpc8301200407/tool/vep.sif"

# Reference FASTA (toplevel) from Ensembl
REF_FASTA_DIR="${OUTPUT_DIR}/FASTA"
REF_FASTA="${REF_FASTA_DIR}/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
SAMTOOLS_BIN="${HOME}/samtools/bin/samtools"

# ------------------------------------------------------------------------------
# 3) Derived paths + checks
# ------------------------------------------------------------------------------
INPUT_VCF="${INPUT_DIR}/chr${CHR}.copied.pvar.vcf"
OUTPUT_VCF="${OUTPUT_DIR}/chr${CHR}.VEP.anno.vcf"

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}" "logs"

[[ -f "${INPUT_VCF}" ]] || { echo "ERROR: Input VCF not found: ${INPUT_VCF}" >&2; exit 1; }
[[ -f "${SIF_IMAGE}" ]] || { echo "ERROR: Singularity image not found: ${SIF_IMAGE}" >&2; exit 1; }
[[ -d "${CACHE_DIR}" ]] || { echo "ERROR: VEP cache dir not found: ${CACHE_DIR}" >&2; exit 1; }
[[ -f "${REF_FASTA}" ]] || { echo "ERROR: Reference FASTA not found: ${REF_FASTA}" >&2; exit 1; }

# FASTA indexes expected by VEP/samtools (common pitfalls)
[[ -f "${REF_FASTA}.fai" ]] || { echo "ERROR: Missing FASTA index: ${REF_FASTA}.fai" >&2; exit 1; }

# This gzi naming matches your original script; keep as-is, but check existence.
GZI_PATH="${REF_FASTA%.gz}.gz.gzi"
[[ -f "${GZI_PATH}" ]] || { echo "ERROR: Missing FASTA gzi index: ${GZI_PATH}" >&2; exit 1; }

# ------------------------------------------------------------------------------
# 4) VEP command (runs inside container)
# ------------------------------------------------------------------------------
LOFTEE_DIR="${CACHE_DIR}/loftee"

VEP_CMD=(
  vep
  -i "${INPUT_VCF}"
  --assembly GRCh38
  --vcf
  --format vcf
  --cache
  --dir_cache /cache
  --fasta /ref.fa.gz
  -o "${OUTPUT_VCF}"
  --plugin "LoF,loftee_path:/plugins,\
human_ancestor_fa:/plugins/vep_data/human_ancestor.fa.gz,\
conservation_file:/plugins/vep_data/loftee.sql,\
gerp_bigwig:/plugins/vep_data/gerp_conservation_scores.homo_sapiens.GRCh38.bw"
  --plugin "dbNSFP,/plugins/dbNSFP5.1a.txt.gz,REVEL_score,CADD_phred"
  --everything
  --offline
  --force_overwrite
)

# ------------------------------------------------------------------------------
# 5) Singularity exec (binds)
# ------------------------------------------------------------------------------
SING_CMD=(
  singularity exec --cleanenv
  -W "${HOME}"
  --bind "${SCRATCH_DIR}"
  --bind "${CACHE_DIR}":/cache
  --bind "${LOFTEE_DIR}":/plugins
  --bind "${REF_FASTA}":/ref.fa.gz
  --bind "${REF_FASTA}.fai":/ref.fa.gz.fai
  --bind "${GZI_PATH}":/ref.fa.gz.gzi
)

# samtools bind only if it exists (keeps script robust for others)
if [[ -x "${SAMTOOLS_BIN}" ]]; then
  SING_CMD+=( --bind "${SAMTOOLS_BIN}":/usr/local/bin/samtools )
fi

SING_CMD+=( "${SIF_IMAGE}" )

# ------------------------------------------------------------------------------
# 6) Run
# ------------------------------------------------------------------------------
echo "[INFO] CHR=${CHR}"
echo "[INFO] INPUT_VCF=${INPUT_VCF}"
echo "[INFO] OUTPUT_VCF=${OUTPUT_VCF}"
echo "[INFO] CACHE_DIR=${CACHE_DIR}"
echo "[INFO] REF_FASTA=${REF_FASTA}"
echo

echo "[INFO] Running command:"
printf "  %q " "${SING_CMD[@]}" "${VEP_CMD[@]}"
echo
echo

"${SING_CMD[@]}" "${VEP_CMD[@]}"

echo "[INFO] Done: ${OUTPUT_VCF}"