# step.4 spatial S-LDSC

set -eo pipefail

## Step 4: running spatial LDSC ##

# --- parse arguments ---
if [[ $# -lt 1 ]]; then
  echo "[ERROR] Usage: $0 <sample_name>" >&2
  exit 1
fi

sample_name="$1"

# --- activate conda env ---
if command -v conda &>/dev/null; then
  eval "$(conda shell.bash hook)"
  conda activate gsMap_new || true
else
  echo "[WARN] conda command not found; assuming gsmap is on PATH" >&2
fi

# number of processes to use
NUM_PROCESSES="${SLURM_CPUS_PER_TASK:-16}"

# --- paths ---
BASE_WD="GWAS_ST/gsMap_pipeline"
LOG_DIR="${BASE_WD}/logs_step4"
working_dir="${BASE_WD}/wd_${sample_name}"
sumstats_file="${BASE_WD}/melanoma_gwas_for_gsmap.tsv"
w_prefix="${working_dir}/${sample_name}/generate_ldscore/w_ld/weights."

# mkdir -p "${LOG_DIR}"
# mkdir -p "${working_dir}"

echo "[INFO] Sample:          ${sample_name}"
echo "[INFO] Working dir:     ${working_dir}"
echo "[INFO] Sumstats file:   ${sumstats_file}"
echo "[INFO] Weight prefix:   ${w_prefix}"
echo "[INFO] Num processes:   ${NUM_PROCESSES}"

# --- run spatial LDSC ---
gsmap run_spatial_ldsc \
  --workdir "${working_dir}" \
  --sample_name "${sample_name}" \
  --trait_name "Melanoma" \
  --sumstats_file "${sumstats_file}" \
  --w_file "${w_prefix}" \
  --num_processes "${NUM_PROCESSES}"

echo "[INFO] Finished Step 4: spatial LDSC for ${sample_name}"