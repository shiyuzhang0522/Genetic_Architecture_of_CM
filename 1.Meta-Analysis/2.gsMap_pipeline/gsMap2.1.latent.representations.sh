#### step.1 find latent representations

#!/bin/bash
# ---- 0) get sample name from command line ----
if [ $# -lt 1 ]; then
  echo "Usage: $0 <sample_name>"
  echo "Example: $0 MEL131HAP_CA21_4_2"
  exit 1
fi

sample_name="$1"

# ---- 1) define paths ----
BASE_WD="/public/home/hpc8301200407/WGS/GWAS_ST/gsMap_pipeline"
working_dir="${BASE_WD}/wd_${sample_name}"
h5ad_file="${BASE_WD}/ST_h5ad/${sample_name}.h5ad"

# make working directory
mkdir -p "${working_dir}"

# sanity check: h5ad exists
if [ ! -f "${h5ad_file}" ]; then
  echo "[ERROR] h5ad file not found: ${h5ad_file}" >&2
  exit 1
fi

echo "[INFO] Sample:       ${sample_name}"
echo "[INFO] Working dir:  ${working_dir}"
echo "[INFO] Input h5ad:   ${h5ad_file}"

# ---- 2) run gsmap latent representation step ----
gsmap run_find_latent_representations \
  --workdir "${working_dir}" \
  --sample_name "${sample_name}" \
  --input_hdf5_path "${h5ad_file}" \
  --data_layer 'count' \
  --epochs 1000

echo "[INFO] Finished run_find_latent_representations for ${sample_name}"

######### submit the scripts for step 1 ############
#!/usr/bin/env bash
set -euo pipefail

# ---- paths ----
PIPE_BASE="/public/home/hpc8301200407/WGS/GWAS_ST/gsMap_pipeline"
STEP1_SCRIPT="${PIPE_BASE}/Step1.latent.representations.sh"
LOG_DIR="${PIPE_BASE}/logs_step1"
mkdir -p "${LOG_DIR}"

# ---- list of samples ----
SAMPLES=(
  MEL109ZSD_CA20_10_9
  MEL126_HTP_CA21_11_12
  MEL131HAP_CA21_4_2
  MEL162YJG_CA21_4_30
  MEL170YFY_CA21_5_21
  MEL176_PDYCA2021_6_4
  MEL201ZLJ_CA21_8_26
  MEL222LBLCA22_9_9
  MEL239CXE_CA22_5_13_2
  MEL256XXD_CA22_7_8
  MEL271TPH_CA22_9_9
  MEL74TXM_CA20_6_19
  MEL75CZL_CA20_7_3
  MEL99LXQ_CA20_9_11_3_MEL222LBL_CA211224
)

# ---- (optional) activate env ----
# source "$HOME/miniconda3/etc/profile.d/conda.sh"
# conda activate gsMap

# ---- loop over samples ----
for s in "${SAMPLES[@]}"; do
  echo "[`date`] >>> Step1 for ${s}"
  out_log="${LOG_DIR}/${s}.step1.out"
  err_log="${LOG_DIR}/${s}.step1.err"

  bash "${STEP1_SCRIPT}" "${s}" >"${out_log}" 2>"${err_log}"

  echo "[`date`] >>> Finished ${s} (logs: ${out_log}, ${err_log})"
  echo
done