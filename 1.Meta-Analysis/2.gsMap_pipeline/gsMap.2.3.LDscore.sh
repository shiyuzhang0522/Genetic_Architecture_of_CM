# step.3 generate LD score

set -eo pipefail

# ---- 0) parse command-line arguments ----
if [[ $# -lt 2 ]]; then
  echo "Usage: sbatch $0 <sample_name> <chr>"
  exit 1
fi

sample_name="$1"   # first input from the command line
chr="$2"           # second input from the command line

echo "[INFO] Generating LDscore for sample ${sample_name}, chr${chr}..."

# ---- 1) define paths ----
BASE_WD="GWAS_ST/gsMap_pipeline"
working_dir="${BASE_WD}/wd_${sample_name}"
ldref_root="GWAS_ST/gsMap_pipeline/LD_ref.panel"

echo "[INFO] Sample:       ${sample_name}"
echo "[INFO] Chromosome:   ${chr}"
echo "[INFO] Working dir:  ${working_dir}"

# ---- 2) activate conda env (your preferred way) ----
if command -v conda &>/dev/null; then
  eval "$(conda shell.bash hook)"
  conda activate gsMap_new || true
else
  echo "[WARN] conda command not found; assuming gsmap is on PATH" >&2
fi

# Use all allocated CPUs if gsmap is multithreaded
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

# ---- 3) run generate ldscore step ----
gsmap run_generate_ldscore \
  --workdir "${working_dir}" \
  --sample_name "${sample_name}" \
  --chrom "${chr}" \
  --bfile_root "${ldref_root}/chr${chr}.UKBB_LD10k_MAF0.05" \
  --gtf_annotation_file "gsMap_resource/gencode.v49.basic.annotation.gtf" \
  --gene_window_size 50000

echo "[INFO] Finished step 3: run_generate_ldscore for ${sample_name} chr${chr}"