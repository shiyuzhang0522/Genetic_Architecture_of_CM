#!/usr/bin/env bash
#SBATCH --job-name=cojo_meta
#SBATCH --qos=high
#SBATCH --partition=GPU1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --array=1-22
#SBATCH --output=/gpfs/hpc/home/zhangsy/Meta-GWAS/GWAS_summary/Meta_analysis/post-GWAS/COJO/log/cojo_%A_%a.out
#SBATCH --error=/gpfs/hpc/home/zhangsy/Meta-GWAS/GWAS_summary/Meta_analysis/post-GWAS/COJO/log/cojo_%A_%a.err

############################################################
# COJO stepwise selection on meta-GWAS (GCTA; chr1–22)
#
# Input:  GWAS_meta_chr${chr}.txt (per-chromosome meta summary stats, formatted for COJO input)
# LD ref: chr${chr}_LD_panel_for_COJO (PLINK bfile; ~5k samples)
# Method: GCTA-COJO stepwise selection (--cojo-slct; --diff-freq 0.3)
# Output: COJO.meta.GWAS.chr${chr}.* (one set per chromosome)
#
# Author: Shelley
############################################################

set -euo pipefail

chr="${SLURM_ARRAY_TASK_ID}"

# -----------------------------
# Paths
# -----------------------------
bfile="/gpfs/hpc/home/zhangsy/0223_DATA_TRANSFER/WGS_data/COJO_LD_Panel/chr${chr}_LD_panel_for_COJO"
sumstats="/gpfs/hpc/home/zhangsy/Meta-GWAS/GWAS_summary/Meta_analysis/post-GWAS/COJO/GWAS_input/GWAS_meta_chr${chr}.txt"
outdir="/gpfs/hpc/home/zhangsy/Meta-GWAS/GWAS_summary/Meta_analysis/post-GWAS/COJO/results"
gcta_exec="/gpfs/hpc/home/zhangsy/software/gcta-1.94.1-linux-kernel-3-x86_64/gcta64"

mkdir -p "${outdir}/log" "${outdir}"

outpref="${outdir}/COJO.meta.GWAS.chr${chr}"

# -----------------------------
# Sanity checks
# -----------------------------
if [ ! -f "${sumstats}" ]; then
  echo "[ERROR] Missing sumstats: ${sumstats}" >&2
  exit 1
fi
if [ ! -f "${bfile}.bed" ] || [ ! -f "${bfile}.bim" ] || [ ! -f "${bfile}.fam" ]; then
  echo "[ERROR] Missing PLINK files for bfile prefix: ${bfile}" >&2
  exit 1
fi
if [ ! -x "${gcta_exec}" ]; then
  echo "[ERROR] GCTA executable not found or not executable: ${gcta_exec}" >&2
  exit 1
fi

# -----------------------------
# Run COJO (stepwise selection)
# -----------------------------
echo "[$(date)] Starting COJO stepwise selection on chr${chr}"
echo "  LD ref bfile: ${bfile}"
echo "  Sumstats:     ${sumstats}"
echo "  Threads:      ${SLURM_CPUS_PER_TASK}"
echo "  Output:       ${outpref}"

"${gcta_exec}" \
  --bfile "${bfile}" \
  --chr "${chr}" \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --cojo-file "${sumstats}" \
  --cojo-slct \
  --diff-freq 0.3 \
  --out "${outpref}"

echo "[$(date)] Finished COJO on chr${chr}"
