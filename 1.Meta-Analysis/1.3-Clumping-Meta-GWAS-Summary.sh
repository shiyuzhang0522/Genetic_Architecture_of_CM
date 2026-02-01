#!/bin/bash
#SBATCH --job-name=clump_meta
#SBATCH --qos=high
#SBATCH --partition=GPU2             
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --array=1-22                 # chromosomes 1–22
#SBATCH --output=/gpfs/hpc/home/zhangsy/Meta-GWAS/GWAS_summary/Meta_analysis/post-GWAS/Clump/log/clump_%A_%a.out
#SBATCH --error=/gpfs/hpc/home/zhangsy/Meta-GWAS/GWAS_summary/Meta_analysis/post-GWAS/Clump/log/clump_%A_%a.err

############################################################
# LD clumping of meta-GWAS summary statistics (PLINK2; chr1–22)
#
# Input:  Meta.GWAS.melanoma1.QCed.txt
# Method: PLINK2 --clump (p1=5e-8, p2=1e-5, r2=0.1, kb=1000)
# Output: clumped.meta.GWAS.chr${chr}.clumped (one per chromosome)
#
# Author: Shelley
############################################################

set -euo pipefail

chr=${SLURM_ARRAY_TASK_ID}

# 1) pick the right reference panel bfile
if [ "$chr" -le 6 ]; then
  bfile="/gpfs/hpc/home/zhangsy/0309_datatrans/chr${chr}/plink_QCed_chr${chr}"
else
  bfile="/gpfs/hpc/home/zhangsy/0223_DATA_TRANSFER/WGS_data/chr${chr}/plink_QCed_chr${chr}"
fi

# 2) meta‐analysis summary statistics
sumstats="/gpfs/hpc/home/zhangsy/Meta-GWAS/GWAS_summary/Meta_analysis/Meta.GWAS.melanoma1.QCed.txt"

# 3) output directory & prefix
outdir="/gpfs/hpc/home/zhangsy/Meta-GWAS/GWAS_summary/Meta_analysis/post-GWAS/Clump"
outpref="${outdir}/clumped.meta.GWAS.chr${chr}"

# 4) PLINK2 executable
plink2_exec="/gpfs/hpc/home/zhangsy/software/plink2"

echo "[$(date)] Starting clumping for chr${chr}"
echo "  Reference BED: $bfile"
echo "  Sumstats:      $sumstats"
echo "  Threads:       ${SLURM_CPUS_PER_TASK}"

"$plink2_exec" \
  --bfile "$bfile" \
  --threads ${SLURM_CPUS_PER_TASK} \
  --clump "$sumstats" \
  --clump-p1 5e-8 --clump-p2 1e-5 --clump-r2 0.1 --clump-kb 1000 --clump-snp-field MarkerName --clump-field P.value \
  --out "$outpref"

echo "[$(date)] Finished clumping chr${chr}; results in ${outpref}.clumped"