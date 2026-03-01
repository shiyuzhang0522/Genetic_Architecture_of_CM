#!/bin/bash
# Script to Calculate ALT allele counts for cases and controls in UKBB WGS data (UKB-RAP via swiss-army-knife)
# Usage:
#   bash get.counts.sh <chr> <case|control>
# Example:
#   bash get.counts.sh 20 case

set -euo pipefail

######################### 0. Inputs ##########################################
chr=${1:-}
group=${2:-}   # "case" or "control"

if [[ -z "${chr}" || -z "${group}" ]]; then
  echo "[ERROR] Usage: bash $0 <chr> <case|control>"
  exit 1
fi

if ! [[ "${chr}" =~ ^([1-9]|1[0-9]|2[0-2])$ ]]; then
  echo "[ERROR] chr must be 1-22. Got: ${chr}"
  exit 1
fi

if [[ "${group}" != "case" && "${group}" != "control" ]]; then
  echo "[ERROR] group must be 'case' or 'control'. Got: ${group}"
  exit 1
fi

echo "[INFO] Running for chromosome ${chr}; sample group = ${group}"

######################### 1. Input files #####################################
# Use a normal path string here; quoting is handled when building the cmd string
pfile_prefix="/mnt/project/Bulk/\"DRAGEN WGS\"/\"DRAGEN population level WGS variants, PLINK format [500k release]\"/ukb24308_c${chr}_b0_v1"

variant_list="/mnt/project/BHR/counts/variant.list/chr${chr}.variant_list.keep.txt"

if [[ "${group}" == "case" ]]; then
  keep_file="/mnt/project/BHR/counts/FID_IID/cases.keep"
else
  keep_file="/mnt/project/BHR/counts/FID_IID/controls.keep"
fi

######################### 2. Output location #################################
# RAP project output folder (DNAnexus)
output_folder="project-GxVyGzQJ5V1qVf14k938ZFf8:/BHR/counts/allele_counts"

# Output prefix for PLINK2 results (.acount/.log)
output_prefix="chr${chr}.${group}.ALT_counts"

######################### 3. PLINK2 command ##################################
# IMPORTANT: quote paths with spaces inside the command string
cmd="plink2 \
  --no-pheno \
  --pfile ${pfile_prefix} \
  --keep ${keep_file} \
  --extract ${variant_list} \
  --freq counts \
  --out ${output_prefix} "

echo "[INFO] Command to run inside swiss-army-knife:"
echo "${cmd}"

######################### 4. Submit job to RAP ###############################
dx run app-swiss-army-knife \
  -icmd="${cmd}" \
  --folder="${output_folder}" \
  --instance-type=mem2_ssd1_v2_x8 \
  --name="ALTcount_chr${chr}_${group}" \
  --priority=high \
  -y




#!/bin/bash
# submit_all_counts.sh
# Submit ALT allele-count jobs for chr 1-22, for both case and control
# Usage: bash submit_all_counts.sh

SCRIPT="./get.counts.sh"   # <- adjust if your script is elsewhere

for chr in $(seq 1 22); do
  for group in case control; do
    echo "[SUBMIT] chr${chr} ${group}"
    bash "${SCRIPT}" "${chr}" "${group}"
  done
done

echo "[DONE] Submitted chr1-22 for case & control"