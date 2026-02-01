############################################################
# Post–meta-GWAS Quality Control Script
#
# Author: Shelley
# Date: 2025-05-10
#
# Purpose:
# Perform post–meta-analysis QC on melanoma GWAS summary statistics.
# Steps include:
#   1. Filter variants with large cross-cohort allele frequency discordance
#   2. Parse marker IDs into genomic coordinates
#   3. Remove malformed entries
#   4. Standardize allele representation and column structure
#   5. Export analysis-ready summary statistics
############################################################

# =========================
# Load required packages
# =========================
library(data.table)
library(dplyr)
library(tidyr)

# =========================
# Read meta-GWAS summary data
# =========================
meta_GWAS <- fread(
  "/gpfs/hpc/home/zhangsy/Meta-GWAS/GWAS_summary/Meta_analysis/Meta.GWAS.melanoma1.txt"
)

# =========================
# 1. Frequency discordance filter (ΔMAF)
# =========================
# ΔMAF = MaxFreq − MinFreq across cohorts
meta_GWAS[, difMAF := MaxFreq - MinFreq]

# Count variants exceeding thresholds
counts <- meta_GWAS[, .(
  n_gt_0.3 = sum(difMAF > 0.3, na.rm = TRUE),
  n_gt_0.5 = sum(difMAF > 0.5, na.rm = TRUE)
)]

cat("Variants with ΔMAF > 0.3:", counts$n_gt_0.3, "\n")
cat("Variants with ΔMAF > 0.5:", counts$n_gt_0.5, "\n")

# Remove variants with large frequency discordance
n_removed <- sum(meta_GWAS$difMAF > 0.3, na.rm = TRUE)
cat("Removing", n_removed, "variants with ΔMAF > 0.3\n")

meta_GWAS_filtered <- meta_GWAS[difMAF <= 0.3]

# =========================
# 2. Parse genomic coordinates
# =========================
# Split MarkerName (chr:pos:ref:alt)
meta_GWAS_filtered[, c("Chr", "Pos", "Ref", "Alt") :=
  tstrsplit(MarkerName, ":", fixed = TRUE)
]

meta_GWAS_filtered[, Pos := as.integer(Pos)]

# Remove entries with missing positions
meta_GWAS_filtered <- meta_GWAS_filtered[!is.na(Pos)]

# =========================
# 3. Count genome-wide significant variants
# =========================
setnames(meta_GWAS_filtered, "P-value", "P.value")
sig_count <- meta_GWAS_filtered[P.value < 5e-08, .N]
cat("Number of genome-wide significant variants (P < 5e-8):", sig_count, "\n")

# =========================
# 4. Select standardized output columns
# =========================
meta_GWAS_selected <- meta_GWAS_filtered[
  , .(
    Chr, Pos, Ref, Alt,
    MarkerName,
    Allele1, Allele2,
    Freq1, FreqSE,
    MinFreq, MaxFreq,
    Effect, StdErr,
    P.value
  )
]

# Replace ":" with "_" in marker IDs
meta_GWAS_selected[, MarkerName := gsub(":", "_", MarkerName)]

# =========================
# 5. Final formatting
# =========================
meta_GWAS_selected[, Chr := as.integer(Chr)]
setorder(meta_GWAS_selected, Chr, Pos)

meta_GWAS_selected[, `:=`(
  Allele1 = toupper(Allele1),
  Allele2 = toupper(Allele2)
)]

# =========================
# 6. Export QCed summary statistics
# =========================
fwrite(
  meta_GWAS_selected,
  file  = "/gpfs/hpc/home/zhangsy/Meta-GWAS/GWAS_summary/Meta_analysis/Meta.GWAS.melanoma1.QCed.txt",
  sep   = "\t",
  quote = FALSE
)
