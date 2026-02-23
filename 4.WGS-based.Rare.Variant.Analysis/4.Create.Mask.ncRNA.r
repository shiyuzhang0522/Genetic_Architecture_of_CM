#!/usr/bin/env Rscript
# Create masks for ncRNA 
############################################################
## R script: Build 3-score variant masks and generate
##           SAIGE-style group/annotation tables for ncRNA
##           exon-only aggregates (subsetted)
##
## Inputs (per chromosome):
##   1) QC-passed variant list (VCF; first 5 columns)
##   2) CADD (SNV + INDEL), GERP, JARVIS score tables
##   3) ncRNA exon aggregates
##
## Outputs:
##   1) var_wide.<...>.txt   (space-delimited, no header)
##   2) mask_wide.<...>.txt  (space-delimited, no header)
##
## Author: Shiyu Zhang (Shelley)
## Date: 2025-06-24
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# ------------------------------------------------------------------------------
# 0) CLI + config
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript build_ncRNA_exon_masks.R <chr> <subset_id>", call. = FALSE)
}

chr <- args[1]
subset_id <- args[2]

if (!grepl("^([1-9]|1[0-9]|2[0-2])$", chr)) {
  stop("chr must be an integer in [1..22].", call. = FALSE)
}
if (!grepl("^[0-9]+$", subset_id)) {
  stop("subset_id must be an integer-like string (e.g., 1, 2, 3...).", call. = FALSE)
}

message("[INFO] Chromosome: ", chr)
message("[INFO] ncRNA exon aggregate subset_id: ", subset_id)

# ---- Base directory (HIDE absolute paths) ----
# Users should set WGS_MELANOMA_BASE once in their shell, e.g.:
#   export WGS_MELANOMA_BASE=/path/to/WGS_melanoma
BASE_DIR <- Sys.getenv("WGS_MELANOMA_BASE", unset = "PATH_TO_PROJECT_ROOT")

# ---- Input directories (relative to BASE_DIR) ----
CADD_DIR   <- file.path(BASE_DIR, "WGS_annotation", "CADD_annotation", "CADD-scripts")
GERP_DIR   <- file.path(BASE_DIR, "WGS_annotation", "GERP_annotation")
JARVIS_DIR <- file.path(BASE_DIR, "WGS_annotation", "JARVIS")
QC_DIR     <- file.path(BASE_DIR, "WGS_annotation", "WGS_pvar", "download_WGS_pvar")

# ncRNA exon aggregate subsets
NCRNA_SUBMASK_DIR <- file.path(BASE_DIR, "RV_aggregates", "ncRNA_submask")

# ---- Output directory ----
OUT_DIR <- file.path(BASE_DIR, "ncRNA_mask", paste0("chr", chr))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1) Inputs: score files + QC-passed variants
# ------------------------------------------------------------------------------

file_CADD_SNV   <- file.path(CADD_DIR,   paste0("CADD.cutoff.20.SNVs.chr", chr, ".tsv"))
file_CADD_INDEL <- file.path(CADD_DIR,   paste0("CADD.cutoff.20.indels.chr", chr, ".tsv"))
file_GERP       <- file.path(GERP_DIR,   paste0("GERP.cutoff.2.chr", chr, ".tsv"))
file_JARVIS     <- file.path(JARVIS_DIR, paste0("jarvis_ranked_deduplicated.", chr, ".cutoff.0.99.hg38.tsv"))

file_WGSQCpassed <- file.path(QC_DIR, paste0("WGS.chr", chr, ".QCpassed.varlist.vcf"))

ncRNA_file <- file.path(
  NCRNA_SUBMASK_DIR,
  paste0("chr", chr),
  paste0("ncRNA.exons_chr", chr, "_part_", subset_id, ".txt")
)

in_files <- c(file_CADD_SNV, file_CADD_INDEL, file_GERP, file_JARVIS, file_WGSQCpassed, ncRNA_file)
missing <- in_files[!file.exists(in_files)]
if (length(missing) > 0) {
  stop("Missing input file(s):\n  ", paste(missing, collapse = "\n  "), call. = FALSE)
}

# ---- Read score tables (no header) ----
CADD_SNV_data   <- fread(file_CADD_SNV, header = FALSE)
CADD_INDEL_data <- fread(file_CADD_INDEL, header = FALSE)
GERP_data       <- fread(file_GERP, header = FALSE)
JARVIS_data     <- fread(file_JARVIS, header = FALSE)

# ---- Read QC-passed variants (first 5 columns only) ----
WGS_QCpassed <- fread(file_WGSQCpassed, skip = "#CHROM", select = 1:5, header = TRUE, sep = "\t")
setnames(WGS_QCpassed, 1, "CHR")

# ------------------------------------------------------------------------------
# 2) Normalize IDs + merge scores into QC-passed variant table
# ------------------------------------------------------------------------------

# CADD tables: build DRAGEN-style ID
CADD_SNV_data[,   ID := paste0("DRAGEN:chr", V1, ":", V2, ":", V3, ":", V4)]
CADD_INDEL_data[, ID := paste0("DRAGEN:chr", V1, ":", V2, ":", V3, ":", V4)]
setnames(CADD_SNV_data,   "V6", "CADD_PHRED")
setnames(CADD_INDEL_data, "V6", "CADD_PHRED")

CADD_data <- rbind(
  CADD_SNV_data[,   .(ID, CADD_PHRED)],
  CADD_INDEL_data[, .(ID, CADD_PHRED)]
)
WGS_QCpassed <- as.data.frame(WGS_QCpassed) %>%
  left_join(as.data.frame(CADD_data), by = "ID")

CADD_percent_above_20 <- WGS_QCpassed %>%
  summarise(
    total_variants = n(),
    count_above_20 = sum(CADD_PHRED >= 20, na.rm = TRUE),
    percentage = count_above_20 / total_variants
  )
message("[QC] CADD PHRED >= 20")
print(CADD_percent_above_20)

# GERP: V3 -> ID, V6 -> GERP
setnames(GERP_data, old = c("V3", "V6"), new = c("ID", "GERP"))
WGS_QCpassed <- WGS_QCpassed %>%
  left_join(as.data.frame(GERP_data[, .(ID, GERP)]), by = "ID")

GERP_percent_above_2 <- WGS_QCpassed %>%
  summarise(
    total_variants = n(),
    count_above_2 = sum(GERP >= 2, na.rm = TRUE),
    percentage = count_above_2 / total_variants
  )
message("[QC] GERP >= 2")
print(GERP_percent_above_2)

# JARVIS: V1->CHR, V2->POS, V3->JARVIS_score
setnames(JARVIS_data, old = c("V1", "V2", "V3"), new = c("CHR", "POS", "JARVIS_score"))
jarvis_df <- as.data.frame(JARVIS_data[, .(CHR, POS, JARVIS_score)])
WGS_QCpassed <- WGS_QCpassed %>%
  left_join(jarvis_df, by = c("CHR", "POS"))

JARVIS_percent_above_0.99 <- WGS_QCpassed %>%
  summarise(
    total_variants = n(),
    count_above_0.99 = sum(JARVIS_score >= 0.99, na.rm = TRUE),
    percentage = count_above_0.99 / total_variants
  )
message("[QC] JARVIS_score >= 0.99")
print(JARVIS_percent_above_0.99)

# ------------------------------------------------------------------------------
# 3) Assign mask per variant
# ------------------------------------------------------------------------------

get_mask <- function(cadd, gerp, jarvis) {
  if (is.na(cadd) & is.na(gerp) & is.na(jarvis)) return("Null")

  mask_parts <- c()
  if (!is.na(cadd) && cadd >= 20)  mask_parts <- c(mask_parts, "CADD")
  if (!is.na(gerp) && gerp >= 2)   mask_parts <- c(mask_parts, "GERP")
  if (!is.na(jarvis) && jarvis >= 0.99) mask_parts <- c(mask_parts, "JARVIS")

  if (length(mask_parts) == 0) return("None")
  paste(mask_parts, collapse = "_")
}

WGS_QCpassed <- WGS_QCpassed %>%
  mutate(mask = mapply(get_mask, CADD_PHRED, GERP, JARVIS_score))

mask_counts <- WGS_QCpassed %>%
  group_by(mask) %>%
  summarise(count = n(), .groups = "drop")
print(mask_counts)

# Consistency checks: cutoff counts vs mask membership
num_cadd_variants   <- sum(WGS_QCpassed$CADD_PHRED >= 20, na.rm = TRUE)
num_gerp_variants   <- sum(WGS_QCpassed$GERP >= 2, na.rm = TRUE)
num_jarvis_variants <- sum(WGS_QCpassed$JARVIS_score >= 0.99, na.rm = TRUE)

num_mask_CADD   <- sum(str_detect(WGS_QCpassed$mask, "CADD"))
num_mask_GERP   <- sum(str_detect(WGS_QCpassed$mask, "GERP"))
num_mask_JARVIS <- sum(str_detect(WGS_QCpassed$mask, "JARVIS"))

message("[CHECK] Cutoff vs mask membership")
cat("  CADD:   score>=20 ", num_cadd_variants, " | mask contains CADD   ", num_mask_CADD, "\n")
cat("  GERP:   score>=2  ", num_gerp_variants, " | mask contains GERP   ", num_mask_GERP, "\n")
cat("  JARVIS: score>=0.99 ", num_jarvis_variants, " | mask contains JARVIS ", num_mask_JARVIS, "\n")

# ------------------------------------------------------------------------------
# 4) Read ncRNA exon aggregates (subset) + generate SAIGE-style wide tables
# ------------------------------------------------------------------------------

ncRNA_aggregates <- fread(ncRNA_file, select = 1:2, header = TRUE)
message("[INFO] Loaded ncRNA exon aggregates: ", nrow(ncRNA_aggregates), " rows")

# Expand pipe-separated variant IDs and join masks
ncRNA_expanded <- ncRNA_aggregates %>%
  separate_rows(variants, sep = "\\|") %>%
  rename(var = variants) %>%
  left_join(WGS_QCpassed %>% select(ID, mask), by = c("var" = "ID")) %>%
  group_by(ncRNA_uniqueID) %>%
  mutate(v_order = row_number()) %>%
  ungroup()

# One row per aggregate: variant IDs
var_wide <- ncRNA_expanded %>%
  select(ncRNA_uniqueID, v_order, var) %>%
  pivot_wider(names_from = v_order, values_from = var, names_prefix = "ID") %>%
  mutate(row_type = "var") %>%
  select(ncRNA_uniqueID, row_type, everything())

# One row per aggregate: mask labels
mask_wide <- ncRNA_expanded %>%
  select(ncRNA_uniqueID, v_order, mask) %>%
  pivot_wider(names_from = v_order, values_from = mask, names_prefix = "Mask") %>%
  mutate(row_type = "anno") %>%
  select(ncRNA_uniqueID, row_type, everything())

# ------------------------------------------------------------------------------
# 5) Write outputs (space-delimited; no header)
# ------------------------------------------------------------------------------

var_wide_file  <- file.path(OUT_DIR,  paste0("var_wide.ncRNA.exons.chr", chr, ".subset", subset_id, ".txt"))
mask_wide_file <- file.path(OUT_DIR, paste0("mask_wide.ncRNA.exons.chr", chr, ".subset", subset_id, ".txt"))

fwrite(var_wide,  file = var_wide_file,  sep = " ", quote = FALSE, col.names = FALSE)
fwrite(mask_wide, file = mask_wide_file, sep = " ", quote = FALSE, col.names = FALSE)

message("✅ ncRNA mask (exons only) created for chr", chr, ", subset ", subset_id)
message("[INFO] Wrote:\n  ", var_wide_file, "\n  ", mask_wide_file)