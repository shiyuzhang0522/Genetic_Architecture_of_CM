#!/usr/bin/env Rscript
# Create mask for melanocyte-specific cCREs 
############################################################
## R script: Build 3-score variant masks and generate
##           SAIGE-style group/annotation tables per aggregate
##
## Inputs (per chromosome):
##   1) QC-passed variant list (VCF; read first 5 columns)
##   2) CADD (SNV + INDEL), GERP, JARVIS score tables
##   3) Aggregate definitions (e.g., melanocyte cCRE aggregates)
##
## Outputs (per chromosome):
##   1) QC-passed variants with mask labels (tab-delimited)
##   2) SAIGE group/anno table for aggregates (space-delimited, no header)
##
## Author: Shiyu Zhang (Shelley)
## Date: 2025-04-12
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
if (length(args) < 1) {
  stop("Usage: Rscript build_masks_and_groups.R <chr> [--base <BASE_DIR>]", call. = FALSE)
}

chr <- args[1]
if (!grepl("^([1-9]|1[0-9]|2[0-2])$", chr)) {
  stop("chr must be an integer in [1..22].", call. = FALSE)
}
message("[INFO] Processing chromosome: ", chr)

# ---- Base directory (HIDE absolute paths) ----
# Users should set this via environment variable or edit once.
BASE_DIR <- Sys.getenv("WGS_MELANOMA_BASE", unset = "PATH_TO_PROJECT_ROOT")

# ---- Input directories (relative to BASE_DIR) ----
CADD_DIR   <- file.path(BASE_DIR, "WGS_annotation", "CADD_annotation", "CADD-scripts")
GERP_DIR   <- file.path(BASE_DIR, "WGS_annotation", "GERP_annotation")
JARVIS_DIR <- file.path(BASE_DIR, "WGS_annotation", "JARVIS")
QC_DIR     <- file.path(BASE_DIR, "WGS_annotation", "WGS_pvar", "download_WGS_pvar")
AGG_DIR    <- file.path(BASE_DIR, "RV_aggregates")

# ---- Output directory ----
OUT_DIR <- file.path(BASE_DIR, "reordered_masks")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1) Inputs: score files + QC-passed variants
# ------------------------------------------------------------------------------

file_CADD_SNV   <- file.path(CADD_DIR,   paste0("CADD.cutoff.20.SNVs.chr", chr, ".tsv"))
file_CADD_INDEL <- file.path(CADD_DIR,   paste0("CADD.cutoff.20.indels.chr", chr, ".tsv"))
file_GERP       <- file.path(GERP_DIR,   paste0("GERP.cutoff.2.chr", chr, ".tsv"))
file_JARVIS     <- file.path(JARVIS_DIR, paste0("jarvis_ranked_deduplicated.", chr, ".cutoff.0.99.hg38.tsv"))

file_WGSQCpassed <- file.path(QC_DIR, paste0("WGS.chr", chr, ".QCpassed.varlist.vcf"))

in_files <- c(file_CADD_SNV, file_CADD_INDEL, file_GERP, file_JARVIS, file_WGSQCpassed)
missing <- in_files[!file.exists(in_files)]
if (length(missing) > 0) stop("Missing input file(s):\n  ", paste(missing, collapse = "\n  "), call. = FALSE)

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
CADD_data_df <- as.data.frame(CADD_data)

WGS_QCpassed <- WGS_QCpassed %>%
  left_join(CADD_data_df, by = "ID")

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
GERP_keep_df <- as.data.frame(GERP_data[, .(ID, GERP)])

WGS_QCpassed <- WGS_QCpassed %>%
  left_join(GERP_keep_df, by = "ID")

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
  # all NA -> Null
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

total_assigned <- sum(mask_counts$count)
total_variants <- nrow(WGS_QCpassed)
if (total_assigned == total_variants) {
  message("[CHECK] Mask counts sum to total variants: ", total_assigned)
} else {
  message("[WARN] Mask count mismatch: ", total_assigned, " vs total ", total_variants)
}

# Consistency checks: mask membership should match cutoff counts
num_cadd_variants   <- sum(WGS_QCpassed$CADD_PHRED >= 20, na.rm = TRUE)
num_gerp_variants   <- sum(WGS_QCpassed$GERP      >= 2,  na.rm = TRUE)
num_jarvis_variants <- sum(WGS_QCpassed$JARVIS_score >= 0.99, na.rm = TRUE)

num_mask_CADD   <- sum(str_detect(WGS_QCpassed$mask, "CADD"))
num_mask_GERP   <- sum(str_detect(WGS_QCpassed$mask, "GERP"))
num_mask_JARVIS <- sum(str_detect(WGS_QCpassed$mask, "JARVIS"))

message("[CHECK] Cutoff vs mask membership")
cat("  CADD:   score>=20 ", num_cadd_variants, " | mask contains CADD   ", num_mask_CADD, "\n")
cat("  GERP:   score>=2  ", num_gerp_variants, " | mask contains GERP   ", num_mask_GERP, "\n")
cat("  JARVIS: score>=0.99 ", num_jarvis_variants, " | mask contains JARVIS ", num_mask_JARVIS, "\n")

# Sort by CHR, POS
setDT(WGS_QCpassed)
setorder(WGS_QCpassed, CHR, POS)

out_mask_table <- file.path(OUT_DIR, paste0("WGS_QCpassed.chr", chr, ".three.score.mask.txt"))
fwrite(WGS_QCpassed, file = out_mask_table, sep = "\t")
message("[INFO] Wrote: ", out_mask_table)

# ------------------------------------------------------------------------------
# 4) Aggregates: melanocyte cCRE -> SAIGE-like wide var/anno
# ------------------------------------------------------------------------------

mel_cCRE_aggregates_file <- file.path(
  AGG_DIR, paste0("chr", chr),
  paste0("mel_cCRE_aggregates.chr", chr, ".txt")
)
if (!file.exists(mel_cCRE_aggregates_file)) {
  stop("Missing aggregate file: ", mel_cCRE_aggregates_file, call. = FALSE)
}

mel_cCRE_aggregates <- fread(mel_cCRE_aggregates_file, select = 1:2, header = TRUE)
original_order <- mel_cCRE_aggregates$mel_cCRE_uniqueID

mel_expanded <- mel_cCRE_aggregates %>%
  separate_rows(variants, sep = "\\|") %>%
  rename(var = variants) %>%
  left_join(as.data.frame(WGS_QCpassed[, .(ID, mask)]), by = c("var" = "ID")) %>%
  mutate(mel_cCRE_uniqueID = factor(mel_cCRE_uniqueID, levels = original_order)) %>%
  group_by(mel_cCRE_uniqueID) %>%
  mutate(v_order = row_number()) %>%
  ungroup() %>%
  arrange(mel_cCRE_uniqueID)

var_wide <- mel_expanded %>%
  select(mel_cCRE_uniqueID, v_order, var) %>%
  pivot_wider(names_from = v_order, values_from = var, names_prefix = "ID") %>%
  mutate(row_type = "var") %>%
  select(mel_cCRE_uniqueID, row_type, everything())

mask_wide <- mel_expanded %>%
  select(mel_cCRE_uniqueID, v_order, mask) %>%
  pivot_wider(names_from = v_order, values_from = mask, names_prefix = "Mask") %>%
  mutate(row_type = "anno") %>%
  select(mel_cCRE_uniqueID, row_type, everything())

var_wide  <- as.data.table(var_wide)
mask_wide <- as.data.table(mask_wide)

rename_common_cols <- function(dt, prefix = "V") {
  setnames(dt, c("mel_cCRE_uniqueID", "row_type", paste0(prefix, seq_len(ncol(dt) - 2))))
  dt
}

combined_wide <- rbind(
  rename_common_cols(var_wide,  prefix = "V"),
  rename_common_cols(mask_wide, prefix = "V")
)

combined_wide$mel_cCRE_uniqueID <- factor(combined_wide$mel_cCRE_uniqueID, levels = original_order)
combined_wide$row_type <- factor(combined_wide$row_type, levels = c("var", "anno"))
combined_wide <- combined_wide[order(combined_wide$mel_cCRE_uniqueID, combined_wide$row_type), ]

out_group <- file.path(OUT_DIR, paste0("combined_mask.mel.cCRE.chr", chr, ".txt"))
fwrite(combined_wide, file = out_group, sep = " ", quote = FALSE, col.names = FALSE)

message("✅ melanocyte cCRE mask has been successfully created for chr", chr)
message("[INFO] Wrote: ", out_group)