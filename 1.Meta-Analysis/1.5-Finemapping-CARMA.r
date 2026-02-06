#!/usr/bin/env Rscript
############### Fine-mapping pipeline —CARMA for meta-GWAS ################
############### Author: Shelley ###########################################
############### Date:   2025-10-04 ########################################
##
## Purpose
##   Run CARMA fine-mapping for a single locus using:
##     (1) per-locus GWAS summary statistics
##     (2) per-locus LD matrix + variant IDs
##
## Key steps
##   1) Read per-locus GWAS + LD inputs
##   2) Align GWAS rows to the LD variant order (exact 1:1 match)
##   3) Filter variants with missing LD values (two-stage NaN filtering)
##   4) Compute Z = Effect / StdErr
##   5) Run CARMA with outlier detection enabled
##   6) Save results + QC outputs for reproducibility
##
## Usage
##   Rscript run_carma.R <locus_id>
##
###############################################################################

##  Load packages
suppressPackageStartupMessages({
  library(data.table)
  library(magrittr)
  library(R.utils)
  library(dplyr)
  library(CARMA)
})

##  args: locus_id from command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: Rscript run_carma.R <locus_id>")
locus_id <- args[[1]]
# check current locus id
cat(sprintf("[INFO] locus_id = %s\n", locus_id))

## Make per-locus working directory and cd into it ----------------------------
base_wd  <- "/public/home/hpc8301200407/WGS/Melanoma_WGS/GWAS_Finemapping/Finemapping-CARMA/4.CARMA"
locus_dir <- file.path(base_wd, paste0("locusA_", locus_id))

dir.create(locus_dir, recursive = TRUE, showWarnings = FALSE)
setwd(locus_dir)

cat(sprintf("[INFO] Working directory set to: %s\n", getwd()))

## paths
locus_gwas <- sprintf("/public/home/hpc8301200407/WGS/Melanoma_WGS/GWAS_Finemapping/Finemapping-CARMA/1.Extract_GWAS_loci/per_locus.filtered/locusA_%s.txt", locus_id)
locus_ld   <- sprintf("/public/home/hpc8301200407/WGS/Melanoma_WGS/GWAS_Finemapping/Finemapping-CARMA/3.LD/locusA_%s.LD.ld", locus_id)
locus_ids  <- sprintf("/public/home/hpc8301200407/WGS/Melanoma_WGS/GWAS_Finemapping/Finemapping-CARMA/3.LD/locusA_%s.LD.tags.list", locus_id)
locus_bim  <- sprintf("/public/home/hpc8301200407/WGS/Melanoma_WGS/GWAS_Finemapping/Finemapping-CARMA/2.subset.GT/BED_subset/locusA_%s.bim", locus_id)

## read (no colname surgery)
gwas_dt <- data.table::fread(locus_gwas)
ld_ids  <- data.table::fread(locus_ids, header = FALSE)$V1
ld_mat  <- as.matrix(data.table::fread(locus_ld, header = FALSE))
bim_dt  <- data.table::fread(locus_bim, header = FALSE)

## quick summaries
cat(sprintf("[OK] GWAS rows: %d\n", nrow(gwas_dt)))
cat(sprintf("[OK] LD ids   : %d\n", length(ld_ids)))
cat(sprintf("[OK] LD matrix: %d x %d\n", nrow(ld_mat), ncol(ld_mat)))
cat(sprintf("[OK] BIM rows : %d\n", nrow(bim_dt)))

# print and remove the first row of ld_ids which is "SNP"
## peek
cat("[HEAD] ld_ids (first 5):", paste(utils::head(ld_ids, 5), collapse = ", "), "\n")
cat("[INFO] first ld_id:", ld_ids[1], "\n")

## if the first token is a header "SNP", drop it
if (identical(as.character(ld_ids[1]), "SNP")) {
  ld_ids <- ld_ids[-1]
  cat("[FIX] Removed leading 'SNP' header from ld_ids\n")
}

## optional: ensure length matches LD matrix dimension
if (length(ld_ids) != nrow(ld_mat)) {
  warning(sprintf("ld_ids length (%d) != LD matrix nrow (%d)", length(ld_ids), nrow(ld_mat)))
} else {
  cat("[OK] ld_ids now matches LD matrix dimension\n")
}

# --- sanity check 1: BIM order vs LD ids -------------------------------------
ids_bim <- as.character(bim_dt$V2)
ids_ld  <- as.character(ld_ids)

cat(sprintf("[INFO] lengths: BIM=%d  LD_ids=%d\n", length(ids_bim), length(ids_ld)))

# quick shape checks
if (nrow(ld_mat) != length(ids_ld)) {
  warning(sprintf("LD matrix nrow (%d) != ld_ids length (%d)", nrow(ld_mat), length(ids_ld)))
}
if (length(ids_bim) != length(ids_ld)) {
  warning("BIM and ld_ids differ in length; order equality cannot hold exactly.")
}

# exact order check
mismatch_idx <- which(ids_bim != ids_ld)
if (length(mismatch_idx) == 0 &&
    length(ids_bim) == length(ids_ld)) {
  cat("[OK] BIM$V2 order matches ld_ids exactly.\n")
} else {
  cat(sprintf("[WARN] Order mismatch: %d positions differ.\n", length(mismatch_idx)))
  # show first few mismatches
  if (length(mismatch_idx) > 0) {
    k <- min(10L, length(mismatch_idx))
    idx <- mismatch_idx[seq_len(k)]
    mis_tb <- data.frame(
      index = idx,
      BIM_V2 = ids_bim[idx],
      LD_ID  = ids_ld[idx],
      stringsAsFactors = FALSE
    )
    print(mis_tb, row.names = FALSE)
  }
  # also check if the sets are the same (permutation vs. true mismatch)
  same_set <- setequal(ids_bim, ids_ld)
  cat(sprintf("[INFO] Same set (ignoring order)? %s\n", ifelse(same_set, "YES", "NO")))
  if (same_set) {
    cat("[HINT] They contain the same variants but in different order.\n")
    cat("       If needed later, you can realign with: ord <- match(ids_ld, ids_bim)\n")
  }
}

# if the sanity check passed, reorder the gwas_dt by its UKBB_WGS_ID column as the order of variants in ids_ld
## Reorder GWAS to LD order (all LD IDs are present) --------------------------
stopifnot("UKBB_WGS_ID" %in% names(gwas_dt))

ids_ld <- as.character(ld_ids)               # authoritative LD order
ids_gw <- as.character(gwas_dt$UKBB_WGS_ID)

# For each LD ID, find its row in GWAS; assert none are missing
ord <- match(ids_ld, ids_gw)
stopifnot(!any(is.na(ord)))                  # all LD IDs must be in GWAS

# Subset GWAS to LD order (no filtering needed since all exist)
gwas_dt_ld <- gwas_dt[ord, , drop = FALSE]

# Sanity: exact 1:1 equality to ids_ld
stopifnot(identical(as.character(gwas_dt_ld$UKBB_WGS_ID), ids_ld))

cat(sprintf("[OK] GWAS aligned to LD: %d variants (all LD IDs found in GWAS)\n",
            nrow(gwas_dt_ld)))

## Build LD submatrix exactly matching the reordered GWAS ---------------------
# Since all LD IDs exist and the LD matrix is in ids_ld order, just take all rows/cols.
keep_ld     <- seq_along(ids_ld)
ld_mat_sub0 <- ld_mat[keep_ld, keep_ld, drop = FALSE]
ids_sub0    <- ids_ld[keep_ld]
gwas_sub0   <- gwas_dt_ld

stopifnot(
  nrow(ld_mat_sub0) == nrow(gwas_sub0),
  ncol(ld_mat_sub0) == nrow(gwas_sub0),
  identical(as.character(gwas_sub0$UKBB_WGS_ID), ids_sub0)
)

## Stage 1: drop variants with >=1% NaN in their row OR column ----------------
frac_nan     <- function(v) mean(is.nan(v))
row_nan_frac <- apply(ld_mat_sub0, 1, frac_nan)
col_nan_frac <- apply(ld_mat_sub0, 2, frac_nan)
bad01        <- (row_nan_frac >= 0.01) | (col_nan_frac >= 0.01)

failed01_idx <- which(bad01)
keep1_idx    <- which(!bad01)
failed01_ids <- ids_sub0[failed01_idx]

if (length(failed01_idx) > 0) {
  cat(sprintf("[FILTER-1%%] dropping %d/%d variants (>=1%% NaN in row/col). Example: %s\n",
              length(failed01_idx), length(ids_sub0),
              paste(utils::head(failed01_ids, 5), collapse = ", ")))
}

# Survivors of stage 1
ld_mat_sub1 <- ld_mat_sub0[keep1_idx, keep1_idx, drop = FALSE]
ids_sub1    <- ids_sub0[keep1_idx]
gwas_sub1   <- gwas_sub0[keep1_idx]

## Stage 2: among survivors, drop any variant with ANY NaN --------------------
row_has_nan1 <- apply(ld_mat_sub1, 1, function(v) any(is.nan(v)))
col_has_nan1 <- apply(ld_mat_sub1, 2, function(v) any(is.nan(v)))
bad_any1     <- row_has_nan1 | col_has_nan1

failedAny_idx <- which(bad_any1)
keep2_idx     <- which(!bad_any1)
failedAny_ids <- ids_sub1[failedAny_idx]

if (length(failedAny_idx) > 0) {
  cat(sprintf("[FILTER-any] dropping %d/%d additional variants (any NaN). Example: %s\n",
              length(failedAny_idx), length(ids_sub1),
              paste(utils::head(failedAny_ids, 5), collapse = ", ")))
}

## Final kept set after both stages -------------------------------------------
fixed_ld_mat     <- ld_mat_sub1[keep2_idx, keep2_idx, drop = FALSE]
fixed_gwas_dt_ld <- gwas_sub1[keep2_idx]
fixed_ids_keep   <- ids_sub1[keep2_idx]

# Collate failed lists (stage-1 + stage-2)
failed_ids <- c(failed01_ids, failedAny_ids)
gwas_failed <- data.table::rbindlist(
  list(
    gwas_sub0[failed01_idx],
    if (length(failedAny_idx) > 0) gwas_sub1[failedAny_idx] else gwas_sub1[0]
  ),
  use.names = TRUE, fill = TRUE
)

cat(sprintf("[FIX] kept %d variants; dropped %d total (stage1=%d, stage2=%d).\n",
            nrow(fixed_ld_mat), length(failed_ids),
            length(failed01_ids), length(failedAny_ids)))

# gwas_pass should match the final kept set:
gwas_pass <- fixed_gwas_dt_ld

## 4) Calculate Z-score and drop non-finite z (keeps alignment)
stopifnot(all(c("Effect","StdErr") %in% names(fixed_gwas_dt_ld)))
fixed_gwas_dt_ld[, z := Effect / StdErr]

# 5) write GWAS split and ID lists to disk (write the z-added passed table)
tag <- if (exists("locus_id")) paste0(".locusA_", locus_id) else ""
out_failed_gwas <- sprintf("GWAS_failed_NaN%s.tsv", tag)
out_passed_gwas <- sprintf("GWAS_passed%s.tsv", tag)
out_failed_ids  <- sprintf("variant_ids_failedNaN%s.txt", tag)
out_passed_ids  <- sprintf("variant_ids_passed%s.txt", tag)

data.table::fwrite(gwas_failed,        out_failed_gwas, sep = "\t", quote = FALSE)
data.table::fwrite(fixed_gwas_dt_ld,   out_passed_gwas, sep = "\t", quote = FALSE)
writeLines(as.character(failed_ids),     con = out_failed_ids)
writeLines(as.character(fixed_ids_keep), con = out_passed_ids)

cat("[OK] wrote outputs:\n")
cat("   - ", out_failed_gwas, " (", nrow(gwas_failed), " rows)\n", sep = "")
cat("   - ", out_passed_gwas, " (", nrow(fixed_gwas_dt_ld),  " rows)\n", sep = "")
cat("   - ", out_failed_ids,  " (", length(failed_ids), " IDs)\n", sep = "")
cat("   - ", out_passed_ids,  " (", length(fixed_ids_keep), " IDs)\n", sep = "")

# 6) final sanity before CARMA
stopifnot(
  is.matrix(fixed_ld_mat),
  nrow(fixed_ld_mat) == ncol(fixed_ld_mat),
  nrow(fixed_ld_mat) == nrow(fixed_gwas_dt_ld),
  length(fixed_ids_keep) == nrow(fixed_ld_mat),
  all(is.finite(diag(fixed_ld_mat)))
)

# 7) run CARMA
z.list <- list()
ld.list <- list()
lambda.list <- list()
# lambda: as suggested in the paper, set as 1

z.list[[1]] <- fixed_gwas_dt_ld$z
ld.list[[1]] <- as.matrix(fixed_ld_mat)
lambda.list[[1]] <- 1

# PIP is set as 0.99 by default
CARMA_results <- CARMA(z.list, ld.list, lambda.list=lambda.list, outlier.switch=T)

# save the CARMA_results as rds out 
out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
tag <- if (exists("locus_id")) paste0(".locusA_", locus_id) else ""
out_rds <- file.path(out_dir, sprintf("CARMA_results%s.rds", tag))

saveRDS(CARMA_results, out_rds)
cat("[OK] CARMA results saved to: ", out_rds, "\n", sep = "")