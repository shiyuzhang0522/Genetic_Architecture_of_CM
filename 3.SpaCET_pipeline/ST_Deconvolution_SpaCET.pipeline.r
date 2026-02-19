############################################################
## SpaCET Pipeline for Spatial Transcriptomics Deconvolution
## Author: Shiyu Zhang
## Last updated: 2026-02-10
############################################################

############################################################
## Step 1) Load packages
############################################################
suppressPackageStartupMessages({
  library(SpaCET)
  library(data.table)
  library(Matrix)
})

############################################################
## Step 2) Read ST data and create SpaCET object (10X Visium)
############################################################
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript run_spacet.R <sample_ST_ID>\n")
  cat("Example: Rscript run_spacet.R MEL162YJG_CA21_4_30\n")
  quit(save = "no", status = 1)
}

sample_ST_ID <- args[1]

base_dir   <- "/public/home/hpc8301200407/MEL.ST.10x.Visium"
visiumPath <- file.path(base_dir, sample_ST_ID)

# ---- output dir (per sample) ----
BASE_OUTDIR <- "/public/home/hpc8301200407/WGS/GWAS_ST/SpaCET_process_updated_20260210"
outdir <- file.path(BASE_OUTDIR, sample_ST_ID)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "tables"), showWarnings = FALSE)
dir.create(file.path(outdir, "rds"), showWarnings = FALSE)
dir.create(file.path(outdir, "log"), showWarnings = FALSE)
dir.create(file.path(outdir, "fig"), showWarnings = FALSE)

# ---- sanity checks ----
if (!dir.exists(visiumPath)) {
  stop("[ERROR] visiumPath not found: ", visiumPath)
}

cat("[INFO] sample_ST_ID: ", sample_ST_ID, "\n", sep = "")
cat("[INFO] visiumPath:   ", visiumPath, "\n", sep = "")
cat("[INFO] outdir:       ", outdir, "\n", sep = "")

# ---- create SpaCET object from 10X Visium folder ----
SpaCET_obj <- create.SpaCET.object.10X(visiumPath = visiumPath)

# ---- robust accessors (matches your SpaCET object structure) ----
stopifnot("input" %in% slotNames(SpaCET_obj))
stopifnot(all(c("counts", "spotCoordinates", "metaData") %in% names(SpaCET_obj@input)))

counts <- SpaCET_obj@input$counts                 # dgCMatrix, genes x spots
coords <- SpaCET_obj@input$spotCoordinates        # data.frame, spots x coord_cols
meta   <- SpaCET_obj@input$metaData               # data.frame, spots x meta_cols

stopifnot(inherits(counts, "dgCMatrix"))
stopifnot(ncol(counts) == nrow(coords))

cat("[INFO] SpaCET object created for sample: ", sample_ST_ID, "\n", sep = "")
cat("[INFO] Genes x Spots: ", nrow(counts), " x ", ncol(counts), "\n", sep = "")

############################################################
## Step 3) Extract per-spot QC metrics: UMI + nGene
############################################################
spot_ids <- colnames(counts)

# UMI (library size) per spot: sum of counts
umi_vec <- Matrix::colSums(counts)

# nGene per spot: number of genes with >=1 UMI
ngene_vec <- Matrix::colSums(counts > 0)

qc_dt <- data.table(
  spot_id   = spot_ids,
  barcode   = if ("barcode" %in% colnames(meta)) meta$barcode else NA_character_,
  UMI       = as.numeric(umi_vec),
  nGene     = as.numeric(ngene_vec),
  log1p_UMI = log1p(as.numeric(umi_vec))
)

# add coordinates (keep original column names from SpaCET)
coord_dt <- data.table::as.data.table(coords, keep.rownames = FALSE)
if (nrow(coord_dt) == nrow(qc_dt)) {
  qc_dt <- cbind(qc_dt, coord_dt)
}

# save
qc_outfile <- file.path(outdir, "tables", paste0(sample_ST_ID, ".spot_QC_metrics.tsv"))
data.table::fwrite(qc_dt, qc_outfile, sep = "\t")

cat("[INFO] Saved spot QC metrics: ", qc_outfile, "\n", sep = "")

############################################################
## Step 4) QC filtering + QC visualization (SpaCET)
############################################################

# ---- 4.1 QC filter: remove low-complexity spots ----
SpaCET_obj <- SpaCET.quality.control(SpaCET_obj, min.genes = 100)

cat("[INFO] QC filtering done (min.genes = 100).\n")

# ---- 4.2 QC plot (UMI + Gene) ----
fig.QC <- SpaCET.visualize.spatialFeature(
  SpaCET_obj,
  spatialType     = "QualityControl",
  spatialFeatures = c("UMI", "Gene"),
  imageBg         = TRUE,
  imageSize = "CaptureArea",
  pointSize = 0.4
)

# ---- 4.3 Save QC figure as Nature-friendly PDF ----
qc_pdf <- file.path(outdir, "fig", paste0(sample_ST_ID, ".QC_UMI_Gene.minGenes100.pdf"))

# A4-friendly small panel (single column style):
# - width ~ 7.0 in, height ~ 3.3–3.8 in works well for 2 panels side-by-side
pdf(qc_pdf, width = 7.2, height = 2.4, useDingbats = FALSE)
print(fig.QC)
dev.off()

cat("[INFO] Saved QC PDF: ", qc_pdf, "\n", sep = "")

############################################################
## Step 5) Deconvolution + save key outputs (prop matrix + signature genes)
############################################################

# ---- 5.1 Run deconvolution ----
SpaCET_obj <- SpaCET.deconvolution(SpaCET_obj, cancerType = "SKCM", coreNo = 8)
cat("[INFO] SpaCET deconvolution finished.\n")

# ---- 5.2 Robustly extract results ----
stopifnot("results" %in% slotNames(SpaCET_obj))
stopifnot("deconvolution" %in% names(SpaCET_obj@results))

deconv <- SpaCET_obj@results$deconvolution
stopifnot("propMat" %in% names(deconv))

propMat <- deconv$propMat  # numeric matrix: celltypes x spots (per your str())

# sanity checks
stopifnot(is.matrix(propMat))
stopifnot(!is.null(rownames(propMat)))
stopifnot(!is.null(colnames(propMat)))

cat("[INFO] propMat dim (celltypes x spots): ",
    nrow(propMat), " x ", ncol(propMat), "\n", sep = "")

# ---- 5.3 Save cell-type proportion matrix ----
# Save as: spots x celltypes (more convenient for joins with qc_dt)
prop_dt <- data.table::as.data.table(t(propMat), keep.rownames = "spot_id")
prop_outfile <- file.path(outdir, "tables", paste0(sample_ST_ID, ".SpaCET.celltype_prop.tsv"))
data.table::fwrite(prop_dt, prop_outfile, sep = "\t")
cat("[INFO] Saved cell-type proportions: ", prop_outfile, "\n", sep = "")

# (optional) also save long/tidy format
prop_long_dt <- melt(
  prop_dt,
  id.vars = "spot_id",
  variable.name = "cell_type",
  value.name = "proportion"
)
prop_long_outfile <- file.path(outdir, "tables", paste0(sample_ST_ID, ".SpaCET.celltype_prop.long.tsv"))
data.table::fwrite(prop_long_dt, prop_long_outfile, sep = "\t")
cat("[INFO] Saved cell-type proportions (long): ", prop_long_outfile, "\n", sep = "")

# ---- 5.4 Save signature genes per cell type ----
# In your object: deconv$Ref$sigGenes is a named list (celltype -> character vector of genes)
sigGenes <- NULL
if ("Ref" %in% names(deconv) && "sigGenes" %in% names(deconv$Ref)) {
  sigGenes <- deconv$Ref$sigGenes
}

if (!is.null(sigGenes)) {
  stopifnot(is.list(sigGenes), !is.null(names(sigGenes)))

  sig_long_dt <- rbindlist(
    lapply(names(sigGenes), function(ct) {
      data.table(cell_type = ct, gene = as.character(sigGenes[[ct]]))
    }),
    use.names = TRUE, fill = TRUE
  )

  sig_outfile <- file.path(outdir, "tables", paste0(sample_ST_ID, ".SpaCET.signature_genes.tsv"))
  data.table::fwrite(sig_long_dt, sig_outfile, sep = "\t")
  cat("[INFO] Saved signature genes: ", sig_outfile, "\n", sep = "")
} else {
  cat("[WARN] Signature genes not found at SpaCET_obj@results$deconvolution$Ref$sigGenes\n")
}

############################################################
## Step 6) Cell–cell interaction (CCI) — colocalization
############################################################

# ---- 6.1 Run CCI colocalization ----
SpaCET_obj <- SpaCET.CCI.colocalization(SpaCET_obj)
cat("[INFO] SpaCET CCI colocalization finished.\n")

# ---- 6.2 Extract and save colocalization table ----
stopifnot("results" %in% slotNames(SpaCET_obj))
stopifnot("CCI" %in% names(SpaCET_obj@results))
stopifnot("colocalization" %in% names(SpaCET_obj@results$CCI))

coloc_df <- SpaCET_obj@results$CCI$colocalization
stopifnot(is.data.frame(coloc_df))

# standardize to data.table + save
coloc_dt <- data.table::as.data.table(coloc_df)

# print positively colocalized significant pairs
print(coloc_dt[fraction_rho > 0 & fraction_pv < 0.05])

# (optional) enforce expected columns if present
expected_cols <- c(
  "cell_type_1", "cell_type_2",
  "fraction_product", "fraction_rho", "fraction_pv",
  "reference_rho", "reference_pv"
)
missing_cols <- setdiff(expected_cols, colnames(coloc_dt))
if (length(missing_cols) > 0) {
  cat("[WARN] Missing expected columns in colocalization table: ",
      paste(missing_cols, collapse = ","), "\n", sep = "")
}

coloc_outfile <- file.path(outdir, "tables", paste0(sample_ST_ID, ".SpaCET.CCI.colocalization.tsv"))
data.table::fwrite(coloc_dt, coloc_outfile, sep = "\t")
cat("[INFO] Saved CCI colocalization table: ", coloc_outfile, "\n", sep = "")

############################################################
## Step 6.3) CCI — L–R network enrichment score per spot
############################################################

# ---- 6.3.1 Run L–R network scoring ----
SpaCET_obj <- SpaCET.CCI.LRNetworkScore(SpaCET_obj, coreNo = 8)
cat("[INFO] SpaCET CCI LRNetworkScore finished.\n")

# ---- 6.3.2 Extract LRNetworkScore matrix ----
stopifnot("results" %in% slotNames(SpaCET_obj))
stopifnot("CCI" %in% names(SpaCET_obj@results))
stopifnot("LRNetworkScore" %in% names(SpaCET_obj@results$CCI))

lr_mat <- SpaCET_obj@results$CCI$LRNetworkScore
stopifnot(is.matrix(lr_mat))
stopifnot(!is.null(rownames(lr_mat)))
stopifnot(!is.null(colnames(lr_mat)))

# per your str(): rows are metrics, cols are spots
# rows: Raw_expr, Network_Score, Network_Score_pv
cat("[INFO] LRNetworkScore dim (metrics x spots): ",
    nrow(lr_mat), " x ", ncol(lr_mat), "\n", sep = "")
cat("[INFO] LRNetworkScore rows: ", paste(rownames(lr_mat), collapse = ","), "\n", sep = "")

# ---- 6.3.3 Convert to tidy dt (spot-wise) and save ----
lr_dt <- data.table(
  spot_id = colnames(lr_mat),
  Raw_expr = as.numeric(lr_mat["Raw_expr", ]),
  Network_Score = as.numeric(lr_mat["Network_Score", ]),
  Network_Score_pv = as.numeric(lr_mat["Network_Score_pv", ])
)

lr_outfile <- file.path(outdir, "tables", paste0(sample_ST_ID, ".SpaCET.CCI.LRNetworkScore.tsv"))
data.table::fwrite(lr_dt, lr_outfile, sep = "\t")
cat("[INFO] Saved LRNetworkScore table: ", lr_outfile, "\n", sep = "")

############################################################
## Step 7) Tumor–stroma interface calling + export table
############################################################

# ---- 7.1 Identify interface (default: Malignant >= 0.5 => Tumor) ----
SpaCET_obj <- SpaCET.identify.interface(SpaCET_obj)
cat("[INFO] SpaCET identify.interface finished.\n")

# ---- 7.2 Extract interface labels and save as dt ----
stopifnot("results" %in% slotNames(SpaCET_obj))
stopifnot("CCI" %in% names(SpaCET_obj@results))
stopifnot("interface" %in% names(SpaCET_obj@results$CCI))

iface_mat <- SpaCET_obj@results$CCI$interface
stopifnot(is.matrix(iface_mat))
stopifnot(nrow(iface_mat) == 1)

iface_dt <- data.table(
  spot_id = colnames(iface_mat),
  interface_label = as.character(iface_mat[1, ])
)

# quick summary
iface_tab <- iface_dt[, .N, by = interface_label][order(-N)]
print(iface_tab)

iface_outfile <- file.path(outdir, "tables", paste0(sample_ST_ID, ".SpaCET.interface_labels.tsv"))
data.table::fwrite(iface_dt, iface_outfile, sep = "\t")
cat("[INFO] Saved interface labels: ", iface_outfile, "\n", sep = "")

############################################################
## Step 8) Gene set scoring per spot (Hallmark / CancerCellState / TLS)
############################################################

# Helper: extract SpaCET GeneSetScore matrix -> spot-wise dt and save
save_geneset_score <- function(SpaCET_obj, geneset_name, outdir, sample_ST_ID) {

  stopifnot("results" %in% slotNames(SpaCET_obj))
  stopifnot("GeneSetScore" %in% names(SpaCET_obj@results))
  stopifnot(is.matrix(SpaCET_obj@results$GeneSetScore))

  gs_mat <- SpaCET_obj@results$GeneSetScore  # matrix: gene_sets x spots (per your str())

  stopifnot(!is.null(rownames(gs_mat)), !is.null(colnames(gs_mat)))

  # spots x gene_sets (easy join with qc_dt)
  gs_dt <- data.table::as.data.table(t(gs_mat), keep.rownames = "spot_id")

  outfile <- file.path(outdir, "tables",
                       paste0(sample_ST_ID, ".SpaCET.GeneSetScore.", geneset_name, ".tsv"))
  data.table::fwrite(gs_dt, outfile, sep = "\t")

  cat("[INFO] Saved GeneSetScore (", geneset_name, "): ", outfile, "\n", sep = "")
  invisible(outfile)
}

############################################################
## Step 8.1) Hallmark scores
############################################################
SpaCET_obj <- SpaCET.GeneSetScore(SpaCET_obj, GeneSets = "Hallmark")
cat("[INFO] GeneSetScore finished: Hallmark\n")
save_geneset_score(SpaCET_obj, geneset_name = "Hallmark", outdir = outdir, sample_ST_ID = sample_ST_ID)

############################################################
## Step 8.2) Cancer cell state module scores
############################################################
SpaCET_obj <- SpaCET.GeneSetScore(SpaCET_obj, GeneSets = "CancerCellState")
cat("[INFO] GeneSetScore finished: CancerCellState\n")
save_geneset_score(SpaCET_obj, geneset_name = "CancerCellState", outdir = outdir, sample_ST_ID = sample_ST_ID)

############################################################
## Step 8.3) TLS score (30-gene TLS signature)
############################################################
SpaCET_obj <- SpaCET.GeneSetScore(SpaCET_obj, GeneSets = "TLS")
cat("[INFO] GeneSetScore finished: TLS\n")
save_geneset_score(SpaCET_obj, geneset_name = "TLS", outdir = outdir, sample_ST_ID = sample_ST_ID)

### step 9. Compute weight matrix 
W <- calWeights(SpaCET_obj, radius=200, sigma=100, diagAsZero=TRUE)

### step 10. find spatially variable genes
suppressPackageStartupMessages(library(future))
options(future.globals.maxSize = 20 * 1024^3)  # 20 GB

# run for genome-wide genes, 50 mins
SpaCET_obj <- SpaCET.SpatialCorrelation(
  SpaCET_obj, 
  mode="univariate", 
  item=NULL, 
  W=W, 
  nPermutation=1000
)

############################################################
## Step 10) Spatially variable genes (SVGs) — summarize + export
############################################################

# ---- 10.1 Sanity check: spatial correlation results exist ----
stopifnot("results" %in% slotNames(SpaCET_obj))
stopifnot("SpatialCorrelation" %in% names(SpaCET_obj@results))
stopifnot("univariate" %in% names(SpaCET_obj@results$SpatialCorrelation))

svg_df <- SpaCET_obj@results$SpatialCorrelation$univariate
stopifnot(is.data.frame(svg_df))
stopifnot(all(c("p.Moran_I", "p.Moran_Z", "p.Moran_P", "p.Moran_Padj") %in% colnames(svg_df)))

# gene names are typically rownames
if (!is.null(rownames(svg_df))) {
  svg_dt <- data.table::as.data.table(svg_df, keep.rownames = "gene")
} else if ("gene" %in% colnames(svg_df)) {
  svg_dt <- data.table::as.data.table(svg_df)
} else {
  stop("[ERROR] Cannot find gene identifiers (no rownames and no 'gene' column).")
}

# ---- 10.2 Rank SVGs (primary: Padj, secondary: Moran_I) ----
setorder(svg_dt, p.Moran_Padj, -p.Moran_I)

# ---- 10.3 Add convenient flags ----
svg_dt[, is_SVG := (p.Moran_Padj < 0.05)]

# ---- 10.4 Quick report ----
cat("[INFO] SVG table rows: ", nrow(svg_dt), "\n", sep = "")
cat("[INFO] Significant SVGs (FDR < 0.05): ", svg_dt[is_SVG == TRUE, .N], "\n", sep = "")

# print top 20
print(svg_dt[1:min(20, .N)])

# ---- 10.5 Save full SVG table ----
svg_outfile <- file.path(outdir, "tables", paste0(sample_ST_ID, ".SpaCET.SVG.univariate.Moran.tsv"))
data.table::fwrite(svg_dt, svg_outfile, sep = "\t")
cat("[INFO] Saved SVG table: ", svg_outfile, "\n", sep = "")

# step.11 Spatially co-expression ligand-receptor interactions
# run for in-house ligand–receptor database, 20 mins
Sys.time()
SpaCET_obj <- SpaCET.SpatialCorrelation(
  SpaCET_obj, 
  mode="bivariate", 
  item=NULL, 
  W=W, 
  nPermutation=1000
)
Sys.time()

############################################################
## Step 11) Spatially co-expressed ligand–receptor pairs
##         (bivariate Moran's I) — summarize + export
############################################################

# ---- 11.1 Get bivariate results ----
co_expression <- SpaCET_obj@results$SpatialCorrelation$bivariate
stopifnot(is.data.frame(co_expression))

# ---- 11.2 Convert to data.table and keep pair IDs ----
co_expression <- data.table::as.data.table(co_expression, keep.rownames = "pair_id")

# split ligand and receptor (your format: LIGAND_RECEPTOR)
co_expression[, c("ligand", "receptor") := tstrsplit(pair_id, "_", fixed = TRUE)]

# ---- 11.3 Rank pairs ----
setorder(co_expression, p.Moran_Padj, -p.Moran_I)

# flag significant
co_expression[, is_sig := (p.Moran_Padj < 0.05)]

# preview
print(co_expression[1:20])

# ---- 11.4 Save full table ----
outfile_all <- file.path(outdir, "tables",
                         paste0(sample_ST_ID, ".SpaCET.co_expression.spatial_bivariate.tsv"))
data.table::fwrite(co_expression, outfile_all, sep = "\t")

############################################################
## Final Step) Save SpaCET object as R data
############################################################

# make sure rds folder exists
dir.create(file.path(outdir, "rds"), showWarnings = FALSE, recursive = TRUE)

saveRDS(
  SpaCET_obj,
  file = file.path(outdir, "rds", paste0(sample_ST_ID, ".SpaCET_obj.final.rds")),
  compress = "xz"   # strong compression, smaller file
)

cat("[INFO] SpaCET object saved as RDS.\n")