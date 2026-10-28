#!/usr/bin/env Rscript

# step_1_harmonize_datasets.R
#
# Purpose
# -------
# Build one atlas-wide sparse count matrix + unified cell metadata from all
# cohorts staged under CLLCMA_alpha/datasets/. Downstream steps assume this
# shared schema (same genes, unique cell_ids, common metadata columns).
#
# Pipeline position
# -----------------
#   step_1  (this script)  -> CLLCMA_alpha/integrated_raw/
#   step_1.5               -> per-sample QC on integrated_raw
#   step_2                 -> Harmony + AIFI label transfer
#   step_3+                -> malignant labels, BayesPrism ref, trajectory
#
# Inputs:  CLLCMA_alpha/datasets/<cohort>/*_counts.rds + *_cell_metadata.tsv
# Outputs: CLLCMA_alpha/integrated_raw/
#   integrated_counts.rds           genes x cells (CsparseMatrix)
#   integrated_cell_metadata.tsv    one row per cell
#   integrated_gene_metadata.tsv    shared gene symbols
#   integrated_summary.tsv          n_cells by dataset/sample/disease/tissue
#
# Design notes
# ------------
# - Cohorts use different metadata column names; each entry in `cohorts`
#   maps those names onto a common schema.
# - Only genes present in EVERY cohort are kept (intersection). This drops
#   cohort-specific genes but makes matrices cbind-compatible.
# - Cell IDs are prefixed with dataset ("GSE295491_...") so barcodes never
#   collide across studies.
# - This step does NOT run QC, normalize, or integrate — only harmonize IDs
#   and stack matrices. Use step_1.5 for QC.

suppressPackageStartupMessages({
  library(Matrix)
})

# Project root: allow override via CLL_TME_DIR for non-default installs.
project_dir <- Sys.getenv(
  "CLL_TME_DIR",
  unset = path.expand("~/research/CLL_TME")
)
project_dir <- normalizePath(project_dir, mustWork = TRUE)
setwd(project_dir)

atlas_dir <- "CLLCMA_alpha"
data_dir <- file.path(atlas_dir, "datasets")
out_dir <- file.path(atlas_dir, "integrated_raw")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------- Cohort definitions ------------------------------
#
# Each list() is one cohort. Fields:
#   dataset          Atlas name (also used as cell_id prefix)
#   counts / meta    Paths to genes x cells RDS and per-cell TSV
#   id_col           Metadata column matching colnames(counts)
#   sample_col       Library / sample label (batch unit for Harmony later)
#   donor_col        Patient / donor ID
#   disease          Fixed disease label for all cells in this cohort, OR
#   disease_col      Metadata column with per-cell disease (overrides `disease`)
#   tissue           Fixed tissue (e.g. "PB"), OR
#   tissue_col       Metadata column with per-cell tissue
#   cell_type_col    Optional prior annotation (used as transfer reference for HC)
#   cell_type_L1_col Optional coarse prior annotation
#   annotation_col   Optional "locked" labels (Nadeu keeps these in step_2/3)
#   gsm_col          Optional GEO sample accession
#
# NULL for a *_col field means "not available in this cohort's metadata".

cohorts <- list(
  # Healthy PBMC reference with AIFI L1/L2 labels (used for transfer in step_2).
  list(
    dataset = "PBMCpedia_HC",
    counts = file.path(data_dir, "PBMCpedia_HC", "PBMCpedia_HC_raw_subsample_counts.rds"),
    meta = file.path(data_dir, "PBMCpedia_HC", "PBMCpedia_HC_cell_metadata.tsv"),
    id_col = "barcode_sample",
    sample_col = "sample_id",
    donor_col = "donor_id",
    disease = "HC",
    tissue = "PB",
    cell_type_col = "AIFI_L2_name",
    cell_type_L1_col = "AIFI_L1_name",
    gsm_col = NULL,
    annotation_col = NULL,
    disease_col = NULL
  ),
  # Early CLL PB day-0 libraries.
  list(
    dataset = "GSE111014",
    counts = file.path(data_dir, "GSE111014_day0", "GSE111014_day0_counts.rds"),
    meta = file.path(data_dir, "GSE111014_day0", "GSE111014_day0_cell_metadata.tsv"),
    id_col = "cell_id",
    sample_col = "sample_label",
    donor_col = "patient",
    disease = "CLL",
    tissue = "PB",
    cell_type_col = NULL,
    cell_type_L1_col = NULL,
    gsm_col = "gsm",
    annotation_col = NULL,
    disease_col = NULL
  ),
  list(
    dataset = "GSE132065",
    counts = file.path(data_dir, "GSE132065", "GSE132065_CLL_counts.rds"),
    meta = file.path(data_dir, "GSE132065", "GSE132065_CLL_cell_metadata.tsv"),
    id_col = "cell_id",
    sample_col = "sample_label",
    donor_col = "patient",
    disease = "CLL",
    tissue = NULL,
    tissue_col = "tissue",
    cell_type_col = NULL,
    cell_type_L1_col = NULL,
    gsm_col = "gsm",
    annotation_col = NULL,
    disease_col = NULL
  ),
  list(
    dataset = "GSE152469",
    counts = file.path(data_dir, "GSE152469", "GSE152469_PB_counts.rds"),
    meta = file.path(data_dir, "GSE152469", "GSE152469_PB_cell_metadata.tsv"),
    id_col = "cell_id",
    sample_col = "sample_label",
    donor_col = "patient",
    disease = "CLL",
    tissue = NULL,
    tissue_col = "tissue",
    cell_type_col = NULL,
    cell_type_L1_col = NULL,
    gsm_col = "gsm",
    annotation_col = NULL,
    disease_col = NULL
  ),
  list(
    dataset = "GSE161610",
    counts = file.path(data_dir, "GSE161610", "GSE161610_PB_counts.rds"),
    meta = file.path(data_dir, "GSE161610", "GSE161610_PB_cell_metadata.tsv"),
    id_col = "cell_id",
    sample_col = "sample_label",
    donor_col = "patient",
    disease = "CLL",
    tissue = NULL,
    tissue_col = "tissue",
    cell_type_col = NULL,
    cell_type_L1_col = NULL,
    gsm_col = "gsm",
    annotation_col = NULL,
    disease_col = NULL
  ),
  # Includes CLL and Richter (RS) samples; disease taken from metadata.
  list(
    dataset = "GSE205704",
    counts = file.path(data_dir, "GSE205704", "GSE205704_counts.rds"),
    meta = file.path(data_dir, "GSE205704", "GSE205704_cell_metadata.tsv"),
    id_col = "cell_id",
    sample_col = "sample_label",
    donor_col = "patient",
    disease = NULL,
    disease_col = "disease",
    tissue = NULL,
    tissue_col = "tissue",
    cell_type_col = NULL,
    cell_type_L1_col = NULL,
    gsm_col = "gsm",
    annotation_col = NULL
  ),
  # Lymph-node CLL (tissue from metadata, usually LN).
  list(
    dataset = "GSE261917",
    counts = file.path(data_dir, "GSE261917", "GSE261917_LN_counts.rds"),
    meta = file.path(data_dir, "GSE261917", "GSE261917_LN_cell_metadata.tsv"),
    id_col = "cell_id",
    sample_col = "sample_label",
    donor_col = "patient",
    disease = "CLL",
    tissue = NULL,
    tissue_col = "tissue",
    cell_type_col = NULL,
    cell_type_L1_col = NULL,
    gsm_col = "gsm",
    annotation_col = NULL,
    disease_col = NULL
  ),
  list(
    dataset = "GSE295489",
    counts = file.path(data_dir, "GSE295489", "GSE295489_PB_counts.rds"),
    meta = file.path(data_dir, "GSE295489", "GSE295489_PB_cell_metadata.tsv"),
    id_col = "cell_id",
    sample_col = "sample_label",
    donor_col = "patient",
    disease = "CLL",
    tissue = NULL,
    tissue_col = "tissue",
    cell_type_col = NULL,
    cell_type_L1_col = NULL,
    gsm_col = "gsm",
    annotation_col = NULL,
    disease_col = NULL
  ),
  # Paired MBL / CLL RNA (disease_stage: CLL, HC-MBL, LC-MBL).
  list(
    dataset = "GSE295491",
    counts = file.path(data_dir, "GSE295491", "GSE295491_PB_counts.rds"),
    meta = file.path(data_dir, "GSE295491", "GSE295491_PB_cell_metadata.tsv"),
    id_col = "cell_id",
    sample_col = "sample_label",
    donor_col = "patient",
    disease = NULL,
    disease_col = "disease_stage",
    tissue = NULL,
    tissue_col = "tissue",
    cell_type_col = NULL,
    cell_type_L1_col = NULL,
    gsm_col = "gsm",
    annotation_col = NULL
  ),
  # Author annotations kept as annotation_final (not overwritten by transfer).
  list(
    dataset = "Nadeu_2022",
    counts = file.path(data_dir, "Nadeu_2022", "Nadeu2022_counts.rds"),
    meta = file.path(data_dir, "Nadeu_2022", "Nadeu2022_cell_metadata.tsv"),
    id_col = "cell_id",
    sample_col = "sample_label",
    donor_col = "patient",
    disease = NULL,
    disease_col = "disease_stage",
    tissue = NULL,
    tissue_col = "tissue",
    cell_type_col = "annotation_final",
    cell_type_L1_col = NULL,
    gsm_col = NULL,
    annotation_col = "annotation_final"
  )
)

# Map alternate disease spellings to the atlas standard ("RS" = Richter).
normalize_disease <- function(x) {
  x <- as.character(x)
  x[x %in% c("Richter_transformation", "Richter", "richter")] <- "RS"
  x
}

# Map alternate tissue spellings to atlas standard ("PB" = peripheral blood).
normalize_tissue <- function(x) {
  x <- as.character(x)
  x[x %in% c("PBMC", "pbmc")] <- "PB"
  x
}

# Reorder / subset metadata rows so they exactly match colnames(counts).
# Stops if any count column lacks a metadata row (would break Seurat later).
align_meta_to_cells <- function(counts, meta, id_col, dataset_label) {
  # PBMCpedia export can duplicate barcode_sample; keep first unique match.
  if (anyDuplicated(meta[[id_col]])) {
    meta <- meta[!duplicated(meta[[id_col]]), , drop = FALSE]
  }
  if (!id_col %in% names(meta)) {
    stop(dataset_label, ": missing id column ", id_col)
  }
  idx <- match(colnames(counts), meta[[id_col]])
  if (anyNA(idx)) {
    stop(
      dataset_label, ": ", sum(is.na(idx)),
      " cells in counts missing from metadata"
    )
  }
  meta[idx, , drop = FALSE]
}

# Safely pull a metadata column; return NA vector if the column is absent/NULL.
get_col_or <- function(df, col, default = NA_character_) {
  if (is.null(col) || !nzchar(col) || !col %in% names(df)) {
    return(rep(default, nrow(df)))
  }
  as.character(df[[col]])
}

# ------------------------------- Load ----------------------------------------
# Read each cohort, standardize metadata, and store matrices for later cbind.

message("Loading count matrices ...")
counts_list <- list()
meta_list <- list()
gene_lists <- list()

for (cfg in cohorts) {
  ds <- cfg$dataset
  message("  ", ds, " ...")
  if (!file.exists(cfg$counts)) {
    stop("Missing counts: ", cfg$counts)
  }
  if (!file.exists(cfg$meta)) {
    stop("Missing metadata: ", cfg$meta)
  }

  # genes x cells sparse matrix; gene symbols must be unique within cohort.
  counts <- readRDS(cfg$counts)
  counts <- as(counts, "CsparseMatrix")
  if (anyDuplicated(rownames(counts))) {
    stop(ds, ": gene names are not unique")
  }

  raw_meta <- read.delim(
    cfg$meta,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  # Some PBMCpedia exports put the ID in the first unnamed column.
  if (ds == "PBMCpedia_HC" && !"barcode_sample" %in% names(raw_meta)) {
    names(raw_meta)[1] <- "barcode_sample"
  }
  raw_meta <- align_meta_to_cells(counts, raw_meta, cfg$id_col, ds)

  # Tissue: prefer a per-cell column; else use the fixed cohort value.
  if (!is.null(cfg[["tissue_col"]])) {
    tissue_vals <- get_col_or(raw_meta, cfg$tissue_col)
  } else if (!is.null(cfg$tissue)) {
    tissue_vals <- rep(as.character(cfg$tissue), nrow(raw_meta))
  } else {
    tissue_vals <- rep(NA_character_, nrow(raw_meta))
  }
  tissue_vals <- normalize_tissue(tissue_vals)

  # Disease: same pattern — column (e.g. MBL/CLL/RS) or fixed label (e.g. HC).
  if (!is.null(cfg$disease_col)) {
    disease_vals <- normalize_disease(get_col_or(raw_meta, cfg$disease_col))
  } else {
    disease_vals <- rep(as.character(cfg$disease), nrow(raw_meta))
  }

  # Optional annotations: filled for HC / Nadeu; NA for unlabeled GSE cohorts
  # (those get labels later via AIFI transfer in step_2).
  cell_type <- get_col_or(raw_meta, cfg$cell_type_col)
  cell_type_L1 <- get_col_or(raw_meta, cfg$cell_type_L1_col)
  annotation_original <- get_col_or(raw_meta, cfg$annotation_col)
  gsm_vals <- get_col_or(raw_meta, cfg$gsm_col)

  # Common metadata schema written for every cell in the atlas.
  meta <- data.frame(
    cell_id_raw = as.character(raw_meta[[cfg$id_col]]),
    dataset = ds,
    sample = as.character(raw_meta[[cfg$sample_col]]),
    donor = as.character(raw_meta[[cfg$donor_col]]),
    disease = disease_vals,
    tissue = tissue_vals,
    cell_type = cell_type,
    cell_type_L1 = cell_type_L1,
    annotation_original = annotation_original,
    gsm = gsm_vals,
    stringsAsFactors = FALSE
  )

  # Globally unique cell IDs: "<dataset>_<original_barcode>".
  new_ids <- paste0(ds, "_", colnames(counts))
  colnames(counts) <- new_ids
  meta$cell_id <- new_ids

  message("    ", ncol(counts), " cells x ", nrow(counts), " genes")
  counts_list[[ds]] <- counts
  meta_list[[ds]] <- meta
  gene_lists[[ds]] <- rownames(counts)
  rm(counts, raw_meta, meta)
  gc(verbose = FALSE)
}

# -------------------- Gene intersection + matrix merge -----------------------
# Keep only genes shared by all cohorts so cbind does not introduce NAs.

shared_genes <- Reduce(intersect, gene_lists)
shared_genes <- sort(shared_genes)
message("Shared genes: ", length(shared_genes))
if (length(shared_genes) < 8000L) {
  # Guardrail: a very small intersection usually means symbol/Ensembl mismatch.
  stop("Shared gene count looks too low; check symbol mapping.")
}

message("Subsetting to shared genes and combining ...")
for (ds in names(counts_list)) {
  counts_list[[ds]] <- counts_list[[ds]][shared_genes, , drop = FALSE]
}

# Column-bind all cohorts into one genes x cells matrix.
counts <- do.call(cbind, counts_list)
cell_meta <- do.call(rbind, meta_list)
rownames(cell_meta) <- NULL

# Sanity checks before writing (hard failures are better than silent misalignment).
stopifnot(identical(colnames(counts), cell_meta$cell_id))
stopifnot(!any(duplicated(cell_meta$cell_id)))
stopifnot(identical(rownames(counts), shared_genes))

# Stable column order for downstream scripts.
cell_meta <- cell_meta[, c(
  "cell_id", "cell_id_raw", "dataset", "sample", "donor",
  "disease", "tissue", "cell_type", "cell_type_L1",
  "annotation_original", "gsm"
)]

gene_meta <- data.frame(
  feature_name = shared_genes,
  stringsAsFactors = FALSE
)

# Compact sample-level census (useful for QC plots and manuscript tables).
summary_tab <- as.data.frame(
  table(
    dataset = cell_meta$dataset,
    sample = cell_meta$sample,
    disease = cell_meta$disease,
    tissue = cell_meta$tissue
  ),
  stringsAsFactors = FALSE
)
names(summary_tab)[names(summary_tab) == "Freq"] <- "n_cells"
summary_tab <- summary_tab[summary_tab$n_cells > 0, ]
summary_tab <- summary_tab[
  order(summary_tab$dataset, summary_tab$sample),
  ,
  drop = FALSE
]
rownames(summary_tab) <- NULL

# ------------------------------- Save ----------------------------------------

message("Writing outputs to ", out_dir, " ...")
# compress=FALSE: faster write/read for large sparse RDS on local disk.
saveRDS(counts, file.path(out_dir, "integrated_counts.rds"), compress = FALSE)
write.table(
  cell_meta,
  file.path(out_dir, "integrated_cell_metadata.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  gene_meta,
  file.path(out_dir, "integrated_gene_metadata.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  summary_tab,
  file.path(out_dir, "integrated_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

message(
  "Done: ", ncol(counts), " cells x ", nrow(counts),
  " genes -> ", normalizePath(out_dir)
)
print(aggregate(n_cells ~ dataset + disease, data = summary_tab, FUN = sum))
