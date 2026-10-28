#!/usr/bin/env Rscript

# Export a CLL single-cell reference for CIBERSORT / BayesPrism from CLLCMA_alpha.
# Uses fine AIFI-L2–like cell types (xCell-style resolution). Malignant cells are
# taken only from disease==CLL patients (no RT, no MBL). TME (T/NK/mono/DC/...)
# may still come from any GSE cohort. Excludes PBMCpedia_HC and Nadeu_2022.
#
# Inputs:  CLLCMA_alpha/integrated_annotated/seu_integrated_final.rds
# Outputs: CLLCMA_alpha/integrated_annotated/bayesprism/
#   sc_counts.rds
#   cell_type_labels.rds
#   cell_state_labels.rds
#   bayesprism_cell_table.tsv
#   bayesprism_label_summary.tsv
#   bayesprism_type_counts.tsv
#   bayesprism_type_map.tsv

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

project_dir <- Sys.getenv("CLL_TME_DIR", unset = path.expand("~/research/CLL_TME"))
setwd(normalizePath(project_dir, mustWork = TRUE))

seu_file <- "CLLCMA_alpha/integrated_annotated/seu_integrated_final.rds"
out_dir <- "CLLCMA_alpha/integrated_annotated/bayesprism"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

max_per_type <- 200L
min_per_type <- 20L
set.seed(42L)

gse_datasets <- c(
  "GSE111014",
  "GSE132065",
  "GSE152469",
  "GSE161610",
  "GSE205704",
  "GSE261917",
  "GSE295489",
  "GSE295491"
)

# Fine cell_type_final -> deconvolution reference label (xCell-like resolution).
# Names are CIBERSORT/BayesPrism-safe (no spaces). Platelets/erythrocytes excluded.
# Malignant compartment: CLL only (no RT, no MBL) — filtered by disease below.
type_map <- c(
  # malignant from CLL-disease patients only (see disease filter)
  "CLL" = "CLL",
  # CD4 T (xCell: CD4 naive/memory/Treg/Th-like)
  "Naive CD4 T cell" = "CD4_T_naive",
  "Memory CD4 T cell" = "CD4_T_memory",
  "Treg" = "Treg",
  "Proliferating T cell" = "T_proliferating",
  # CD8 / unconventional T (xCell: CD8 naive/memory, gdT, NKT-like)
  "Naive CD8 T cell" = "CD8_T_naive",
  "Memory CD8 T cell" = "CD8_T_memory",
  "CD8aa" = "CD8aa",
  "gdT" = "gdT",
  "DN T cell" = "DN_T",
  "MAIT" = "MAIT",
  # NK (xCell: NK cell)
  "CD56dim NK cell" = "NK_CD56dim",
  "CD56bright NK cell" = "NK_CD56bright",
  "Proliferating NK cell" = "NK_proliferating",
  # monocytes (xCell: Monocyte / macrophage-adjacent)
  "CD14 monocyte" = "Monocyte_CD14",
  "CD16 monocyte" = "Monocyte_CD16",
  "Intermediate monocyte" = "Monocyte_intermediate",
  # dendritic (xCell: mDC / pDC / activated)
  "cDC1" = "cDC1",
  "cDC2" = "cDC2",
  "ASDC" = "ASDC",
  "pDC" = "pDC",
  # innate / progenitor (xCell: ILC, CLP/CMP/HSC-like)
  "ILC" = "ILC",
  "Progenitor cell" = "Progenitor"
  # excluded malignant: RT, MBL (Richter disease==RS samples)
  # excluded other: Platelet, Erythrocyte, residual B subtypes
)

malignant_types <- c("CLL")
cll_disease_only <- c("CLL")

if (!file.exists(seu_file)) {
  stop("Missing: ", seu_file)
}

message("Loading ", seu_file)
seu <- readRDS(seu_file)
meta <- as.data.table(seu[[]], keep.rownames = "cell_id")

need <- c("dataset", "sample", "donor", "disease", "cell_type_final")
miss <- setdiff(need, names(meta))
if (length(miss)) {
  stop("Missing metadata columns: ", paste(miss, collapse = ", "))
}

meta[, fine_label := as.character(cell_type_final)]
meta[, cell_type := unname(type_map[fine_label])]

unmapped <- sort(unique(meta[is.na(cell_type), fine_label]))
if (length(unmapped)) {
  message(
    "Unmapped labels (excluded from reference): ",
    paste(unmapped, collapse = ", ")
  )
}

meta <- meta[!is.na(cell_type)]
meta <- meta[dataset %in% gse_datasets]

# Malignant cells: CLL disease patients only (drop RT / HC-MBL / LC-MBL / MBL)
n_malig_before <- meta[cell_type %in% malignant_types, .N]
meta <- meta[
  !(cell_type %in% malignant_types) |
    (cell_type == "CLL" & disease %in% cll_disease_only)
]
n_malig_after <- meta[cell_type == "CLL", .N]
message(
  "Malignant CLL cells kept (disease==CLL only): ",
  n_malig_after, " / ", n_malig_before,
  " (excluded RT/MBL malignant)"
)

if (!nrow(meta)) {
  stop("No GSE cells left after type mapping / dataset / disease filter")
}
if (!"CLL" %in% meta$cell_type) {
  stop("No CLL malignant cells left after disease==CLL filter")
}

meta[, patient := as.character(donor)]
meta[, cell_state := fifelse(
  cell_type %in% malignant_types,
  paste0(cell_type, "_", patient),
  cell_type
)]

message("Eligible cells by fine type (pre-subsample):")
print(meta[, .N, by = cell_type][order(-N)])
message("N fine types: ", uniqueN(meta$cell_type))

# Subsample up to max_per_type cells per fine cell_type
keep <- meta[
  ,
  {
    if (.N <= max_per_type) {
      .SD
    } else {
      .SD[sample(.N, max_per_type)]
    }
  },
  by = cell_type
]

type_n <- keep[, .N, by = cell_type]
drop_types <- type_n[N < min_per_type, cell_type]
if (length(drop_types)) {
  message("Dropping rare types (<", min_per_type, "): ", paste(drop_types, collapse = ", "))
  keep <- keep[!cell_type %in% drop_types]
}

setorder(keep, cell_type, dataset, sample, cell_id)

message("Reference cells after subsample (max ", max_per_type, " / type):")
print(keep[, .N, by = cell_type][order(-N)])
message("By dataset:")
print(keep[, .N, by = dataset][order(-N)])
message("Final n types: ", uniqueN(keep$cell_type))

counts <- LayerData(seu, assay = "RNA", layer = "counts")
counts <- counts[, keep$cell_id, drop = FALSE]
counts <- counts[Matrix::rowSums(counts) > 0, , drop = FALSE]
sc_dat <- t(counts) # cells x genes
sc_dat <- as(sc_dat, "CsparseMatrix")

cell_type_labels <- keep$cell_type
cell_state_labels <- keep$cell_state
names(cell_type_labels) <- keep$cell_id
names(cell_state_labels) <- keep$cell_id

stopifnot(
  identical(rownames(sc_dat), keep$cell_id),
  "CLL" %in% cell_type_labels,
  all(sc_dat@x >= 0)
)

saveRDS(sc_dat, file.path(out_dir, "sc_counts.rds"), compress = FALSE)
saveRDS(unname(cell_type_labels), file.path(out_dir, "cell_type_labels.rds"))
saveRDS(unname(cell_state_labels), file.path(out_dir, "cell_state_labels.rds"))

fwrite(
  keep[, .(
    cell_id, dataset, sample, patient, disease, tissue,
    fine_label, cell_type, cell_state
  )],
  file.path(out_dir, "bayesprism_cell_table.tsv"),
  sep = "\t"
)
fwrite(
  keep[, .(
    n_cells = .N,
    n_patients = uniqueN(patient),
    n_datasets = uniqueN(dataset)
  ), by = .(cell_type, cell_state)],
  file.path(out_dir, "bayesprism_label_summary.tsv"),
  sep = "\t"
)
fwrite(
  keep[, .(n_cells = .N), by = cell_type][order(-n_cells)],
  file.path(out_dir, "bayesprism_type_counts.tsv"),
  sep = "\t"
)
fwrite(
  data.table(
    fine_label = names(type_map),
    cell_type = unname(type_map)
  ),
  file.path(out_dir, "bayesprism_type_map.tsv"),
  sep = "\t"
)

message(
  "Wrote ", nrow(sc_dat), " cells x ", ncol(sc_dat),
  " genes, ", uniqueN(keep$cell_type), " fine types -> ",
  normalizePath(out_dir)
)
