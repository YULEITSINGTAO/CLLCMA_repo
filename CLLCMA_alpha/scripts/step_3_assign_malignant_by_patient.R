#!/usr/bin/env Rscript

# Assign patient B-lineage cells as malignant CLL/RT/MBL for unlabeled cohorts.
# Nadeu_2022 is excluded: cell_type_final stays annotation_original.
#
# Inputs:  CLLCMA_alpha/integrated_annotated/seu_integrated.rds
# Outputs: CLLCMA_alpha/integrated_annotated/seu_integrated_final.rds + tables/UMAPs

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
})

project_dir <- Sys.getenv(
  "CLL_TME_DIR",
  unset = path.expand("~/research/CLL_TME")
)
project_dir <- normalizePath(project_dir, mustWork = TRUE)
setwd(project_dir)

out_dir <- "CLLCMA_alpha/integrated_annotated"
seu_file <- file.path(out_dir, "seu_integrated.rds")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

b_as_malignant <- c(
  "Naive B cell",
  "Memory B cell",
  "Transitional B cell",
  "Effector B cell",
  "Plasma cell",
  "pDC"
)

skip_datasets <- c("PBMCpedia_HC", "Nadeu_2022")

# -------------------------------- Load ---------------------------------------

if (!file.exists(seu_file)) {
  stop("Missing: ", seu_file)
}

message("Loading ", seu_file)
seu <- readRDS(seu_file)

need <- c("dataset", "donor", "disease", "sample", "cell_type_final")
miss <- setdiff(need, colnames(seu[[]]))
if (length(miss)) {
  stop("Missing metadata: ", paste(miss, collapse = ", "))
}

seu$patient <- as.character(seu$donor)
seu$cell_type_final_pre_malignant <- as.character(seu$cell_type_final)

pre_label <- as.character(seu$cell_type_final_pre_malignant)

# --------------------- Assign malignant in CLL/RT cohorts --------------------

is_assignable <- !(seu$dataset %in% skip_datasets)
is_b <- pre_label %in% b_as_malignant
assign_idx <- which(is_assignable & is_b)

seu$is_malignant_cll <- FALSE
seu$malignant_call_method <- NA_character_

seu$cell_type_final <- pre_label
if (length(assign_idx)) {
  disease_at <- as.character(seu$disease[assign_idx])
  # Richter patients (disease RS) -> RT; MBL patients -> MBL; else -> CLL
  new_label <- ifelse(
    disease_at == "RS", "RT",
    ifelse(disease_at %in% c("HC-MBL", "LC-MBL", "MBL"), "MBL", "CLL")
  )
  seu$cell_type_final[assign_idx] <- new_label
  seu$is_malignant_cll[assign_idx] <- TRUE
  seu$malignant_call_method[assign_idx] <- "all_patient_B_as_CLL_RT_or_MBL"
}

# Nadeu: lock original annotations and cell_state
nadeu_idx <- which(seu$dataset == "Nadeu_2022")
if (length(nadeu_idx)) {
  nadeu_label <- as.character(seu$annotation_original[nadeu_idx])
  seu$cell_type_final[nadeu_idx] <- nadeu_label
  seu$is_malignant_cll[nadeu_idx] <- FALSE
  seu$malignant_call_method[nadeu_idx] <- NA_character_
}

seu$cell_state <- as.character(seu$cell_type_final)
malig_idx <- which(seu$cell_type_final %in% c("CLL", "RT", "MBL"))
if (length(malig_idx)) {
  seu$cell_state[malig_idx] <- paste0(
    seu$cell_type_final[malig_idx],
    "_",
    seu$patient[malig_idx]
  )
}
if (length(nadeu_idx)) {
  seu$cell_state[nadeu_idx] <- as.character(seu$annotation_original[nadeu_idx])
}

message(
  "Assigned ", length(assign_idx), " non-Nadeu patient B-lineage cells as CLL/RT (",
  sprintf("%.1f", 100 * length(assign_idx) / max(1, sum(is_assignable))),
  "% of assignable cells)"
)

# ------------------------------ Summaries ------------------------------------

meta <- as.data.table(seu[[]], keep.rownames = "cell_id")

summ_ds <- meta[
  ,
  .(
    n_cells = .N,
    n_CLL = sum(cell_type_final == "CLL"),
    n_MBL = sum(cell_type_final == "MBL"),
    n_RT = sum(cell_type_final == "RT"),
    pct_malignant = 100 * mean(is_malignant_cll)
  ),
  by = dataset
]

summ_sample <- meta[
  !(dataset %in% skip_datasets),
  .(
    n_cells = .N,
    n_CLL = sum(cell_type_final == "CLL"),
    n_MBL = sum(cell_type_final == "MBL"),
    n_RT = sum(cell_type_final == "RT"),
    pct_malignant = 100 * mean(is_malignant_cll),
    n_B_pre = sum(cell_type_final_pre_malignant %in% b_as_malignant)
  ),
  by = .(dataset, sample, patient, disease)
]

summ_type <- meta[
  ,
  .N,
  by = .(dataset, cell_type_final)
]
setorder(summ_type, dataset, -N)

fwrite(summ_sample, file.path(out_dir, "malignant_assignment_summary.tsv"), sep = "\t")
fwrite(summ_ds, file.path(out_dir, "malignant_assignment_by_dataset.tsv"), sep = "\t")
message("Per-dataset:")
print(summ_ds)

# ------------------------------ Save -----------------------------------------

anno <- data.frame(
  cell_id = colnames(seu),
  dataset = seu$dataset,
  sample = seu$sample,
  donor = seu$donor,
  patient = seu$patient,
  disease = seu$disease,
  tissue = seu$tissue,
  annotation_original = if ("annotation_original" %in% colnames(seu[[]])) {
    seu$annotation_original
  } else {
    NA_character_
  },
  cell_type_final_pre_malignant = seu$cell_type_final_pre_malignant,
  cell_type_final = seu$cell_type_final,
  cell_state = seu$cell_state,
  is_malignant_cll = seu$is_malignant_cll,
  malignant_call_method = seu$malignant_call_method,
  prediction_score = if ("prediction_score" %in% colnames(seu[[]])) {
    seu$prediction_score
  } else {
    NA_real_
  },
  stringsAsFactors = FALSE
)
fwrite(anno, file.path(out_dir, "cell_annotations_final.tsv"), sep = "\t")
saveRDS(seu, file.path(out_dir, "seu_integrated_final.rds"), compress = FALSE)

message("Wrote:\n  ", file.path(out_dir, "seu_integrated_final.rds"))
message("  ", file.path(out_dir, "cell_annotations_final.tsv"))
message("  ", file.path(out_dir, "malignant_assignment_summary.tsv"))

# ------------------------------ UMAPs ----------------------------------------

if (!"umap" %in% Reductions(seu)) {
  stop("No UMAP in object; re-run step 2.")
}

p1 <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "cell_type_final",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  ggtitle("UMAP cell_type_final (patient B → CLL/RT/MBL; Nadeu original)") +
  theme(legend.text = element_text(size = 6))
ggsave(
  file.path(out_dir, "umap_cell_type_final.pdf"),
  p1, width = 14, height = 9
)

p2 <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "dataset",
  shuffle = TRUE,
  raster = TRUE
) +
  ggtitle("UMAP by dataset")
ggsave(
  file.path(out_dir, "umap_dataset_final.pdf"),
  p2, width = 10, height = 7
)

p3 <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "is_malignant_cll",
  raster = TRUE
) +
  ggtitle("UMAP is_malignant_cll")
ggsave(
  file.path(out_dir, "umap_is_malignant_cll.pdf"),
  p3, width = 8, height = 6
)

p4 <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "disease",
  shuffle = TRUE,
  raster = TRUE
) +
  ggtitle("UMAP by disease")
ggsave(
  file.path(out_dir, "umap_disease_final.pdf"),
  p4, width = 8, height = 6
)

message("Wrote UMAPs to ", normalizePath(out_dir))
message("Done.")
