#!/usr/bin/env Rscript

# Integrate all CLLCMA_alpha cohorts with Harmony, then transfer PBMCpedia
# AIFI_L2 labels to unlabeled CLL/RS/MBL cohorts. Nadeu_2022 keeps original
# annotation_final as cell_type_final (no transfer).
#
# Prefer running step_1.5_qc_per_sample.R after step_1 and before this script.
#
# Memory-safe transfer: after Harmony, predictions are transferred one query
# dataset at a time (avoids OOM on FindTransferAnchors with ~500k cells).
# Harmony object is checkpointed so transfer can resume after a crash.
#
# Inputs:  CLLCMA_alpha/integrated_raw/
# Outputs: CLLCMA_alpha/integrated_annotated/

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(harmony)
})

options(future.globals.maxSize = 50 * 1024^3)

project_dir <- Sys.getenv(
  "CLL_TME_DIR",
  unset = path.expand("~/research/CLL_TME")
)
project_dir <- normalizePath(project_dir, mustWork = TRUE)
setwd(project_dir)

in_dir <- "CLLCMA_alpha/integrated_raw"
out_dir <- "CLLCMA_alpha/integrated_annotated"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_features <- 3000L
n_pcs <- 30L
batch_key <- "sample"
ref_dataset <- "PBMCpedia_HC"
skip_transfer_datasets <- c("PBMCpedia_HC", "Nadeu_2022")

checkpoint_file <- file.path(out_dir, "seu_harmony_checkpoint.rds")
force_reharmony <- identical(toupper(Sys.getenv("CLLCMA_FORCE_REHARMONY", "")), "TRUE")

# ------------------------------- Helpers -------------------------------------

diet_for_transfer <- function(obj) {
  # Drop dense scale.data / unused layers before subsetting huge objects.
  if (inherits(obj[["RNA"]], "Assay5")) {
    tryCatch(
      {
        layers <- Layers(obj[["RNA"]])
        drop_layers <- grep("^scale\\.data", layers, value = TRUE)
        for (ly in drop_layers) {
          try(obj[["RNA"]][[ly]] <- NULL, silent = TRUE)
        }
      },
      error = function(e) invisible(NULL)
    )
  }
  tryCatch(
    DietSeurat(
      obj,
      assays = "RNA",
      dimreducs = intersect(c("pca", "harmony", "umap"), names(obj@reductions)),
      graphs = names(obj@graphs),
      counts = TRUE,
      data = TRUE,
      scale.data = FALSE
    ),
    error = function(e) obj
  )
}

transfer_one_query <- function(ref, query_cells, seu, n_features, n_pcs,
                               max_query_cells = 80000L) {
  # Split huge datasets into chunks to keep FindTransferAnchors memory bounded.
  if (length(query_cells) > max_query_cells) {
    set.seed(1L)
    query_cells <- sample(query_cells)
    starts <- seq(1L, length(query_cells), by = max_query_cells)
    parts <- lapply(starts, function(i) {
      j <- min(i + max_query_cells - 1L, length(query_cells))
      message("  chunk cells ", i, "-", j, " / ", length(query_cells))
      transfer_one_query(
        ref, query_cells[i:j], seu, n_features, n_pcs,
        max_query_cells = max_query_cells
      )
    })
    return(do.call(rbind, parts))
  }

  message("  query cells: ", length(query_cells))
  query <- subset(seu, cells = query_cells)

  query <- NormalizeData(query, verbose = FALSE)
  query <- FindVariableFeatures(query, nfeatures = n_features, verbose = FALSE)
  # Scale only shared variable features to limit dense matrix size
  feats <- intersect(VariableFeatures(ref), rownames(query))
  if (length(feats) < 500L) {
    feats <- VariableFeatures(ref)
  }
  query <- ScaleData(query, features = feats, verbose = FALSE)
  query <- RunPCA(query, features = feats, npcs = 50L, verbose = FALSE)

  anchors <- FindTransferAnchors(
    reference = ref,
    query = query,
    reference.reduction = "pca",
    dims = 1:n_pcs,
    features = feats,
    normalization.method = "LogNormalize",
    k.filter = NA,
    verbose = TRUE
  )

  pred <- TransferData(
    anchorset = anchors,
    refdata = ref$cell_type_ref,
    dims = 1:n_pcs,
    verbose = FALSE
  )

  pred_l1 <- TransferData(
    anchorset = anchors,
    refdata = ref$cell_type_L1,
    dims = 1:n_pcs,
    verbose = FALSE
  )

  out <- data.frame(
    cell_id = colnames(query),
    predicted_cell_type = as.character(pred$predicted.id),
    prediction_score = as.numeric(pred$prediction.score.max),
    predicted_L1 = as.character(pred_l1$predicted.id),
    prediction_score_L1 = as.numeric(pred_l1$prediction.score.max),
    stringsAsFactors = FALSE
  )
  
  rm(query, anchors, pred, pred_l1)
  gc(verbose = FALSE)
  out
}

# ------------------------------- Load / Harmony ------------------------------

if (file.exists(checkpoint_file) && !force_reharmony) {
  message("Loading Harmony checkpoint: ", checkpoint_file)
  seu <- readRDS(checkpoint_file)
} else {
  message("Loading harmonized counts + metadata ...")
  counts <- readRDS(file.path(in_dir, "integrated_counts.rds"))
  meta <- read.delim(
    file.path(in_dir, "integrated_cell_metadata.tsv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  stopifnot(identical(colnames(counts), meta$cell_id))
  rownames(meta) <- meta$cell_id

  seu <- CreateSeuratObject(
    counts = counts,
    meta.data = meta,
    project = "CLLCMA_alpha",
    min.cells = 0,
    min.features = 0
  )
  rm(counts, meta)
  gc(verbose = FALSE)

  seu$cell_type_ref <- ifelse(
    seu$dataset == ref_dataset,
    as.character(seu$cell_type),
    NA_character_
  )

  message(
    "Seurat object: ", ncol(seu), " cells x ", nrow(seu), " genes; datasets:\n",
    paste(capture.output(print(table(seu$dataset))), collapse = "\n")
  )

  message("NormalizeData / FindVariableFeatures / ScaleData / RunPCA ...")
  seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e4)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = n_features)
  seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)
  seu <- RunPCA(seu, npcs = 50L, verbose = FALSE)

  message("Harmony integration on batch key: ", batch_key, " ...")
  DefaultAssay(seu) <- "RNA"
  if (inherits(seu[["RNA"]], "Assay5")) {
    seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
  }

  seu <- RunHarmony(
    seu,
    group.by.vars = batch_key,
    reduction.use = "pca",
    dims.use = 1:50,
    reduction.save = "harmony",
    project.dim = FALSE
  )

  seu <- RunUMAP(seu, reduction = "harmony", dims = 1:n_pcs, verbose = FALSE)
  seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:n_pcs, verbose = FALSE)
  seu <- FindClusters(seu, resolution = 0.6, verbose = FALSE)

  message("Dieting object + writing Harmony checkpoint ...")
  seu <- diet_for_transfer(seu)
  saveRDS(seu, checkpoint_file, compress = FALSE)
  message("Checkpoint saved: ", checkpoint_file)
}

# -------------------- Label transfer (per dataset) ---------------------------

all_datasets <- unique(as.character(seu$dataset))
query_datasets <- setdiff(all_datasets, skip_transfer_datasets)
message(
  "Transfer query datasets (", length(query_datasets), "): ",
  paste(query_datasets, collapse = ", ")
)

ref_cells <- colnames(seu)[seu$dataset == ref_dataset]
nadeu_cells <- colnames(seu)[seu$dataset == "Nadeu_2022"]

message("Preparing reference (", length(ref_cells), " HC cells) ...")
ref <- subset(seu, cells = ref_cells)
ref <- NormalizeData(ref, verbose = FALSE)
ref <- FindVariableFeatures(ref, nfeatures = n_features, verbose = FALSE)
ref <- ScaleData(ref, features = VariableFeatures(ref), verbose = FALSE)
ref <- RunPCA(ref, npcs = 50L, verbose = FALSE)

seu$predicted_cell_type <- NA_character_
seu$prediction_score <- NA_real_
seu$predicted_L1 <- NA_character_
seu$prediction_score_L1 <- NA_real_

# HC truth
seu$predicted_cell_type[ref_cells] <- as.character(seu$cell_type_ref[ref_cells])
seu$predicted_L1[ref_cells] <- as.character(seu$cell_type_L1[ref_cells])

pred_rows <- list()
for (ds in query_datasets) {
  message("FindTransferAnchors for dataset: ", ds, " ...")
  q_cells <- colnames(seu)[seu$dataset == ds]
  if (!length(q_cells)) {
    next
  }
  pred_rows[[ds]] <- transfer_one_query(ref, q_cells, seu, n_features, n_pcs)
}

rm(ref)
gc(verbose = FALSE)

if (length(pred_rows)) {
  pred_all <- do.call(rbind, pred_rows)
  rownames(pred_all) <- NULL
  stopifnot(!anyDuplicated(pred_all$cell_id))
  idx <- match(pred_all$cell_id, colnames(seu))
  stopifnot(!anyNA(idx))
  seu$predicted_cell_type[idx] <- pred_all$predicted_cell_type
  seu$prediction_score[idx] <- pred_all$prediction_score
  seu$predicted_L1[idx] <- pred_all$predicted_L1
  seu$prediction_score_L1[idx] <- pred_all$prediction_score_L1
}

# Nadeu keeps original annotations
if (length(nadeu_cells)) {
  seu$predicted_cell_type[nadeu_cells] <- as.character(
    seu$annotation_original[nadeu_cells]
  )
}

seu$cell_type_final <- seu$predicted_cell_type

# -------------------- Save tables / object -----------------------------------

anno <- data.frame(
  cell_id = colnames(seu),
  dataset = seu$dataset,
  sample = seu$sample,
  donor = seu$donor,
  disease = seu$disease,
  tissue = seu$tissue,
  cell_type_ref = seu$cell_type_ref,
  annotation_original = seu$annotation_original,
  predicted_cell_type = seu$predicted_cell_type,
  prediction_score = seu$prediction_score,
  predicted_L1 = seu$predicted_L1,
  prediction_score_L1 = seu$prediction_score_L1,
  cell_type_final = seu$cell_type_final,
  seurat_clusters = as.character(seu$seurat_clusters),
  stringsAsFactors = FALSE
)

write.table(
  anno,
  file.path(out_dir, "cell_annotations.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

summary_tab <- as.data.frame(table(
  dataset = anno$dataset,
  cell_type_final = anno$cell_type_final
))
names(summary_tab)[3] <- "n_cells"
summary_tab <- summary_tab[summary_tab$n_cells > 0, ]
summary_tab <- summary_tab[order(summary_tab$dataset, -summary_tab$n_cells), ]
write.table(
  summary_tab,
  file.path(out_dir, "transfer_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

message("Saving Seurat object (large) ...")
saveRDS(seu, file.path(out_dir, "seu_integrated.rds"), compress = FALSE)

message("Writing UMAP PDFs ...")
p1 <- DimPlot(seu, reduction = "umap", group.by = "dataset", shuffle = TRUE) +
  ggtitle("Harmony UMAP by dataset")
ggsave(file.path(out_dir, "umap_dataset.pdf"), p1, width = 10, height = 7)

p2 <- DimPlot(
  seu, reduction = "umap", group.by = "cell_type_final",
  label = TRUE, repel = TRUE, raster = TRUE
) +
  ggtitle("cell_type_final (HC truth + transferred + Nadeu original)") +
  theme(legend.text = element_text(size = 6))
ggsave(file.path(out_dir, "umap_predicted_cell_type.pdf"), p2, width = 14, height = 9)

p3 <- FeaturePlot(seu, features = "prediction_score", reduction = "umap", raster = TRUE) +
  ggtitle("Transfer prediction score (NA on HC and Nadeu)")
ggsave(file.path(out_dir, "umap_prediction_score.pdf"), p3, width = 8, height = 6)

message("Done -> ", normalizePath(out_dir))
