#!/usr/bin/env Rscript
# Replot final UMAPs using the same palette as cell_type pie charts.

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
source(file.path(project_dir, "CLLCMA_alpha", "scripts", "helpers_cell_type_palette.R"))

anno_dir <- file.path(project_dir, "CLLCMA_alpha", "integrated_annotated")
seu_file <- file.path(anno_dir, "seu_integrated_final.rds")
out_celltype <- file.path(anno_dir, "umap_cell_type_final.pdf")
out_dataset <- file.path(anno_dir, "umap_dataset_final.pdf")

ref_types <- fread(
  file.path(anno_dir, "cell_annotations_final.tsv"),
  select = "cell_type_final"
)$cell_type_final
pal <- cell_type_palette(ref_types)

if (!file.exists(seu_file)) {
  stop("Seurat object not found: ", seu_file, " — run steps 1–3 first.")
}

message("Loading Seurat object (large) ...")
seu <- readRDS(seu_file)

if (!"umap" %in% Reductions(seu)) {
  stop("No UMAP in object; re-run step 2.")
}

missing_types <- setdiff(unique(as.character(seu$cell_type_final)), names(pal))
if (length(missing_types)) {
  extra <- cell_type_palette(missing_types)
  pal <- c(pal, extra)
}

save_umap <- function(plot, path, width, height) {
  if (capabilities("cairo") && tools::file_ext(path) == "pdf") {
    ggsave(path, plot, width = width, height = height, device = cairo_pdf, bg = "white")
  } else {
    ggsave(path, plot, width = width, height = height, bg = "white")
  }
}

p_celltype <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "cell_type_final",
  cols = pal,
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  ggtitle("cell_type_final (patient B malignant: CLL / RT / MBL; Nadeu original)") +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 8, face = "bold"),
    legend.title = element_text(face = "bold")
  )

save_umap(p_celltype, out_celltype, width = 14, height = 9)
message("Wrote: ", out_celltype)

p_dataset <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "dataset",
  shuffle = TRUE,
  raster = TRUE
) +
  ggtitle("UMAP by dataset") +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 8, face = "bold"),
    legend.title = element_text(face = "bold")
  )

save_umap(p_dataset, out_dataset, width = 10, height = 7)
message("Wrote: ", out_dataset)
