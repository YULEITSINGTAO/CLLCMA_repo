#!/usr/bin/env Rscript

# Marker-gene heatmap for Normal B, MBL, CLL, RT from saved pseudobulk / DESeq2.
#
# A marker of a stage is:
#   - upregulated vs the other three (one-vs-rest padj < 0.05, LFC >= 1)
#   - highest mean expression in that stage
#   - not MT / ribosomal / IG-variable / hemoglobin
#
# Inputs:  CLLCMA_alpha/results/bcell_stage_deseq2/
# Outputs: same folder
#   stage_marker_genes.tsv
#   heatmap_stage_markers.pdf          (samples, columns ordered by stage)
#   heatmap_stage_markers_mean.pdf     (four stage means)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pheatmap)
})

project_dir <- Sys.getenv(
  "CLL_TME_DIR",
  unset = path.expand("~/research/CLL_TME")
)
setwd(normalizePath(project_dir, mustWork = TRUE))

out_dir <- "CLLCMA_alpha/results/bcell_stage_deseq2"
n_per_stage <- 15L
padj_cut <- 0.05
lfc_cut <- 1
max_show_per_stage <- 20L

stage_levels <- c("Normal_B", "MBL", "CLL", "RT")
stage_colors <- c(
  Normal_B = "#4C78A8",
  MBL = "#54A24B",
  CLL = "#F58518",
  RT = "#E45756"
)
ovr_files <- c(
  Normal_B = "deseq2_one_vs_rest_HC_B.tsv",
  MBL = "deseq2_one_vs_rest_MBL.tsv",
  CLL = "deseq2_one_vs_rest_CLL.tsv",
  RT = "deseq2_one_vs_rest_RT.tsv"
)

is_nonmarker <- function(g) {
  grepl(
    paste(
      "^MT-",
      "^RP[SL][0-9]",
      "^HB[AB]",
      "^IGH[VDJ]",
      "^IGK[VDJ]",
      "^IGL[VDJ]",
      "^HIST",
      "^AC[0-9]",
      "^AL[0-9]",
      "^LINC",
      sep = "|"
    ),
    g
  )
}

pb_mat <- readRDS(file.path(out_dir, "pseudobulk_counts.rds"))
pb_meta <- fread(file.path(out_dir, "pseudobulk_colData.tsv"))
stopifnot(identical(colnames(pb_mat), pb_meta$pb_id))

pb_meta[, stage := fifelse(b_stage == "HC_B", "Normal_B", as.character(b_stage))]
keep <- pb_meta$stage %in% stage_levels &
  !(pb_meta$b_stage == "RT" & pb_meta$disease != "RS")
pb_mat <- pb_mat[, keep, drop = FALSE]
pb_meta <- pb_meta[keep]
pb_meta[, stage := factor(stage, levels = stage_levels)]

lib <- pmax(colSums(pb_mat), 1)
logcpm <- log1p(t(t(pb_mat) / lib * 1e6))

mean_by_stage <- sapply(stage_levels, function(st) {
  rowMeans(logcpm[, pb_meta$stage == st, drop = FALSE])
})
best <- stage_levels[max.col(mean_by_stage, ties.method = "first")]
names(best) <- rownames(mean_by_stage)
second <- apply(mean_by_stage, 1L, function(x) sort(x, decreasing = TRUE)[2L])
delta <- mean_by_stage[cbind(seq_len(nrow(mean_by_stage)), match(best, stage_levels))] - second

marker_list <- list()
for (st in stage_levels) {
  ovr <- fread(file.path(out_dir, ovr_files[[st]]))
  up <- ovr[
    is.finite(padj) & padj < padj_cut & log2FoldChange >= lfc_cut &
      !is_nonmarker(gene)
  ]
  up[, highest_in := best[gene]]
  up[, delta := delta[gene]]
  up[, mean_logcpm := mean_by_stage[gene, st]]
  up <- up[highest_in == st]
  setorder(up, -delta, -log2FoldChange)
  up <- head(up, n_per_stage)
  up[, marker_of := st]
  up[, rank_in_stage := seq_len(.N)]
  marker_list[[st]] <- up
}

markers <- rbindlist(marker_list)
markers <- markers[
  ,
  .(
    marker_of, rank_in_stage, gene, baseMean, log2FoldChange, padj,
    delta, mean_logcpm
  )
]
fwrite(markers, file.path(out_dir, "stage_marker_genes.tsv"), sep = "\t")
message("Markers per stage:")
print(markers[, .N, by = marker_of])
print(markers[, .(marker_of, rank_in_stage, gene, log2FoldChange, padj)])

genes <- unique(markers$gene)
genes <- genes[genes %in% rownames(logcpm)]

# Mean heatmap (4 columns) — clearest marker view
mean_mat <- mean_by_stage[genes, stage_levels, drop = FALSE]
row_ann <- data.frame(
  marker_of = factor(markers$marker_of[match(genes, markers$gene)], levels = stage_levels),
  row.names = genes
)
ann_col_stage <- data.frame(
  stage = factor(stage_levels, levels = stage_levels),
  row.names = stage_levels
)
ann_colors <- list(
  stage = stage_colors,
  marker_of = stage_colors
)

pdf(file.path(out_dir, "heatmap_stage_markers_mean.pdf"), width = 6.5, height = 11)
pheatmap(
  mean_mat,
  scale = "row",
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  annotation_row = row_ann,
  annotation_col = ann_col_stage,
  annotation_colors = ann_colors,
  fontsize_row = 8,
  fontsize_col = 11,
  border_color = NA,
  main = "Stage marker genes (mean log1p CPM, row z-score)"
)
dev.off()

# Sample heatmap: up to max_show_per_stage donors per stage, ordered
set.seed(42L)
show_ids <- unlist(lapply(stage_levels, function(st) {
  ids <- pb_meta[stage == st, pb_id]
  if (length(ids) > max_show_per_stage) sample(ids, max_show_per_stage) else ids
}), use.names = FALSE)
show_meta <- pb_meta[match(show_ids, pb_id)]
setorder(show_meta, stage)
show_ids <- show_meta$pb_id
sample_mat <- logcpm[genes, show_ids, drop = FALSE]
ann_col <- data.frame(
  stage = show_meta$stage,
  row.names = show_ids
)

z <- t(scale(t(sample_mat)))
z[z > 3] <- 3
z[z < -3] <- -3
z_breaks <- seq(-3, 3, length.out = 101)
z_cols <- colorRampPalette(c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"))(100)

pdf(file.path(out_dir, "heatmap_stage_markers.pdf"), width = 11, height = 11)
pheatmap(
  z,
  scale = "none",
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  annotation_row = row_ann,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  color = z_cols,
  breaks = z_breaks,
  show_colnames = FALSE,
  fontsize_row = 8,
  border_color = NA,
  gaps_col = cumsum(as.integer(table(show_meta$stage)))[-length(stage_levels)],
  gaps_row = cumsum(as.integer(table(row_ann$marker_of)))[-length(stage_levels)],
  main = "Stage marker genes across donor pseudobulks (row z-score, -3 to 3)"
)
dev.off()

message("Wrote heatmap_stage_markers.pdf and heatmap_stage_markers_mean.pdf")
