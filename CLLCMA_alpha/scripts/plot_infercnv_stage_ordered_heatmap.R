#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

project_dir <- Sys.getenv("CLL_TME_DIR", unset = path.expand("~/research/CLL_TME"))
setwd(normalizePath(project_dir, mustWork = TRUE))

out_dir <- "CLLCMA_alpha/results/bcell_stage_deseq2/publication_figures/infercnv_bcell_stage"
fig_dir <- file.path(out_dir, "stage_ordered_panels")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

stage_order <- c("Normal_B", "MBL", "CLL", "RT")
stage_cols <- c(Normal_B = "#4C78A8", MBL = "#54A24B", CLL = "#F58518", RT = "#E45756")
chr_levels <- paste0("chr", 1:22)
set.seed(42L)

obj <- readRDS(file.path(out_dir, "infercnv_obj_final.rds"))
expr <- obj@expr.data
anno <- fread(file.path(out_dir, "infercnv_cell_annotations.tsv"), header = FALSE)
setnames(anno, c("cell_id", "stage"))
anno <- anno[cell_id %in% colnames(expr)]
anno[, stage := factor(stage, levels = stage_order)]
setorder(anno, stage, cell_id)
expr <- expr[, anno$cell_id, drop = FALSE]

gene_order <- fread(file.path(out_dir, "infercnv_gene_order.tsv"), header = FALSE)
setnames(gene_order, c("gene", "chr", "start", "end"))
gene_order <- gene_order[gene %in% rownames(expr) & chr %in% chr_levels]
gene_order[, chr := factor(chr, levels = chr_levels)]
setorder(gene_order, chr, start, end, gene)
expr <- expr[gene_order$gene, , drop = FALSE]

# inferCNV final expr.data is fold-change-like and centered near 1.
# Convert to log2 deviation for symmetric visualization, then bin by genomic order.
log_dev <- log2(pmax(expr, 1e-6))
log_dev[log_dev > 0.5] <- 0.5
log_dev[log_dev < -0.5] <- -0.5

bin_size <- 80L
gene_order[, bin := paste0(as.character(chr), ":", floor((seq_len(.N) - 1L) / bin_size) + 1L), by = chr]
bin_info <- gene_order[, .(
  chr = as.character(chr[1]),
  start = min(start),
  end = max(end),
  n_genes = .N
), by = bin]
bin_info[, chr_order := factor(chr, levels = chr_levels)]
setorder(bin_info, chr_order, start)

bin_mat_cells <- sapply(bin_info$bin, function(b) {
  genes <- gene_order[bin == b, gene]
  colMeans(log_dev[genes, , drop = FALSE], na.rm = TRUE)
})
rownames(bin_mat_cells) <- colnames(log_dev)
bin_mat_cells <- bin_mat_cells[anno$cell_id, , drop = FALSE]

# A compact per-cell view: balanced subsampling from the inferCNV-processed matrix.
max_cells_per_stage_plot <- as.integer(Sys.getenv("INFERCNV_STAGE_ORDERED_MAX_CELLS", "500"))
plot_anno <- rbindlist(lapply(stage_order, function(st) {
  dt <- anno[stage == st]
  if (nrow(dt) > max_cells_per_stage_plot) dt <- dt[sample.int(nrow(dt), max_cells_per_stage_plot)]
  dt
}))
setorder(plot_anno, stage, cell_id)
cell_plot <- bin_mat_cells[plot_anno$cell_id, , drop = FALSE]
cell_plot[cell_plot > 0.35] <- 0.35
cell_plot[cell_plot < -0.35] <- -0.35

# Stage mean view using all selected cells.
stage_mean <- sapply(stage_order, function(st) {
  ids <- anno[stage == st, cell_id]
  colMeans(bin_mat_cells[ids, , drop = FALSE], na.rm = TRUE)
})
stage_mean <- t(stage_mean)
stage_mean[stage_mean > 0.25] <- 0.25
stage_mean[stage_mean < -0.25] <- -0.25

fwrite(as.data.table(stage_mean, keep.rownames = "stage"), file.path(fig_dir, "infercnv_stage_mean_binned_log2_deviation.tsv"), sep = "\t")
fwrite(bin_info[, !"chr_order"], file.path(fig_dir, "infercnv_stage_ordered_bins.tsv"), sep = "\t")

col_fun_cells <- colorRamp2(c(-0.35, 0, 0.35), c("#2166AC", "white", "#B2182B"))
col_fun_mean <- colorRamp2(c(-0.25, 0, 0.25), c("#2166AC", "white", "#B2182B"))
col_split <- factor(bin_info$chr, levels = chr_levels)

ha_cells <- rowAnnotation(
  Stage = plot_anno$stage,
  col = list(Stage = stage_cols),
  simple_anno_size = unit(3, "mm"),
  annotation_name_gp = gpar(fontsize = 8)
)
ht_cells <- Heatmap(
  cell_plot,
  name = "inferCNV\nlog2 dev.",
  col = col_fun_cells,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = plot_anno$stage,
  column_split = col_split,
  show_row_names = FALSE,
  show_column_names = FALSE,
  left_annotation = ha_cells,
  column_title = "inferCNV-processed single-cell signal ordered by genomic position",
  row_title = NULL
)

pdf(file.path(fig_dir, "infercnv_stage_ordered_cells_binned.pdf"), width = 14, height = 7)
draw(ht_cells, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
png(file.path(fig_dir, "infercnv_stage_ordered_cells_binned.png"), width = 4200, height = 2100, res = 300)
draw(ht_cells, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

ht_mean <- Heatmap(
  stage_mean,
  name = "Mean inferCNV\nlog2 dev.",
  col = col_fun_mean,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = col_split,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  show_column_names = FALSE,
  column_title = "Stage-mean inferCNV signal: Normal_B -> MBL -> CLL -> RT",
  row_title = NULL
)

pdf(file.path(fig_dir, "infercnv_stage_ordered_stage_means_binned.pdf"), width = 14, height = 3.5)
draw(ht_mean, heatmap_legend_side = "right")
dev.off()
png(file.path(fig_dir, "infercnv_stage_ordered_stage_means_binned.png"), width = 4200, height = 1050, res = 300)
draw(ht_mean, heatmap_legend_side = "right")
dev.off()

caption <- c(
  "Figure X. Stage-ordered inferCNV heatmap of B-cell disease progression.",
  "inferCNV was run on raw single-cell counts from donor-balanced B-lineage cells using Normal_B cells as the reference group and MBL, CLL, and RT cells as observations.",
  "The displayed panel uses the final inferCNV-processed expression matrix, converted to log2 deviation from the diploid-centered reference and smoothed into 80-gene genomic bins for readability.",
  paste0("The per-cell panel shows up to ", max_cells_per_stage_plot, " cells per stage, ordered Normal_B, MBL, CLL, and RT; the stage-mean panel summarizes all inferCNV-selected cells."),
  "Red indicates inferred relative gain and blue indicates inferred relative loss compared with the Normal_B reference. These are RNA-based inferCNV estimates and should be interpreted with matched genomic CNA data when available."
)
writeLines(caption, file.path(fig_dir, "infercnv_stage_ordered_caption.txt"))
message("Wrote stage-ordered inferCNV panels to: ", normalizePath(fig_dir))
