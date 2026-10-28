#!/usr/bin/env Rscript

# Rebuild the stage-characteristic pseudobulk heatmap after excluding
# immunoglobulin variable-region genes. The original heatmap and tables remain
# unchanged; revised outputs use the suffix "_no_IGV".

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

project_dir <- Sys.getenv("CLL_TME_DIR", unset = path.expand("~/research/CLL_TME"))
setwd(normalizePath(project_dir, mustWork = TRUE))

out_dir <- "CLLCMA_alpha/results/bcell_stage_deseq2"
pb_counts_file <- file.path(out_dir, "pseudobulk_counts.rds")
pb_meta_file <- file.path(out_dir, "pseudobulk_colData.tsv")

stage_levels <- c("HC_B", "MBL", "CLL", "RT")
stage_display <- c(HC_B = "B", MBL = "MBL", CLL = "CLL", RT = "RT")
stage4 <- c("B", "MBL", "CLL", "RT")
stage_colors <- c(
  B = "#0072B2",
  MBL = "#009E73",
  CLL = "#E69F00",
  RT = "#D55E00"
)
min_count <- 10
min_samples_expr <- 3L
n_per_stage <- 50L
igv_pattern <- "^(IGHV|IGKV|IGLV)"

if (!file.exists(pb_counts_file) || !file.exists(pb_meta_file)) {
  stop("Missing pseudobulk input files under ", out_dir)
}

message("Loading saved pseudobulk counts and metadata")
pb_mat <- readRDS(pb_counts_file)
pb_meta <- fread(pb_meta_file)
if ("stage" %in% names(pb_meta) && !"b_stage" %in% names(pb_meta)) {
  pb_meta[, b_stage := stage]
}
pb_meta[b_stage %in% c("Normal_B", "B"), b_stage := "HC_B"]
pb_meta[, b_stage := factor(as.character(b_stage), levels = stage_levels)]
stopifnot(identical(colnames(pb_mat), pb_meta$pb_id))

# Match the four-group sample set used by step_5_bcell_stage_deseq2.R.
pb_meta <- pb_meta[
  b_stage %in% stage_levels &
    !(b_stage == "RT" & disease != "RS")
]
pb_mat <- pb_mat[, pb_meta$pb_id, drop = FALSE]
keep_gene <- rowSums(pb_mat >= min_count) >= min_samples_expr
mat <- pb_mat[keep_gene, , drop = FALSE]

cd <- as.data.frame(pb_meta)
rownames(cd) <- cd$pb_id
cd$stage <- factor(
  unname(stage_display[as.character(cd$b_stage)]),
  levels = stage4
)

message("Computing VST for ", nrow(mat), " genes and ", ncol(mat), " pseudobulks")
dds <- DESeqDataSetFromMatrix(
  countData = round(mat),
  colData = cd,
  design = ~ stage
)
dds <- dds[rowSums(counts(dds)) > 0, ]
vsd <- vst(dds, blind = TRUE)

# Stage-characteristic ranking matches the original heatmap method:
# stage mean log1p(CPM), assign each gene to its maximum stage, then rank by
# separation from the second-highest stage. Exclude IG variable genes before
# ranking so every stage is replenished to 50 non-IGV genes.
lib_hm <- pmax(colSums(mat), 1)
logcpm_hm <- log1p(t(t(mat) / lib_hm * 1e6))
mean_by_stage <- sapply(stage4, function(st) {
  rowMeans(logcpm_hm[, cd$stage == st, drop = FALSE])
})
best <- stage4[max.col(mean_by_stage, ties.method = "first")]
names(best) <- rownames(mean_by_stage)
second <- apply(mean_by_stage, 1L, function(x) sort(x, decreasing = TRUE)[2L])
delta <- mean_by_stage[
  cbind(seq_len(nrow(mean_by_stage)), match(best, stage4))
] - second

eligible <- !grepl(igv_pattern, names(best), ignore.case = TRUE)
marker_list <- lapply(stage4, function(st) {
  idx <- which(best == st & eligible)
  dt <- data.table(
    gene = names(best)[idx],
    marker_of = st,
    mean_logcpm = mean_by_stage[idx, st],
    delta = delta[idx]
  )
  setorder(dt, -delta, -mean_logcpm, gene)
  head(dt, n_per_stage)
})
markers <- rbindlist(marker_list)
markers[, rank_in_stage := seq_len(.N), by = marker_of]

if (markers[, any(grepl(igv_pattern, gene, ignore.case = TRUE))]) {
  stop("IG variable-region genes remain after filtering")
}
if (any(markers[, .N, by = marker_of]$N < n_per_stage)) {
  warning("At least one stage has fewer than ", n_per_stage, " eligible markers")
}

genes <- markers$gene
genes <- genes[genes %in% rownames(assay(vsd))]
mat_hm <- assay(vsd)[genes, , drop = FALSE]
z <- t(scale(t(mat_hm)))
z[!is.finite(z)] <- 0
z[z > 2.5] <- 2.5
z[z < -2.5] <- -2.5

ord <- order(cd$stage, cd$dataset, rownames(cd))
z <- z[, ord, drop = FALSE]
stage_ord <- cd$stage[ord]
marker_stage <- factor(
  markers$marker_of[match(rownames(z), markers$gene)],
  levels = stage4
)

fwrite(
  markers,
  file.path(out_dir, "heatmap_stage_top50_no_IGV_genes.tsv"),
  sep = "\t"
)
fwrite(
  as.data.table(z, keep.rownames = "gene"),
  file.path(out_dir, "heatmap_stage_top50_no_IGV_zscore.tsv"),
  sep = "\t"
)

col_fun <- colorRamp2(
  c(-2.5, 0, 2.5),
  c("#0072B2", "#F7F7F7", "#D55E00")
)
ha_col <- HeatmapAnnotation(
  Stage = stage_ord,
  col = list(Stage = stage_colors),
  annotation_name_side = "left",
  simple_anno_size = unit(3, "mm")
)

ht <- Heatmap(
  z,
  name = "VST\nz-score",
  col = col_fun,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 5),
  top_annotation = ha_col,
  column_split = factor(stage_ord, levels = stage4),
  row_split = marker_stage,
  cluster_column_slices = FALSE,
  cluster_row_slices = FALSE,
  column_title = NULL,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  border = FALSE,
  heatmap_legend_param = list(
    at = c(-2.5, 0, 2.5),
    title = "VST z-score"
  ),
  column_gap = unit(1.5, "mm"),
  row_gap = unit(1.5, "mm"),
  width = unit(11.5, "cm")
)

draw_heatmap <- function() {
  draw(
    ht,
    column_title = "Top 50 stage-characteristic genes per group (IG variable genes excluded)",
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    merge_legend = TRUE,
    padding = unit(c(2, 2, 2, 2), "mm")
  )
}

pdf(
  file.path(out_dir, "heatmap_stage_top50_no_IGV.pdf"),
  width = 8.2, height = 14,
  useDingbats = FALSE
)
draw_heatmap()
dev.off()

png(
  file.path(out_dir, "heatmap_stage_top50_no_IGV.png"),
  width = 2460, height = 4200, res = 300,
  bg = "white"
)
draw_heatmap()
dev.off()

counts <- markers[, .N, by = marker_of]
writeLines(
  c(
    "Figure X. Stage-characteristic pseudobulk expression heatmap excluding immunoglobulin variable-region genes.",
    "Genes beginning IGHV, IGKV, or IGLV were excluded before marker ranking.",
    "The top 50 remaining genes per stage were selected by the original max-in-stage and delta-from-second-stage method.",
    "Rows show gene symbols only; full gene names and functional annotations were removed.",
    "Columns are donor-stage pseudobulks ordered by stage and dataset.",
    "Values are row z-scores of DESeq2 VST expression, clipped to [-2.5, 2.5] for display.",
    paste0(
      "Genes displayed by stage: ",
      paste(paste0(counts$marker_of, "=", counts$N), collapse = ", "),
      "."
    ),
    "The original heatmap and source tables were preserved."
  ),
  file.path(out_dir, "heatmap_stage_top50_no_IGV_caption.txt")
)

message("Wrote simplified no-IGV heatmap to ", normalizePath(out_dir))
