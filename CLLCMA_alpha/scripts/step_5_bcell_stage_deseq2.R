#!/usr/bin/env Rscript

# Pseudobulk DESeq2 for B-lineage disease stages (four-group).
#
# Inputs:  CLLCMA_alpha/integrated_annotated/seu_integrated_final.rds
# Outputs: CLLCMA_alpha/results/bcell_stage_deseq2/
#
# Single DEG model: ~ stage
# Groups: B (normal B / HC), MBL, CLL, RT
# Unit: donor x stage pseudobulk. Not donor-blocked.
# Note: stages are partly cohort-confounded (B ~ PBMCpedia; MBL ~ GSE295491).

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(DESeq2)
})

project_dir <- Sys.getenv(
  "CLL_TME_DIR",
  unset = path.expand("~/research/CLL_TME")
)
project_dir <- normalizePath(project_dir, mustWork = TRUE)
setwd(project_dir)

seu_file <- "CLLCMA_alpha/integrated_annotated/seu_integrated_final.rds"
out_dir <- "CLLCMA_alpha/results/bcell_stage_deseq2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(42L)
min_cells_pb <- 30L
min_count <- 10
min_samples_expr <- 3L
min_reps_per_group <- 2L
padj_cut <- 0.05
lfc_cut <- 1

# Internal labels from cell assignment; display as B / MBL / CLL / RT
stage_levels <- c("HC_B", "MBL", "CLL", "RT")
stage_display <- c(HC_B = "B", MBL = "MBL", CLL = "CLL", RT = "RT")
stage4 <- c("B", "MBL", "CLL", "RT")

# Okabe-Ito (Wong, Nat Methods 2011); color + shape used together on PCA
stage_colors <- c(
  B = "#0072B2",
  MBL = "#009E73",
  CLL = "#E69F00",
  RT = "#D55E00"
)
stage_shapes <- c(B = 16, MBL = 17, CLL = 15, RT = 18)

hc_b_labels <- c(
  "Naive B cell",
  "Memory B cell",
  "Transitional B cell",
  "Effector B cell"
)
nadeu_rs_like <- c(
  "RS", "RT", "RT proliferative", "RT quiescent",
  "TP53INP1hi RT", "MIR155HGhi RT", "CCND2hi RT", "CCND2lo RT",
  "CXCR4loCD27hi RT"
)
nadeu_cll_like <- c(
  "CLL", "CXCR4hiCD27lo", "CXCR4loCD27hi",
  "MIR155HGhi", "MIR155HGhi CLL",
  "CD83hiMIR155HGhi", "CD83loMIR155HGhi",
  "MZB1hiIGHMhiXBP1hi"
)

# ------------------------------ Helpers --------------------------------------

assign_b_stage <- function(cell_type_final, disease, dataset) {
  ct <- as.character(cell_type_final)
  dis <- as.character(disease)
  ds <- as.character(dataset)
  stage <- rep(NA_character_, length(ct))
  stage[ct %in% hc_b_labels & dis == "HC"] <- "HC_B"
  stage[ct == "MBL"] <- "MBL"
  stage[ct == "CLL"] <- "CLL"
  stage[ct == "RT"] <- "RT"
  stage[ct %in% nadeu_cll_like & dis == "CLL"] <- "CLL"
  stage[ct %in% nadeu_rs_like] <- "RT"
  stage[dis == "RS" & ct %in% nadeu_cll_like] <- "RT"
  stage[ds == "PBMCpedia_HC" & !ct %in% hc_b_labels] <- NA_character_
  stage
}

get_counts <- function(seu) {
  tryCatch(
    LayerData(seu, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(seu, assay = "RNA", slot = "counts")
  )
}

pseudobulk_sum <- function(counts, group) {
  group <- droplevels(factor(group))
  mm <- sparseMatrix(
    i = seq_along(group),
    j = as.integer(group),
    x = 1,
    dims = c(length(group), nlevels(group))
  )
  colnames(mm) <- levels(group)
  as.matrix(counts %*% mm)
}

filter_genes <- function(mat) {
  keep <- rowSums(mat >= min_count) >= min_samples_expr
  mat[keep, , drop = FALSE]
}

deseq_table <- function(res) {
  dt <- as.data.table(as.data.frame(res), keep.rownames = "gene")
  setorder(dt, padj, na.last = TRUE)
  dt[, sig := !is.na(padj) & padj < padj_cut & abs(log2FoldChange) >= lfc_cut]
  dt
}

theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
      legend.key = element_blank(),
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(size = base_size - 1, colour = "grey30")
    )
}

volcano_plot <- function(tab, title, file) {
  plot_dt <- copy(tab)
  plot_dt <- plot_dt[is.finite(log2FoldChange) & is.finite(pvalue)]
  plot_dt[, neglog10p := pmin(-log10(pmax(pvalue, 1e-300)), 50)]
  plot_dt[, col := fifelse(
    is.na(padj) | padj >= padj_cut,
    "ns",
    fifelse(log2FoldChange >= lfc_cut, "up",
            fifelse(log2FoldChange <= -lfc_cut, "down", "ns"))
  )]
  plot_dt[, col := factor(col, levels = c("up", "down", "ns"))]
  label_dt <- rbindlist(list(
    plot_dt[col == "up"][order(padj)][seq_len(min(10L, .N))],
    plot_dt[col == "down"][order(padj)][seq_len(min(10L, .N))]
  ), fill = TRUE)

  n_up <- plot_dt[col == "up", .N]
  n_down <- plot_dt[col == "down", .N]
  lab_map <- c(
    up = paste0("Up (", n_up, ")"),
    down = paste0("Down (", n_down, ")"),
    ns = "NS"
  )

  lab_layer <- NULL
  if (requireNamespace("ggrepel", quietly = TRUE) && nrow(label_dt)) {
    lab_layer <- ggrepel::geom_text_repel(
      data = label_dt,
      aes(label = gene),
      size = 2.6,
      max.overlaps = 30,
      min.segment.length = 0,
      colour = "grey20",
      show.legend = FALSE
    )
  }

  p <- ggplot(plot_dt, aes(log2FoldChange, neglog10p, color = col, shape = col)) +
    geom_point(size = 1.1, alpha = 0.75) +
    lab_layer +
    scale_color_manual(
      name = NULL,
      values = c(up = "#D55E00", down = "#0072B2", ns = "grey75"),
      labels = lab_map
    ) +
    scale_shape_manual(
      name = NULL,
      values = c(up = 17, down = 15, ns = 16),
      labels = lab_map
    ) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 2, linewidth = 0.35) +
    geom_hline(yintercept = -log10(padj_cut), linetype = 2, linewidth = 0.35) +
    theme_pub() +
    labs(
      title = title,
      subtitle = paste0(
        "padj < ", padj_cut, ", |log2FC| >= ", lfc_cut, "; design ~ stage"
      ),
      x = expression(log[2] ~ "fold change"),
      y = expression(-log[10] ~ p)
    )
  ggsave(file, p, width = 7.2, height = 5.4, device = cairo_pdf)
}

write_contrast <- function(tab, stem, title) {
  fwrite(tab, file.path(out_dir, paste0(stem, ".tsv")), sep = "\t")
  n_sig <- tab[sig == TRUE, .N]
  message("  ", stem, " significant: ", n_sig)
  volcano_plot(tab, title, file.path(out_dir, paste0(stem, "_volcano.pdf")))
  n_sig
}

# ------------------------------ Load / pseudobulk ----------------------------

pb_counts_file <- file.path(out_dir, "pseudobulk_counts.rds")
pb_meta_file <- file.path(out_dir, "pseudobulk_colData.tsv")
reuse_pb <- identical(Sys.getenv("REUSE_PSEUDOBULK"), "1") &&
  file.exists(pb_counts_file) && file.exists(pb_meta_file)

if (reuse_pb) {
  message("Reusing saved pseudobulk matrices")
  pb_mat <- readRDS(pb_counts_file)
  pb_meta <- fread(pb_meta_file)
  if ("stage" %in% names(pb_meta) && !"b_stage" %in% names(pb_meta)) {
    pb_meta[, b_stage := stage]
  }
  pb_meta[b_stage %in% c("Normal_B", "B"), b_stage := "HC_B"]
  pb_meta[, b_stage := factor(as.character(b_stage), levels = stage_levels)]
  stopifnot(identical(colnames(pb_mat), pb_meta$pb_id))
} else {
  if (!file.exists(seu_file)) {
    stop("Missing: ", seu_file, " — run steps 1–3 first.")
  }

  message("Loading ", seu_file, " ...")
  seu <- readRDS(seu_file)

  need <- c("dataset", "sample", "donor", "disease", "cell_type_final")
  miss <- setdiff(need, colnames(seu[[]]))
  if (length(miss)) {
    stop("Missing metadata: ", paste(miss, collapse = ", "))
  }

  meta <- as.data.table(seu[[]], keep.rownames = "cell_id")
  meta[, b_stage := assign_b_stage(cell_type_final, disease, dataset)]
  meta[, b_stage := factor(b_stage, levels = stage_levels)]
  keep <- meta[!is.na(b_stage)]
  keep[, pb_id := paste(dataset, donor, b_stage, sep = "|")]
  message("B-lineage cells by stage:")
  print(keep[, .N, by = b_stage][order(b_stage)])

  if (!nrow(keep)) {
    stop("No cells assigned to HC_B / MBL / CLL / RT")
  }

  message("Extracting counts for ", nrow(keep), " cells ...")
  counts <- get_counts(seu)
  miss_cells <- setdiff(keep$cell_id, colnames(counts))
  if (length(miss_cells)) {
    stop("Cells in metadata missing from counts: ", length(miss_cells))
  }
  counts <- counts[, keep$cell_id, drop = FALSE]
  rm(seu)
  gc(verbose = FALSE)

  message("Building donor x stage pseudobulk ...")
  pb_mat <- pseudobulk_sum(counts, keep$pb_id)
  rm(counts)
  gc(verbose = FALSE)

  pb_meta <- keep[
    ,
    .(
      n_cells = .N,
      n_samples = uniqueN(sample),
      dataset = dataset[1],
      donor = as.character(donor[1]),
      b_stage = as.character(b_stage[1]),
      disease = names(sort(table(disease), decreasing = TRUE))[1]
    ),
    by = pb_id
  ]
  pb_meta <- pb_meta[match(colnames(pb_mat), pb_id)]
  pb_meta[, b_stage := factor(b_stage, levels = stage_levels)]

  n_drop <- pb_meta[n_cells < min_cells_pb, .N]
  if (n_drop) {
    message("Dropping ", n_drop, " pseudobulks with <", min_cells_pb, " cells")
  }
  keep_pb <- pb_meta$n_cells >= min_cells_pb
  pb_mat <- pb_mat[, keep_pb, drop = FALSE]
  pb_meta <- pb_meta[keep_pb]
  stopifnot(identical(colnames(pb_mat), pb_meta$pb_id))

  fwrite(pb_meta, pb_meta_file, sep = "\t")
  saveRDS(pb_mat, pb_counts_file, compress = FALSE)
  rm(meta, keep)
  gc(verbose = FALSE)
}

# -------------------- Four-group model: ~ stage ------------------------------

pb_meta <- pb_meta[
  b_stage %in% stage_levels &
    !(b_stage == "RT" & disease != "RS")
]
pb_mat <- pb_mat[, pb_meta$pb_id, drop = FALSE]

mat <- filter_genes(pb_mat)
cd <- as.data.frame(pb_meta)
rownames(cd) <- cd$pb_id
cd$stage <- factor(
  unname(stage_display[as.character(cd$b_stage)]),
  levels = stage4
)

message("Four-group sample counts (stage x dataset):")
print(table(cd$stage, cd$dataset))
if (any(table(cd$stage) < min_reps_per_group)) {
  stop("At least one stage has <", min_reps_per_group, " pseudobulks")
}

fwrite(
  as.data.table(cd)[, .(pb_id, dataset, donor, stage, b_stage, disease, n_cells)],
  file.path(out_dir, "deseq2_four_group_samples.tsv"),
  sep = "\t"
)
fwrite(
  as.data.table(cd)[
    ,
    .(n_pseudobulk = .N, n_donors = uniqueN(donor), n_cells = sum(n_cells)),
    by = .(stage, dataset)
  ][order(stage, dataset)],
  file.path(out_dir, "deseq2_four_group_summary.tsv"),
  sep = "\t"
)

message("DESeq2: ~ stage  n=", ncol(mat), " samples, ", nrow(mat), " genes")
dds <- DESeqDataSetFromMatrix(
  countData = round(mat),
  colData = cd,
  design = ~ stage
)
dds <- dds[rowSums(counts(dds)) > 0, ]
dds <- DESeq(dds, quiet = TRUE)

pairwise <- list(
  MBL_vs_B = c("MBL", "B"),
  CLL_vs_B = c("CLL", "B"),
  RT_vs_B = c("RT", "B"),
  CLL_vs_MBL = c("CLL", "MBL"),
  RT_vs_MBL = c("RT", "MBL"),
  RT_vs_CLL = c("RT", "CLL")
)

contrast_n <- list()
for (nm in names(pairwise)) {
  res <- results(
    dds,
    contrast = c("stage", pairwise[[nm]][1], pairwise[[nm]][2]),
    alpha = padj_cut
  )
  tab <- deseq_table(res)
  title <- paste0(pairwise[[nm]][1], " vs ", pairwise[[nm]][2])
  contrast_n[[nm]] <- write_contrast(tab, paste0("deseq2_", nm), title)
}

message("LRT for stage (~ stage vs ~ 1) ...")
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1, quiet = TRUE)
lrt_tab <- deseq_table(results(dds_lrt, alpha = padj_cut))
lrt_tab[, sig := !is.na(padj) & padj < padj_cut]
fwrite(lrt_tab, file.path(out_dir, "deseq2_four_group_LRT.tsv"), sep = "\t")
contrast_n$four_group_LRT <- lrt_tab[sig == TRUE, .N]
message("  LRT significant (padj<", padj_cut, "): ", contrast_n$four_group_LRT)

vsd <- vst(dds, blind = TRUE)
pca_obj <- plotPCA(vsd, intgroup = c("stage", "dataset"), returnData = TRUE)
pca_dt <- as.data.table(pca_obj)
pct <- attr(pca_obj, "percentVar")
pct1 <- round(100 * pct[1], 1)
pct2 <- round(100 * pct[2], 1)

p_pca <- ggplot(pca_dt, aes(PC1, PC2, color = stage, shape = stage)) +
  geom_point(size = 2.6, alpha = 0.9) +
  scale_color_manual(values = stage_colors) +
  scale_shape_manual(values = stage_shapes) +
  theme_pub() +
  labs(
    title = "Pseudobulk PCA (VST)",
    subtitle = "Four groups; design ~ stage (not donor-blocked)",
    x = paste0("PC1 (", pct1, "%)"),
    y = paste0("PC2 (", pct2, "%)"),
    color = "Stage",
    shape = "Stage"
  )
ggsave(
  file.path(out_dir, "pca_four_group.pdf"),
  p_pca,
  width = 6.8,
  height = 5.2,
  device = cairo_pdf
)

if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
    requireNamespace("circlize", quietly = TRUE)) {
  n_per_stage <- 50L
  # Mean log1p CPM by stage; pick genes whose max is that stage
  lib_hm <- pmax(colSums(mat), 1)
  logcpm_hm <- log1p(t(t(mat) / lib_hm * 1e6))
  mean_by_stage <- sapply(stage4, function(st) {
    rowMeans(logcpm_hm[, cd$stage == st, drop = FALSE])
  })
  best <- stage4[max.col(mean_by_stage, ties.method = "first")]
  names(best) <- rownames(mean_by_stage)
  second <- apply(mean_by_stage, 1L, function(x) sort(x, decreasing = TRUE)[2L])
  delta <- mean_by_stage[cbind(seq_len(nrow(mean_by_stage)), match(best, stage4))] - second

  marker_list <- lapply(stage4, function(st) {
    idx <- which(best == st)
    dt <- data.table(
      gene = names(best)[idx],
      marker_of = st,
      mean_logcpm = mean_by_stage[idx, st],
      delta = delta[idx]
    )
    setorder(dt, -delta, -mean_logcpm)
    head(dt, n_per_stage)
  })
  markers <- rbindlist(marker_list)
  markers[, rank_in_stage := seq_len(.N), by = marker_of]
  fwrite(markers, file.path(out_dir, "heatmap_stage_top50_genes.tsv"), sep = "\t")
  message(
    "Stage top-", n_per_stage, " genes: ",
    paste(markers[, .N, by = marker_of]$N, collapse = "/")
  )

  genes <- markers$gene
  genes <- genes[genes %in% rownames(assay(vsd))]
  mat_hm <- assay(vsd)[genes, , drop = FALSE]
  z <- t(scale(t(mat_hm)))
  z[!is.finite(z)] <- 0
  z[z > 2.5] <- 2.5
  z[z < -2.5] <- -2.5

  ord <- order(cd$stage, cd$dataset)
  z <- z[, ord, drop = FALSE]
  stage_ord <- cd$stage[ord]

  gene_ann <- data.table(
    gene = rownames(z),
    marker_of = markers$marker_of[match(rownames(z), markers$gene)],
    description = NA_character_
  )
  go_sets <- list(
    "Immune / lymphocyte" = c(
      "GO:0002376", "GO:0006955", "GO:0045087", "GO:0002250", "GO:0006959"
    ),
    "Cell cycle / proliferation" = c(
      "GO:0007049", "GO:0022402", "GO:0051301", "GO:0000278"
    ),
    "Stress / apoptosis" = c(
      "GO:0006915", "GO:0012501", "GO:0006950", "GO:0034976"
    ),
    "Transcription / chromatin" = c(
      "GO:0006351", "GO:0006355", "GO:0016569", "GO:0006325", "GO:0006338"
    ),
    "Signaling" = c("GO:0007165", "GO:0023052", "GO:0007154"),
    "RNA / ribosome" = c(
      "GO:0006412", "GO:0006396", "GO:0008380", "GO:0016071"
    ),
    "Cytoskeleton / adhesion" = c(
      "GO:0007010", "GO:0007155", "GO:0030030", "GO:0030155"
    ),
    "Metabolism" = c(
      "GO:0006091", "GO:0006629", "GO:0005975", "GO:0009058", "GO:0006119"
    )
  )
  class_levels <- c(names(go_sets), "Other / unknown")

  if (requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
      requireNamespace("AnnotationDbi", quietly = TRUE)) {
    map <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = gene_ann$gene,
      columns = c("SYMBOL", "GENENAME"),
      keytype = "SYMBOL"
    )
    map <- as.data.table(map)[!is.na(SYMBOL)]
    map <- map[!duplicated(SYMBOL)]
    gene_ann[map, description := i.GENENAME, on = c(gene = "SYMBOL")]

    goall <- as.data.table(AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = gene_ann$gene,
      columns = c("GOALL", "ONTOLOGY"),
      keytype = "SYMBOL"
    ))
    goall <- goall[!is.na(GOALL) & ONTOLOGY == "BP"]
    go_by_gene <- split(goall$GOALL, goall$SYMBOL)
    assign_from_go <- function(g) {
      terms <- go_by_gene[[g]]
      if (is.null(terms)) {
        return(NA_character_)
      }
      for (cl in names(go_sets)) {
        if (any(go_sets[[cl]] %in% terms)) {
          return(cl)
        }
      }
      NA_character_
    }
    gene_ann[, fn_go := vapply(gene, assign_from_go, character(1))]
  } else {
    gene_ann[, fn_go := NA_character_]
  }
  gene_ann[is.na(description) | description == "", description := "unannotated"]

  assign_from_kw <- function(gene, desc) {
    g <- toupper(gene)
    d <- tolower(desc)
    if (grepl("^IGH|^IGK|^IGL|^HLA-|^CD[0-9]", g) ||
        grepl("immun|lymphocyte|cytokine|chemokine|antigen|leukocyte|neutrophil", d)) {
      return("Immune / lymphocyte")
    }
    if (grepl("^HIST|^H[1-4][-_]", g) ||
        grepl("transcri|chromatin|histone", d)) {
      return("Transcription / chromatin")
    }
    if (grepl("cell cycle|mitos|cyclin|prolif", d)) {
      return("Cell cycle / proliferation")
    }
    if (grepl("kinase|phosphatase|signal|receptor", d)) {
      return("Signaling")
    }
    if (grepl("ribosom|translation|mrna|splic|helicase", d)) {
      return("RNA / ribosome")
    }
    if (grepl("mitochond|metabol|lipid|glycol", d)) {
      return("Metabolism")
    }
    if (grepl("apopt|stress|caspase", d)) {
      return("Stress / apoptosis")
    }
    if (grepl("cytoskeleton|actin|adhesion|integrin", d)) {
      return("Cytoskeleton / adhesion")
    }
    "Other / unknown"
  }
  gene_ann[, fn_kw := mapply(assign_from_kw, gene, description)]
  gene_ann[, fn_class := fifelse(!is.na(fn_go), fn_go, fn_kw)]
  gene_ann[, fn_class := factor(fn_class, levels = class_levels)]
  gene_ann[, desc_short := fifelse(
    nchar(description) > 42L,
    paste0(substr(description, 1L, 39L), "..."),
    description
  )]

  fwrite(
    gene_ann[, .(gene, marker_of, fn_class, description)],
    file.path(out_dir, "heatmap_stage_top50_gene_annotation.tsv"),
    sep = "\t"
  )
  fwrite(
    as.data.table(z, keep.rownames = "gene"),
    file.path(out_dir, "heatmap_stage_top50_zscore.tsv"),
    sep = "\t"
  )

  fn_cols <- c(
    "Immune / lymphocyte" = "#0072B2",
    "Cell cycle / proliferation" = "#D55E00",
    "Transcription / chromatin" = "#CC79A7",
    "Signaling" = "#E69F00",
    "Metabolism" = "#009E73",
    "RNA / ribosome" = "#56B4E9",
    "Stress / apoptosis" = "#F0E442",
    "Cytoskeleton / adhesion" = "#000000",
    "Other / unknown" = "#999999"
  )

  col_fun <- circlize::colorRamp2(
    c(-2.5, 0, 2.5),
    c("#0072B2", "#F7F7F7", "#D55E00")
  )
  ha_col <- ComplexHeatmap::HeatmapAnnotation(
    stage = stage_ord,
    col = list(stage = stage_colors),
    annotation_name_side = "left",
    simple_anno_size = grid::unit(3, "mm")
  )
  ha_row <- ComplexHeatmap::rowAnnotation(
    Stage = gene_ann$marker_of,
    Function = gene_ann$fn_class,
    Description = ComplexHeatmap::anno_text(
      gene_ann$desc_short,
      gp = grid::gpar(fontsize = 3.8),
      location = 0,
      just = "left",
      width = grid::unit(4.8, "cm")
    ),
    col = list(Stage = stage_colors, Function = fn_cols),
    annotation_name_gp = grid::gpar(fontsize = 8),
    gap = grid::unit(1, "mm")
  )

  ht <- ComplexHeatmap::Heatmap(
    z,
    name = "VST\nz-score",
    col = col_fun,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 3.8),
    top_annotation = ha_col,
    right_annotation = ha_row,
    column_split = factor(stage_ord, levels = stage4),
    row_split = factor(gene_ann$marker_of, levels = stage4),
    cluster_column_slices = FALSE,
    cluster_row_slices = FALSE,
    column_title = NULL,
    row_title_rot = 0,
    row_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    border = FALSE,
    heatmap_legend_param = list(
      at = c(-2.5, 0, 2.5),
      title = "VST z-score"
    ),
    column_gap = grid::unit(1.5, "mm"),
    row_gap = grid::unit(1.5, "mm"),
    width = grid::unit(10, "cm")
  )

  pdf(file.path(out_dir, "heatmap_stage_top50.pdf"), width = 12, height = 10)
  ComplexHeatmap::draw(
    ht,
    column_title = "Top 50 stage-characteristic genes per group (B / MBL / CLL / RT)",
    column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    merge_legend = TRUE,
    padding = grid::unit(c(2, 2, 2, 2), "mm")
  )
  dev.off()
}


writeLines(
  c(
    "# B-lineage stage pseudobulk DESeq2 (four-group)",
    "",
    "Model: `~ stage`",
    "Groups: B (normal B / HC), MBL, CLL, RT",
    "Unit: donor x stage pseudobulk. Not donor-blocked.",
    "RT restricted to disease == RS.",
    "",
    "Caveat: stages are partly cohort-confounded (B mostly PBMCpedia;",
    "MBL mostly GSE295491). Treat B vs disease contrasts as exploratory.",
    "",
    paste0("Min cells / pseudobulk: ", min_cells_pb),
    paste0(
      "Significance: padj < ", padj_cut,
      ", |log2FC| >= ", lfc_cut, " (Wald); LRT uses padj only."
    ),
    "",
    "Palette: Okabe-Ito (Wong, Nat Methods 2011). PCA uses color + shape.",
    "Heatmap values: row z-score of VST, clipped to [-2.5, 2.5] for display.",
    "Heatmap genes: top 50 per stage (max-in-stage + delta). Right-side Function uses GOALL BP ancestors (priority) with keyword/symbol fallback;",
    "Description is NCBI GENENAME (truncated on plot; full text in annotation TSV).",
    "Publisher formatting: provisional (journal not specified).",
    "",
    "Outputs:",
    "- `deseq2_four_group_samples.tsv`, `deseq2_four_group_summary.tsv`",
    "- `deseq2_MBL_vs_B.tsv`, `deseq2_CLL_vs_B.tsv`, `deseq2_RT_vs_B.tsv`",
    "- `deseq2_CLL_vs_MBL.tsv`, `deseq2_RT_vs_MBL.tsv`, `deseq2_RT_vs_CLL.tsv`",
    "- `deseq2_four_group_LRT.tsv`",
    "- `pca_four_group.pdf`, `heatmap_stage_top50.pdf`",
    "- `heatmap_stage_top50_gene_annotation.tsv`",
    "",
    "Significant gene counts:",
    paste(sprintf("- %s: %s", names(contrast_n), unlist(contrast_n)), collapse = "\n")
  ),
  file.path(out_dir, "README.md")
)

message("Done -> ", normalizePath(out_dir))
