#!/usr/bin/env Rscript

# Plots for splicing_outlier_tme_association.py

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

project_dir <- Sys.getenv(
  "CLL_TME_DIR",
  unset = path.expand("~/research/CLL_TME")
)
project_dir <- normalizePath(project_dir, mustWork = TRUE)
setwd(project_dir)

deconv_method <- tolower(Sys.getenv("DECONV_METHOD", unset = "bayesprism"))
method_label <- switch(
  deconv_method,
  cibersort = "CIBERSORT",
  bayesprism = "BayesPrism",
  mcp_counter = "MCP-counter",
  quantiseq = "quanTIseq",
  xcell = "xCell",
  deconv_method
)

out_dir <- file.path(
  "CLLCMA_alpha", "results", "deconvolution", "cll_broad_2022",
  "splicing_outlier_tme_association", deconv_method
)
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
message("Method: ", deconv_method, " -> ", normalizePath(out_dir, mustWork = FALSE))

counts <- fread(file.path(out_dir, "outlier_gene_counts_per_sample.tsv"))
assoc <- fread(file.path(out_dir, "association_all.tsv"))
micro <- fread(file.path(
  "CLLCMA_alpha", "results", "deconvolution", "cll_broad_2022",
  paste0(deconv_method, "_sample_by_feature.tsv")
))
merged <- merge(counts, micro, by = "SAMPLE_ID")

as_types <- c(
  "alt_3prime", "alt_5prime", "exon_skip", "intron_retention",
  "mult_exon_skip", "mutex_exons", "all_types_union"
)
as_types <- intersect(as_types, names(counts))

# ---- Count distributions ----
long_c <- melt(
  counts,
  id.vars = "SAMPLE_ID",
  measure.vars = as_types,
  variable.name = "AS_TYPE",
  value.name = "N_OUTLIER_GENES"
)
p_hist <- ggplot(long_c, aes(x = N_OUTLIER_GENES, fill = AS_TYPE)) +
  geom_histogram(bins = 40, color = NA, alpha = 0.9) +
  facet_wrap(~AS_TYPE, scales = "free", ncol = 3) +
  scale_fill_brewer(palette = "Dark2", guide = "none") +
  labs(
    title = "Outlier splicing genes per sample",
    subtitle = sprintf("%s | gene-collapsed two-sided P10/P90 + |PSI-median|>0.10", method_label),
    x = "Outlier gene count",
    y = "Samples"
  ) +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 8))
ggsave(
  file.path(out_dir, "outlier_gene_count_distributions.pdf"),
  p_hist, width = 10, height = 7
)

# ---- Association heatmap (adjusted beta) ----
heat <- assoc[BURDEN %in% as_types]
heat[, BURDEN := factor(BURDEN, levels = as_types)]
heat[, LABEL := fifelse(
  !is.na(FDR_ADJ) & FDR_ADJ < 0.05, "*",
  fifelse(!is.na(FDR_ADJ) & FDR_ADJ < 0.10, "\u00b7", "")
)]
# cap color
lim <- quantile(abs(heat$BETA_COUNT_Z), 0.98, na.rm = TRUE)
lim <- max(lim, 0.01)
heat[, BETA_PLOT := pmax(pmin(BETA_COUNT_Z, lim), -lim)]

p_heat <- ggplot(heat, aes(x = CELL_TYPE, y = BURDEN, fill = BETA_PLOT)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = LABEL), size = 3.2, vjust = 0.45) +
  scale_fill_gradient2(
    low = "#0072B2", mid = "white", high = "#D55E00",
    midpoint = 0, name = "Beta\n(z-count)"
  ) +
  labs(
    title = "Splicing-outlier gene burden vs TME fractions",
    subtitle = sprintf(
      "%s | LM: fraction ~ z(count) + IGHV + sex + age/10; * FDR<0.05; · FDR<0.10",
      method_label
    ),
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 8),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank()
  )
ggsave(
  file.path(out_dir, "association_heatmap_adjusted.pdf"),
  p_heat, width = 12, height = 5.5
)

# ---- Scatter for top adjusted hits ----
top <- assoc[is.finite(P_ADJ)][order(P_ADJ)][, head(.SD, 12L)]
n_plotted <- 0L
for (i in seq_len(nrow(top))) {
  r <- top[i]
  if (!r$BURDEN %in% names(merged) || !r$CELL_TYPE %in% names(merged)) {
    next
  }
  dt <- merged[, c("SAMPLE_ID", r$BURDEN, r$CELL_TYPE), with = FALSE]
  setnames(dt, c("SAMPLE_ID", "COUNT", "FRACTION"))
  dt <- dt[is.finite(COUNT) & is.finite(FRACTION)]
  p <- ggplot(dt, aes(x = COUNT, y = FRACTION)) +
    geom_point(alpha = 0.45, size = 1.6, color = "#333333") +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, color = "#D55E00") +
    labs(
      title = sprintf("%s · %s", r$BURDEN, r$CELL_TYPE),
      subtitle = sprintf(
        "rho=%+.2f | beta_z=%+.3g | adj-p=%.2g FDR=%s",
        r$SPEARMAN_RHO, r$BETA_COUNT_Z, r$P_ADJ,
        ifelse(is.na(r$FDR_ADJ), "NA", sprintf("%.2g", r$FDR_ADJ))
      ),
      x = "Outlier gene count",
      y = "TME fraction"
    ) +
    theme_bw(base_size = 11)
  ggsave(
    file.path(
      plot_dir,
      sprintf("scatter__%s__%s.pdf", gsub("[^A-Za-z0-9._-]+", "_", r$BURDEN),
              gsub("[^A-Za-z0-9._-]+", "_", r$CELL_TYPE))
    ),
    p, width = 4.2, height = 4.0
  )
  n_plotted <- n_plotted + 1L
}

message("Plotted ", n_plotted, " scatters -> ", normalizePath(plot_dir))

# ---- Per-gene volcano (faceted by AS type; color by cell type for top hits) ----
per_gene_file <- file.path(out_dir, "association_per_gene_all.tsv.gz")
if (file.exists(per_gene_file)) {
  pg <- fread(per_gene_file)
  pg <- pg[is.finite(P_ADJ) & is.finite(BETA_ADJ)]
  pg[, NEGLOG10_P := -log10(pmax(P_ADJ, 1e-300))]
  pg[, SIG := !is.na(FDR_ADJ) & FDR_ADJ < 0.05]
  # Cap extreme for display
  pg[, NEGLOG10_P_PLOT := pmin(NEGLOG10_P, quantile(NEGLOG10_P, 0.999, na.rm = TRUE))]

  p_vol <- ggplot(pg, aes(x = BETA_ADJ, y = NEGLOG10_P_PLOT)) +
    geom_point(
      data = pg[SIG == FALSE],
      color = "#B0B0B0", alpha = 0.15, size = 0.4
    ) +
    geom_point(
      data = pg[SIG == TRUE],
      aes(color = CELL_TYPE), alpha = 0.55, size = 0.7
    ) +
    facet_wrap(~AS_TYPE, scales = "free_y", ncol = 3) +
    labs(
      title = "Per-gene splicing-outlier associations with TME",
      subtitle = sprintf(
        "%s | LM: fraction ~ outlier + IGHV + sex + age/10; colored points FDR<0.05",
        method_label
      ),
      x = "Adjusted beta (outlier vs not)",
      y = expression(-log[10](P[adj])),
      color = "Cell type"
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(0.35, "cm"),
      strip.text = element_text(face = "bold", size = 8)
    ) +
    guides(color = guide_legend(nrow = 3, override.aes = list(size = 2, alpha = 1)))
  ggsave(
    file.path(out_dir, "association_per_gene_volcano.pdf"),
    p_vol, width = 12, height = 8
  )

  # Hits by cell type bar
  if (nrow(pg[SIG == TRUE])) {
    hit_ct <- pg[SIG == TRUE, .N, by = .(AS_TYPE, CELL_TYPE)]
    p_bar <- ggplot(hit_ct, aes(x = CELL_TYPE, y = N, fill = AS_TYPE)) +
      geom_col(position = "dodge") +
      labs(
        title = "Significant per-gene outlier-TME hits (FDR<0.05)",
        x = NULL, y = "No. significant gene associations", fill = "AS type"
      ) +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 8))
    ggsave(
      file.path(out_dir, "association_per_gene_hits_by_celltype.pdf"),
      p_bar, width = 11, height = 5
    )
  }

  # ---- Heatmap: top N genes × cell types ----
  n_top <- as.integer(Sys.getenv("TOP_N_GENES", unset = "100"))
  if (is.na(n_top) || n_top < 1L) {
    n_top <- 100L
  }
  # Prefer a single AS matrix for coherent rows; fall back to best-AS per gene
  heat_as <- Sys.getenv("HEATMAP_AS_TYPE", unset = "all_types_union")
  if (!heat_as %in% pg$AS_TYPE) {
    heat_as <- "all_types_union"
  }

  # Rank genes by best adjusted p across all AS types
  best_any <- pg[, .SD[which.min(P_ADJ)], by = GENE][order(P_ADJ)]
  top_genes <- best_any[seq_len(min(n_top, .N)), GENE]

  heat_src <- pg[AS_TYPE == heat_as & GENE %in% top_genes]
  # Genes missing from preferred AS type: fill from that gene's best AS type
  missing <- setdiff(top_genes, unique(heat_src$GENE))
  if (length(missing)) {
    fill_as <- best_any[GENE %in% missing, .(GENE, AS_TYPE)]
    extra <- merge(pg, fill_as, by = c("GENE", "AS_TYPE"))
    heat_src <- rbind(heat_src, extra, use.names = TRUE, fill = TRUE)
  }

  # Gene symbols
  gene_ids <- unique(heat_src$GENE)
  sym_map <- data.table(GENE = gene_ids, SYMBOL = gene_ids)
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    mapped <- tryCatch(
      as.data.table(AnnotationDbi::select(
        org.Hs.eg.db,
        keys = gene_ids,
        keytype = "ENSEMBL",
        columns = "SYMBOL"
      )),
      error = function(e) NULL
    )
    if (!is.null(mapped) && nrow(mapped)) {
      mapped <- unique(mapped[!is.na(SYMBOL) & SYMBOL != "", .(GENE = ENSEMBL, SYMBOL)])
      sym_map <- merge(sym_map[, .(GENE)], mapped, by = "GENE", all.x = TRUE)
      sym_map[is.na(SYMBOL) | SYMBOL == "", SYMBOL := GENE]
    }
  }
  # Order genes by best overall p (top at top of heatmap)
  gene_order <- data.table(GENE = top_genes)
  gene_order[, GENE_ORD := seq_len(.N)]
  heat_src <- merge(heat_src, sym_map, by = "GENE", all.x = TRUE)
  heat_src <- merge(heat_src, gene_order, by = "GENE", all.x = TRUE)
  heat_src[is.na(SYMBOL) | SYMBOL == "", SYMBOL := GENE]
  dup_sym <- heat_src[, .(n_gene = uniqueN(GENE)), by = SYMBOL][n_gene > 1L, SYMBOL]
  heat_src[, LABEL := fifelse(
    SYMBOL %in% dup_sym,
    sprintf("%s (%s)", SYMBOL, GENE),
    SYMBOL
  )]
  label_order <- unique(heat_src[, .(GENE, LABEL, GENE_ORD)])
  setorder(label_order, GENE_ORD)
  heat_src[, LABEL := factor(LABEL, levels = rev(label_order$LABEL))]

  ct_order <- c(
    "CLL", "RS",
    "CD4_T_naive", "CD4_T_memory", "CD8_T_naive", "CD8_T_memory", "CD8aa", "DN_T",
    "T_proliferating", "Treg", "gdT",
    "NK_CD56dim", "NK_CD56bright", "NK_proliferating",
    "Monocyte_CD14", "Monocyte_CD16", "Monocyte_intermediate",
    "cDC1", "cDC2", "ASDC", "ILC", "Progenitor"
  )
  ct_present <- intersect(ct_order, unique(heat_src$CELL_TYPE))
  ct_present <- c(ct_present, setdiff(unique(heat_src$CELL_TYPE), ct_present))
  heat_src[, CELL_TYPE := factor(CELL_TYPE, levels = ct_present)]

  heat_src[, MARK := fifelse(
    !is.na(FDR_ADJ) & FDR_ADJ < 0.05, "*",
    fifelse(!is.na(FDR_ADJ) & FDR_ADJ < 0.10, "\u00b7", "")
  )]
  # Color by Cliff's delta (bounded); fall back to beta if needed
  fill_col <- if (all(is.finite(heat_src$CLIFFS_DELTA))) "CLIFFS_DELTA" else "BETA_ADJ"
  lim <- if (fill_col == "CLIFFS_DELTA") {
    0.85
  } else {
    as.numeric(quantile(abs(heat_src[[fill_col]]), 0.98, na.rm = TRUE))
  }
  heat_src[, FILL := pmax(pmin(get(fill_col), lim), -lim)]

  p_gene_heat <- ggplot(heat_src, aes(x = CELL_TYPE, y = LABEL, fill = FILL)) +
    geom_tile(color = "white", linewidth = 0.15) +
    geom_text(aes(label = MARK), size = 2.0, vjust = 0.4, color = "#222222") +
    scale_fill_gradient2(
      low = "#0072B2", mid = "white", high = "#D55E00",
      midpoint = 0, limits = c(-lim, lim),
      name = if (fill_col == "CLIFFS_DELTA") "Cliff's delta" else "Beta"
    ) +
    labs(
      title = sprintf("Top %d splicing-outlier genes associated with TME", n_top),
      subtitle = sprintf(
        "%s | matrix=%s; genes ranked by best adj-p; * FDR<0.05; \u00b7 FDR<0.10",
        method_label, heat_as
      ),
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 7.5),
      axis.text.y = element_text(size = 6.5, face = "italic"),
      panel.grid = element_blank(),
      legend.position = "right"
    )

  out_heat <- file.path(out_dir, sprintf("association_per_gene_top%d_heatmap.pdf", n_top))
  ggsave(out_heat, p_gene_heat, width = 11.5, height = max(10, 0.11 * n_top + 2.5))
  ggsave(
    sub("\\.pdf$", ".png", out_heat),
    p_gene_heat,
    width = 11.5,
    height = max(10, 0.11 * n_top + 2.5),
    dpi = 300,
    bg = "white"
  )

  # Sidecar table of genes shown
  gene_tab <- merge(
    best_any[GENE %in% top_genes, .(
      GENE, BEST_AS_TYPE = AS_TYPE, BEST_CELL_TYPE = CELL_TYPE,
      BEST_P_ADJ = P_ADJ, BEST_FDR_ADJ = FDR_ADJ,
      BEST_BETA = BETA_ADJ, BEST_DELTA = CLIFFS_DELTA
    )],
    sym_map,
    by = "GENE",
    all.x = TRUE
  )
  gene_tab <- merge(gene_tab, gene_order, by = "GENE")
  setorder(gene_tab, GENE_ORD)
  fwrite(
    gene_tab,
    file.path(out_dir, sprintf("association_per_gene_top%d_genes.tsv", n_top)),
    sep = "\t"
  )
  message("Top-", n_top, " gene heatmap -> ", normalizePath(out_heat))
  message("Per-gene plots written")
}

message("Done -> ", normalizePath(out_dir))
