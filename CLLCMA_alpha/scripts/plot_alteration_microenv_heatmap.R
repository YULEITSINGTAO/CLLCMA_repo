#!/usr/bin/env Rscript

# Publication polish for alteration × microenvironment heatmap + n-sig bar.
#
# Journal / article type / submission phase: not specified.
# This is a provisional general figure; publisher-specific rules remain pending.
#
# Inputs:
#   alteration_microenv_association/<method>/association_all.tsv
# Outputs:
#   alteration_microenv_heatmap.pdf
#   alteration_microenv_heatmap.png
#   alteration_nsig_celltypes.tsv

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
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

fdr_star <- 0.05
fdr_dot <- 0.10

out_dir <- file.path(
  "CLLCMA_alpha", "results", "deconvolution", "cll_broad_2022",
  "alteration_microenv_association", deconv_method
)
assoc_file <- file.path(out_dir, "association_all.tsv")
if (!file.exists(assoc_file)) {
  stop("Missing: ", assoc_file)
}

assoc <- fread(assoc_file)

# Logical cell-type order for reading
cell_order <- c(
  "CLL",
  "Progenitor",
  "CD4_T_naive", "CD4_T_memory", "Treg",
  "CD8_T_naive", "CD8_T_memory", "CD8aa", "DN_T", "gdT", "MAIT",
  "T_proliferating",
  "NK_CD56bright", "NK_CD56dim", "NK_proliferating",
  "ILC",
  "Monocyte_CD14", "Monocyte_intermediate", "Monocyte_CD16",
  "cDC1", "cDC2", "ASDC"
)
cell_order <- cell_order[cell_order %in% unique(assoc$CELL_TYPE)]
extra_ct <- setdiff(unique(assoc$CELL_TYPE), cell_order)
cell_order <- c(cell_order, sort(extra_ct))

cell_lab <- c(
  CLL = "CLL",
  Progenitor = "Progenitor",
  CD4_T_naive = "CD4 naive",
  CD4_T_memory = "CD4 memory",
  Treg = "Treg",
  CD8_T_naive = "CD8 naive",
  CD8_T_memory = "CD8 memory",
  CD8aa = "CD8aa",
  DN_T = "DN T",
  gdT = "gdT",
  MAIT = "MAIT",
  T_proliferating = "T prolif.",
  NK_CD56bright = "NK bright",
  NK_CD56dim = "NK dim",
  NK_proliferating = "NK prolif.",
  ILC = "ILC",
  Monocyte_CD14 = "Mono CD14",
  Monocyte_intermediate = "Mono int.",
  Monocyte_CD16 = "Mono CD16",
  cDC1 = "cDC1",
  cDC2 = "cDC2",
  ASDC = "ASDC"
)

pretty_cnv_label <- function(x) {
  x <- as.character(x)
  x <- sub("^arm_", "", x)
  x <- gsub("tri\\(12\\)", "+12", x)
  x <- sub("_altered$", "", x)
  x <- sub("_Gain$", " gain", x)
  x <- sub("_Loss$", " loss", x)
  x
}

heat <- assoc[
  (ALTERATION_CLASS == "mutation") |
    (ALTERATION_CLASS == "arm_cna" & grepl("_altered$", ALTERATION)) |
    (ALTERATION_CLASS == "focal_cna" & grepl("_altered$", ALTERATION))
]
heat[, CELL_TYPE := factor(CELL_TYPE, levels = cell_order)]
heat[, CELL_LABEL := factor(
  fifelse(
    as.character(CELL_TYPE) %in% names(cell_lab),
    unname(cell_lab[as.character(CELL_TYPE)]),
    as.character(CELL_TYPE)
  ),
  levels = unname(cell_lab[cell_order])
)]
heat[, ALTERATION_TYPE := factor(
  fcase(
    ALTERATION_CLASS == "mutation", "Mutation",
    ALTERATION_CLASS == "arm_cna", "Arm CNV",
    ALTERATION_CLASS == "focal_cna", "Focal CNV"
  ),
  levels = c("Mutation", "Arm CNV", "Focal CNV")
)]
heat[, ALTERATION_LABEL := fifelse(
  ALTERATION_CLASS == "mutation",
  as.character(ALTERATION),
  pretty_cnv_label(ALTERATION)
)]
heat[, SIG_MARK := fcase(
  FDR < fdr_star, "*",
  FDR < fdr_dot, "\u00b7",
  default = ""
)]

# Row order: within class, by number of FDR<0.05 hits then min FDR
row_rank <- heat[
  ,
  .(
    n_sig = sum(FDR < fdr_star, na.rm = TRUE),
    min_fdr = min(FDR, na.rm = TRUE),
    max_abs = max(abs(CLIFFS_DELTA), na.rm = TRUE)
  ),
  by = .(ALTERATION_TYPE, ALTERATION_LABEL, ALTERATION_CLASS, ALTERATION)
][order(ALTERATION_TYPE, -n_sig, min_fdr, -max_abs)]
row_levels <- rev(unique(as.character(row_rank$ALTERATION_LABEL)))
heat[, ALTERATION_LABEL := factor(ALTERATION_LABEL, levels = row_levels)]
row_rank[, ALTERATION_LABEL := factor(ALTERATION_LABEL, levels = row_levels)]

fwrite(
  row_rank[, .(
    ALTERATION_TYPE, ALTERATION_LABEL, ALTERATION_CLASS, ALTERATION,
    N_SIG_CELLTYPES = n_sig, MIN_FDR = min_fdr
  )][order(ALTERATION_TYPE, -N_SIG_CELLTYPES)],
  file.path(out_dir, "alteration_nsig_celltypes.tsv"),
  sep = "\t"
)

n_rows <- length(row_levels)
n_ct <- length(cell_order)
delta_max <- max(abs(heat$CLIFFS_DELTA), na.rm = TRUE)
delta_lim <- min(0.85, max(0.5, ceiling(delta_max * 20) / 20))

# Okabe-Ito-friendly category fills (high contrast on white)
type_fill <- c(
  Mutation = "#0072B2",
  `Arm CNV` = "#D55E00",
  `Focal CNV` = "#009E73"
)

# Diverging map centered at 0 (blue–gray–vermillion; color + * / · marks)
div_low <- "#0072B2"
div_mid <- "#F7F7F7"
div_high <- "#D55E00"

theme_pub <- function(base_size = 9) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.5, color = "#333333", hjust = 0),
      axis.title = element_text(size = base_size),
      axis.text = element_text(color = "#111111"),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "#666666", fill = NA, linewidth = 0.4),
      strip.background = element_rect(fill = "#ECECEC", color = NA),
      legend.title = element_text(size = base_size - 0.5),
      legend.text = element_text(size = base_size - 1),
      plot.margin = margin(4, 4, 4, 4)
    )
}

p_heat <- ggplot(
  heat,
  aes(x = CELL_LABEL, y = ALTERATION_LABEL, fill = CLIFFS_DELTA)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(
    aes(label = SIG_MARK),
    size = 2.6,
    color = "#111111",
    vjust = 0.55,
    fontface = "bold"
  ) +
  facet_grid(
    ALTERATION_TYPE ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_fill_gradient2(
    low = div_low,
    mid = div_mid,
    high = div_high,
    midpoint = 0,
    limits = c(-delta_lim, delta_lim),
    oob = squish,
    breaks = seq(-delta_lim, delta_lim, length.out = 5),
    labels = function(x) sprintf("%.2f", x),
    name = "Cliff's \u03b4\n(altered vs WT)"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    title = "Genomic alterations associated with circulating microenvironment composition",
    subtitle = sprintf(
      "%s fractions; Cliff's \u03b4 (descriptive); Storey adaptive FDR from IGHV/sex/age-adjusted LM within class. * FDR < %.2f; \u00b7 FDR < %.2f",
      method_label, fdr_star, fdr_dot
    ),
    x = NULL,
    y = NULL
  ) +
  theme_pub(9) +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 7.5),
    axis.text.y = element_text(size = 7.5, face = "plain"),
    strip.placement = "outside",
    strip.text.y.left = element_text(
      angle = 0, face = "bold", size = 8, hjust = 1
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = unit(1.1, "cm"),
    legend.key.height = unit(0.28, "cm"),
    plot.title.position = "plot"
  )

p_bar <- ggplot(
  row_rank,
  aes(x = n_sig, y = ALTERATION_LABEL, fill = ALTERATION_TYPE)
) +
  geom_col(width = 0.72, color = NA) +
  geom_text(
    data = row_rank[n_sig > 0],
    aes(label = n_sig),
    hjust = -0.15,
    size = 2.2,
    color = "#222222"
  ) +
  facet_grid(
    ALTERATION_TYPE ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_fill_manual(values = type_fill, guide = "none") +
  scale_x_continuous(
    limits = c(0, n_ct),
    breaks = pretty_breaks(n = 4),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    x = sprintf("No. cell types\n(FDR < %.2f)", fdr_star),
    y = NULL
  ) +
  theme_pub(9) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_blank(),
    strip.background = element_blank(),
    plot.margin = margin(4, 8, 4, 0),
    # leave room so x label aligns under heatmap x labels
    axis.title.x = element_text(size = 8, margin = margin(t = 2))
  )

# Match bottom space so facet panels align (heatmap has angled x labels)
p_heat_aligned <- p_heat + theme(
  plot.margin = margin(4, 2, 4, 4)
)
p_bar_aligned <- p_bar + theme(
  plot.margin = margin(4, 10, 28, 0)
)

fig <- p_heat_aligned + p_bar_aligned +
  plot_layout(widths = c(4.6, 1.05), guides = "collect") &
  theme(legend.position = "bottom")

fig_h <- max(8.2, 0.22 * n_rows + 3.4)
fig_w <- 12.2

out_pdf <- file.path(out_dir, "alteration_microenv_heatmap.pdf")
out_png <- file.path(out_dir, "alteration_microenv_heatmap.png")
ggsave(out_pdf, fig, width = fig_w, height = fig_h, device = cairo_pdf)
ggsave(out_png, fig, width = fig_w, height = fig_h, dpi = 400, bg = "white")

# Small provenance sidecar
writeLines(
  c(
    "Provisional general figure (journal/phase not specified).",
    paste0("Method: ", method_label),
    paste0("Source: ", normalizePath(assoc_file)),
    "Encoding: Cliff's delta (altered vs WT, descriptive); * adaptive FDR<0.05; · adaptive FDR<0.10",
    "FDR from LM: fraction ~ alteration + IGHV + sex + age/10 (complete cases)",
    "Right bars: count of cell types with adaptive FDR<0.05 per alteration",
    "Diverging colors: Okabe-Ito blue (#0072B2) / vermillion (#D55E00); midpoint 0",
    "Category bars: Okabe-Ito on-white subset (blue / vermillion / bluish green)",
    paste0("Exported: ", Sys.time())
  ),
  file.path(out_dir, "alteration_microenv_heatmap_provenance.txt")
)

message("Wrote ", normalizePath(out_pdf))
message("Wrote ", normalizePath(out_png))
message("Rows: ", n_rows, " | cell types: ", n_ct, " | delta limit: \u00b1", delta_lim)
