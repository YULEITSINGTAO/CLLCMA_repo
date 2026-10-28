#!/usr/bin/env Rscript

# Three transition-focused volcano plots for B -> MBL -> CLL -> RT.
# Uses the existing DESeq2 result tables; no DEG values are recomputed.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

project_dir <- Sys.getenv("CLL_TME_DIR", unset = path.expand("~/research/CLL_TME"))
setwd(normalizePath(project_dir, mustWork = TRUE))

in_dir <- "CLLCMA_alpha/results/bcell_stage_deseq2"
out_dir <- file.path(in_dir, "publication_figures", "stage_progression")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

padj_cut <- 0.05
lfc_cut <- 1
x_cap <- 10
y_cap <- 50

contrasts <- data.table(
  stem = c("MBL_vs_B", "CLL_vs_MBL", "RT_vs_CLL"),
  title = c("MBL vs B", "CLL vs MBL", "RT vs CLL"),
  transition = c("B \u2192 MBL", "MBL \u2192 CLL", "CLL \u2192 RT")
)

prepare_volcano <- function(stem, title, transition) {
  tab <- fread(file.path(in_dir, paste0("deseq2_", stem, ".tsv")))
  tab <- tab[
    is.finite(log2FoldChange) & !is.na(padj)
  ]
  tab[, direction := fcase(
    padj < padj_cut & log2FoldChange >= lfc_cut, "Up",
    padj < padj_cut & log2FoldChange <= -lfc_cut, "Down",
    default = "Not significant"
  )]
  tab[, direction := factor(direction, levels = c("Up", "Down", "Not significant"))]
  tab[, x_plot := pmax(-x_cap, pmin(x_cap, log2FoldChange))]
  tab[, y_raw := -log10(pmax(padj, .Machine$double.xmin))]
  tab[, y_plot := pmin(y_raw, y_cap)]

  labels <- rbindlist(list(
    tab[direction == "Up"][order(padj, -abs(log2FoldChange))][seq_len(min(5L, .N))],
    tab[direction == "Down"][order(padj, -abs(log2FoldChange))][seq_len(min(5L, .N))]
  ))

  n_up <- tab[direction == "Up", .N]
  n_down <- tab[direction == "Down", .N]
  n_x_clip <- tab[abs(log2FoldChange) > x_cap, .N]
  n_y_clip <- tab[y_raw > y_cap, .N]

  p <- ggplot(tab, aes(x_plot, y_plot, color = direction, shape = direction)) +
    geom_point(size = 0.75, alpha = 0.65, stroke = 0.15) +
    geom_vline(
      xintercept = c(-lfc_cut, lfc_cut),
      linetype = 2, linewidth = 0.3, color = "grey35"
    ) +
    geom_hline(
      yintercept = -log10(padj_cut),
      linetype = 2, linewidth = 0.3, color = "grey35"
    ) +
    geom_text_repel(
      data = labels,
      aes(label = gene),
      size = 2.4,
      min.segment.length = 0,
      segment.size = 0.2,
      max.overlaps = Inf,
      box.padding = 0.25,
      seed = 42,
      fontface = "bold",
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(Up = "#D55E00", Down = "#0072B2", `Not significant` = "#B3B3B3"),
      labels = c(Up = "Up", Down = "Down", `Not significant` = "Not significant")
    ) +
    scale_shape_manual(
      values = c(Up = 24, Down = 25, `Not significant` = 16),
      labels = c(Up = "Up", Down = "Down", `Not significant` = "Not significant")
    ) +
    scale_x_continuous(
      limits = c(-x_cap, x_cap),
      breaks = seq(-10, 10, 5),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      limits = c(0, y_cap),
      breaks = seq(0, y_cap, 10),
      expand = expansion(mult = c(0, 0.03))
    ) +
    labs(
      title = title,
      subtitle = paste0(
        transition, "; padj<", padj_cut, ", |log2FC|\u2265", lfc_cut,
        "; Up=", n_up, ", Down=", n_down,
        "\nBoundary-capped: ", n_x_clip, " x values; ", n_y_clip, " y values"
      ),
      x = expression(log[2] * " fold change"),
      y = expression(-log[10] * " adjusted p-value"),
      color = NULL,
      shape = NULL
    ) +
    theme_classic(base_size = 9) +
    theme(
      text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.key.width = grid::unit(3, "mm"),
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 7.5, color = "grey25"),
      axis.title = element_text(size = 8.5)
    )

  summary <- data.table(
    contrast = stem,
    transition = transition,
    n_tested = nrow(tab),
    n_up = n_up,
    n_down = n_down,
    n_x_boundary_capped = n_x_clip,
    n_y_boundary_capped = n_y_clip
  )
  list(plot = p, summary = summary)
}

made <- lapply(seq_len(nrow(contrasts)), function(i) {
  prepare_volcano(contrasts$stem[i], contrasts$title[i], contrasts$transition[i])
})

summary_dt <- rbindlist(lapply(made, `[[`, "summary"))
fwrite(summary_dt, file.path(out_dir, "transition_volcano_summary.tsv"), sep = "\t")

combined <- wrap_plots(lapply(made, `[[`, "plot"), nrow = 3, guides = "collect") +
  plot_annotation(
    title = "Differential expression across ordered B-cell disease transitions",
    subtitle = "DESeq2 pseudobulk contrasts; positive log2FC indicates higher expression in the later stage",
    caption = paste(
      "Points at \u00b110 log2FC or 50 on the -log10(padj) axis exceed the displayed range and are boundary-capped.",
      "Stages are ordered disease states, not uniformly sampled longitudinal time.",
      "B-vs-disease contrasts are partly cohort-confounded; interpret MBL vs B as exploratory.",
      sep = "\n"
    ),
    theme = theme(
      text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(face = "bold", size = 9),
      plot.caption = element_text(face = "bold", size = 7, hjust = 0)
    )
  ) &
  theme(
    text = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  file.path(out_dir, "transition_volcanoes.pdf"),
  combined, width = 7.2, height = 15, device = cairo_pdf
)
ggsave(
  file.path(out_dir, "transition_volcanoes.png"),
  combined, width = 7.2, height = 15, dpi = 400, bg = "white"
)

writeLines(
  c(
    "Figure X. Differential expression across ordered B-cell disease transitions.",
    "Three volcano plots show existing pseudobulk DESeq2 results for MBL versus normal B, CLL versus MBL, and RT versus CLL.",
    "Positive log2 fold change indicates higher expression in the later disease stage.",
    "Genes are called significant at adjusted p<0.05 and |log2 fold change|\u22651.",
    "For readability, values outside \u00b110 log2 fold change and -log10(adjusted p)=50 are drawn at the corresponding boundary; counts are reported in transition_volcano_summary.tsv and each panel subtitle.",
    "Color and point shape redundantly encode direction. Labels show the five most significant up- and downregulated genes per contrast.",
    "The stages represent an ordered disease progression rather than uniformly sampled longitudinal time.",
    "The MBL-versus-B contrast is partly cohort-confounded and should be interpreted as exploratory."
  ),
  file.path(out_dir, "transition_volcanoes_caption.txt")
)

message("Wrote volcano figure to ", normalizePath(out_dir))
