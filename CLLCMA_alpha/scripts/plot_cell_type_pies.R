#!/usr/bin/env Rscript
# Pie charts of cell_type_final: individual PDFs + polished 2-row panel figure.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
})

project_dir <- Sys.getenv(
  "CLL_TME_DIR",
  unset = path.expand("~/research/CLL_TME")
)
project_dir <- normalizePath(project_dir, mustWork = TRUE)
source(file.path(project_dir, "CLLCMA_alpha", "scripts", "helpers_cell_type_palette.R"))

in_file <- file.path(
  project_dir, "CLLCMA_alpha", "integrated_annotated", "cell_annotations_final.tsv"
)
anno_dir <- file.path(project_dir, "CLLCMA_alpha", "integrated_annotated")
out_dir <- file.path(anno_dir, "cell_type_pie_charts")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dt <- fread(in_file, select = c("dataset", "disease", "cell_type_final"))
counts <- dt[, .N, by = .(dataset, cell_type_final)]
counts <- counts[N > 0]

mbl_datasets <- unique(dt[disease %in% c("HC-MBL", "LC-MBL", "MBL"), dataset])
cll_datasets <- c(
  "GSE111014", "GSE132065", "GSE152469", "GSE161610",
  "GSE205704", "GSE261917", "GSE295489"
)
rs_datasets <- "Nadeu_2022"

dataset_order <- c("PBMCpedia_HC", mbl_datasets, cll_datasets, rs_datasets)
missing <- setdiff(unique(counts$dataset), dataset_order)
if (length(missing)) {
  stop("Unassigned datasets: ", paste(missing, collapse = ", "))
}

counts[, dataset := factor(dataset, levels = dataset_order)]
counts[, pct := N / sum(N), by = dataset]

all_types <- sort(unique(counts$cell_type_final))
pal <- cell_type_palette(all_types)
malig_types <- c("CLL", "RT", "MBL")

prepare_pie_data <- function(data) {
  data <- copy(data)
  setorder(data, -N)
  data[, cell_type_final := factor(
    cell_type_final,
    levels = rev(as.character(cell_type_final))
  )]
  data[, label := ""]
  data[, label_color := NA_character_]
  data[
    cell_type_final %in% malig_types,
    label := paste0(as.character(cell_type_final), "\n", scales::percent(pct, accuracy = 1))
  ]
  data[
    cell_type_final %in% malig_types,
    label_color := label_text_color(as.character(pal[as.character(cell_type_final)]))
  ]
  data[is.na(label_color), label_color := "#1A1A1A"]
  data
}

pie_plot <- function(
    data,
    label_size = 3.4,
    title_size = 13,
    subtitle_size = 9,
    base_size = 13,
    show_legend = TRUE,
    legend_ncol = 1,
    full_palette = FALSE) {
  ds <- as.character(data$dataset[1L])
  n_total <- sum(data$N)
  data <- prepare_pie_data(data)

  p <- ggplot(data, aes(x = 2, y = N, fill = cell_type_final)) +
    geom_col(width = 1, color = "white", linewidth = 0.45) +
    geom_text(
      aes(label = label, color = label_color),
      position = position_stack(vjust = 0.5),
      size = label_size,
      fontface = "bold",
      lineheight = 0.92,
      show.legend = FALSE
    ) +
    coord_polar(theta = "y", clip = "off", start = 0) +
    xlim(0.5, 2.5) +
    scale_fill_manual(
      values = pal,
      drop = !full_palette,
      name = "Cell type"
    ) +
    scale_color_identity() +
    labs(
      title = ds,
      subtitle = paste0("n = ", format(n_total, big.mark = ","), " cells")
    ) +
    theme_void(base_size = base_size) +
    theme(
      aspect.ratio = 1,
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(
        size = title_size, face = "bold", hjust = 0.5, margin = margin(b = 2)
      ),
      plot.subtitle = element_text(
        size = subtitle_size, face = "bold", hjust = 0.5,
        color = "#444444", margin = margin(b = 3)
      ),
      legend.text = element_text(size = 9, face = "bold"),
      legend.title = element_text(size = 11, face = "bold"),
      legend.key.height = unit(0.42, "cm"),
      legend.key.width = unit(0.42, "cm"),
      plot.margin = margin(6, 6, 6, 6)
    ) +
    guides(fill = guide_legend(ncol = legend_ncol, byrow = TRUE))

  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }

  p
}

make_shared_legend <- function(ncol = 8L, text_size = 8.5, title_size = 10) {
  dummy <- data.table(cell_type_final = factor(all_types, levels = all_types))
  p <- ggplot(dummy, aes(x = 1, y = 1, fill = cell_type_final)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_manual(values = pal, drop = FALSE, name = "Cell type") +
    guides(fill = guide_legend(ncol = ncol, byrow = TRUE)) +
    theme(
      legend.text = element_text(size = text_size, face = "bold"),
      legend.title = element_text(size = title_size, face = "bold"),
      legend.key.height = unit(0.42, "cm"),
      legend.key.width = unit(0.42, "cm"),
      legend.box = "horizontal",
      legend.position = "bottom",
      legend.background = element_rect(fill = "white", color = NA)
    ) +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank()
    )
  get_legend(p)
}

save_plot <- function(path, plot, width, height, dpi = 300) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "png") {
    ggsave(path, plot, width = width, height = height, dpi = dpi, bg = "white")
  } else if (capabilities("cairo") && ext == "pdf") {
    ggsave(path, plot, width = width, height = height, device = cairo_pdf, bg = "white")
  } else {
    ggsave(path, plot, width = width, height = height, bg = "white")
  }
}

# A4 landscape (297 mm x 210 mm)
A4_LANDSCAPE_W <- 11.69
A4_LANDSCAPE_H <- 8.27
PANEL_DPI <- 300L

for (ds in dataset_order) {
  sub <- counts[dataset == ds]
  if (!nrow(sub)) {
    next
  }
  out_file <- file.path(out_dir, paste0("cell_type_pie_", ds, ".pdf"))
  save_plot(out_file, pie_plot(sub, label_size = 3.6), width = 8, height = 8)
  message("Wrote: ", out_file)
}

row_plots <- lapply(dataset_order, function(ds) {
  sub <- counts[dataset == ds]
  pie_plot(
    sub,
    label_size = 3.4,
    title_size = 12,
    subtitle_size = 10,
    base_size = 14,
    show_legend = FALSE,
    full_palette = TRUE
  )
})

shared_legend <- make_shared_legend(ncol = 6L, text_size = 10, title_size = 12)
if (is.null(shared_legend)) {
  stop("Failed to build shared color legend")
}

panel_combined <- wrap_plots(row_plots, nrow = 2) /
  wrap_elements(shared_legend) +
  plot_layout(heights = c(1, 0.42))

panel_pdf <- file.path(anno_dir, "cell_type_pie_charts_row.pdf")
panel_png <- file.path(anno_dir, "cell_type_pie_charts_row.png")
save_plot(panel_pdf, panel_combined, A4_LANDSCAPE_W, A4_LANDSCAPE_H)
save_plot(panel_png, panel_combined, A4_LANDSCAPE_W, A4_LANDSCAPE_H, dpi = PANEL_DPI)
message("Wrote: ", panel_pdf)
message("Wrote: ", panel_png)
message("Done: ", length(dataset_order), " individual + A4 panel (PDF + PNG)")
