#!/usr/bin/env Rscript
# Bar plots of cell counts from integrated label-transfer summary.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

project_dir <- Sys.getenv(
  "CLL_TME_DIR",
  unset = path.expand("~/research/CLL_TME")
)
project_dir <- normalizePath(project_dir, mustWork = TRUE)

in_file <- file.path(
  project_dir, "CLLCMA_alpha", "integrated_annotated", "transfer_summary.tsv"
)
out_dir <- dirname(in_file)

dt <- fread(in_file)
setnames(dt, c("dataset", "cell_type_final", "n_cells"))

tot <- dt[, .(total = sum(n_cells)), by = dataset]
setorder(tot, -total)
dt[, dataset := factor(dataset, levels = rev(tot$dataset))]

dt[, cell_type_final := reorder(cell_type_final, n_cells)]

p_total <- ggplot(tot, aes(reorder(dataset, total), total)) +
  geom_col(fill = "#4575b4", width = 0.75) +
  geom_text(aes(label = scales::comma(total)), hjust = -0.05, size = 4.5, fontface = "bold") +
  coord_flip(clip = "off") +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.12))) +
  labs(
    x = NULL,
    y = "Total cells",
    title = "Integrated reference + query datasets",
    subtitle = "Cell counts after Harmony integration and label transfer"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12, face = "bold"),
    axis.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(5.5, 40, 5.5, 5.5)
  )

p_types <- ggplot(dt, aes(n_cells, cell_type_final)) +
  geom_col(fill = "#74add1", width = 0.75) +
  geom_text(aes(label = scales::comma(n_cells)), hjust = -0.05, size = 3.5, fontface = "bold") +
  scale_x_continuous(
    labels = scales::comma,
    trans = "log10",
    expand = expansion(mult = c(0, 0.15))
  ) +
  facet_wrap(~ dataset, scales = "free_y", nrow = 2) +
  labs(
    x = "Cell count (log10 scale)",
    y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    panel.spacing = unit(1.0, "lines")
  )

ggsave(
  file.path(out_dir, "transfer_summary_total_cells.pdf"),
  p_total,
  width = 8,
  height = 5
)

ggsave(
  file.path(out_dir, "transfer_summary_cell_counts.pdf"),
  p_types,
  width = 24,
  height = 12
)

message("Wrote: ", file.path(out_dir, "transfer_summary_total_cells.pdf"))
message("Wrote: ", file.path(out_dir, "transfer_summary_cell_counts.pdf"))
