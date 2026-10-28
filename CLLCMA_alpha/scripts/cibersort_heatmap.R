#!/usr/bin/env Rscript

# Heatmaps of CLLCMA deconvolution cell fractions for the Broad CLL cohort
# (CIBERSORT and BayesPrism).
#
# Rows = RNA samples; columns = cell-type features. Values are column-wise
# z-scores (fractions vary widely across cell types; CLL dominates raw scale).
# Samples are split by IGHV status and clustered within each group.
#
# Inputs:
#   CLLCMA_alpha/results/deconvolution/cll_broad_2022/{cibersort,bayesprism}_sample_by_feature.tsv
#   CLLCMA_alpha/bulk/cll_broad_2022/data_clinical_{sample,patient}.txt
#
# Outputs (per method):
#   CLLCMA_alpha/results/deconvolution/cll_broad_2022/{cibersort,bayesprism}_heatmap/
#     {method}_heatmap_zscore.pdf
#     {method}_heatmap_raw.pdf
#     {method}_heatmap_sample_order.tsv
#     {method}_heatmap_matrix_zscore.tsv

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

project_dir <- Sys.getenv(
  "CLL_TME_DIR",
  unset = path.expand("~/research/CLL_TME")
)
project_dir <- normalizePath(project_dir, mustWork = TRUE)
setwd(project_dir)

deconv_dir <- file.path(
  "CLLCMA_alpha", "results", "deconvolution", "cll_broad_2022"
)
clinical_sample_file <- file.path(
  "CLLCMA_alpha", "bulk", "cll_broad_2022", "data_clinical_sample.txt"
)
clinical_patient_file <- file.path(
  "CLLCMA_alpha", "bulk", "cll_broad_2022", "data_clinical_patient.txt"
)

z_score_limit <- 3
ighv_row_order <- c("mutated", "unmutated", "unknown")

methods <- list(
  cibersort = list(
    input = file.path(deconv_dir, "cibersort_sample_by_feature.tsv"),
    output_dir = file.path(deconv_dir, "cibersort_heatmap"),
    title = "CIBERSORT fractions · Broad CLL (CLLCMA reference)",
    preferred_order = c(
      "CLL", "RS", "CD4_T", "CD8_T", "NK", "Monocyte", "DC", "ILC", "Progenitor"
    )
  ),
  bayesprism = list(
    input = file.path(deconv_dir, "bayesprism_sample_by_feature.tsv"),
    output_dir = file.path(deconv_dir, "bayesprism_heatmap"),
    title = "BayesPrism fractions · Broad CLL (CLLCMA reference)",
    preferred_order = c(
      "CLL", "RS",
      "CD4_T_naive", "CD4_T_memory",
      "CD8_T_naive", "CD8_T_memory", "CD8aa", "DN_T", "gdT",
      "Treg", "T_proliferating", "MAIT",
      "NK_CD56dim", "NK_CD56bright", "NK_proliferating",
      "Monocyte_CD14", "Monocyte_CD16", "Monocyte_intermediate",
      "cDC1", "cDC2", "ASDC",
      "ILC", "Progenitor"
    )
  )
)

# ------------------------------ Helpers --------------------------------------

read_cbioportal_clinical <- function(path) {
  fread(
    path,
    skip = 4L,
    na.strings = c("NA", "", "[Not Available]", "[Not Applicable]")
  )
}

safe_z_score <- function(x) {
  x <- as.numeric(x)
  spread <- sd(x, na.rm = TRUE)
  if (!is.finite(spread) || spread == 0) {
    return(rep(0, length(x)))
  }
  z <- (x - mean(x, na.rm = TRUE)) / spread
  z[!is.finite(z)] <- 0
  z
}

clean_annotation <- function(x, unknown = "Unknown") {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- unknown
  x
}

map_ighv_group <- function(x) {
  x <- tolower(clean_annotation(x, unknown = "unknown"))
  out <- rep("unknown", length(x))
  out[grepl("unmut", x)] <- "unmutated"
  out[grepl("^mut", x) | x == "mutated"] <- "mutated"
  out[x %in% c("m", "m-cll", "mutated ighv")] <- "mutated"
  out[x %in% c("u", "u-cll", "unmutated ighv")] <- "unmutated"
  factor(out, levels = ighv_row_order)
}

draw_heatmap <- function(
    mat,
    row_meta,
    col_fun,
    legend_title,
    file_path,
    column_title,
    cluster_columns = FALSE
) {
  ighv <- factor(row_meta$IGHV_GROUP, levels = ighv_row_order)
  sex <- factor(row_meta$SEX, levels = c("Female", "Male", "Unknown"))

  row_anno <- rowAnnotation(
    IGHV = ighv,
    Sex = sex,
    col = list(
      IGHV = c(
        mutated = "#2E7D4F",
        unmutated = "#C44E52",
        unknown = "#9E9E9E"
      ),
      Sex = c(
        Female = "#7A9E9F",
        Male = "#405D72",
        Unknown = "#B0B0B0"
      )
    ),
    annotation_name_side = "top",
    simple_anno_size = unit(3.5, "mm"),
    gap = unit(1, "mm")
  )

  ht <- Heatmap(
    mat,
    name = legend_title,
    col = col_fun,
    left_annotation = row_anno,
    cluster_rows = TRUE,
    cluster_columns = cluster_columns,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 9),
    row_split = ighv,
    cluster_row_slices = FALSE,
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 10),
    border = TRUE,
    use_raster = TRUE,
    raster_quality = 5,
    heatmap_legend_param = list(
      title = legend_title,
      direction = "vertical"
    ),
    column_title = column_title,
    column_title_gp = gpar(fontsize = 12, fontface = "bold")
  )

  pdf(file_path, width = 9.5, height = 11)
  drawn <- draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legend = TRUE
  )
  dev.off()

  row_idx <- unlist(row_order(drawn), use.names = FALSE)
  list(drawn = drawn, row_order = rownames(mat)[row_idx])
}

# --------------------------- Clinical metadata -------------------------------

clinical_sample <- read_cbioportal_clinical(clinical_sample_file)
clinical_patient <- read_cbioportal_clinical(clinical_patient_file)

meta <- merge(
  clinical_sample[, .(SAMPLE_ID, PATIENT_ID)],
  clinical_patient[, .(
    PATIENT_ID,
    SEX,
    IGHV_MUTATION_STATUS
  )],
  by = "PATIENT_ID",
  all.x = TRUE
)

# ------------------------------ Per method ----------------------------------

for (method in names(methods)) {
  cfg <- methods[[method]]
  if (!file.exists(cfg$input)) {
    warning("Skipping ", method, ": missing ", cfg$input)
    next
  }

  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
  message("=== ", method, " ===")

  cs <- fread(cfg$input)
  feature_cols <- setdiff(names(cs), "SAMPLE_ID")
  feature_cols <- c(
    intersect(cfg$preferred_order, feature_cols),
    setdiff(feature_cols, cfg$preferred_order)
  )

  dat <- merge(cs, meta, by = "SAMPLE_ID", all.x = TRUE)
  dat[, SEX := clean_annotation(SEX)]
  dat[, SEX := fifelse(SEX %in% c("Female", "Male"), SEX, "Unknown")]
  dat[, IGHV_GROUP := map_ighv_group(IGHV_MUTATION_STATUS)]
  setorder(dat, IGHV_GROUP, SAMPLE_ID)

  raw_mat <- as.matrix(dat[, ..feature_cols])
  storage.mode(raw_mat) <- "numeric"
  rownames(raw_mat) <- dat$SAMPLE_ID

  z_mat <- apply(raw_mat, 2L, safe_z_score)
  z_mat <- pmax(pmin(z_mat, z_score_limit), -z_score_limit)
  rownames(z_mat) <- rownames(raw_mat)

  row_meta <- dat[, .(SAMPLE_ID, PATIENT_ID, SEX, IGHV_GROUP)]

  message(
    "Heatmap matrix: ", nrow(raw_mat), " samples x ", ncol(raw_mat),
    " features"
  )
  print(table(row_meta$IGHV_GROUP))

  z_col <- colorRamp2(
    c(-z_score_limit, 0, z_score_limit),
    c("#2166AC", "#F7F7F7", "#B2182B")
  )
  raw_max <- quantile(raw_mat, 0.99, na.rm = TRUE)
  raw_col <- colorRamp2(
    c(0, raw_max / 2, raw_max),
    c("#F7FBFF", "#6BAED6", "#08306B")
  )

  z_out <- draw_heatmap(
    mat = z_mat,
    row_meta = row_meta,
    col_fun = z_col,
    legend_title = "z-score",
    file_path = file.path(cfg$output_dir, paste0(method, "_heatmap_zscore.pdf")),
    column_title = cfg$title,
    cluster_columns = FALSE
  )

  raw_out <- draw_heatmap(
    mat = raw_mat,
    row_meta = row_meta,
    col_fun = raw_col,
    legend_title = "fraction",
    file_path = file.path(cfg$output_dir, paste0(method, "_heatmap_raw.pdf")),
    column_title = cfg$title,
    cluster_columns = FALSE
  )

  fwrite(
    data.table(
      SAMPLE_ORDER = seq_along(z_out$row_order),
      SAMPLE_ID = z_out$row_order
    )[
      row_meta, on = "SAMPLE_ID"
    ],
    file.path(cfg$output_dir, paste0(method, "_heatmap_sample_order.tsv")),
    sep = "\t"
  )

  fwrite(
    as.data.table(z_mat, keep.rownames = "SAMPLE_ID"),
    file.path(cfg$output_dir, paste0(method, "_heatmap_matrix_zscore.tsv")),
    sep = "\t"
  )

  message("Wrote heatmaps to ", normalizePath(cfg$output_dir))
}

message("Done")
