#!/usr/bin/env Rscript

# Per-sample QC on the harmonized CLLCMA_alpha count matrix.
#
# Run AFTER step_1_harmonize_datasets.R and BEFORE step_2_annotate_integrated.R.
# Numbered 1.5 because it filters the merged object that step_1 wrote to
# integrated_raw/. Downstream Harmony / label transfer assume this matrix is
# already cleaned.
#
# Why QC here (not inside Seurat CreateSeuratObject):
#   The atlas is a concatenated sparse matrix from many cohorts. A single
#   global nCount/nFeature cutoff would over-filter shallow 10x libraries and
#   under-filter deeper ones. Thresholds that depend on library size (MAD
#   upper bounds) are therefore computed *within each sample*. Hard floors
#   (min genes / min UMIs / mito %) are still global so empty droplets and
#   dying cells are removed everywhere.
#
# Inputs:  CLLCMA_alpha/integrated_raw/
#   integrated_counts.rds          genes x cells CsparseMatrix
#   integrated_cell_metadata.tsv   one row per cell; must include sample, dataset
#   integrated_gene_metadata.tsv   feature names aligned to count rownames
#
# Outputs:
#   CLLCMA_alpha/integrated_raw/qc/
#     cell_qc_metrics.tsv          every cell, with pass/fail flags and reason
#     qc_summary_by_sample.tsv
#     qc_summary_by_dataset.tsv
#     qc_filter_params.tsv         thresholds used for this run
#     qc_sample_mad_bounds.tsv     per-sample upper cutoffs actually applied
#     plots/qc_violin_nFeature_by_dataset.pdf
#     plots/qc_violin_percent_mt_by_dataset.pdf
#     plots/qc_scatter_nCount_nFeature.pdf
#     plots/qc_pct_removed_by_sample.pdf
#   Filtered matrix is written back onto integrated_raw/ (counts + cell meta +
#   gene meta + summary). A pre-QC snapshot is archived as
#   integrated_raw_pre_qc_<timestamp>/ unless CLLCMA_QC_NO_ARCHIVE=TRUE.
#
# Filters (AND-ed; a cell must pass all of them):
#   - nFeature >= min_features (default 200)
#       Low gene counts usually mean empty droplets or ambient RNA.
#   - nCount   >= min_counts   (default 500)
#       Low UMI counts are the same failure mode, complementary to nFeature.
#   - percent.mt <= max_mt
#       Default 20% in peripheral blood, 25% in LN/BM. Lymphoid tissues often
#       run a bit higher mitochondrial fraction without being dying cells.
#   - nCount   <= median + nmads * MAD  (default nmads = 5)
#   - nFeature <= median + nmads * MAD
#       Upper tails catch doublets / multiplets and rare debris with inflated
#       complexity. MAD is robust to the heavy right tail of 10x libraries.
#       Samples with < min_cells_per_sample_for_mad cells skip the upper cut
#       so a handful of cells cannot define a nonsense bound.
#
# Env:
#   CLLCMA_QC_DRY_RUN=TRUE    metrics and plots only; do not rewrite integrated_raw
#   CLLCMA_QC_NO_ARCHIVE=TRUE overwrite integrated_raw without a pre-QC snapshot

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
})

project_dir <- Sys.getenv(
  "CLL_TME_DIR",
  unset = path.expand("~/research/CLL_TME")
)
project_dir <- normalizePath(project_dir, mustWork = TRUE)
setwd(project_dir)

atlas_dir <- "CLLCMA_alpha"
in_dir <- file.path(atlas_dir, "integrated_raw")
qc_dir <- file.path(in_dir, "qc")
plot_dir <- file.path(qc_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

dry_run <- identical(toupper(Sys.getenv("CLLCMA_QC_DRY_RUN", "")), "TRUE")
no_archive <- identical(toupper(Sys.getenv("CLLCMA_QC_NO_ARCHIVE", "")), "TRUE")

# ----------------------------- Parameters ------------------------------------
# These are the only numbers that change typical QC behavior. They are also
# written to qc_filter_params.tsv so a later run can be reproduced.

min_features <- 200L   # drop empty / ambient droplets
min_counts <- 500L     # drop ultra-low UMI barcodes
max_mt_pb <- 20        # PB dying-cell / stressed-cell ceiling (%)
max_mt_other <- 25     # slightly looser for LN / BM
n_mads_upper <- 5      # how far above the sample median is still "a singlet"
min_cells_per_sample_for_mad <- 50L  # skip MAD upper cut on tiny samples

params <- data.frame(
  param = c(
    "min_features", "min_counts", "max_mt_pb", "max_mt_other",
    "n_mads_upper", "min_cells_per_sample_for_mad", "mito_pattern"
  ),
  value = c(
    min_features, min_counts, max_mt_pb, max_mt_other,
    n_mads_upper, min_cells_per_sample_for_mad, "^MT-"
  ),
  stringsAsFactors = FALSE
)

# ----------------------------- Load ------------------------------------------
# Expect the step_1 layout: counts colnames == meta$cell_id, and gene names
# either already on the matrix rownames or in genes$feature_name.

counts_file <- file.path(in_dir, "integrated_counts.rds")
meta_file <- file.path(in_dir, "integrated_cell_metadata.tsv")
gene_file <- file.path(in_dir, "integrated_gene_metadata.tsv")

if (!file.exists(counts_file)) {
  stop("Missing ", counts_file, " — run step_1 first.")
}

message("Loading ", counts_file, " ...")
counts <- readRDS(counts_file)
# Column sums on a dgCMatrix are much faster than on triplet (Tsparse) form.
counts <- as(counts, "CsparseMatrix")
meta <- read.delim(meta_file, stringsAsFactors = FALSE, check.names = FALSE)
genes <- read.delim(gene_file, stringsAsFactors = FALSE, check.names = FALSE)

stopifnot(identical(colnames(counts), meta$cell_id))
if (!"sample" %in% names(meta)) {
  stop("integrated_cell_metadata.tsv needs a 'sample' column")
}
if (!"dataset" %in% names(meta)) {
  stop("integrated_cell_metadata.tsv needs a 'dataset' column")
}

# Prefer the human-readable gene symbol for MT- matching. Step_1 may have
# stored Ensembl IDs on the matrix and symbols in feature_name.
feature_names <- if ("feature_name" %in% names(genes)) {
  genes$feature_name
} else {
  rownames(counts)
}
if (!identical(rownames(counts), feature_names) &&
    length(feature_names) == nrow(counts)) {
  rownames(counts) <- feature_names
}

# Human mitochondrial genes are named MT-ND1, MT-CO1, ... (not mt- as in mouse).
mito_genes <- grep("^MT-", rownames(counts), value = TRUE)
message(
  "Cells=", ncol(counts), " genes=", nrow(counts),
  " samples=", length(unique(meta$sample)),
  " mito_genes=", length(mito_genes)
)
if (!length(mito_genes)) {
  warning("No MT- genes found; percent.mt will be 0")
}

# ----------------------------- Metrics ---------------------------------------
# Three standard scRNA-seq QC metrics, computed once for every barcode:
#   nCount_RNA   total UMIs (library size)
#   nFeature_RNA number of detected genes (complexity)
#   percent_mt   mitochondrial UMI share (stress / dying cells)

message("Computing per-cell QC metrics ...")
nCount <- Matrix::colSums(counts)
nFeature <- Matrix::colSums(counts > 0)
if (length(mito_genes)) {
  # pmax(..., 1) avoids 0/0 in empty barcodes; those fail min_counts anyway.
  percent_mt <- 100 * Matrix::colSums(counts[mito_genes, , drop = FALSE]) / pmax(nCount, 1)
} else {
  percent_mt <- rep(0, ncol(counts))
}

qc <- data.frame(
  cell_id = meta$cell_id,
  dataset = meta$dataset,
  sample = meta$sample,
  donor = if ("donor" %in% names(meta)) meta$donor else NA_character_,
  disease = if ("disease" %in% names(meta)) meta$disease else NA_character_,
  tissue = if ("tissue" %in% names(meta)) meta$tissue else NA_character_,
  nCount_RNA = as.numeric(nCount),
  nFeature_RNA = as.integer(nFeature),
  percent_mt = as.numeric(percent_mt),
  stringsAsFactors = FALSE
)

# Tissue-aware mito ceiling: PB keeps the stricter 20%; LN/BM allow 25%.
# Unlabeled / other tissues inherit the PB cutoff (safer default).
tissue_u <- toupper(trimws(as.character(qc$tissue)))
qc$max_mt_allowed <- ifelse(
  tissue_u %in% c("LN", "BM", "LYMPH_NODE", "BONE_MARROW"),
  max_mt_other,
  max_mt_pb
)

# Median absolute deviation of a vector, with two fallbacks:
#   1) mad(..., constant = 1) is the raw MAD (no 1.4826 Gaussian scaling).
#   2) If that is 0 (many identical values), try the scaled MAD.
#   3) If still 0 / non-finite, return NA so the caller can use a 98th
#      percentile instead of a zero-width bound that would keep everything.
mad_safe <- function(x) {
  m <- mad(x, constant = 1, na.rm = TRUE)
  if (!is.finite(m) || m == 0) {
    m <- mad(x, na.rm = TRUE)
  }
  if (!is.finite(m) || m == 0) {
    # fall back to 98th percentile span
    return(NA_real_)
  }
  m
}

# One upper bound per sample. Using the sample (not dataset or whole atlas)
# keeps 3' 10x v2 libraries from being judged against v3 / 5' depth.
sample_bounds <- lapply(split(seq_len(nrow(qc)), qc$sample), function(ii) {
  n <- length(ii)
  nf <- qc$nFeature_RNA[ii]
  nc <- qc$nCount_RNA[ii]
  if (n < min_cells_per_sample_for_mad) {
    # Tiny samples: do not invent an upper cutoff from <50 cells.
    return(data.frame(
      sample = qc$sample[ii[1]],
      n_cells_pre = n,
      use_mad = FALSE,
      max_nFeature = Inf,
      max_nCount = Inf,
      stringsAsFactors = FALSE
    ))
  }
  med_f <- median(nf, na.rm = TRUE)
  med_c <- median(nc, na.rm = TRUE)
  mad_f <- mad_safe(nf)
  mad_c <- mad_safe(nc)
  max_f <- if (is.finite(mad_f)) med_f + n_mads_upper * mad_f else as.numeric(quantile(nf, 0.98, na.rm = TRUE))
  max_c <- if (is.finite(mad_c)) med_c + n_mads_upper * mad_c else as.numeric(quantile(nc, 0.98, na.rm = TRUE))
  data.frame(
    sample = qc$sample[ii[1]],
    n_cells_pre = n,
    use_mad = TRUE,
    max_nFeature = max_f,
    max_nCount = max_c,
    stringsAsFactors = FALSE
  )
})
bounds <- do.call(rbind, sample_bounds)
rownames(bounds) <- NULL

qc <- merge(qc, bounds[, c("sample", "use_mad", "max_nFeature", "max_nCount")],
            by = "sample", all.x = TRUE, sort = FALSE)
# merge() reorders rows; put cells back in matrix column order before slicing.
qc <- qc[match(meta$cell_id, qc$cell_id), ]
stopifnot(identical(qc$cell_id, meta$cell_id))

qc$pass_min_features <- qc$nFeature_RNA >= min_features
qc$pass_min_counts <- qc$nCount_RNA >= min_counts
qc$pass_mt <- qc$percent_mt <= qc$max_mt_allowed
qc$pass_max_features <- qc$nFeature_RNA <= qc$max_nFeature
qc$pass_max_counts <- qc$nCount_RNA <= qc$max_nCount
qc$pass_qc <- (
  qc$pass_min_features &
    qc$pass_min_counts &
    qc$pass_mt &
    qc$pass_max_features &
    qc$pass_max_counts
)

# Single fail_reason for summaries. Order is the typical failure cascade:
# empty droplet (low genes) -> low UMIs -> dying (mito) -> doublet (high).
# Cells that fail several checks keep only the first matching reason.
qc$fail_reason <- NA_character_
qc$fail_reason[!qc$pass_min_features] <- "low_nFeature"
qc$fail_reason[!qc$pass_min_counts & is.na(qc$fail_reason)] <- "low_nCount"
qc$fail_reason[!qc$pass_mt & is.na(qc$fail_reason)] <- "high_percent_mt"
qc$fail_reason[!qc$pass_max_features & is.na(qc$fail_reason)] <- "high_nFeature"
qc$fail_reason[!qc$pass_max_counts & is.na(qc$fail_reason)] <- "high_nCount"

# ----------------------------- Summaries -------------------------------------
# Sample-level table is the one to inspect for a library that lost too many
# cells (e.g. a failed 10x channel). Dataset-level is the atlas overview.

summ_sample <- do.call(rbind, lapply(split(qc, list(qc$dataset, qc$sample), drop = TRUE), function(d) {
  data.frame(
    dataset = d$dataset[1],
    sample = d$sample[1],
    disease = d$disease[1],
    tissue = d$tissue[1],
    n_cells_pre = nrow(d),
    n_cells_post = sum(d$pass_qc),
    n_removed = sum(!d$pass_qc),
    pct_removed = 100 * mean(!d$pass_qc),
    median_nFeature_pre = median(d$nFeature_RNA),
    median_nFeature_post = if (any(d$pass_qc)) median(d$nFeature_RNA[d$pass_qc]) else NA_real_,
    median_nCount_pre = median(d$nCount_RNA),
    median_percent_mt_pre = median(d$percent_mt),
    median_percent_mt_post = if (any(d$pass_qc)) median(d$percent_mt[d$pass_qc]) else NA_real_,
    max_mt_allowed = d$max_mt_allowed[1],
    use_mad = d$use_mad[1],
    stringsAsFactors = FALSE
  )
}))
rownames(summ_sample) <- NULL
summ_sample <- summ_sample[order(summ_sample$dataset, -summ_sample$n_cells_pre), ]

summ_dataset <- do.call(rbind, lapply(split(qc, qc$dataset), function(d) {
  data.frame(
    dataset = d$dataset[1],
    n_samples = length(unique(d$sample)),
    n_cells_pre = nrow(d),
    n_cells_post = sum(d$pass_qc),
    n_removed = sum(!d$pass_qc),
    pct_removed = 100 * mean(!d$pass_qc),
    stringsAsFactors = FALSE
  )
}))
rownames(summ_dataset) <- NULL
summ_dataset <- summ_dataset[order(-summ_dataset$n_cells_pre), ]

message("QC summary by dataset:")
print(summ_dataset)
message(
  "Total kept: ", sum(qc$pass_qc), " / ", nrow(qc),
  " (", sprintf("%.1f", 100 * mean(qc$pass_qc)), "%)"
)

write.table(params, file.path(qc_dir, "qc_filter_params.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(qc, file.path(qc_dir, "cell_qc_metrics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summ_sample, file.path(qc_dir, "qc_summary_by_sample.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summ_dataset, file.path(qc_dir, "qc_summary_by_dataset.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bounds, file.path(qc_dir, "qc_sample_mad_bounds.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ----------------------------- Plots -----------------------------------------
# Violins use all cells (pass vs fail fill). The nCount–nFeature scatter is
# subsampled so a ~hundred-thousand-cell atlas does not produce a multi-GB PDF.

message("Writing QC plots ...")
# subsample points for scatter readability
set.seed(1L)
plot_qc <- qc
if (nrow(plot_qc) > 150000L) {
  plot_qc <- plot_qc[sample.int(nrow(plot_qc), 150000L), ]
}

p1 <- ggplot(qc, aes(x = dataset, y = nFeature_RNA, fill = pass_qc)) +
  geom_violin(scale = "width", linewidth = 0.2) +
  scale_y_log10() +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(title = "nFeature_RNA by dataset", y = "nFeature_RNA (log10)")
ggsave(file.path(plot_dir, "qc_violin_nFeature_by_dataset.pdf"), p1, width = 11, height = 5)

p2 <- ggplot(qc, aes(x = dataset, y = percent_mt, fill = pass_qc)) +
  geom_violin(scale = "width", linewidth = 0.2) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(title = "percent.mt by dataset")
ggsave(file.path(plot_dir, "qc_violin_percent_mt_by_dataset.pdf"), p2, width = 11, height = 5)

p3 <- ggplot(plot_qc, aes(x = nCount_RNA, y = nFeature_RNA, color = pass_qc)) +
  geom_point(size = 0.15, alpha = 0.35) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~dataset) +
  theme_bw(base_size = 10) +
  labs(title = "nCount vs nFeature (subsampled points)")
ggsave(file.path(plot_dir, "qc_scatter_nCount_nFeature.pdf"), p3, width = 12, height = 9)

p4 <- ggplot(summ_sample, aes(x = reorder(sample, -pct_removed), y = pct_removed, fill = dataset)) +
  geom_col(width = 0.8) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(title = "% cells removed per sample", x = "sample", y = "% removed")
ggsave(file.path(plot_dir, "qc_pct_removed_by_sample.pdf"), p4, width = 12, height = 4.5)

# ----------------------------- Filter + write --------------------------------
# Default path rewrites integrated_raw in place so step_2 never has to know
# about QC. Dry-run stops here. Otherwise we snapshot the pre-QC files first
# (hardlink when the filesystem allows) so the unfiltered matrix can be
# recovered without re-running step_1.

keep <- qc$pass_qc
if (!any(keep)) {
  stop("All cells failed QC — check thresholds")
}

if (dry_run) {
  message("CLLCMA_QC_DRY_RUN=TRUE — metrics/plots written; integrated_raw NOT modified.")
  message("QC reports in: ", normalizePath(qc_dir))
  quit(save = "no", status = 0)
}

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
if (!no_archive) {
  archive_dir <- file.path(atlas_dir, paste0("integrated_raw_pre_qc_", stamp))
  message("Archiving pre-QC integrated_raw -> ", archive_dir)
  dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
  for (f in c(
    "integrated_counts.rds",
    "integrated_cell_metadata.tsv",
    "integrated_gene_metadata.tsv",
    "integrated_summary.tsv"
  )) {
    src <- file.path(in_dir, f)
    if (file.exists(src)) {
      # hardlink when possible to save space; fall back to copy
      ok <- suppressWarnings(file.link(src, file.path(archive_dir, f)))
      if (!isTRUE(ok) && !file.exists(file.path(archive_dir, f))) {
        file.copy(src, file.path(archive_dir, f), overwrite = FALSE)
      }
    }
  }
  writeLines(
    c(
      paste0("Pre-QC snapshot of integrated_raw at ", stamp),
      paste0("n_cells_pre=", nrow(qc)),
      paste0("n_cells_post=", sum(keep)),
      "Created by step_1.5_qc_per_sample.R"
    ),
    file.path(archive_dir, "README.txt")
  )
}

message("Writing filtered integrated_raw (", sum(keep), " cells) ...")
counts_f <- counts[, keep, drop = FALSE]
# Genes with zero counts in every remaining cell cannot affect PCA / Harmony
# and would only inflate the feature space.
keep_genes <- Matrix::rowSums(counts_f) > 0
counts_f <- counts_f[keep_genes, , drop = FALSE]
meta_f <- meta[keep, , drop = FALSE]
# attach QC metrics used for filtering
meta_f$nCount_RNA <- qc$nCount_RNA[keep]
meta_f$nFeature_RNA <- qc$nFeature_RNA[keep]
meta_f$percent_mt <- qc$percent_mt[keep]

genes_f <- data.frame(feature_name = rownames(counts_f), stringsAsFactors = FALSE)

summary_f <- as.data.frame(table(
  dataset = meta_f$dataset,
  sample = meta_f$sample,
  disease = meta_f$disease,
  tissue = meta_f$tissue
))
names(summary_f)[5] <- "n_cells"
summary_f <- summary_f[summary_f$n_cells > 0, ]
summary_f <- summary_f[order(summary_f$dataset, summary_f$sample), ]

saveRDS(counts_f, counts_file, compress = FALSE)
write.table(meta_f, meta_file, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(genes_f, gene_file, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(
  summary_f,
  file.path(in_dir, "integrated_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Keep a pointer in VERSIONS-style note
ver_note <- file.path(atlas_dir, "VERSIONS.md")
extra <- paste0(
  "\n## pre_qc_", stamp, "\n\n",
  "Snapshot before step_1.5 QC: `integrated_raw_pre_qc_", stamp, "/`\n",
  "- cells pre/post: ", nrow(qc), " / ", sum(keep), "\n"
)
if (file.exists(ver_note)) {
  cat(extra, file = ver_note, append = TRUE)
}

message(
  "Done. Kept ", sum(keep), " / ", nrow(qc),
  " cells (", nrow(counts_f), " genes). QC dir: ", normalizePath(qc_dir)
)
if (!no_archive) {
  message("Pre-QC archive: ", normalizePath(file.path(atlas_dir, paste0("integrated_raw_pre_qc_", stamp))))
}
