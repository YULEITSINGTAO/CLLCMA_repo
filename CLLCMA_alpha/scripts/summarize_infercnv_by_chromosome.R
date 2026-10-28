#!/usr/bin/env Rscript

# Donor-level chromosome summaries from the final inferCNV matrix.
# HMM was not run; this does not call discrete CNVs.
#
# Estimator: for each cell, chromosome-mean log2(inferCNV expr.data),
# then the donor-stage mean of those cell values (unit of replication = donor).
# Focal companion: 80-gene genomic bins, then the bin with largest |mean|
# within each donor x chromosome.

suppressPackageStartupMessages({
  library(data.table)
})

project_dir <- Sys.getenv("CLL_TME_DIR", unset = path.expand("~/research/CLL_TME"))
setwd(normalizePath(project_dir, mustWork = TRUE))

out_dir <- "CLLCMA_alpha/results/bcell_stage_deseq2/publication_figures/infercnv_bcell_stage"
fig_dir <- file.path(out_dir, "chromosome_summary")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

stage_order <- c("Normal_B", "MBL", "CLL", "RT")
chr_levels <- paste0("chr", 1:22)
arm_levels <- as.vector(rbind(paste0(1:22, "p"), paste0(1:22, "q")))
min_cells_per_donor <- 5L
min_genes_per_arm <- 10L
bin_size <- 80L

obj <- readRDS(file.path(out_dir, "infercnv_obj_final.rds"))
expr <- obj@expr.data

anno <- fread(file.path(out_dir, "infercnv_cell_annotations.tsv"), header = FALSE)
setnames(anno, c("cell_id", "stage"))
anno <- anno[cell_id %in% colnames(expr)]
anno[, stage := factor(stage, levels = stage_order)]

meta <- fread(file.path(out_dir, "infercnv_selected_cells_metadata.tsv"))
id_cols <- names(meta)[names(meta) %in% c("cell_id", "cell_id.1", "cell_id_raw")]
match_col <- NULL
for (col in id_cols) {
  if (any(meta[[col]] %in% colnames(expr))) {
    match_col <- col
    break
  }
}
if (is.null(match_col)) {
  stop("Could not match metadata cell IDs to inferCNV colnames")
}
meta[, cell_id := as.character(get(match_col))]
meta <- unique(meta[, .(cell_id, donor, sample, dataset, infercnv_group, b_stage)])
meta <- meta[cell_id %in% colnames(expr)]
anno <- merge(anno, meta[, .(cell_id, donor, sample, dataset)], by = "cell_id", all.x = TRUE)
if (anyNA(anno$donor)) {
  stop("Missing donor labels for ", sum(is.na(anno$donor)), " cells")
}
# Donor labels are not globally unique across source datasets. Qualify them so
# unrelated donors with reused labels remain separate biological replicates.
anno[, donor_original := donor]
anno[, donor := paste(dataset, donor_original, sep = "::")]

gene_order <- fread(file.path(out_dir, "infercnv_gene_order.tsv"), header = FALSE)
setnames(gene_order, c("gene", "chr", "start", "end"))
gene_order <- gene_order[gene %in% rownames(expr) & chr %in% chr_levels]
gene_order[, chr := factor(chr, levels = chr_levels)]
setorder(gene_order, chr, start, end, gene)
expr <- expr[gene_order$gene, anno$cell_id, drop = FALSE]

# inferCNV final expr.data is fold-change-like and centered near 1.
log_dev <- log2(pmax(expr, 1e-6))

chr_n_genes <- gene_order[, .N, by = chr]
fwrite(chr_n_genes, file.path(fig_dir, "chromosome_gene_counts.tsv"), sep = "\t")

chr_cell <- sapply(chr_levels, function(ch) {
  genes <- gene_order[as.character(chr) == ch, gene]
  colMeans(log_dev[genes, , drop = FALSE], na.rm = TRUE)
})
chr_cell <- as.data.table(chr_cell, keep.rownames = "cell_id")
chr_cell <- melt(
  chr_cell,
  id.vars = "cell_id",
  variable.name = "chr",
  value.name = "log2_dev",
  variable.factor = FALSE
)
chr_cell <- merge(chr_cell, anno, by = "cell_id")
chr_cell[, chr := factor(chr, levels = chr_levels)]

donor_chr <- chr_cell[, .(
  n_cells = .N,
  mean_log2_dev = mean(log2_dev, na.rm = TRUE),
  sd_cells = sd(log2_dev, na.rm = TRUE),
  median_log2_dev = median(log2_dev, na.rm = TRUE)
), by = .(donor, stage, dataset, chr)]
donor_chr <- donor_chr[n_cells >= min_cells_per_donor]
donor_chr[, chr := factor(chr, levels = chr_levels)]
setorder(donor_chr, stage, donor, chr)
fwrite(donor_chr, file.path(fig_dir, "donor_chromosome_mean_log2_deviation.tsv"), sep = "\t")

# ------------------------- Chromosome-arm summaries --------------------------

arm_boundaries <- fread("CLLCMA_alpha/references/hg38_chromosome_arm_boundaries.tsv")
stopifnot(
  identical(as.character(arm_boundaries$chr), chr_levels),
  all(arm_boundaries$p_end == arm_boundaries$q_start)
)
gene_order[, midpoint := floor((start + end) / 2)]
gene_order <- merge(
  gene_order,
  arm_boundaries[, .(chr, q_start)],
  by = "chr",
  all.x = TRUE,
  sort = FALSE
)
gene_order[, arm := paste0(
  sub("^chr", "", as.character(chr)),
  fifelse(midpoint < q_start, "p", "q")
)]
gene_order[, arm := factor(arm, levels = arm_levels)]
setorder(gene_order, arm, start, end, gene)

arm_n_genes <- merge(
  data.table(arm = factor(arm_levels, levels = arm_levels)),
  gene_order[, .N, by = arm],
  by = "arm",
  all.x = TRUE,
  sort = FALSE
)
arm_n_genes[is.na(N), N := 0L]
arm_n_genes[, usable := N >= min_genes_per_arm]
fwrite(
  arm_n_genes,
  file.path(fig_dir, "chromosome_arm_gene_counts.tsv"),
  sep = "\t"
)
usable_arms <- as.character(arm_n_genes[usable == TRUE, arm])

arm_cell <- sapply(usable_arms, function(a) {
  genes <- gene_order[as.character(arm) == a, gene]
  colMeans(log_dev[genes, , drop = FALSE], na.rm = TRUE)
})
arm_cell <- as.data.table(arm_cell, keep.rownames = "cell_id")
arm_cell <- melt(
  arm_cell,
  id.vars = "cell_id",
  variable.name = "arm",
  value.name = "log2_dev",
  variable.factor = FALSE
)
arm_cell <- merge(arm_cell, anno, by = "cell_id")
arm_cell[, arm := factor(arm, levels = arm_levels)]

donor_arm <- arm_cell[, .(
  n_cells = .N,
  mean_log2_dev = mean(log2_dev, na.rm = TRUE),
  sd_cells = sd(log2_dev, na.rm = TRUE),
  median_log2_dev = median(log2_dev, na.rm = TRUE)
), by = .(donor, stage, dataset, arm)]
donor_arm <- donor_arm[n_cells >= min_cells_per_donor]
setorder(donor_arm, stage, donor, arm)
fwrite(
  donor_arm,
  file.path(fig_dir, "donor_chromosome_arm_mean_log2_deviation.tsv"),
  sep = "\t"
)

stage_arm <- donor_arm[, .(
  n_donors = .N,
  n_cells_total = sum(n_cells),
  mean_of_donor_means = mean(mean_log2_dev),
  sd_of_donor_means = sd(mean_log2_dev),
  median_of_donor_means = median(mean_log2_dev),
  q25 = quantile(mean_log2_dev, 0.25, names = FALSE),
  q75 = quantile(mean_log2_dev, 0.75, names = FALSE),
  min_donor_mean = min(mean_log2_dev),
  max_donor_mean = max(mean_log2_dev)
), by = .(stage, arm)]
arm_ref <- stage_arm[
  stage == "Normal_B",
  .(arm, ref_mean = mean_of_donor_means)
]
stage_arm <- merge(stage_arm, arm_ref, by = "arm")
stage_arm[, mean_minus_normal_b := mean_of_donor_means - ref_mean]
stage_arm[, arm := factor(arm, levels = arm_levels)]
stage_arm[, stage := factor(stage, levels = stage_order)]
setorder(stage_arm, arm, stage)
fwrite(
  stage_arm,
  file.path(fig_dir, "chromosome_arm_stage_summary.tsv"),
  sep = "\t"
)

obs_arm <- donor_arm[stage != "Normal_B"]
rank_arm <- obs_arm[, .(
  n_observation_donors = .N,
  mad_observation_donors = mad(mean_log2_dev),
  sd_observation_donors = sd(mean_log2_dev),
  mean_observation = mean(mean_log2_dev)
), by = arm]
rank_arm <- merge(rank_arm, arm_ref, by = "arm")
rank_arm[, mean_minus_normal_b := mean_observation - ref_mean]
setorder(rank_arm, -mad_observation_donors)
fwrite(
  rank_arm,
  file.path(fig_dir, "chromosome_arm_ranked_by_donor_mad.tsv"),
  sep = "\t"
)

# 80-gene bins to keep focal events from being averaged away.
gene_order[, bin := paste0(as.character(chr), ":", floor((seq_len(.N) - 1L) / bin_size) + 1L), by = chr]
bin_info <- gene_order[, .(
  chr = as.character(chr[1]),
  start = min(start),
  end = max(end),
  n_genes = .N
), by = bin]
bin_info[, chr := factor(chr, levels = chr_levels)]
setorder(bin_info, chr, start)
fwrite(bin_info, file.path(fig_dir, "chromosome_bins.tsv"), sep = "\t")

bin_cell <- sapply(bin_info$bin, function(b) {
  genes <- gene_order[bin == b, gene]
  colMeans(log_dev[genes, , drop = FALSE], na.rm = TRUE)
})
rownames(bin_cell) <- colnames(log_dev)
bin_donor <- sapply(bin_info$bin, function(b) {
  v <- bin_cell[, b]
  tapply(v[anno$cell_id], paste(anno$donor, anno$stage, sep = "||"), mean)
})
bin_donor <- as.data.table(bin_donor, keep.rownames = "donor_stage")
bin_donor[, c("donor", "stage") := tstrsplit(donor_stage, "||", fixed = TRUE)]
bin_donor[, donor_stage := NULL]
bin_long <- melt(
  bin_donor,
  id.vars = c("donor", "stage"),
  variable.name = "bin",
  value.name = "mean_log2_dev",
  variable.factor = FALSE
)
bin_long <- merge(bin_long, bin_info, by = "bin")
keep_donors <- unique(donor_chr[, .(donor, stage)])
bin_long <- merge(bin_long, keep_donors, by = c("donor", "stage"))
bin_long[, chr := factor(chr, levels = chr_levels)]
bin_long[, abs_dev := abs(mean_log2_dev)]
focal <- bin_long[bin_long[, .I[which.max(abs_dev)], by = .(donor, stage, chr)]$V1]
focal <- focal[, .(
  donor, stage, chr, bin,
  focal_start = start,
  focal_end = end,
  focal_n_genes = n_genes,
  focal_bin_mean_log2_dev = mean_log2_dev
)]
fwrite(focal, file.path(fig_dir, "donor_chromosome_focal_bin.tsv"), sep = "\t")

stage_chr <- donor_chr[, .(
  n_donors = .N,
  n_cells_total = sum(n_cells),
  mean_of_donor_means = mean(mean_log2_dev),
  sd_of_donor_means = sd(mean_log2_dev),
  median_of_donor_means = median(mean_log2_dev),
  q25 = quantile(mean_log2_dev, 0.25, names = FALSE),
  q75 = quantile(mean_log2_dev, 0.75, names = FALSE),
  min_donor_mean = min(mean_log2_dev),
  max_donor_mean = max(mean_log2_dev)
), by = .(stage, chr)]
ref <- stage_chr[stage == "Normal_B", .(chr, ref_mean = mean_of_donor_means)]
stage_chr <- merge(stage_chr, ref, by = "chr")
stage_chr[, mean_minus_normal_b := mean_of_donor_means - ref_mean]
stage_chr[, chr := factor(chr, levels = chr_levels)]
stage_chr[, stage := factor(stage, levels = stage_order)]
setorder(stage_chr, chr, stage)
fwrite(stage_chr, file.path(fig_dir, "chromosome_stage_summary.tsv"), sep = "\t")

obs <- donor_chr[stage != "Normal_B"]
rank_chr <- obs[, .(
  n_observation_donors = .N,
  mad_observation_donors = mad(mean_log2_dev),
  sd_observation_donors = sd(mean_log2_dev),
  mean_observation = mean(mean_log2_dev)
), by = chr]
rank_chr <- merge(rank_chr, ref, by = "chr")
rank_chr[, mean_minus_normal_b := mean_observation - ref_mean]
setorder(rank_chr, -mad_observation_donors)
fwrite(rank_chr, file.path(fig_dir, "chromosome_ranked_by_donor_mad.tsv"), sep = "\t")

donor_n <- unique(donor_chr[, .(donor, stage, n_cells)])
donor_n <- donor_n[, .(n_donors = .N, median_cells = median(n_cells), min_cells = min(n_cells), max_cells = max(n_cells)), by = stage]
fwrite(donor_n, file.path(fig_dir, "donor_counts_by_stage.tsv"), sep = "\t")

message("Donors retained (>=", min_cells_per_donor, " cells):")
print(donor_n)
message("Wrote ", normalizePath(fig_dir))
