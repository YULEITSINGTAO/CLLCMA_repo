#!/usr/bin/env Rscript

# Mutation–tumor microenvironment association (Broad CLL × deconvolution method).
#
# Uses non-silent mutations from data_mutations.txt and deconvolution fractions.
# For each gene × cell type:
#   - Wilcoxon (altered vs WT) + Cliff's delta (unadjusted descriptive effect)
#   - Primary association: linear model controlling for IGHV, sex, and age
#     at sampling (per decade): fraction ~ mutation + IGHV + sex + age10
#
# FDR:
#   - FDR_ALL / FDR_DRIVER : Storey adaptive q-values on covariate-adjusted p
#     (qvalue; estimates pi0). FDR_*_BH columns keep standard BH for comparison.
#
# Inputs:
#   CLLCMA_alpha/bulk/cll_broad_2022/data_mutations.txt
#   CLLCMA_alpha/bulk/cll_broad_2022/data_clinical_{sample,patient}.txt
#   CLLCMA_alpha/results/deconvolution/cll_broad_2022/<method>_sample_by_feature.tsv
#
# Outputs (shared with genomic_microenv_association.R):
#   CLLCMA_alpha/results/deconvolution/cll_broad_2022/alteration_microenv_association/<method>/
# Mutation-specific tables / driver FDR live here.
# The publication heatmap is written by plot_alteration_microenv_heatmap.R.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(qvalue)
})

project_dir <- Sys.getenv(
  "CLL_TME_DIR",
  unset = path.expand("~/research/CLL_TME")
)
project_dir <- normalizePath(project_dir, mustWork = TRUE)
setwd(project_dir)

deconv_method <- tolower(Sys.getenv("DECONV_METHOD", unset = "cibersort"))
method_label <- switch(
  deconv_method,
  cibersort = "CIBERSORT",
  bayesprism = "BayesPrism",
  mcp_counter = "MCP-counter",
  quantiseq = "quanTIseq",
  xcell = "xCell",
  deconv_method
)

bulk_dir <- file.path("CLLCMA_alpha", "bulk", "cll_broad_2022")
deconv_dir <- file.path(
  "CLLCMA_alpha", "results", "deconvolution", "cll_broad_2022"
)
deconv_file <- file.path(deconv_dir, paste0(deconv_method, "_sample_by_feature.tsv"))
if (!file.exists(deconv_file)) {
  stop("Missing deconvolution matrix: ", deconv_file)
}
output_dir <- file.path(deconv_dir, "alteration_microenv_association", deconv_method)
plot_dir <- file.path(output_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
message("Method: ", deconv_method, " -> ", normalizePath(output_dir, mustWork = FALSE))

min_mutated_samples <- 15L
fdr_cutoff <- 0.05
nominal_p_cutoff <- 0.05

cll_driver_genes <- c(
  "TP53", "ATM", "SF3B1", "NOTCH1", "POT1", "BIRC3", "MYD88", "XPO1",
  "CHD2", "KLHL6", "MGA", "BRAF", "KRAS", "NRAS", "MED12", "FBXW7",
  "EGR2", "BCOR", "DDX3X", "ZMYM3", "NFKBIE", "RPS15", "IKZF3",
  "ARID1A", "ZNF292", "IGLL5"
)

silent_like <- c(
  "Silent", "Intron", "3'UTR", "5'UTR", "3'Flank", "5'Flank",
  "IGR", "RNA", "lincRNA", "Targeted_Region"
)

cell_type_order <- c(
  "CLL", "RS", "CD4_T", "CD8_T", "NK", "Monocyte", "DC", "ILC", "Progenitor"
)

# ------------------------------ Helpers --------------------------------------

sanitize_filename <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
}

# Storey adaptive q-values (estimates pi0). Falls back to BH if qvalue fails.
adaptive_qvalue <- function(p) {
  p <- as.numeric(p)
  n <- length(p)
  q <- rep(NA_real_, n)
  pi0 <- NA_real_
  ok <- which(is.finite(p) & p >= 0 & p <= 1)
  if (!length(ok)) {
    return(list(qvalues = q, pi0 = pi0))
  }
  if (length(ok) < 5L) {
    q[ok] <- p.adjust(p[ok], method = "BH")
    pi0 <- 1
    return(list(qvalues = q, pi0 = pi0))
  }
  qobj <- tryCatch(
    qvalue(p[ok]),
    error = function(e) NULL
  )
  if (is.null(qobj)) {
    q[ok] <- p.adjust(p[ok], method = "BH")
    pi0 <- 1
  } else {
    q[ok] <- qobj$qvalues
    pi0 <- qobj$pi0
  }
  list(qvalues = q, pi0 = pi0)
}

source(file.path(
  project_dir, "CLLCMA_alpha", "scripts", "helpers_microenv_dotplot.R"
))
source(file.path(
  project_dir, "CLLCMA_alpha", "scripts", "helpers_microenv_covariates.R"
))

read_cbioportal_clinical <- function(path) {
  fread(
    path,
    skip = 4L,
    na.strings = c("NA", "", "[Not Available]", "[Not Applicable]")
  )
}

safe_wilcox <- function(a, b) {
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  if (length(a) < 3L || length(b) < 3L) {
    return(NA_real_)
  }
  tryCatch(
    wilcox.test(a, b, alternative = "two.sided", exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

rank_biserial <- function(altered, wt) {
  altered <- altered[is.finite(altered)]
  wt <- wt[is.finite(wt)]
  n1 <- length(altered)
  n0 <- length(wt)
  if (!n1 || !n0) {
    return(NA_real_)
  }
  ranks <- rank(c(altered, wt), ties.method = "average")
  u1 <- sum(ranks[seq_len(n1)]) - n1 * (n1 + 1) / 2
  2 * u1 / (n1 * n0) - 1
}

# ------------------------------- Load data -----------------------------------

cs <- fread(deconv_file)
feature_cols <- intersect(cell_type_order, setdiff(names(cs), "SAMPLE_ID"))
feature_cols <- c(
  feature_cols,
  setdiff(setdiff(names(cs), "SAMPLE_ID"), feature_cols)
)

clinical_sample <- read_cbioportal_clinical(
  file.path(bulk_dir, "data_clinical_sample.txt")
)
clinical_patient <- read_cbioportal_clinical(
  file.path(bulk_dir, "data_clinical_patient.txt")
)

micro <- attach_clinical_covariates(
  cs[, c("SAMPLE_ID", feature_cols), with = FALSE],
  clinical_sample,
  clinical_patient
)
sample_universe <- unique(micro$SAMPLE_ID)
message(sprintf(
  "Covariate coverage: IGHV=%d sex=%d age=%d / %d samples",
  sum(is.finite(micro$IGHV_MUTATED)),
  sum(is.finite(micro$SEX_MALE)),
  sum(is.finite(micro$AGE_YEARS)),
  nrow(micro)
))

mut <- fread(file.path(bulk_dir, "data_mutations.txt"), skip = 2L)
mut <- mut[
  !is.na(Hugo_Symbol) & Hugo_Symbol != "" &
    Tumor_Sample_Barcode %in% sample_universe &
    !Variant_Classification %in% silent_like
]
mut[, SAMPLE_ID := as.character(Tumor_Sample_Barcode)]

# Samples with any mutation call (proxy for sequenced / MAF-covered)
assayed_samples <- unique(mut$SAMPLE_ID)
micro_g <- micro[SAMPLE_ID %in% assayed_samples]
message(
  "Mutation-assayed ", method_label, " samples: ", nrow(micro_g),
  " / ", length(sample_universe)
)

# Gene panel
gene_n <- mut[, .(N_MUTATED = uniqueN(SAMPLE_ID)), by = Hugo_Symbol]
keep_genes <- sort(unique(c(
  gene_n[N_MUTATED >= min_mutated_samples, Hugo_Symbol],
  intersect(cll_driver_genes, gene_n$Hugo_Symbol)
)))
driver_genes <- intersect(cll_driver_genes, keep_genes)
message(
  "Genes tested: ", length(keep_genes),
  " (recurrent >=", min_mutated_samples, " or drivers); drivers: ",
  length(driver_genes)
)

mut_by_gene <- unique(mut[Hugo_Symbol %in% keep_genes, .(SAMPLE_ID, Hugo_Symbol)])
mut_map <- split(mut_by_gene$SAMPLE_ID, mut_by_gene$Hugo_Symbol)

# --------------------------- Gene × cell-type tests --------------------------

rows <- list()
for (gene in keep_genes) {
  alt_ids <- mut_map[[gene]]
  altered <- as.integer(micro_g$SAMPLE_ID %in% alt_ids)
  n_alt <- sum(altered == 1L)
  n_wt <- sum(altered == 0L)
  if (n_alt < 5L || n_wt < 5L) {
    next
  }
  for (ct in feature_cols) {
    vals <- micro_g[[ct]]
    a <- vals[altered == 1L]
    w <- vals[altered == 0L]
    p_w <- safe_wilcox(a, w)
    delta <- rank_biserial(a, w)
    adj <- fit_alt_covariate_adjusted(
      y = vals,
      altered = altered,
      ighv_mutated = micro_g$IGHV_MUTATED,
      sex_male = micro_g$SEX_MALE,
      age_years = micro_g$AGE_YEARS
    )
    rows[[length(rows) + 1L]] <- data.table(
      GENE = gene,
      IS_DRIVER = gene %in% driver_genes,
      CELL_TYPE = ct,
      N_MUTATED = n_alt,
      N_WT = n_wt,
      MEDIAN_MUT = median(a, na.rm = TRUE),
      MEDIAN_WT = median(w, na.rm = TRUE),
      MEAN_MUT = mean(a, na.rm = TRUE),
      MEAN_WT = mean(w, na.rm = TRUE),
      CLIFFS_DELTA = delta,
      P_WILCOX = p_w,
      BETA_ADJ = adj$beta,
      SE_ADJ = adj$se,
      P_ADJ = adj$p,
      N_ADJ = adj$n
    )
  }
}

assoc <- rbindlist(rows, use.names = TRUE)
# Primary FDR on covariate-adjusted p-values
assoc[, FDR_ALL_BH := p.adjust(P_ADJ, method = "BH")]
aq_all <- adaptive_qvalue(assoc$P_ADJ)
assoc[, FDR_ALL := aq_all$qvalues]
message(sprintf("Adaptive FDR_ALL (covariate-adjusted) pi0=%.3f", aq_all$pi0))

assoc[, FDR_DRIVER := NA_real_]
assoc[, FDR_DRIVER_BH := NA_real_]
driver_idx <- which(assoc$IS_DRIVER == TRUE)
if (length(driver_idx)) {
  assoc$FDR_DRIVER_BH[driver_idx] <- p.adjust(assoc$P_ADJ[driver_idx], method = "BH")
  aq_drv <- adaptive_qvalue(assoc$P_ADJ[driver_idx])
  assoc$FDR_DRIVER[driver_idx] <- aq_drv$qvalues
  message(sprintf("Adaptive FDR_DRIVER (covariate-adjusted) pi0=%.3f", aq_drv$pi0))
}
assoc[, DIRECTION := fifelse(
  is.na(CLIFFS_DELTA), NA_character_,
  fifelse(CLIFFS_DELTA > 0, "higher_in_mutated", "lower_in_mutated")
)]
assoc[, DIRECTION_ADJ := fifelse(
  is.na(BETA_ADJ), NA_character_,
  fifelse(BETA_ADJ > 0, "higher_in_mutated", "lower_in_mutated")
)]
setorder(assoc, P_ADJ)

fwrite(assoc, file.path(output_dir, "mutation_celltype_association.tsv"), sep = "\t")

sig_all <- assoc[FDR_ALL < fdr_cutoff]
sig_driver <- assoc[IS_DRIVER == TRUE & FDR_DRIVER < fdr_cutoff]
nominal_driver <- assoc[IS_DRIVER == TRUE & P_ADJ < nominal_p_cutoff]

fwrite(sig_all, file.path(output_dir, "significant_fdr_all.tsv"), sep = "\t")
fwrite(sig_driver, file.path(output_dir, "significant_fdr_driver.tsv"), sep = "\t")
# Compatibility alias for prior IGHV-only filename (now full covariates)
fwrite(
  sig_driver,
  file.path(output_dir, "significant_fdr_driver_covariate_adjusted.tsv"),
  sep = "\t"
)
fwrite(
  nominal_driver[order(P_ADJ)],
  file.path(output_dir, "driver_nominal_p05.tsv"),
  sep = "\t"
)
fwrite(
  assoc[order(P_ADJ)][, head(.SD, 50L)],
  file.path(output_dir, "top50_by_pvalue.tsv"),
  sep = "\t"
)

# Best cell type per gene (by adjusted p)
best <- assoc[, .SD[which.min(P_ADJ)], by = GENE][order(P_ADJ)]
fwrite(best, file.path(output_dir, "best_celltype_per_gene.tsv"), sep = "\t")

# ------------------------------ Key findings ---------------------------------

key <- rbindlist(list(
  if (nrow(sig_driver)) {
    copy(sig_driver)[, EVIDENCE := "driver_FDR_covariate_adjusted"]
  } else {
    NULL
  },
  if (nrow(sig_all)) {
    copy(sig_all)[, EVIDENCE := "all_genes_FDR_covariate_adjusted"]
  } else {
    NULL
  },
  copy(nominal_driver)[, EVIDENCE := "driver_nominal_p05_covariate_adjusted"]
), use.names = TRUE, fill = TRUE)
if (nrow(key)) {
  key <- unique(key, by = c("GENE", "CELL_TYPE", "EVIDENCE"))
  setorder(key, EVIDENCE, P_ADJ)
}
fwrite(key, file.path(output_dir, "key_mutation_microenv_hits.tsv"), sep = "\t")

lines <- c(
  "MUTATION–TUMOR MICROENVIRONMENT ASSOCIATION",
  "Primary model: fraction ~ mutation + IGHV + sex + age_sampling/10",
  sprintf("Assayed samples with deconvolution + mutations: %d", nrow(micro_g)),
  sprintf("Genes tested: %d (min mutated samples=%d)", length(keep_genes), min_mutated_samples),
  sprintf("Driver genes in panel: %s", paste(driver_genes, collapse = ", ")),
  "",
  sprintf("FDR_ALL (covariate-adj) < %.2f hits: %d", fdr_cutoff, nrow(sig_all)),
  sprintf("FDR_DRIVER (covariate-adj) < %.2f hits: %d", fdr_cutoff, nrow(sig_driver)),
  sprintf("Driver nominal adj-p < %.2f hits: %d", nominal_p_cutoff, nrow(nominal_driver)),
  "",
  "=== Top driver associations (by covariate-adjusted p) ==="
)
top_d <- assoc[IS_DRIVER == TRUE][order(P_ADJ)]
for (i in seq_len(min(25L, nrow(top_d)))) {
  r <- top_d[i]
  lines <- c(lines, sprintf(
    "%s ~ %s | n_mut=%d | delta=%+.2f | beta_adj=%+.3g | adj-p=%.3g FDR_driver=%s FDR_all=%s | wilcox-p=%s",
    r$GENE, r$CELL_TYPE, r$N_MUTATED, r$CLIFFS_DELTA,
    ifelse(is.na(r$BETA_ADJ), NA, r$BETA_ADJ),
    ifelse(is.na(r$P_ADJ), NA, r$P_ADJ),
    ifelse(is.na(r$FDR_DRIVER), "NA", sprintf("%.3g", r$FDR_DRIVER)),
    ifelse(is.na(r$FDR_ALL), "NA", sprintf("%.3g", r$FDR_ALL)),
    ifelse(is.na(r$P_WILCOX), "NA", sprintf("%.3g", r$P_WILCOX))
  ))
}
writeLines(lines, file.path(output_dir, "KEY_FINDINGS.txt"))

message(paste(lines[1:min(20L, length(lines))], collapse = "\n"))

# --------------------------------- Plots -------------------------------------

# Mutation rate vs mean cell-type fraction in mutated samples
plot_rate_vs_mean_celltype(
  dt = assoc,
  gene_col = "GENE",
  n_alt_col = "N_MUTATED",
  n_wt_col = "N_WT",
  mean_alt_col = "MEAN_MUT",
  p_col = "P_ADJ",
  fdr_col = "FDR_ALL",
  out_file = file.path(output_dir, "mutation_rate_vs_celltype_mean.pdf"),
  title = "Mutation rate vs mean microenvironment percentage",
  subtitle = "Each point is a gene; significance from IGHV/sex/age-adjusted model",
  highlight_genes = c("TP53", "ATM", "SF3B1", "NOTCH1", "POT1", "MYD88"),
  cell_levels = feature_cols
)

# Main heatmap is produced by plot_alteration_microenv_heatmap.R

# Boxplots for strongest driver hits
plot_set <- unique(rbindlist(list(
  sig_driver,
  assoc[IS_DRIVER == TRUE][order(P_ADJ)][, head(.SD, 20L)]
), use.names = TRUE, fill = TRUE), by = c("GENE", "CELL_TYPE"))

for (i in seq_len(nrow(plot_set))) {
  r <- plot_set[i]
  dt <- copy(micro_g[, c("SAMPLE_ID", r$CELL_TYPE), with = FALSE])
  setnames(dt, r$CELL_TYPE, "FRACTION")
  dt[, GROUP := fifelse(SAMPLE_ID %in% mut_map[[r$GENE]], "Mutated", "WT")]
  dt <- dt[is.finite(FRACTION)]
  dt[, GROUP := factor(GROUP, levels = c("WT", "Mutated"))]
  p <- ggplot(dt, aes(x = GROUP, y = FRACTION, fill = GROUP)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.85) +
    geom_jitter(width = 0.12, size = 0.9, alpha = 0.35) +
    scale_fill_manual(values = c(WT = "#9E9E9E", Mutated = "#C44E52")) +
    labs(
      title = sprintf("%s mutation · %s", r$GENE, r$CELL_TYPE),
      subtitle = sprintf(
        "n_mut=%d n_wt=%d | delta=%+.2f | adj-p=%.2g | FDR_driver=%s",
        r$N_MUTATED, r$N_WT, r$CLIFFS_DELTA, r$P_ADJ,
        ifelse(is.na(r$FDR_DRIVER), "NA", sprintf("%.2g", r$FDR_DRIVER))
      ),
      x = NULL,
      y = paste0(method_label, " fraction / score")
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))
  ggsave(
    file.path(
      plot_dir,
      sprintf("%s__%s.pdf", sanitize_filename(r$GENE), sanitize_filename(r$CELL_TYPE))
    ),
    p,
    width = 4.2,
    height = 4.2
  )
}

message("Done -> ", normalizePath(output_dir))
message(
  "Driver FDR hits (covariate-adj): ", nrow(sig_driver),
  "; all-gene FDR hits: ", nrow(sig_all)
)
