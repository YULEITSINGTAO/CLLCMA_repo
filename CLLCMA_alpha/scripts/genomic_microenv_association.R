#!/usr/bin/env Rscript

# Associate Broad CLL genomic alterations with deconvolution microenvironment.
#
# Alteration classes:
#   mutation  — non-silent gene-level (mutated vs WT)
#   focal_cna — focal locus Gain/Loss vs Unchanged (and Gain-only / Loss-only)
#   arm_cna   — arm-level Gain/Loss vs Unchanged
#   sv        — gene involved in SV (Site1 or Site2) vs not
#
# Primary association: linear model controlling for IGHV, sex, and age at
# sampling (per decade): fraction ~ alteration + IGHV + sex + age10.
# Also report unadjusted Wilcoxon p and Cliff's delta for description.
# Adaptive FDR (Storey qvalue) within each ALTERATION_CLASS on adjusted p.
# Standard BH is also stored as FDR_BH for comparison.
#
# Inputs:
#   CLLCMA_alpha/results/deconvolution/cll_broad_2022/<method>_sample_by_feature.tsv
#   CLLCMA_alpha/bulk/cll_broad_2022/data_{mutations,focal_cna,armlevel_cna,sv}.txt
#   CLLCMA_alpha/bulk/cll_broad_2022/data_clinical_{sample,patient}.txt
#
# Outputs (shared with mutation_microenv_association.R):
#   CLLCMA_alpha/results/deconvolution/cll_broad_2022/alteration_microenv_association/<method>/
#   association_*.tsv — mutation / focal_cna / arm_cna / sv tests
#   alteration_microenv_heatmap.pdf via plot_alteration_microenv_heatmap.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
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

min_altered <- 20L
fdr_cutoff <- 0.05
# Relaxed FDR for the main heatmap / key-alteration tables (report more associations)
heatmap_fdr_cutoff <- 0.10
max_plots_per_class <- 40L

# Always keep these CLL drivers in the mutation panel if present at all.
cll_driver_genes <- c(
  "TP53", "ATM", "SF3B1", "NOTCH1", "POT1", "BIRC3", "MYD88", "XPO1",
  "CHD2", "KLHL6", "MGA", "BRAF", "KRAS", "NRAS", "MED12", "FBXW7",
  "EGR2", "BCOR", "DDX3X", "ZMYM3", "NFKBIE", "RPS15", "IKZF3"
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
  project_dir, "CLLCMA_alpha", "scripts", "helpers_microenv_covariates.R"
))

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

# Cliff's delta / rank-biserial: +1 = altered group always higher.
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

read_cbioportal_clinical <- function(path) {
  fread(
    path,
    skip = 4L,
    na.strings = c("NA", "", "[Not Available]", "[Not Applicable]")
  )
}

build_binary_long <- function(sample_ids, altered_map, class_name) {
  # altered_map: named list or data.table of ALTERATION -> character vector of
  # SAMPLE_IDs that are altered.
  rows <- lapply(names(altered_map), function(alt) {
    alt_ids <- unique(as.character(altered_map[[alt]]))
    alt_ids <- intersect(alt_ids, sample_ids)
    data.table(
      SAMPLE_ID = sample_ids,
      ALTERATION_CLASS = class_name,
      ALTERATION = alt,
      ALTERED = as.integer(sample_ids %in% alt_ids)
    )
  })
  rbindlist(rows, use.names = TRUE)
}

# --------------------------- Load microenvironment ---------------------------

cs <- fread(deconv_file)
feature_cols <- intersect(cell_type_order, setdiff(names(cs), "SAMPLE_ID"))
feature_cols <- c(feature_cols, setdiff(setdiff(names(cs), "SAMPLE_ID"), feature_cols))

clinical_sample <- read_cbioportal_clinical(
  file.path(bulk_dir, "data_clinical_sample.txt")
)
clinical_patient <- read_cbioportal_clinical(
  file.path(bulk_dir, "data_clinical_patient.txt")
)

# Genomics are sample-keyed; use sample-level fractions matched to SAMPLE_ID.
micro <- attach_clinical_covariates(
  cs[, c("SAMPLE_ID", feature_cols), with = FALSE],
  clinical_sample,
  clinical_patient
)

sample_universe <- unique(micro$SAMPLE_ID)
message(method_label, " samples: ", length(sample_universe))
message(sprintf(
  "Covariate coverage: IGHV=%d sex=%d age=%d / %d samples",
  sum(is.finite(micro$IGHV_MUTATED)),
  sum(is.finite(micro$SEX_MALE)),
  sum(is.finite(micro$AGE_YEARS)),
  nrow(micro)
))

# ------------------------------- Mutations -----------------------------------

mut <- fread(file.path(bulk_dir, "data_mutations.txt"), skip = 2L)
mut <- mut[
  !is.na(Hugo_Symbol) & Hugo_Symbol != "" &
    !is.na(Tumor_Sample_Barcode) &
    Tumor_Sample_Barcode %in% sample_universe &
    !Variant_Classification %in% silent_like
]
mut[, SAMPLE_ID := as.character(Tumor_Sample_Barcode)]
mut_by_gene <- mut[, .(SAMPLE_ID = unique(SAMPLE_ID)), by = Hugo_Symbol]
gene_n <- mut_by_gene[, .(N_ALTERED = .N), by = Hugo_Symbol]
keep_genes <- unique(c(
  gene_n[N_ALTERED >= min_altered, Hugo_Symbol],
  intersect(cll_driver_genes, gene_n$Hugo_Symbol)
))
mut_by_gene <- mut_by_gene[Hugo_Symbol %in% keep_genes]
mut_map <- split(mut_by_gene$SAMPLE_ID, mut_by_gene$Hugo_Symbol)
message(
  "Mutation genes tested: ", length(mut_map),
  " (min altered=", min_altered, " or CLL drivers)"
)

# ------------------------------- Focal CNA -----------------------------------

focal <- fread(file.path(bulk_dir, "data_focal_cna.txt"))
focal_samples <- intersect(
  setdiff(names(focal), c("ENTITY_STABLE_ID", "NAME", "DESCRIPTION")),
  sample_universe
)
focal_long <- melt(
  focal,
  id.vars = c("ENTITY_STABLE_ID", "NAME"),
  measure.vars = focal_samples,
  variable.name = "SAMPLE_ID",
  value.name = "STATUS",
  variable.factor = FALSE
)
focal_long[, SAMPLE_ID := as.character(SAMPLE_ID)]
focal_long[, STATUS := trimws(as.character(STATUS))]
focal_long <- focal_long[STATUS %in% c("Gain", "Loss", "Unchanged")]

# Binary: any alteration (Gain or Loss)
focal_any <- focal_long[
  STATUS %in% c("Gain", "Loss"),
  .(SAMPLE_ID = unique(SAMPLE_ID)),
  by = .(ALTERATION = paste0(NAME, "_altered"))
]
# Separate Gain / Loss when recurrent
focal_gain <- focal_long[
  STATUS == "Gain",
  .(SAMPLE_ID = unique(SAMPLE_ID)),
  by = .(ALTERATION = paste0(NAME, "_Gain"))
]
focal_loss <- focal_long[
  STATUS == "Loss",
  .(SAMPLE_ID = unique(SAMPLE_ID)),
  by = .(ALTERATION = paste0(NAME, "_Loss"))
]
focal_events <- rbindlist(list(focal_any, focal_gain, focal_loss), use.names = TRUE)
focal_n <- focal_events[, .(N_ALTERED = uniqueN(SAMPLE_ID)), by = ALTERATION]
focal_keep <- focal_n[N_ALTERED >= min_altered, ALTERATION]
focal_events <- focal_events[ALTERATION %in% focal_keep]
focal_map <- split(focal_events$SAMPLE_ID, focal_events$ALTERATION)
message("Focal CNA events tested: ", length(focal_map))

# ------------------------------- Arm CNA -------------------------------------

arm <- fread(file.path(bulk_dir, "data_armlevel_cna.txt"))
arm_samples <- intersect(
  setdiff(names(arm), c("ENTITY_STABLE_ID", "NAME", "DESCRIPTION")),
  sample_universe
)
arm_long <- melt(
  arm,
  id.vars = c("ENTITY_STABLE_ID", "NAME"),
  measure.vars = arm_samples,
  variable.name = "SAMPLE_ID",
  value.name = "STATUS",
  variable.factor = FALSE
)
arm_long[, SAMPLE_ID := as.character(SAMPLE_ID)]
arm_long[, STATUS := trimws(as.character(STATUS))]
arm_any <- arm_long[
  STATUS %in% c("Gain", "Loss"),
  .(SAMPLE_ID = unique(SAMPLE_ID)),
  by = .(ALTERATION = paste0("arm_", NAME, "_altered"))
]
arm_gain <- arm_long[
  STATUS == "Gain",
  .(SAMPLE_ID = unique(SAMPLE_ID)),
  by = .(ALTERATION = paste0("arm_", NAME, "_Gain"))
]
arm_loss <- arm_long[
  STATUS == "Loss",
  .(SAMPLE_ID = unique(SAMPLE_ID)),
  by = .(ALTERATION = paste0("arm_", NAME, "_Loss"))
]
arm_events <- rbindlist(list(arm_any, arm_gain, arm_loss), use.names = TRUE)
arm_n <- arm_events[, .(N_ALTERED = uniqueN(SAMPLE_ID)), by = ALTERATION]
arm_keep <- arm_n[N_ALTERED >= min_altered, ALTERATION]
arm_events <- arm_events[ALTERATION %in% arm_keep]
arm_map <- split(arm_events$SAMPLE_ID, arm_events$ALTERATION)
message("Arm CNA events tested: ", length(arm_map))

# ----------------------------------- SV --------------------------------------

sv <- fread(file.path(bulk_dir, "data_sv.txt"))
sv <- sv[Sample_Id %in% sample_universe]
sv_genes <- rbindlist(list(
  sv[!is.na(Site1_Hugo_Symbol) & Site1_Hugo_Symbol != "",
     .(SAMPLE_ID = Sample_Id, gene = Site1_Hugo_Symbol)],
  sv[!is.na(Site2_Hugo_Symbol) & Site2_Hugo_Symbol != "",
     .(SAMPLE_ID = Sample_Id, gene = Site2_Hugo_Symbol)]
))
sv_by_gene <- unique(sv_genes)[, .(SAMPLE_ID = unique(SAMPLE_ID)), by = gene]
sv_n <- sv_by_gene[, .(N_ALTERED = .N), by = gene]
# SV is rarer — lower threshold, plus drivers
sv_keep <- unique(c(
  sv_n[N_ALTERED >= 10L, gene],
  intersect(cll_driver_genes, sv_n$gene)
))
sv_by_gene <- sv_by_gene[gene %in% sv_keep]
sv_map <- split(sv_by_gene$SAMPLE_ID, paste0("SV_", sv_by_gene$gene))
message("SV gene events tested: ", length(sv_map))

# --------------------------- Association tests -------------------------------

alteration_tables <- list(
  mutation = build_binary_long(sample_universe, mut_map, "mutation"),
  focal_cna = build_binary_long(sample_universe, focal_map, "focal_cna"),
  arm_cna = build_binary_long(sample_universe, arm_map, "arm_cna"),
  sv = build_binary_long(sample_universe, sv_map, "sv")
)
# Drop empty classes
alteration_tables <- alteration_tables[vapply(alteration_tables, nrow, 1L) > 0]

assoc_rows <- list()

for (class_name in names(alteration_tables)) {
  alt_dt <- alteration_tables[[class_name]]
  # Restrict to samples that have data for this class (appeared in that file)
  if (class_name == "mutation") {
    assayed <- intersect(sample_universe, unique(mut$SAMPLE_ID))
  } else if (class_name == "focal_cna") {
    assayed <- focal_samples
  } else if (class_name == "arm_cna") {
    assayed <- arm_samples
  } else {
    assayed <- intersect(sample_universe, unique(sv$Sample_Id))
    # For SV, WT = no SV involving that gene among samples with any SV call
    # OR among all CIBERSORT samples with SV file coverage.
    # Use samples present in SV file OR all micro samples — SV file is
    # event-only, so absence ≠ assayed. Use all CIBERSORT samples that also
    # have WES-like clinical sample overlap; treat unlisted as WT only among
    # samples that appear at least once in SV OR mutation/CNA panels.
    # Practical approach: WT = not altered among all CIBERSORT samples that
    # have mutation data (proxy for genomically profiled).
    assayed <- intersect(sample_universe, unique(c(mut$SAMPLE_ID, focal_samples)))
  }
  alt_dt <- alt_dt[SAMPLE_ID %in% assayed]

  for (alt in unique(alt_dt$ALTERATION)) {
    sub <- alt_dt[ALTERATION == alt]
    n_alt <- sum(sub$ALTERED == 1L)
    n_wt <- sum(sub$ALTERED == 0L)
    if (n_alt < 5L || n_wt < 5L) {
      next
    }
    merged <- merge(
      micro[, c(
        "SAMPLE_ID", feature_cols,
        "IGHV_MUTATED", "SEX_MALE", "AGE_YEARS"
      ), with = FALSE],
      sub[, .(SAMPLE_ID, ALTERED)],
      by = "SAMPLE_ID",
      all = FALSE
    )
    for (ct in feature_cols) {
      vals <- merged[[ct]]
      a <- vals[merged$ALTERED == 1L]
      w <- vals[merged$ALTERED == 0L]
      p_w <- safe_wilcox(a, w)
      delta <- rank_biserial(a, w)
      adj <- fit_alt_covariate_adjusted(
        y = vals,
        altered = merged$ALTERED,
        ighv_mutated = merged$IGHV_MUTATED,
        sex_male = merged$SEX_MALE,
        age_years = merged$AGE_YEARS
      )
      assoc_rows[[length(assoc_rows) + 1L]] <- data.table(
        ALTERATION_CLASS = class_name,
        ALTERATION = alt,
        CELL_TYPE = ct,
        N_ALTERED = length(a[is.finite(a)]),
        N_WT = length(w[is.finite(w)]),
        MEDIAN_ALTERED = median(a, na.rm = TRUE),
        MEDIAN_WT = median(w, na.rm = TRUE),
        MEAN_ALTERED = mean(a, na.rm = TRUE),
        MEAN_WT = mean(w, na.rm = TRUE),
        CLIFFS_DELTA = delta,
        P_WILCOX = p_w,
        BETA_ADJ = adj$beta,
        SE_ADJ = adj$se,
        P_ADJ = adj$p,
        N_ADJ = adj$n,
        P_VALUE = adj$p
      )
    }
  }
}

assoc <- rbindlist(assoc_rows, use.names = TRUE, fill = TRUE)
if (!nrow(assoc)) {
  stop("No association tests were produced")
}

# Primary FDR on covariate-adjusted p (P_VALUE == P_ADJ)
assoc[, FDR_BH := p.adjust(P_VALUE, method = "BH"), by = ALTERATION_CLASS]
assoc[, FDR := NA_real_]

pi0_by_class <- list()
for (cls in unique(assoc$ALTERATION_CLASS)) {
  idx <- which(assoc$ALTERATION_CLASS == cls)
  aq <- adaptive_qvalue(assoc$P_VALUE[idx])
  set(assoc, i = idx, j = "FDR", value = aq$qvalues)
  pi0_by_class[[cls]] <- aq$pi0
}
assoc[, DIRECTION := fifelse(
  is.na(CLIFFS_DELTA),
  NA_character_,
  fifelse(CLIFFS_DELTA > 0, "higher_in_altered", "lower_in_altered")
)]
assoc[, DIRECTION_ADJ := fifelse(
  is.na(BETA_ADJ),
  NA_character_,
  fifelse(BETA_ADJ > 0, "higher_in_altered", "lower_in_altered")
)]
setorder(assoc, ALTERATION_CLASS, P_VALUE)

fwrite(assoc, file.path(output_dir, "association_all.tsv"), sep = "\t")

sig <- assoc[FDR < fdr_cutoff]
fwrite(sig, file.path(output_dir, "association_significant_fdr05.tsv"), sep = "\t")

# Key mutations / arms for reporting + heatmap (relaxed FDR)
key_mut <- assoc[
  ALTERATION_CLASS == "mutation" & FDR < heatmap_fdr_cutoff
][order(FDR, -abs(CLIFFS_DELTA))]
fwrite(key_mut, file.path(output_dir, "key_mutations_microenv.tsv"), sep = "\t")

# Key arm CNAs for reporting: prefer "_altered" to avoid Gain/Loss duplicates
key_arm <- assoc[
  ALTERATION_CLASS == "arm_cna" &
    grepl("_altered$", ALTERATION) &
    FDR < heatmap_fdr_cutoff
][order(FDR, -abs(CLIFFS_DELTA))]
fwrite(key_arm, file.path(output_dir, "key_arm_cna_microenv.tsv"), sep = "\t")

# Key focal CNAs (_altered)
key_focal <- assoc[
  ALTERATION_CLASS == "focal_cna" &
    grepl("_altered$", ALTERATION) &
    FDR < heatmap_fdr_cutoff
][order(FDR, -abs(CLIFFS_DELTA))]
fwrite(key_focal, file.path(output_dir, "key_focal_cna_microenv.tsv"), sep = "\t")

# Combined mutation + CNV set for key-alteration table
key_heatmap <- rbindlist(
  list(key_mut, key_arm, key_focal),
  use.names = TRUE,
  fill = TRUE
)[order(ALTERATION_CLASS, FDR, -abs(CLIFFS_DELTA))]
fwrite(key_heatmap, file.path(output_dir, "key_alterations_microenv.tsv"), sep = "\t")

# Also top mutation hits by covariate-adjusted p even if FDR soft (exploratory)
top_mut_raw <- assoc[
  ALTERATION_CLASS == "mutation"
][order(P_VALUE)][, head(.SD, 50L)]
fwrite(top_mut_raw, file.path(output_dir, "top50_mutations_by_pvalue.tsv"), sep = "\t")

# Summary: best cell-type hit per alteration
best_per_alt <- assoc[
  ,
  .SD[which.min(P_VALUE)],
  by = .(ALTERATION_CLASS, ALTERATION)
][order(ALTERATION_CLASS, P_VALUE)]
fwrite(best_per_alt, file.path(output_dir, "best_hit_per_alteration.tsv"), sep = "\t")

# Class-level summary counts
summary_tab <- assoc[
  ,
  .(
    N_TESTS = .N,
    PI0 = NA_real_,
    N_SIG_FDR05 = sum(FDR < fdr_cutoff, na.rm = TRUE),
    N_SIG_FDR10 = sum(FDR < heatmap_fdr_cutoff, na.rm = TRUE),
    N_SIG_BH05 = sum(FDR_BH < fdr_cutoff, na.rm = TRUE),
    N_SIG_P001 = sum(P_VALUE < 0.001, na.rm = TRUE)
  ),
  by = ALTERATION_CLASS
]
for (i in seq_len(nrow(summary_tab))) {
  cls <- summary_tab$ALTERATION_CLASS[i]
  summary_tab$PI0[i] <- pi0_by_class[[cls]]
}
fwrite(summary_tab, file.path(output_dir, "association_summary.tsv"), sep = "\t")

message("Association summary (IGHV/sex/age-adjusted; adaptive FDR / qvalue):")
print(summary_tab)
message(sprintf(
  "Heatmap mutation hits (covariate-adj adaptive FDR<%.2f): %d",
  heatmap_fdr_cutoff, nrow(key_mut)
))
if (nrow(key_mut)) {
  print(key_mut[, .(
    ALTERATION, CELL_TYPE, N_ALTERED, CLIFFS_DELTA, BETA_ADJ, P_ADJ, FDR
  )][1:min(30L, .N)])
}
message(sprintf(
  "Heatmap arm-CNA (_altered) hits (covariate-adj adaptive FDR<%.2f): %d",
  heatmap_fdr_cutoff, nrow(key_arm)
))
if (nrow(key_arm)) {
  print(key_arm[, .(
    ALTERATION, CELL_TYPE, N_ALTERED, CLIFFS_DELTA, BETA_ADJ, P_ADJ, FDR
  )])
}
message(sprintf(
  "Heatmap focal-CNA (_altered) hits (covariate-adj adaptive FDR<%.2f): %d",
  heatmap_fdr_cutoff, nrow(key_focal)
))
if (nrow(key_focal)) {
  print(key_focal[, .(
    ALTERATION, CELL_TYPE, N_ALTERED, CLIFFS_DELTA, BETA_ADJ, P_ADJ, FDR
  )][1:min(20L, .N)])
}

# --------------------------------- Plots -------------------------------------

# Mutation-rate scatter is written by mutation_microenv_association.R into the
# same alteration_microenv_association/<method>/ folder.

plot_one <- function(row, micro_dt, alt_samples) {
  ct <- row$CELL_TYPE
  alt <- row$ALTERATION
  class_name <- row$ALTERATION_CLASS
  dt <- copy(micro_dt[, c("SAMPLE_ID", ct), with = FALSE])
  setnames(dt, ct, "FRACTION")
  dt[, GROUP := fifelse(SAMPLE_ID %in% alt_samples, "Altered", "WT")]
  dt <- dt[is.finite(FRACTION)]
  dt[, GROUP := factor(GROUP, levels = c("WT", "Altered"))]

  ggplot(dt, aes(x = GROUP, y = FRACTION, fill = GROUP)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.85) +
    geom_jitter(width = 0.12, height = 0, size = 0.9, alpha = 0.35) +
    scale_fill_manual(values = c(WT = "#9E9E9E", Altered = "#C44E52")) +
    labs(
      title = sprintf("%s · %s", alt, ct),
      subtitle = sprintf(
        "%s | n_alt=%d n_wt=%d | delta=%.2f | adj-p=%.2g FDR=%.2g",
        class_name, row$N_ALTERED, row$N_WT,
        row$CLIFFS_DELTA, row$P_VALUE, row$FDR
      ),
      x = NULL,
      y = "deconvolution fraction/score"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8)
    )
}

# Build quick lookup of altered samples per alteration
all_maps <- c(mut_map, focal_map, arm_map, sv_map)

plot_candidates <- rbindlist(list(
  if (nrow(sig)) sig[order(FDR)][, head(.SD, max_plots_per_class), by = ALTERATION_CLASS] else NULL,
  # Always plot top mutation hits for interpretability
  assoc[ALTERATION_CLASS == "mutation"][order(P_VALUE)][, head(.SD, 20L)]
), use.names = TRUE, fill = TRUE)
if (nrow(plot_candidates)) {
  plot_candidates <- unique(plot_candidates, by = c("ALTERATION_CLASS", "ALTERATION", "CELL_TYPE"))
}

n_plotted <- 0L
for (i in seq_len(nrow(plot_candidates))) {
  row <- plot_candidates[i]
  alt_ids <- all_maps[[row$ALTERATION]]
  if (is.null(alt_ids)) {
    next
  }
  p <- tryCatch(
    plot_one(row, micro, alt_ids),
    error = function(e) NULL
  )
  if (is.null(p)) {
    next
  }
  out <- file.path(
    plot_dir,
    sprintf(
      "%s__%s__%s.pdf",
      sanitize_filename(row$ALTERATION_CLASS),
      sanitize_filename(row$ALTERATION),
      sanitize_filename(row$CELL_TYPE)
    )
  )
  ggsave(out, p, width = 4.2, height = 4.2)
  n_plotted <- n_plotted + 1L
}

# Publication heatmap (+ right bar) is rendered from association_all.tsv
plot_heat_script <- file.path(
  project_dir, "CLLCMA_alpha", "scripts", "plot_alteration_microenv_heatmap.R"
)
status <- system2(
  file.path(R.home("bin"), "Rscript"),
  args = plot_heat_script,
  env = c(
    sprintf("CLL_TME_DIR=%s", project_dir),
    sprintf("DECONV_METHOD=%s", deconv_method)
  )
)
if (!identical(as.integer(status), 0L)) {
  warning("plot_alteration_microenv_heatmap.R failed with status ", status)
}

message("Plotted ", n_plotted, " boxplots -> ", normalizePath(plot_dir))
message("Done -> ", normalizePath(output_dir))
