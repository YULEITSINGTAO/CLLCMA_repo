#!/usr/bin/env Rscript

# Kaplan-Meier OS / FFS for TME-associated genomic alterations.
#
# For each mutation or arm CNV with adaptive FDR < cutoff vs microenvironment,
# split patients into Variated vs Non-variated and plot OS and FFS KM curves
# with log-rank p-values.
#
# Inputs:
#   alteration_microenv_association/<method>/association_significant_fdr05.tsv
#   bulk/cll_broad_2022/data_{mutations,armlevel_cna,clinical_*}.txt
#
# Outputs:
#   alteration_microenv_association/<method>/survival_km/
#     km_stats.tsv
#     km_<alteration>_OS.pdf / _FFS.pdf / _OS_FFS.pdf
#     km_all_alterations_OS.pdf / km_all_alterations_FFS.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(survival)
  library(patchwork)
})

project_dir <- Sys.getenv(
  "CLL_TME_DIR",
  unset = path.expand("~/research/CLL_TME")
)
project_dir <- normalizePath(project_dir, mustWork = TRUE)
setwd(project_dir)

deconv_method <- tolower(Sys.getenv("DECONV_METHOD", unset = "bayesprism"))
fdr_cutoff <- as.numeric(Sys.getenv("KM_FDR_CUTOFF", unset = "0.05"))

bulk_dir <- file.path("CLLCMA_alpha", "bulk", "cll_broad_2022")
assoc_dir <- file.path(
  "CLLCMA_alpha", "results", "deconvolution", "cll_broad_2022",
  "alteration_microenv_association", deconv_method
)
sig_file <- file.path(assoc_dir, "association_significant_fdr05.tsv")
if (!file.exists(sig_file)) {
  stop("Missing association table: ", sig_file)
}
output_dir <- file.path(assoc_dir, "survival_km")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
message("Method: ", deconv_method, " -> ", normalizePath(output_dir, mustWork = FALSE))

silent_like <- c(
  "Silent", "Intron", "3'UTR", "5'UTR", "3'Flank", "5'Flank",
  "IGR", "RNA", "lincRNA", "Targeted_Region"
)

group_colors <- c(
  `Non-variated` = "#4C78A8",
  Variated = "#E45756"
)

# ------------------------------ Helpers --------------------------------------

read_cbioportal_clinical <- function(path) {
  fread(
    path,
    skip = 4L,
    na.strings = c("NA", "", "[Not Available]", "[Not Applicable]")
  )
}

sanitize_filename <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
}

pretty_label <- function(alt, class_name) {
  if (class_name == "arm_cna") {
    x <- sub("^arm_", "", alt)
    x <- gsub("tri\\(12\\)", "+12", x)
    x <- sub("_altered$", "", x)
    return(paste0(x, " CNV"))
  }
  paste0(alt, " mutation")
}

event_from_status <- function(status, endpoint) {
  status <- toupper(trimws(as.character(status)))
  if (endpoint == "OS") {
    return(fifelse(
      is.na(status),
      NA_integer_,
      fifelse(status == "DECEASED", 1L,
        fifelse(status == "LIVING", 0L, NA_integer_)
      )
    ))
  }
  fifelse(
    is.na(status),
    NA_integer_,
    fifelse(grepl("^1(:|$)|FAILURE$", status), 1L,
      fifelse(grepl("^0(:|$)|FAILUREFREE", status), 0L, NA_integer_)
    )
  )
}

make_km_plot <- function(surv_data, title) {
  surv_data <- copy(surv_data)
  surv_data[, GROUP := factor(GROUP, levels = c("Non-variated", "Variated"))]
  fit <- survfit(Surv(TIME, EVENT) ~ GROUP, data = surv_data)
  sdiff <- tryCatch(
    survdiff(Surv(TIME, EVENT) ~ GROUP, data = surv_data),
    error = function(e) NULL
  )
  p_value <- if (!is.null(sdiff)) {
    pchisq(sdiff$chisq, df = length(sdiff$n) - 1L, lower.tail = FALSE)
  } else {
    NA_real_
  }

  strata_names <- gsub(".*=", "", names(fit$strata))
  if (!length(strata_names)) {
    strata_names <- levels(surv_data$GROUP)
  }
  strata_vec <- rep(strata_names, times = as.integer(fit$strata))
  km_dt <- data.table(
    TIME = fit$time,
    SURV = fit$surv,
    GROUP = factor(strata_vec, levels = c("Non-variated", "Variated"))
  )
  zero <- data.table(
    TIME = 0,
    SURV = 1,
    GROUP = factor(
      c("Non-variated", "Variated"),
      levels = c("Non-variated", "Variated")
    )
  )
  km_dt <- rbind(zero, km_dt)
  setorder(km_dt, GROUP, TIME)

  n_tab <- surv_data[, .(N = .N, N_EVENT = sum(EVENT)), by = GROUP]
  setorder(n_tab, GROUP)
  subtitle <- paste0(
    paste(
      sprintf("%s n=%d (ev=%d)", n_tab$GROUP, n_tab$N, n_tab$N_EVENT),
      collapse = "; "
    ),
    sprintf("; log-rank p=%.3g", p_value)
  )

  p <- ggplot(km_dt, aes(x = TIME, y = SURV, color = GROUP)) +
    geom_step(linewidth = 0.9) +
    scale_color_manual(values = group_colors, name = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Months from sampling",
      y = "Survival probability"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

  list(
    plot = p,
    p_value = p_value,
    chisq = if (!is.null(sdiff)) unname(sdiff$chisq) else NA_real_,
    n_tab = n_tab
  )
}

# ----------------------------- Clinical meta ---------------------------------

clinical_sample <- read_cbioportal_clinical(
  file.path(bulk_dir, "data_clinical_sample.txt")
)
clinical_patient <- read_cbioportal_clinical(
  file.path(bulk_dir, "data_clinical_patient.txt")
)
meta <- merge(
  clinical_sample[, .(SAMPLE_ID, PATIENT_ID)],
  clinical_patient[, .(
    PATIENT_ID, OS_STATUS, OS_MONTHS, FFS_STATUS, FFS_MONTHS
  )],
  by = "PATIENT_ID",
  all.x = TRUE
)
meta[, SAMPLE_ID := as.character(SAMPLE_ID)]
meta[, PATIENT_ID := as.character(PATIENT_ID)]
meta[, OS_EVENT := event_from_status(OS_STATUS, "OS")]
meta[, FFS_EVENT := event_from_status(FFS_STATUS, "FFS")]
meta[, OS_MONTHS := suppressWarnings(as.numeric(OS_MONTHS))]
meta[, FFS_MONTHS := suppressWarnings(as.numeric(FFS_MONTHS))]

patient_surv <- meta[
  ,
  .(
    OS_MONTHS = OS_MONTHS[1L],
    OS_EVENT = OS_EVENT[1L],
    FFS_MONTHS = FFS_MONTHS[1L],
    FFS_EVENT = FFS_EVENT[1L]
  ),
  by = PATIENT_ID
]

# ------------------------ TME-associated alterations -------------------------

sig <- fread(sig_file)
alts <- unique(sig[
  (ALTERATION_CLASS == "mutation") |
    (ALTERATION_CLASS == "arm_cna" & grepl("_altered$", ALTERATION)),
  .(ALTERATION_CLASS, ALTERATION, FDR)
][
  ,
  .(MIN_FDR = min(FDR, na.rm = TRUE)),
  by = .(ALTERATION_CLASS, ALTERATION)
][MIN_FDR < fdr_cutoff][order(ALTERATION_CLASS, MIN_FDR)])

if (!nrow(alts)) {
  stop("No mutation/arm alterations with adaptive FDR < ", fdr_cutoff)
}
message(
  "Alterations for KM (adaptive FDR < ", fdr_cutoff, "): ",
  nrow(alts)
)
print(alts)

# ------------------------------- Mutation map --------------------------------

mut <- fread(file.path(bulk_dir, "data_mutations.txt"), skip = 2L)
mut <- mut[
  !is.na(Hugo_Symbol) & Hugo_Symbol != "" &
    !is.na(Tumor_Sample_Barcode) &
    !Variant_Classification %in% silent_like
]
mut[, SAMPLE_ID := as.character(Tumor_Sample_Barcode)]
mut <- merge(mut[, .(SAMPLE_ID, Hugo_Symbol)], meta[, .(SAMPLE_ID, PATIENT_ID)],
             by = "SAMPLE_ID", all = FALSE)
mut_pat <- unique(mut[, .(PATIENT_ID, ALTERATION = Hugo_Symbol)])
mut_assayed_patients <- unique(mut_pat$PATIENT_ID)

# -------------------------------- Arm CNV map --------------------------------

arm <- fread(file.path(bulk_dir, "data_armlevel_cna.txt"))
arm_samples <- setdiff(
  names(arm),
  c("ENTITY_STABLE_ID", "NAME", "DESCRIPTION")
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
arm_long <- merge(
  arm_long,
  meta[, .(SAMPLE_ID, PATIENT_ID)],
  by = "SAMPLE_ID",
  all = FALSE
)
arm_assayed_patients <- unique(arm_long$PATIENT_ID)
arm_var <- arm_long[
  STATUS %in% c("Gain", "Loss"),
  .(PATIENT_ID = unique(PATIENT_ID)),
  by = .(ALTERATION = paste0("arm_", NAME, "_altered"))
]

# --------------------------------- KM loop -----------------------------------

stats_rows <- list()
os_plots <- list()
ffs_plots <- list()

for (i in seq_len(nrow(alts))) {
  row <- alts[i]
  alt <- row$ALTERATION
  cls <- row$ALTERATION_CLASS
  label <- pretty_label(alt, cls)
  message("KM: ", label)

  if (cls == "mutation") {
    assayed <- mut_assayed_patients
    variated <- unique(mut_pat[ALTERATION == alt, PATIENT_ID])
  } else {
    assayed <- arm_assayed_patients
    variated <- unique(arm_var[ALTERATION == alt, PATIENT_ID])
  }
  assayed <- intersect(assayed, patient_surv$PATIENT_ID)
  if (length(assayed) < 20L || length(variated) < 5L) {
    message("  skip (too few patients): n_assay=", length(assayed),
            " n_var=", length(variated))
    next
  }

  dt <- patient_surv[PATIENT_ID %in% assayed]
  dt[, GROUP := fifelse(
    PATIENT_ID %in% variated, "Variated", "Non-variated"
  )]

  for (endpoint in c("OS", "FFS")) {
    time_col <- paste0(endpoint, "_MONTHS")
    event_col <- paste0(endpoint, "_EVENT")
    surv_data <- dt[
      is.finite(get(time_col)) & get(time_col) > 0 &
        !is.na(get(event_col)),
      .(
        PATIENT_ID,
        TIME = get(time_col),
        EVENT = get(event_col),
        GROUP
      )
    ]
    # need both groups
    if (uniqueN(surv_data$GROUP) < 2L || nrow(surv_data) < 20L) {
      message("  skip ", endpoint, ": insufficient complete cases")
      next
    }

    ep_title <- if (endpoint == "OS") {
      "Overall survival"
    } else {
      "Failure-free survival"
    }
    res <- make_km_plot(
      surv_data,
      title = paste0(label, " · ", ep_title)
    )

    n_var <- surv_data[GROUP == "Variated", .N]
    n_wt <- surv_data[GROUP == "Non-variated", .N]
    ev_var <- surv_data[GROUP == "Variated", sum(EVENT)]
    ev_wt <- surv_data[GROUP == "Non-variated", sum(EVENT)]

    stats_rows[[paste(alt, endpoint)]] <- data.table(
      ALTERATION_CLASS = cls,
      ALTERATION = alt,
      LABEL = label,
      ENDPOINT = endpoint,
      N_VARIATED = n_var,
      N_NON_VARIATED = n_wt,
      EVENTS_VARIATED = ev_var,
      EVENTS_NON_VARIATED = ev_wt,
      LOG_RANK_P = res$p_value,
      LOG_RANK_CHISQ = res$chisq,
      TME_MIN_FDR = row$MIN_FDR
    )

    out_stem <- file.path(
      output_dir,
      sprintf("km_%s_%s", sanitize_filename(alt), endpoint)
    )
    ggsave(paste0(out_stem, ".pdf"), res$plot, width = 5.2, height = 4.4)

    if (endpoint == "OS") {
      os_plots[[alt]] <- res$plot
    } else {
      ffs_plots[[alt]] <- res$plot
    }
  }

  if (!is.null(os_plots[[alt]]) && !is.null(ffs_plots[[alt]])) {
    combo <- os_plots[[alt]] + ffs_plots[[alt]] +
      plot_layout(ncol = 2, guides = "collect") &
      theme(legend.position = "bottom")
    ggsave(
      file.path(
        output_dir,
        sprintf("km_%s_OS_FFS.pdf", sanitize_filename(alt))
      ),
      combo,
      width = 10,
      height = 4.6
    )
  }
}

stats <- rbindlist(stats_rows, use.names = TRUE, fill = TRUE)
if (!nrow(stats)) {
  stop("No KM statistics were produced")
}
stats[, FDR_LOG_RANK := p.adjust(LOG_RANK_P, method = "BH"), by = ENDPOINT]
setorder(stats, ENDPOINT, LOG_RANK_P)
fwrite(stats, file.path(output_dir, "km_stats.tsv"), sep = "\t")

# Combined multi-page PDFs
if (length(os_plots)) {
  pdf(file.path(output_dir, "km_all_alterations_OS.pdf"), width = 5.2, height = 4.4)
  for (p in os_plots) print(p)
  dev.off()
}
if (length(ffs_plots)) {
  pdf(file.path(output_dir, "km_all_alterations_FFS.pdf"), width = 5.2, height = 4.4)
  for (p in ffs_plots) print(p)
  dev.off()
}

message("Wrote KM stats for ", nrow(stats), " alteration × endpoint tests")
print(stats[, .(LABEL, ENDPOINT, N_VARIATED, N_NON_VARIATED, LOG_RANK_P, FDR_LOG_RANK)])
message("Done -> ", normalizePath(output_dir))
