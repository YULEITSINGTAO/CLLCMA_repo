#!/usr/bin/env Rscript

# Median-split Kaplan–Meier survival for CLLCMA deconvolution fractions
# (CIBERSORT and BayesPrism).
#
# For each feature, split patients into High / Low by the feature median,
# then plot OS and FFS KM curves with log-rank p-values.
#
# Inputs:
#   CLLCMA_alpha/results/deconvolution/cll_broad_2022/{cibersort,bayesprism}_sample_by_feature.tsv
#   CLLCMA_alpha/bulk/cll_broad_2022/data_clinical_{sample,patient}.txt
#
# Outputs (per method):
#   CLLCMA_alpha/results/deconvolution/cll_broad_2022/{cibersort,bayesprism}_km/
#     km_stats.tsv
#     km_<feature>_OS.pdf
#     km_<feature>_FFS.pdf
#     km_<feature>_OS_FFS.pdf
#     km_all_features_OS.pdf
#     km_all_features_FFS.pdf

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

deconv_dir <- file.path(
  "CLLCMA_alpha", "results", "deconvolution", "cll_broad_2022"
)
clinical_sample_file <- file.path(
  "CLLCMA_alpha", "bulk", "cll_broad_2022", "data_clinical_sample.txt"
)
clinical_patient_file <- file.path(
  "CLLCMA_alpha", "bulk", "cll_broad_2022", "data_clinical_patient.txt"
)

methods <- list(
  cibersort = list(
    input = file.path(deconv_dir, "cibersort_sample_by_feature.tsv"),
    output_dir = file.path(deconv_dir, "cibersort_km"),
    label = "CIBERSORT"
  ),
  bayesprism = list(
    input = file.path(deconv_dir, "bayesprism_sample_by_feature.tsv"),
    output_dir = file.path(deconv_dir, "bayesprism_km"),
    label = "BayesPrism"
  )
)

group_colors <- c(
  Low = "#4C78A8",
  High = "#E45756"
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
  surv_data[, GROUP := factor(GROUP, levels = c("Low", "High"))]
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
    GROUP = factor(strata_vec, levels = c("Low", "High"))
  )
  zero <- data.table(
    TIME = 0,
    SURV = 1,
    GROUP = factor(c("Low", "High"), levels = c("Low", "High"))
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
    n_tab = n_tab,
    median_surv = summary(fit)$table
  )
}

run_method_km <- function(method, cfg, meta) {
  if (!file.exists(cfg$input)) {
    warning("Skipping ", method, ": missing ", cfg$input)
    return(invisible(NULL))
  }

  output_dir <- cfg$output_dir
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("=== ", method, " ===")

  cs <- fread(cfg$input)
  feature_cols <- setdiff(names(cs), "SAMPLE_ID")
  if (!length(feature_cols)) {
    warning("Skipping ", method, ": no features found")
    return(invisible(NULL))
  }

  dat <- merge(cs, meta, by = "SAMPLE_ID", all = FALSE)
  if (!nrow(dat)) {
    warning("Skipping ", method, ": no SAMPLE_ID overlap with clinical")
    return(invisible(NULL))
  }

  # One row per patient: average fractions across RNA samples.
  patient_frac <- dat[
    ,
    c(
      .(
        N_SAMPLES = .N,
        OS_STATUS = OS_STATUS[1L],
        OS_MONTHS = OS_MONTHS[1L],
        FFS_STATUS = FFS_STATUS[1L],
        FFS_MONTHS = FFS_MONTHS[1L]
      ),
      lapply(.SD, function(x) mean(as.numeric(x), na.rm = TRUE))
    ),
    by = PATIENT_ID,
    .SDcols = feature_cols
  ]

  message(
    "Patients with ", cfg$label, " + clinical: ", nrow(patient_frac),
    " (from ", nrow(dat), " samples); features: ", length(feature_cols)
  )

  stats_rows <- list()
  os_plots <- list()
  ffs_plots <- list()

  for (feature in feature_cols) {
    values <- patient_frac[[feature]]
    med <- median(values, na.rm = TRUE)
    if (!is.finite(med) || length(unique(values[is.finite(values)])) < 2L) {
      message("Skipping ", feature, ": no variation")
      next
    }

    # Ties at the median: assign == median to Low so groups are exhaustive.
    patient_frac[, GROUP := fifelse(get(feature) > med, "High", "Low")]

    for (endpoint in c("OS", "FFS")) {
      time_col <- paste0(endpoint, "_MONTHS")
      status_col <- paste0(endpoint, "_STATUS")

      surv_data <- patient_frac[
        ,
        .(
          PATIENT_ID,
          FEATURE = feature,
          FEATURE_VALUE = get(feature),
          GROUP,
          MEDIAN_CUTOFF = med,
          TIME = suppressWarnings(as.numeric(get(time_col))),
          EVENT = event_from_status(get(status_col), endpoint)
        )
      ]
      surv_data <- surv_data[
        is.finite(TIME) & TIME > 0 & !is.na(EVENT) & !is.na(GROUP)
      ]

      if (nrow(surv_data) < 20L || uniqueN(surv_data$GROUP) < 2L) {
        message("Skipping ", feature, " ", endpoint, ": insufficient data")
        next
      }

      title <- sprintf(
        "%s %s · %s (median=%.4g)",
        cfg$label, feature, endpoint, med
      )
      km <- make_km_plot(surv_data, title)

      out_pdf <- file.path(
        output_dir,
        sprintf("km_%s_%s.pdf", sanitize_filename(feature), endpoint)
      )
      ggsave(out_pdf, km$plot, width = 5.5, height = 4.5)

      n_low <- km$n_tab[GROUP == "Low", N]
      n_high <- km$n_tab[GROUP == "High", N]
      ev_low <- km$n_tab[GROUP == "Low", N_EVENT]
      ev_high <- km$n_tab[GROUP == "High", N_EVENT]
      if (!length(n_low)) n_low <- 0L
      if (!length(n_high)) n_high <- 0L
      if (!length(ev_low)) ev_low <- 0L
      if (!length(ev_high)) ev_high <- 0L

      stats_rows[[paste(feature, endpoint, sep = "_")]] <- data.table(
        METHOD = method,
        FEATURE = feature,
        ENDPOINT = endpoint,
        MEDIAN_CUTOFF = med,
        N_PATIENTS = nrow(surv_data),
        N_LOW = as.integer(n_low),
        N_HIGH = as.integer(n_high),
        EVENTS_LOW = as.integer(ev_low),
        EVENTS_HIGH = as.integer(ev_high),
        LOG_RANK_CHISQ = km$chisq,
        LOG_RANK_P = km$p_value
      )

      if (endpoint == "OS") {
        os_plots[[feature]] <- km$plot
      } else {
        ffs_plots[[feature]] <- km$plot
      }

      message(
        feature, " ", endpoint, ": n=", nrow(surv_data),
        " median=", signif(med, 4), " p=", signif(km$p_value, 3)
      )
    }

    if (!is.null(os_plots[[feature]]) && !is.null(ffs_plots[[feature]])) {
      combo <- os_plots[[feature]] + ffs_plots[[feature]] +
        plot_layout(ncol = 2) +
        plot_annotation(
          title = paste0(cfg$label, " median-split KM: ", feature)
        )
      ggsave(
        file.path(
          output_dir,
          sprintf("km_%s_OS_FFS.pdf", sanitize_filename(feature))
        ),
        combo,
        width = 10,
        height = 4.5
      )
    }
  }

  stats <- rbindlist(stats_rows, use.names = TRUE, fill = TRUE)
  if (!nrow(stats)) {
    warning("No KM statistics produced for ", method)
    return(invisible(NULL))
  }

  stats[, FDR := p.adjust(LOG_RANK_P, method = "BH"), by = ENDPOINT]
  setorder(stats, ENDPOINT, LOG_RANK_P)
  fwrite(stats, file.path(output_dir, "km_stats.tsv"), sep = "\t")

  if (length(os_plots)) {
    os_grid <- wrap_plots(os_plots, ncol = 3) +
      plot_annotation(
        title = paste0(cfg$label, " median-split KM · Overall survival")
      )
    ggsave(
      file.path(output_dir, "km_all_features_OS.pdf"),
      os_grid,
      width = 12,
      height = 4 * ceiling(length(os_plots) / 3)
    )
  }
  if (length(ffs_plots)) {
    ffs_grid <- wrap_plots(ffs_plots, ncol = 3) +
      plot_annotation(
        title = paste0(cfg$label, " median-split KM · Failure-free survival")
      )
    ggsave(
      file.path(output_dir, "km_all_features_FFS.pdf"),
      ffs_grid,
      width = 12,
      height = 4 * ceiling(length(ffs_plots) / 3)
    )
  }

  message("Wrote KM results to ", normalizePath(output_dir))
  print(stats[1:min(10L, .N)])
  invisible(stats)
}

# ------------------------------- Load clinical -------------------------------

clinical_sample <- read_cbioportal_clinical(clinical_sample_file)
clinical_patient <- read_cbioportal_clinical(clinical_patient_file)

need_sample <- c("SAMPLE_ID", "PATIENT_ID")
need_patient <- c(
  "PATIENT_ID", "OS_STATUS", "OS_MONTHS", "FFS_STATUS", "FFS_MONTHS"
)
stopifnot(
  all(need_sample %in% names(clinical_sample)),
  all(need_patient %in% names(clinical_patient))
)

meta <- merge(
  clinical_sample[, ..need_sample],
  clinical_patient[, ..need_patient],
  by = "PATIENT_ID",
  all.x = TRUE
)

# ------------------------------ Per method ----------------------------------

all_stats <- lapply(names(methods), function(method) {
  run_method_km(method, methods[[method]], meta)
})
all_stats <- rbindlist(all_stats, use.names = TRUE, fill = TRUE)

if (nrow(all_stats)) {
  fwrite(
    all_stats,
    file.path(deconv_dir, "km_stats_cibersort_bayesprism.tsv"),
    sep = "\t"
  )
  message(
    "Combined stats -> ",
    normalizePath(file.path(deconv_dir, "km_stats_cibersort_bayesprism.tsv"))
  )
}

message("Done")
