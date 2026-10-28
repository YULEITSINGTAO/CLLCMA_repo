#!/usr/bin/env Rscript

# Forest plots of Cox HRs for deconvolution cell types (CIBERSORT / BayesPrism)
# adjusted for IGHV, sex, and age.
#
# For each cell type and endpoint (OS, FFS):
#   Surv ~ scale(fraction) + IGHV_unmutated + SexMale + Age_per_10y
# Clinical-only rows come from:
#   Surv ~ IGHV_unmutated + SexMale + Age_per_10y
#
# Inputs:
#   CLLCMA_alpha/results/deconvolution/cll_broad_2022/{cibersort,bayesprism}_sample_by_feature.tsv
#   CLLCMA_alpha/bulk/cll_broad_2022/data_clinical_{sample,patient}.txt
#
# Outputs (per method):
#   CLLCMA_alpha/results/deconvolution/cll_broad_2022/{method}_forest/
#     forest_OS.pdf
#     forest_FFS.pdf
#     forest_stats.tsv
# Combined:
#   .../forest_stats_cibersort_bayesprism.tsv
#   .../forest_OS_cibersort_bayesprism.pdf
#   .../forest_FFS_cibersort_bayesprism.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(survival)
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
    output_dir = file.path(deconv_dir, "cibersort_forest"),
    label = "CIBERSORT"
  ),
  bayesprism = list(
    input = file.path(deconv_dir, "bayesprism_sample_by_feature.tsv"),
    output_dir = file.path(deconv_dir, "bayesprism_forest"),
    label = "BayesPrism"
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
      fifelse(grepl("^0(:|$)|FAILUREFREE$", status), 0L, NA_integer_)
    )
  )
}

map_ighv <- function(x) {
  x <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(x))
  out[grepl("unmut", x) | x %in% c("u", "u-cll", "unmutated ighv")] <- "unmutated"
  out[is.na(out) & (grepl("^mut", x) | x %in% c("m", "m-cll", "mutated", "mutated ighv"))] <- "mutated"
  factor(out, levels = c("mutated", "unmutated"))
}

map_sex <- function(x) {
  x <- trimws(as.character(x))
  out <- rep(NA_character_, length(x))
  out[x %in% c("Female", "F", "female")] <- "Female"
  out[x %in% c("Male", "M", "male")] <- "Male"
  factor(out, levels = c("Female", "Male"))
}

safe_cox <- function(form, data) {
  tryCatch(
    coxph(form, data = data, ties = "efron", model = FALSE, x = FALSE, y = FALSE),
    error = function(e) NULL
  )
}

cox_row <- function(fit, term, label, group) {
  if (is.null(fit) || !(term %in% names(coef(fit)))) {
    return(data.table(
      GROUP = group,
      TERM = term,
      LABEL = label,
      HR = NA_real_,
      HR_LOW = NA_real_,
      HR_HIGH = NA_real_,
      P = NA_real_,
      N = NA_integer_,
      N_EVENT = NA_integer_
    ))
  }
  s <- summary(fit)
  ci <- confint(fit)
  data.table(
    GROUP = group,
    TERM = term,
    LABEL = label,
    HR = unname(exp(coef(fit)[term])),
    HR_LOW = unname(exp(ci[term, 1])),
    HR_HIGH = unname(exp(ci[term, 2])),
    P = unname(s$coefficients[term, "Pr(>|z|)"]),
    N = as.integer(fit$n),
    N_EVENT = as.integer(fit$nevent)
  )
}

plot_forest <- function(stats, title, file_path) {
  dt <- copy(stats)
  dt <- dt[is.finite(HR) & is.finite(HR_LOW) & is.finite(HR_HIGH)]
  if (!nrow(dt)) {
    warning("No finite HRs for ", title)
    return(invisible(NULL))
  }

  # Order: clinical first, then cell types by HR
  clin <- dt[GROUP == "Clinical"]
  cells <- dt[GROUP == "Cell type"]
  setorder(cells, HR)
  dt <- rbind(clin, cells)
  dt[, LABEL := factor(LABEL, levels = rev(unique(LABEL)))]
  dt[, SIG := fifelse(!is.na(P) & P < 0.05, "p < 0.05", "NS")]

  # Cap display CI for plotting only (keep table unclipped)
  dt[, HR_LOW_PLOT := pmax(HR_LOW, 0.05)]
  dt[, HR_HIGH_PLOT := pmin(HR_HIGH, 20)]
  dt[, HR_PLOT := pmin(pmax(HR, 0.05), 20)]

  p <- ggplot(dt, aes(x = HR_PLOT, y = LABEL, color = SIG)) +
    geom_vline(xintercept = 1, linetype = 2, linewidth = 0.4, colour = "grey40") +
    geom_errorbar(
      aes(xmin = HR_LOW_PLOT, xmax = HR_HIGH_PLOT),
      orientation = "y",
      height = 0.25,
      linewidth = 0.55
    ) +
    geom_point(size = 2.2) +
    scale_x_log10(
      breaks = c(0.25, 0.5, 1, 2, 4, 8),
      labels = c("0.25", "0.5", "1", "2", "4", "8")
    ) +
    scale_color_manual(
      values = c("p < 0.05" = "#D55E00", NS = "#0072B2"),
      name = NULL
    ) +
    labs(
      title = title,
      subtitle = paste0(
        "Cell types: Cox adjusted for IGHV + sex + age (per SD of fraction). ",
        "Clinical: model without cell type. HR > 1 = higher risk."
      ),
      x = "Hazard ratio (95% CI, log scale)",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 8, colour = "grey30"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )

  # Facet clinical vs cell type for readability
  if (uniqueN(dt$GROUP) > 1L) {
    p <- p + facet_grid(GROUP ~ ., scales = "free_y", space = "free_y")
  }

  n_rows <- nrow(dt)
  h <- max(5, 0.28 * n_rows + 2.2)
  ggsave(file_path, p, width = 8.5, height = h, limitsize = FALSE)
  invisible(p)
}

build_patient_table <- function(frac_file, meta) {
  cs <- fread(frac_file)
  feature_cols <- setdiff(names(cs), "SAMPLE_ID")
  dat <- merge(cs, meta, by = "SAMPLE_ID", all = FALSE)
  if (!nrow(dat)) {
    return(NULL)
  }

  patient <- dat[
    ,
    c(
      .(
        N_SAMPLES = .N,
        AGE = {
          a <- suppressWarnings(as.numeric(AGE_SAMPLING[1L]))
          if (!is.finite(a)) {
            a <- suppressWarnings(as.numeric(AGE[1L]))
          }
          a
        },
        SEX = as.character(SEX[1L]),
        IGHV_MUTATION_STATUS = as.character(IGHV_MUTATION_STATUS[1L]),
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
  patient[, IGHV := map_ighv(IGHV_MUTATION_STATUS)]
  patient[, SEX := map_sex(SEX)]
  patient[, AGE10 := AGE / 10]
  list(patient = patient, feature_cols = feature_cols)
}

fit_endpoint_stats <- function(patient, feature_cols, endpoint, method, label) {
  time_col <- paste0(endpoint, "_MONTHS")
  status_col <- paste0(endpoint, "_STATUS")

  df <- copy(patient)
  df[, TIME := suppressWarnings(as.numeric(get(time_col)))]
  df[, EVENT := event_from_status(get(status_col), endpoint)]
  df <- df[
    is.finite(TIME) & TIME > 0 & !is.na(EVENT) &
      !is.na(IGHV) & !is.na(SEX) & is.finite(AGE10)
  ]

  rows <- list()

  # Clinical-only model
  fit_clin <- safe_cox(Surv(TIME, EVENT) ~ IGHV + SEX + AGE10, data = df)
  rows[["clin_ighv"]] <- cox_row(
    fit_clin, "IGHVunmutated", "IGHV unmutated (vs mutated)", "Clinical"
  )
  rows[["clin_sex"]] <- cox_row(
    fit_clin, "SEXMale", "Sex Male (vs Female)", "Clinical"
  )
  rows[["clin_age"]] <- cox_row(
    fit_clin, "AGE10", "Age (per 10 years)", "Clinical"
  )

  for (feature in feature_cols) {
    dfi <- copy(df)
    dfi[, CELL := as.numeric(get(feature))]
    if (length(unique(dfi$CELL[is.finite(dfi$CELL)])) < 2L) {
      next
    }
    dfi[, CELL_Z := as.numeric(scale(CELL))]
    fit <- safe_cox(
      Surv(TIME, EVENT) ~ CELL_Z + IGHV + SEX + AGE10,
      data = dfi
    )
    rows[[feature]] <- cox_row(
      fit,
      "CELL_Z",
      paste0(feature, " (per SD)"),
      "Cell type"
    )
  }

  out <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  out[, METHOD := method]
  out[, METHOD_LABEL := label]
  out[, ENDPOINT := endpoint]
  out
}

# ------------------------------- Clinical ------------------------------------

clinical_sample <- read_cbioportal_clinical(clinical_sample_file)
clinical_patient <- read_cbioportal_clinical(clinical_patient_file)

need_sample <- c("SAMPLE_ID", "PATIENT_ID")
need_patient <- c(
  "PATIENT_ID", "AGE", "AGE_SAMPLING", "SEX", "IGHV_MUTATION_STATUS",
  "OS_STATUS", "OS_MONTHS", "FFS_STATUS", "FFS_MONTHS"
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

all_stats <- list()

for (method in names(methods)) {
  cfg <- methods[[method]]
  if (!file.exists(cfg$input)) {
    warning("Skipping ", method, ": missing ", cfg$input)
    next
  }

  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
  message("=== ", method, " ===")

  built <- build_patient_table(cfg$input, meta)
  if (is.null(built)) {
    warning("Skipping ", method, ": no SAMPLE_ID overlap")
    next
  }
  patient <- built$patient
  feature_cols <- built$feature_cols
  message(
    "Patients: ", nrow(patient),
    "; features: ", length(feature_cols),
    "; with complete covariates: ",
    patient[!is.na(IGHV) & !is.na(SEX) & is.finite(AGE10), .N]
  )

  method_stats <- rbindlist(lapply(c("OS", "FFS"), function(endpoint) {
    fit_endpoint_stats(patient, feature_cols, endpoint, method, cfg$label)
  }), use.names = TRUE, fill = TRUE)

  method_stats[, FDR_CELL := NA_real_]
  method_stats[
    GROUP == "Cell type",
    FDR_CELL := p.adjust(P, method = "BH"),
    by = ENDPOINT
  ]
  setorder(method_stats, ENDPOINT, GROUP, HR)
  fwrite(method_stats, file.path(cfg$output_dir, "forest_stats.tsv"), sep = "\t")

  for (endpoint in c("OS", "FFS")) {
    plot_forest(
      method_stats[ENDPOINT == endpoint],
      title = paste0(cfg$label, " · ", endpoint, " multivariable Cox forest"),
      file_path = file.path(cfg$output_dir, paste0("forest_", endpoint, ".pdf"))
    )
  }

  all_stats[[method]] <- method_stats
  message("Wrote forest plots to ", normalizePath(cfg$output_dir))
}

combined <- rbindlist(all_stats, use.names = TRUE, fill = TRUE)
if (nrow(combined)) {
  fwrite(
    combined,
    file.path(deconv_dir, "forest_stats_cibersort_bayesprism.tsv"),
    sep = "\t"
  )

  # Side-by-side method comparison forests (cell types + clinical)
  for (endpoint in c("OS", "FFS")) {
    dt <- combined[ENDPOINT == endpoint & is.finite(HR)]
    if (!nrow(dt)) {
      next
    }
    clin <- dt[GROUP == "Clinical"][, .SD[1], by = LABEL]
    clin[, METHOD_LABEL := "Clinical"]
    cells <- dt[GROUP == "Cell type"]
    plot_dt <- rbind(
      clin[, .(METHOD_LABEL, GROUP, LABEL, HR, HR_LOW, HR_HIGH, P)],
      cells[, .(METHOD_LABEL, GROUP, LABEL, HR, HR_LOW, HR_HIGH, P)],
      fill = TRUE
    )
    plot_dt[, SIG := fifelse(!is.na(P) & P < 0.05, "p < 0.05", "NS")]
    plot_dt[, HR_LOW_PLOT := pmax(HR_LOW, 0.05)]
    plot_dt[, HR_HIGH_PLOT := pmin(HR_HIGH, 20)]
    plot_dt[, HR_PLOT := pmin(pmax(HR, 0.05), 20)]

    # Shared label order: clinical first, then cell types by mean HR across methods
    cell_order <- cells[
      ,
      .(mean_hr = exp(mean(log(HR), na.rm = TRUE))),
      by = LABEL
    ]
    setorder(cell_order, mean_hr)
    lab_levels <- c(
      rev(unique(as.character(clin$LABEL))),
      rev(cell_order$LABEL)
    )
    plot_dt[, LABEL := factor(LABEL, levels = lab_levels)]
    plot_dt[, FACET := fifelse(GROUP == "Clinical", "Clinical", as.character(METHOD_LABEL))]
    plot_dt[, FACET := factor(
      FACET,
      levels = c("Clinical", "CIBERSORT", "BayesPrism")
    )]

    p <- ggplot(plot_dt, aes(x = HR_PLOT, y = LABEL, color = SIG)) +
      geom_vline(xintercept = 1, linetype = 2, linewidth = 0.4, colour = "grey40") +
      geom_errorbar(
        aes(xmin = HR_LOW_PLOT, xmax = HR_HIGH_PLOT),
        orientation = "y",
        height = 0.25,
        linewidth = 0.5
      ) +
      geom_point(size = 1.8) +
      facet_wrap(~ FACET, nrow = 1) +
      scale_x_log10(
        breaks = c(0.25, 0.5, 1, 2, 4, 8),
        labels = c("0.25", "0.5", "1", "2", "4", "8")
      ) +
      scale_color_manual(
        values = c("p < 0.05" = "#D55E00", NS = "#0072B2"),
        name = NULL
      ) +
      labs(
        title = paste0("Deconvolution Cox forest · ", endpoint),
        subtitle = "Cell-type HR per SD, adjusted for IGHV + sex + age",
        x = "Hazard ratio (95% CI, log scale)",
        y = NULL
      ) +
      theme_bw(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey95", colour = NA)
      )

    ggsave(
      file.path(deconv_dir, paste0("forest_", endpoint, "_cibersort_bayesprism.pdf")),
      p,
      width = 12,
      height = max(6, 0.22 * uniqueN(plot_dt$LABEL) + 2.5),
      limitsize = FALSE
    )
  }

  message(
    "Combined stats -> ",
    normalizePath(file.path(deconv_dir, "forest_stats_cibersort_bayesprism.tsv"))
  )
}

message("Done")
