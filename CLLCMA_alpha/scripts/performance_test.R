#!/usr/bin/env Rscript

# performance_test.R
#
# Benchmark deconvolution methods for predicting OS and FFS.
# Features = cell-type / score columns from each method. Evaluate time-dependent
# AUC at 1, 3, and 5 years using 10-fold CV with Cox, Lasso Cox, Cox ridge,
# random survival forest (ranger), and DeepSurv (Cox-PH MLP). Folds are fixed
# once per endpoint on the common patient set and reused across methods/learners.
#
# Inputs:
#   CLLCMA_alpha/results/deconvolution/cll_broad_2022/all_methods_long.tsv
#   CLLCMA_alpha/bulk/cll_broad_2022/data_clinical_{sample,patient}.txt
#
# Outputs:
#   CLLCMA_alpha/results/deconvolution/cll_broad_2022/deconvolution_method_benchmark/
#   td AUC, Harrell C, Uno C_tau, integrated AUC, IPCW Brier / IBS / IPA
#
# Set METRICS_FROM_OOF=1 to recompute metrics from oof_risk_scores.tsv without CV.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(survival)
  library(glmnet)
  library(ranger)
})

project_dir <- Sys.getenv(
  "CLL_TME_DIR",
  unset = path.expand("~/research/CLL_TME")
)
project_dir <- normalizePath(project_dir, mustWork = TRUE)
setwd(project_dir)

estimate_file <- file.path(
  "CLLCMA_alpha", "results", "deconvolution", "cll_broad_2022",
  "all_methods_long.tsv"
)
clinical_sample_file <- file.path(
  "CLLCMA_alpha", "bulk", "cll_broad_2022", "data_clinical_sample.txt"
)
clinical_patient_file <- file.path(
  "CLLCMA_alpha", "bulk", "cll_broad_2022", "data_clinical_patient.txt"
)
output_dir <- file.path(
  "CLLCMA_alpha", "results", "deconvolution", "cll_broad_2022",
  "deconvolution_method_benchmark"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_folds <- 10L
horizons_months <- c(12, 36, 60)
horizon_labels <- c("1-year", "3-year", "5-year")
names(horizon_labels) <- as.character(horizons_months)
set.seed(123L)

# Learners compared on the same features + fixed folds.
# Override with LEARNERS=cox,cox_lasso (comma-separated) to run a subset.
all_learners <- c(
  "cox", "cox_lasso", "cox_ridge", "random_forest", "deepsurv"
)
learner_order <- all_learners
learners_env <- Sys.getenv("LEARNERS", unset = "")
if (nzchar(learners_env)) {
  requested <- trimws(strsplit(learners_env, ",", fixed = TRUE)[[1]])
  learner_order <- intersect(all_learners, requested)
  if (!length(learner_order)) {
    stop("LEARNERS env selected no valid learners")
  }
  message("Running learners: ", paste(learner_order, collapse = ", "))
}
rf_num_trees <- 500L
deepsurv_python <- Sys.getenv(
  "DEEPSURV_PYTHON",
  unset = path.expand("~/miniconda3/envs/singlecell/bin/python")
)
deepsurv_script <- file.path(
  project_dir, "CLLCMA_alpha", "scripts", "deepsurv_cv.py"
)

method_order <- c(
  "bayesprism", "cibersort", "quantiseq", "mcp_counter", "xcell"
)

# xCell composite scores kept; they can still be predictive.
xcell_drop <- character()

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
      is.na(status), NA_integer_,
      fifelse(status == "DECEASED", 1L,
        fifelse(status == "LIVING", 0L, NA_integer_))
    ))
  }
  fifelse(
    is.na(status), NA_integer_,
    fifelse(grepl("^1(:|$)|FAILURE$", status), 1L,
      fifelse(grepl("^0(:|$)|FAILUREFREE", status), 0L, NA_integer_))
  )
}

zscore_matrix <- function(x_train, x_test) {
  mu <- colMeans(x_train, na.rm = TRUE)
  sdv <- apply(x_train, 2L, sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  train_z <- sweep(sweep(x_train, 2L, mu, "-"), 2L, sdv, "/")
  test_z <- sweep(sweep(x_test, 2L, mu, "-"), 2L, sdv, "/")
  train_z[!is.finite(train_z)] <- 0
  test_z[!is.finite(test_z)] <- 0
  list(train = train_z, test = test_z)
}

# Censoring distribution G(t): KM treating event as censoring indicator flipped.
censoring_surv <- function(time, event) {
  # Censoring indicator: 1 if censored
  cens <- as.integer(event == 0L)
  # If no censored observations, G=1
  if (sum(cens) == 0L) {
    return(function(t) rep(1, length(t)))
  }
  fit <- survfit(Surv(time, cens) ~ 1)
  function(t) {
    # Left-continuous G(t-): use step function just before t
    out <- sapply(t, function(tt) {
      if (!is.finite(tt) || tt <= 0) {
        return(1)
      }
      idx <- max(which(fit$time < tt), 0L)
      if (idx == 0L) 1 else as.numeric(fit$surv[idx])
    })
    pmax(out, 1e-6)
  }
}

# Cumulative/dynamic time-dependent AUC with IPCW (Uno-style weights).
td_auc_ipcw <- function(time, event, score, t_horizon) {
  ok <- is.finite(time) & is.finite(event) & is.finite(score) & time > 0
  time <- time[ok]
  event <- as.integer(event[ok])
  score <- as.numeric(score[ok])
  n <- length(time)
  if (n < 20L) {
    return(NA_real_)
  }

  case <- event == 1L & time <= t_horizon
  control <- time > t_horizon
  if (sum(case) < 5L || sum(control) < 5L) {
    return(NA_real_)
  }

  G <- censoring_surv(time, event)
  # Weights: cases 1/G(T_i), controls 1/G(t)
  w_case <- 1 / G(time[case])
  w_ctrl <- rep(1 / G(t_horizon), sum(control))
  s_case <- score[case]
  s_ctrl <- score[control]

  # Weighted Mann-Whitney: P(score_case > score_ctrl)
  # Efficient double loop avoided via ranking approximation with weights
  auc_num <- 0
  auc_den <- 0
  for (i in seq_along(s_case)) {
    gt <- as.numeric(s_case[i] > s_ctrl)
    eq <- as.numeric(s_case[i] == s_ctrl)
    contrib <- sum(w_ctrl * (gt + 0.5 * eq))
    auc_num <- auc_num + w_case[i] * contrib
    auc_den <- auc_den + w_case[i] * sum(w_ctrl)
  }
  if (auc_den <= 0) {
    return(NA_real_)
  }
  as.numeric(auc_num / auc_den)
}

# Uno (2011) IPCW C-index. Higher score = higher risk. tau in the same time unit.
uno_cindex <- function(time, event, score, tau = Inf) {
  ok <- is.finite(time) & is.finite(event) & is.finite(score) & time > 0
  time <- time[ok]
  event <- as.integer(event[ok])
  score <- as.numeric(score[ok])
  n <- length(time)
  if (n < 20L || sum(event == 1L) < 5L) {
    return(NA_real_)
  }
  G <- censoring_surv(time, event)
  num <- 0
  den <- 0
  for (i in seq_len(n)) {
    if (event[i] != 1L || time[i] >= tau) {
      next
    }
    w <- 1 / (G(time[i])^2)
    later <- time > time[i]
    n_later <- sum(later)
    if (n_later < 1L) {
      next
    }
    s_i <- score[i]
    s_j <- score[later]
    den <- den + w * n_later
    num <- num + w * sum((s_i > s_j) + 0.5 * (s_i == s_j))
  }
  if (den <= 0) {
    return(NA_real_)
  }
  as.numeric(num / den)
}

harrell_cindex <- function(time, event, score) {
  # survival::concordance treats larger x as longer survival; flip risk scores.
  cfit <- tryCatch(
    concordance(Surv(time, event) ~ I(-score)),
    error = function(e) NULL
  )
  if (is.null(cfit)) {
    return(list(c = NA_real_, se = NA_real_))
  }
  list(
    c = as.numeric(cfit$concordance),
    se = as.numeric(sqrt(cfit$var))
  )
}

km_surv_fun <- function(time, event) {
  fit <- survfit(Surv(time, event) ~ 1)
  function(t) {
    vapply(t, function(tt) {
      if (!is.finite(tt) || tt <= 0) {
        return(1)
      }
      idx <- max(which(fit$time <= tt), 0L)
      if (idx == 0L) 1 else as.numeric(fit$surv[idx])
    }, numeric(1))
  }
}

# KM-weighted and uniform integrated C/D AUC on [t_min, t_max].
integrated_auc <- function(
    time, event, score, t_min = 12, t_max = 60, step = 3) {
  grid <- seq(t_min, t_max, by = step)
  aucs <- vapply(grid, function(tt) {
    td_auc_ipcw(time, event, score, tt)
  }, numeric(1))
  ok <- is.finite(aucs)
  if (sum(ok) < 3L) {
    return(list(iauc_km = NA_real_, iauc_uniform = NA_real_))
  }
  g <- grid[ok]
  a <- aucs[ok]
  iauc_uniform <- sum(diff(g) * (a[-1] + a[-length(a)]) / 2) / (max(g) - min(g))
  S <- km_surv_fun(time, event)
  sv <- S(g)
  w <- -diff(sv)
  w[w < 0] <- 0
  iauc_km <- if (sum(w) <= 0) {
    iauc_uniform
  } else {
    sum(((a[-1] + a[-length(a)]) / 2) * w) / sum(w)
  }
  list(iauc_km = as.numeric(iauc_km), iauc_uniform = as.numeric(iauc_uniform))
}

breslow_cumhaz <- function(time, event, lp) {
  ok <- is.finite(time) & is.finite(event) & is.finite(lp) & time > 0
  time <- time[ok]
  event <- as.integer(event[ok])
  rsk <- exp(pmin(pmax(lp[ok], -20), 20))
  evt_times <- sort(unique(time[event == 1L]))
  if (!length(evt_times)) {
    return(list(time = numeric(), H = numeric()))
  }
  H <- numeric(length(evt_times))
  h <- 0
  for (k in seq_along(evt_times)) {
    tk <- evt_times[k]
    d <- sum(event == 1L & abs(time - tk) < 1e-12)
    den <- sum(rsk[time >= tk - 1e-12])
    h <- h + d / max(den, 1e-12)
    H[k] <- h
  }
  list(time = evt_times, H = H)
}

H0_at <- function(bh, t) {
  if (!length(bh$time)) {
    return(rep(0, length(t)))
  }
  vapply(t, function(tt) {
    idx <- max(which(bh$time <= tt), 0L)
    if (idx == 0L) 0 else bh$H[idx]
  }, numeric(1))
}

# Fold-wise Breslow Ŝ(t) from OOF linear predictors (higher = higher risk).
oof_surv_matrix <- function(dt, times) {
  m <- matrix(NA_real_, nrow(dt), length(times))
  folds <- sort(unique(dt$FOLD))
  for (f in folds) {
    train <- dt$FOLD != f
    test <- which(dt$FOLD == f)
    if (sum(train) < 20L || !length(test)) {
      next
    }
    lp_c <- mean(dt$RISK_SCORE[train], na.rm = TRUE)
    bh <- breslow_cumhaz(
      dt$TIME[train],
      dt$EVENT[train],
      dt$RISK_SCORE[train] - lp_c
    )
    ht <- H0_at(bh, times)
    lp_te <- dt$RISK_SCORE[test] - lp_c
    m[test, ] <- exp(-outer(exp(pmin(pmax(lp_te, -20), 20)), ht, `*`))
  }
  pmin(pmax(m, 1e-6), 1)
}

oof_km_null_matrix <- function(dt, times) {
  m <- matrix(NA_real_, nrow(dt), length(times))
  folds <- sort(unique(dt$FOLD))
  for (f in folds) {
    train <- dt$FOLD != f
    test <- which(dt$FOLD == f)
    if (sum(train) < 20L || !length(test)) {
      next
    }
    S <- km_surv_fun(dt$TIME[train], dt$EVENT[train])
    sv <- S(times)
    m[test, ] <- matrix(sv, nrow = length(test), ncol = length(times), byrow = TRUE)
  }
  pmin(pmax(m, 1e-6), 1)
}

brier_ipcw_from_surv <- function(time, event, surv_hat, t_horizon) {
  ok <- is.finite(time) & is.finite(event) & is.finite(surv_hat) & time > 0
  time <- time[ok]
  event <- as.integer(event[ok])
  surv_hat <- pmin(pmax(surv_hat[ok], 0), 1)
  n <- length(time)
  if (n < 20L) {
    return(NA_real_)
  }
  G <- censoring_surv(time, event)
  contrib <- numeric(n)
  evented <- event == 1L & time <= t_horizon
  alive <- time > t_horizon
  if (any(evented)) {
    contrib[evented] <- (0 - surv_hat[evented])^2 / G(time[evented])
  }
  if (any(alive)) {
    contrib[alive] <- (1 - surv_hat[alive])^2 / G(t_horizon)
  }
  mean(contrib)
}

ibs_from_brier_curve <- function(times, bs, tau) {
  ok <- is.finite(times) & is.finite(bs) & times > 0 & times <= tau
  if (sum(ok) < 2L) {
    return(NA_real_)
  }
  t <- c(0, times[ok])
  y <- c(0, bs[ok])
  keep <- t <= tau
  t <- t[keep]
  y <- y[keep]
  if (length(t) < 2L) {
    return(NA_real_)
  }
  if (max(t) < tau) {
    t <- c(t, tau)
    y <- c(y, y[length(y)])
  }
  as.numeric(sum(diff(t) * (y[-1] + y[-length(y)]) / 2) / tau)
}

summarise_oof_metrics <- function(
    oof_dt,
    nfeat_map,
    horizons_months,
    horizon_labels,
    iauc_tmin = 12,
    iauc_tmax = 60,
    iauc_step = 3,
    uno_tau = 60,
    ibs_tau = 60,
    ibs_step = 3) {
  keys <- unique(oof_dt[, .(METHOD, LEARNER, ENDPOINT)])
  auc_rows <- list()
  fold_auc_rows <- list()
  cindex_rows <- list()
  uno_rows <- list()
  iauc_rows <- list()
  brier_rows <- list()
  ibs_grid <- seq(ibs_step, ibs_tau, by = ibs_step)
  brier_times <- sort(unique(c(horizons_months, ibs_grid)))

  for (i in seq_len(nrow(keys))) {
    key <- keys[i]
    dt <- oof_dt[
      METHOD == key$METHOD & LEARNER == key$LEARNER & ENDPOINT == key$ENDPOINT
    ]
    dt <- dt[is.finite(TIME) & is.finite(EVENT) & is.finite(RISK_SCORE)]
    nfeat <- nfeat_map[METHOD == key$METHOD, N_FEATURES]
    if (!length(nfeat) || is.na(nfeat[1])) {
      nfeat <- NA_integer_
    } else {
      nfeat <- as.integer(nfeat[1])
    }
    n_pat <- nrow(dt)
    n_evt <- sum(dt$EVENT)
    if (n_pat < 30L || n_evt < 5L) {
      next
    }

    for (h in horizons_months) {
      auc_rows[[length(auc_rows) + 1L]] <- data.table(
        METHOD = key$METHOD,
        LEARNER = key$LEARNER,
        ENDPOINT = key$ENDPOINT,
        HORIZON_MONTHS = h,
        HORIZON = horizon_labels[[as.character(h)]],
        AUC = td_auc_ipcw(dt$TIME, dt$EVENT, dt$RISK_SCORE, h),
        N_PATIENTS = n_pat,
        N_EVENTS = n_evt,
        N_FEATURES = nfeat
      )
    }
    for (fold in sort(unique(dt$FOLD))) {
      test <- dt[FOLD == fold]
      if (nrow(test) < 5L) {
        next
      }
      for (h in horizons_months) {
        fold_auc_rows[[length(fold_auc_rows) + 1L]] <- data.table(
          METHOD = key$METHOD,
          LEARNER = key$LEARNER,
          ENDPOINT = key$ENDPOINT,
          FOLD = fold,
          HORIZON_MONTHS = h,
          HORIZON = horizon_labels[[as.character(h)]],
          AUC = td_auc_ipcw(test$TIME, test$EVENT, test$RISK_SCORE, h),
          N_TEST = nrow(test),
          N_EVENTS_TEST = sum(test$EVENT)
        )
      }
    }

    hc <- harrell_cindex(dt$TIME, dt$EVENT, dt$RISK_SCORE)
    cindex_rows[[length(cindex_rows) + 1L]] <- data.table(
      METHOD = key$METHOD,
      LEARNER = key$LEARNER,
      ENDPOINT = key$ENDPOINT,
      C_INDEX = hc$c,
      C_SE = hc$se,
      N_PATIENTS = n_pat,
      N_EVENTS = n_evt,
      N_FEATURES = nfeat
    )
    uno_rows[[length(uno_rows) + 1L]] <- data.table(
      METHOD = key$METHOD,
      LEARNER = key$LEARNER,
      ENDPOINT = key$ENDPOINT,
      UNO_C = uno_cindex(dt$TIME, dt$EVENT, dt$RISK_SCORE, tau = uno_tau),
      UNO_TAU_MONTHS = uno_tau,
      N_PATIENTS = n_pat,
      N_EVENTS = n_evt,
      N_FEATURES = nfeat
    )
    ia <- integrated_auc(
      dt$TIME, dt$EVENT, dt$RISK_SCORE,
      t_min = iauc_tmin, t_max = iauc_tmax, step = iauc_step
    )
    iauc_rows[[length(iauc_rows) + 1L]] <- data.table(
      METHOD = key$METHOD,
      LEARNER = key$LEARNER,
      ENDPOINT = key$ENDPOINT,
      IAUC_KM = ia$iauc_km,
      IAUC_UNIFORM = ia$iauc_uniform,
      T_MIN = iauc_tmin,
      T_MAX = iauc_tmax,
      N_PATIENTS = n_pat,
      N_EVENTS = n_evt,
      N_FEATURES = nfeat
    )

    s_hat <- oof_surv_matrix(dt, brier_times)
    s_null <- oof_km_null_matrix(dt, brier_times)
    bs_curve <- vapply(seq_along(brier_times), function(j) {
      brier_ipcw_from_surv(dt$TIME, dt$EVENT, s_hat[, j], brier_times[j])
    }, numeric(1))
    bs_null_curve <- vapply(seq_along(brier_times), function(j) {
      brier_ipcw_from_surv(dt$TIME, dt$EVENT, s_null[, j], brier_times[j])
    }, numeric(1))
    for (h in horizons_months) {
      j <- match(h, brier_times)
      bs <- bs_curve[j]
      bs0 <- bs_null_curve[j]
      brier_rows[[length(brier_rows) + 1L]] <- data.table(
        METHOD = key$METHOD,
        LEARNER = key$LEARNER,
        ENDPOINT = key$ENDPOINT,
        HORIZON_MONTHS = h,
        HORIZON = horizon_labels[[as.character(h)]],
        BRIER = bs,
        BRIER_NULL = bs0,
        IPA = if (is.finite(bs) && is.finite(bs0) && bs0 > 1e-8) 1 - bs / bs0 else NA_real_,
        N_PATIENTS = n_pat,
        N_EVENTS = n_evt,
        N_FEATURES = nfeat
      )
    }
    ibs <- ibs_from_brier_curve(brier_times, bs_curve, ibs_tau)
    ibs0 <- ibs_from_brier_curve(brier_times, bs_null_curve, ibs_tau)
    brier_rows[[length(brier_rows) + 1L]] <- data.table(
      METHOD = key$METHOD,
      LEARNER = key$LEARNER,
      ENDPOINT = key$ENDPOINT,
      HORIZON_MONTHS = ibs_tau,
      HORIZON = "IBS",
      BRIER = ibs,
      BRIER_NULL = ibs0,
      IPA = if (is.finite(ibs) && is.finite(ibs0) && ibs0 > 1e-8) 1 - ibs / ibs0 else NA_real_,
      N_PATIENTS = n_pat,
      N_EVENTS = n_evt,
      N_FEATURES = nfeat
    )
    message(
      key$METHOD, " ", key$LEARNER, " ", key$ENDPOINT,
      ": UnoC=", sprintf("%.3f", uno_rows[[length(uno_rows)]]$UNO_C),
      " iAUC=", sprintf("%.3f", ia$iauc_km),
      " IBS=", sprintf("%.3f", ibs)
    )
  }

  list(
    auc = rbindlist(auc_rows, use.names = TRUE, fill = TRUE),
    fold_auc = rbindlist(fold_auc_rows, use.names = TRUE, fill = TRUE),
    cindex = rbindlist(cindex_rows, use.names = TRUE, fill = TRUE),
    uno = rbindlist(uno_rows, use.names = TRUE, fill = TRUE),
    iauc = rbindlist(iauc_rows, use.names = TRUE, fill = TRUE),
    brier = rbindlist(brier_rows, use.names = TRUE, fill = TRUE)
  )
}

fit_cox_glmnet <- function(x, time, event, alpha = 0) {
  y <- Surv(time, event)
  # Prefer cv.glmnet; fall back to a small fixed lambda path if CV fails.
  fit <- tryCatch(
    {
      cv <- cv.glmnet(
        x = x,
        y = y,
        family = "cox",
        alpha = alpha,
        nfolds = min(5L, max(3L, nrow(x) %/% 10L)),
        standardize = FALSE,
        maxit = 1e5
      )
      list(model = cv$glmnet.fit, lambda = cv$lambda.min, type = "cv")
    },
    error = function(e) NULL
  )
  if (is.null(fit)) {
    path <- tryCatch(
      glmnet(
        x = x, y = y, family = "cox", alpha = alpha,
        standardize = FALSE, maxit = 1e5
      ),
      error = function(e) NULL
    )
    if (is.null(path)) {
      return(NULL)
    }
    lam <- path$lambda[max(1L, length(path$lambda) %/% 2L)]
    fit <- list(model = path, lambda = lam, type = "fixed")
  }
  fit
}

predict_cox_glmnet <- function(fit, x) {
  as.numeric(predict(fit$model, newx = x, s = fit$lambda, type = "link"))
}

fit_cox <- function(x, time, event) {
  df <- as.data.frame(x)
  colnames(df) <- make.names(colnames(df), unique = TRUE)
  df$TIME <- as.numeric(time)
  df$EVENT <- as.integer(event)
  # Drop near-constant columns that break coxph
  feat_names <- setdiff(colnames(df), c("TIME", "EVENT"))
  keep <- vapply(feat_names, function(nm) {
    v <- df[[nm]]
    is.finite(sd(v, na.rm = TRUE)) && sd(v, na.rm = TRUE) > 1e-8
  }, logical(1))
  feat_names <- feat_names[keep]
  if (!length(feat_names)) {
    return(NULL)
  }
  form <- as.formula(
    paste("Surv(TIME, EVENT) ~", paste(feat_names, collapse = " + "))
  )
  suppressWarnings(
    tryCatch(
      coxph(form, data = df, ties = "efron", model = FALSE, x = FALSE, y = FALSE),
      error = function(e) NULL
    )
  )
}

predict_cox <- function(fit, x) {
  if (is.null(fit)) {
    return(rep(NA_real_, nrow(x)))
  }
  df <- as.data.frame(x)
  colnames(df) <- make.names(colnames(df), unique = TRUE)
  needed <- attr(terms(fit), "term.labels")
  for (nm in setdiff(needed, colnames(df))) {
    df[[nm]] <- 0
  }
  tryCatch(
    as.numeric(predict(fit, newdata = df, type = "lp")),
    error = function(e) rep(NA_real_, nrow(df))
  )
}

fit_random_forest <- function(x, time, event, seed = 123L) {
  df <- as.data.frame(x)
  # ranger requires syntactically valid, unique names
  colnames(df) <- make.names(colnames(df), unique = TRUE)
  df$TIME <- as.numeric(time)
  df$EVENT <- as.integer(event)
  tryCatch(
    ranger(
      Surv(TIME, EVENT) ~ .,
      data = df,
      num.trees = rf_num_trees,
      min.node.size = 15L,
      importance = "none",
      seed = seed,
      verbose = FALSE
    ),
    error = function(e) NULL
  )
}

predict_random_forest <- function(fit, x) {
  if (is.null(fit)) {
    return(rep(NA_real_, nrow(x)))
  }
  df <- as.data.frame(x)
  colnames(df) <- make.names(colnames(df), unique = TRUE)
  pred <- tryCatch(
    predict(fit, data = df),
    error = function(e) NULL
  )
  if (is.null(pred) || is.null(pred$chf)) {
    return(rep(NA_real_, nrow(x)))
  }
  # Mean cumulative hazard: higher = higher risk (same direction as Cox LP)
  as.numeric(rowMeans(pred$chf))
}

predict_deepsurv_oof <- function(dt, features, seed = 123L) {
  n <- nrow(dt)
  if (!file.exists(deepsurv_python)) {
    warning("DeepSurv python not found: ", deepsurv_python)
    return(rep(NA_real_, n))
  }
  if (!file.exists(deepsurv_script)) {
    warning("DeepSurv script not found: ", deepsurv_script)
    return(rep(NA_real_, n))
  }
  in_file <- tempfile(pattern = "deepsurv_in_", fileext = ".tsv")
  out_file <- tempfile(pattern = "deepsurv_out_", fileext = ".tsv")
  log_file <- tempfile(pattern = "deepsurv_log_", fileext = ".txt")
  on.exit(unlink(c(in_file, out_file, log_file)), add = TRUE)
  keep <- c("PATIENT_ID", "FOLD", "TIME", "EVENT", features)
  fwrite(dt[, ..keep], in_file, sep = "\t")
  status <- tryCatch(
    system2(
      deepsurv_python,
      c(
        deepsurv_script,
        "--input", in_file,
        "--output", out_file,
        "--seed", as.character(seed),
        "--device", "cpu"
      ),
      stdout = log_file,
      stderr = log_file
    ),
    error = function(e) {
      writeLines(conditionMessage(e), log_file)
      1L
    }
  )
  if (!identical(status, 0L) || !file.exists(out_file) || file.info(out_file)$size < 10) {
    msg <- if (file.exists(log_file)) paste(readLines(log_file, warn = FALSE), collapse = "\n") else ""
    warning("DeepSurv failed (status=", status, "): ", msg)
    return(rep(NA_real_, n))
  }
  pred <- fread(out_file)
  as.numeric(pred$RISK_SCORE[match(dt$PATIENT_ID, pred$PATIENT_ID)])
}

make_folds <- function(n, event, k = 10L) {
  # Stratify approximately by event status.
  idx <- seq_len(n)
  folds <- integer(n)
  for (ev in c(0L, 1L)) {
    ids <- idx[event == ev]
    if (!length(ids)) {
      next
    }
    ids <- sample(ids)
    fold_ids <- rep_len(seq_len(k), length(ids))
    folds[ids] <- fold_ids
  }
  # Assign any leftovers
  if (any(folds == 0L)) {
    folds[folds == 0L] <- sample(seq_len(k), sum(folds == 0L), replace = TRUE)
  }
  folds
}

# ------------------------------- Load data -----------------------------------

estimates <- fread(estimate_file)
clinical_sample <- read_cbioportal_clinical(clinical_sample_file)
clinical_patient <- read_cbioportal_clinical(clinical_patient_file)

need_patient <- c(
  "PATIENT_ID", "OS_STATUS", "OS_MONTHS", "FFS_STATUS", "FFS_MONTHS"
)
meta <- merge(
  clinical_sample[, .(SAMPLE_ID, PATIENT_ID)],
  clinical_patient[, ..need_patient],
  by = "PATIENT_ID",
  all.x = TRUE
)

estimates <- estimates[METHOD %in% method_order]
if ("xcell" %in% estimates$METHOD && length(xcell_drop)) {
  estimates <- estimates[
    !(METHOD == "xcell" & CELL_TYPE %in% xcell_drop)
  ]
}

# Wide sample-level matrices per method, then average to patient.
method_patient <- list()
for (method in intersect(method_order, unique(estimates$METHOD))) {
  long <- estimates[METHOD == method]
  wide <- dcast(
    long,
    SAMPLE_ID ~ CELL_TYPE,
    value.var = "ESTIMATE",
    fun.aggregate = mean
  )
  wide <- merge(wide, meta, by = "SAMPLE_ID", all.x = TRUE)
  feature_cols <- setdiff(
    names(wide),
    c("SAMPLE_ID", "PATIENT_ID", need_patient)
  )
  # Patient-level mean of features; survival from first clinical row
  pat <- wide[
    !is.na(PATIENT_ID),
    c(
      .(
        N_SAMPLES = .N,
        OS_STATUS = OS_STATUS[1L],
        OS_MONTHS = OS_MONTHS[1L],
        FFS_STATUS = FFS_STATUS[1L],
        FFS_MONTHS = FFS_MONTHS[1L]
      ),
      lapply(.SD, function(z) mean(as.numeric(z), na.rm = TRUE))
    ),
    by = PATIENT_ID,
    .SDcols = feature_cols
  ]
  # Drop all-NA / zero-variance features
  keep <- feature_cols[
    vapply(feature_cols, function(f) {
      v <- pat[[f]]
      is.finite(sd(v, na.rm = TRUE)) && sd(v, na.rm = TRUE) > 0
    }, logical(1))
  ]
  pat <- pat[, c("PATIENT_ID", "N_SAMPLES", need_patient[-1], keep), with = FALSE]
  method_patient[[method]] <- list(data = pat, features = keep)
  message(method, ": ", nrow(pat), " patients, ", length(keep), " features")
}

# --------------------------- Fixed folds per endpoint ------------------------

# Define the analysis cohort once from the intersection of patients present in
# every method (with complete survival), then assign folds once per endpoint
# and reuse the same PATIENT_ID → fold map across methods.
common_patients <- Reduce(
  intersect,
  lapply(method_patient, function(x) x$data$PATIENT_ID)
)
message("Common patients across methods: ", length(common_patients))

# Survival labels are identical across methods; take from the first method.
surv_ref <- method_patient[[1]]$data[
  PATIENT_ID %in% common_patients,
  .(
    PATIENT_ID,
    OS_STATUS, OS_MONTHS, FFS_STATUS, FFS_MONTHS
  )
]
setorder(surv_ref, PATIENT_ID)

fold_maps <- list()
for (endpoint in c("OS", "FFS")) {
  time_col <- paste0(endpoint, "_MONTHS")
  status_col <- paste0(endpoint, "_STATUS")
  surv_ep <- surv_ref[
    ,
    .(
      PATIENT_ID,
      TIME = suppressWarnings(as.numeric(get(time_col))),
      EVENT = event_from_status(get(status_col), endpoint)
    )
  ]
  surv_ep <- surv_ep[is.finite(TIME) & TIME > 0 & !is.na(EVENT)]
  setorder(surv_ep, PATIENT_ID)
  if (nrow(surv_ep) < 50L || sum(surv_ep$EVENT) < 10L) {
    stop("Insufficient survival data for endpoint ", endpoint)
  }
  set.seed(123L + ifelse(endpoint == "OS", 0L, 1L))
  folds <- make_folds(nrow(surv_ep), surv_ep$EVENT, k = n_folds)
  fold_maps[[endpoint]] <- data.table(
    PATIENT_ID = surv_ep$PATIENT_ID,
    TIME = surv_ep$TIME,
    EVENT = surv_ep$EVENT,
    FOLD = folds
  )
  message(
    endpoint, " fixed folds: n=", nrow(fold_maps[[endpoint]]),
    " events=", sum(fold_maps[[endpoint]]$EVENT)
  )
}
fwrite(
  rbindlist(lapply(names(fold_maps), function(ep) {
    copy(fold_maps[[ep]])[, ENDPOINT := ep]
  }), use.names = TRUE),
  file.path(output_dir, "fixed_folds_by_endpoint.tsv"),
  sep = "\t"
)

# --------------------------- CV learners + time AUC --------------------------

metrics_from_oof <- identical(Sys.getenv("METRICS_FROM_OOF"), "1")
if (metrics_from_oof) {
  oof_file <- file.path(output_dir, "oof_risk_scores.tsv")
  if (!file.exists(oof_file)) {
    stop("METRICS_FROM_OOF=1 but missing: ", oof_file)
  }
  message("Recomputing metrics from existing OOF scores: ", oof_file)
  oof_dt <- fread(oof_file)
} else {

oof_rows <- list()
auc_rows <- list()
cindex_rows <- list()
fold_auc_rows <- list()

for (endpoint in c("OS", "FFS")) {
  fold_map <- fold_maps[[endpoint]]

  for (method in names(method_patient)) {
    pat <- copy(method_patient[[method]]$data)
    features <- method_patient[[method]]$features

    dt <- merge(
      fold_map,
      pat[, c("PATIENT_ID", features), with = FALSE],
      by = "PATIENT_ID",
      all = FALSE
    )
    setorder(dt, PATIENT_ID)
    feat_ok <- complete.cases(dt[, ..features])
    dt <- dt[feat_ok]
    if (nrow(dt) < 50L || sum(dt$EVENT) < 10L) {
      message("Skipping ", method, " ", endpoint, ": insufficient events")
      next
    }
    if (!identical(sort(dt$PATIENT_ID), sort(fold_map$PATIENT_ID))) {
      warning(
        method, " ", endpoint,
        ": patient set differs from fixed fold map after complete-case filter; ",
        "using intersection only."
      )
    }

    for (learner in learner_order) {
      risk_oof <- rep(NA_real_, nrow(dt))

      if (learner == "deepsurv") {
        message(method, " deepsurv ", endpoint, ": training Cox-PH MLP on 10 folds")
        risk_oof <- predict_deepsurv_oof(
          dt,
          features,
          seed = 123L + ifelse(endpoint == "OS", 0L, 100L)
        )
      } else {
        for (fold in seq_len(n_folds)) {
          test_idx <- which(dt$FOLD == fold)
          train_idx <- which(dt$FOLD != fold)
          if (length(test_idx) < 5L || length(train_idx) < 30L) {
            next
          }
          if (sum(dt$EVENT[train_idx]) < 5L) {
            next
          }

          x_train <- as.matrix(dt[train_idx, ..features])
          x_test <- as.matrix(dt[test_idx, ..features])
          storage.mode(x_train) <- "numeric"
          storage.mode(x_test) <- "numeric"

          if (learner == "cox") {
            zs <- zscore_matrix(x_train, x_test)
            fit <- fit_cox(
              zs$train,
              dt$TIME[train_idx],
              dt$EVENT[train_idx]
            )
            if (is.null(fit)) {
              next
            }
            risk_oof[test_idx] <- predict_cox(fit, zs$test)
          } else if (learner %in% c("cox_lasso", "cox_ridge")) {
            zs <- zscore_matrix(x_train, x_test)
            alpha <- if (learner == "cox_lasso") 1 else 0
            fit <- fit_cox_glmnet(
              zs$train,
              dt$TIME[train_idx],
              dt$EVENT[train_idx],
              alpha = alpha
            )
            if (is.null(fit)) {
              next
            }
            risk_oof[test_idx] <- predict_cox_glmnet(fit, zs$test)
          } else if (learner == "random_forest") {
            fit <- fit_random_forest(
              x_train,
              dt$TIME[train_idx],
              dt$EVENT[train_idx],
              seed = 123L + fold
            )
            if (is.null(fit)) {
              next
            }
            risk_oof[test_idx] <- predict_random_forest(fit, x_test)
          } else {
            stop("Unknown learner: ", learner)
          }
        }
      }

      for (fold in seq_len(n_folds)) {
        test_idx <- which(dt$FOLD == fold)
        if (!length(test_idx)) {
          next
        }
        for (h in horizons_months) {
          auc_f <- td_auc_ipcw(
            dt$TIME[test_idx],
            dt$EVENT[test_idx],
            risk_oof[test_idx],
            h
          )
          fold_auc_rows[[length(fold_auc_rows) + 1L]] <- data.table(
            METHOD = method,
            LEARNER = learner,
            ENDPOINT = endpoint,
            FOLD = fold,
            HORIZON_MONTHS = h,
            HORIZON = horizon_labels[[as.character(h)]],
            AUC = auc_f,
            N_TEST = length(test_idx),
            N_EVENTS_TEST = sum(dt$EVENT[test_idx])
          )
        }
      }

      ok_risk <- is.finite(risk_oof)
      if (sum(ok_risk) < 30L) {
        message(
          "Skipping AUCs for ", method, " ", learner, " ", endpoint,
          ": too few OOF scores"
        )
        next
      }

      oof_rows[[length(oof_rows) + 1L]] <- data.table(
        METHOD = method,
        LEARNER = learner,
        ENDPOINT = endpoint,
        PATIENT_ID = dt$PATIENT_ID,
        FOLD = dt$FOLD,
        TIME = dt$TIME,
        EVENT = dt$EVENT,
        RISK_SCORE = risk_oof
      )

      for (h in horizons_months) {
        auc <- td_auc_ipcw(dt$TIME, dt$EVENT, risk_oof, h)
        auc_rows[[length(auc_rows) + 1L]] <- data.table(
          METHOD = method,
          LEARNER = learner,
          ENDPOINT = endpoint,
          HORIZON_MONTHS = h,
          HORIZON = horizon_labels[[as.character(h)]],
          AUC = auc,
          N_PATIENTS = sum(ok_risk),
          N_EVENTS = sum(dt$EVENT[ok_risk]),
          N_FEATURES = length(features)
        )
      }

      cfit <- tryCatch(
        concordance(Surv(dt$TIME, dt$EVENT) ~ risk_oof),
        error = function(e) NULL
      )
      cindex_rows[[length(cindex_rows) + 1L]] <- data.table(
        METHOD = method,
        LEARNER = learner,
        ENDPOINT = endpoint,
        C_INDEX = if (!is.null(cfit)) as.numeric(cfit$concordance) else NA_real_,
        C_SE = if (!is.null(cfit)) as.numeric(sqrt(cfit$var)) else NA_real_,
        N_PATIENTS = sum(ok_risk),
        N_EVENTS = sum(dt$EVENT[ok_risk]),
        N_FEATURES = length(features)
      )

      message(
        method, " ", learner, " ", endpoint, ": n=", sum(ok_risk),
        " events=", sum(dt$EVENT[ok_risk]),
        " AUC1/3/5=",
        paste(sprintf("%.3f", sapply(horizons_months, function(h) {
          td_auc_ipcw(dt$TIME, dt$EVENT, risk_oof, h)
        })), collapse = "/")
      )
    }
  }
}

}

if (!metrics_from_oof) {
  oof_dt <- rbindlist(oof_rows, use.names = TRUE, fill = TRUE)
  if (Sys.getenv("MERGE_PREVIOUS", unset = "1") == "1") {
    merge_keep <- function(new_dt, path) {
      if (!file.exists(path) || !nrow(new_dt)) {
        return(new_dt)
      }
      old <- tryCatch(fread(path), error = function(e) NULL)
      if (is.null(old) || !nrow(old) || !"LEARNER" %in% names(old)) {
        return(new_dt)
      }
      ran <- unique(as.character(new_dt$LEARNER))
      old <- old[!as.character(LEARNER) %in% ran]
      if (!nrow(old)) {
        return(new_dt)
      }
      rbindlist(list(old, new_dt), use.names = TRUE, fill = TRUE)
    }
    oof_dt <- merge_keep(oof_dt, file.path(output_dir, "oof_risk_scores.tsv"))
  }
}

if (!nrow(oof_dt)) {
  stop("No OOF risk scores were produced")
}

nfeat_map <- rbindlist(lapply(names(method_patient), function(m) {
  data.table(METHOD = m, N_FEATURES = length(method_patient[[m]]$features))
}))
message("Summarising pooled OOF metrics (AUC, Harrell C, Uno C, iAUC, Brier) ...")
summ <- summarise_oof_metrics(
  oof_dt,
  nfeat_map,
  horizons_months,
  horizon_labels
)
auc_dt <- summ$auc
fold_auc_dt <- summ$fold_auc
cindex_dt <- summ$cindex
uno_dt <- summ$uno
iauc_dt <- summ$iauc
brier_dt <- summ$brier

if (!nrow(auc_dt)) {
  stop("No AUC results were produced")
}

# Facet order always shows all known learners present in the table.
plot_learner_order <- all_learners
set_plot_factors <- function(dt) {
  dt[, METHOD := factor(METHOD, levels = intersect(method_order, unique(METHOD)))]
  dt[, LEARNER := factor(LEARNER, levels = intersect(plot_learner_order, unique(LEARNER)))]
  dt[, ENDPOINT := factor(ENDPOINT, levels = c("OS", "FFS"))]
  dt
}
auc_dt <- set_plot_factors(auc_dt)
auc_dt[, HORIZON := factor(HORIZON, levels = horizon_labels)]
fold_auc_dt <- set_plot_factors(fold_auc_dt)
fold_auc_dt[, HORIZON := factor(HORIZON, levels = horizon_labels)]
cindex_dt <- set_plot_factors(cindex_dt)
uno_dt <- set_plot_factors(uno_dt)
iauc_dt <- set_plot_factors(iauc_dt)
brier_dt <- set_plot_factors(brier_dt)
brier_dt[, HORIZON := factor(HORIZON, levels = c(horizon_labels, "IBS"))]

# Fold-level mean ± SE
fold_summary <- fold_auc_dt[
  ,
  .(
    AUC_MEAN = mean(AUC, na.rm = TRUE),
    AUC_SD = sd(AUC, na.rm = TRUE),
    AUC_SE = sd(AUC, na.rm = TRUE) / sqrt(sum(is.finite(AUC))),
    N_FOLDS = sum(is.finite(AUC))
  ),
  by = .(METHOD, LEARNER, ENDPOINT, HORIZON_MONTHS, HORIZON)
]
fold_summary[, METHOD := factor(METHOD, levels = levels(auc_dt$METHOD))]
fold_summary[, LEARNER := factor(LEARNER, levels = levels(auc_dt$LEARNER))]
fold_summary[, ENDPOINT := factor(ENDPOINT, levels = c("OS", "FFS"))]
fold_summary[, HORIZON := factor(HORIZON, levels = horizon_labels)]

fwrite(auc_dt, file.path(output_dir, "td_auc_pooled_oof.tsv"), sep = "\t")
fwrite(fold_auc_dt, file.path(output_dir, "td_auc_per_fold.tsv"), sep = "\t")
fwrite(fold_summary, file.path(output_dir, "td_auc_fold_summary.tsv"), sep = "\t")
fwrite(cindex_dt, file.path(output_dir, "cindex_oof.tsv"), sep = "\t")
fwrite(uno_dt, file.path(output_dir, "uno_cindex_oof.tsv"), sep = "\t")
fwrite(iauc_dt, file.path(output_dir, "iauc_oof.tsv"), sep = "\t")
fwrite(brier_dt, file.path(output_dir, "brier_oof.tsv"), sep = "\t")
fwrite(oof_dt, file.path(output_dir, "oof_risk_scores.tsv"), sep = "\t")

auc_wide <- dcast(
  auc_dt,
  METHOD + LEARNER + ENDPOINT + N_PATIENTS + N_EVENTS + N_FEATURES ~ HORIZON,
  value.var = "AUC"
)
fwrite(auc_wide, file.path(output_dir, "td_auc_wide.tsv"), sep = "\t")

brier_wide <- dcast(
  brier_dt,
  METHOD + LEARNER + ENDPOINT + N_PATIENTS + N_EVENTS + N_FEATURES ~ HORIZON,
  value.var = "BRIER"
)
fwrite(brier_wide, file.path(output_dir, "brier_wide.tsv"), sep = "\t")

extra_wide <- merge(
  merge(
    cindex_dt[, .(METHOD, LEARNER, ENDPOINT, HARRELL_C = C_INDEX, HARRELL_C_SE = C_SE)],
    uno_dt[, .(METHOD, LEARNER, ENDPOINT, UNO_C, UNO_TAU_MONTHS)],
    by = c("METHOD", "LEARNER", "ENDPOINT")
  ),
  iauc_dt[, .(METHOD, LEARNER, ENDPOINT, IAUC_KM, IAUC_UNIFORM, T_MIN, T_MAX)],
  by = c("METHOD", "LEARNER", "ENDPOINT")
)
ibs_wide <- dcast(
  brier_dt,
  METHOD + LEARNER + ENDPOINT ~ HORIZON,
  value.var = c("BRIER", "IPA")
)
extra_wide <- merge(extra_wide, ibs_wide, by = c("METHOD", "LEARNER", "ENDPOINT"))
fwrite(extra_wide, file.path(output_dir, "extra_metrics_wide.tsv"), sep = "\t")

# --------------------------------- Plots -------------------------------------

method_colors <- c(
  bayesprism = "#4C78A8",
  cibersort = "#E45756",
  quantiseq = "#72B7B2",
  mcp_counter = "#F58518",
  xcell = "#54A24B"
)
horizon_colors <- c(
  "1-year" = "#4C78A8",
  "3-year" = "#F58518",
  "5-year" = "#E45756"
)
theme_bench <- theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "bottom"
  )

p1 <- ggplot(auc_dt, aes(x = METHOD, y = AUC, fill = HORIZON)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "grey40") +
  facet_grid(LEARNER ~ ENDPOINT) +
  coord_cartesian(ylim = c(0.4, 0.85)) +
  scale_fill_manual(values = horizon_colors) +
  labs(
    title = "Deconvolution methods: time-dependent AUC (10-fold CV, pooled OOF)",
    subtitle = "Cox / Lasso Cox / Cox ridge / RSF / DeepSurv; IPCW cumulative/dynamic AUC",
    x = NULL,
    y = "AUC",
    fill = "Horizon"
  ) +
  theme_bench
ggsave(file.path(output_dir, "td_auc_pooled_by_method.pdf"), p1, width = 10, height = 14)

p2 <- ggplot(fold_summary, aes(x = METHOD, y = AUC_MEAN, fill = HORIZON)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = AUC_MEAN - AUC_SE, ymax = AUC_MEAN + AUC_SE),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "grey40") +
  facet_grid(LEARNER ~ ENDPOINT) +
  coord_cartesian(ylim = c(0.35, 0.9)) +
  scale_fill_manual(values = horizon_colors) +
  labs(
    title = "Mean per-fold time-dependent AUC (± SE)",
    x = NULL,
    y = "AUC",
    fill = "Horizon"
  ) +
  theme_bench
ggsave(file.path(output_dir, "td_auc_foldmean_by_method.pdf"), p2, width = 10, height = 14)

p3 <- ggplot(auc_dt, aes(x = HORIZON, y = AUC, color = METHOD, group = METHOD)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "grey40") +
  facet_grid(LEARNER ~ ENDPOINT) +
  coord_cartesian(ylim = c(0.4, 0.85)) +
  scale_color_manual(values = method_colors) +
  labs(
    title = "Time-dependent AUC across survival horizons",
    x = NULL,
    y = "AUC",
    color = "Method"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "bottom"
  )
ggsave(file.path(output_dir, "td_auc_lines_by_horizon.pdf"), p3, width = 9, height = 12)

p4 <- ggplot(cindex_dt, aes(x = METHOD, y = C_INDEX, fill = ENDPOINT)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_errorbar(
    aes(ymin = C_INDEX - C_SE, ymax = C_INDEX + C_SE),
    position = position_dodge(width = 0.75),
    width = 0.2
  ) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "grey40") +
  facet_wrap(~LEARNER, nrow = 1) +
  coord_cartesian(ylim = c(0.45, 0.8)) +
  scale_fill_manual(values = c(OS = "#4C78A8", FFS = "#E45756")) +
  labs(
    title = "Harrell C-index from 10-fold OOF risk scores",
    subtitle = "Higher score = higher risk (concordance on -risk)",
    x = NULL,
    y = "C-index",
    fill = "Endpoint"
  ) +
  theme_bench
ggsave(file.path(output_dir, "cindex_by_method.pdf"), p4, width = 14, height = 4.8)

p5 <- ggplot(uno_dt, aes(x = METHOD, y = UNO_C, fill = ENDPOINT)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "grey40") +
  facet_wrap(~LEARNER, nrow = 1) +
  coord_cartesian(ylim = c(0.45, 0.8)) +
  scale_fill_manual(values = c(OS = "#4C78A8", FFS = "#E45756")) +
  labs(
    title = "Uno's C-index (IPCW, tau = 60 months)",
    subtitle = "Higher score = higher risk; inverse-probability-of-censoring weights",
    x = NULL,
    y = "Uno C",
    fill = "Endpoint"
  ) +
  theme_bench
ggsave(file.path(output_dir, "uno_cindex_by_method.pdf"), p5, width = 14, height = 4.8)

p6 <- ggplot(iauc_dt, aes(x = METHOD, y = IAUC_KM, fill = ENDPOINT)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "grey40") +
  facet_wrap(~LEARNER, nrow = 1) +
  coord_cartesian(ylim = c(0.45, 0.8)) +
  scale_fill_manual(values = c(OS = "#4C78A8", FFS = "#E45756")) +
  labs(
    title = "Integrated AUC (KM-weighted, 12-60 months)",
    subtitle = "Trapezoid of cumulative/dynamic IPCW AUC, weighted by KM survival drop",
    x = NULL,
    y = "iAUC",
    fill = "Endpoint"
  ) +
  theme_bench
ggsave(file.path(output_dir, "iauc_by_method.pdf"), p6, width = 14, height = 4.8)

brier_plot <- brier_dt[HORIZON %in% horizon_labels]
p7 <- ggplot(brier_plot, aes(x = METHOD, y = BRIER, fill = HORIZON)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_grid(LEARNER ~ ENDPOINT) +
  scale_fill_manual(values = horizon_colors) +
  labs(
    title = "IPCW Brier score (lower is better)",
    subtitle = "Fold-wise Breslow survival from OOF risk; Graf inverse-probability weights",
    x = NULL,
    y = "Brier score",
    fill = "Horizon"
  ) +
  theme_bench
ggsave(file.path(output_dir, "brier_by_method.pdf"), p7, width = 10, height = 14)

ibs_plot <- brier_dt[HORIZON == "IBS"]
p8 <- ggplot(ibs_plot, aes(x = METHOD, y = BRIER, fill = ENDPOINT)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  facet_wrap(~LEARNER, nrow = 1) +
  scale_fill_manual(values = c(OS = "#4C78A8", FFS = "#E45756")) +
  labs(
    title = "Integrated Brier score (0-60 months; lower is better)",
    x = NULL,
    y = "IBS",
    fill = "Endpoint"
  ) +
  theme_bench
ggsave(file.path(output_dir, "ibs_by_method.pdf"), p8, width = 14, height = 4.8)

p9 <- ggplot(brier_plot, aes(x = METHOD, y = IPA, fill = HORIZON)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
  facet_grid(LEARNER ~ ENDPOINT) +
  scale_fill_manual(values = horizon_colors) +
  labs(
    title = "Index of prediction accuracy (1 - Brier / Brier_KM)",
    subtitle = "0 = KM null; higher is better. Negative = worse than KM.",
    x = NULL,
    y = "IPA",
    fill = "Horizon"
  ) +
  theme_bench
ggsave(file.path(output_dir, "ipa_by_method.pdf"), p9, width = 10, height = 14)

rank_dt <- copy(auc_dt)
rank_dt[, RANK := frank(-AUC, ties.method = "min"), by = .(LEARNER, ENDPOINT, HORIZON)]
fwrite(
  rank_dt[order(LEARNER, ENDPOINT, HORIZON, RANK)],
  file.path(output_dir, "td_auc_ranks.tsv"),
  sep = "\t"
)

message("Wrote benchmark results to ", normalizePath(output_dir))
print(auc_wide[order(LEARNER, ENDPOINT, METHOD)])
print(extra_wide[order(LEARNER, ENDPOINT, -UNO_C)])
