# Shared clinical covariates and adjusted association model for
# mutation_microenv_association.R and genomic_microenv_association.R.
#
# Primary association model:
#   fraction ~ alteration + IGHV + sex + age_at_sampling/10
#
# Returns the alteration coefficient (beta), SE, p-value, and complete-case n.

encode_ighv_mutated <- function(status) {
  status <- as.character(status)
  out <- rep(NA_real_, length(status))
  out[grepl("unmut", status, ignore.case = TRUE)] <- 0
  out[
    grepl("^mut|mutated", status, ignore.case = TRUE) &
      !grepl("unmut", status, ignore.case = TRUE)
  ] <- 1
  out
}

encode_sex_male <- function(sex) {
  sex <- toupper(trimws(as.character(sex)))
  fifelse(sex %in% "MALE", 1, fifelse(sex %in% "FEMALE", 0, NA_real_))
}

attach_clinical_covariates <- function(micro_dt, clinical_sample, clinical_patient) {
  meta <- merge(
    clinical_sample[, .(SAMPLE_ID, PATIENT_ID)],
    clinical_patient[, .(
      PATIENT_ID,
      SEX,
      AGE,
      AGE_SAMPLING,
      IGHV_MUTATION_STATUS
    )],
    by = "PATIENT_ID",
    all.x = TRUE
  )
  meta[, IGHV_MUTATED := encode_ighv_mutated(IGHV_MUTATION_STATUS)]
  meta[, SEX_MALE := encode_sex_male(SEX)]
  meta[, AGE_YEARS := suppressWarnings(as.numeric(AGE_SAMPLING))]
  meta[is.na(AGE_YEARS), AGE_YEARS := suppressWarnings(as.numeric(AGE))]
  merge(micro_dt, meta, by = "SAMPLE_ID", all.x = TRUE)
}

# Linear model controlling for IGHV, sex, and age (per decade).
fit_alt_covariate_adjusted <- function(
    y,
    altered,
    ighv_mutated,
    sex_male,
    age_years,
    min_n = 40L
) {
  ok <- is.finite(y) &
    !is.na(altered) &
    is.finite(ighv_mutated) &
    is.finite(sex_male) &
    is.finite(age_years)
  n <- sum(ok)
  empty <- list(
    beta = NA_real_, se = NA_real_, p = NA_real_, n = n,
    n_alt = NA_integer_, n_wt = NA_integer_
  )
  if (n < min_n) {
    return(empty)
  }
  alt <- as.numeric(altered[ok])
  if (length(unique(alt)) < 2L) {
    return(empty)
  }
  df <- data.frame(
    y = as.numeric(y[ok]),
    alt = alt,
    ighv = as.numeric(ighv_mutated[ok]),
    sex = as.numeric(sex_male[ok]),
    age10 = as.numeric(age_years[ok]) / 10
  )
  # Require both sexes and both IGHV classes when possible; still fit if one
  # covariate is nearly constant (lm drops collinear terms).
  fit <- tryCatch(
    lm(y ~ alt + ighv + sex + age10, data = df),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(empty)
  }
  s <- summary(fit)$coefficients
  if (!"alt" %in% rownames(s)) {
    return(empty)
  }
  list(
    beta = unname(s["alt", "Estimate"]),
    se = unname(s["alt", "Std. Error"]),
    p = unname(s["alt", "Pr(>|t|)"]),
    n = n,
    n_alt = as.integer(sum(alt == 1)),
    n_wt = as.integer(sum(alt == 0))
  )
}
