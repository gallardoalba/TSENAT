## Helper for mixed-model retries using different optimizers
.tsenat_try_lmer <- function(formula, data, REML = FALSE, verbose = FALSE, progress = FALSE) {
  mm_suppress_pattern <- "failed to converge|convergence|singular|boundary (.*) singular"
  control_list <- list(lmerControl(optimizer = "bobyqa"), lmerControl(optimizer = "Nelder_M"))
  # primary try
  fit <- try(lme4::lmer(formula, data = data, REML = REML), silent = TRUE)
  if (!inherits(fit, "try-error")) {
    return(fit)
  }
  # try alternate optimizers
  for (ctrl in control_list) {
    fit2 <- try(suppressWarnings(lme4::lmer(formula, data = data, REML = REML, control = ctrl)), silent = TRUE)
    if (!inherits(fit2, "try-error")) return(fit2)
  }
  return(fit)
}
## Helpers for calculate_lm_interaction

# Try lme4::lmer with multiple optimizers and controlled warnings.
.tsenat_try_lmer <- function(formula, data, suppress_lme4_warnings = TRUE, verbose = FALSE,
    mm_suppress_pattern = "boundary \\\\(singular\\\\) fit|Computed variance-covariance matrix problem|not a positive definite matrix") {
    if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("Package 'lme4' is required for mixed-model fitting")
    }
    opts <- list(list(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+05)),
        list(optimizer = "nloptwrap", optCtrl = list(maxfun = 5e+05)))
    for (o in opts) {
        ctrl <- lme4::lmerControl(optimizer = o$optimizer, optCtrl = o$optCtrl)
        muffle_cond <- suppress_lme4_warnings || (!verbose)
        fit_try <- withCallingHandlers(try(lme4::lmer(formula, data = data, REML = FALSE, control = ctrl),
            silent = TRUE), warning = function(w) {
            if (muffle_cond && grepl(mm_suppress_pattern, conditionMessage(w), ignore.case = TRUE))
                invokeRestart("muffleWarning")
        }, message = function(m) {
            if (muffle_cond && grepl(mm_suppress_pattern, conditionMessage(m), ignore.case = TRUE))
                invokeRestart("muffleMessage")
        })
        if (!inherits(fit_try, "try-error")) {
            # check singularity if function available
            is_sing <- FALSE
            if (exists("isSingular", where = asNamespace("lme4"), inherits = FALSE)) {
                is_sing <- tryCatch(lme4::isSingular(fit_try, tol = 1e-04), error = function(e) FALSE)
            }
            attr(fit_try, "singular") <- is_sing
            return(fit_try)
        }
    }
    # all attempts failed
    return(structure("error", class = "try-error"))
}

## Helpers for calculate_lm_interaction fallbacks, LRT and Satterthwaite p-values
.tsenat_try_lm_fallbacks <- function(df, verbose = FALSE) {
  # try lm with subject as fixed effect
  fit0_lm <- try(stats::lm(entropy ~ q + group + subject, data = df), silent = TRUE)
  fit1_lm <- try(stats::lm(entropy ~ q * group + subject, data = df), silent = TRUE)
  if (!inherits(fit0_lm, "try-error") && !inherits(fit1_lm, "try-error")) {
    return(list(fit0 = fit0_lm, fit1 = fit1_lm, method = "lm_subject"))
  }
  # last resort: drop subject, plain lm
  fit0_lm2 <- try(stats::lm(entropy ~ q + group, data = df), silent = TRUE)
  fit1_lm2 <- try(stats::lm(entropy ~ q * group, data = df), silent = TRUE)
  if (!inherits(fit0_lm2, "try-error") && !inherits(fit1_lm2, "try-error")) {
    return(list(fit0 = fit0_lm2, fit1 = fit1_lm2, method = "lm_nosubject"))
  }
  return(NULL)
}

.tsenat_extract_lrt_p <- function(fit0, fit1) {
  an <- try(stats::anova(fit0, fit1), silent = TRUE)
  if (!inherits(an, "try-error") && nrow(an) >= 2) {
    pcol <- grep("Pr\\(>F\\)|Pr\\(>Chisq\\)|Pr\\(>Chi\\)", colnames(an), value = TRUE)
    if (length(pcol) == 0) return(as.numeric(an[2, ncol(an)])) else return(as.numeric(an[2, pcol[1]]))
  }
  return(NA_real_)
}

.tsenat_extract_satterthwaite_p <- function(fit1, fallback_lm = NULL) {
  # if lmerTest available and fit1 is not singular use that
  if (!is.null(fallback_lm)) {
    coefs <- try(summary(fallback_lm$fit1)$coefficients, silent = TRUE)
    if (!inherits(coefs, "try-error")) {
      ia_idx <- grep("^q:group", rownames(coefs))
      if (length(ia_idx) > 0) return(coefs[ia_idx[1], "Pr(>|t|)"])
    }
    return(NA_real_)
  }
  if (requireNamespace("lmerTest", quietly = TRUE) && inherits(fit1, "lmerMod") && !isTRUE(attr(fit1, "singular"))) {
    fit_lt <- try(lmerTest::lmer(stats::formula(fit1), data = stats::model.frame(fit1), REML = FALSE), silent = TRUE)
    if (!inherits(fit_lt, "try-error")) {
      coefs <- summary(fit_lt)$coefficients
      ia_idx <- grep("^q:group", rownames(coefs))
      if (length(ia_idx) > 0) return(coefs[ia_idx[1], "Pr(>|t|)"])
    }
  }
  return(NA_real_)
}

.tsenat_prepare_fpca_matrix <- function(mat, min_frac = 0.01) {
  # prepare matrix for FPCA: center, scale, and drop near-constant rows
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  row_vars <- apply(mat, 1, stats::var, na.rm = TRUE)
  keep <- row_vars > (min_frac * max(row_vars, na.rm = TRUE))
  if (sum(keep) == 0) keep <- rep(TRUE, nrow(mat))
  m2 <- mat[keep, , drop = FALSE]
  m2 <- t(scale(t(m2)))
  return(list(mat = m2, keep = keep))
}
## Additional helpers for calculate_lm_interaction

# Try simple lm fallbacks when lmer fails. Returns list(fit0, fit1, method)
.tsenat_try_lm_fallbacks <- function(df, verbose = FALSE) {
    # try lm with subject as fixed effect
    fit0_lm <- try(stats::lm(entropy ~ q + group + subject, data = df), silent = TRUE)
    fit1_lm <- try(stats::lm(entropy ~ q * group + subject, data = df), silent = TRUE)
    if (!inherits(fit0_lm, "try-error") && !inherits(fit1_lm, "try-error")) {
        return(list(fit0 = fit0_lm, fit1 = fit1_lm, method = "lm_subject"))
    }
    # last resort: drop subject
    fit0_lm2 <- try(stats::lm(entropy ~ q + group, data = df), silent = TRUE)
    fit1_lm2 <- try(stats::lm(entropy ~ q * group, data = df), silent = TRUE)
    if (!inherits(fit0_lm2, "try-error") && !inherits(fit1_lm2, "try-error")) {
        return(list(fit0 = fit0_lm2, fit1 = fit1_lm2, method = "lm_nosubject"))
    }
    # indicate failure
    return(NULL)
}

# Extract LRT p-value from two fitted models (lm or lmer). Returns NA_real_ on failure.
.tsenat_extract_lrt_p <- function(fit0, fit1) {
    an <- try(stats::anova(fit0, fit1), silent = TRUE)
    if (inherits(an, "try-error") || nrow(an) < 2) return(NA_real_)
    pcol <- grep("Pr\\(>F\\)|Pr\\(>Chisq\\)|Pr\\(>Chi\\)", colnames(an), value = TRUE)
    if (length(pcol) == 0) {
        return(as.numeric(an[2, ncol(an)]))
    }
    return(as.numeric(an[2, pcol[1]]))
}

# Extract Satterthwaite p-value for q:group interaction from lmer/lmerTest or fallback lm.
.tsenat_extract_satterthwaite_p <- function(fit1, fallback_lm = NULL) {
    # prefer lmerTest when available
    if (!is.null(fallback_lm)) {
        coefs <- summary(fallback_lm$fit1)$coefficients
        ia_idx <- grep("^q:group", rownames(coefs))
        if (length(ia_idx) > 0) return(coefs[ia_idx[1], "Pr(>|t|)"])
        return(NA_real_)
    }
    if (requireNamespace("lmerTest", quietly = TRUE) && inherits(fit1, "lmerMod")) {
        fit_lt <- try(lmerTest::lmer(as.formula(formula(fit1)), data = model.frame(fit1), REML = FALSE), silent = TRUE)
        if (!inherits(fit_lt, "try-error")) {
            coefs <- summary(fit_lt)$coefficients
            ia_idx <- grep("^q:group", rownames(coefs))
            if (length(ia_idx) > 0) return(coefs[ia_idx[1], "Pr(>|t|)"])
        }
    }
    return(NA_real_)
}

# Helper for FPCA-style preprocessing used in calculate_lm_interaction fpca method.
# Builds curve_mat, filters good rows, imputes column means, and returns list(mat_sub, used_samples)
.tsenat_prepare_fpca_matrix <- function(mat, sample_names, q_vals, min_obs = 10) {
    uq <- sort(unique(q_vals))
    samples_u <- unique(sample_names)
    curve_mat <- matrix(NA_real_, nrow = length(samples_u), ncol = length(uq))
    if (length(samples_u) > 0) rownames(curve_mat) <- samples_u
    for (i in seq_along(sample_names)) {
        s <- sample_names[i]
        qv <- q_vals[i]
        qi <- match(qv, uq)
        if (is.na(qi)) next
        if (s %in% rownames(curve_mat)) curve_mat[s, qi] <- as.numeric(mat[, i])
    }
    # keep samples with at least half of q points present
    good_rows <- which(rowSums(!is.na(curve_mat)) >= max(2, ceiling(ncol(curve_mat)/2)))
    if (length(good_rows) < min_obs) return(NULL)
    mat_sub <- curve_mat[good_rows, , drop = FALSE]
    col_means <- apply(mat_sub, 2, function(col) mean(col, na.rm = TRUE))
    for (r in seq_len(nrow(mat_sub))) mat_sub[r, is.na(mat_sub[r, ])] <- col_means[is.na(mat_sub[r, ])]
    list(mat_sub = mat_sub, used_samples = rownames(mat_sub))
}
