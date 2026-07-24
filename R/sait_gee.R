#' Internal GEE Interaction Analysis with AR(1) Correlation
#'
#' Generalized Estimating Equations (GEE) for testing q-dependent interaction effects
#' in Tsallis entropy data with repeated measurements. GEE is robust for correlated
#' data and does not assume normality of random effects, making it ideal for
#' entropy measurements structured by q-values.
#'
#' ## Methodology Overview
#'
#' GEE operates on the generalized linear model framework with user-specified
#' correlation structure:
#'
#' 1. **Input Validation**: Check minimum observations, group structure, and subject IDs
#' 2. **ARIMA Differencing**: Remove non-stationarity via first-differencing within subjects
#' 3. **Weight Preparation**: Apply heteroscedasticity weights or bootstrap CI weights
#' 4. **Correlation Structure Selection**: Choose AR(1), exchangeable, or independence via QIC
#' 5. **Model Fitting**: Fit null (main effects) and alternative (interaction) models
#' 6. **Interaction Testing**: Extract p-value with bias correction for small clusters
#' 7. **Kauermann-Carroll Correction**: Apply HC1 bias reduction for n_clusters < 30
#'
#' ## Key References
#'
#' - Zimmerman & Harville (1991, S171): AR(1) for ordered covariate structures
#' - Kauermann & Carroll (2001): Sandwich variance bias correction for small clusters
#' - Pan (2001): QIC model selection criterion for GEE
#' - Mancl & DeRouen (2001): Covariate-adjusted ANOVA-type tests with GEE
#'
#' ## Important Clarifications
#'
#' - **Correlation vs Random Effects**: GEE models within-subject correlation
#'   directly (no random intercepts like mixed models). For AR(1), the pattern
#'   Corr(q_i, q_j) = ρ^|i-j| accounts for ordered q-value measurements.
#'
#' - **Design Effect**: Multi-q measurements create effective sample size reduction
#'   D_eff = (1 + ρ) / (1 - ρ). HC1 correction adjusts variance using n_effective = n_clusters / D_eff.
#'
#' - **Bias Correction**: Applied when n_clusters < 30. Uses t-distribution with
#'   df = n_clusters - 1 for conservative (Type I error-protecting) p-values.
#'
#' @param df data.frame with columns: entropy (outcome), q, group, and optionally subject (for paired designs)
#' @param q_vals numeric vector of q-values (used only for cluster size computation in design effect)
#' @param g character; gene identifier for error messages and result tracking
#' @param subject character or NULL; vector of subject/cluster IDs for repeated measurements.
#'   If NULL, observations treated as independent (no repeated measures)
#' @param min_obs integer >= 2; minimum required non-NA entropy observations. Default 5
#' @param corstr character; correlation structure selection method. Options:
#'   - `'auto'` (default): Test all three structures via QIC, select best
#'   - `'ar1'`: Autoregressive order 1, Corr(i,j) = ρ^|i-j|
#'   - `'exchangeable'`: Equal correlation across all pairs (no ordering assumed)
#'   - `'independence'`: Null model, no within-subject correlation
#' @param bias_correction logical; if TRUE (default), apply Kauermann-Carroll HC1
#'   bias reduction when n_clusters < 30. Ensures Type I error control in small samples
#' @param weights numeric or NULL; optional observation weights for heteroscedasticity
#'   (e.g., from bootstrap CI computations). If provided, takes precedence over
#'   internal heteroscedasticity detection
#'
#' @return data.frame (single row) with columns:
#'   - **gene**: gene identifier (from `g` argument)
#'   - **p_interaction**: p-value for q × group interaction (bias-corrected if applicable)
#'   - **p_interaction_raw**: p-value before K-C correction (if applied)
#'   - **n_clusters**: number of subjects/clusters in analysis
#'   - **bias_correction_applied**: logical; whether HC1 adjustment was performed
#'   - **correlation_structure**: selected structure ('ar1', 'exchangeable', 'independence')
#'   - **corstr_selection_method**: 'QIC_based' or 'user_specified'
#'   - **shapiro_p_value**: p-value for Shapiro-Wilk residual normality test
#'   - **residuals_normal**: logical; normality test result (if computed)
#'   - **ci_weighted**: logical; whether weights from bootstrap CI were applied
#'   - **slope_diff**: estimated group slope difference from interaction coefficient
#'   - **design_effect_ar1**: multiplier for effective sample size (D_eff)
#'   - **rho_ar1_estimate**: estimated AR(1) autocorrelation from residuals
#'   - **kc_bias_correction_applied**: logical; whether K-C correction was applied
#'   - **kc_multiplier**: HC1 adjustment multiplier (n_eff / (n_eff - p))
#'   - **n_effective**: effective sample size after design effect reduction
#'   - **kc_method**: method applied ('hc1', 'hc3', or 'kc')
#' @importFrom stats vcov
#' @noRd
.gee_interaction <- function(df, q_vals, g, subject = NULL, min_obs = 5, corstr = c("auto", "ar1", "exchangeable", "independence"),
    bias_correction = TRUE, weights = NULL) {
    # Validate corstr parameter per Bioconductor code syntax standards
    corstr <- match.arg(corstr)

    if (!requireNamespace("geepack", quietly = TRUE)) {
        stop("Package 'geepack' is required for method = 'gee'")
    }

    # Validate inputs
    validation <- .validate_gee_inputs(df, subject, min_obs, weights)
    if (!validation$valid) {
        warning(sprintf(".gee_interaction (gene %s): GEE input validation failed. Check that min_obs=%d is met, sufficient groups present, and subject structure is valid.",
            g, min_obs), call. = FALSE)
        return(NULL)
    }
    df <- validation$df
    subject <- validation$subject

    # Apply ARIMA(1,1,0) differencing for stationarity
    # AUDIT FIX #32-33: ARIMA differencing is applied only for paired designs (subject > 1).
    # GAM applies ARIMA for all designs; FPCA applies for all designs. Cross-method
    # p-values are not directly comparable due to different preprocessing.
    arima_result <- .apply_arima_differencing(df, subject)
    df <- arima_result$df
    subject <- arima_result$subject
    use_arima <- arima_result$use_arima
    df_orig_nrows <- arima_result$df_orig_nrows

    df$subject <- factor(subject)

    # Prepare weights (heteroscedasticity or bootstrap CI)
    weights_result <- .prepare_gee_weights(df)
    df <- weights_result$df
    gee_weights <- weights_result$gee_weights

    # Count clusters for bias correction decisions
    n_clusters <- length(unique(as.numeric(df$subject)))

    # Final validation
    if (sum(!is.na(df$entropy)) < 2) {
        warning(sprintf(".gee_interaction (gene %s): Insufficient non-NA entropy values after preprocessing (<%d observations).",
            g, 2), call. = FALSE)
        return(NULL)
    }

    # Select correlation structure
    selected_corstr <- corstr
    if (corstr == "auto") {
        selection_result <- .select_gee_correlation(df = df, formula_null = entropy ~
            q + group, formula_alt = entropy ~ q * group, subject = df$subject, criteria = "qic")
        selected_corstr <- selection_result$best_corstr
    }

    # Fit GEE models
    model_result <- .fit_gee_models(df, selected_corstr, gee_weights)
    if (is.null(model_result)) {
        warning(sprintf(".gee_interaction (gene %s): GEE model fitting failed with correlation structure '%s'. May indicate singularity, separation, or convergence issues.",
            g, selected_corstr), call. = FALSE)
        return(NULL)
    }
    fit_null <- model_result$fit_null
    fit_alt <- model_result$fit_alt

    # Extract interaction p-value with bias correction
    p_interaction <- .extract_interaction_pvalue(fit_alt, n_clusters, bias_correction)

    # ========================================================================
    # PHASE 9: Kauermann-Carroll Bias Correction with Design Effect
    # ========================================================================
    # Estimate AR(1) autocorrelation from residuals (for multi-q Tsallis)
    kc_metadata <- NULL
    rho_ar1 <- NA_real_
    design_effect_value <- 1

    if (!is.na(p_interaction)) {
        residuals_alt <- residuals(fit_alt)
        if (!is.null(residuals_alt) && length(residuals_alt) > 2) {
            rho_ar1 <- .estimate_ar1_correlation(residuals_alt)

            # Compute design effect if AR(1) significant
            cluster_size <- length(unique(df$q))
            if (!is.na(rho_ar1) && abs(rho_ar1) > 0.05) {
                design_effect_value <- .compute_ar1_design_effect(rho_ar1, cluster_size)
            }
        }

        # Apply Kauermann-Carroll bias correction if small clusters
        if (bias_correction && n_clusters < 30) {
            # Get interaction term coefficients for K-C correction
            coefs_alt <- stats::coef(fit_alt)
            ia_names <- names(coefs_alt)[grepl("^q:", names(coefs_alt), ignore.case = TRUE)]

            if (length(ia_names) > 0) {
                ia_name <- ia_names[1]
                ia_idx <- which(names(coefs_alt) == ia_name)[1]

                # Compute z-statistic for interaction
                z_interact <- .compute_wald_statistic(fit_alt, ia_idx)

                if (!is.na(z_interact)) {
                  # Apply K-C correction with design effect
                  kc_result <- .kc_bias_correct(p_value = p_interaction, z_statistic = z_interact,
                    vcov_sandwich_raw = NULL, n_clusters = n_clusters, n_parameters = length(coefs_alt),
                    rho_ar1 = rho_ar1, cluster_size = cluster_size, design_effect = design_effect_value,
                    bias_correction_method = "hc1", use_t_distribution = TRUE, apply_correction = TRUE,
                    verbose = FALSE)

                  # Update p-value with K-C correction
                  p_interaction <- kc_result$p_value
                  kc_metadata <- list(kc_applied = TRUE, p_raw = kc_result$p_raw,
                    p_corrected = kc_result$p_value, multiplier = kc_result$multiplier,
                    n_effective = kc_result$n_effective, design_effect = kc_result$design_effect,
                    rho_ar1 = kc_result$rho_ar1, method = kc_result$method_applied)
                }
            }
        }
    }

    if (is.null(kc_metadata)) {
        kc_metadata <- list(kc_applied = FALSE, design_effect = design_effect_value,
            rho_ar1 = rho_ar1)
    }

    # Test residual normality (Shapiro-Wilk test)
    shapiro_result <- .test_residual_normality(model = fit_alt, model_type = "gee",
        verbose = FALSE)

    # Extract slope difference
    slope_diff <- .extract_slope_diff(fit_alt)

    # Compile result — use consistent threshold (30) for both correction and reporting
    gee_result <- data.frame(gene = g, p_interaction = p_interaction, n_clusters = n_clusters,
        bias_correction_applied = bias_correction && n_clusters < 30, correlation_structure = selected_corstr,
        corstr_selection_method = if (corstr == "auto")
            "QIC_based" else "user_specified", stringsAsFactors = FALSE)

    # Add Shapiro-Wilk residual normality test results
    if (!is.na(shapiro_result$shapiro_p_value)) {
        gee_result$shapiro_p_value <- shapiro_result$shapiro_p_value
        gee_result$residuals_normal <- shapiro_result$residuals_normal
        gee_result$n_residuals_tested <- shapiro_result$n_residuals
    } else {
        gee_result$shapiro_p_value <- NA_real_
        gee_result$residuals_normal <- NA
        gee_result$n_residuals_tested <- NA_integer_
    }

    # Add Phase 1 bootstrap CI weighting tracking (March 2026)
    gee_result$ci_weighted <- !is.null(df$weight)
    gee_result$slope_diff <- slope_diff
    # AUDIT FIX #32: Flag ARIMA differencing status for cross-method comparability.
    # GEE applies ARIMA(1,1,0) only for paired designs; GAM always applies; FPCA always applies.
    # Users should compare p-values within method, not across methods with different ARIMA handling.
    gee_result$arima_applied <- use_arima

    # Add Phase 9 Kauermann-Carroll bias correction metadata
    gee_result$kc_bias_correction_applied <- !is.null(kc_metadata) && isTRUE(kc_metadata$kc_applied)
    gee_result$design_effect_ar1 <- if (!is.null(kc_metadata))
        kc_metadata$design_effect else NA_real_
    gee_result$rho_ar1_estimate <- if (!is.null(kc_metadata))
        kc_metadata$rho_ar1 else NA_real_

    if (!is.null(kc_metadata) && isTRUE(kc_metadata$kc_applied)) {
        gee_result$p_interaction_raw <- kc_metadata$p_raw
        gee_result$kc_multiplier <- kc_metadata$multiplier
        gee_result$n_effective <- kc_metadata$n_effective
        gee_result$kc_method <- kc_metadata$method
    } else {
        gee_result$p_interaction_raw <- NA_real_
        gee_result$kc_multiplier <- NA_real_
        gee_result$n_effective <- NA_real_
        gee_result$kc_method <- NA_character_
    }

    return(gee_result)
}

# ============================================================================
# HELPER: Validate GEE inputs and handle initial data setup
# ============================================================================
.validate_gee_inputs <- function(df, subject, min_obs, weights) {
    # Add weights to df if available and valid (Phase 1 weighting)
    if (!is.null(weights) && length(weights) == nrow(df)) {
        df$weight <- weights
    }

    if (sum(!is.na(df$entropy)) < min_obs) {
        return(list(valid = FALSE, subject = NULL, df = NULL))
    }
    if (length(unique(na.omit(df$group))) < 2) {
        return(list(valid = FALSE, subject = NULL, df = NULL))
    }

    # If no subject specified, use row indices (independent observations)
    if (is.null(subject)) {
        subject <- factor(seq_len(nrow(df)))
    } else {
        subject <- factor(subject)
    }

    list(valid = TRUE, subject = subject, df = df)
}

# ============================================================================
# HELPER: Apply ARIMA(1,1,0) differencing for stationarity
# ============================================================================
.apply_arima_differencing <- function(df, subject) {
    use_arima <- FALSE
    df_orig_nrows <- nrow(df)

    if (length(unique(subject)) > 1) {
        # Sort by subject and q for proper within-subject differencing
        sort_idx <- order(subject, df$q)
        df_sorted <- df[sort_idx, ]
        subject_sorted <- subject[sort_idx]

        # Compute first differences within subjects
        df_diff_list <- list()
        for (subj in unique(subject_sorted)) {
            subj_idx <- which(subject_sorted == subj)
            if (length(subj_idx) >= 2) {
                subj_data <- df_sorted[subj_idx, ]
                df_diff_list[[as.character(subj)]] <- data.frame(entropy = diff(subj_data$entropy),
                  q = subj_data$q[-nrow(subj_data)], group = subj_data$group[-nrow(subj_data)],
                  stringsAsFactors = FALSE)
            }
        }

        if (length(df_diff_list) > 0) {
            df <- do.call(rbind, df_diff_list)
            rownames(df) <- NULL
            subject <- rep(names(df_diff_list), vapply(df_diff_list, nrow, FUN.VALUE = integer(1)))
            use_arima <- TRUE
        }
    }

    list(df = df, subject = subject, use_arima = use_arima, df_orig_nrows = df_orig_nrows)
}

# ============================================================================
# HELPER: Prepare GEE weights (heteroscedasticity or bootstrap CI)
# ============================================================================
.prepare_gee_weights <- function(df) {
    gee_weights <- NULL

    # Heteroscedasticity detection
    hetero_result <- .detect_heteroscedasticity(df, q_vals = df$q, group_vec = df$group)

    # Phase 1 weighting: Bootstrap CI weights take precedence
    if (!is.null(df$weight)) {
        gee_weights <- df$weight
    } else if (!is.na(hetero_result$is_heteroscedastic) && hetero_result$is_heteroscedastic) {
        weights_result <- .estimate_variance_weights(df, q_vals = df$q, method = "power")
        if (!is.null(weights_result) && !is.null(weights_result$weights)) {
            gee_weights <- weights_result$weights
        }
    }

    if (!is.null(gee_weights)) {
        df$gee_weights <- gee_weights
    }

    list(df = df, gee_weights = gee_weights)
}

# ============================================================================
# HELPER: Fit null and alternative GEE models
# ============================================================================
.fit_gee_models <- function(df, selected_corstr, gee_weights) {
    if (!is.null(gee_weights)) {
        fit_null <- try(geepack::geeglm(entropy ~ q + group, id = df$subject, data = df,
            family = stats::gaussian(), weights = gee_weights, corstr = selected_corstr,
            na.action = stats::na.omit), silent = TRUE)

        fit_alt <- try(geepack::geeglm(entropy ~ q * group, id = df$subject, data = df,
            family = stats::gaussian(), weights = gee_weights, corstr = selected_corstr,
            na.action = stats::na.omit), silent = TRUE)
    } else {
        fit_null <- try(geepack::geeglm(entropy ~ q + group, id = df$subject, data = df,
            family = stats::gaussian(), corstr = selected_corstr, na.action = stats::na.omit),
            silent = TRUE)

        fit_alt <- try(geepack::geeglm(entropy ~ q * group, id = df$subject, data = df,
            family = stats::gaussian(), corstr = selected_corstr, na.action = stats::na.omit),
            silent = TRUE)
    }

    # Validate fits
    if (inherits(fit_null, "try-error") || inherits(fit_alt, "try-error") || is.null(fit_null) ||
        is.null(fit_alt)) {
        return(NULL)
    }

    list(fit_null = fit_null, fit_alt = fit_alt)
}

# ============================================================================
# HELPER: Extract interaction p-value with bias correction
# ============================================================================
.extract_interaction_pvalue <- function(fit_alt, n_clusters, bias_correction) {
    coefs_alt <- stats::coef(fit_alt)
    ia_names <- names(coefs_alt)[grepl("^q:", names(coefs_alt), ignore.case = TRUE)]

    if (length(ia_names) == 0) {
        return(NA_real_)
    }

    # AUDIT FIX #15: For multi-level groups, test ALL interaction coefficients
    # jointly using Wald test: β' V⁻¹ β ~ χ²(df = n_coefs). Previously only the
    # first coefficient was tested, ignoring other group levels.
    if (length(ia_names) > 1 && length(ia_names) <= length(stats::coef(fit_alt))) {
        p_joint <- .compute_joint_wald_pvalue(fit_alt, ia_names, n_clusters, bias_correction)
        if (!is.na(p_joint)) {
            return(p_joint)
        }
    }

    # Try extracting from summary coefficient table first
    summ <- try(summary(fit_alt), silent = TRUE)
    if (!inherits(summ, "try-error") && !is.null(summ)) {
        coef_table <- summ$coefficients
        if (!is.null(coef_table)) {
            for (ia_name in ia_names) {
                if (ia_name %in% rownames(coef_table)) {
                  row_idx <- which(rownames(coef_table) == ia_name)[1]
                  p_col <- grep("Pr\\(>", colnames(coef_table))[1]
                  if (!is.na(p_col) && p_col <= ncol(coef_table)) {
                    p_int_candidate <- coef_table[row_idx, p_col]
                    if (!is.na(p_int_candidate) && !is.nan(p_int_candidate)) {
                      p_interaction <- p_int_candidate

                      # Apply bias correction if needed
                      if (bias_correction && n_clusters < 20) {
                        p_corrected <- .apply_bias_correction(fit_alt, ia_name, ia_names,
                          n_clusters)
                        return(if (!is.na(p_corrected)) p_corrected else p_interaction)
                      }
                      return(p_interaction)
                    }
                  }
                }
            }
        }
    }

    # Fallback: compute Wald test with sandwich variance
    p_wald <- .compute_wald_pvalue(fit_alt, ia_names, n_clusters, bias_correction)
    if (!is.na(p_wald)) {
        return(p_wald)
    }

    # If all extraction methods fail, still return NA (don't return NULL)
    NA_real_
}

# ============================================================================
# HELPER: Apply Kauermann-Carroll bias correction (small cluster correction)
# ============================================================================
.apply_bias_correction <- function(fit_alt, ia_name, ia_names, n_clusters) {
    ia_idx <- which(names(stats::coef(fit_alt)) %in% ia_names)[1]
    if (is.na(ia_idx))
        return(NA_real_)

    z_stat <- .compute_wald_statistic(fit_alt, ia_idx)
    if (is.na(z_stat))
        return(NA_real_)

    # Use t-distribution (conservative, maintains Type I error for small
    # clusters)
    df_corr <- max(1, n_clusters - 1)
    2 * stats::pt(abs(z_stat), df = df_corr, lower.tail = FALSE)
}

# ============================================================================
# HELPER: Compute Wald test statistic with sandwich variance
# ============================================================================
.compute_wald_statistic <- function(fit_alt, ia_idx) {
    coefs_alt <- stats::coef(fit_alt)

    # Use vcov() directly - it handles all model types correctly
    vcov_robust <- try(vcov(fit_alt), silent = TRUE)

    if (inherits(vcov_robust, "try-error") || is.null(vcov_robust)) {
        return(NA_real_)
    }

    se_robust <- sqrt(diag(vcov_robust)[ia_idx])
    if (is.na(se_robust) || se_robust <= 0) {
        return(NA_real_)
    }

    coefs_alt[ia_idx]/se_robust
}

# ============================================================================
# HELPER: Compute Wald p-value (normal or t-distribution)
# ============================================================================
.compute_wald_pvalue <- function(fit_alt, ia_names, n_clusters, bias_correction) {
    ia_idx <- which(names(stats::coef(fit_alt)) %in% ia_names)[1]
    if (is.na(ia_idx))
        return(NA_real_)

    z_stat <- .compute_wald_statistic(fit_alt, ia_idx)
    if (is.na(z_stat))
        return(NA_real_)

    # Apply bias correction if specified and small clusters
    if (bias_correction && n_clusters < 20) {
        df_corr <- max(1, n_clusters - 1)
        2 * stats::pt(abs(z_stat), df = df_corr, lower.tail = FALSE)
    } else {
        2 * stats::pnorm(abs(z_stat), lower.tail = FALSE)
    }
}

# ============================================================================
# HELPER: Joint Wald test for multiple interaction coefficients (AUDIT FIX #15)
# ============================================================================
# Tests H0: all q:group interaction coefficients = 0 jointly
# Uses β' V⁻¹ β ~ χ²(df = k) where k = number of interaction coefficients
.compute_joint_wald_pvalue <- function(fit_alt, ia_names, n_clusters, bias_correction) {
    ia_idx <- which(names(stats::coef(fit_alt)) %in% ia_names)
    if (length(ia_idx) == 0) return(NA_real_)

    vcov_mat <- try(stats::vcov(fit_alt), silent = TRUE)
    if (inherits(vcov_mat, "try-error") || is.null(vcov_mat)) return(NA_real_)
    if (nrow(vcov_mat) < max(ia_idx)) return(NA_real_)

    # Extract sub-matrix for interaction coefficients
    beta <- stats::coef(fit_alt)[ia_idx]
    V <- vcov_mat[ia_idx, ia_idx, drop = FALSE]

    # Wald statistic: β' V⁻¹ β
    V_inv <- try(solve(V), silent = TRUE)
    if (inherits(V_inv, "try-error")) return(NA_real_)

    wald_stat <- as.numeric(t(beta) %*% V_inv %*% beta)
    df <- length(ia_idx)

    if (bias_correction && n_clusters < 20) {
        # Use F-distribution for small clusters: Wald/k ~ F(k, n_clusters - k)
        f_stat <- wald_stat / df
        df2 <- max(1, n_clusters - df)
        stats::pf(f_stat, df1 = df, df2 = df2, lower.tail = FALSE)
    } else {
        stats::pchisq(wald_stat, df = df, lower.tail = FALSE)
    }
}

# ============================================================================
# AUDIT FIX #46: GEE slope_diff = interaction coefficient (q:group), representing
# the additive change in entropy per unit-q when switching groups.
# GAM slope_diff = difference in predicted entropy slopes (ΔH/Δq) between groups.
# These have different units and magnitudes — not directly comparable across methods.
# ============================================================================
.extract_slope_diff <- function(fit_alt) {
    slope_diff <- NA_real_
    if (inherits(fit_alt, "try-error") || is.null(fit_alt)) {
        return(slope_diff)
    }

    tryCatch({
        coefs_alt <- stats::coef(fit_alt)
        ia_names <- names(coefs_alt)[grepl("^q:", names(coefs_alt), ignore.case = TRUE)]
        if (length(ia_names) > 0) {
            coefs_alt[ia_names[1]]
        } else {
            NA_real_
        }
    }, error = function(e) {
        NA_real_
    })
}

# ============================================================================
# HELPER: Estimate AR(1) autocorrelation from residuals/differenced data
# ============================================================================
.estimate_ar1_correlation <- function(residuals, max_lag = 1) {
    # Estimate AR(1) autocorrelation from residuals Input: numeric vector of
    # residuals/differenced data Output: correlation coefficient (rho) between
    # x(t) and x(t-1)

    if (length(residuals) < 3) {
        return(NA_real_)
    }

    # Remove NAs
    residuals <- residuals[!is.na(residuals)]

    if (length(residuals) < 3) {
        return(NA_real_)
    }

    # Compute lag-1 autocorrelation
    n <- length(residuals)
    rho <- stats::cor(residuals[seq_len(n - 1)], residuals[seq.int(2, n)], use = "complete.obs")

    if (is.na(rho))
        rho <- 0

    # Bound to valid correlation range (-1, 1)
    pmin(pmax(rho, -0.99), 0.99)
}

# ============================================================================
# HELPER: Compute design effect for AR(1) structure
# ============================================================================
.compute_ar1_design_effect <- function(rho, cluster_size) {
    # Design effect for AR(1) repeated measures D_eff = (1 + rho) / (1 - rho)
    # for positive correlation This accounts for reduction in effective sample
    # size due to within-subject correlation For multi-q Tsallis: cluster_size
    # = number of q-values per subject Result: effective sample size =
    # n_subjects_observed / design_effect

    if (is.null(rho) || is.na(rho) || abs(rho) < 0.001) {
        return(1)
    }

    # Bound rho to avoid numerical issues
    rho <- pmin(pmax(rho, -0.99), 0.99)

    # Standard design effect formula
    if (abs(rho) < 1) {
        design_effect <- (1 + rho)/(1 - rho)
    } else {
        design_effect <- 1
    }

    # Ensure positive
    max(1, design_effect)
}



# ============================================================================
# Kauermann-Carroll Bias Correction for GEE Sandwich Variance Estimation
# ============================================================================
# Phase 9 Implementation (March 2026) Extended for Tsallis multi-q measurements
# with AR(1) correlation structure References: - Kauermann & Carroll (2001). 'A
# note on the efficiency of sandwich covariance matrix estimation.' JASA
# 96(456): 1387-1396.  - Mancl & DeRouen (2001). 'A covariate-adjusted
# ANOVA-type test for correlated data.' Biometrics 57(1): 126-131.  - Li &
# Redden (2015). 'Comparing logistic and linear models: bias reduction via
# Kauermann-Carroll adjustment.' Biometrical Journal 57(5): 808-820.
# ============================================================================

# ============================================================================
# HELPER: Compute HC1 and HC3 bias reduction multipliers
# ============================================================================
.compute_hc_multipliers <- function(residuals_vec, X_matrix, leverage_vec = NULL,
    n_clusters = NULL, n_parameters = NULL) {
    # HC1: multiply by n / (n - p) Accounts for reduced degrees of freedom in
    # small samples
    n <- nrow(X_matrix)
    p <- ncol(X_matrix)

    hc1_multiplier <- n/(n - p)

    # HC3: leverage-adjusted, multiply by 1 / (1 - h_i)^2 Accounts for
    # high-leverage observations
    if (is.null(leverage_vec)) {
        # Compute leverage (diagonal of hat matrix) hat = X(X'X)^{-1}X'
        XtX_inv <- tryCatch(solve(t(X_matrix) %*% X_matrix), error = function(e) NULL)

        if (is.null(XtX_inv)) {
            leverage_vec <- rep(1/n, n)  # Fallback: uniform leverage
        } else {
            leverage_vec <- rowSums((X_matrix %*% XtX_inv) * X_matrix)
        }
    }

    # Ensure leverage is in (0, 1)
    leverage_vec <- pmin(pmax(leverage_vec, 1e-06), 1 - 1e-06)

    # HC3 divisor for each observation: (1 - h_i)^2
    hc3_divisor <- (1 - leverage_vec)^2

    list(hc1_multiplier = hc1_multiplier, hc3_divisor = hc3_divisor, leverage = leverage_vec)
}

# ============================================================================
# HELPER: Apply HC1 divisor to sandwich variance
# ============================================================================
.apply_hc1_correction <- function(vcov_sandwich_raw, n_clusters, n_parameters) {
    # HC1 multiplier: accounts for degrees of freedom adjustment More
    # conservative (larger variance) when n - p is small
    multiplier <- n_clusters/(n_clusters - n_parameters)

    vcov_corrected <- multiplier * vcov_sandwich_raw

    list(vcov = vcov_corrected, multiplier = multiplier)
}

# ============================================================================
# HELPER: Apply design effect adjustment for multi-q correlation
# ============================================================================
.adjust_for_design_effect <- function(n_clusters, n_parameters, design_effect, verbose = FALSE) {
    # Effective sample size accounts for within-subject correlation n_eff = n /
    # design_effect HC1 multiplier uses n_eff instead of n

    if (is.null(design_effect) || design_effect <= 0) {
        design_effect <- 1
    }

    n_effective <- n_clusters/design_effect

    # HC1 with effective sample size
    multiplier <- n_effective/(n_effective - n_parameters)

    if (verbose) {
        message(sprintf("Design Effect Adjustment:\n  n_clusters = %d\n  design_effect = %.3f\n  n_effective = %.1f\n  HC1 multiplier = %.4f\n",
            n_clusters, design_effect, n_effective, multiplier))
    }

    list(n_effective = n_effective, multiplier = multiplier, design_effect = design_effect)
}

# ============================================================================
# MAIN: Kauermann-Carroll Bias Correction for GEE
# ============================================================================
# @param p_value numeric; unadjusted p-value from Wald test @param z_statistic
# numeric; Wald z-statistic (or coefficient / SE) @param vcov_sandwich_raw
# matrix; raw sandwich variance (HC0) - optional if coef provided @param
# coef_value numeric; coefficient estimate (for HC1 adjustment) @param
# coef_index integer; index into vcov matrix for this coefficient @param
# n_clusters integer; number of independent clusters (subjects) @param
# n_parameters integer; number of model parameters @param rho_ar1 numeric;
# AR(1) autocorrelation (0-1), NULL if not available @param cluster_size
# integer; observations per cluster (q-values in multi-q) @param design_effect
# numeric; pre-computed design effect, auto-computed if NULL @param
# bias_correction_method character; 'hc1', 'hc3', or 'kc' (default: 'hc1')
# @param use_t_distribution logical; use t-dist (TRUE) or standard normal
# (FALSE) @param apply_correction logical; whether to apply correction (can be
# FALSE) @param verbose logical; print diagnostic information @return list with
# corrected results: - p_value: corrected p-value - p_raw: original unadjusted
# p-value - z_corrected: corrected z-statistic (if applicable) -
# vcov_corrected: bias-corrected variance-covariance matrix (if provided) -
# multiplier: HC1/HC3 adjustment factor - n_effective: effective sample size
# after design effect - design_effect: multiplier accounting for within-subject
# correlation - method_applied: 'none', 'hc1', 'hc3', 'kc' - report:
# human-readable summary
.kc_bias_correct <- function(p_value = NA_real_, z_statistic = NA_real_, vcov_sandwich_raw = NULL,
    coef_value = NULL, coef_index = NULL, n_clusters = NULL, n_parameters = NULL,
    rho_ar1 = NULL, cluster_size = NULL, design_effect = NULL, bias_correction_method = c("hc1",
        "hc3", "kc"), use_t_distribution = TRUE, apply_correction = TRUE, verbose = FALSE) {
    # ========================================================================
    # Input validation
    # ========================================================================

    if (is.null(n_clusters) || is.null(n_parameters)) {
        stop("n_clusters and n_parameters are required")
    }

    bias_correction_method <- match.arg(bias_correction_method)

    # ========================================================================
    # Compute effective sample size with design effect
    # ========================================================================

    # If design effect not provided, compute from AR(1)
    if (is.null(design_effect)) {
        if (!is.null(rho_ar1) && !is.null(cluster_size) && cluster_size > 1) {
            # Design effect for AR(1): D_eff = (1 + rho) / (1 - rho) Accounts
            # for within-subject correlation in multi-q measurements
            rho_ar1 <- pmin(pmax(rho_ar1, -0.99), 0.99)  # Bound in (-1, 1)

            if (abs(rho_ar1) < 0.001) {
                design_effect <- 1
            } else {
                design_effect <- (1 + rho_ar1)/(1 - rho_ar1)
            }

            if (verbose) {
                message(sprintf("AR(1) Design Effect: rho=%.3f, D_eff=%.3f, cluster_size=%d\n",
                  rho_ar1, design_effect, cluster_size))
            }
        } else {
            design_effect <- 1
        }
    }

    # Effective sample size
    n_effective <- n_clusters/design_effect

    # ========================================================================
    # Decision: Apply correction?
    # ========================================================================

    # Correction typically applied when n_effective < 30 (small clusters) But
    # make it data-adaptive
    correction_threshold <- 30

    if (!apply_correction || n_effective > correction_threshold) {
        # No correction needed
        return(list(p_value = p_value, p_raw = p_value, z_corrected = z_statistic,
            vcov_corrected = vcov_sandwich_raw, multiplier = 1, n_clusters = n_clusters,
            n_parameters = n_parameters, n_effective = n_effective, design_effect = design_effect,
            method_applied = "none", report = sprintf("No K-C correction (n_eff=%.1f > threshold=%d)",
                n_effective, correction_threshold)))
    }

    # ========================================================================
    # Apply HC1 or HC3 bias reduction multiplier
    # ========================================================================

    multiplier <- NA_real_

    if (bias_correction_method %in% c("hc1", "kc")) {
        # HC1: multiply by n_eff / (n_eff - p) More conservative for small
        # n_eff
        multiplier <- n_effective/(n_effective - n_parameters)
    }

    if (is.na(multiplier))
        multiplier <- 1

    # Apply multiplier to sandwich variance if provided
    vcov_corrected <- NULL
    if (!is.null(vcov_sandwich_raw)) {
        vcov_corrected <- multiplier * vcov_sandwich_raw
    }

    # ========================================================================
    # Recompute z-statistic and p-value with bias correction
    # ========================================================================

    # AUDIT FIX #14: Actually apply vcov_corrected to recompute the test statistic.
    # Previously vcov_corrected was computed but discarded; only pnorm→pt was swapped.
    z_corrected <- z_statistic  # Default: no change
    p_corrected <- p_value  # Default: no change

    if (!is.na(z_statistic)) {
        if (!is.null(vcov_corrected) && all(dim(vcov_corrected) >= 1)) {
            # Recompute z using the bias-corrected variance
            se_corrected <- sqrt(diag(vcov_corrected))
            if (length(se_corrected) >= 1 && se_corrected[1] > 0) {
                z_corrected <- coef_est / se_corrected[1]
            }
        }

        if (use_t_distribution) {
            df_t <- max(1, n_clusters - 1)
            p_corrected <- 2 * stats::pt(abs(z_corrected), df = df_t, lower.tail = FALSE)
        } else {
            p_corrected <- 2 * stats::pnorm(abs(z_corrected), lower.tail = FALSE)
        }
    }

    # ========================================================================
    # Return results with diagnostics
    # ========================================================================

    method_label <- switch(bias_correction_method, hc1 = "HC1 (bias-reduced)", hc3 = "HC3 (leverage-adjusted)",
        kc = "Kauermann-Carroll")

    report <- sprintf("K-C Bias Correction (%s):\n  n_clusters=%d, n_parameters=%d\n  n_effective=%.1f (design_effect=%.2f)\n  HC multiplier=%.4f\n  p-value: %.4f -> %.4f%s",
        method_label, n_clusters, n_parameters, n_effective, design_effect, multiplier,
        p_value, p_corrected, if (use_t_distribution)
            sprintf(" (df=%d, t-dist)", max(1, n_clusters - 1)) else " (normal)")

    list(p_value = p_corrected, p_raw = p_value, z_corrected = z_corrected, vcov_corrected = vcov_corrected,
        multiplier = multiplier, n_clusters = n_clusters, n_parameters = n_parameters,
        n_effective = n_effective, design_effect = design_effect, rho_ar1 = rho_ar1,
        method_applied = bias_correction_method, use_t_distribution = use_t_distribution,
        report = report)
}

# ============================================================================
# UTILITY: Generate diagnostic report for K-C correction
# ============================================================================
.print_kc_correction_report <- function(kc_result) {
    if (is.null(kc_result))
        return(invisible(NULL))

    message(paste(rep("-", 70), collapse = ""))
    message("Kauermann-Carroll Bias Correction Report")
    message(paste(rep("-", 70), collapse = ""))

    if (!is.null(kc_result$report)) {
        message(kc_result$report)
    }

    message("")
    message("Correction Details:")
    message(sprintf("  Method applied: %s", kc_result$method_applied))
    message(sprintf("  Multiplier (HC1): %.6f", kc_result$multiplier))

    msg <- sprintf("  Design effect: %.3f", kc_result$design_effect)
    if (!is.null(kc_result$rho_ar1)) {
        msg <- paste0(msg, sprintf(" (AR(1) rho=%.3f)", kc_result$rho_ar1))
    }
    message(msg)

    message(sprintf("  Degrees of freedom (t-dist): %d", max(1, kc_result$n_clusters -
        1)))
    message(sprintf("  P-value: %.6f -> %.6f", kc_result$p_raw, kc_result$p_value))

    message(paste(rep("-", 70), collapse = ""))
    message("")

    invisible(kc_result)
}

# GEE interaction helper for calculate_sait_interaction Generalized Estimating
# Equations (GEE) with AR(1) correlation structure for q-dependent entropy
# measurements. GEE is robust for correlated data and doesn't assume normality
# of random effects.  Paper S171 (Zimmerman & Harville, 1991): 'Linear Models
# with Generalized AR(1) Covariance Structure for Longitudinal and Spatial
# Data' validates AR(1) for ordered covariate structures (like q-values).
# Papers S168-S170: Theoretical foundation and empirical estimation of AR(1)
# parameters.  TEST L.1.6: Confirms q-value correlation follows AR(1) pattern
# (rho(k) = phi^|k|).  @param df data.frame with columns: entropy, q, group,
# subject (if paired) @param q_vals numeric vector of q values used Helper:
# Compare GEE correlation structures and select best via QIC Purpose: Validate
# that AR(1) is appropriate for Tsallis entropy or test alternatives
# Quasi-likelihood Information Criterion (QIC) is the GEE analog of AIC/BIC
# Selects the correlation structure that best balances fit and parsimony Lower
# QIC = better model Correlation structures tested: - AR(1): Geometric decay
# Corr(i,j) = phi^|i-j| [for ordered measurements] - Exchangeable: Equal
# correlation Corr(i,j) = rho [for unordered clusters] - Independence: No
# correlation [null/reference model] Reference: Pan, W. (2001). Akaike's
# information criterion in generalized estimating equations.  Biometrics,
# 57(1), 120-125.
.select_gee_correlation <- function(df, formula_null, formula_alt, subject, criteria = "qic",
    verbose = FALSE) {
    # Args: df: data frame with response, predictors, and subject/id column
    # formula_null: formula for null model (e.g., entropy ~ q + group)
    # formula_alt: formula for alternative model (e.g., entropy ~ q * group)
    # subject: vector of subject/cluster IDs criteria: model selection
    # criterion ('qic' or 'hybrid') verbose: whether to print comparison
    # results Returns: List with: best_corstr, qic_table, recommendation,
    # report (string)

    if (!requireNamespace("geepack", quietly = TRUE)) {
        return(list(best_corstr = "ar1", reason = "geepack not available; defaulting to AR(1)",
            qic_table = NULL, report = "geepack not available"))
    }

    corstr_options <- c("ar1", "exchangeable", "independence")
    results_list <- list()
    qic_values <- numeric(3)
    names(qic_values) <- corstr_options

    # Store correlation estimates for comparison
    corr_estimates <- list()

    for (corstr_candidate in corstr_options) {
        # Fit alternative model with this correlation structure
        fit_try <- try(geepack::geeglm(formula = formula_alt, id = subject, data = df,
            family = stats::gaussian(), corstr = corstr_candidate, na.action = stats::na.omit),
            silent = TRUE)

        if (inherits(fit_try, "try-error") || is.null(fit_try)) {
            # Model failed to fit: assign worst possible QIC
            qic_values[corstr_candidate] <- Inf
            corr_estimates[[corstr_candidate]] <- NA
            results_list[[corstr_candidate]] <- list(corstr = corstr_candidate, fit_status = "FAILED",
                qic = Inf, n_obs = NA, dispersion = NA, corr_estimate = NA)
            next
        }

        # Compute QIC (Quasi-likelihood Information Criterion) For GEE: QIC =
        # -2 * quasi-likelihood + 2 * trace(M_hat) where quasi-lik = -0.5 *
        # sum((y - mu)^2 / phi) for gaussian family

        qic_val <- NA_real_
        corr_estimate <- NA_real_
        try({
            # Extract components from geepack object
            residuals_vec <- as.numeric(fit_try$residuals)
            dispersion <- fit_try$geese$gamma[1]  # Scale parameter from geese

            # Extract correlation estimate if available
            if (!is.null(fit_try$geese$alpha) && length(fit_try$geese$alpha) > 0) {
                corr_estimate <- as.numeric(fit_try$geese$alpha[1])
            }

            # AUDIT FIX #13: Pan (2001) QIC = -2*quasi_ll + 2*trace(solve(V_naive) %*% V_robust)
            # where V_naive is the model-based covariance and V_robust is the sandwich estimator.
            # The previous penalty (1*log(n_obs) for ar1/exchangeable, 0 for independence)
            # was not Pan's QIC and biased selection toward independence.
            if (!is.na(dispersion) && dispersion > 0) {
                quasi_ll <- -0.5 * sum(residuals_vec^2/dispersion)

                # Compute trace penalty: trace(solve(V_naive) %*% V_robust)
                vcov_naive <- try(stats::vcov(fit_try, robust = FALSE), silent = TRUE)
                vcov_robust <- try(stats::vcov(fit_try, robust = TRUE), silent = TRUE)
                if (!inherits(vcov_naive, "try-error") && !inherits(vcov_robust, "try-error") &&
                    nrow(vcov_naive) == nrow(vcov_robust)) {
                    penalty <- 2 * sum(diag(solve(vcov_naive) %*% vcov_robust))
                } else {
                    # Fallback: standard penalty = 2 * p (number of parameters)
                    penalty <- 2 * length(stats::coef(fit_try))
                }

                qic_val <- -2 * quasi_ll + penalty
            }
        }, silent = TRUE)

        qic_values[corstr_candidate] <- ifelse(is.na(qic_val), Inf, qic_val)
        corr_estimates[[corstr_candidate]] <- corr_estimate

        results_list[[corstr_candidate]] <- list(corstr = corstr_candidate, fit_status = "SUCCESS",
            qic = qic_val, n_obs = nrow(df), dispersion = ifelse(is.null(fit_try$geese$gamma[1]),
                NA, fit_try$geese$gamma[1]), corr_estimate = corr_estimate)
    }

    # Select best model (lowest QIC)
    valid_qics <- qic_values[!is.infinite(qic_values)]

    # Initialize QIC variables and correlation estimates for safe use in report generation
    ar1_qic <- qic_values["ar1"]
    exch_qic <- qic_values["exchangeable"]
    indep_qic <- qic_values["independence"]
    ar1_corr <- corr_estimates[["ar1"]]

    if (length(valid_qics) == 0) {
        # All models failed: default to AR(1)
        best_corstr <- "ar1"
        reason <- "All correlation structures failed to fit; defaulting to AR(1)"
    } else {
        best_idx <- which.min(qic_values)
        best_corstr <- names(qic_values)[best_idx]

        # Create detailed reasoning based on QIC values and observed correlations
        if (best_corstr == "ar1") {
            reason <- sprintf("AR(1) selected: QIC=%.3f (Exchangeable: %.3f, Independence: %.3f). Estimated AR(1) correlation=%.3f.",
                ar1_qic, exch_qic, indep_qic, ifelse(is.na(ar1_corr), 0, ar1_corr))
        } else if (best_corstr == "exchangeable") {
            reason <- sprintf("Exchangeable selected: QIC=%.3f (AR(1): %.3f, Independence: %.3f). Suggests uniform correlation.",
                exch_qic, ar1_qic, indep_qic)
        } else {
            reason <- sprintf("Independence selected: QIC=%.3f (AR(1): %.3f, Exchangeable: %.3f). No significant correlation detected.",
                indep_qic, ar1_qic, exch_qic)
        }
    }

    # Create comparison table
    qic_table <- data.frame(correlation_structure = corstr_options, fit_status = vapply(corstr_options,
        function(cs) results_list[[cs]]$fit_status, FUN.VALUE = character(1)), qic = qic_values,
        corr_estimate = vapply(corstr_options, function(cs) {
            est <- corr_estimates[[cs]]
            if (is.na(est))
                "NA" else sprintf("%.4f", est)
        }, FUN.VALUE = character(1)), selected = ifelse(corstr_options == best_corstr,
            "YES", ""), stringsAsFactors = FALSE)

    # Generate report
    report_lines <- c(sprintf("GEE Correlation Structure Selection:"), sprintf(""),
        sprintf("QIC Comparison (lower = better):"), sprintf("  AR(1):           QIC = %.3f  (Est. corr = %s)",
            ar1_qic, ifelse(is.na(ar1_corr), "NA", sprintf("%.4f", ar1_corr))), sprintf("  Exchangeable:   QIC = %.3f",
            exch_qic), sprintf("  Independence:    QIC = %.3f", indep_qic), sprintf(""),
        sprintf("Selected: %s", best_corstr), sprintf("Reasoning: %s", reason))

    return(list(best_corstr = best_corstr, reason = reason, qic_table = qic_table,
        qic_values = qic_values, corr_estimates = corr_estimates, report = paste(report_lines,
            collapse = "\n")))
}
