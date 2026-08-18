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
#' 2. **No ARIMA Differencing**: the test targets the functional
#'    interaction on the original H(q) curve
#' 3. **Weight Preparation**: Apply heteroscedasticity weights or bootstrap CI weights
#' 4. **Correlation Structure Selection**: Choose AR(1), exchangeable, or independence via QIC
#' 5. **Model Fitting**: Fit null (main effects) and alternative (interaction) models
#' 6. **Interaction Testing**: Extract p-value with bias correction for small clusters
#' 7. **Small-Cluster Sandwich Correction**: HC1 bias reduction for n_clusters < 30
#'    (empirical correction inspired by Kauermann & Carroll 2001, combined
#'    with a t-reference; NOT a literal KC estimator)
#'
#' ## Key References
#'
#' - Zimmerman & Harville (1991): AR(1) for ordered covariate structures
#' - Kauermann & Carroll (2001): inspiration for small-cluster sandwich
#'   variance bias reduction (the implementation below is an EMPIRICAL
#'   small-cluster correction, not the literal KC estimator)
#' - Pan (2001): QIC model selection criterion for GEE
#' - Mancl & DeRouen (2001): Covariate-adjusted ANOVA-type tests with GEE
#'
#' ## Important Clarifications
#'
#' - **Correlation vs Random Effects**: GEE models within-subject correlation
#'   directly (no random intercepts like mixed models). For AR(1), the pattern
#'   Corr(q_i, q_j) = ρ^|i-j| accounts for ordered q-value measurements.
#'
#' - **Design Effect**: Multi-q measurements create within-cluster correlation,
#'   quantified descriptively as D_eff (finite-m AR(1) form). It
#'   is NOT used to rescale the sandwich variance (which would double-count the
#'   dependence); HC1 uses the number of clusters.
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
#' @param bias_correction logical; if TRUE (default), apply an empirical
#'   small-cluster sandwich (HC1) bias reduction when n_clusters < 30, with a
#'   t-reference (df = n_clusters - p). Inspired by Kauermann & Carroll (2001)
#'   but NOT a literal KC estimator. Ensures Type I error control in small samples
#' @param weights numeric or NULL; optional observation weights for heteroscedasticity
#'   (e.g., from bootstrap CI computations). If provided, takes precedence over
#'   internal heteroscedasticity detection
#'
#' @return data.frame (single row) with columns:
#'   - **gene**: gene identifier (from `g` argument)
#'   - **p_interaction**: p-value for q × group interaction (bias-corrected if applicable)
#'   - **p_interaction_raw**: p-value before small-cluster correction (if applied)
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
#'   - **kc_bias_correction_applied**: logical; whether the empirical
#'     small-cluster sandwich correction was applied (column name kept for
#'     backwards compatibility)
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

    # ========================================================================
    # REMOVED ARIMA(1,1,0) differencing.
    # ========================================================================
    # Differencing rows sorted by (subject, q) interleaves condition and q
    # (A(q1), B(q1), A(q2), B(q2), ...), so diff() mixes within-q condition
    # contrasts with between-q transitions. The estimand stops being the
    # difference between entropy curves and becomes a mixture of DeltaH terms.
    # The GEE is therefore fitted on the ORIGINAL entropy H(q). For corstr=
    # 'ar1', rows are ordered by subject -> condition (group) -> q so that
    # adjacent observations within a cluster are q-adjacent within a
    # condition, never interleaved across conditions.
    df$subject <- factor(subject)
    df <- df[order(as.character(df$subject), as.character(df$group), df$q), , drop = FALSE]
    use_arima <- FALSE

    # Prepare weights (heteroscedasticity or bootstrap CI)
    weights_result <- .prepare_gee_weights(df)
    df <- weights_result$df
    gee_weights <- weights_result$gee_weights

    # Count clusters for bias correction decisions
    n_clusters <- nlevels(df$subject)

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
    # PHASE 7: Kauermann-Carroll Bias Correction
    # When the interaction test is JOINT over several q:group
    # coefficients, the joint Wald test already carries its own small-sample
    # correction (.compute_joint_wald_pvalue uses an F approximation).
    # Overwriting the joint p-value with a K-C correction computed on a SINGLE
    # coefficient destroys the joint-test interpretation. K-C is only applied
    # when the model has exactly one interaction coefficient.
    # ========================================================================
    n_ia_coefs <- .count_interaction_coefs(fit_alt)
    if (n_ia_coefs == 1L) {
        kc_result <- .gee_apply_kc_correction(
            fit_alt = fit_alt, df = df,
            p_interaction = p_interaction,
            n_clusters = n_clusters,
            bias_correction = bias_correction
        )
        p_interaction <- kc_result$p_interaction
        kc_metadata <- kc_result$kc_metadata
    } else {
        kc_metadata <- list(kc_applied = FALSE, design_effect = 1, rho_ar1 = NA_real_,
            note = if (n_ia_coefs > 1L)
                "joint Wald test used; per-coefficient K-C skipped" else
                "no interaction coefficients")
    }

    # ========================================================================
    # PHASE 8: Assemble result row
    # ========================================================================
    gee_result <- .gee_assemble_result_row(
        g = g, p_interaction = p_interaction,
        n_clusters = n_clusters,
        bias_correction = bias_correction,
        selected_corstr = selected_corstr,
        corstr = corstr,
        fit_alt = fit_alt,
        df = df,
        use_arima = use_arima,
        kc_metadata = kc_metadata
    )

    return(gee_result)
}

# ============================================================================
# HELPER: Count q:group interaction coefficients in fitted GEE model
# ============================================================================
.count_interaction_coefs <- function(fit_alt) {
    if (inherits(fit_alt, "try-error") || is.null(fit_alt)) {
        return(0L)
    }
    coef_names <- names(stats::coef(fit_alt))
    length(coef_names[grepl("^q:", coef_names, ignore.case = TRUE)])
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
# Apply Kauermann-Carroll bias correction (July 2026 refactoring)
# ============================================================================
# Extracted from .gee_interaction() to reduce cyclomatic complexity (75→~35).
# Handles AR(1) rho estimation, design effect computation, and K-C correction.
#
#' @param fit_alt Fitted GEE alternative model (geeglm object)
#' @param df Data frame with subject column
#' @param p_interaction Raw interaction p-value before correction
#' @param n_clusters Number of clusters/subjects
#' @param bias_correction Logical; whether to apply correction
#' @return List with p_interaction (possibly corrected) and kc_metadata
#' @noRd
.gee_apply_kc_correction <- function(fit_alt, df, p_interaction, n_clusters,
    bias_correction) {
    kc_metadata <- NULL
    rho_ar1 <- NA_real_
    design_effect_value <- 1

    if (is.na(p_interaction)) {
        return(list(p_interaction = p_interaction,
            kc_metadata = list(kc_applied = FALSE, design_effect = 1, rho_ar1 = NA_real_)))
    }

    residuals_alt <- residuals(fit_alt)
    if (!is.null(residuals_alt) && length(residuals_alt) > 2) {
        # Estimate AR(1) rho within each cluster (subject) and pool via
        # Fisher z-transform
        subject_levels <- levels(df$subject)
        if (length(subject_levels) > 0 && length(residuals_alt) == nrow(df)) {
            rho_per_subject <- vapply(subject_levels, function(s) {
                idx <- which(df$subject == s)
                if (length(idx) >= 3) {
                    .estimate_ar1_correlation(residuals_alt[idx])
                } else NA_real_
            }, FUN.VALUE = numeric(1))
            rho_valid <- rho_per_subject[!is.na(rho_per_subject)]
            if (length(rho_valid) > 0) {
                z_vals <- atanh(pmin(pmax(rho_valid, -0.99), 0.99))
                rho_ar1 <- tanh(mean(z_vals))
            }
        } else {
            rho_ar1 <- .estimate_ar1_correlation(residuals_alt)
        }

        # Compute design effect if AR(1) significant
        cluster_size <- length(unique(df$q))
        if (!is.na(rho_ar1) && abs(rho_ar1) > 0.05) {
            design_effect_value <- .compute_ar1_design_effect(rho_ar1, cluster_size)
        }
    }

    # Apply Kauermann-Carroll bias correction if small clusters
    if (bias_correction && n_clusters < 30) {
        coefs_alt <- stats::coef(fit_alt)
        ia_names <- names(coefs_alt)[grepl("^q:", names(coefs_alt), ignore.case = TRUE)]

        if (length(ia_names) > 0) {
            ia_name <- ia_names[1]
            ia_idx <- which(names(coefs_alt) == ia_name)[1]
            z_interact <- .compute_wald_statistic(fit_alt, ia_idx)

            if (!is.na(z_interact)) {
                vcov_sandwich <- try(vcov(fit_alt), silent = TRUE)
                if (inherits(vcov_sandwich, "try-error")) vcov_sandwich <- NULL
                if (!is.null(vcov_sandwich) &&
                    (!is.matrix(vcov_sandwich) || any(!is.finite(vcov_sandwich)))) {
                    vcov_sandwich <- NULL
                }
                coef_value <- coefs_alt[ia_idx]

                kc_result <- .kc_bias_correct(p_value = p_interaction,
                    z_statistic = z_interact, vcov_sandwich_raw = vcov_sandwich,
                    coef_value = coef_value, coef_index = ia_idx,
                    n_clusters = n_clusters, n_parameters = length(coefs_alt),
                    rho_ar1 = rho_ar1, cluster_size = cluster_size,
                    design_effect = design_effect_value,
                    bias_correction_method = "hc1", use_t_distribution = TRUE,
                    apply_correction = TRUE, verbose = FALSE)

                p_interaction <- kc_result$p_value
                kc_metadata <- list(kc_applied = TRUE, p_raw = kc_result$p_raw,
                    p_corrected = kc_result$p_value, multiplier = kc_result$multiplier,
                    n_effective = kc_result$n_effective,
                    design_effect = kc_result$design_effect,
                    rho_ar1 = kc_result$rho_ar1, method = kc_result$method_applied)
            }
        }
    }

    if (is.null(kc_metadata)) {
        kc_metadata <- list(kc_applied = FALSE, design_effect = design_effect_value,
            rho_ar1 = rho_ar1)
    }

    list(p_interaction = p_interaction, kc_metadata = kc_metadata)
}

# ============================================================================
# Assemble GEE result row (July 2026 refactoring)
# ============================================================================
# Extracted from .gee_interaction() to reduce cyclomatic complexity (75→~35).
# Assembles all computed components into the final single-row data.frame.
#
#' @noRd
.gee_assemble_result_row <- function(g, p_interaction, n_clusters, bias_correction,
    selected_corstr, corstr, fit_alt, df, use_arima, kc_metadata) {
    
    # Shapiro-Wilk residual normality test
    shapiro_result <- .test_residual_normality(model = fit_alt, model_type = "gee",
        verbose = FALSE)
    
    # Extract slope difference
    slope_diff <- .extract_slope_diff(fit_alt)
    
    # Compile base result
    gee_result <- data.frame(
        gene = g,
        p_interaction = p_interaction,
        n_clusters = n_clusters,
        bias_correction_applied = bias_correction && n_clusters < 30,
        correlation_structure = selected_corstr,
        corstr_selection_method = if (corstr == "auto") "QIC_based" else "user_specified",
        stringsAsFactors = FALSE
    )
    
    # Shapiro-Wilk results
    if (!is.na(shapiro_result$shapiro_p_value)) {
        gee_result$shapiro_p_value <- shapiro_result$shapiro_p_value
        gee_result$residuals_normal <- shapiro_result$residuals_normal
        gee_result$n_residuals_tested <- shapiro_result$n_residuals
    } else {
        gee_result$shapiro_p_value <- NA_real_
        gee_result$residuals_normal <- NA
        gee_result$n_residuals_tested <- NA_integer_
    }
    
    # Bootstrap CI weighting and slope
    gee_result$ci_weighted <- !is.null(df$weight)
    gee_result$slope_diff <- slope_diff
    gee_result$arima_applied <- use_arima
    
    # K-C correction metadata
    gee_result$kc_bias_correction_applied <- !is.null(kc_metadata) &&
        isTRUE(kc_metadata$kc_applied)
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
    
    gee_result
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

    # For multi-level groups, test ALL interaction coefficients
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

                      # Apply bias correction if needed (consistent threshold with caller)
                      if (bias_correction && n_clusters < 30) {
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
    # clusters). df = n_clusters - p (estimated coefficients).
    df_corr <- max(1, n_clusters - length(stats::coef(fit_alt)))
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

    # Apply bias correction if specified and small clusters (threshold: 30, consistent)
    if (bias_correction && n_clusters < 30) {
        # Use df = n_clusters - p (estimated coefficients) instead
        # of n_clusters - 1 (anti-conservative for few clusters)
        df_corr <- max(1, n_clusters - length(stats::coef(fit_alt)))
        2 * stats::pt(abs(z_stat), df = df_corr, lower.tail = FALSE)
    } else {
        2 * stats::pnorm(abs(z_stat), lower.tail = FALSE)
    }
}

# ============================================================================
# HELPER: Joint Wald test for multiple interaction coefficients
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

    if (bias_correction && n_clusters < 30) {
        # Use F-distribution for small clusters: Wald/k ~ F(k, n_clusters - p)
        # where p = number of estimated coefficients. Validation (
        # tests/testthat/test-gee-small-sample.R) showed that
        # df2 = n_clusters - k is anti-conservative for few clusters; the
        # n_clusters - p reference is conservative but controls the type I.
        f_stat <- wald_stat / df
        p_model <- length(stats::coef(fit_alt))
        df2 <- max(1, n_clusters - p_model)
        stats::pf(f_stat, df1 = df, df2 = df2, lower.tail = FALSE)
    } else {
        stats::pchisq(wald_stat, df = df, lower.tail = FALSE)
    }
}

# ============================================================================
# GEE slope_diff = interaction coefficient (q:group), representing
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
# HELPER: Compute design effect for AR(1) structure (finite-m form)
# ============================================================================
# Uses the exact finite-m formula (Diggle et al. 2002, Crowder 1995):
#   D_eff = 1 + 2 * Σ_{k=1}^{m-1} (1 - k/m) * ρ^k
# This correctly accounts for edge effects in small clusters (few q-values),
# unlike the asymptotic (1+ρ)/(1-ρ) which overstates the design effect.
.compute_ar1_design_effect <- function(rho, cluster_size) {
    if (is.null(rho) || is.na(rho) || rho <= 0 || cluster_size <= 1) {
        return(1)
    }

    if (rho >= 1) {
        return(as.numeric(cluster_size))
    }

    # Bound rho to avoid numerical issues
    rho <- pmin(pmax(rho, -0.99), 0.99)

    # Finite-m AR(1) design effect with edge-effect weights
    summed <- 0
    for (k in seq_len(cluster_size - 1)) {
        lambda_k <- 1 - k / cluster_size
        summed <- summed + lambda_k * (rho^k)
    }

    d_eff <- 1 + 2 * summed
    d_eff <- pmax(1, d_eff)  # Ensure D_eff >= 1

    return(d_eff)
}



# ============================================================================
# Empirical Small-Cluster Sandwich Correction for GEE Variance Estimation
# ============================================================================
# Phase 9 Implementation (March 2026) Extended for Tsallis multi-q measurements
# with AR(1) correlation structure.
#
# AUDIT R4: this is NOT a literal Kauermann-Carroll estimator. It combines an
# HC1-style sandwich variance multiplier, a t-reference with degrees of
# freedom based on the number of clusters, and a descriptive AR(1) design
# effect. It should be described as an EMPIRICAL small-cluster sandwich
# correction (validated by Monte Carlo in test-gee-small-sample.R), not as
# "the KC correction". References: - Kauermann & Carroll (2001). 'A note on
# the efficiency of sandwich covariance matrix estimation.' JASA 96(456):
# 1387-1396.  - Mancl & DeRouen (2001). 'A covariate-adjusted ANOVA-type test
# for correlated data.' Biometrics 57(1): 126-131.  - Li & Redden (2015).
# 'Comparing logistic and linear models: bias reduction via
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
# MAIN: Empirical Small-Cluster Sandwich Correction for GEE
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

    # The GEE sandwich estimator already accounts for
    # within-cluster dependence in the variance. Scaling the variance again by
    # a design-effect-reduced sample size (n_eff = n_clusters / D_eff) counts
    # the dependence TWICE. The HC1 multiplier therefore uses the number of
    # CLUSTERS (empirical small-cluster correction, inspired by Kauermann &
    # Carroll 2001 but not the literal KC estimator). The design effect is
    # retained as descriptive metadata only and never enters the variance
    # multiplier.
    n_effective <- n_clusters

    # ========================================================================
    # Decision: Apply correction?
    # ========================================================================

    # Correction typically applied when n_clusters < 30 (small clusters) But
    # make it data-adaptive
    correction_threshold <- 30

    if (!apply_correction || n_effective > correction_threshold) {
        # No correction needed
        return(list(p_value = p_value, p_raw = p_value, z_corrected = z_statistic,
            vcov_corrected = vcov_sandwich_raw, multiplier = 1, n_clusters = n_clusters,
            n_parameters = n_parameters, n_effective = n_effective, design_effect = design_effect,
            method_applied = "none", report = sprintf("No small-cluster correction (n_eff=%.1f > threshold=%d)",
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
        # Guard against NaN/Inf in sandwich variance (singular fits)
        if (any(!is.finite(vcov_corrected))) {
            vcov_corrected <- NULL
        }
    }

    # ========================================================================
    # Recompute z-statistic and p-value with bias correction
    # ========================================================================

    # Actually apply vcov_corrected to recompute the test statistic.
    # Previously vcov_corrected was computed but discarded; only pnorm→pt was swapped.
    z_corrected <- z_statistic  # Default: no change
    p_corrected <- p_value  # Default: no change

    if (!is.na(z_statistic)) {
        if (!is.null(vcov_corrected) && all(dim(vcov_corrected) >= 1)) {
            # Recompute z using the bias-corrected variance
            se_corrected <- sqrt(diag(vcov_corrected))
            if (length(se_corrected) >= 1 && is.finite(se_corrected[1]) && se_corrected[1] > 0) {
                z_new <- coef_value / se_corrected[1]
                if (is.finite(z_new)) z_corrected <- z_new
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
        kc = "empirical small-cluster (KC-style)")

    report <- sprintf("Empirical small-cluster sandwich correction (%s):\n  n_clusters=%d, n_parameters=%d\n  n_effective=%.1f (design_effect=%.2f)\n  HC multiplier=%.4f\n  p-value: %.4f -> %.4f%s",
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
# of random effects.  Zimmerman & Harville (1991): 'Linear Models
# with Generalized AR(1) Covariance Structure for Longitudinal and Spatial
# Data' validates AR(1) for ordered covariate structures (like q-values).
# Grunwald, Hyndman & Tedesco (2000): theoretical foundation and empirical
# estimation of AR(1) parameters.  TEST L.1.6: Confirms q-value correlation follows AR(1) pattern
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

            # Pan (2001) QIC = -2*quasi_ll + 2*trace(solve(V_naive) %*% V_robust)
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
