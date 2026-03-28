# ===============================================================================
# HELPER FUNCTION: Setup and validate GAM data
# ===============================================================================
# Ensures 'group' is a factor and validates package dependencies
.setup_gam_data <- function(df) {
    # Ensure 'group' is a factor for 'by' argument in GAM/GAMM smooths
    # This is required for s(q, by = group, ...) to work correctly
    if (!requireNamespace("mgcv", quietly = TRUE)) {
        stop("Package 'mgcv' is required for method = 'gam'")
    }
    
    if (!is.null(df$group)) {
        df$group <- factor(df$group)
    }
    
    return(df)
}

# ===============================================================================
# HELPER FUNCTION: Select appropriate GAM family and link
# ===============================================================================
# Selects family based on bounded support and subject info
# GAMM does NOT support extended families, so falls back to gaussian for paired designs
.select_gam_family <- function(bounded_result, subject) {
    # GAMM COMPATIBILITY FIX (March 2026): mgcv::gamm() does NOT support extended families
    # (beta, gamma, Tweedie, etc.). For paired designs (subject != NULL -> uses gamm),
    # fall back to gaussian family instead of extended families.
    # Reference: C042/C043 (GAMM Tutorial, mgcv Documentation)
    
    use_bounded_family <- bounded_result$use_gamma
    family_gam <- bounded_result$family_obj
    inverse_link_fn <- bounded_result$inverse_link
    
    if (!is.null(subject)) {
        # For paired/mixed designs using gamm(), force gaussian family
        # CRITICAL FIX (March 2026): mgcv::gamm() does NOT support extended families
        # Warn user that bounded family selection is being overridden for reproducibility
        if (use_bounded_family) {
            warning(
                "[calculate_lm_interaction] GAMM with paired design detected. ",
                "mgcv::gamm() does not support extended families. ",
                "Forcing gaussian family. Results may be less accurate for bounded data.",
                call. = FALSE
            )
        }
        family_gam <- stats::gaussian()
        inverse_link_fn <- function(eta) eta
    }
    
    # If beta regression is selected, use stabilized entropy
    if (bounded_result$use_beta && !is.null(bounded_result$stabilized_df)) {
        df <- bounded_result$stabilized_df
    }
    
    return(list(use_bounded_family = use_bounded_family,
                family_gam = family_gam,
                inverse_link_fn = inverse_link_fn))
}

# ===============================================================================
# HELPER FUNCTION: Prepare weights for GAM/GAMM fitting
# ===============================================================================
# Handles heteroscedasticity detection and ARIMA differencing
.prepare_gam_weights <- function(df, q_vals, weights_input, hetero_result, subject) {
    # PHASE 1 WEIGHTING (March 2026): Use bootstrap CI weights if provided
    # These take precedence over heteroscedasticity-estimated weights
    gam_weights_original <- NULL
    
    if (!is.null(weights_input) && length(weights_input) == nrow(df)) {
        gam_weights_original <- weights_input  # Bootstrap CI weights for Phase 1
    } else if (!is.na(hetero_result$is_heteroscedastic) && hetero_result$is_heteroscedastic) {
        weights_result <- .estimate_variance_weights(df, q_vals, method = "power")
        if (!is.null(weights_result)) {
            gam_weights_original <- weights_result$weights  # Store original weights
        }
    }
    
    return(gam_weights_original)
}

# ===============================================================================
# HELPER FUNCTION: Handle ARIMA transformations and weight updates
# ===============================================================================
# Applies ARIMA(1,1,0) differencing and updates weights accordingly
.handle_arima_and_weights <- function(df, q_vals, subject, gam_weights_original) {
    # ARIMA(1,1,0) IMPLEMENTATION: Compute first differences for stationarity
    # Differencing removes monotone trend from Tsallis entropy, enabling valid AR(1) inference
    # This is applied when subject information is available (paired design)
    # CRITICAL: ARIMA is applied AFTER heteroscedasticity detection on original data
    
    use_arima <- FALSE
    gam_weights <- gam_weights_original
    n_samples_original <- nrow(df)
    
    if (!is.null(subject)) {
        # Only apply ARIMA differencing for paired designs (has subject info)
        arima_result <- .compute_arima_differences(df, q_vals, df$group, factor(subject))
        
        if (!is.null(arima_result) && nrow(arima_result$df) >= 3) {
            # Sufficient data for ARIMA differencing
            df <- arima_result$df
            use_arima <- TRUE
            
            # FIX: When ARIMA is applied, DON'T use weights computed on original data
            # Reason: Differencing changes the variance structure, weights would be invalid
            # Conservative approach: Better to lose efficiency than introduce bias
            gam_weights <- NULL
        }
    }
    
    # PHASE 1 WEIGHTING (March 2026): Set df$weight for ci_weighted flag tracking
    if (!is.null(gam_weights)) {
        df$weight <- gam_weights
    }
    
    return(list(df = df, use_arima = use_arima, gam_weights = gam_weights, n_samples = nrow(df)))
}

# ===============================================================================
# HELPER FUNCTION: Compute adaptive spline knots
# ===============================================================================
# Selects knot parameters based on sample size and data complexity
.compute_adaptive_knots <- function(df, q_vals, adaptive_knots) {
    uq_len <- length(unique(na.omit(q_vals)))
    
    # Adaptive knot selection: compute k based on gene's entropy curve complexity
    if (adaptive_knots) {
        k_q <- .adaptive_spline_knots_memo(entropy_vals = df$entropy, q_vals = q_vals,
                                                 n_q_unique = uq_len, min_k = 2, max_k = 10)
    } else {
        # Fallback to static knot selection
        k_q <- max(2, min(10, uq_len - 1))
    }
    
    return(list(k_q = k_q, uq_len = uq_len))
}

# ===============================================================================
# HELPER FUNCTION: Fit single GAMM with AR(1) model
# ===============================================================================
# Fits null and alternative GAMM models with AR(1) correlation
.fit_gamm_ar1_single <- function(formula, df, family_gam, gam_weights) {
    # Helper to fit GAMM with AR(1) correlation
    # Used internally by .fit_gamm_ar1
    
    if (!is.null(gam_weights)) {
        fit <- try(
            mgcv::gamm(formula = formula,
                      random = list(subject = ~1), 
                      correlation = nlme::corAR1(form = ~obs_seq|subject),
                      family = family_gam,
                      weights = gam_weights,
                      data = df),
            silent = TRUE
        )
    } else {
        fit <- try(
            mgcv::gamm(formula = formula,
                      random = list(subject = ~1), 
                      correlation = nlme::corAR1(form = ~obs_seq|subject),
                      family = family_gam,
                      data = df),
            silent = TRUE
        )
    }
    
    return(fit)
}

# ===============================================================================
# HELPER FUNCTION: Fit GAMM with AR(1) correlation (both models)
# ===============================================================================
# Handles GAMM fitting with AR(1) - both null and alternative models
.fit_gamm_ar1 <- function(df, family_gam, k_q_marginal, k_q_interaction, gam_weights) {
    # BUG FIX: Use smooth splines s() instead of poly() for actual GAM fitting
    # Adaptive spline basis with thin-plate (tp) for flexible curve fitting
    # PRIORITY 1: Try GAMM with AR(1) correlation (full autocorrelation model)
    
    # Prepare data if weights provided
    if (!is.null(gam_weights)) {
        df$gam_weights <- gam_weights
        df$group_numeric <- as.numeric(df$group)
    }
    
    # Fit null model: entropy ~ group + s(q)
    fit_null <- .fit_gamm_ar1_single(
        entropy ~ group + s(q, bs="tp", k=k_q_marginal),
        df, family_gam, gam_weights
    )
    
    # Fit alternative model: entropy ~ group + s(q, by=group)
    fit_alt <- .fit_gamm_ar1_single(
        entropy ~ group + s(q, bs="tp", k=k_q_interaction, by=group),
        df, family_gam, gam_weights
    )
    
    return(list(fit_null = fit_null, fit_alt = fit_alt, use_ar1 = TRUE))
}

# ===============================================================================
# HELPER FUNCTION: Fit single GAMM without correlation
# ===============================================================================
# Fits GAMM with random intercept only (no AR(1))
.fit_gamm_nocorr_single <- function(formula, df, family_gam, gam_weights) {
    # Helper to fit GAMM without correlation structure
    # Used when AR(1) convergence fails
    
    if (!is.null(gam_weights)) {
        fit <- try(
            mgcv::gamm(formula = formula,
                      random = list(subject = ~1),
                      family = family_gam,
                      weights = gam_weights,
                      data = df),
            silent = TRUE
        )
    } else {
        fit <- try(
            mgcv::gamm(formula = formula,
                      random = list(subject = ~1),
                      family = family_gam,
                      data = df),
            silent = TRUE
        )
    }
    
    return(fit)
}

# ===============================================================================
# HELPER FUNCTION: Fit GAMM without correlation (fallback 2)
# ===============================================================================
# PRIORITY 2: Falls back from AR(1) GAMM to simple GAMM
.fit_gamm_fallback <- function(df, family_gam, k_q_marginal, k_q_interaction, gam_weights) {
    # Try GAMM without correlation structure (fallback from AR(1))
    
    fit_null <- .fit_gamm_nocorr_single(
        entropy ~ group + s(q, bs="tp", k=k_q_marginal),
        df, family_gam, gam_weights
    )
    
    fit_alt <- .fit_gamm_nocorr_single(
        entropy ~ group + s(q, bs="tp", k=k_q_interaction, by=group),
        df, family_gam, gam_weights
    )
    
    return(list(fit_null = fit_null, fit_alt = fit_alt))
}

# ===============================================================================
# HELPER FUNCTION: Fit standard GAM (fallback 3)
# ===============================================================================
# PRIORITY 3: Falls back from GAMM to standard GAM

# ===============================================================================
# HELPER FUNCTION: Compare GAMM/GAM models and extract p-value
# ===============================================================================
# Performs model comparison and extracts interaction p-value
.compare_gam_models <- function(fit_null, fit_alt) {
    p_interaction <- NA_real_
    
    old_warn <- options(warn = -1)
    
    if (!is.null(fit_null$lme) && !is.null(fit_alt$lme)) {
        # GAMM comparison via LME component
        an <- try(anova(fit_null$lme, fit_alt$lme), silent = TRUE)
        if (!inherits(an, "try-error") && nrow(an) >= 2) {
            if ("p-value" %in% colnames(an)) {
                p_interaction <- an[2, "p-value"]
            } else if ("Pr(>Chisq)" %in% colnames(an)) {
                p_interaction <- an[2, "Pr(>Chisq)"]
            } else if ("Pr(>F)" %in% colnames(an)) {
                p_interaction <- an[2, "Pr(>F)"]
            }
        }
    } else {
        # Standard GAM comparison via anova.gam()
        an <- try(anova(fit_null, fit_alt, test = "Chisq"), silent = TRUE)
        if (!inherits(an, "try-error") && nrow(an) >= 2) {
            if ("p-value" %in% colnames(an)) {
                p_interaction <- an[2, "p-value"]
            } else if ("Pr(>Chi)" %in% colnames(an)) {
                p_interaction <- an[2, "Pr(>Chi)"]
            } else if ("p-value" %in% tolower(colnames(an))) {
                col_idx <- grep("p-value", tolower(colnames(an)))[1]
                if (!is.na(col_idx) && nrow(an) >= 2) {
                    p_interaction <- an[2, col_idx]
                }
            }
        }
    }
    
    options(old_warn)
    return(list(p_interaction = p_interaction, anova_result = an))
}

# ===============================================================================
# HELPER FUNCTION: Fit standard GAM (unpaired design)
# ===============================================================================
# Fitting GAM without subject/paired structure
.fit_standard_gam <- function(df, family_gam, k_q_marginal, k_q_interaction, gam_weights) {
    # BUG FIX: Use smooth splines s() instead of poly() for actual GAM fitting
    # Adaptive spline basis with thin-plate (tp) for flexible curve fitting
    
    if (!is.null(gam_weights)) {
        # Fitting with weights
        df$gam_weights <- gam_weights
        fit_null <- try(
            mgcv::gam(entropy ~ group + s(q, bs="tp", k=k_q_marginal), 
                     family = family_gam,
                     weights = gam_weights,
                     data = df), 
            silent = TRUE
        )
        # Use group-specific smooth for interaction testing
        fit_alt <- try(
            mgcv::gam(entropy ~ group + s(q, bs="tp", k=k_q_interaction, by=group),
                     family = family_gam,
                     weights = gam_weights,
                     data = df),
            silent = TRUE
        )
    } else {
        fit_null <- try(
            mgcv::gam(entropy ~ group + s(q, bs="tp", k=k_q_marginal), 
                     family = family_gam,
                     data = df), 
            silent = TRUE
        )
        # Use group-specific smooth for interaction testing
        fit_alt <- try(
            mgcv::gam(entropy ~ group + s(q, bs="tp", k=k_q_interaction, by=group),
                     family = family_gam,
                     data = df),
            silent = TRUE
        )
    }
    
    return(list(fit_null = fit_null, fit_alt = fit_alt))
}

# ===============================================================================
# HELPER FUNCTION: Fit standard GAM fallback (last resort after GAMM fails)
# ===============================================================================
# Standard GAM used as fallback when GAMM fitting fails
.fit_gam_fallback <- function(df, family_gam, k_q_marginal, k_q_interaction, gam_weights) {
    # PRIORITY 3: Fall back from GAMM to standard GAM (independence assumption)
    # This is used when both AR(1) GAMM and simple GAMM fail
    # Same implementation as .fit_standard_gam
    
    .fit_standard_gam(df, family_gam, k_q_marginal, k_q_interaction, gam_weights)
}

# ===============================================================================
# HELPER FUNCTION: Extract effect size from model summary
# ===============================================================================
# Extracts dev.expl or r.sq depending on model type
.extract_effect_size <- function(gam_summary, is_gamm) {
    # For standard GAM: use dev.expl (deviance explained)
    # For GAMM: dev.expl may be NA due to random effects, use r.sq instead
    effect_size <- NA_real_
    
    if (!is.null(gam_summary$dev.expl) && length(gam_summary$dev.expl) > 0 && is.finite(gam_summary$dev.expl)) {
        effect_size <- as.numeric(gam_summary$dev.expl)[1]
    } else if (!is_gamm && !is.null(gam_summary$r.sq) && length(gam_summary$r.sq) > 0 && is.finite(gam_summary$r.sq)) {
        effect_size <- as.numeric(gam_summary$r.sq)[1]
    } else if (is_gamm && !is.null(gam_summary$r.sq) && length(gam_summary$r.sq) > 0 && is.finite(gam_summary$r.sq)) {
        effect_size <- as.numeric(gam_summary$r.sq)[1]
    }
    
    if (!is.finite(effect_size)) {
        effect_size <- NA_real_
    }
    
    return(effect_size)
}

# ===============================================================================
# HELPER FUNCTION: Extract test statistic from anova results
# ===============================================================================
# Extracts F, L.Ratio, or Chisq from anova results
.extract_test_statistic <- function(anova_result) {
    # Extract F-statistic or likelihood ratio if anova results are available
    test_statistic <- NA_real_
    
    if (!is.null(anova_result) && nrow(anova_result) >= 2 && !inherits(anova_result, "try-error")) {
        tryCatch({
            if ("L.Ratio" %in% colnames(anova_result)) {
                test_statistic <- as.numeric(anova_result[2, "L.Ratio"])[1]
            } else if ("F" %in% colnames(anova_result)) {
                test_statistic <- as.numeric(anova_result[2, "F"])[1]
            } else if ("Chisq" %in% colnames(anova_result)) {
                test_statistic <- as.numeric(anova_result[2, "Chisq"])[1]
            }
            if (!is.finite(test_statistic)) {
                test_statistic <- NA_real_
            }
        }, error = function(e) { NULL })
    }
    
    return(test_statistic)
}

# ===============================================================================
# HELPER FUNCTION: Extract GAM/GAMM statistics
# ===============================================================================
# Extracts effect size, test statistic, and residual df from fitted model
.extract_gam_statistics <- function(fit_alt, anova_result) {
    test_statistic <- NA_real_
    effect_size <- NA_real_
    df_residual <- NA_real_
    model_converged <- !inherits(fit_alt, "try-error")
    
    if (!inherits(fit_alt, "try-error") && !is.null(fit_alt)) {
        # Extract effect size (deviance explained / R-squared equivalent)
        # Handle both GAMM (list with $gam component) and GAM (gam object directly)
        is_gamm <- is.list(fit_alt) && !is.null(fit_alt$gam)
        gam_obj <- if (is_gamm) fit_alt$gam else fit_alt
        
        # Extract summary and compute effect size
        if (!is.null(gam_obj)) {
            tryCatch({
                gam_summary <- summary(gam_obj)
                if (!is.null(gam_summary)) {
                    effect_size <- .extract_effect_size(gam_summary, is_gamm)
                    
                    # Residual df
                    if (!is.null(gam_summary$residual.df)) {
                        df_residual <- as.numeric(gam_summary$residual.df)[1]
                        if (!is.finite(df_residual)) {
                            df_residual <- NA_real_
                        }
                    }
                }
            }, error = function(e) { NULL })
        }
        
        # Extract test statistic from anova results
        test_statistic <- .extract_test_statistic(anova_result)
    }
    
    return(list(test_statistic = test_statistic,
                effect_size = effect_size,
                df_residual = df_residual,
                model_converged = model_converged))
}

# ===============================================================================
# HELPER FUNCTION: Compute slope difference between groups
# ===============================================================================
# Extracts slope_diff from GAM by computing predicted slopes for each group
.compute_slope_diff <- function(fit_alt, df, q_vals, subject) {
    slope_diff <- NA_real_
    
    if (!inherits(fit_alt, "try-error") && !is.null(fit_alt)) {
        tryCatch({
            # Get the GAM/GAMM object (may be nested in list for GAMM)
            gam_obj <- if (is.list(fit_alt) && !is.null(fit_alt$gam)) fit_alt$gam else fit_alt
            
            # Create prediction grid at min and max q for each group
            q_range <- range(df$q, na.rm = TRUE)
            unique_groups <- unique(na.omit(as.character(df$group)))
            
            if (length(unique_groups) == 2 && is.finite(q_range[1]) && is.finite(q_range[2])) {
                pred_slopes <- numeric(2)
                for (g_idx in seq_along(unique_groups)) {
                    pred_grid <- data.frame(
                        q = c(q_range[1], q_range[2]),
                        group = factor(rep(unique_groups[g_idx], 2), levels = levels(df$group))
                    )
                    
                    # Add subject if needed for GAMM
                    if (!is.null(subject) && "subject" %in% colnames(df)) {
                        pred_grid$subject <- df$subject[1]  # Use first subject as reference
                    }
                    
                    preds <- tryCatch(
                        predict(gam_obj, newdata = pred_grid, type = "response", se.fit = FALSE),
                        error = function(e) NULL
                    )
                    
                    if (!is.null(preds) && length(preds) == 2 && all(is.finite(preds))) {
                        pred_slopes[g_idx] <- (preds[2] - preds[1]) / (q_range[2] - q_range[1])
                    }
                }
                
                # Compute slope_diff if we got both slopes
                if (all(is.finite(pred_slopes))) {
                    slope_diff <- pred_slopes[2] - pred_slopes[1]
                }
            }
        }, error = function(e) { NULL })
    }
    
    return(slope_diff)
}

# ===============================================================================
# HELPER FUNCTION: Compile final results and metadata
# ===============================================================================
# Creates result data frame with all statistics and metadata
.compile_gam_results <- function(g, bc_result, test_statistic, effect_size, df_residual,
                                        model_converged, slope_diff, fit_alt, df, bounded_result,
                                        use_arima, subject) {
    # Return result with bias correction information
    result <- data.frame(
        gene = g, 
        p_interaction = bc_result$p_value,
        p_raw = bc_result$p_raw,  # Always include for reference
        n_observations = bc_result$n_observations,  # Total observations
        n_subjects = bc_result$n_subjects,  # Independent sample units
        n_effective = bc_result$n_effective,  # Effective sample size with AR(1) correlation
        rho_ar1 = bc_result$rho_estimate,  # AR(1) rho estimate
        stringsAsFactors = FALSE
    )
    
    # Explicitly add effect size, test statistic, and df columns
    result$test_statistic <- test_statistic
    result$effect_size <- effect_size
    result$df_residual <- df_residual
    result$model_converged <- model_converged
    result$slope_diff <- slope_diff
    
    # Add bias correction metadata if applied
    if (bc_result$bias_correction_applied) {
        result$bias_correction_applied <- TRUE
        result$correction_method <- bc_result$correction_method
        result$adjustment_factor <- bc_result$adjustment_factor
        result$rho_method <- if(bc_result$rho_data_driven) "data-driven" else "default"
    }
    
    # Add ARIMA transformation flag
    result$arima_transformation <- use_arima
    
    # Add bounded support handling flag and family selection information
    result$bounded_support_model <- bounded_result$use_gamma
    result$family_used <- if (bounded_result$use_gamma) "Gamma" else "Gaussian"
    result$heteroscedasticity_detected <- bounded_result$family_info$heteroscedastic
    result$variance_ratio_q <- bounded_result$family_info$var_ratio_q
    
    # Add fit method tag
    result$fit_method <- ifelse(use_arima, "mgcv::gamm_arima(1,1,0)", "mgcv::gamm")
    
    # Add residual normality testing results
    shapiro_result <- .test_residual_normality(
        model = fit_alt,
        model_type = if (!is.null(subject)) "gamm" else "gam",
        verbose = FALSE
    )
    
    if (!is.na(shapiro_result$shapiro_p_value)) {
        result$shapiro_p_value <- shapiro_result$shapiro_p_value
        result$residuals_normal <- shapiro_result$residuals_normal
        result$n_residuals_tested <- shapiro_result$n_residuals
    } else {
        result$shapiro_p_value <- NA_real_
        result$residuals_normal <- NA
        result$n_residuals_tested <- NA_integer_
    }
    
    # Add Phase 1 bootstrap CI weighting tracking
    result$ci_weighted <- !is.null(df$weight)
    
    # VALIDATION: Ensure p_interaction is always present
    if (!"p_interaction" %in% colnames(result)) {
        stop(sprintf("[.gam_interaction] CRITICAL: p_interaction missing from result for gene %s",
                     g))
    }
    
    return(result)
}

# ===============================================================================
# MAIN GAM INTERACTION FUNCTION
# ===============================================================================
# GAM interaction helper - enhanced with regularization and bias correction support
.gam_interaction <- function(df, q_vals, g, min_obs = 10, subject = NULL,
                                   regularization = c("pca", "gamsel", "spline"),
                                   bias_correction = TRUE, adaptive_knots = TRUE, weights = NULL) {
    # GAMM with ARIMA(1,1,0) covariance structure for q-dependent entropy measurements
    # Paper S171 (Zimmerman & Harville, 1991): Validates generalized covariance structures.
    # Papers S168-S170: Theoretical foundation for time series correlation patterns.
    # TEST L.1.6: Confirms first differences DeltaH_q follow AR(1) pattern via ARIMA(1,1,0).
    #
    # Adaptive knot selection (ENHANCED - March 2026):
    # - automatic_knots = TRUE: per-gene adaptivity based on entropy curve complexity
    # - Measures curve roughness via coefficient of variation of slopes
    # - Allocates more basis functions (k) to complex curves, fewer to simple curves
    # - Improves model fit efficiency and reduces overfitting/underfitting trade-off
    #
    regularization <- match.arg(regularization)
    
    # Setup and validate data
    df <- .setup_gam_data(df)
    
    # ===============================================================================
    # BOUNDED SUPPORT HANDLING (CRITICAL FIX - March 2026)
    # ===============================================================================
    # MUST be done BEFORE ARIMA differencing to:
    # 1. Detect bounds on ORIGINAL entropy (not differenced)
    # 2. Stabilize entropy before any transformations
    # 3. Preserve ARIMA structure for differenced data
    #
    bounded_result <- .handle_bounded_support(df, q_vals, group_vec = df$group, verbose = FALSE)
    
    # Select appropriate family (gaussian for GAMM paired designs)
    family_result <- .select_gam_family(bounded_result, subject)
    family_gam <- family_result$family_gam
    if (family_result$use_bounded_family && !is.null(bounded_result$stabilized_df)) {
        df <- bounded_result$stabilized_df
    }
    
    # ===============================================================================
    # HETEROSCEDASTICITY DETECTION AND VARIANCE WEIGHTING (FIXED - March 2026)
    # ===============================================================================
    # CRITICAL FIX: Detect heteroscedasticity on ORIGINAL entropy BEFORE ARIMA differencing
    hetero_result <- .detect_heteroscedasticity(df, q_vals, df$group)
    gam_weights_original <- .prepare_gam_weights(df, q_vals, weights, hetero_result, subject)
    
    # Handle ARIMA and weight updates
    arima_weights_result <- .handle_arima_and_weights(df, q_vals, subject, gam_weights_original)
    df <- arima_weights_result$df
    use_arima <- arima_weights_result$use_arima
    gam_weights <- arima_weights_result$gam_weights
    n_samples <- arima_weights_result$n_samples
    
    # Compute adaptive knots based on sample size and complexity
    knot_result <- .compute_adaptive_knots(df, q_vals, adaptive_knots)
    k_q <- knot_result$k_q
    uq_len <- knot_result$uq_len
    
    # Apply regularization for variable selection if requested (not "pca")
    reg_result <- NULL
    if (regularization != "pca") {
        reg_result <- .gam_regularization(entropy_vals = df$entropy, 
                                                 q_vals = q_vals,
                                                 group_vec = df$group,
                                                 regularization = regularization)
    }
    
    # Fit model based on presence of subject (paired design)
    # Initialize model comparison results
    p_interaction <- NA_real_
    anova_result <- NULL
    fit_null <- NULL
    fit_alt <- NULL
    
    if (!is.null(subject)) {
        # ===================================================================
        # PAIRED DESIGN: GAMM with subject random intercept (ARIMA-ready)
        # ===================================================================
        # If ARIMA differencing was not applied, add subject factor to df
        # (ARIMA result already has subject as a column)
        if (!use_arima) {
            df$subject <- factor(subject)
        } else {
            # After ARIMA differencing, subject is already in df, ensure it's a factor
            df$subject <- factor(df$subject)
        }
        
        # CRITICAL FIX: Ensure data is sorted by q|subject for corAR1 correlation structure
        # (nlme::corAR1 assumes observations are ordered by the within-group ordering variable)
        # This is especially important after adding new columns like gam_weights or group_numeric
        df <- df[order(df$subject, df$q), ]
        rownames(df) <- NULL
        
        # FIX: Create observation sequence within each subject for corAR1 ordering
        # corAR1 requires unique values within each group; use sequence index
        # This respects the q-ordering (data is sorted by q within subject) while providing unique IDs
        df$obs_seq <- unlist(lapply(rle(as.numeric(df$subject))$lengths, seq_len))
        
        # Check if we have at least 2 subjects
        if (length(unique(na.omit(df$subject))) < 2) {
            return(NULL)
        }
        
        # OPTIMIZATION (March 2026): Adaptive knot selection for smooth terms
        # Issue #4 & #8: Replace hardcoded/aggressive k parameters with data-driven selection
        # FIX (April 2026): mgcv thin-plate splines s(x, bs="tp") have minimum basis dimension k=3
        # mgcv automatically increases k if specified value is too low, causing warning
        # Previous code set min_k_adaptive=2 for small samples, triggering mgcv auto-increase
        min_k_adaptive <- 3L  # Thin-plate spline minimum basis dimension
        # Marginal smooth for q: Use conservative knots to avoid overfitting
        # Formula: at least 3 knots, but not more than available unique q-values - 1
        k_q_marginal <- as.integer(max(min_k_adaptive, min(k_q, max(min_k_adaptive, nrow(df) / 15))))
        # Interaction smooth: Use fewer knots since it must accommodate both q and group
        # Also enforce minimum of 3 to prevent mgcv auto-increase warnings
        k_q_interaction <- as.integer(max(3L, min(k_q / 2, 4L)))  # Min 3 for tp splines, cap at 4
        
        # PRIORITY 1: Try GAMM with AR(1) correlation
        fit_result <- .fit_gamm_ar1(df, family_gam, k_q_marginal, k_q_interaction, gam_weights)
        fit_null <- fit_result$fit_null
        fit_alt <- fit_result$fit_alt
        
        # PRIORITY 2: If AR(1) convergence failed, try GAMM without correlation structure
        if (inherits(fit_null, "try-error") || inherits(fit_alt, "try-error")) {
            fit_result <- .fit_gamm_fallback(df, family_gam, 
                                                    k_q_marginal, k_q_interaction, gam_weights)
            fit_null <- fit_result$fit_null
            fit_alt <- fit_result$fit_alt
        }
        
        # PRIORITY 3: If GAMM fails entirely, fall back to standard GAM
        if (inherits(fit_null, "try-error") || inherits(fit_alt, "try-error")) {
            fit_result <- .fit_gam_fallback(df, family_gam, 
                                                   k_q_marginal, k_q_interaction, gam_weights)
            fit_null <- fit_result$fit_null
            fit_alt <- fit_result$fit_alt
        }
        
        if (inherits(fit_null, "try-error") && inherits(fit_alt, "try-error")) {
            return(NULL)
        }
        
        # Compare models and extract p-value
        compare_result <- .compare_gam_models(fit_null, fit_alt)
        p_interaction <- compare_result$p_interaction
        anova_result <- compare_result$anova_result
        
    } else {
        # ===================================================================
        # UNPAIRED DESIGN: Standard GAM with independence assumption
        # ===================================================================
        # Determine k values for standard GAM (non-paired design)
        # Adaptive minimum k based on sample size: mgcv needs ~15-20 points per basis function
        min_k_adaptive <- if (nrow(df) < 25) 2L else if (nrow(df) < 50) 3L else 4L
        k_q_marginal <- as.integer(max(min_k_adaptive, min(k_q, max(min_k_adaptive, nrow(df) / 25))))
        k_q_interaction <- as.integer(max(2L, min(k_q, max(2L, nrow(df) / 20))))
        
        fit_result <- .fit_standard_gam(df, family_gam, k_q_marginal, k_q_interaction, gam_weights)
        fit_null <- fit_result$fit_null
        fit_alt <- fit_result$fit_alt
        
        # Allow model fitting errors to proceed with NA values
        # This ensures genes with convergence issues still appear in output
        if (inherits(fit_null, "try-error") && inherits(fit_alt, "try-error")) {
            return(NULL)
        }
        
        # Suppress NaN warnings from anova.gam F-test with small samples or edge cases
        old_warn <- options(warn = -1)
        anova_result <- try(mgcv::anova.gam(fit_null, fit_alt, test = "F"), silent = TRUE)
        options(old_warn)
        
        # If anova fails, proceed with NA p-value instead of returning NULL
        if (inherits(anova_result, "try-error")) {
            p_interaction <- NA_real_
        } else {
            p_interaction <- NA_real_
            if (nrow(anova_result) >= 2) {
                if ("Pr(F)" %in% colnames(anova_result)) {
                    p_interaction <- anova_result[2, "Pr(F)"]
                } else if ("Pr(>F)" %in% colnames(anova_result)) {
                    p_interaction <- anova_result[2, "Pr(>F)"]
                } else if ("p-value" %in% colnames(anova_result)) {
                    p_interaction <- anova_result[2, "p-value"]
                }
            }
        }
    }
    
    # Apply GAM-specific bias correction for small samples (C071)
    # Account for ARIMA(1,1,0) correlation structure in Tsallis entropy measurements
    # n_observations = total data points; n_subjects = independent observational units
    n_subjects_bc <- if (!is.null(subject)) length(unique(na.omit(subject))) else NULL
    bc_result <- .gam_bias_correct(p_interaction, n_observations = n_samples,
                                          n_subjects = n_subjects_bc,
                                          ar1_correlation = TRUE,
                                          bias_correction = bias_correction,
                                          entropy_data = df$entropy,
                                          subject_data = df$subject)
    
    # Extract test statistic and effect size from models
    stats_result <- .extract_gam_statistics(fit_alt, anova_result)
    
    # Compute slope difference between groups
    slope_diff <- .compute_slope_diff(fit_alt, df, q_vals, subject)
    
    # Compile final results with all metadata
    result <- .compile_gam_results(g, bc_result, stats_result$test_statistic,
                                         stats_result$effect_size, stats_result$df_residual,
                                         stats_result$model_converged, slope_diff, fit_alt,
                                         df, bounded_result, use_arima, subject)
    
    return(result)
}
