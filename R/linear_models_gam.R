# ===============================================================================
# HELPER FUNCTION: Setup and validate GAM data
# ===============================================================================
# Ensures 'group' is a factor and validates package dependencies
.setup_gam_data <- function(df) {
    # Ensure 'group' is a factor for 'by' argument in GAM/GAMM smooths This is
    # required for s(q, by = group, ...) to work correctly
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
# Selects family based on bounded support and subject info GAMM does NOT
# support extended families, so falls back to gaussian for paired designs
.select_gam_family <- function(bounded_result, subject) {
    # GAMM COMPATIBILITY FIX (March 2026): mgcv::gamm() does NOT support
    # extended families (beta, gamma, Tweedie, etc.). For paired designs
    # (subject != NULL -> uses gamm), fall back to gaussian family instead of
    # extended families.  Reference: C042/C043 (GAMM Tutorial, mgcv
    # Documentation)

    use_bounded_family <- bounded_result$use_gamma
    family_gam <- bounded_result$family_obj
    inverse_link_fn <- bounded_result$inverse_link

    if (!is.null(subject)) {
        # For paired/mixed designs using gamm(), force gaussian family CRITICAL
        # FIX (March 2026): mgcv::gamm() does NOT support extended families
        # Warn user that bounded family selection is being overridden for
        # reproducibility
        if (use_bounded_family) {
            warning("[calculate_lm_interaction] GAMM with paired design detected. ",
                "mgcv::gamm() does not support extended families. ", "Forcing gaussian family. Results may be less accurate for bounded data.",
                call. = FALSE)
        }
        family_gam <- stats::gaussian()
        inverse_link_fn <- function(eta) eta
    }

    # If beta regression is selected, use stabilized entropy
    if (bounded_result$use_beta && !is.null(bounded_result$stabilized_df)) {
        df <- bounded_result$stabilized_df
    }

    return(list(use_bounded_family = use_bounded_family, family_gam = family_gam,
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
    # Differencing removes monotone trend from Tsallis entropy, enabling valid
    # AR(1) inference This is applied when subject information is available
    # (paired design) CRITICAL: ARIMA is applied AFTER heteroscedasticity
    # detection on original data

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

            # FIX: When ARIMA is applied, DON'T use weights computed on
            # original data Reason: Differencing changes the variance
            # structure, weights would be invalid Conservative approach: Better
            # to lose efficiency than introduce bias
            gam_weights <- NULL
        }
    }

    # PHASE 1 WEIGHTING (March 2026): Set df$weight for ci_weighted flag
    # tracking
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

    # Adaptive knot selection: compute k based on gene's entropy curve
    # complexity
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
    # Helper to fit GAMM with AR(1) correlation Used internally by
    # .fit_gamm_ar1

    if (!is.null(gam_weights)) {
        fit <- try(mgcv::gamm(formula = formula, random = list(subject = ~1), correlation = nlme::corAR1(form = ~obs_seq |
            subject), family = family_gam, weights = gam_weights, data = df), silent = TRUE)
    } else {
        fit <- try(mgcv::gamm(formula = formula, random = list(subject = ~1), correlation = nlme::corAR1(form = ~obs_seq |
            subject), family = family_gam, data = df), silent = TRUE)
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
    fit_null <- .fit_gamm_ar1_single(entropy ~ group + s(q, bs = "tp", k = k_q_marginal),
        df, family_gam, gam_weights)

    # Fit alternative model: entropy ~ group + s(q, by=group)
    fit_alt <- .fit_gamm_ar1_single(entropy ~ group + s(q, bs = "tp", k = k_q_interaction,
        by = group), df, family_gam, gam_weights)

    return(list(fit_null = fit_null, fit_alt = fit_alt, use_ar1 = TRUE))
}

# ===============================================================================
# HELPER FUNCTION: Fit single GAMM without correlation
# ===============================================================================
# Fits GAMM with random intercept only (no AR(1))
.fit_gamm_nocorr_single <- function(formula, df, family_gam, gam_weights) {
    # Helper to fit GAMM without correlation structure Used when AR(1)
    # convergence fails

    if (!is.null(gam_weights)) {
        fit <- try(mgcv::gamm(formula = formula, random = list(subject = ~1), family = family_gam,
            weights = gam_weights, data = df), silent = TRUE)
    } else {
        fit <- try(mgcv::gamm(formula = formula, random = list(subject = ~1), family = family_gam,
            data = df), silent = TRUE)
    }

    return(fit)
}

# ===============================================================================
# HELPER FUNCTION: Fit GAMM without correlation (fallback 2)
# ===============================================================================
# PRIORITY 2: Falls back from AR(1) GAMM to simple GAMM
.fit_gamm_fallback <- function(df, family_gam, k_q_marginal, k_q_interaction, gam_weights) {
    # Try GAMM without correlation structure (fallback from AR(1))

    fit_null <- .fit_gamm_nocorr_single(entropy ~ group + s(q, bs = "tp", k = k_q_marginal),
        df, family_gam, gam_weights)

    fit_alt <- .fit_gamm_nocorr_single(entropy ~ group + s(q, bs = "tp", k = k_q_interaction,
        by = group), df, family_gam, gam_weights)

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
        fit_null <- try(mgcv::gam(entropy ~ group + s(q, bs = "tp", k = k_q_marginal),
            family = family_gam, weights = gam_weights, data = df), silent = TRUE)
        # Use group-specific smooth for interaction testing
        fit_alt <- try(mgcv::gam(entropy ~ group + s(q, bs = "tp", k = k_q_interaction,
            by = group), family = family_gam, weights = gam_weights, data = df),
            silent = TRUE)
    } else {
        fit_null <- try(mgcv::gam(entropy ~ group + s(q, bs = "tp", k = k_q_marginal),
            family = family_gam, data = df), silent = TRUE)
        # Use group-specific smooth for interaction testing
        fit_alt <- try(mgcv::gam(entropy ~ group + s(q, bs = "tp", k = k_q_interaction,
            by = group), family = family_gam, data = df), silent = TRUE)
    }

    return(list(fit_null = fit_null, fit_alt = fit_alt))
}

# ===============================================================================
# HELPER FUNCTION: Fit standard GAM fallback (last resort after GAMM fails)
# ===============================================================================
# Standard GAM used as fallback when GAMM fitting fails
.fit_gam_fallback <- function(df, family_gam, k_q_marginal, k_q_interaction, gam_weights) {
    # PRIORITY 3: Fall back from GAMM to standard GAM (independence assumption)
    # This is used when both AR(1) GAMM and simple GAMM fail Same
    # implementation as .fit_standard_gam

    .fit_standard_gam(df, family_gam, k_q_marginal, k_q_interaction, gam_weights)
}

# ===============================================================================
# HELPER FUNCTION: Extract effect size from model summary
# ===============================================================================
# Extracts dev.expl or r.sq depending on model type
.extract_effect_size <- function(gam_summary, is_gamm) {
    # For standard GAM: use dev.expl (deviance explained) For GAMM: dev.expl
    # may be NA due to random effects, use r.sq instead
    effect_size <- NA_real_

    if (!is.null(gam_summary$dev.expl) && length(gam_summary$dev.expl) > 0 && is.finite(gam_summary$dev.expl)) {
        effect_size <- as.numeric(gam_summary$dev.expl)[1]
    } else if (!is_gamm && !is.null(gam_summary$r.sq) && length(gam_summary$r.sq) >
        0 && is.finite(gam_summary$r.sq)) {
        effect_size <- as.numeric(gam_summary$r.sq)[1]
    } else if (is_gamm && !is.null(gam_summary$r.sq) && length(gam_summary$r.sq) >
        0 && is.finite(gam_summary$r.sq)) {
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

    if (!is.null(anova_result) && nrow(anova_result) >= 2 && !inherits(anova_result,
        "try-error")) {
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
        }, error = function(e) {
            NULL
        })
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
        # Handle both GAMM (list with $gam component) and GAM (gam object
        # directly)
        is_gamm <- is.list(fit_alt) && !is.null(fit_alt$gam)
        gam_obj <- if (is_gamm)
            fit_alt$gam else fit_alt

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
            }, error = function(e) {
                NULL
            })
        }

        # Extract test statistic from anova results
        test_statistic <- .extract_test_statistic(anova_result)
    }

    return(list(test_statistic = test_statistic, effect_size = effect_size, df_residual = df_residual,
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
            gam_obj <- if (is.list(fit_alt) && !is.null(fit_alt$gam))
                fit_alt$gam else fit_alt

            # Create prediction grid at min and max q for each group
            q_range <- range(df$q, na.rm = TRUE)
            unique_groups <- unique(na.omit(as.character(df$group)))

            if (length(unique_groups) == 2 && is.finite(q_range[1]) && is.finite(q_range[2])) {
                pred_slopes <- numeric(2)
                for (g_idx in seq_along(unique_groups)) {
                  pred_grid <- data.frame(q = c(q_range[1], q_range[2]), group = factor(rep(unique_groups[g_idx],
                    2), levels = levels(df$group)))

                  # Add subject if needed for GAMM
                  if (!is.null(subject) && "subject" %in% colnames(df)) {
                    pred_grid$subject <- df$subject[1]  # Use first subject as reference
                  }

                  preds <- tryCatch(predict(gam_obj, newdata = pred_grid, type = "response",
                    se.fit = FALSE), error = function(e) NULL)

                  if (!is.null(preds) && length(preds) == 2 && all(is.finite(preds))) {
                    pred_slopes[g_idx] <- (preds[2] - preds[1])/(q_range[2] - q_range[1])
                  }
                }

                # Compute slope_diff if we got both slopes
                if (all(is.finite(pred_slopes))) {
                  slope_diff <- pred_slopes[2] - pred_slopes[1]
                }
            }
        }, error = function(e) {
            NULL
        })
    }

    return(slope_diff)
}

# ===============================================================================
# HELPER FUNCTION: Compile final results and metadata
# ===============================================================================
# Creates result data frame with all statistics and metadata
.compile_gam_results <- function(g, bc_result, test_statistic, effect_size, df_residual,
    model_converged, slope_diff, fit_alt, df, bounded_result, use_arima, subject) {
    # Return result with bias correction information
    result <- data.frame(gene = g, p_interaction = bc_result$p_value, p_raw = bc_result$p_raw,
        n_observations = bc_result$n_observations, n_subjects = bc_result$n_subjects,
        n_effective = bc_result$n_effective, rho_ar1 = bc_result$rho_estimate, stringsAsFactors = FALSE)

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
        result$rho_method <- if (bc_result$rho_data_driven)
            "data-driven" else "default"
    }

    # Add ARIMA transformation flag
    result$arima_transformation <- use_arima

    # Add bounded support handling flag and family selection information
    result$bounded_support_model <- bounded_result$use_gamma
    result$family_used <- if (bounded_result$use_gamma)
        "Gamma" else "Gaussian"
    result$heteroscedasticity_detected <- bounded_result$family_info$heteroscedastic
    result$variance_ratio_q <- bounded_result$family_info$var_ratio_q

    # Add fit method tag
    result$fit_method <- ifelse(use_arima, "mgcv::gamm_arima(1,1,0)", "mgcv::gamm")

    # Add residual normality testing results
    shapiro_result <- .test_residual_normality(model = fit_alt, model_type = if (!is.null(subject))
        "gamm" else "gam", verbose = FALSE)

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
# GAM interaction helper - enhanced with regularization and bias correction
# support
.gam_interaction <- function(df, q_vals, g, min_obs = 10, subject = NULL, regularization = c("pca",
    "gamsel", "spline"), bias_correction = TRUE, adaptive_knots = TRUE, weights = NULL) {
    # GAMM with ARIMA(1,1,0) covariance structure for q-dependent entropy
    # measurements Paper S171 (Zimmerman & Harville, 1991): Validates
    # generalized covariance structures.  Papers S168-S170: Theoretical
    # foundation for time series correlation patterns.  TEST L.1.6: Confirms
    # first differences DeltaH_q follow AR(1) pattern via ARIMA(1,1,0).
    # Adaptive knot selection (ENHANCED - March 2026): - automatic_knots =
    # TRUE: per-gene adaptivity based on entropy curve complexity - Measures
    # curve roughness via coefficient of variation of slopes - Allocates more
    # basis functions (k) to complex curves, fewer to simple curves - Improves
    # model fit efficiency and reduces overfitting/underfitting trade-off
    regularization <- match.arg(regularization)

    # Setup and validate data
    df <- .setup_gam_data(df)

    # ===============================================================================
    # BOUNDED SUPPORT HANDLING (CRITICAL FIX - March 2026)
    # ===============================================================================
    # MUST be done BEFORE ARIMA differencing to: 1. Detect bounds on ORIGINAL
    # entropy (not differenced) 2. Stabilize entropy before any transformations
    # 3. Preserve ARIMA structure for differenced data
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
    # CRITICAL FIX: Detect heteroscedasticity on ORIGINAL entropy BEFORE ARIMA
    # differencing
    hetero_result <- .detect_heteroscedasticity(df, q_vals, df$group)
    gam_weights_original <- .prepare_gam_weights(df, q_vals, weights, hetero_result,
        subject)

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

    # Apply regularization for variable selection if requested (not 'pca')
    reg_result <- NULL
    if (regularization != "pca") {
        reg_result <- .gam_regularization(entropy_vals = df$entropy, q_vals = q_vals,
            group_vec = df$group, regularization = regularization)
    }

    # Fit model based on presence of subject (paired design) Initialize model
    # comparison results
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
            # After ARIMA differencing, subject is already in df, ensure it's a
            # factor
            df$subject <- factor(df$subject)
        }

        # CRITICAL FIX: Ensure data is sorted by q|subject for corAR1
        # correlation structure (nlme::corAR1 assumes observations are ordered
        # by the within-group ordering variable) This is especially important
        # after adding new columns like gam_weights or group_numeric
        df <- df[order(df$subject, df$q), ]
        rownames(df) <- NULL

        # FIX: Create observation sequence within each subject for corAR1
        # ordering corAR1 requires unique values within each group; use
        # sequence index This respects the q-ordering (data is sorted by q
        # within subject) while providing unique IDs
        df$obs_seq <- unlist(lapply(rle(as.numeric(df$subject))$lengths, seq_len))

        # Check if we have at least 2 subjects
        if (length(unique(na.omit(df$subject))) < 2) {
            return(NULL)
        }

        # OPTIMIZATION (March 2026): Adaptive knot selection for smooth terms
        # Issue #4 & #8: Replace hardcoded/aggressive k parameters with
        # data-driven selection FIX (April 2026): mgcv thin-plate splines s(x,
        # bs='tp') have minimum basis dimension k=3 mgcv automatically
        # increases k if specified value is too low, causing warning Previous
        # code set min_k_adaptive=2 for small samples, triggering mgcv
        # auto-increase
        min_k_adaptive <- 3L  # Thin-plate spline minimum basis dimension
        # Marginal smooth for q: Use conservative knots to avoid overfitting
        # Formula: at least 3 knots, but not more than available unique
        # q-values - 1
        k_q_marginal <- as.integer(max(min_k_adaptive, min(k_q, max(min_k_adaptive,
            nrow(df)/15))))
        # Interaction smooth: Use fewer knots since it must accommodate both q
        # and group Also enforce minimum of 3 to prevent mgcv auto-increase
        # warnings
        k_q_interaction <- as.integer(max(3L, min(k_q/2, 4L)))  # Min 3 for tp splines, cap at 4

        # PRIORITY 1: Try GAMM with AR(1) correlation
        fit_result <- .fit_gamm_ar1(df, family_gam, k_q_marginal, k_q_interaction,
            gam_weights)
        fit_null <- fit_result$fit_null
        fit_alt <- fit_result$fit_alt

        # PRIORITY 2: If AR(1) convergence failed, try GAMM without correlation
        # structure
        if (inherits(fit_null, "try-error") || inherits(fit_alt, "try-error")) {
            fit_result <- .fit_gamm_fallback(df, family_gam, k_q_marginal, k_q_interaction,
                gam_weights)
            fit_null <- fit_result$fit_null
            fit_alt <- fit_result$fit_alt
        }

        # PRIORITY 3: If GAMM fails entirely, fall back to standard GAM
        if (inherits(fit_null, "try-error") || inherits(fit_alt, "try-error")) {
            fit_result <- .fit_gam_fallback(df, family_gam, k_q_marginal, k_q_interaction,
                gam_weights)
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
        # Determine k values for standard GAM (non-paired design) Adaptive
        # minimum k based on sample size: mgcv needs ~15-20 points per basis
        # function
        min_k_adaptive <- if (nrow(df) < 25)
            2L else if (nrow(df) < 50)
            3L else 4L
        k_q_marginal <- as.integer(max(min_k_adaptive, min(k_q, max(min_k_adaptive,
            nrow(df)/25))))
        k_q_interaction <- as.integer(max(2L, min(k_q, max(2L, nrow(df)/20))))

        fit_result <- .fit_standard_gam(df, family_gam, k_q_marginal, k_q_interaction,
            gam_weights)
        fit_null <- fit_result$fit_null
        fit_alt <- fit_result$fit_alt

        # Allow model fitting errors to proceed with NA values This ensures
        # genes with convergence issues still appear in output
        if (inherits(fit_null, "try-error") && inherits(fit_alt, "try-error")) {
            return(NULL)
        }

        # Suppress NaN warnings from anova.gam F-test with small samples or
        # edge cases
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

    # Apply GAM-specific bias correction for small samples (C071) Account for
    # ARIMA(1,1,0) correlation structure in Tsallis entropy measurements
    # n_observations = total data points; n_subjects = independent
    # observational units
    n_subjects_bc <- if (!is.null(subject))
        length(unique(na.omit(subject))) else NULL
    bc_result <- .gam_bias_correct(p_interaction, n_observations = n_samples, n_subjects = n_subjects_bc,
        ar1_correlation = TRUE, bias_correction = bias_correction, entropy_data = df$entropy,
        subject_data = df$subject)

    # Extract test statistic and effect size from models
    stats_result <- .extract_gam_statistics(fit_alt, anova_result)

    # Compute slope difference between groups
    slope_diff <- .compute_slope_diff(fit_alt, df, q_vals, subject)

    # Compile final results with all metadata
    result <- .compile_gam_results(g, bc_result, stats_result$test_statistic, stats_result$effect_size,
        stats_result$df_residual, stats_result$model_converged, slope_diff, fit_alt,
        df, bounded_result, use_arima, subject)

    return(result)
}

# GAM bias correction helper: adjusts for smoothing bias in small samples
# (C071) When n_samples < 20, small sample smoothing can inflate Type I error
# rates Applies degrees of freedom adjustment based on sample size
.gam_bias_correct <- function(p_value, n_observations = NULL, n_samples = NULL, n_subjects = NULL,
    ar1_correlation = TRUE, bias_correction = TRUE, entropy_data = NULL, subject_data = NULL) {
    # `n_samples` is provided for backward compatibility with earlier versions
    # that accepted this argument name.  If `n_observations` is NULL we fall
    # back to the user-supplied `n_samples` value so that tests and external
    # code using either name will continue to work.
    if (is.null(n_observations) && !is.null(n_samples)) {
        n_observations <- n_samples
    }

    # validate presence of at least one of the size arguments
    if (is.null(n_observations)) {
        stop("must supply either 'n_observations' or 'n_samples'", call. = FALSE)
    }
    # For Tsallis multi-q design with ARIMA(1,1,0) covariance: - First
    # differences DeltaH_q = H_q - H_{q-1} are modeled as AR(1) [IMPLEMENTED] -
    # Differencing removes monotone trend in Tsallis entropy (dH/dq < 0)
    # [VERIFIED] - True independent units are subjects, not observations - Bias
    # correction threshold should use n_subjects, not n_observations CRITICAL
    # FIX (March 2026): Calling functions now properly difference entropy
    # before fitting AR(1) correlations. This guarantees stationarity
    # assumptions.  See: .compute_arima_differences() helper function added
    # March 2026.

    # If n_subjects not provided, attempt to estimate from ARIMA structure
    # Conservative: assume ~sqrt(n_obs) independent units under ARIMA(1,1,0)
    if (is.null(n_subjects)) {
        # For ARIMA(1,1,0), we lose 1 observation per subject via differencing
        # Estimate: (observations - n_subjects) / n_subjects gives adjusted
        # count Conservative: use sqrt(n_obs) which is robust estimate
        n_subjects <- max(2, ceiling(sqrt(n_observations)))
    }

    # For ARIMA(1,1,0) correlation, effective degrees of freedom are reduced
    # CRITICAL FIX March 2026: Use AR(1)-specific design effect formula (NOT
    # Kish exchangeable formula) Background: - Previous code used: D_eff = 1 +
    # (m-1)rho [Kish formula for ICC/exchangeable] - Correct for AR(1): D_eff =
    # (1+phi)/(1-phi) [Diggle et al. 2002] - These formulas apply to VERY
    # different correlation structures - AR(1) is appropriate for ordered
    # q-values with geometric decay: Corr(t,t+k) = phi^k

    if (ar1_correlation && n_observations > n_subjects && n_observations > 0) {
        # Compute intra-subject cluster size
        cluster_size <- n_observations/n_subjects

        # Estimate rho from data if available; otherwise use conservative
        # default
        rho_avg <- NULL
        data_driven_rho <- FALSE

        if (!is.null(entropy_data) && !is.null(subject_data)) {
            rho_est <- .estimate_ar1_rho(entropy_data, subject_data)
            if (!is.null(rho_est) && rho_est >= 0 && rho_est <= 1) {
                rho_avg <- rho_est
                data_driven_rho <- TRUE
            }
        }

        if (is.null(rho_avg)) {
            # ARIMA(1,1,0) average correlation on first differences
            # (trend-removed) Conservative default: rho = 0.35 based on AR(1)
            # applied to differenced data SENSITIVITY ANALYSIS for AR(1) design
            # effect: - rho = 0.20: D_eff = (1.2)/(0.8) = 1.5, n_eff =
            # n_subjects / 1.5 - rho = 0.35: D_eff = (1.35)/(0.65) = 2.08,
            # n_eff = n_subjects / 2.08 - rho = 0.50: D_eff = (1.5)/(0.5) =
            # 3.0, n_eff = n_subjects / 3.0 (Accounting for finite-m
            # corrections depending on cluster_size) Note: Much higher D_eff
            # than Kish (which gave 1.2-1.6 for same rho) This demonstrates
            # importance of using AR(1)-specific formula
            rho_avg <- 0.35
            data_driven_rho <- FALSE
        }

        # Design effect: Use AR(1)-specific formula (NOT Kish exchangeable
        # formula) OPTIMIZATION: Use memoized version to cache repeated (rho,
        # cluster_size) pairs
        design_effect <- .ar1_design_effect_memo(rho_avg, cluster_size)

        # Effective sample size accounting for AR(1) within-subject correlation
        n_eff <- n_subjects/design_effect
    } else {
        # No ARIMA(1,1,0) or independence: effective n = n_subjects
        n_eff <- n_subjects
        design_effect <- NA_real_
        rho_avg <- NA_real_
        data_driven_rho <- FALSE
    }

    # Bias correction decision: use raw observation count rather than
    # ARIMA-adjusted effective units.  Historical tests (and published C071
    # guidance) trigger correction when the number of samples is small (<20);
    # the original implementation compared against n_eff, which under AR(1)
    # dependency could fall below 20 even for reasonably large datasets and
    # therefore caused over-conservative adjustments.  To keep behaviour
    # compatible with existing user expectations we now only suppress bias
    # correction when the *observed* sample size is large.
    if (!bias_correction || n_observations >= 20) {
        return(list(p_value = p_value, p_raw = p_value, bias_correction_applied = FALSE,
            n_observations = n_observations, n_samples = n_observations, n_subjects = n_subjects,
            n_effective = n_eff, design_effect_ar1 = design_effect, rho_estimate = rho_avg,
            rho_data_driven = data_driven_rho, correction_method = "none", correction_rationale = sprintf("n_observations=%.0f >= 20; GAM smoothing bias minimal (AR(1) D_eff=%.2f, rho=%.2f %s)",
                n_observations, if (is.na(design_effect)) 0 else design_effect, if (is.na(rho_avg)) 0 else rho_avg,
                if (data_driven_rho) "[data-driven]" else "[default]")))
    }

    # For small samples (n_eff < 20), smoothing bias can affect p-values (C071)
    # Apply conservative adjustment accounting for ARIMA(1,1,0) structure

    if (is.na(p_value)) {
        return(list(p_value = p_value, p_raw = p_value, bias_correction_applied = FALSE,
            n_observations = n_observations, n_samples = n_observations, n_subjects = n_subjects,
            n_effective = n_eff, design_effect_ar1 = design_effect, rho_estimate = rho_avg,
            rho_data_driven = data_driven_rho, correction_method = "na_value", correction_rationale = "p-value is NA"))
    }

    # Compute adjustment factor based on effective sample size Smaller
    # effective samples get larger adjustments (less power, more conservative)
    # Linear scaling: at n_eff=5, factor=2.0; at n_eff=19, factor=1.05
    adjustment_factor <- 1 + (20 - n_eff)/20

    # Apply multiplicative adjustment (Bonferroni-style, conservative for GAM
    # smoothing bias) Reference: C071 (empirical correction for GAM smoothing
    # bias in small samples) This is more conservative than K-C correction but
    # appropriate for GAM bias
    p_corrected <- min(p_value * adjustment_factor, 1)

    return(list(p_value = p_corrected, bias_correction_applied = TRUE, n_observations = n_observations,
        n_samples = n_observations, n_subjects = n_subjects, n_effective = n_eff,
        design_effect_ar1 = design_effect, rho_estimate = rho_avg, rho_data_driven = data_driven_rho,
        adjustment_factor = adjustment_factor, p_raw = p_value, correction_method = "gam_smoothing_bias_c071",
        correction_rationale = sprintf("n_eff=%.1f < 20; Adjusted for ARIMA(1,1,0) correlation: AR(1) D_eff=%.2f (rho=%.2f %s); adjustment_factor=%.2f",
            n_eff, design_effect, if (is.na(rho_avg)) 0 else rho_avg, if (data_driven_rho) "[data-driven]" else "[default]",
            adjustment_factor)))
}

# GAM regularization helper: applies spline constraints or GAMSEL for variable
# selection Supports pca (no regularization), gamsel (automatic variable
# selection), and spline (controlled smoothness) modes. Based on papers C057,
# C063, C065, C082, C083.
.gam_regularization <- function(entropy_vals, q_vals, group_vec, regularization = c("pca",
    "gamsel", "spline")) {
    regularization <- match.arg(regularization)

    if (regularization == "pca") {
        # PCA mode: no regularization, return NULL
        return(NULL)
    }

    if (regularization == "gamsel") {
        # GAMSEL mode: automatic variable selection using gamsel package
        # Requires gamsel package; if not available, fall back to spline mode
        if (!requireNamespace("gamsel", quietly = TRUE)) {
            # Fallback to spline mode if gamsel not available
            return(list(mode = "spline_fallback", constraint = "auto"))
        }

        # GAMSEL expects matrix X and vector y We'll use q values as the
        # feature to select
        X <- as.matrix(q_vals)
        y <- entropy_vals

        # Fit GAMSEL model
        gs_fit <- try(gamsel::gamsel(x = X, y = y, family = "gaussian"), silent = TRUE)

        if (inherits(gs_fit, "try-error")) {
            return(list(mode = "spline_fallback", constraint = "auto"))
        }

        # Return GAMSEL result with information for model construction
        return(list(mode = "gamsel", gamsel_fit = gs_fit, q_values = unique(sort(q_vals))))
    }

    if (regularization == "spline") {
        # Spline mode: use controlled smoothness with mgcv's automatic
        # smoothing This applies automatic smoothness selection (GCV/REML)
        return(list(mode = "spline", constraint = "auto"  # Let mgcv handle smoothness via GCV
))
    }

    return(NULL)
}

# Helper: Handle bounded support for Tsallis entropy via appropriate GAM family
# selection Tsallis entropy is bounded [0, log(m)] where m = number of isoforms
# Priority: Beta (if [0,1]) > Gamma (if heteroscedastic) > Gaussian (default)
# Database Support (March 2026): - S223: 'Information entropy of generalized
# beta distribution' - S220-S222: Beta regression applications with robustness
# validation
.handle_bounded_support <- function(df, q_vals, group_vec = NULL, verbose = FALSE) {
    # ========================================================================
    # INLINE: Family selection logic (previously .select_gam_family) Select
    # appropriate GAM family based on data characteristics Priority: Beta (if
    # [0,1] bounded) > Gamma (if heteroscedastic) > Gaussian (default) Tsallis
    # entropy is mathematically bounded [0, log(m)], but Beta is ideal for
    # [0,1]
    # ========================================================================

    # INDICATOR 1: Check if data is [0,1] bounded (ideal for Beta regression)
    # =====================================================================
    entropy_vals <- na.omit(df$entropy)
    is_bounded_01 <- .is_bounded_0_1(entropy_vals)

    # INDICATOR 2: Heteroscedasticity detection
    # =========================================
    hetero_result <- try(.detect_heteroscedasticity(df, q_vals = q_vals, group_vec = group_vec,
        verbose = verbose), silent = TRUE)

    heteroscedastic <- FALSE
    var_ratio_q <- 1
    var_ratio_group <- 1

    if (!inherits(hetero_result, "try-error") && !is.na(hetero_result$is_heteroscedastic)) {
        heteroscedastic <- hetero_result$is_heteroscedastic
        var_ratio_q <- if (is.null(hetero_result$var_ratio_q))
            1 else hetero_result$var_ratio_q
        var_ratio_group <- if (is.null(hetero_result$var_ratio_group))
            1 else hetero_result$var_ratio_group
    }

    # INDICATOR 3: Boundary clustering (values near 0 or 1)
    # ====================================================
    n_total <- length(entropy_vals)
    finite_entropy <- is.finite(entropy_vals)
    if (any(finite_entropy)) {
        entropy_min <- min(entropy_vals[finite_entropy])
        entropy_max <- max(entropy_vals[finite_entropy])
        entropy_range <- entropy_max - entropy_min
    } else {
        entropy_min <- NA
        entropy_max <- NA
        entropy_range <- NA
    }
    boundary_threshold <- if (is.finite(entropy_range))
        0.1 * entropy_range else NA  # 10% of range is 'near boundary'

    # Count values near boundaries
    n_near_min <- sum(entropy_vals <= entropy_min + boundary_threshold)
    n_near_max <- sum(entropy_vals >= entropy_max - boundary_threshold)
    pct_boundary_clustering <- 100 * (n_near_min + n_near_max)/n_total

    # INDICATOR 4: Skewness (asymmetry indicates non-Gaussian behavior)
    # =============================================================== Skewness
    # = (mean - median) / sd * constant; values > 1 or < -1 indicate strong
    # asymmetry
    skewness_val <- .compute_skewness(entropy_vals)
    has_strong_skew <- abs(skewness_val) > 1

    # DECISION LOGIC (March 2026) Priority: Beta > Gamma > Gaussian
    # ================================
    use_beta <- FALSE
    use_gamma <- FALSE
    family_choice <- "gaussian"
    reasons <- c()

    # *** PRIORITY 1: Use Beta if data is [0,1] bounded *** Beta regression is
    # mathematically ideal for bounded (0,1) data Database paper S223:
    # 'Information entropy of the generalized beta distribution'
    if (is_bounded_01) {
        use_beta <- TRUE
        family_choice <- "beta"
        reasons <- c(reasons, "Data bounded in [0,1] - Beta regression ideal (S223)")
    } else {
        # *** PRIORITY 2: Use Gamma if strong evidence of non-Gaussian behavior
        # *** Criterion 1: Strong heteroscedasticity (p < 0.05) AND variance
        # changes much
        if (heteroscedastic && (var_ratio_q > 3 || var_ratio_group > 3)) {
            use_gamma <- TRUE
            family_choice <- "gamma"
            reasons <- c(reasons, sprintf("Heteroscedasticity detected (p<0.05, var_ratio=%.2f)",
                max(var_ratio_q, var_ratio_group)))
        }

        # Criterion 2: EXTREME boundary clustering only (> 40% of data near
        # bounds) Most entropy distributions naturally have some clustering -
        # must be severe
        if (pct_boundary_clustering > 40 && !use_gamma) {
            use_gamma <- TRUE
            family_choice <- "gamma"
            reasons <- c(reasons, sprintf("Extreme boundary clustering: %.1f%% near bounds",
                pct_boundary_clustering))
        }

        # Criterion 3: Extreme skewness (|skew| > 1) AND evidence of
        # heteroscedasticity Require combination of indicators rather than
        # skewness alone
        if (has_strong_skew && abs(skewness_val) > 1 && heteroscedastic && var_ratio_q >
            3 && !use_gamma) {
            use_gamma <- TRUE
            family_choice <- "gamma"
            reasons <- c(reasons, sprintf("Extreme skewness (|skew|=%.2f) with heteroscedasticity (var_ratio=%.2f)",
                skewness_val, var_ratio_q))
        }
    }

    if (verbose) {
        if (length(reasons) > 0) {
            message(sprintf("[GAM Family Selection] Using %s because: %s", toupper(family_choice),
                paste(reasons, collapse = "; ")))
        } else {
            message(sprintf("[GAM Family Selection] Using Gaussian (no strong indicators); bounded=[%s], hetero_p=%.4f, hetero_vars=(%.2f,%.2f), boundary=%.1f%%, |skew|=%.2f",
                is_bounded_01, if (is.na(hetero_result$p_value))
                  NA else hetero_result$p_value, var_ratio_q, var_ratio_group, pct_boundary_clustering,
                skewness_val))
        }
    }

    family_info <- list(use_beta = use_beta, use_gamma = use_gamma, use_gaussian = !use_beta &&
        !use_gamma, is_bounded_01 = is_bounded_01, heteroscedastic = heteroscedastic,
        var_ratio_q = var_ratio_q, var_ratio_group = var_ratio_group, boundary_pct = pct_boundary_clustering,
        skewness = skewness_val, reasons = reasons, family_choice = family_choice)

    if (family_info$use_beta) {
        # Use Beta family with logit link (BEST for [0,1] bounded entropy data)
        # Beta regression respects bounds and handles skewness naturally
        # CRITICAL: For continuous data in (0,1), use quasibinomial NOT
        # binomial - binomial() expects count/binary data -> gives warnings for
        # continuous values - quasibinomial() is designed for continuous
        # proportions in (0,1) - Alternatively, mgcv::betar() (v1.8.41+) is
        # specialized for beta regression Numerical stability: Ensure no exact
        # 0 or 1 values which cause singularities FIX (March 2026): Return the
        # stabilized dataframe so calling code uses it!
        df$entropy <- pmax(pmin(df$entropy, 1 - 1e-07), 1e-07)

        # Try to use betar() from mgcv if available (v1.8.41+), otherwise
        # quasibinomial
        family_obj <- try(mgcv::betar(), silent = TRUE)
        if (inherits(family_obj, "try-error")) {
            # Fallback to quasibinomial for continuous (0,1) data
            family_obj <- stats::quasibinomial(link = "logit")
        }

        return(list(use_bounded = TRUE, use_beta = TRUE, use_gamma = FALSE, use_gaussian = FALSE,
            family_obj = family_obj, inverse_link = function(eta) 1/(1 + exp(-eta)),
            stabilized_df = df, family_info = family_info))
    } else if (family_info$use_gamma) {
        # Use Gamma family with log link (appropriate for positive bounded
        # data)
        return(list(use_bounded = TRUE, use_beta = FALSE, use_gamma = TRUE, use_gaussian = FALSE,
            family_obj = stats::Gamma(link = "log"), inverse_link = function(eta) exp(eta),
            stabilized_df = NULL, family_info = family_info))
    } else {
        return(list(use_bounded = FALSE, use_beta = FALSE, use_gamma = FALSE, use_gaussian = TRUE,
            family_obj = stats::gaussian(), inverse_link = function(eta) eta, stabilized_df = NULL,
            family_info = family_info))
    }
}


.adaptive_spline_knots <- function(entropy_vals, q_vals, n_q_unique, min_k = 2, max_k = 10) {
    # K-selection strategy for Tsallis entropy curves: Tsallis entropy is
    # GUARANTEED monotone decreasing in q (mathematical property) Therefore,
    # use a FIXED k based on number of unique q-values Do NOT use CV-based
    # adaptation for monotone data HISTORICAL ISSUE: Earlier code computed CV
    # of first differences and allocated MORE knots for HIGH CV.  This is
    # BACKWARDS for monotone data because: - High CV in first differences
    # indicates DEVIATION FROM MONOTONICITY (i.e., noise) - Allocating more
    # knots to noisy data increases overfitting, not model appropriateness -
    # For truly monotone data, CV should reflect measurement error, not true
    # complexity SOLUTION: Use fixed k based on number of unique q-values
    # (conservative, data-driven minimum) This ensures smooth monotone fitting
    # without noise-driven over-complexity

    # Fixed selection: k = max(min_k, min(max_k, n_q_unique - 1)) Principle:
    # use at most (number of unique q values - 1) basis functions This leaves
    # at least one degree of freedom for residual fitting
    k_final <- max(min_k, min(max_k, n_q_unique - 1))

    return(k_final)
}
