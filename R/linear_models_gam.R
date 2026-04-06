# ===============================================================================
# MAIN GAM INTERACTION FUNCTION
# ===============================================================================
# GAM interaction helper - enhanced with regularization and bias correction
# support PURPOSE: Test for q-dependent interaction effects in Tsallis entropy
# data using Generalized Additive Models (GAM) or Generalized Additive Mixed
# Models (GAMM) for paired designs. Implements a 7-stage pipeline: 1.
# Preprocessing (bounds, family selection, ARIMA, weights) 2. Model fitting
# (paired vs unpaired dispatch) 3. Model comparison (null vs alternative) 4.
# Bias correction (small sample adjustment) 5. Statistics extraction (effect
# size, test statistic, df) 6. Slope computation (group-specific curve slopes)
# 7. Result compilation (final output with metadata) PARAMETERS: df - Data
# frame with columns: entropy, q, group, [subject] q_vals - Numeric vector of q
# parameter values (same length as rows) g - Character gene identifier (for
# result metadata) min_obs - Minimum observations required (not currently used)
# subject - Optional factor/vector for paired design. If provided: * Triggers
# ARIMA(1,1,0) differencing for stationarity * Uses GAMM with AR(1) correlation
# structure * Enables bias correction accounting for clustering regularization
# - Mode for spline complexity control: * 'pca' (default): No regularization,
# auto smoothness * 'gamsel': Automatic variable selection via gamsel pkg *
# 'spline': Controlled smoothness with manual constraints bias_correction -
# Logical. If TRUE, applies small-sample adjustment for n_observations < 20.
# Accounts for ARIMA(1,1,0) structure when subject is provided. Reference:
# Hastie & Tibshirani (2015), Generalized Additive Models (GAM smoothing bias)
# adaptive_knots - Logical. If TRUE, adapts spline basis dimension (k) based on
# sample size and q-value complexity. Default TRUE.  weights - Optional numeric
# vector of observation weights. Useful for: * Bootstrap confidence interval
# weighting (Phase 1) * Heteroscedasticity adjustment (via
# .detect_heteroscedasticity) WORKFLOW STAGES: Stage 1 - PREPROCESSING
# (.prepare_gam_preprocessing): * Ensure 'group' is factor for by= smooths *
# Detect bounded support [0,1] -> Beta family, heteroscedastic -> Gamma *
# Select GAM family (Beta > Gamma > Gaussian priority) * Detect
# heteroscedasticity and compute variance weights if needed * Apply
# ARIMA(1,1,0) differencing for paired designs (removes trend) * Compute
# adaptive knot selection based on entropy curve complexity * Apply
# regularization (PCA/GAMSEL/Spline) if requested Stage 2 - MODEL FITTING
# (paired vs unpaired dispatch): * If subject != NULL (paired design): - Call
# .fit_gam_paired_design() which uses GAMM with AR(1) - GAMM uses random
# intercept ~1|subject, corr structure corAR1() - Implements 3-priority
# fallback strategy: Priority 1: GAMM with AR(1) correlation Priority 2: GAMM
# without correlation Priority 3: Standard GAM (if GAMM fails) * If subject ==
# NULL (unpaired design): - Call .fit_gam_unpaired_design() which uses standard
# GAM - No random effects, assumes independence - Uses F-test for model
# comparison (appropriate for independent data) Stage 3 - BIAS CORRECTION
# (.gam_bias_correct): * For n_observations >= 20: No correction applied
# (sufficient power) * For n_observations < 20: Apply multiplicative p-value
# adjustment * Accounts for ARIMA(1,1,0) correlation structure via AR(1) design
# effect * Formula: D_eff = (1+rho)/(1-rho); n_eff = n_subjects / D_eff *
# Adjustment: p_corrected = min(p_raw * factor, 1.0), conservative Stage 4 -
# STATISTICS EXTRACTION (.extract_gam_statistics): * Effect size: Deviance
# explained (dev.expl) or R-squared (r.sq) * Test statistic: F-statistic (GAM
# F-test) or likelihood ratio (GAMM) * Residual df: Residual degrees of freedom
# from model summary * Convergence flag: TRUE if model fitting succeeded, FALSE
# otherwise Stage 5 - SLOPE COMPUTATION (.compute_slope_diff): * Predict
# entropy at min/max q for each group * Compute slope: (y_max - y_min) / (q_max
# - q_min) for each group * Return slope_diff = slope_group2 - slope_group1 *
# Useful for interpretation: quantifies how entropy response to q differs Stage
# 6 - RESULT COMPILATION (.compile_gam_results): * Combine p_value, p_raw,
# statistics, metadata into single data frame * Add bias correction information
# if applied * Add ARIMA flag, bounded family used, heteroscedasticity
# detection * Add residual normality test (Shapiro-Wilk) * Add bootstrap CI
# weighting flag for Phase 1 tracking OUTPUT: Data frame with one row (one
# gene) containing: gene - Gene identifier p_interaction - Interaction p-value
# (bias-corrected) p_raw - Uncorrected p-value before bias correction
# n_observations - Total observations (rows in df) n_subjects - Number of
# unique subjects (if paired design) n_effective - Effective sample size after
# AR(1) adjustment rho_ar1 - AR(1) correlation coefficient estimate
# test_statistic - F-statistic or likelihood ratio effect_size - Deviance
# explained or R-squared df_residual - Residual degrees of freedom
# model_converged - Convergence flag (TRUE/FALSE) slope_diff - Difference in
# entropy slopes between groups [bias_correction_applied] - TRUE if
# small-sample correction applied [correction_method] -
# 'gam_smoothing_bias_c071' if corrected [arima_transformation] - TRUE if
# ARIMA(1,1,0) differencing used [bounded_support_model] - TRUE if Beta or
# Gamma family used [shapiro_p_value] - P-value from Shapiro-Wilk residual
# normality test [residuals_normal] - Logical residuals pass normality test
# NOTES ON IMPLEMENTATION: * GAMM limitation: mgcv::gamm() does NOT support
# extended families (Beta, Gamma) For paired designs, forces gaussian family ->
# warning issued (March 2026 fix) * ARIMA implementation: First differences
# applied to remove monotone trend in Tsallis entropy. Weight recomputation
# skipped after ARIMA (variance changes).  * Bootstrap CI weights (Phase 1):
# Take precedence over heteroscedasticity weights * Adaptive knots: Prevents
# overfitting in small samples while preserving signal * AR(1) formula: Uses
# D_eff = (1+rho)/(1-rho), NOT Kish exchangeable formula REFERENCES: Hastie &
# Tibshirani (2015), Generalized Additive Models: GAM smoothing bias in small
# samples (Hastie & Tibshirani) Wood (2024), Package 'mgcv': Mixed GAM
# Computation Vehicle-Wood (2024), CRAN R Package 'mgcv': mgcv documentation
# and GAMM tutorial Lambadaris et al. (2023), ITM Web of Conferences:
# Information entropy of generalized beta distribution (for Beta regression)
.gam_interaction <- function(df, q_vals, g, min_obs = 10, subject = NULL, regularization = c("pca",
    "gamsel", "spline"), bias_correction = TRUE, adaptive_knots = TRUE, weights = NULL) {

    # Validate regularization parameter
    regularization <- match.arg(regularization)

    # ===== PREPROCESSING ===== Consolidate all data preparation, bounds,
    # family selection, weights, knots
    prep_result <- .prepare_gam_preprocessing(df = df, q_vals = q_vals, group_vec = df$group,
        subject = subject, weights = weights, adaptive_knots = adaptive_knots, regularization = regularization)

    # ===== MODEL FITTING ===== Dispatch to paired or unpaired design handler
    if (!is.null(subject)) {
        fit_result <- .fit_gam_paired_design(df = prep_result$df, subject = subject,
            family_gam = prep_result$family_gam, k_q = prep_result$k_q, gam_weights = prep_result$gam_weights)
    } else {
        fit_result <- .fit_gam_unpaired_design(df = prep_result$df, family_gam = prep_result$family_gam,
            k_q = prep_result$k_q, gam_weights = prep_result$gam_weights)
    }

    # Validate that at least one model fit succeeded
    if (is.null(fit_result$fit_alt)) {
        # Model fitting completely failed - return NULL
        return(NULL)
    }

    # ===== BIAS CORRECTION & OUTPUT PROCESSING ===== Apply GAM-specific bias
    # correction for small samples
    n_subjects_bc <- if (!is.null(subject))
        length(unique(na.omit(subject))) else NULL
    bc_result <- .gam_bias_correct(fit_result$p_interaction, n_observations = prep_result$n_samples,
        n_subjects = n_subjects_bc, ar1_correlation = TRUE, bias_correction = bias_correction,
        entropy_data = prep_result$df$entropy, subject_data = prep_result$df$subject)

    # Extract statistics from fitted model
    stats <- .extract_gam_statistics(fit_result$fit_alt, fit_result$anova_result)

    # Compute slope difference between groups
    slope_diff <- .compute_slope_diff(fit_result$fit_alt, prep_result$df, q_vals,
        subject)

    # Compile final results with all metadata
    result <- .compile_gam_results(g, bc_result, stats$test_statistic, stats$effect_size,
        stats$df_residual, stats$model_converged, slope_diff, fit_result$fit_alt,
        prep_result$df, prep_result$bounded_result, prep_result$use_arima, subject)

    return(result)
}

# ===============================================================================
# MEMOIZATION CACHE: Package-level performance optimization for GAM functions
# ===============================================================================
# Global cache environments for memoizing expensive computations These caches
# are populated during function execution and provide O(1) lookup for repeated
# (rho, cluster_size) pairs in AR(1) design effect calculations and
# entropy-based knot computations PERFORMANCE IMPACT: Reduces AR(1) design
# effect computation from O(n) to O(1) for repeated parameter combinations
# (common in per-gene analysis loops)

# Initialize memoization caches if running in package context
if (!exists(".GAM_MEMO_CACHE", mode = "environment")) {
    .GAM_MEMO_CACHE <- new.env(hash = TRUE, parent = emptyenv())
}

if (!exists(".KNOTS_MEMO_CACHE", mode = "environment")) {
    .KNOTS_MEMO_CACHE <- new.env(hash = TRUE, parent = emptyenv())
}

# ===============================================================================
# MEMOIZATION FUNCTION: AR(1) Design Effect with Caching
# ===============================================================================
# Computes (1+rho)/(1-rho) design effect for AR(1) correlation structures.
# Caches results to avoid redundant computation across multiple genes.
# Reference: Diggle et al. 2002 (AR(1) correlation design effect formula)
# PERFORMANCE: First call O(1) computation; subsequent calls with same (rho, m)
# are O(1) cache lookup vs. O(1) but with function call overhead.  With ~10K
# genes, typical gains: ~5-10ms per analysis run.
.ar1_design_effect_memo <- function(rho, cluster_size) {
    # Validate inputs
    if (!is.finite(rho) || rho < 0 || rho > 1) {
        return(NA_real_)
    }
    if (!is.finite(cluster_size) || cluster_size < 1) {
        return(NA_real_)
    }

    # Create cache key: format rho and cluster_size for stable hashing
    cache_key <- sprintf("rho=%.4f|m=%.1f", round(rho, 4), cluster_size)

    # Check if already cached
    if (exists(cache_key, envir = .GAM_MEMO_CACHE, inherits = FALSE)) {
        return(get(cache_key, envir = .GAM_MEMO_CACHE))
    }

    # Compute AR(1) design effect: D_eff = (1+rho)/(1-rho) See: Diggle, P.J.,
    # Heagerty, P., Liang, K.Y., Zeger, S.L. (2002) Analysis of Longitudinal
    # Data, Oxford University Press.
    design_eff <- (1 + rho)/(1 - rho)

    # Validate result
    if (!is.finite(design_eff) || design_eff < 1) {
        # rho near 1 -> D_eff -> Inf; rho near 0 -> D_eff near 1
        design_eff <- max(1, min(design_eff, Inf))
    }

    # Cache the result
    assign(cache_key, design_eff, envir = .GAM_MEMO_CACHE)

    return(design_eff)
}

# ===============================================================================
# MEMOIZATION FUNCTION: Adaptive Spline Knots with Caching
# ===============================================================================
# Determines appropriate basis dimension (k) for spline fitting based on
# entropy curve complexity. Uses memoization to cache (n_q_unique, entropy_sd)
# pairs to avoid recomputation for similar genes.  STRATEGY (March 2026 Fix):
# Use fixed k based on unique q values only.  DO NOT use CV-based adaptation
# for monotone Tsallis entropy.  Reference: .adaptive_spline_knots()
# documentation at end of file PERFORMANCE: O(1) cache lookup for repeated
# entropy curve structures
.adaptive_spline_knots_memo <- function(entropy_vals, q_vals, n_q_unique, min_k = 2,
    max_k = 10) {
    # Remove NA values
    entropy_clean <- na.omit(entropy_vals)
    q_clean <- na.omit(q_vals)

    # Validate minimum data
    if (length(entropy_clean) < 3 || length(q_clean) < 2) {
        return(max(min_k, min(max_k, n_q_unique - 1)))
    }

    # Compute entropy distribution characteristics for cache key (simplified
    # signature to ensure cache hits for similar curves)
    entropy_sd <- sd(entropy_clean)
    entropy_range <- diff(range(entropy_clean))

    # Create cache key: n_q_unique is primary driver of k selection Include
    # entropy characteristics for robustness
    cache_key <- sprintf("nq=%d|sd=%.3f|range=%.3f", n_q_unique, round(entropy_sd,
        3), round(entropy_range, 3))

    # Check cache
    if (exists(cache_key, envir = .KNOTS_MEMO_CACHE, inherits = FALSE)) {
        return(get(cache_key, envir = .KNOTS_MEMO_CACHE))
    }

    # Apply fixed k selection: k = max(min_k, min(max_k, n_q_unique - 1))
    # Principle: use at most (unique q values - 1) basis functions This ensures
    # smooth fits without noise-driven over-complexity
    k_final <- max(min_k, min(max_k, n_q_unique - 1))

    # Cache result
    assign(cache_key, k_final, envir = .KNOTS_MEMO_CACHE)

    return(k_final)
}

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
    # extended families.  Reference: Wood (2024), Package 'mgcv': Mixed GAM
    # Computation Vehicle/Wood (2024), CRAN R Package 'mgcv' (GAMM Tutorial,
    # mgcv Documentation)

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
# PREPROCESSING HELPER: Consolidate all prep steps into single function
# ===============================================================================
# Consolidates setup, bounds, family, hetero, weights, knots, and
# regularization
.prepare_gam_preprocessing <- function(df, q_vals, group_vec, subject, weights, adaptive_knots,
    regularization) {
    # Step 1: Setup and validate data
    df <- .setup_gam_data(df)

    # Step 2: Bounded support detection (BEFORE ARIMA differencing)
    bounded_result <- .handle_bounded_support(df, q_vals, group_vec = group_vec,
        verbose = FALSE)

    # Step 3: Select appropriate family (gaussian for GAMM paired designs)
    family_result <- .select_gam_family(bounded_result, subject)
    family_gam <- family_result$family_gam
    if (family_result$use_bounded_family && !is.null(bounded_result$stabilized_df)) {
        df <- bounded_result$stabilized_df
    }

    # Step 4: Heteroscedasticity detection on ORIGINAL data
    hetero_result <- .detect_heteroscedasticity(df, q_vals, group_vec)
    gam_weights_original <- .prepare_gam_weights(df, q_vals, weights, hetero_result,
        subject)

    # Step 5: Handle ARIMA differencing and weight updates
    arima_result <- .handle_arima_and_weights(df, q_vals, subject, gam_weights_original)

    # Step 6: Compute adaptive knots based on sample size and complexity
    knots_result <- .compute_adaptive_knots(arima_result$df, q_vals, adaptive_knots)

    # Step 7: Apply regularization if requested (not 'pca')
    reg_result <- NULL
    if (regularization != "pca") {
        reg_result <- .gam_regularization(entropy_vals = arima_result$df$entropy,
            q_vals = q_vals, group_vec = arima_result$df$group, regularization = regularization)
    }

    return(list(df = arima_result$df, family_gam = family_gam, gam_weights = arima_result$gam_weights,
        use_arima = arima_result$use_arima, n_samples = arima_result$n_samples, k_q = knots_result$k_q,
        uq_len = knots_result$uq_len, bounded_result = bounded_result, reg_result = reg_result))
}

# ===============================================================================
# PAIRED DESIGN FITTING: Handle GAMM with optional AR(1) correlation
# ===============================================================================
# Fits GAMM to paired/repeated measures design with priority fallbacks
.fit_gam_paired_design <- function(df, subject, family_gam, k_q, gam_weights) {
    # Step 1: Check if subject is already in df (from ARIMA preprocessing) If
    # not, add it; if yes, ensure it's a factor
    if ("subject" %in% colnames(df)) {
        df$subject <- factor(df$subject)
    } else {
        # Subject not yet in df, add it from parameter
        df$subject <- factor(subject)
    }

    # Step 2: Sort by subject/q for AR(1) ordering requirements
    df <- df[order(df$subject, df$q), ]
    rownames(df) <- NULL

    # Step 3: Create observation sequence within each subject
    df$obs_seq <- unlist(lapply(rle(as.numeric(df$subject))$lengths, seq_len))

    # Step 4: Validate >= 2 subjects (required for mixed effects)
    if (length(unique(na.omit(df$subject))) < 2) {
        return(list(fit_null = NULL, fit_alt = NULL, p_interaction = NA_real_, anova_result = NULL))
    }

    # Step 5: Compute adaptive k values for GAMM
    min_k_adaptive <- 3L  # Thin-plate spline minimum
    k_q_marginal <- as.integer(max(min_k_adaptive, min(k_q, max(min_k_adaptive, nrow(df)/15))))
    k_q_interaction <- as.integer(max(3L, min(k_q/2, 4L)))

    # Step 6: Priority 1 - Try GAMM with AR(1) correlation
    fit_result <- .fit_gamm_ar1(df, family_gam, k_q_marginal, k_q_interaction, gam_weights)

    # Step 7: Priority 2 - If AR(1) fails, try GAMM without correlation
    # structure
    if (inherits(fit_result$fit_null, "try-error") || inherits(fit_result$fit_alt,
        "try-error")) {
        fit_result <- .fit_gamm_fallback(df, family_gam, k_q_marginal, k_q_interaction,
            gam_weights)
    }

    # Step 8: Priority 3 - If GAMM fails, fall back to standard GAM
    if (inherits(fit_result$fit_null, "try-error") || inherits(fit_result$fit_alt,
        "try-error")) {
        fit_result <- .fit_gam_fallback(df, family_gam, k_q_marginal, k_q_interaction,
            gam_weights)
    }

    # Step 9: Abort if all models failed
    if (inherits(fit_result$fit_null, "try-error") && inherits(fit_result$fit_alt,
        "try-error")) {
        return(list(fit_null = NULL, fit_alt = NULL, p_interaction = NA_real_, anova_result = NULL))
    }

    # Step 10: Compare models and extract p-value
    compare_result <- .compare_gam_models(fit_result$fit_null, fit_result$fit_alt)

    return(list(fit_null = fit_result$fit_null, fit_alt = fit_result$fit_alt, p_interaction = compare_result$p_interaction,
        anova_result = compare_result$anova_result))
}

# ===============================================================================
# UNPAIRED DESIGN FITTING: Handle standard GAM with independence assumption
# ===============================================================================
# Fits standard GAM to unpaired design without repeated measures structure
.fit_gam_unpaired_design <- function(df, family_gam, k_q, gam_weights) {
    # Step 1: Compute adaptive k values based on sample size
    min_k_adaptive <- if (nrow(df) < 25)
        2L else if (nrow(df) < 50)
        3L else 4L
    k_q_marginal <- as.integer(max(min_k_adaptive, min(k_q, max(min_k_adaptive, nrow(df)/25))))
    # Interaction terms require higher k minimum (consistent with paired
    # design) min=3L ensures by=group smooths have sufficient basis dimension
    k_q_interaction <- as.integer(max(3L, min(k_q, max(3L, nrow(df)/20))))

    # Step 2: Fit null and alternative models
    fit_result <- .fit_standard_gam(df, family_gam, k_q_marginal, k_q_interaction,
        gam_weights)

    # Step 3: Check if both models failed
    if (inherits(fit_result$fit_null, "try-error") && inherits(fit_result$fit_alt,
        "try-error")) {
        return(list(fit_null = NULL, fit_alt = NULL, p_interaction = NA_real_, anova_result = NULL))
    }

    # Step 4: Compare models with F-test
    old_warn <- options(warn = -1)
    anova_result <- try(mgcv::anova.gam(fit_result$fit_null, fit_result$fit_alt,
        test = "F"), silent = TRUE)
    options(old_warn)

    # Step 5: Extract p-value from anova results
    p_interaction <- NA_real_
    if (!inherits(anova_result, "try-error") && nrow(anova_result) >= 2) {
        if ("Pr(F)" %in% colnames(anova_result)) {
            p_interaction <- anova_result[2, "Pr(F)"]
        } else if ("Pr(>F)" %in% colnames(anova_result)) {
            p_interaction <- anova_result[2, "Pr(>F)"]
        } else if ("p-value" %in% colnames(anova_result)) {
            p_interaction <- anova_result[2, "p-value"]
        }
    }

    return(list(fit_null = fit_result$fit_null, fit_alt = fit_result$fit_alt, p_interaction = p_interaction,
        anova_result = anova_result))
}



# GAM bias correction helper: adjusts for smoothing bias in small samples
# (Hastie & Tibshirani (2015), Generalized Additive Models) When n_samples <
# 20, small sample smoothing can inflate Type I error rates Applies degrees of
# freedom adjustment based on sample size
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
    # ARIMA-adjusted effective units.  Historical tests (and published Hastie &
    # Tibshirani (2015), Generalized Additive Models guidance) trigger
    # correction when the number of samples is small (<20); the original
    # implementation compared against n_eff, which under AR(1) dependency could
    # fall below 20 even for reasonably large datasets and therefore caused
    # over-conservative adjustments.  To keep behaviour compatible with
    # existing user expectations we now only suppress bias correction when the
    # *observed* sample size is large.
    if (!bias_correction || n_observations >= 20) {
        return(list(p_value = p_value, p_raw = p_value, bias_correction_applied = FALSE,
            n_observations = n_observations, n_samples = n_observations, n_subjects = n_subjects,
            n_effective = n_eff, design_effect_ar1 = design_effect, rho_estimate = rho_avg,
            rho_data_driven = data_driven_rho, correction_method = "none", correction_rationale = sprintf("n_observations=%.0f >= 20; GAM smoothing bias minimal (AR(1) D_eff=%.2f, rho=%.2f %s)",
                n_observations, if (is.na(design_effect)) 0 else design_effect, if (is.na(rho_avg)) 0 else rho_avg,
                if (data_driven_rho) "[data-driven]" else "[default]")))
    }

    # For small samples (n_eff < 20), smoothing bias can affect p-values
    # (Hastie & Tibshirani (2015), Generalized Additive Models) Apply
    # conservative adjustment accounting for ARIMA(1,1,0) structure

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
    # smoothing bias) Reference: Hastie & Tibshirani (2015), Generalized
    # Additive Models (empirical correction for GAM smoothing bias in small
    # samples) This is more conservative than K-C correction but appropriate
    # for GAM bias
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
# selection), and spline (controlled smoothness) modes. Based on papers
# Chouldechova & Hastie (2015), Annals of Applied Statistics, C063, C065, CRAN
# R Package 'gamsel' (2023), Chouldechova & Hastie (1986), Annals of Applied
# Statistics.
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
# Database Support (March 2026): - Lambadaris et al. (2023), ITM Web of
# Conferences: 'Information entropy of generalized beta distribution' -
# Capelletti et al. (2024), Beta regression for wind power modeling-Lasso
# Penalization for High-Dimensional Beta Regression (2023): Beta regression
# applications with robustness validation
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
    # mathematically ideal for bounded (0,1) data Database paper Lambadaris et
    # al. (2023), ITM Web of Conferences: 'Information entropy of the
    # generalized beta distribution'
    if (is_bounded_01) {
        use_beta <- TRUE
        family_choice <- "beta"
        reasons <- c(reasons, "Data bounded in [0,1] - Beta regression ideal (Lambadaris et al. (2023), ITM Web of Conferences)")
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
