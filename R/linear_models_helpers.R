# ================================================================================
# BIAS CORRECTION STRATEGY IN TSENAT - CURRENT STATUS
# ================================================================================
#
# IMPLEMENTATION SUMMARY (Phase 8-14):
#
# [PLANNED] GEE: K-C bias correction for sandwich variance (Future implementation)
#   * Intended location: R/effect_size.R
#   * Purpose: Address sandwich variance underestimation in small samples with GEE
#   * Proposed method: Kenward-Roger/HC1 correction
#   * Literature support: 34 papers recommend K-C for GEE sandwich variance estimation
#   * Target use case: GEE models with small clusters (n < 30)
#   * Status: Not yet implemented
#
# [OK] GAM: IMPLEMENTED - Smoothing bias correction for small samples
#   * Location: R/calc_lm_helpers.R (function .gam_bias_correct, line 281)
#   * Called from: R/linear_models_gam.R (line 766)
#   * Purpose: Address Type I error inflation from smoothing splines when n_samples < 20
#   * Method: Conservative p-value adjustment (factor = 1 + (n_eff-20)/n_eff), accounts for AR(1) correlation
#   * Literature support: 10+ papers discuss smoothing bias in GAM; paper C071
#   * When applied: GAM method with n_observations < 20 and bias_correction=TRUE
#   * Test file: test-integration-linear_models.R (lines 1686-1754+) - 10+ tests passing
#
# [X] LMM: ENHANCED (Phase 14) - AR(1) correlation & sample size tracking
#   * Analysis: Papers S88, S90, S92, S93, S135, S147 confirm LMM best practices
#   * New Feature 1: AR(1) correlation structure option (.try_lmm_ar1)
#     - Implements corAR1(form = ~time_idx | subject) for ordered q-values
#     - Appropriate for entropy curves with autocorrelated differences
#     - Tries AR(1) FIRST in fallback hierarchy
#   * New Feature 2: Sample size documentation (.extract_lrt_p returns list)
#     - Returns: p_value, n_subjects, small_sample_flag (TRUE if n < 5)
#     - Results include n_subjects column for transparency
#     - Flag indicates LOW POWER region (requires caution)
#
# LITERATURE BASIS:
# Papers confirming proposed strategy:
#   * S160: Selection bias in linear mixed models (bias in parameter estimation context)
#   * S161: Bias Correction in GLMM (focus: variance component and coefficient estimation)
#   * S163-S164: Bias correction for parameter estimation, not test statistics
#   * Phase 11 validation: Database analysis confirming GEE/GAM bias corrections would be needed
#
# KEY DISTINCTION:
# * Parameter Estimation Bias: Average value of estimator differs from true parameter
#   -> Discussed in papers S160-S164 for LMM/GLMM
#   -> Could affect confidence intervals if severe
#   -> Not currently affecting HYPOTHESIS TESTING in TSENAT (Satterthwaite is sufficient)
#
# * Hypothesis Testing Bias: Type I error rate differs from nominal alpha
#   -> Proposed for GEE (sandwich variance bias) - PLANNED (ref: effect_size.R not yet created)
#   -> Implemented for GAM (smoothing bias) - VERIFIED (ref: .gam_bias_correct in calc_lm_helpers.R)
#   -> NOT needed for LMM (Satterthwaite inherently conservative) - VERIFIED
#
# ================================================================================

# Summary reporting helper
.report_fit_summary <- function(res, verbose = TRUE) {
    if (verbose && "fit_method" %in% colnames(res)) {
        total_genes <- nrow(res)
        fallback_mask <- !is.na(res$fit_method) & res$fit_method != "lmer"
        n_fallback <- sum(fallback_mask)
        
        # Report model convergence and fit quality
        message(sprintf("[calculate_lm_interaction] Analyzed %d genes | Primary fits: %d | Alternative method: %d",
                        total_genes, total_genes - n_fallback, n_fallback))
        
        if (n_fallback > 0) {
            tab <- table(res$fit_method[fallback_mask])
            tab_str <- paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse = ", ")
            message(sprintf("[calculate_lm_interaction]   Methods: %s", tab_str))
        }
        
        # Report numerical/convergence issues
        if ("singular" %in% colnames(res)) {
            n_sing <- sum(as.logical(res$singular), na.rm = TRUE)
            if (n_sing > 0) {
                message(sprintf("[calculate_lm_interaction]   Singular fits (collinear effects): %d genes", n_sing))
            }
        }
        
        # Report effect size range (q-interaction magnitude)
        if ("f_statistic" %in% colnames(res)) {
            f_vals <- res$f_statistic[!is.na(res$f_statistic)]
            if (length(f_vals) > 0) {
                message(sprintf("[calculate_lm_interaction] q x condition interaction strength: F-stat range [%.2f, %.2f]",
                                min(f_vals), max(f_vals)))
            }
        }
        
        # Report significance summary
        if ("adj_p_interaction" %in% colnames(res)) {
            sig_p <- sum(res$adj_p_interaction < 0.05, na.rm = TRUE)
            message(sprintf("[calculate_lm_interaction] Significant results (adj.p < 0.05): %d/%d genes (%.1f%%)",
                            sig_p, total_genes, 100 * sig_p / total_genes))
        }
    }
}

# Helper: Estimate autocorrelation rho from differenced entropy data
# Used to compute design effect for bias correction
# Helper: Compute AR(1) design effect (NOT Kish formula for exchangeable ICC)
# 
# CRITICAL CORRECTION March 2026: Previous code mistakenly used Kish formula
# designed for exchangeable correlation (ICC), but TSENAT uses AR(1) correlation
# after ARIMA(1,1,0) differencing. These are fundamentally different.
# 
# For AR(1) correlation with autocorrelation coefficient phi:
#   D_eff = (1 + phi) / (1 - phi)  [for large m: m >> 1]
#   
# For moderate m (typical in multi-q designs where m = q-values per subject):
#   D_eff = (1 + phi) / (1 - phi) * [1 - phi^m] / [m - (m-1)phi^m]
#   
# This formula assumes:
#   - rho applied to DIFFERENCED entropy (ARIMA(1,1,0) applied first)
#   - m = cluster_size = observations per subject (typically q-values)
#   - phi = autocorrelation on differenced data (0 < phi < 1)
#
# Reference:
#   Diggle et al. (2002) "Analysis of Longitudinal Data" Section 4.3
#   Crowder (1995) "Generalised Estimating Equations for repeated measurements"
#   Liang & Zeger (1986) "Longitudinal data analysis using GEE"
#
# COMPARISON: Why AR(1) formula differs from Kish:
#   Kish formula 1 + (m-1)rho: Assumes exchangeable correlation (ICC)
#     - All pairs equally correlated with ICC rho
#     - Appropriate for clusters with homogeneous correlation
#   AR(1) formula (1+phi)/(1-phi): Assumes geometric correlation decay
#     - Correlation decreases as lag k increases: Corr(t, t+k) = phi^k
#     - Appropriate for ordered measurements (like q-values)
#     - Applied to differenced data (ARIMA(1,1,0) stationarity)
#

.ar1_design_effect <- function(rho, cluster_size) {
    # Compute design effect for AR(1) correlation
    # Args:
    #   rho: autocorrelation coefficient phi on differenced data (0 <= phi <= 1)
    #   cluster_size: m = observations per subject (e.g., number of q-values)
    # Returns:
    #   D_eff = design effect to adjust effective sample size as n_eff = n_subjects / D_eff
    
    if (is.null(rho) || is.na(rho) || rho <= 0 || cluster_size <= 1) {
        # No correlation or invalid input: D_eff = 1 (independence)
        return(1.0)
    }
    
    if (rho >= 1) {
        # Perfect correlation: D_eff = m (complete non-independence)
        return(as.numeric(cluster_size))
    }
    
    # AR(1) design effect: D_eff = 1 + 2*Sum_{k=1}^{m-1} (1 - k/m)*rho^k
    # This accounts for geometric correlation decay and edge effects
    # Reference: Diggle et al. (2002), Crowder (1995)
    
    # Compute sum efficiently
    summed <- 0
    for (k in seq_len(cluster_size - 1)) {
        lambda_k <- 1 - k / cluster_size  # Edge effect weight
        summed <- summed + lambda_k * (rho ^ k)
    }
    
    d_eff <- 1 + 2 * summed
    d_eff <- pmax(1.0, d_eff)  # Ensure D_eff >= 1
    
    return(d_eff)
}

.estimate_ar1_rho <- function(entropy_diff, subject_vec = NULL) {
    # Estimate first-order autocorrelation rho from differenced entropy
    # Input: entropy_diff = first-differenced entropy values DeltaH_q = H_q - H_{q-1}
    # Returns: rho estimate in [0, 1], or NULL if insufficient data
    # Issues warning if rho is very high (GAMM convergence risk)
    
    if (is.null(entropy_diff) || length(na.omit(entropy_diff)) < 3) {
        return(NULL)
    }
    
    # Remove NA values
    entropy_clean <- na.omit(entropy_diff)
    
    if (length(entropy_clean) < 3) {
        return(NULL)
    }
    
    # OPTIMIZATION (March 2026): Use stats::acf() for numerical stability
    # Previous: Manual computation (divides by n instead of n-1, less stable)
    # New: Built-in acf() for better stability and standard handling
    acf_result <- tryCatch({
        stats::acf(entropy_clean, lag.max = 1, plot = FALSE, demean = TRUE)
    }, error = function(e) NULL)
    
    if (!is.null(acf_result)) {
        rho_est <- as.numeric(acf_result$acf[2, 1, 1])
    } else {
        # Fallback to manual calculation if acf fails
        n <- length(entropy_clean)
        mean_x <- mean(entropy_clean, na.rm = TRUE)
        var_x <- sum((entropy_clean - mean_x)^2, na.rm = TRUE) / n
        
        if (var_x < 1e-10) {
            return(NULL)
        }
        
        x_t <- entropy_clean[-n]
        x_t1 <- entropy_clean[-1]
        cov_lag1 <- sum((x_t - mean_x) * (x_t1 - mean_x), na.rm = TRUE) / n
        rho_est <- cov_lag1 / var_x
    }
    
    # Ensure rho is in [0, 1] (sometimes numerical errors give slight negative values)
    rho_est <- max(0, min(1, rho_est))
    
    # OPTIMIZATION (March 2026): Add tolerance checks for edge cases
    # Issue #9: Missing tolerance checks from CODE_REVIEW_BUGS_FOUND.md
    if (rho_est > 0.95) {
        # Very high autocorrelation - warn about GAMM convergence risk
        warning(sprintf("AR(1) autocorrelation very high (rho=%.3f). GAMM may fail to converge. Consider reducing q-values or checking data for trends.", rho_est), call. = FALSE)
    }
    
    if (rho_est < 0.01) {
        # Very small autocorrelation - independence assumption near valid
        # Return NULL to suggest simpler model without AR(1)
        return(NULL)
    }
    
    return(rho_est)
}

# Helper: Knot selection for Tsallis entropy curve fitting
#
# Tsallis entropy is MATHEMATICALLY GUARANTEED to be monotone decreasing in q
# Therefore, k-selection uses a simple fixed formula based on number of unique q-values
# This ensures adequate smoothing without noise-driven over-complexity
#
# Historical note: Earlier versions attempted CV-based adaptation to detect
# curve complexity, but this was backwards for monotone data (high CV indicates noise)
# 
# Current approach: k = max(min_k, min(max_k, n_q_unique - 1))
# Principle: Use at most (# unique q values - 1) basis functions
# This provides data-driven parsimony while ensuring sufficient flexibility
#
# Reference: Wood (2006) Generalized Additive Models; enforced via fixed formula
# for mathematical monotonicity property of Tsallis entropy
#
# STATIONARITY VALIDATION FRAMEWORK FOR TSALLIS ENTROPY MODELING
# ================================================================
# 
# MATHEMATICAL JUSTIFICATION:
# Tsallis entropy H_q is monotone DECREASING in q parameter (proven in Tsallis 1988)
# This monotonicity makes the series NON-STATIONARY (systematic/deterministic trend)
# 
# AR(1) models assume stationarity (constant mean/variance around trend)
# Solution: ARIMA(1,1,0) = Apply AR(1) to FIRST DIFFERENCES DeltaH_q = H_q - H_{q-1}
# This removes the trend (differencing) allowing AR(1) to model residual correlation
#
# VALIDATION FRAMEWORK: Four complementary tests to validate assumptions
# ========================================================================
# 
# TEST 1: Monotonicity Check (Visual)
#   Purpose: Detect q-value ordering issues or data quality problems
#   Method: Count decreasing vs increasing pairs in ordered q-values
#   Expected for Tsallis: >95% pairs should be decreasing (monotone)
#   
# TEST 2: Augmented Dickey-Fuller (ADF) Test for Unit Root
#   Null Hypothesis (H0): Series has unit root (non-stationary)
#   Alternative (H1): Series is stationary
#   Expected for raw Tsallis entropy: FAIL to reject H0 (non-stationary with unit root)
#   Expected for differenced data: REJECT H0 (stationary, no unit root)
#   Reference: Dickey & Fuller (1979, 1981); MacKinnon (1996) for critical values
#   
# TEST 3: KPSS Test (Reverse of ADF)
#   Null Hypothesis (H0): Series IS stationary
#   Alternative (H1): Series is NON-stationary
#   Expected for raw Tsallis entropy: REJECT H0 (non-stationary)
#   Expected for differenced data: FAIL to reject H0 (stationary)
#   Reference: Kwiatkowski, Phillips, Schmidt & Shin (1992)
#   
# TEST 4: Integration Order Validation
#   Apply ADF/KPSS to differenced data
#   Confirms that ARIMA(1,1,0) with integration order d=1 is appropriate
#   Expected: Differenced data should be I(0) - integrated of order 0 (stationary)
#
# DECISION LOGIC FOR ARIMA(1,1,0):
#   Use ARIMA(1,1,0) if ALL conditions met:
#   [OK] Raw data is NON-monotone (>5% violations) OR fails stationarity tests
#   [OK] ADF test FAILs to reject H0 on raw data (has unit root)
#   [OK] KPSS test REJECTs H0 on raw data (non-stationary)
#   [OK] ADF test REJECTs H0 on differenced data (stationary)
#   [OK] KPSS test FAILs to reject H0 on differenced data (stationary)
#
# REFERENCES:
#   Dickey, D. A., & Fuller, W. A. (1979). Distribution of the estimators for
#     autoregressive time series with a unit root. Journal of the American
#     Statistical Association, 74(366), 427-431.
#   Kwiatkowski, D., Phillips, P. C., Schmidt, P., & Shin, Y. (1992). Testing
#     the null hypothesis of stationarity against the alternative of a unit root.
#     Journal of Econometrics, 54(1-3), 159-178.
#   MacKinnon, J. G. (1996). Numerical distribution functions for unit root
#     and cointegration tests. Journal of Applied Econometrics, 11(6), 601-618.
#   Tsallis, C. (1988). Possible generalization of Boltzmann-Gibbs statistics.
#     Journal of Statistical Physics, 52(1), 479-487.

# STATIONARITY VALIDATION FRAMEWORK
# ===================================
# Validates core mathematical assumptions before applying AR(1) models
# Reference: Null hypothesis tests for time series stationarity (Dickey-Fuller, KPSS)
# 
# For Tsallis entropy in ordered q-values:
# - Raw data: H_q is monotone decreasing -> SHOULD BE NON-STATIONARY
# - Differenced: DeltaH_q = H_q - H_{q-1} -> SHOULD BE STATIONARY
# 
# These tests VALIDATE the ARIMA(1,1,0) modeling approach

# ===============================================================================
# RESIDUAL DIAGNOSTICS: Shapiro-Wilk Normality Testing
# ===============================================================================
# 
# DATABASE EVIDENCE (March 2026):
# * B001 (2001) - Foundations of Systems Biology
# * B004 (2008) - LINEAR MODELS IN [Systems Biology]
# * C017 (2006) - Springer Handbook of Statistical Methods
#
# Purpose: Verify that residuals from GAM/LMM/GEE models satisfy normality assumption
# Method: Shapiro-Wilk test on model residuals (tests H0: residuals are normal)
# Standard Practice: Applied universally in statistical modeling literature
# Interpretation:
#   * p > 0.05: Fail to reject H0 -> Residuals appear normal [OK]
#   * p <= 0.05: Reject H0 -> Residuals show significant departure from normality ?
#
# Implementation: Extract residuals from fitted model, apply shapiro.test()

.test_residual_normality <- function(model, model_type = c("gam", "gamm", "lme", "gee"),
                                            verbose = FALSE) {
    # Args:
    #   model: fitted model object (GAM, GAMM, lme, or geeglm)
    #   model_type: character - type of model for residual extraction
    #   verbose: if TRUE, print diagnostic messages
    # Returns:
    #   List with components:
    #   - shapiro_p_value: p-value from Shapiro-Wilk test (NA if test fails)
    #   - residuals_normal: logical - TRUE if p > 0.05 (residuals appear normal)
    #   - n_residuals: number of residuals tested
    #   - test_status: character - "pass", "fail", or "error"
    #   - report: character - human-readable summary
    
    if (is.null(model) || inherits(model, "try-error")) {
        return(list(
            shapiro_p_value = NA_real_,
            residuals_normal = NA,
            n_residuals = 0,
            test_status = "error",
            report = "Model object is NULL or error class"
        ))
    }
    
    model_type <- match.arg(model_type)
    residuals_vec <- NULL
    
    # Extract residuals based on model type
    tryCatch({
        if (model_type == "gam") {
            # Standard GAM: check family to determine residual type
            # Pearson residuals available for non-gaussian families (gamma, beta, etc.)
            # Gaussian family: use deviance residuals directly
            family_name <- ifelse(!is.null(model$family), model$family$family, "gaussian")
            if (family_name == "gaussian") {
                residuals_vec <- residuals(model, type = "deviance")
            } else {
                # For non-gaussian families, use Pearson residuals
                residuals_vec <- tryCatch(
                    residuals(model, type = "pearson"),
                    error = function(e) residuals(model, type = "deviance")
                )
            }
        } else if (model_type == "gamm") {
            # GAMM: extract residuals from $gam component, check family first
            if (!is.null(model$gam)) {
                family_name <- ifelse(!is.null(model$gam$family), model$gam$family$family, "gaussian")
                if (family_name == "gaussian") {
                    residuals_vec <- residuals(model$gam, type = "deviance")
                } else {
                    residuals_vec <- tryCatch(
                        residuals(model$gam, type = "pearson"),
                        error = function(e) residuals(model$gam, type = "deviance")
                    )
                }
            } else {
                family_name <- ifelse(!is.null(model$family), model$family$family, "gaussian")
                if (family_name == "gaussian") {
                    residuals_vec <- residuals(model, type = "deviance")
                } else {
                    residuals_vec <- tryCatch(
                        residuals(model, type = "pearson"),
                        error = function(e) residuals(model, type = "deviance")
                    )
                }
            }
        } else if (model_type == "lme") {
            # nlme::lme model: use residuals() generic
            residuals_vec <- residuals(model, type = "normalized")
        } else if (model_type == "gee") {
            # geeglm: use residuals() generic
            residuals_vec <- tryCatch(
                residuals(model, type = "pearson"),
                error = function(e) residuals(model, type = "deviance")
            )
        }
    }, error = function(e) {
        if (verbose) {
            message("[.test_residual_normality] Could not extract residuals: ", e$message)
        }
    })
    
    # Check if residuals were extracted successfully
    if (is.null(residuals_vec) || length(residuals_vec) == 0) {
        return(list(
            shapiro_p_value = NA_real_,
            residuals_normal = NA,
            n_residuals = 0,
            test_status = "error",
            report = "Could not extract residuals from model"
        ))
    }
    
    # Remove missing values
    residuals_clean <- as.numeric(na.omit(residuals_vec))
    n_res <- length(residuals_clean)
    
    # Shapiro-Wilk test requires at least 3 observations
    if (n_res < 3) {
        return(list(
            shapiro_p_value = NA_real_,
            residuals_normal = NA,
            n_residuals = n_res,
            test_status = "error",
            report = sprintf("Insufficient residuals for Shapiro-Wilk test (n=%d, need >=3)", n_res)
        ))
    }
    
    # Run Shapiro-Wilk test
    test_result <- tryCatch({
        stats::shapiro.test(residuals_clean)
    }, error = function(e) {
        return(NULL)
    })
    
    if (is.null(test_result)) {
        return(list(
            shapiro_p_value = NA_real_,
            residuals_normal = NA,
            n_residuals = n_res,
            test_status = "error",
            report = "Shapiro-Wilk test execution failed"
        ))
    }
    
    # Extract test statistics
    p_value <- test_result$p.value
    is_normal <- p_value > 0.05  # Fail to reject H0 at alpha=0.05
    
    if (verbose) {
        status_text <- if (is_normal) "PASS [OK]" else "FAIL ?"
        message(sprintf("[.test_residual_normality] %s (p=%.4f, n=%d residuals)",
                       status_text, p_value, n_res))
    }
    
    return(list(
        shapiro_p_value = p_value,
        residuals_normal = is_normal,
        n_residuals = n_res,
        test_status = if (is_normal) "pass" else "fail",
        report = sprintf(
            "Shapiro-Wilk test: p=%.4f, %s normal (n=%d residuals)",
            p_value, if(is_normal) "residuals appear" else "residuals NOT",
            n_res
        )
    ))
}

# Helper: Check visual monotonicity of entropy values
# Purpose: Detect ordering issues or data quality problems before statistical testing
.check_monotonicity <- function(entropy_vals, q_vals, tolerance = 0.05) {
    # Args:
    #   entropy_vals: numeric vector of entropy values
    #   q_vals: numeric vector of q-values (should match entropy_vals length)
    #   tolerance: proportion of non-monotone pairs tolerated (default 5%)
    # Returns:
    #   List with: is_monotone (logical), n_violations, violation_indices, report (string)
    
    if (is.null(entropy_vals) || length(entropy_vals) < 2) {
        return(list(
            is_monotone = NA,
            n_values = 0,
            n_violations = 0,
            violation_indices = integer(0),
            report = "Insufficient data for monotonicity check"
        ))
    }
    
    # Sort by q-values to ensure proper ordering
    order_idx <- order(q_vals)
    entropy_sorted <- entropy_vals[order_idx]
    q_sorted <- q_vals[order_idx]
    
    # Compute differences: d_i = H_{q_{i+1}} - H_{q_i}
    # For monotone decreasing: all d_i < 0
    diffs <- diff(entropy_sorted)
    
    # Identify violations (positive differences indicate increase instead of decrease)
    # Use numerical tolerance to avoid false positives from floating-point errors
    violation_tolerance <- 1e-10
    violations <- which(diffs >= violation_tolerance)
    n_violations <- length(violations)
    n_total_pairs <- length(diffs)
    violation_rate <- n_violations / n_total_pairs
    
    # Determine if monotonicity holds (with tolerance)
    is_monotone <- violation_rate <= tolerance
    
    return(list(
        is_monotone = is_monotone,
        n_values = length(entropy_vals),
        n_violations = n_violations,
        n_pairs = n_total_pairs,
        violation_rate = violation_rate,
        violation_indices = violations + 1,  # Convert to 1-based indexing
        report = sprintf(
            "Monotonicity check: %d/%d pairs decreasing (%.1f%% violations). %s monotone.",
            n_total_pairs - n_violations, n_total_pairs, 100 * violation_rate,
            if(is_monotone) "PASS:" else "FAIL:"
        )
    ))
}

# Helper: Augmented Dickey-Fuller (ADF) test for unit root
# Simple implementation without external package dependencies
.adf_test <- function(time_series, max_lag = 3, alpha = 0.05) {
    # Args:
    #   time_series: numeric vector (observations)
    #   max_lag: maximum lag order for augmentation (default 3)
    #   alpha: significance level (default 0.05)
    # Returns:
    #   List with: test_stat, p_value, lag_used, conclusion, report (string)
    # 
    # H0: Unit root present (non-stationary)
    # Reject H0 -> series is stationary
    # Fail to reject H0 -> series is non-stationary (may have unit root)
    
    if (is.null(time_series) || length(na.omit(time_series)) < 5) {
        return(list(
            test_stat = NA_real_,
            p_value = NA_real_,
            lag_used = NA_integer_,
            stationary = NA,
            conclusion = "INSUFFICIENT_DATA",
            report = "Insufficient data for ADF test"
        ))
    }
    
    # Remove NAs
    ts_clean <- na.omit(as.numeric(time_series))
    
    if (length(ts_clean) < 5) {
        return(list(
            test_stat = NA_real_,
            p_value = NA_real_,
            lag_used = NA_integer_,
            stationary = NA,
            conclusion = "INSUFFICIENT_DATA",
            report = "Insufficient non-NA data for ADF test"
        ))
    }
    
    # Simplified ADF: regress Deltay_t on y_{t-1} and lagged differences
    # y_t = c + beta*y_{t-1} + Sum alpha_i*Deltay_{t-i} + ?_t
    # Test: H0: beta = 0 (unit root, non-stationary)
    
    n <- length(ts_clean)
    y <- ts_clean
    dy <- diff(y)  # First differences
    
    # Use lag order 1 (balance between flexibility and power)
    # In practice, lag selection would use AIC/BIC
    lag_order <- min(max_lag, max(1, floor(sqrt(n))))
    
    # Build regression matrix
    # Dependent variable: dy[2:n] (Deltay_t for t=2,...,n)
    # Predictor 1: y[seq_len(n-1)] (y_{t-1})
    # Predictor 2+: lagged differences dy[seq_len(n-lag_order-1)], etc.
    
    y_lag1 <- y[seq_len(n - 1)]
    dy_response <- dy[2:length(dy)]  # Deltay_t for t=2
    y_lag1_response <- y_lag1[2:length(y_lag1)]  # y_{t-1} aligned with Deltay_t
    
    # Simple regression: just use y_{t-1} without augmentation for stability
    valid_idx <- !is.na(dy_response) & !is.na(y_lag1_response)
    if (sum(valid_idx) < 3) {
        return(list(
            test_stat = NA_real_,
            p_value = NA_real_,
            lag_used = lag_order,
            stationary = NA,
            conclusion = "REGRESSION_FAILED",
            report = "Insufficient valid data for regression"
        ))
    }
    
    dy_model <- dy_response[valid_idx]
    y_lag_model <- y_lag1_response[valid_idx]
    
    # Fit: Deltay_t = beta * y_{t-1} + ?_t
    fit <- try(lm(dy_model ~ y_lag_model), silent = TRUE)
    if (inherits(fit, "try-error")) {
        return(list(
            test_stat = NA_real_,
            p_value = NA_real_,
            lag_used = lag_order,
            stationary = NA,
            conclusion = "REGRESSION_FAILED",
            report = "Regression failed"
        ))
    }
    
    # Extract t-statistic for beta (coefficient on y_lag_model)
    coef_table <- tryCatch(coef(summary(fit)), error = function(e) NULL)
    if (is.null(coef_table) || nrow(coef_table) < 2) {
        return(list(
            test_stat = NA_real_,
            p_value = NA_real_,
            lag_used = lag_order,
            stationary = NA,
            conclusion = "STAT_EXTRACTION_FAILED",
            report = "Could not extract test statistic"
        ))
    }
    
    # Critical values for ADF test (MacKinnon 1996, 5% level)
    # These are approximate; exact values depend on regression specification
    critical_values_5pct <- c(-2.86, -2.57, -2.57)  # For n~25, 50, 100+
    critical_value <- -2.86  # Conservative for small n
    
    t_stat <- as.numeric(coef_table[2, 3])  # t-statistic for y_lag_model coefficient
    
    # Approximate p-value based on t-statistic comparison to critical value
    # This is a simplified approximation; exact p-values require special distribution
    if (is.na(t_stat)) {
        p_value <- NA_real_
        stationary <- NA
    } else {
        # If t_stat < critical_value: REJECT H0 (stationary)
        # If t_stat > critical_value: FAIL TO REJECT H0 (non-stationary)
        stationary <- t_stat < critical_value
        # Approximate p-value (crude)
        p_value <- 2 * pt(t_stat, df = length(dy_model) - 2)  # Two-tailed
        p_value <- max(0.001, min(0.999, p_value))  # Bound to [0.001, 0.999]
    }
    
    return(list(
        test_stat = t_stat,
        p_value = p_value,
        lag_used = lag_order,
        critical_value = critical_value,
        stationary = stationary,
        conclusion = if(is.na(stationary)) "FAILED" else if(stationary) "REJECT_H0:_STATIONARY" else "FAIL_REJECT_H0:_NON-STATIONARY",
        report = sprintf(
            "ADF test (lag=%d): t=%.3f, crit=%.3f. %s -> %s",
            lag_order, t_stat, critical_value,
            if(is.na(stationary)) "FAILED" else if(t_stat < critical_value) "REJECT H0" else "FAIL REJECT H0",
            if(is.na(stationary)) "inconclusive" else if(stationary) "STATIONARY (rejects unit root)" else "NON-STATIONARY (has unit root)"
        )
    ))
}

# Helper: KPSS Test for stationarity (reverse of ADF)
# H0: Series IS stationary
.kpss_test <- function(time_series, trend = "constant", alpha = 0.05) {
    # Args:
    #   time_series: numeric vector
    #   trend: "constant" or "ct" (constant + time trend)
    #   alpha: significance level
    # Returns:
    #   List with: test_stat, p_value, conclusion, report (string)
    # 
    # H0: Series is stationary
    # Reject H0 -> series is NON-stationary
    # Fail to reject H0 -> series is stationary
    
    if (is.null(time_series) || length(na.omit(time_series)) < 5) {
        return(list(
            test_stat = NA_real_,
            p_value = NA_real_,
            stationary = NA,
            conclusion = "INSUFFICIENT_DATA",
            report = "Insufficient data for KPSS test"
        ))
    }
    
    ts_clean <- na.omit(as.numeric(time_series))
    
    if (length(ts_clean) < 5) {
        return(list(
            test_stat = NA_real_,
            p_value = NA_real_,
            stationary = NA,
            conclusion = "INSUFFICIENT_DATA",
            report = "Insufficient non-NA data for KPSS test"
        ))
    }
    
    # Simplified KPSS: Compute cumulative residuals variance ratio
    n <- length(ts_clean)
    y <- ts_clean
    
    # Demean if trend="constant", detrend if trend="ct"
    if (trend == "constant") {
        y_residual <- y - mean(y, na.rm = TRUE)
    } else {
        time_idx <- seq_len(n)
        fit <- try(lm(y ~ time_idx), silent = TRUE)
        if (inherits(fit, "try-error")) {
            y_residual <- y - mean(y, na.rm = TRUE)
        } else {
            y_residual <- residuals(fit)
        }
    }
    
    # Cumulative sum of residuals
    S_t <- cumsum(y_residual)
    
    # Long-run variance estimate (Newey-West with lag=1)
    s2 <- mean(y_residual^2)  # Short-run variance
    
    # Autocovariance at lag 1
    if (n > 1) {
        gamma1 <- sum(y_residual[-n] * y_residual[-1]) / n
    } else {
        gamma1 <- 0
    }
    
    # Long-run variance (Newey-West)
    sigma2_lr <- s2 + 2 * gamma1 * (1 - 1 / n)
    sigma2_lr <- max(s2 * 0.1, sigma2_lr)  # Ensure positive
    
    # KPSS statistic
    kpss_stat <- sum(S_t^2) / (n^2 * sigma2_lr)
    
    # Critical values for KPSS (Kwiatkowski et al. 1992)
    # trend="constant": 10%, 5%, 2.5%, 1% are 0.347, 0.463, 0.574, 0.739
    # With trend: 0.119, 0.146, 0.176, 0.216
    if (trend == "constant") {
        crit_5pct <- 0.463
    } else {
        crit_5pct <- 0.146
    }
    
    # Reject H0 (stationarity) if KPSS stat > critical value
    reject_h0 <- kpss_stat > crit_5pct
    stationary <- !reject_h0
    
    # Approximate p-value
    if (kpss_stat < 0.347) {
        p_value <- 0.10
    } else if (kpss_stat < 0.463) {
        p_value <- 0.05
    } else if (kpss_stat < 0.574) {
        p_value <- 0.025
    } else if (kpss_stat < 0.739) {
        p_value <- 0.01
    } else {
        p_value <- 0.001
    }
    if (!reject_h0) p_value <- 1 - p_value
    
    return(list(
        test_stat = kpss_stat,
        p_value = p_value,
        critical_value = crit_5pct,
        stationary = stationary,
        conclusion = if(reject_h0) "REJECT_H0:_NON-STATIONARY" else "FAIL_REJECT_H0:_STATIONARY",
        report = sprintf(
            "KPSS test (trend=%s): LM=%.3f, crit=%.3f. %s -> %s",
            trend, kpss_stat, crit_5pct,
            if(reject_h0) "REJECT H0" else "FAIL REJECT H0",
            if(stationary) "STATIONARY" else "NON-STATIONARY"
        )
    ))
}

# Comprehensive stationarity validation
# Returns diagnostic report comparing raw and differenced data
.validate_stationarity <- function(entropy_vals, q_vals, subject_vec = NULL, gene_name = NULL) {
    # Args:
    #   entropy_vals: raw entropy values
    #   q_vals: corresponding q-values
    #   subject_vec: (optional) subject identifiers for within-subject validation
    #   gene_name: (optional) identifier for reporting
    # Returns:
    #   List with full diagnostic report including all tests
    
    if (is.null(gene_name)) gene_name <- "Unknown"
    
    # Test 1: Monotonicity
    mono_check <- .check_monotonicity(entropy_vals, q_vals)
    
    # Test 2-3: ADF and KPSS on raw data
    adf_raw <- .adf_test(entropy_vals)
    kpss_raw <- .kpss_test(entropy_vals, trend = "constant")
    
    # Test 4-5: ADF and KPSS on first differences
    if (length(entropy_vals) > 1) {
        # Sort by q first
        sort_idx <- order(q_vals)
        entropy_sorted <- entropy_vals[sort_idx]
        entropy_diff <- diff(entropy_sorted)
        
        adf_diff <- .adf_test(entropy_diff)
        kpss_diff <- .kpss_test(entropy_diff, trend = "constant")
    } else {
        adf_diff <- list(test_stat = NA, p_value = NA, stationary = NA, conclusion = "NO_DATA")
        kpss_diff <- list(test_stat = NA, p_value = NA, stationary = NA, conclusion = "NO_DATA")
        entropy_diff <- NULL
    }
    
    # Summary: Check if ARIMA(1,1,0) is justified
    arima_justified <- !mono_check$is_monotone &&  # Raw is non-monotone
                       !isTRUE(adf_raw$stationary) &&  # Raw is non-stationary (ADF fails to reject H0)
                       isTRUE(kpss_raw$stationary) == FALSE &&  # Raw is non-stationary (KPSS rejects H0)
                       isTRUE(adf_diff$stationary)  # Differences are stationary (ADF rejects H0)
    
    return(list(
        gene = gene_name,
        n_values = length(entropy_vals),
        n_q_values = length(unique(q_vals)),
        monotonicity = mono_check,
        raw_data_tests = list(
            adf = adf_raw,
            kpss = kpss_raw,
            interpretation = sprintf(
                "Raw entropy: ADF=%s, KPSS=%s. %s",
                adf_raw$conclusion, kpss_raw$conclusion,
                if(!isTRUE(adf_raw$stationary) && isTRUE(kpss_raw$stationary) == FALSE)
                    "CONFIRMED non-stationary (unit root likely)" else "QUESTIONABLE stationarity"
            )
        ),
        differenced_data_tests = list(
            adf = adf_diff,
            kpss = kpss_diff,
            interpretation = sprintf(
                "Differenced entropy: ADF=%s, KPSS=%s. %s",
                adf_diff$conclusion, kpss_diff$conclusion,
                if(isTRUE(adf_diff$stationary) && isTRUE(kpss_diff$stationary) == FALSE)
                    "CONFIRMED stationary (differencing effective)" else "QUESTIONABLE stationarity after differencing"
            )
        ),
        arima_justified = arima_justified,
        recommendation = if(arima_justified)
            "[OK] Use ARIMA(1,1,0): differencing removes trend, AR(1) appropriate for residuals" else
            "? REVIEW: Stationarity assumptions may not hold, consider alternative modeling",
        report = sprintf(
            "STATIONARITY VALIDATION for %s:\n%s\nRaw: %s, %s\nDiff: %s, %s\n%s",
            gene_name,
            mono_check$report,
            adf_raw$report, kpss_raw$report,
            adf_diff$report, kpss_diff$report,
            if(arima_justified) "[OK] ARIMA(1,1,0) assumptions validated" else "? Issues detected"
        )
    ))
}

.compute_arima_differences <- function(df, q_vals, group_vec, subject_vec = NULL) {
    # ARIMA(1,1,0) implementation: compute first differences of entropy
    # 
    # Background:
    # - Tsallis entropy H_q is monotone decreasing in q (non-stationary)
    # - AR(1) assumes stationarity (constant mean, variance)
    # - Solution: Apply AR(1) to first differences DeltaH_q = H_q - H_{q-1}
    # - Result: ARIMA(1,1,0) = Integrated AR(1) = AR(1) on differenced data
    #
    # Implementation notes:
    # 1. Order data by q-values to ensure proper differencing
    # 2. Compute differences within each subject (not across subjects)
    # 3. Return data frame with differenced entropy, q values, group, subject
    # 4. Note: Loses 1 observation per subject (trade-off for stationarity)
    #
    # Returns: list(df_diff, n_lost_obs) or NULL if insufficient data
    
    if (is.null(df) || nrow(df) == 0) {
        return(NULL)
    }
    
    # Add q and group to data frame for sorting
    df_full <- data.frame(
        entropy = as.numeric(df$entropy),
        q = as.numeric(q_vals),
        group = factor(group_vec),
        subject = if (!is.null(subject_vec)) factor(subject_vec) else factor(seq_len(nrow(df))),
        stringsAsFactors = FALSE
    )
    
    # Remove NA entropy values
    df_full <- df_full[!is.na(df_full$entropy), ]
    
    if (nrow(df_full) < 2) {
        return(NULL)
    }
    
    # Sort by subject and q to ensure proper differencing within subjects
    df_full <- df_full[order(df_full$subject, df_full$q), ]
    
    # Compute first differences within each subject
    df_diff_list <- list()
    n_lost <- 0
    
    for (subj in levels(df_full$subject)) {
        subj_idx <- which(df_full$subject == subj)
        
        if (length(subj_idx) < 2) {
            # Skip subjects with < 2 observations (can't compute difference)
            n_lost <- n_lost + length(subj_idx)
            next
        }
        
        # Extract subject data (should already be sorted by q)
        subj_data <- df_full[subj_idx, ]
        
        # Compute differences: DeltaH_q = H_q - H_{q-1}
        n_diff <- nrow(subj_data) - 1
        
        df_diff_list[[subj]] <- data.frame(
            entropy_diff = diff(subj_data$entropy),  # DeltaH_q
            q = subj_data$q[-1],                      # q indices for differences
            q_prev = subj_data$q[-nrow(subj_data)],   # q_{q-1} for reference
            group = subj_data$group[-nrow(subj_data)],  # Group for first q in pair
            subject = rep(subj, n_diff),
            stringsAsFactors = FALSE
        )
    }
    
    if (length(df_diff_list) == 0) {
        return(NULL)
    }
    
    # Combine all subject differences
    df_diff <- do.call(rbind, df_diff_list)
    rownames(df_diff) <- NULL
    
    if (nrow(df_diff) == 0) {
        return(NULL)
    }
    
    # Rename entropy_diff to entropy for compatibility with model fitting
    names(df_diff)[names(df_diff) == "entropy_diff"] <- "entropy"
    
    return(list(
        df = df_diff,
        n_observations_original = nrow(df_full) + n_lost,
        n_observations_differenced = nrow(df_diff),
        n_observations_lost = n_lost,
        transformation = "ARIMA(1,1,0): First differences"
    ))
}

# Helper: Check if entropy data is truly bounded in [0, 1]
# Returns TRUE if data appears normalized/proportional
.is_bounded_0_1 <- function(entropy_vals) {
    entropy_clean <- na.omit(entropy_vals)
    if (length(entropy_clean) == 0) return(FALSE)
    finite_clean <- is.finite(entropy_clean)
    if (!any(finite_clean)) return(FALSE)
    min_val <- min(entropy_clean[finite_clean])
    max_val <- max(entropy_clean[finite_clean])
    tolerance <- 0.01
    bounds_check <- min_val >= -tolerance && max_val <= 1 + tolerance
    if (!bounds_check) return(FALSE)
    approaches_lower_bound <- min_val <= 0.1
    approaches_upper_bound <- max_val >= 0.9
    return(approaches_lower_bound || approaches_upper_bound)
}

# ============================================================================
# MEMOIZATION: Cache expensive computations
# ============================================================================
# Memoization reduces redundant calculations in multi-q iterative analyses
# Expected speedup: 10-30% for typical multi-q analyses

if (getOption("TSENAT.memoization", TRUE)) {
    # Cache knot selection: input = (entropy, q_vals, n_unique, ...)
    # Avoids recomputing knots for same entropy data across iterations
    .adaptive_spline_knots_memo <- memoise::memoise(
        .adaptive_spline_knots,
        cache = memoise::cache_memory()
    )
    
    # Cache design effect: input = (rho, cluster_size)
    # Avoids recomputing design effect for repeated (rho, cluster_size) pairs
    .ar1_design_effect_memo <- memoise::memoise(
        .ar1_design_effect,
        cache = memoise::cache_memory()
    )
} else {
    # Fallback: no memoization if disabled globally
    .adaptive_spline_knots_memo <- .adaptive_spline_knots
    .ar1_design_effect_memo <- .ar1_design_effect
}

#' Internal: Clear memoization cache
#' @description Invalidates cached results for new dataset processing
#' @noRd
.clear_lm_helper_cache <- function() {
    if (getOption("TSENAT.memoization", TRUE)) {
        memoise::forget(.adaptive_spline_knots_memo)
        memoise::forget(.ar1_design_effect_memo)
    }
}



# Helper: Compute skewness of a vector
# Positive skew: right tail longer (mode < median < mean)
# Negative skew: left tail longer (mean < median < mode)
.compute_skewness <- function(x, na.rm = TRUE) {
    if (na.rm) x <- na.omit(x)
    if (length(x) < 3) return(NA)
    
    m <- mean(x)
    s <- sd(x)
    n <- length(x)
    
    if (s == 0) return(0)
    
    # Unbiased skewness estimate
    skew <- (sum((x - m)^3) / n) / (s^3)
    return(skew)
}

# ================================================================================
# HETEROSCEDASTICITY DETECTION AND VARIANCE WEIGHTING (March 2026)
# ================================================================================
# Tsallis entropy often exhibits variance that depends on:
#   1. Mean entropy level (mean-variance relationship)
#   2. q-value (variance changes across diversity orders)
#   3. Group/condition (treatment vs control may have different variance)
#
# Consequence: Tests can be biased with inflated Type I error if heteroscedasticity ignored
# Solution: Detect heteroscedasticity and apply appropriate variance adjustment/weighting

# Detect heteroscedasticity using Breusch-Pagan test
.detect_heteroscedasticity <- function(df, q_vals, group_vec, verbose = FALSE) {
    # Fit OLS to get residuals
    fit_ols <- try(
        lm(entropy ~ q + group, data = df),
        silent = TRUE
    )
    
    if (inherits(fit_ols, "try-error")) {
        return(list(
            is_heteroscedastic = NA,
            bp_stat = NA,
            p_value = NA,
            var_ratio_q = NA,
            var_ratio_group = NA
        ))
    }
    
    residuals_sq <- residuals(fit_ols)^2
    
    # Breusch-Pagan auxiliary regression: log(residuals^2) ~ q + group
    fit_aux <- try(
        lm(log(residuals_sq + 1e-8) ~ q + factor(group), data = df),
        silent = TRUE
    )
    
    if (inherits(fit_aux, "try-error")) {
        return(list(
            is_heteroscedastic = NA,
            bp_stat = NA,
            p_value = NA,
            var_ratio_q = NA,
            var_ratio_group = NA
        ))
    }
    
    # BP statistic = RSS from auxiliary model / (2 * RSS from original model)
    rss_aux <- sum(residuals(fit_aux)^2)
    # OPTIMIZATION (March 2026): Cache computation to avoid redundant calculation
    fitted_sq_sum <- sum((fitted(fit_aux) - mean(fitted(fit_aux)))^2)
    tss_aux <- fitted_sq_sum + rss_aux
    
    bp_stat <- (fitted_sq_sum / tss_aux * nrow(df))
    # Compute correct degrees of freedom: number of predictors in auxiliary regression
    # BUG FIX: Was hardcoded to 2, but should be ncol(X) - 1 where X is model.matrix
    df_bp <- ncol(model.matrix(fit_aux)) - 1
    p_value <- 1 - pchisq(bp_stat, df = df_bp)
    
    # Compute variance ratios
    # OPTIMIZATION (March 2026): Use tapply() instead of sapply + subsetting (2-3x faster)
    residuals_vec <- residuals(fit_ols)
    var_by_q <- tapply(residuals_vec, df$q, var)
    finite_q <- is.finite(var_by_q)
    if (any(finite_q)) {
        max_q <- max(var_by_q[finite_q])
        min_q <- min(var_by_q[finite_q])
        var_ratio_q <- max_q / (min_q + 1e-8)
    } else {
        var_ratio_q <- NA
    }

    var_by_group <- tapply(residuals_vec, df$group, var)
    finite_g <- is.finite(var_by_group)
    if (any(finite_g)) {
        max_g <- max(var_by_group[finite_g])
        min_g <- min(var_by_group[finite_g])
        var_ratio_group <- max_g / (min_g + 1e-8)
    } else {
        var_ratio_group <- NA
    }
    
    if (verbose) {
        message(sprintf("[Heteroscedasticity] BP p-value: %.4f, Var ratio (q): %.2f, Var ratio (group): %.2f",
                        p_value, var_ratio_q, var_ratio_group))
    }
    
    return(list(
        is_heteroscedastic = p_value < 0.05,
        bp_stat = bp_stat,
        p_value = p_value,
        var_ratio_q = var_ratio_q,
        var_ratio_group = var_ratio_group
    ))
}

# Estimate variance weights for heteroscedasticity adjustment
.estimate_variance_weights <- function(df, q_vals, method = "power", verbose = FALSE) {
    # Estimate weights to model variance heterogeneity
    # method = "power": Model Var ~ q^?, compute weights w_i = q_i^(-?)
    # method = "residual": Use residual variance from OLS as observation weights
    
    if (method == "power") {
        # Estimate power parameter ? via regression: log(residuals_sq) ~ q
        # First, fit OLS to get residuals
        fit_ols <- try(
            lm(entropy ~ q + group, data = df),
            silent = TRUE
        )
        
        # If OLS with group fails, try just q
        if (inherits(fit_ols, "try-error")) {
            fit_ols <- try(
                lm(entropy ~ q, data = df),
                silent = TRUE
            )
        }
        
        if (inherits(fit_ols, "try-error")) {
            return(NULL)
        }
        
        residuals_ols <- residuals(fit_ols)
        residuals_sq <- residuals_ols^2
        
        # This gives: log(Var) = log(sigma2) + ? * log(q)
        # So: Var ~ sigma2 * q^?
        # Weights: w_i = 1 / (sigma2 * q_i^?) ? q_i^(-?)
        
        w <- 1 / (residuals_sq + 1e-8)
        wfit <- try(
            lm(log(residuals_sq + 1e-8) ~ log(df$q + 1e-8), weights = w),
            silent = TRUE
        )
        
        if (!inherits(wfit, "try-error")) {
            theta_est <- coef(wfit)[2]
            if (!is.na(theta_est)) {
                # Ensure positive weighting
                theta_est <- max(theta_est, 0.01)
                weights <- 1 / (df$q^theta_est + 1e-8)
                weights <- weights / mean(weights, na.rm = TRUE)  # Standardize
                
                if (verbose) {
                    message(sprintf("[Variance Weighting] Estimated power parameter ? = %.3f", theta_est))
                }
                
                return(list(
                    weights = weights,
                    power_param = theta_est,
                    method = "power"
                ))
            }
        }
    }
    
    if (method == "residual") {
        # Use inverse variance as weights
        # Try OLS fit to estimate residual variance
        ols_fit <- NULL
        
        # Try with group if available, otherwise just q
        if ("group" %in% colnames(df)) {
            ols_fit <- try(
                lm(entropy ~ q + group, data = df),
                silent = TRUE
            )
        } else {
            ols_fit <- try(
                lm(entropy ~ q, data = df),
                silent = TRUE
            )
        }
        
        residuals_sq <- NA_real_
        if (!is.null(ols_fit) && !inherits(ols_fit, "try-error")) {
            residuals_sq <- residuals(ols_fit)^2
        }
        
        # If OLS succeeded and we have residuals, compute weights
        if (!all(is.na(residuals_sq))) {
            weights <- 1 / (residuals_sq + 1e-8)
            weights <- weights / mean(weights, na.rm = TRUE)  # Standardize
            
            return(list(
                weights = weights,
                method = "residual"
            ))
        }
    }
    
    # Fallback: uniform weights
    return(list(
        weights = rep(1, nrow(df)),
        method = "uniform"
    ))
}




# Fit function extracted from calculate_lm_interaction
.fit_one_interaction <- function(g, se, mat, q_vals, sample_names, group_vec,
    method, pvalue, subject_col, paired, min_obs, verbose, suppress_lme4_warnings,
    progress, bias_correction = TRUE, regularization = c("pca", "lasso", "elasticnet", "gamsel", "spline"),
    corstr = c("ar1", "exchangeable", "independence"), adaptive_knots = TRUE, weights = NULL) {
    regularization <- match.arg(regularization)
    corstr <- match.arg(corstr)
    vals <- as.numeric(mat[g, ])
    df <- data.frame(entropy = vals, q = q_vals, group = factor(group_vec))
    
    
    # Add inverse-variance weights if provided (Phase 1: Bootstrap CI weighting)
    if (!is.null(weights) && length(weights) == nrow(df)) {
        df$weight <- weights
        if (verbose) {
            message(sprintf("[.fit_one_interaction] Gene '%s': weights applied (n=%d, mean=%.4f, min=%.4f, max=%.4f)",
                           g, length(weights), mean(weights, na.rm=TRUE), min(weights, na.rm=TRUE), max(weights, na.rm=TRUE)))
        }
    } else {
        if (verbose && !is.null(weights)) {
            message(sprintf("[.fit_one_interaction] Gene '%s': weights NOT applied - length mismatch (weights=%d, df rows=%d)",
                           g, length(weights), nrow(df)))
        }
    }
    

    if (method == "lmm") {
        # LMM with AR(1) covariance structure for q-dependent entropy measurements
        # Paper S171 (Zimmerman & Harville, 1991): "Linear Models with Generalized AR(1) 
        # Covariance Structure for Longitudinal and Spatial Data" validates this approach.
        # Papers S168-S170: Theoretical foundation and empirical estimation of AR(1) parameters.
        # TEST L.1.6: Confirms q-value correlation follows AR(1) pattern (rho(k) = phi^|k|).
        #
        # CRITICAL FIX (March 2026): Implements true ARIMA(1,1,0) by:
        # 1. Computing first differences DeltaH_q = H_q - H_{q-1} within each subject
        # 2. Ensuring stationarity: DeltaH_q has constant mean (unlike monotone H_q)
        # 3. Applying AR(1) to differenced data (not raw entropy)
        # 4. Fitting all models (null and alt) on differenced entropy
        #
        if (!requireNamespace("nlme", quietly = TRUE)) {
            stop("Package 'nlme' is required for method = 'lmm' (AR(1) covariance support)")
        }
        # determine subject IDs for random effect
        subject <- NULL
        if (!is.null(subject_col)) {
            if (!(subject_col %in% colnames(SummarizedExperiment::colData(se)))) {
                stop(sprintf("subject_col '%s' not found in colData(se)", subject_col))
            }
            subj_full <- as.character(SummarizedExperiment::colData(se)[, subject_col])
            # CRITICAL FIX: When colData doesn't have 'samples' column, use colnames(mat)
            # but strip "_q=" suffixes to match how sample_names were parsed in calculate_lm_interaction
            col_names_for_indexing <- SummarizedExperiment::colData(se)$samples
            if (is.null(col_names_for_indexing)) {
                # Strip "_q=....." suffixes from column names for proper indexing
                col_names_for_indexing <- sub("_q=.*", "", colnames(mat))
            }
            names(subj_full) <- col_names_for_indexing
            subject <- unname(subj_full[sample_names])
        } else if (paired) {
            coldata <- SummarizedExperiment::colData(se)
            coldata_cols <- colnames(coldata)
            
            # Look for 'paired_samples' or 'sample_base' columns created by map_metadata()
            subject_col_name <- NULL
            if ("paired_samples" %in% coldata_cols) {
                subject_col_name <- "paired_samples"
            } else if ("sample_base" %in% coldata_cols) {
                subject_col_name <- "sample_base"
            } else if (length(coldata_cols) >= 3) {
                subject_col_name <- coldata_cols[3]  # fallback
            }
            
            if (!is.null(subject_col_name)) {
                subject_ids <- as.character(coldata[, subject_col_name])
                # Map sample names by setting names to base names (without _q= suffix)
                sample_names_expanded <- sub("_q=.*", "", rownames(coldata))
                names(subject_ids) <- sample_names_expanded
                subject <- unname(subject_ids[sample_names])  # sample_names passed from calculate_lm_interaction

            } else {
                stop("paired = TRUE requires 'paired_samples' or 'sample_base' column in colData; supply subject_col explicitly or use map_metadata(...)")
            }
        } else {
            # fallback to sample base derived from column names
            subject <- sample_names
        }
        df$subject <- factor(subject)
        # require at least two subjects and at least two groups represented
        n_subjects <- length(unique(na.omit(df$subject)))
        
        # Check minimum observation requirement
        if (nrow(df) < min_obs) {
            return(NULL)
        }
        
        # Check minimum subject requirement for mixed models
        if (n_subjects < 2) {
            return(NULL)
        }
                
        # ARIMA(1,1,0) IMPLEMENTATION: Compute first differences for stationarity
        # Differencing removes monotone trend from Tsallis entropy, enabling valid AR(1) inference
        # BONUS: First differencing of bounded [0, log(m)] data helps normalize distribution
        # (bounded support becomes approximately normal after differencing in many cases)
        arima_result <- .compute_arima_differences(df, q_vals, df$group, df$subject)
        
        if (is.null(arima_result) || nrow(arima_result$df) < 3) {
            # Insufficient data for ARIMA differencing; fall back to raw data with warning
            if (verbose) {
                message("[calculate_lm_interaction] ARIMA(1,1,0) differencing lost too many observations; using raw entropy")
            }
            df_model <- df
            use_arima <- FALSE
        } else {
            df_model <- arima_result$df
            use_arima <- TRUE
            if (verbose) {
                message(sprintf("[calculate_lm_interaction] ARIMA(1,1,0): %d observations -> %d after differencing", 
                    arima_result$n_observations_original, arima_result$n_observations_differenced))
            }
        }
        
        # fit null (no interaction) and alternative (with q:group interaction)
        # Using nlme::lme() for AR(1) covariance structure support (instead of lme4::lmer)
        mm_suppress_pattern <- "boundary \\(singular\\) fit|Computed variance-covariance matrix problem|not a positive definite matrix"

        # Apply regularization for feature selection if requested (not "pca")
        fs_result <- NULL
        # nlme formula syntax: fixed effects ~ random intercept
        formula_null <- entropy ~ q + group
        formula_alt <- entropy ~ q * group
        
        if (regularization != "pca") {
            fs_result <- .lmm_regularization(q_vals = df_model$q, entropy_vals = df_model$entropy,
                                                    group_vec = df_model$group, subject_vec = df_model$subject,
                                                    regularization = regularization)
            if (!is.null(fs_result)) {
                # If we got feature selection results, modify formulas to use only selected interaction terms
                # For simplicity, we construct a reduced model with selected q-value ranges
                # Note: This keeps main effects (q, group) but regularizes the interaction
                uq_levels <- length(fs_result$q_values)
                if (uq_levels > 2 && length(fs_result$selected_features) < (uq_levels - 1)) {
                    # Regularization reduced model complexity - use main effects plus interaction
                    # (full formula still, but now justified by regularization path)
                    if (verbose) {
                        message("[calculate_lmm_interaction] regularization retained ", 
                                length(fs_result$selected_features), " of ", 
                                uq_levels - 1, " q-interaction features")
                    }
                }
            }
        }
        
        # ===============================================================================
        # HETEROSCEDASTICITY DETECTION AND VARIANCE WEIGHTING (NEW - March 2026)
        # ===============================================================================
        # Detect q-dependent and group-dependent variance heterogeneity
        # Apply nlme::varPower() to model variance heterogeneity if detected
        hetero_result <- .detect_heteroscedasticity(df_model, df_model$q, df_model$group)
        use_var_structure <- FALSE
        
        if (!is.na(hetero_result$is_heteroscedastic) && hetero_result$is_heteroscedastic) {
            use_var_structure <- TRUE
            if (verbose) {
                message(sprintf("[calculate_lmm_interaction] Heteroscedasticity detected (BP p = %.4f); applying varPower",
                               hetero_result$p_value))
            }
        }
        
        # Fit nlme models with AR(1) covariance structure for q-measurements within subjects
        # AR(1) model: Cov(Y_t, Y_s) = sigma^2 phi^|t-s| where t,s are q-ordered indices
        # *** CRITICAL: Now applied to differenced entropy DeltaH_q, not raw H_q ***
        # HETEROSCEDASTICITY: Add varPower() structure if heteroscedasticity detected
        
        if (use_var_structure) {
            # Include variance power model: Var(Y) ~ q^?
            fit0 <- try(
                nlme::lme(formula_null, random = ~1 | subject, data = df_model, method = "ML",
                         correlation = nlme::corAR1(form = ~1 | subject),
                         weights = nlme::varPower(form = ~ q)),
                silent = TRUE
            )
            fit1 <- try(
                nlme::lme(formula_alt, random = ~1 | subject, data = df_model, method = "ML",
                         correlation = nlme::corAR1(form = ~1 | subject),
                         weights = nlme::varPower(form = ~ q)),
                silent = TRUE
            )
        } else {
            fit0 <- try(
                nlme::lme(formula_null, random = ~1 | subject, data = df_model, method = "ML",
                         correlation = nlme::corAR1(form = ~1 | subject)),
                silent = TRUE
            )
            fit1 <- try(
                nlme::lme(formula_alt, random = ~1 | subject, data = df_model, method = "ML",
                         correlation = nlme::corAR1(form = ~1 | subject)),
                silent = TRUE
            )
        }

        # If either fit failed, try falling back to simpler approaches
        fallback_lm <- NULL
        used_fit_method <- "nlme::lme"
        used_singular <- FALSE
        if (inherits(fit0, "try-error") || inherits(fit1, "try-error")) {
            used_fit_method <- "fallback"
            used_singular <- FALSE
            if ((verbose && progress) || (!verbose && progress)) {
                message("[calculate_lm_interaction] mixed model failed; trying simpler fixed-effects fallback")
            }
            fb <- .try_lm_fallbacks(df_model, verbose = verbose)
            if (!is.null(fb)) {
                fallback_lm <- fb
                used_fit_method <- fb$method
            }
        } else {
            used_fit_method <- if (use_arima) "nlme::lme_arima(1,1,0)" else "nlme::lme_ar1_raw"
        }

        lrt_result <- list(p_value = NA_real_, n_subjects = NA_integer_, small_sample_flag = FALSE)
        msg <- NULL
        if (!is.null(fallback_lm)) {
            lrt_result <- .extract_lrt_p(fallback_lm$fit0, fallback_lm$fit1, df = df_model)
            # If glmmTMB fallback failed due to convergence, propagate message
            if (!is.null(fallback_lm$message)) {
                msg <- fallback_lm$message
            }
        } else {
            lrt_result <- .extract_lrt_p(fit0, fit1, df = df_model)
        }

        # nlme models use LRT for hypothesis testing (not Satterthwaite)
        # pvalue argument is ignored for nlme method
        satter_p <- NA_real_
        lrt_p <- lrt_result$p_value
        n_subj_lmm <- lrt_result$n_subjects
        small_sample_lmm <- lrt_result$small_sample_flag

        # nlme always uses LRT for hypothesis testing
        p_interaction <- lrt_p
        
        # Extract interaction coefficient (slope_diff) from fitted model
        slope_diff <- NA_real_
        if (!is.null(fallback_lm) && !is.null(fallback_lm$fit1)) {
            # For fallback lm/glm models
            coefs <- tryCatch(coef(fallback_lm$fit1), error = function(e) NULL)
            if (!is.null(coefs)) {
                # Look for q:group interaction term
                interaction_idx <- grep("q:group|group:q", names(coefs), ignore.case = FALSE)
                if (length(interaction_idx) > 0) {
                    slope_diff <- coefs[interaction_idx[1]]
                }
            }
        } else if (!inherits(fit1, "try-error")) {
            # For nlme models
            coefs <- tryCatch(nlme::fixef(fit1), error = function(e) NULL)
            if (!is.null(coefs)) {
                # Look for q:group interaction term
                interaction_idx <- grep("q:group|group:q", names(coefs), ignore.case = FALSE)
                if (length(interaction_idx) > 0) {
                    slope_diff <- coefs[interaction_idx[1]]
                }
            }
        }
        
        # Add weighting information to results (Phase 1)
        has_weights <- !is.null(df$weight)

        # PHASE 14 ENHANCEMENT: Document sample size and power flags
        # Report n_subjects for transparency and small_sample_flag for caution
        res <- data.frame(
            gene = g, 
            p_interaction = p_interaction, 
            p_lrt = lrt_p,
            p_satterthwaite = NA_real_, 
            slope_diff = slope_diff, 
            fit_method = used_fit_method, 
            singular = used_singular, 
            arima_transformation = use_arima, 
            ci_weighted = has_weights,
            n_subjects = n_subj_lmm,
            small_sample_flag = small_sample_lmm,
            stringsAsFactors = FALSE
        )
        if (!is.null(msg)) res$message <- msg
        return(res)
    }

    if (method == "gam") {
        # Extract subject info for paired/repeated measures (same as LMM)

        subject <- NULL
        
        if (!is.null(subject_col)) {
            if (!(subject_col %in% colnames(SummarizedExperiment::colData(se)))) {
                stop(sprintf("subject_col '%s' not found in colData(se)", subject_col))
            }
            subj_full <- as.character(SummarizedExperiment::colData(se)[, subject_col])
            # CRITICAL FIX: When colData doesn't have 'samples' column, use colnames(mat)
            # but strip "_q=" suffixes to match how sample_names were parsed in calculate_lm_interaction
            col_names_for_indexing <- SummarizedExperiment::colData(se)$samples
            if (is.null(col_names_for_indexing)) {
                # Strip "_q=....." suffixes from column names for proper indexing
                col_names_for_indexing <- sub("_q=.*", "", colnames(mat))
            }
            names(subj_full) <- col_names_for_indexing
            subject <- unname(subj_full[sample_names])
            
        } else if (paired) {
        
            coldata <- SummarizedExperiment::colData(se)
            coldata_cols <- colnames(coldata)
            
            # Look for 'paired_samples' or 'sample_base' columns created by map_metadata()
            subject_col_name <- NULL
            if ("paired_samples" %in% coldata_cols) {
                subject_col_name <- "paired_samples"
            } else if ("sample_base" %in% coldata_cols) {
                subject_col_name <- "sample_base"
            } else if (length(coldata_cols) >= 3) {
                subject_col_name <- coldata_cols[3]  # fallback
            }
            
            if (!is.null(subject_col_name)) {
                subject_ids <- as.character(coldata[, subject_col_name])
                # Map sample names by setting names to base names (without _q= suffix)
                sample_names_expanded <- sub("_q=.*", "", rownames(coldata))
                names(subject_ids) <- sample_names_expanded
                subject <- unname(subject_ids[sample_names])  # sample_names passed from calculate_lm_interaction

            } else {
                stop("paired = TRUE requires 'paired_samples' or 'sample_base' column in colData; supply subject_col explicitly or use map_metadata(...)")
            }
        }

        # Pass subject info, regularization, bias correction, adaptive knots parameters, and weights to GAM
        return(.gam_interaction(df, q_vals, g, min_obs = min_obs, subject = subject,
                                       regularization = regularization, bias_correction = bias_correction,
                                       adaptive_knots = adaptive_knots, weights = weights))
    }

    if (method == "fpca") {
        # Extract subject info for paired samples (same as LMM and GAMM)
        subject <- NULL
        if (!is.null(subject_col)) {
            if (!(subject_col %in% colnames(SummarizedExperiment::colData(se)))) {
                stop(sprintf("subject_col '%s' not found in colData(se)", subject_col))
            }
            subj_full <- as.character(SummarizedExperiment::colData(se)[, subject_col])
            # CRITICAL FIX: When colData doesn't have 'samples' column, use colnames(mat)
            # but strip "_q=" suffixes to match how sample_names were parsed in calculate_lm_interaction
            col_names_for_indexing <- SummarizedExperiment::colData(se)$samples
            if (is.null(col_names_for_indexing)) {
                # Strip "_q=....." suffixes from column names for proper indexing
                col_names_for_indexing <- sub("_q=.*", "", colnames(mat))
            }
            names(subj_full) <- col_names_for_indexing
            subject <- unname(subj_full[sample_names])
        } else if (paired) {
            coldata <- SummarizedExperiment::colData(se)
            coldata_cols <- colnames(coldata)
            
            # Look for 'paired_samples' or 'sample_base' columns created by map_metadata()
            subject_col_name <- NULL
            if ("paired_samples" %in% coldata_cols) {
                subject_col_name <- "paired_samples"
            } else if ("sample_base" %in% coldata_cols) {
                subject_col_name <- "sample_base"
            } else if (length(coldata_cols) >= 3) {
                subject_col_name <- coldata_cols[3]  # fallback
            }
            
            if (!is.null(subject_col_name)) {
                subject_ids <- as.character(coldata[, subject_col_name])
                # Map sample names by setting names to base names (without _q= suffix)
                sample_names_expanded <- sub("_q=.*", "", rownames(coldata))
                names(subject_ids) <- sample_names_expanded
                subject <- unname(subject_ids[sample_names])  # sample_names passed from calculate_lm_interaction
            } else {
                stop("paired = TRUE requires 'paired_samples' or 'sample_base' column in colData; supply subject_col explicitly or use map_metadata(...)")
            }
        }
        # FPCA respects q-value ordering: creates curve matrix with q-values as columns (ordered),
        # then applies PCA which naturally captures smooth functional dependence structure (S168-S171).
        # This implicitly models AR(1) correlation: rho(k) = phi^|k| across ordered q-values.
        # Test L.1.6 validates this AR(1) pattern for entropy across q-values.
        return(.fpca_interaction(mat, q_vals, sample_names, group_vec, g,
            min_obs = min_obs, subject = subject, regularization = regularization, weights = weights))
    }

    if (method == "gee") {
        # Extract subject info for GEE (clustered/repeated measures)
        subject <- NULL
        # Note: weights are passed directly to the GEE helper, no need to prepare locally
        if (!is.null(subject_col)) {
            if (!(subject_col %in% colnames(SummarizedExperiment::colData(se)))) {
                stop(sprintf("subject_col '%s' not found in colData(se)", subject_col))
            }
            subj_full <- as.character(SummarizedExperiment::colData(se)[, subject_col])
            # CRITICAL FIX: When colData doesn't have 'samples' column, use colnames(mat)
            # but strip "_q=" suffixes to match how sample_names were parsed in calculate_lm_interaction
            col_names_for_indexing <- SummarizedExperiment::colData(se)$samples
            if (is.null(col_names_for_indexing)) {
                # Strip "_q=....." suffixes from column names for proper indexing
                col_names_for_indexing <- sub("_q=.*", "", colnames(mat))
            }
            names(subj_full) <- col_names_for_indexing
            subject <- unname(subj_full[sample_names])
        } else if (paired) {
            coldata <- SummarizedExperiment::colData(se)
            coldata_cols <- colnames(coldata)
            
            # Look for 'paired_samples' or 'sample_base' columns created by map_metadata()
            subject_col_name <- NULL
            if ("paired_samples" %in% coldata_cols) {
                subject_col_name <- "paired_samples"
            } else if ("sample_base" %in% coldata_cols) {
                subject_col_name <- "sample_base"
            } else if (length(coldata_cols) >= 3) {
                subject_col_name <- coldata_cols[3]  # fallback
            }
            
            if (!is.null(subject_col_name)) {
                subject_ids <- as.character(coldata[, subject_col_name])
                # Map sample names by setting names to base names (without _q= suffix)
                sample_names_expanded <- sub("_q=.*", "", rownames(coldata))
                names(subject_ids) <- sample_names_expanded
                subject <- unname(subject_ids[sample_names])  # sample_names passed from calculate_lm_interaction
            } else {
                stop("paired = TRUE requires 'paired_samples' or 'sample_base' column in colData; supply subject_col explicitly or use map_metadata(...)")
            }
        } else {
            # Use sample names as cluster IDs (each sample is independent)
            subject <- sample_names
        }
        # Pass subject info to GEE helper with AR(1) correlation structure (default)
        return(.gee_interaction(df, q_vals, g, subject = subject, min_obs = min_obs, 
                                       corstr = corstr, bias_correction = bias_correction, weights = weights))
    }

    return(NULL)
}

# Helper functions for .calculate_lm_interaction()
# These internal functions decompose the main function logic into
# focused, testable components that each handle a single responsibility.

#' @title Validate Input Parameters for LM Interaction Testing
#'
#' @description
#' Internal helper that consolidates parameter validation for
#' \code{.calculate_lm_interaction()}. Checks argument types, values,
#' and inter-dependencies to ensure valid model fitting.
#'
#' @param method Character; modeling method (matched from user input)
#' @param pvalue Character; p-value type specification
#' @param corstr Character; correlation structure
#' @param regularization Character; dimensionality reduction method
#' @param multicorr Character; multi-q correction method
#' @param pcorr Character; legacy p-value correction method
#' @param storey Logical; whether to apply Storey correction
#' @param wy_randomizations Integer; number of permutations
#' @param paired Logical; whether design is paired
#' @param subject_col Character or NULL; subject column name
#' @param se SummarizedExperiment object
#' @param verbose Logical; print diagnostic messages
#'
#' @return List with validated and normalized parameters:
#'   \itemize{
#'     \item method: Validated method name
#'     \item pvalue: Validated p-value type
#'     \item corstr: Validated correlation structure
#'     \item regularization: Validated regularization method
#'     \item multicorr: Validated multicorr method
#'     \item pcorr: Validated legacy pcorr
#'     \item subject_col: Auto-detected or user-provided subject column
#'   }
#'

#' @noRd
.validate_lm_interaction_input <- function(
    method,
    pvalue,
    corstr,
    regularization,
    multicorr,
    pcorr,
    storey,
    wy_randomizations,
    paired,
    subject_col,
    se,
    verbose
) {
    # Validate storey parameter
    if (!is.logical(storey)) {
        stop("storey must be TRUE or FALSE", call. = FALSE)
    }

    # Validate wy_randomizations
    if (!is.numeric(wy_randomizations) || wy_randomizations < 1) {
        stop("wy_randomizations must be numeric and >= 1", call. = FALSE)
    }
    if (wy_randomizations < 100) {
        warning(
            "wy_randomizations < 100 may give unreliable p-values; ",
            "recommend >= 100",
            call. = FALSE
        )
    }

    # Auto-detect subject_col from colData if paired=TRUE and subject_col=NULL
    # Prioritize 'paired_samples' or 'sample_base' columns
    if (paired && is.null(subject_col)) {
        cd_colnames <- colnames(SummarizedExperiment::colData(se))

        # Check for paired_samples or sample_base columns
        if ("paired_samples" %in% cd_colnames) {
            subject_col <- "paired_samples"
            if (verbose) {
                message(
                    "[calculate_lm_interaction] paired=TRUE detected; ",
                    "auto-using subject_col='paired_samples'"
                )
            }
        } else if ("sample_base" %in% cd_colnames) {
            subject_col <- "sample_base"
            if (verbose) {
                message(
                    "[calculate_lm_interaction] paired=TRUE detected; ",
                    "auto-using subject_col='sample_base'"
                )
            }
        } else {
            # Error if paired=TRUE but no recognized pairing column found
            stop(
                "paired=TRUE requires either 'paired_samples' or ",
                "'sample_base' column in colData. Available columns: ",
                paste(cd_colnames, collapse = ", "),
                ". Ensure .calculate_diversity() or map_metadata() was ",
                "called with appropriate metadata.",
                call. = FALSE
            )
        }
    }

    return(list(
        method = method,
        pvalue = pvalue,
        corstr = corstr,
        regularization = regularization,
        multicorr = multicorr,
        pcorr = pcorr,
        subject_col = subject_col
    ))
}

#' @title Parse Sample Metadata from SummarizedExperiment
#'
#' @description
#' Internal helper that extracts sample names, q-values, and group
#' assignments from the diversity assay column names and colData.
#'
#' @param se SummarizedExperiment object
#' @param condition_col Character; colData column with group assignments
#' @param assay_name Character; name of diversity assay
#' @param verbose Logical; print diagnostic messages
#'
#' @return List containing:
#'   \itemize{
#'     \item sample_q: Full column names with q= values
#'     \item sample_names: Unique sample identifiers
#'     \item q_vals: Parsed q-value parameters
#'     \item group_vec: Group assignment for each observation
#'     \item has_q: Logical vector indicating cols with q=
#'   }
#'

#' @noRd
.parse_sample_metadata <- function(
    se,
    condition_col,
    assay_name,
    verbose
) {
    mat <- SummarizedExperiment::assay(se, assay_name)
    if (is.null(mat)) {
        stop(sprintf("Assay '%s' not found in SummarizedExperiment",
                     assay_name))
    }

    sample_q <- colnames(mat)
    if (is.null(sample_q) || length(sample_q) == 0) {
        stop("No column names found on diversity assay")
    }

    # Parse sample names and q values from column names like
    # 'Sample_q=0.01'
    sample_names <- sub("_q=.*", "", sample_q)
    has_q <- grepl("_q=", sample_q)
    if (!any(has_q)) {
        stop(
            "Could not parse q values; expected '_q=' in column names",
            call. = FALSE
        )
    }
    if (!all(has_q)) {
        stop(
            "Some column names are missing '_q='; ensure all diversity ",
            "columns include a q value",
            call. = FALSE
        )
    }
    q_vals <- as.numeric(sub(".*_q=", "", sample_q))

    # Determine group for each sample
    condition_in_coldata <- !is.null(condition_col) &&
        condition_col %in% colnames(SummarizedExperiment::colData(se))
    if (condition_in_coldata) {
        st <- as.character(
            SummarizedExperiment::colData(se)[, condition_col]
        )
        names(st) <- rownames(SummarizedExperiment::colData(se))
        # Index by the FULL column names (sample_q), not by sample_names
        group_vec <- unname(st[sample_q])
    } else {
        stop(
            "No sample grouping found: please supply `condition_col` ",
            "or map sample types into `colData(se)` before calling ",
            ".calculate_lm_interaction().",
            call. = FALSE
        )
    }

    if (verbose) {
        message(
            "[calculate_lm_interaction] parsed samples and groups: ",
            length(unique(sample_names)), " samples, ",
            length(unique(q_vals)), " q-values"
        )
    }

    return(list(
        sample_q = sample_q,
        sample_names = sample_names,
        q_vals = q_vals,
        group_vec = group_vec,
        has_q = has_q
    ))
}

#' @title Fit Linear Models for All Genes
#'
#' @description
#' Internal helper that orchestrates parallel or sequential fitting
#' of models to all genes in the diversity matrix. Consolidates the
#' fitting loop and result collection logic.
#'
#' @param mat Matrix; diversity assay data
#' @param se SummarizedExperiment object
#' @param metadata List; output from .parse_sample_metadata()
#' @param method Character; modeling method
#' @param pvalue Character; p-value type
#' @param subject_col Character or NULL; subject column
#' @param paired Logical; whether design is paired
#' @param min_obs Integer; minimum observations per gene
#' @param nthreads Integer; number of parallel threads
#' @param verbose Logical; print diagnostics
#' @param bias_correction Logical; apply KC bias correction (GEE)
#' @param regularization Character; dimensionality reduction method
#' @param corstr Character; correlation structure
#' @param adaptive_knots Logical; adaptive knot selection (GAM)
#'
#' @return Data.frame with fitted model results for all genes
#'

#' @noRd
.fit_all_genes <- function(
    mat,
    se,
    metadata,
    method,
    pvalue,
    subject_col,
    paired,
    min_obs,
    nthreads,
    verbose,
    bias_correction,
    regularization,
    corstr,
    adaptive_knots
) {
    suppress_lme4_warnings <- TRUE
    progress <- FALSE
    gene_weights <- NULL

    fit_one <- function(g, group_vec_override = NULL) {
        # Use override group_vec if provided (for permutation testing),
        # otherwise use outer scope
        gv <- if (!is.null(group_vec_override)) {
            group_vec_override
        } else {
            metadata$group_vec
        }

        .fit_one_interaction(
            g = g,
            se = se,
            mat = mat,
            q_vals = metadata$q_vals,
            sample_names = metadata$sample_names,
            group_vec = gv,
            method = method,
            pvalue = pvalue,
            subject_col = subject_col,
            paired = paired,
            min_obs = min_obs,
            verbose = verbose,
            suppress_lme4_warnings = suppress_lme4_warnings,
            progress = progress,
            bias_correction = bias_correction,
            regularization = regularization,
            corstr = corstr,
            adaptive_knots = adaptive_knots,
            weights = gene_weights
        )
    }

    if (nthreads > 1) {
        res_list <- .bplapply(rownames(mat), fit_one,
                                      nthreads = nthreads)
    } else {
        res_list <- lapply(rownames(mat), fit_one)
    }
    all_results <- Filter(Negate(is.null), res_list)

    if (length(all_results) == 0) {
        return(data.frame())
    }
    res <- do.call(rbind, all_results)

    # VALIDATION: Ensure critical columns exist after rbind
    if (nrow(res) == 0) {
        warning(
            "[calculate_lm_interaction] No genes analyzed (all filtered out)",
            call. = FALSE
        )
        return(res)
    }

    critical_cols <- c("p_interaction", "gene")
    missing_cols <- setdiff(critical_cols, colnames(res))
    if (length(missing_cols) > 0) {
        stop(
            "[calculate_lm_interaction] CRITICAL: Missing columns in ",
            "results for ", method, " method: ",
            paste(missing_cols, collapse = ", "),
            "\nAvailable columns: ",
            paste(colnames(res), collapse = ", "),
            call. = FALSE
        )
    }

    # Ensure Shapiro-Wilk columns exist for methods that add them
    if (method %in% c("gam", "gee")) {
        if (!"shapiro_p_value" %in% colnames(res)) {
            res$shapiro_p_value <- NA_real_
        }
        if (!"residuals_normal" %in% colnames(res)) {
            res$residuals_normal <- NA
        }
        if (!"n_residuals_tested" %in% colnames(res)) {
            res$n_residuals_tested <- NA_integer_
        }
    }

    # Ensure ci_weighted column exists (Phase 1 tracking)
    if (!"ci_weighted" %in% colnames(res)) {
        res$ci_weighted <- NA  # Fallback
        if (verbose) {
            warning(
                "[calculate_lm_interaction] ci_weighted column was ",
                "missing; added as NAs. This suggests a method helper ",
                "did not properly set ci_weighted.",
                call. = FALSE
            )
        }
    }

    return(res)
}

#' @title Adjust P-Values for Multiple Q-Values
#'
#' @description
#' Internal router function that applies the specified primary multi-q
#' p-value correction method (Hochberg, Westfall-Young, or
#' Benjamini-Yekutieli), then optionally applies Storey adaptive FDR
#' enhancement.
#'
#' @param p_values Numeric vector; raw p-values to adjust
#' @param multicorr Character; primary correction method
#' @param wy_randomizations Integer; number of permutations (WY only)
#' @param fit_one_fn Function; function to refit models (WY only)
#' @param metadata List; output from .parse_sample_metadata()
#' @param mat Matrix; diversity assay data (WY only)
#' @param rownames_mat Character; rownames of matrix (WY only)
#' @param se SummarizedExperiment object (WY only)
#' @param assay_name Character; assay name (WY only)
#' @param method Character; modeling method (WY only)
#' @param pvalue Character; p-value type (WY only)
#' @param subject_col Character or NULL; subject column (WY only)
#' @param paired Logical; paired design (WY only)
#' @param min_obs Integer; min observations (WY only)
#' @param nthreads Integer; parallel threads (WY only)
#' @param verbose Logical; print diagnostics
#' @param bias_correction Logical; KC bias correction (WY/GEE)
#' @param regularization Character; dimensionality reduction (WY/FPCA)
#' @param corstr Character; correlation structure (WY/GEE)
#' @param adaptive_knots Logical; adaptive knots (WY/GAM)
#' @param storey Logical; apply Storey after primary correction
#'
#' @return Adjusted p-values vector
#'

#' @noRd
.adjust_pvalues_multicorr <- function(
    p_values,
    multicorr,
    wy_randomizations,
    fit_one_fn = NULL,
    metadata = NULL,
    mat = NULL,
    rownames_mat = NULL,
    se = NULL,
    assay_name = "diversity",
    method = NULL,
    pvalue = NULL,
    subject_col = NULL,
    paired = FALSE,
    min_obs = 10,
    nthreads = 1,
    verbose = FALSE,
    bias_correction = TRUE,
    regularization = NULL,
    corstr = "ar1",
    adaptive_knots = TRUE,
    storey = FALSE
) {
    if (multicorr == "hochberg") {
        adj_p <- .hochberg_stepup(p_values)
        if (verbose) {
            message(
                "[calculate_lm_interaction] Applied Hochberg stepup ",
                "adjustment for multi-q correlation"
            )
        }
    } else if (multicorr == "westfall-young") {
        if (verbose) {
            message(
                "[calculate_lm_interaction] Computing true ",
                "Westfall-Young via ", wy_randomizations,
                " permutations (may be slow)..."
            )
        }

        # Save original group vector for safe restoration
        group_vec_orig <- metadata$group_vec

        # Run Westfall-Young permutation
        perm_result <- .westfall_young_permutation(
            n_genes = length(p_values),
            wy_randomizations = wy_randomizations,
            permute_fn = function() {
                # Shuffle group labels separately within each q-level
                q_unique <- unique(metadata$q_vals)
                perm_assignment <- group_vec_orig
                for (q_val in q_unique) {
                    q_idx <- which(metadata$q_vals == q_val)
                    perm_assignment[q_idx] <-
                        sample(group_vec_orig[q_idx])
                }
                return(perm_assignment)
            },
            refit_fn = function(perm_assignment) {
                # Refit all genes with permuted group assignment
                perm_pvalues <- numeric(length(p_values))
                for (g_idx in seq_along(rownames_mat)) {
                    gene_name <- rownames_mat[g_idx]
                    tryCatch({
                        gene_result <- fit_one_fn(
                            gene_name,
                            group_vec_override = perm_assignment
                        )
                        if (!is.null(gene_result) &&
                            !is.na(gene_result$p_interaction)) {
                            perm_pvalues[g_idx] <-
                                gene_result$p_interaction
                        }
                    }, error = function(e) { NULL })
                }
                return(perm_pvalues)
            },
            nthreads = nthreads,
            verbose = verbose
        )

        # Adjust p-values based on permutation distribution
        adj_p <- vapply(p_values, function(p_obs) {
            pmin(1.0,
                 (sum(perm_result$perm_minima <= p_obs) + 1) /
                     (wy_randomizations + 1))
        }, FUN.VALUE = numeric(1))

        if (verbose) {
            message(
                "[calculate_lm_interaction] Applied true ",
                "Westfall-Young (permutation) adjustment"
            )
        }
    } else if (multicorr == "benjamini-yekutieli") {
        adj_p <- .benjamini_yekutieli(p_values)
        if (verbose) {
            message(
                "[calculate_lm_interaction] Applied ",
                "Benjamini-Yekutieli adjustment for dependent tests"
            )
        }
    } else {
        stop("Unknown multicorr method: ", multicorr, call. = FALSE)
    }

    # Apply optional Storey adaptive FDR enhancement layer
    if (storey) {
        if (requireNamespace("fdrtool", quietly = TRUE)) {
            tryCatch({
                adj_p <- .compute_storey_qvalues(adj_p)
                if (verbose) {
                    message(
                        "[calculate_lm_interaction] Applied Storey ",
                        "adaptive FDR pi0 correction to ",
                        multicorr, " p-values"
                    )
                }
            }, error = function(e) {
                if (verbose) {
                    message(
                        "[calculate_lm_interaction] Storey adjustment ",
                        "failed: ", conditionMessage(e)
                    )
                }
            })
        } else if (verbose) {
            message(
                "[calculate_lm_interaction] fdrtool package not ",
                "available for Storey (install with: ",
                "install.packages('fdrtool'))"
            )
        }
    }

    return(adj_p)
}

#' @title Map Gene Identifiers to Annotations
#'
#' @description
#' Internal helper that maps gene rownames to gene_id and gene_name
#' columns using rowData from the SummarizedExperiment. Ensures
#' consistent gene annotation across downstream analyses.
#'
#' @param res Data.frame; results with gene column (rownames)
#' @param se SummarizedExperiment object
#' @param verbose Logical; print diagnostic messages
#'
#' @return Modified data.frame with gene_id and gene_name columns added
#'

#' @noRd
.map_gene_annotations <- function(
    res,
    se,
    verbose
) {
    rd <- SummarizedExperiment::rowData(se)

    # Look for gene_name column from calculate_diversity or build_se
    gene_name_col <- if ("gene_name" %in% colnames(rd)) {
        "gene_name"
    } else {
        NULL
    }

    if (verbose) {
        message(
            "[calculate_lm_interaction] Gene annotations: ",
            paste(colnames(rd), collapse = ", ")
        )
    }

    if (!is.null(gene_name_col)) {
        # Determine gene_id column if it exists
        id_col <- if ("genes" %in% colnames(rd)) {
            "genes"
        } else if ("gene_id" %in% colnames(rd)) {
            "gene_id"
        } else {
            NA  # rownames will be used as ID
        }

        # Build lookup tables: rowname -> gene_id and rowname -> gene_name
        if (is.na(id_col)) {
            # rownames ARE the gene IDs
            rowname_to_id <- setNames(
                as.character(rownames(rd)),
                as.character(rownames(rd))
            )
        } else {
            # gene IDs are in a column
            rowname_to_id <- setNames(
                as.character(rd[[id_col]]),
                as.character(rownames(rd))
            )
        }

        rowname_to_name <- setNames(
            as.character(rd[[gene_name_col]]),
            as.character(rownames(rd))
        )

        # Vectorized lookup: map res$gene to gene_id and gene_name
        res$gene_id <- unname(rowname_to_id[as.character(res$gene)])
        res$gene_name <- unname(rowname_to_name[as.character(res$gene)])

        # For any unmapped genes, use gene column as fallback
        unmapped_idx <- is.na(res$gene_name)
        n_mapped <- sum(!unmapped_idx)
        n_unmapped <- sum(unmapped_idx)

        if (any(unmapped_idx)) {
            res$gene_name[unmapped_idx] <- res$gene[unmapped_idx]
        }

        if (verbose && n_unmapped > 0) {
            message(
                "[calculate_lm_interaction] Gene mapping: ",
                n_mapped, " mapped, ", n_unmapped,
                " used ID as fallback"
            )
        }
    } else if (verbose) {
        message(
            "[calculate_lm_interaction] gene_name column not found ",
            "in rowData - using gene ID as fallback"
        )
    }

    # Ensure gene_name column is always present and populated
    if (is.null(res$gene_name) || !"gene_name" %in% colnames(res)) {
        res$gene_name <- res$gene
    }

    # Ensure gene_id column is always present and populated
    if (is.null(res$gene_id) || !"gene_id" %in% colnames(res)) {
        res$gene_id <- res$gene
    }

    return(res)
}

#' @title Assemble Model Metadata for Diagnostics
#'
#' @description
#' Internal helper that builds the comprehensive metadata list returned
#' when \code{return_model_data = TRUE}. Contains method info, q-values,
#' per-group statistics, and test configuration for downstream visualization.
#'
#' @param se SummarizedExperiment object
#' @param res Data.frame; fitted model results
#' @param mat Matrix; diversity assay data
#' @param metadata List; output from .parse_sample_metadata()
#' @param method Character; modeling method name
#' @param pvalue Character; p-value type
#' @param multicorr Character; multi-q correction method
#' @param assay_name Character; assay name
#' @param bias_correction Logical; KC bias correction setting
#' @param regularization Character; dimensionality reduction method
#' @param corstr Character; correlation structure
#' @param adaptive_knots Logical; adaptive knot selection setting
#'
#' @return List with comprehensive model metadata:
#'   \itemize{
#'     \item method: Modeling method used
#'     \item n_genes: Number of genes analyzed
#'     \item n_q_values: Number of q parameters
#'     \item q_values: The specific q-value vector
#'     \item sample_names: Unique sample identifiers
#'     \item group_levels: Group factor levels
#'     \item per_group_statistics: Summary statistics by group
#'     \item test_configuration: Complete test settings
#'     \item genes_analyzed: Vector of gene identifiers
#'     \item call_time: Timestamp of analysis
#'     \item notes: Usage information
#'   }
#'

#' @noRd
.assemble_model_metadata <- function(
    se,
    res,
    mat,
    metadata,
    method,
    pvalue,
    multicorr,
    assay_name = "diversity",
    bias_correction = TRUE,
    regularization = "pca",
    corstr = "ar1",
    adaptive_knots = TRUE
) {
    # Extract per-group statistics from SE
    per_group_stats <- list()

    for (gr in unique(metadata$group_vec)) {
        gr_idx <- which(metadata$group_vec == gr)
        gr_mat <- mat[, gr_idx, drop = FALSE]

        per_group_stats[[gr]] <- list(
            group = gr,
            n_samples = length(unique(metadata$sample_names[gr_idx])),
            n_observations = ncol(gr_mat),
            entropy_mean = mean(as.numeric(gr_mat), na.rm = TRUE),
            entropy_sd = sd(as.numeric(gr_mat), na.rm = TRUE),
            entropy_min = min(as.numeric(gr_mat), na.rm = TRUE),
            entropy_max = max(as.numeric(gr_mat), na.rm = TRUE),
            entropy_median = median(as.numeric(gr_mat), na.rm = TRUE),
            n_na = sum(is.na(gr_mat))
        )
    }

    model_data <- list(
        method = method,
        n_genes = nrow(res),
        n_q_values = length(unique(metadata$q_vals)),
        q_values = sort(unique(metadata$q_vals)),
        sample_names = unique(metadata$sample_names),
        group_levels = levels(factor(metadata$group_vec)),
        per_group_statistics = per_group_stats,
        test_configuration = list(
            method = method,
            pvalue_method = pvalue,
            multicorr = multicorr,
            bias_correction = bias_correction,
            regularization = regularization,
            corstr = corstr,
            adaptive_knots = adaptive_knots
        ),
        genes_analyzed = res$gene,
        call_time = Sys.time(),
        notes = paste(
            "Use this model_data with plotting functions to ",
            "visualize model fits and diagnostics. See ",
            "per_group_statistics for condition-specific ",
            "entropy summaries."
        )
    )

    return(model_data)
}






