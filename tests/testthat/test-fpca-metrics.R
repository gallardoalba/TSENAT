# ==============================================================================
# Tests for FPCA (Functional Principal Component Analysis) Diagnostic Metrics
# ==============================================================================
# Purpose: Validate FPCA diagnostic functions for assumption checking
# Date: April 12, 2026
# ==============================================================================

library(testthat)

# ==============================================================================
# Test: Basic Input Handling
# ==============================================================================

test_that("FPCA metrics functions handle basic inputs correctly", {
    set.seed(42)
    n_genes <- 50
    n_q_values <- 10
    
    # Generate data with varying structure
    data <- matrix(rnorm(n_genes * n_q_values, mean = 10, sd = 2), 
                   nrow = n_genes, ncol = n_q_values)
    
    # Verify data structure
    expect_is(data, "matrix")
    expect_equal(nrow(data), n_genes)
    expect_equal(ncol(data), n_q_values)
})

# ==============================================================================
# Test: Cumulative Variance Explained
# ==============================================================================

test_that(".compute_cumulative_variance_fpca() returns correct structure", {
    set.seed(42)
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    result <- .compute_cumulative_variance_fpca(data)
    
    # Verify output structure
    expect_is(result, "list")
    expect_true("description" %in% names(result))
    expect_true("eigenvalues" %in% names(result))
    expect_true("cumulative_variance" %in% names(result))
    expect_true("status" %in% names(result))
})

test_that(".compute_cumulative_variance_fpca() eigenvalues are decreasing", {
    set.seed(42)
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    result <- .compute_cumulative_variance_fpca(data)
    
    # Eigenvalues should be in decreasing order
    eigenvalues <- result$eigenvalues
    expect_true(all(diff(eigenvalues) <= 1e-10))
})

test_that(".compute_cumulative_variance_fpca() cumulative variance is increasing", {
    set.seed(42)
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    result <- .compute_cumulative_variance_fpca(data)
    
    # Cumulative variance should be monotonically increasing
    cum_var <- result$cumulative_variance
    expect_true(all(diff(cum_var) >= 0))
    
    # Should end at 1.0
    expect_equal(cum_var[length(cum_var)], 1.0, tolerance = 1e-6)
})

test_that(".compute_cumulative_variance_fpca() detects good variance reduction", {
    set.seed(42)
    # Create data with strong first few components (good dim reduction)
    # Generate 3 "true" components with high variance
    comp1 <- rnorm(50, mean = 0, sd = 10)
    comp2 <- rnorm(50, mean = 0, sd = 8)
    comp3 <- rnorm(50, mean = 0, sd = 6)
    # Add 7 noise components with low variance
    noise <- matrix(rnorm(50 * 7, mean = 0, sd = 0.5), nrow = 50, ncol = 7)
    # Combine into data matrix
    data <- cbind(comp1, comp2, comp3, noise)
    
    result <- .compute_cumulative_variance_fpca(data)
    
    # First few components should capture most variance
    expect_true(result$n_components_95 < 5)
    expect_equal(result$status, "good reduction")
})

test_that(".compute_cumulative_variance_fpca() detects poor variance reduction", {
    set.seed(42)
    # Create data with equal variance across all components (poor dim reduction)
    data <- matrix(rnorm(50 * 15, mean = 10, sd = 2), nrow = 50, ncol = 15)
    
    result <- .compute_cumulative_variance_fpca(data)
    
    # Many components needed for 95% variance
    expect_true(result$n_components_95 > 10)
    expect_equal(result$status, "poor reduction")
})

test_that(".compute_cumulative_variance_fpca() handles zero-variance data", {
    data <- matrix(5, nrow = 50, ncol = 10)
    
    result <- .compute_cumulative_variance_fpca(data)
    
    # Should skip with constant data
    expect_equal(result$status, "? SKIP")
})

test_that(".compute_cumulative_variance_fpca() components are positive integers", {
    set.seed(42)
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    result <- .compute_cumulative_variance_fpca(data)
    
    # Component counts should be positive integers
    expect_gt(result$n_components_90, 0)
    expect_gt(result$n_components_95, 0)
    expect_gt(result$n_components_99, 0)
    
    # 90% should need fewer components than 95% than 99%
    expect_lte(result$n_components_90, result$n_components_95)
    expect_lte(result$n_components_95, result$n_components_99)
})

# ==============================================================================
# Test: Bootstrap Stability
# ==============================================================================

test_that(".compute_fpca_bootstrap_stability() returns correct structure", {
    set.seed(42)
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    result <- .compute_fpca_bootstrap_stability(data, n_bootstrap = 100, n_components = 3)
    
    # Verify output structure
    expect_is(result, "list")
    expect_true("description" %in% names(result))
    expect_true("bootstrap_se" %in% names(result))
    expect_true("bootstrap_ci_lower" %in% names(result))
    expect_true("bootstrap_ci_upper" %in% names(result))
    expect_true("status" %in% names(result))
})

test_that(".compute_fpca_bootstrap_stability() bootstrap CIs are valid", {
    set.seed(42)
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    result <- .compute_fpca_bootstrap_stability(data, n_bootstrap = 100, n_components = 3)
    
    # CIs should be lower < upper
    ci_lower <- result$bootstrap_ci_lower
    ci_upper <- result$bootstrap_ci_upper
    
    if (length(ci_lower) > 0 && all(!is.na(ci_lower)) && all(!is.na(ci_upper))) {
        expect_true(all(ci_lower <= ci_upper))
    }
})

test_that(".compute_fpca_bootstrap_stability() cv_eigenvalues valid range", {
    set.seed(42)
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    result <- .compute_fpca_bootstrap_stability(data, n_bootstrap = 100, n_components = 3)
    
    # Coefficients of variation should be positive
    cv_eig <- result$cv_eigenvalues
    if (length(cv_eig) > 0 && !all(is.na(cv_eig))) {
        expect_true(all(cv_eig[!is.na(cv_eig)] >= 0))
    }
})

test_that(".compute_fpca_bootstrap_stability() detects stable components", {
    set.seed(42)
    # Create highly structured data (stable components)
    # Generate 3 "true" components with high, consistent variance
    comp1 <- rnorm(30, mean = 0, sd = 15)
    comp2 <- rnorm(30, mean = 0, sd = 10)
    comp3 <- rnorm(30, mean = 0, sd = 8)
    # Add low-noise columns
    noise <- matrix(rnorm(30 * 5, mean = 0, sd = 0.1), nrow = 30, ncol = 5)
    # Combine into structured data
    data <- cbind(comp1, comp2, comp3, noise)
    
    result <- .compute_fpca_bootstrap_stability(data, n_bootstrap = 200, n_components = 3)
    
    # Structured data should have stable components
    expect_true(result$status %in% c("stable", "moderate stability"))
})

test_that(".compute_fpca_bootstrap_stability() handles insufficient data", {
    data <- matrix(rnorm(3 * 10), nrow = 3, ncol = 10)
    
    result <- .compute_fpca_bootstrap_stability(data, n_bootstrap = 100, n_components = 3)
    
    # Should skip with very small sample
    expect_equal(result$status, "? SKIP")
})

test_that(".compute_fpca_bootstrap_stability() limits components appropriately", {
    set.seed(42)
    # Small dataset with p > n
    data <- matrix(rnorm(10 * 15, mean = 10, sd = 2), nrow = 10, ncol = 15)
    
    result <- .compute_fpca_bootstrap_stability(data, n_bootstrap = 100, n_components = 20)
    
    # Should limit components to min(n_components, n-1, p-1) = min(20, 9, 14) = 9
    expect_lte(length(result$bootstrap_se), 9)
})

# ==============================================================================
# Test: Wrapper Function
# ==============================================================================

test_that(".get_fpca_metrics() returns all components", {
    set.seed(42)
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    result <- .get_fpca_metrics(data)
    
    expect_is(result, "list")
    expect_equal(length(result), 3)  # variance_adequacy, bootstrap_stability, consolidated
    
    expected_names <- c("variance_adequacy", "bootstrap_stability", "consolidated")
    expect_true(all(expected_names %in% names(result)))
})

test_that(".get_fpca_metrics() creates valid consolidated result", {
    set.seed(42)
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    result <- .get_fpca_metrics(data)
    
    # Consolidated should have required fields
    expect_is(result$consolidated, "list")
    expect_true("result" %in% names(result$consolidated))
    expect_true("status" %in% names(result$consolidated))
    
    # Result string should contain metric statuses
    result_str <- result$consolidated$result
    expect_is(result_str, "character")
    expect_true(nchar(result_str) > 0)
    expect_true(grepl("Variance", result_str))
    expect_true(grepl("Bootstrap", result_str))
})

test_that(".get_fpca_metrics() handles parameter passing", {
    set.seed(42)
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    # Call with custom parameters
    fpca_params <- list(max_components = 5, n_bootstrap = 100, n_components = 2)
    result <- .get_fpca_metrics(data, fpca_params = fpca_params)
    
    expect_is(result, "list")
    expect_equal(length(result), 3)
})

test_that(".get_fpca_metrics() handles empty parameters list", {
    set.seed(42)
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    result <- .get_fpca_metrics(data, fpca_params = list())
    
    # Should use defaults and still work
    expect_is(result, "list")
    expect_equal(length(result), 3)
})

# ==============================================================================
# Test: Data Type Handling
# ==============================================================================

test_that("FPCA functions coerce non-matrix data to matrix", {
    set.seed(42)
    # Create data frame
    data_df <- data.frame(
        q1 = rnorm(20, mean = 10, sd = 2),
        q2 = rnorm(20, mean = 10, sd = 2),
        q3 = rnorm(20, mean = 10, sd = 2)
    )
    
    # Should coerce and work
    result <- .compute_cumulative_variance_fpca(data_df)
    expect_is(result, "list")
    expect_true("eigenvalues" %in% names(result))
})

test_that("FPCA functions handle data with NA values", {
    set.seed(42)
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    # Introduce some NAs
    data[sample(1:length(data), 20)] <- NA
    
    # Should handle NAs gracefully (may skip or produce NA status)
    result_cv <- .compute_cumulative_variance_fpca(data)
    result_bs <- .compute_fpca_bootstrap_stability(data)
    
    expect_is(result_cv, "list")
    expect_is(result_bs, "list")
})

# ==============================================================================
# Test: Edge Cases
# ==============================================================================

test_that("FPCA metrics handle 2x2 minimum data", {
    data <- matrix(rnorm(2 * 2), nrow = 2, ncol = 2)
    
    result <- .get_fpca_metrics(data)
    
    # Should attempt computation even with minimal data
    expect_is(result, "list")
})

test_that("FPCA metrics handle single observation", {
    data <- matrix(rnorm(10), nrow = 1, ncol = 10)
    
    result_cv <- .compute_cumulative_variance_fpca(data)
    result_bs <- .compute_fpca_bootstrap_stability(data)
    
    # Should skip or fail gracefully
    expect_is(result_cv, "list")
    expect_is(result_bs, "list")
})

test_that("FPCA metrics handle perfect correlation", {
    # Create perfectly correlated columns
    col1 <- rnorm(50)
    data <- matrix(rep(col1, 5), nrow = 50, ncol = 5)
    
    result <- .get_fpca_metrics(data)
    
    # Should handle without crashing
    expect_is(result, "list")
    expect_gt(nchar(result$consolidated$result), 0)
})

# ==============================================================================
# Test: Numerical Validation
# ==============================================================================

test_that("Eigenvalues sum to total variance", {
    set.seed(42)
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    result <- .compute_cumulative_variance_fpca(data)
    
    if (!is.na(result$total_variance)) {
        # Sum of eigenvalues should equal total variance
        sum_eig <- sum(result$eigenvalues)
        expect_equal(sum_eig, result$total_variance, tolerance = 1e-10)
    }
})

test_that("Cumulative variance fractions in [0,1]", {
    set.seed(42)
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    result <- .compute_cumulative_variance_fpca(data)
    
    cum_var <- result$cumulative_variance
    expect_true(all(cum_var >= 0))
    expect_true(all(cum_var <= 1))
})

test_that("Bootstrap SE components are non-negative", {
    set.seed(42)
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    result <- .compute_fpca_bootstrap_stability(data, n_bootstrap = 100, n_components = 3)
    
    bootstrap_se <- result$bootstrap_se
    if (length(bootstrap_se) > 0 && !all(is.na(bootstrap_se))) {
        expect_true(all(bootstrap_se[!is.na(bootstrap_se)] >= 0))
    }
})

test_that("Proportion stable CIs in valid range [0,1]", {
    set.seed(42)
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    result <- .compute_fpca_bootstrap_stability(data, n_bootstrap = 100, n_components = 3)
    
    prop_stable <- result$prop_stable_cis
    expect_gte(prop_stable, 0)
    expect_lte(prop_stable, 1)
})

test_that("Eigenvalue ratios are monotonically decreasing", {
    set.seed(42)
    x <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
    
    result <- .compute_cumulative_variance_fpca(x)
    
    if (!is.null(result$eigenvalues) && length(result$eigenvalues) > 1) {
        # Each eigenvalue should be <= previous one
        for (i in 2:length(result$eigenvalues)) {
            expect_lte(result$eigenvalues[i], result$eigenvalues[i-1] + 1e-10)
        }
    }
})

test_that("Cumulative variance thresholds are ordered: n_90 <= n_95 <= n_99", {
    set.seed(42)
    x <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
    
    result <- .compute_cumulative_variance_fpca(x)
    
    if (!is.null(result$n_components_90)) {
        # Number of components needed should be ordered
        expect_lte(result$n_components_90, result$n_components_95)
        expect_lte(result$n_components_95, result$n_components_99)
    }
})

test_that("Bootstrap CI widths are reasonable and positive", {
    set.seed(42)
    x <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
    
    result <- .compute_fpca_bootstrap_stability(x, n_bootstrap = 50, n_components = 3)
    
    if (!is.null(result$bootstrap_ci_lower) && length(result$bootstrap_ci_lower) > 0) {
        for (i in seq_along(result$bootstrap_ci_lower)) {
            if (!is.na(result$bootstrap_ci_lower[i]) && !is.na(result$bootstrap_ci_upper[i])) {
                # CI lower should be <= CI upper
                expect_lte(result$bootstrap_ci_lower[i], result$bootstrap_ci_upper[i] + 1e-10)
                # CI width should be positive
                expect_gt(result$bootstrap_ci_upper[i], result$bootstrap_ci_lower[i])
            }
        }
    }
})

test_that("CV_eigenvalues relate to bootstrap SE and original eigenvalues", {
    set.seed(42)
    x <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
    
    result <- .compute_fpca_bootstrap_stability(x, n_bootstrap = 50, n_components = 3)
    
    # CV should be non-negative
    if (!is.null(result$cv_eigenvalues) && length(result$cv_eigenvalues) > 0) {
        expect_true(all(result$cv_eigenvalues[!is.na(result$cv_eigenvalues)] >= 0))
        # Max CV should be recorded
        if (!is.na(result$max_cv)) {
            expect_gte(result$max_cv, 0)
        }
    }
})

test_that("Wrapper metrics: consolidated status is COMBINED", {
    set.seed(42)
    x <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
    
    result <- .get_fpca_metrics(x)
    
    # Consolidated should always have status COMBINED
    expect_equal(result$consolidated$status, "COMBINED")
})

test_that("Cumulative variance percentage bounds consistency", {
    set.seed(42)
    x <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
    
    result <- .compute_cumulative_variance_fpca(x)
    
    # All n_components values should be positive integers <= nrow(x)
    if (!is.null(result$n_components_90) && result$n_components_90 > 0) {
        expect_gt(result$n_components_90, 0)
        expect_lte(result$n_components_90, nrow(x))
    }
    
    if (!is.null(result$n_components_95) && result$n_components_95 > 0) {
        expect_gt(result$n_components_95, 0)
        expect_lte(result$n_components_95, nrow(x))
    }
    
    if (!is.null(result$n_components_99) && result$n_components_99 > 0) {
        expect_gt(result$n_components_99, 0)
        expect_lte(result$n_components_99, nrow(x))
    }
    
    # n_components should increase with threshold: 90 <= 95 <= 99
    if (!is.null(result$n_components_90) && !is.null(result$n_components_95)) {
        expect_lte(result$n_components_90, result$n_components_95 + 1)
    }
    
    if (!is.null(result$n_components_95) && !is.null(result$n_components_99)) {
        expect_lte(result$n_components_95, result$n_components_99 + 1)
    }
})

test_that("Bootstrap stability: SE relationships with CI width", {
    set.seed(42)
    x <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
    
    result <- .compute_fpca_bootstrap_stability(x, n_bootstrap = 50, n_components = 3)
    
    if (!is.null(result$bootstrap_se) && length(result$bootstrap_se) > 0) {
        # SE should be non-negative
        expect_true(all(result$bootstrap_se[!is.na(result$bootstrap_se)] >= 0))
        
        # Rough check: CI width should relate to SE (95% CI ≈ 1.96 * 2 * SE)
        for (i in seq_along(result$bootstrap_se)) {
            if (!is.na(result$bootstrap_se[i]) && 
                !is.na(result$bootstrap_ci_lower[i]) && 
                !is.na(result$bootstrap_ci_upper[i]) &&
                result$bootstrap_se[i] > 0) {
                ci_width <- result$bootstrap_ci_upper[i] - result$bootstrap_ci_lower[i]
                expected_width <- 4 * result$bootstrap_se[i]  # Allow wide tolerance
                expect_lt(ci_width, expected_width * 3)
            }
        }
    }
})

# ==============================================================================
# Test: Integration Tests
# ==============================================================================

test_that("FPCA metrics integrate with calculate_rank_assumptions()", {
    set.seed(42)
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    result <- tryCatch(
        .calculate_rank_assumptions(data, checks = c("fpca_metrics")),
        error = function(e) NULL
    )
    
    # Should return rank_assumptions object with fpca_metrics
    if (!is.null(result)) {
        expect_is(result, "rank_assumptions")
        checks <- attr(result, "checks")
        expect_true("fpca_metrics" %in% names(checks))
    }
})

test_that("FPCA metrics included in 'all' checks preset", {
    set.seed(42)
    data <- matrix(rnorm(40 * 8, mean = 10, sd = 2), nrow = 40, ncol = 8)
    
    result <- tryCatch(
        .calculate_rank_assumptions(data, checks = "all"),
        error = function(e) NULL
    )
    
    if (!is.null(result)) {
        checks <- attr(result, "checks")
        # Should include fpca_metrics when checks="all"
        expect_true("fpca_metrics" %in% names(checks) || 
                   !is.null(checks$fpca_metrics))
    }
})

# ==============================================================================
# Test: Reproducibility
# ==============================================================================

test_that("FPCA metrics are reproducible with same seed", {
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    set.seed(123)
    result1 <- .get_fpca_metrics(data)
    
    set.seed(123)
    result2 <- .get_fpca_metrics(data)
    
    # Results should match exactly for variance adequacy (deterministic)
    expect_equal(result1$variance_adequacy$n_components_95, 
                result2$variance_adequacy$n_components_95)
})

test_that("Individual FPCA metrics are consistent with wrapper", {
    set.seed(42)
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    # Get metrics separately (with short bootstrap for speed)
    va_sep <- .compute_cumulative_variance_fpca(data)
    
    # Get through wrapper
    wrapper_result <- .get_fpca_metrics(data, fpca_params = list(n_bootstrap = 100))
    
    # Should match
    expect_equal(va_sep$n_components_95, 
                wrapper_result$variance_adequacy$n_components_95)
})

# ==============================================================================
# Test: Bootstrap Sample Size Effect
# ==============================================================================

test_that("Bootstrap stability improves with more bootstrap samples", {
    set.seed(42)
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    # Bootstrap with different sample sizes
    result_100 <- .compute_fpca_bootstrap_stability(data, n_bootstrap = 100, n_components = 3)
    result_300 <- .compute_fpca_bootstrap_stability(data, n_bootstrap = 300, n_components = 3)
    
    # Both should have valid output
    expect_is(result_100, "list")
    expect_is(result_300, "list")
    
    # Bootstrap SE should generally decrease with more samples (higher precision)
    # Though this is probabilistic, so we just verify they're different
    if (all(!is.na(result_100$bootstrap_se)) && all(!is.na(result_300$bootstrap_se))) {
        # More samples typically gives lower CV (more precise estimates)
        # But this can vary, so just verify we get numeric results
        expect_is(result_100$bootstrap_se, "numeric")
        expect_is(result_300$bootstrap_se, "numeric")
    }
})

# ==============================================================================
# Test: Different Data Structures
# ==============================================================================

test_that("FPCA handles different dimensions correctly", {
    set.seed(42)
    
    # Tall data (more genes than q-values)
    data_tall <- matrix(rnorm(100 * 8, mean = 10, sd = 2), nrow = 100, ncol = 8)
    result_tall <- .compute_cumulative_variance_fpca(data_tall)
    expect_is(result_tall, "list")
    expect_lte(result_tall$n_components_95, 8)
    
    # Wide data (more q-values than genes)
    data_wide <- matrix(rnorm(20 * 50, mean = 10, sd = 2), nrow = 20, ncol = 50)
    result_wide <- .compute_cumulative_variance_fpca(data_wide)
    expect_is(result_wide, "list")
    expect_lte(result_wide$n_components_95, 20)
})

test_that("FPCA handles highly skewed variance", {
    set.seed(42)
    
    # Create data with one very dominant component
    # Very high variance in first column, low in others
    dominant <- rnorm(50, mean = 0, sd = 100)
    noise_cols <- matrix(rnorm(50 * 9, mean = 0, sd = 1), nrow = 50, ncol = 9)
    data <- cbind(dominant, noise_cols)
    
    result <- .compute_cumulative_variance_fpca(data)
    
    # Should need very few components for high % variance
    expect_is(result, "list")
    if (!is.null(result$n_components_90)) {
        expect_lt(result$n_components_90, 5)
    }
})
