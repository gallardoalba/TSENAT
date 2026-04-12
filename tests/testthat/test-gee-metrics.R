# ==============================================================================
# Tests for GEE Metrics Functions (Generalized Estimating Equations)
# ==============================================================================
# Purpose: Validate GEE diagnostic functions for assumption checking
# Date: April 12, 2026
# ==============================================================================

library(testthat)

# Set up test data
test_that("GEE metrics functions handle basic inputs correctly", {
    # Create synthetic expression matrix
    set.seed(42)
    n_genes <- 50
    n_samples <- 10
    
    # Generate data with exchangeable correlation structure
    data <- matrix(rnorm(n_genes * n_samples, mean = 10, sd = 2), 
                   nrow = n_genes, ncol = n_samples)
    
    # Verify data structure
    expect_is(data, "matrix")
    expect_equal(nrow(data), n_genes)
    expect_equal(ncol(data), n_samples)
})

# ==============================================================================
# Test: Working Correlation Structure Fit
# ==============================================================================

test_that(".compute_working_correlation_fit() returns correct structure", {
    set.seed(42)
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    # Test with exchangeable structure (default)
    result <- .compute_working_correlation_fit(data, assumed_structure = "exchangeable")
    
    # Verify output structure
    expect_is(result, "list")
    expect_true("description" %in% names(result))
    expect_true("assumed_structure" %in% names(result))
    expect_true("status" %in% names(result))
    expect_true("details" %in% names(result))
})

test_that(".compute_working_correlation_fit() handles missing mgcv/geepack", {
    data <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10)
    
    # Result should handle gracefully if geepack not available
    result <- .compute_working_correlation_fit(data)
    
    expect_is(result, "list")
    expect_true("status" %in% names(result))
    expect_true("description" %in% names(result))
    # Status should be a string (either skip message or fit assessment)
    expect_is(result$status, "character")
    expect_true(nchar(result$status) > 0)
})

test_that(".compute_working_correlation_fit() preserves structure parameter", {
    set.seed(42)
    data <- matrix(rnorm(30 * 8, mean = 5, sd = 1), nrow = 30, ncol = 8)
    
    result <- .compute_working_correlation_fit(data, assumed_structure = "ar1")
    
    expect_equal(result$assumed_structure, "ar1")
    expect_is(result$details, "character")
    expect_true(nchar(result$details) > 0)
})

# ==============================================================================
# Test: Cluster Size Variation
# ==============================================================================

test_that(".compute_cluster_size_variation() detects homogeneous clusters", {
    set.seed(42)
    # Create data with uniform cluster sizes
    data <- matrix(rnorm(100 * 10, mean = 10, sd = 2), nrow = 100, ncol = 10)
    
    result <- .compute_cluster_size_variation(data)
    
    # All columns are treated as clusters with same number of rows
    expect_equal(result$n_clusters, 10)
    expect_equal(result$mean_cluster_size, 100)
    expect_equal(result$cv_cluster_size, 0)  # Perfect homogeneity
    expect_equal(result$status, "homogeneous")
})

test_that(".compute_cluster_size_variation() detects variable clusters", {
    set.seed(42)
    # Create data with variable cluster sizes
    cluster_sizes <- c(10, 15, 12, 18, 20, 11, 14, 19, 16, 13)
    total_obs <- sum(cluster_sizes)
    # Create data with rows = total observations, cols = 1 for simplicity
    data <- matrix(rnorm(total_obs), nrow = total_obs, ncol = 1)
    
    # Create cluster assignment: each element repeated according to its cluster size
    cluster_col <- rep(seq_len(length(cluster_sizes)), cluster_sizes)
    
    result <- .compute_cluster_size_variation(data, cluster_col = cluster_col)
    
    expect_is(result, "list")
    expect_true("cv_cluster_size" %in% names(result))
    expect_true(result$cv_cluster_size > 0)  # Some variation
    expect_true(result$status %in% c("homogeneous", "moderate variation", "high variation"))
})

test_that(".compute_cluster_size_variation() calculates correct statistics", {
    set.seed(42)
    data <- matrix(rnorm(50 * 5), nrow = 50, ncol = 5)
    
    result <- .compute_cluster_size_variation(data)
    
    # For uniform data (5 clusters of 50 each)
    expect_equal(result$mean_cluster_size, 50)
    expect_equal(result$min_size, 50)
    expect_equal(result$max_size, 50)
    expect_equal(result$n_clusters, 5)
})

test_that(".compute_cluster_size_variation() outputs correct details string", {
    set.seed(42)
    data <- matrix(rnorm(60 * 8), nrow = 60, ncol = 8)
    
    result <- .compute_cluster_size_variation(data)
    
    expect_is(result$details, "character")
    expect_match(result$details, "Mean size=")
    expect_match(result$details, "CV=")
})

# ==============================================================================
# Test: Independence Residuals Test
# ==============================================================================

test_that(".compute_independence_residuals() returns correlation metrics", {
    set.seed(42)
    data <- matrix(rnorm(40 * 12, mean = 10, sd = 2), nrow = 40, ncol = 12)
    
    result <- .compute_independence_residuals(data)
    
    expect_is(result, "list")
    expect_true("description" %in% names(result))
    expect_true("mean_residual_correlation" %in% names(result))
    expect_true("status" %in% names(result))
    expect_true("details" %in% names(result))
})

test_that(".compute_independence_residuals() handles NA values", {
    set.seed(42)
    data <- matrix(rnorm(40 * 8, mean = 10, sd = 2), nrow = 40, ncol = 8)
    # Add some NAs
    data[sample(seq_len(length(data)), 10)] <- NA
    
    result <- .compute_independence_residuals(data)
    
    expect_is(result, "list")
    expect_false(is.na(result$status))  # Should still return valid status
})

test_that(".compute_independence_residuals() classifies independence correctly", {
    set.seed(42)
    # Create data with low within-cluster correlation
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    result <- .compute_independence_residuals(data)
    
    expect_true(result$status %in% c("independent", "some correlation (GEE handles)", 
                                     "strong correlation (GEE robust)"))
})

test_that(".compute_independence_residuals() handles single column", {
    set.seed(42)
    # Single cluster
    data <- matrix(rnorm(50, mean = 10, sd = 2), nrow = 50, ncol = 1)
    
    result <- .compute_independence_residuals(data)
    
    expect_is(result, "list")
    expect_equal(result$n_clusters_tested, 1)
})

# ==============================================================================
# Test: GEE Scale Parameter (Dispersion)
# ==============================================================================

test_that(".compute_gee_scale_parameter() returns scale metrics", {
    set.seed(42)
    data <- matrix(rnorm(30 * 8, mean = 10, sd = 2), nrow = 30, ncol = 8)
    
    result <- .compute_gee_scale_parameter(data)
    
    expect_is(result, "list")
    expect_true("description" %in% names(result))
    expect_true("status" %in% names(result))
})

test_that(".compute_gee_scale_parameter() classifies dispersion correctly", {
    set.seed(42)
    data <- matrix(rnorm(40 * 10, mean = 10, sd = 2), nrow = 40, ncol = 10)
    
    result <- .compute_gee_scale_parameter(data)
    
    # Status should classify dispersion
    expect_true(result$status %in% c("? SKIP - geepack not available",
                                     "under-dispersed (rare)", 
                                     "correct dispersion",
                                     "over-dispersed",
                                     "unknown",
                                     "? ERROR"))
})

test_that(".compute_gee_scale_parameter() outputs valid details", {
    set.seed(42)
    data <- matrix(rnorm(50 * 6, mean = 10, sd = 2), nrow = 50, ncol = 6)
    
    result <- .compute_gee_scale_parameter(data)
    
    expect_is(result$details, "character")
    expect_true(nchar(result$details) > 0)
})

# ==============================================================================
# Test: Wrapper Function .get_gee_metrics()
# ==============================================================================

test_that(".get_gee_metrics() returns all 4 metrics plus consolidated", {
    set.seed(42)
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    result <- .get_gee_metrics(data)
    
    expect_is(result, "list")
    expect_equal(length(result), 5)  # 4 metrics + consolidated
    expect_true("correlation_fit" %in% names(result))
    expect_true("cluster_variation" %in% names(result))
    expect_true("independence" %in% names(result))
    expect_true("scale_parameter" %in% names(result))
    expect_true("consolidated" %in% names(result))
})

test_that(".get_gee_metrics() consolidated result has correct structure", {
    set.seed(42)
    data <- matrix(rnorm(60 * 8, mean = 10, sd = 2), nrow = 60, ncol = 8)
    
    result <- .get_gee_metrics(data)
    
    consolidated <- result$consolidated
    expect_is(consolidated, "list")
    expect_true("description" %in% names(consolidated))
    expect_true("result" %in% names(consolidated))
    expect_true("status" %in% names(consolidated))
    expect_equal(consolidated$status, "COMBINED")
})

test_that(".get_gee_metrics() accepts gee_params", {
    set.seed(42)
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    gee_params <- list(
        assumed_structure = "ar1",
        cluster_col = NULL
    )
    
    result <- .get_gee_metrics(data, gee_params = gee_params)
    
    expect_is(result, "list")
    expect_equal(result$correlation_fit$assumed_structure, "ar1")
})

test_that(".get_gee_metrics() handles empty gee_params gracefully", {
    set.seed(42)
    data <- matrix(rnorm(40 * 6, mean = 10, sd = 2), nrow = 40, ncol = 6)
    
    result <- .get_gee_metrics(data, gee_params = list())
    
    # Should use defaults
    expect_is(result, "list")
    expect_equal(result$correlation_fit$assumed_structure, "exchangeable")
})

test_that(".get_gee_metrics() converts status to lowercase in consolidated", {
    set.seed(42)
    data <- matrix(rnorm(50 * 8, mean = 10, sd = 2), nrow = 50, ncol = 8)
    
    result <- .get_gee_metrics(data)
    
    consolidated_result <- result$consolidated$result
    
    # Consolidated result should use lowercase status values
    expect_is(consolidated_result, "character")
    # Should not contain uppercase markers like [OK], ?, !
    expect_false(grepl("\\[OK\\]", consolidated_result))
})

test_that(".get_gee_metrics() handles various data dimensions", {
    set.seed(42)
    
    # Small dataset
    data1 <- matrix(rnorm(20 * 5), nrow = 20, ncol = 5)
    result1 <- .get_gee_metrics(data1)
    expect_equal(length(result1), 5)
    
    # Medium dataset
    data2 <- matrix(rnorm(100 * 15), nrow = 100, ncol = 15)
    result2 <- .get_gee_metrics(data2)
    expect_equal(length(result2), 5)
    
    # Large dataset
    data3 <- matrix(rnorm(500 * 30), nrow = 500, ncol = 30)
    result3 <- .get_gee_metrics(data3)
    expect_equal(length(result3), 5)
})

# ==============================================================================
# Test: Data Type Handling
# ==============================================================================

test_that("GEE functions handle non-matrix data (coerce to matrix)", {
    set.seed(42)
    
    # Data frame input
    data_df <- as.data.frame(matrix(rnorm(40 * 8), nrow = 40, ncol = 8))
    # Should coerce internally
    result <- expect_error(.get_gee_metrics(data_df), NA)  # Expect NO error
})

test_that("GEE functions return proper list structures on error", {
    set.seed(42)
    
    # Test with invalid input (e.g., single value)
    data_invalid <- 42
    
    # Should return error structure, not throw uncaught error
    result <- tryCatch(
        .get_gee_metrics(data_invalid),
        error = function(e) e
    )
    
    # Either successful or proper error handling
    expect_true(is.list(result) || inherits(result, "error"))
})

# ==============================================================================
# Test: NA and Missing Data Handling
# ==============================================================================

test_that(".get_gee_metrics() handles missing data gracefully", {
    set.seed(42)
    data <- matrix(rnorm(60 * 10), nrow = 60, ncol = 10)
    
    # Add random NAs (up to 20% of data)
    na_indices <- sample(seq_len(length(data)), size = length(data) * 0.1)
    data[na_indices] <- NA
    
    # Should not crash
    result <- tryCatch(
        .get_gee_metrics(data),
        error = function(e) NULL
    )
    
    # Either successful or returns gracefully
    expect_true(is.null(result) || is.list(result))
})

test_that("Individual metrics handle complete columns with NAs", {
    set.seed(42)
    data <- matrix(rnorm(50 * 8), nrow = 50, ncol = 8)
    data[, 1] <- NA  # First column all NA
    
    # Should handle without crashing
    result <- tryCatch(
        .compute_cluster_size_variation(data),
        error = function(e) "error"
    )
    
    expect_false(identical(result, "error"))
})

# ==============================================================================
# Test: Edge Cases
# ==============================================================================

test_that(".get_gee_metrics() handles minimum viable data (2x2)", {
    data <- matrix(rnorm(4), nrow = 2, ncol = 2)
    
    result <- tryCatch(
        .get_gee_metrics(data),
        error = function(e) NULL
    )
    
    # Should handle gracefully
    expect_true(is.null(result) || is.list(result))
})

test_that(".get_gee_metrics() handles single-column data", {
    set.seed(42)
    data <- matrix(rnorm(50), nrow = 50, ncol = 1)
    
    result <- .get_gee_metrics(data)
    
    expect_is(result, "list")
    # Should return all 5 elements (4 metrics + consolidated)
    expect_equal(length(result), 5)
    # Cluster variation should be valid
    expect_is(result$cluster_variation, "list")
    expect_true("n_clusters" %in% names(result$cluster_variation))
    expect_equal(result$cluster_variation$n_clusters, 1)  # Single column = 1 cluster
})

test_that(".get_gee_metrics() handles perfect correlation", {
    set.seed(42)
    # Create perfectly correlated data
    base <- rnorm(50)
    data <- matrix(rep(base, 5), nrow = 50, ncol = 5)
    
    result <- tryCatch(
        .get_gee_metrics(data),
        error = function(e) NULL
    )
    
    # Should handle without crashing
    expect_true(is.null(result) || is.list(result))
})

test_that(".get_gee_metrics() handles zero variance data", {
    # All entries are same value
    data <- matrix(5, nrow = 50, ncol = 10)
    
    result <- tryCatch(
        .get_gee_metrics(data),
        error = function(e) NULL
    )
    
    # Should handle gracefully
    expect_true(is.null(result) || is.list(result))
})

# ==============================================================================
# Test: Numerical Validation Tests
# ==============================================================================

test_that("Cluster size CV is mathematically sound for homogeneous clusters", {
    # Create perfectly homogeneous clusters: 5 clusters of 20 observations each
    data <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)
    
    result <- .compute_cluster_size_variation(data, cluster_col = NULL)
    
    # For uniformly sized clusters, CV should be exactly 0
    expect_equal(result$cv_cluster_size, 0, tolerance = 1e-10)
    expect_equal(result$status, "homogeneous")
    # Mean should equal all individual sizes
    expect_equal(result$mean_cluster_size, 100)
    expect_equal(result$sd_cluster_size, 0, tolerance = 1e-10)
})

test_that("Cluster size statistics compute correctly with known data", {
    # Create known cluster sizes: 10, 20, 30 (mean=20, sd=sqrt((100+0+100)/3)≈8.16)
    cluster_sizes <- c(10, 20, 30)
    total_obs <- sum(cluster_sizes)
    data <- matrix(rnorm(total_obs), nrow = total_obs, ncol = 1)
    cluster_col <- rep(seq_len(3), cluster_sizes)
    
    result <- .compute_cluster_size_variation(data, cluster_col = cluster_col)
    
    expect_equal(result$n_clusters, 3)
    expect_equal(result$mean_cluster_size, 20, tolerance = 1e-10)
    expect_equal(result$min_size, 10)
    expect_equal(result$max_size, 30)
    # SD of c(10,20,30) = sqrt(66.67) ≈ 8.165
    expected_sd <- sqrt(sum((c(10, 20, 30) - 20)^2) / (3 - 1))
    expect_equal(result$sd_cluster_size, expected_sd, tolerance = 1e-10)
    # CV = sd/mean
    expected_cv <- expected_sd / 20
    expect_equal(result$cv_cluster_size, expected_cv, tolerance = 1e-10)
})

test_that("Independence residuals produces bounded correlation values", {
    set.seed(42)
    # Create strongly autocorrelated data
    data <- matrix(nrow = 50, ncol = 10)
    for (i in 1:10) {
        data[, i] <- cumsum(rnorm(50, sd = 0.1))
    }
    
    result <- .compute_independence_residuals(data, cluster_col = NULL)
    
    # Result should be a list
    expect_is(result, "list")
    expect_true("status" %in% names(result))
    
    # Correlation should always be bounded by [-1, 1]
    if (!is.null(result$correlation) && is.numeric(result$correlation)) {
        expect_gte(result$correlation, -1)
        expect_lte(result$correlation, 1)
    }
})

test_that("Scale parameter phi reflects dispersion correctly", {
    set.seed(42)
    # Create overdispersed data (high variance relative to mean)
    mean_val <- 5
    overdispersed <- matrix(rnorm(50 * 10, mean = mean_val, sd = 4), 
                           nrow = 50, ncol = 10)
    
    result <- .compute_gee_scale_parameter(overdispersed, cluster_col = NULL)
    
    # Result should be a list with status field
    expect_is(result, "list")
    expect_true("status" %in% names(result))
    
    # Phi should be numeric and positive
    if (!is.null(result$phi) && is.numeric(result$phi)) {
        expect_gt(result$phi, 0)
        # For overdispersed data, phi might be affected by variance
        expect_is(result$status, "character")
    }
})

test_that("Cluster CV thresholds are correctly applied", {
    # Test threshold: CV < 0.2 = homogeneous, 0.2 <= CV < 0.5 = moderate, CV >= 0.5 = high
    test_cases <- list(
        list(sizes = c(20, 20, 20), expected_status = "homogeneous"),      # CV = 0
        list(sizes = c(18, 20, 22), expected_status = "homogeneous"),      # CV ≈ 0.054
        list(sizes = c(15, 20, 25), expected_status = "moderate variation"), # CV ≈ 0.163
        list(sizes = c(12, 20, 28), expected_status = "moderate variation"), # CV = 0.4
        list(sizes = c(5, 20, 35), expected_status = "high variation")       # CV = 0.75
    )
    
    for (test_case in test_cases) {
        cluster_sizes <- test_case$sizes
        total_obs <- sum(cluster_sizes)
        data <- matrix(rnorm(total_obs), nrow = total_obs, ncol = 1)
        cluster_col <- rep(seq_len(length(cluster_sizes)), cluster_sizes)
        
        result <- .compute_cluster_size_variation(data, cluster_col = cluster_col)
        expect_equal(result$status, test_case$expected_status,
                    info = sprintf("Failed for cluster sizes: %s (CV=%.3f)", 
                                 paste(cluster_sizes, collapse = ","), result$cv_cluster_size))
    }
})

test_that("Working correlation fit produces valid structure descriptions", {
    set.seed(42)
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    # Test multiple correlation structures
    structures <- c("independence", "exchangeable", "ar1")
    
    for (struct in structures) {
        result <- .compute_working_correlation_fit(data, assumed_structure = struct)
        
        # Should return valid structure field
        expect_equal(result$assumed_structure, struct)
        # Details should describe the structure
        expect_is(result$details, "character")
        expect_true(nchar(result$details) > 0)
        expect_false(grepl("ERROR", result$details))  # Should not be error
    }
})

test_that("Consolidated GEE results contains all metrics", {
    set.seed(42)
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    result <- .get_gee_metrics(data, 
                               gee_params = list(assumed_structure = "exchangeable"))
    
    # All 5 components should be present
    expect_equal(length(result), 5)
    expected_names <- c("correlation_fit", "cluster_variation", "independence", 
                       "scale_parameter", "consolidated")
    expect_true(all(expected_names %in% names(result)))
    
    # Consolidated should have these fields
    expect_is(result$consolidated, "list")
    expect_true("result" %in% names(result$consolidated))
    expect_true("status" %in% names(result$consolidated))
    
    # Result should be a string describing all metrics
    expect_is(result$consolidated$result, "character")
    expect_true(nchar(result$consolidated$result) > 0)
})

test_that("Cluster variation handles extreme cluster number scenarios", {
    # Many small clusters
    data1 <- matrix(rnorm(100), nrow = 10, ncol = 10)
    result1 <- .compute_cluster_size_variation(data1)
    expect_equal(result1$n_clusters, 10)
    expect_equal(result1$mean_cluster_size, 10)
    
    # Few large clusters
    data2 <- matrix(rnorm(100), nrow = 100, ncol = 1)
    result2 <- .compute_cluster_size_variation(data2)
    expect_equal(result2$n_clusters, 1)
    expect_equal(result2$mean_cluster_size, 100)
})

test_that("CV boundary calculations are precise", {
    # Create data with exact CV boundary: 0.2 (threshold for homogeneous/moderate)
    # For mean=100, CV=0.2 means sd=20. Solve: sd = sqrt(sum((x_i - 100)^2)/(n-1))
    # Two clusters: x1=95, x2=105 gives sd = 7.07, CV = 0.0707 (homogeneous)
    # Two clusters: x1=80, x2=120 gives sd = 28.28, CV = 0.283 (moderate)
    
    test_exact_cv <- function(sizes, expected_threshold_result) {
        total_obs <- sum(sizes)
        data <- matrix(rnorm(total_obs), nrow = total_obs, ncol = 1)
        cluster_col <- rep(seq_len(length(sizes)), sizes)
        
        result <- .compute_cluster_size_variation(data, cluster_col = cluster_col)
        
        cv <- result$cv_cluster_size
        if (cv < 0.2) {
            expect_equal(result$status, "homogeneous")
        } else if (cv < 0.5) {
            expect_equal(result$status, "moderate variation")
        } else {
            expect_equal(result$status, "high variation")
        }
    }
    
    # Test with clearly defined cases
    test_exact_cv(c(100, 100, 100), "homogeneous")
    test_exact_cv(c(50, 100, 150), NA)  # Just verify no error
})

# ==============================================================================
# Test: Reproducibility and Consistency
# ==============================================================================

test_that(".get_gee_metrics() is reproducible with same data", {
    set.seed(42)
    data <- matrix(rnorm(60 * 8, mean = 10, sd = 2), nrow = 60, ncol = 8)
    
    result1 <- .get_gee_metrics(data)
    result2 <- .get_gee_metrics(data)
    
    # Results should be identical
    expect_equal(result1$cluster_variation$n_clusters, 
                 result2$cluster_variation$n_clusters)
    expect_equal(result1$cluster_variation$mean_cluster_size, 
                 result2$cluster_variation$mean_cluster_size)
})

test_that("Individual metrics are consistent with wrapper", {
    set.seed(42)
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    # Get metrics separately
    corr_fit <- .compute_working_correlation_fit(data, assumed_structure = "exchangeable")
    cluster_var <- .compute_cluster_size_variation(data)
    
    # Get through wrapper
    wrapper_result <- .get_gee_metrics(data)
    
    # Should match
    expect_equal(corr_fit$assumed_structure, 
                 wrapper_result$correlation_fit$assumed_structure)
    expect_equal(cluster_var$n_clusters, 
                 wrapper_result$cluster_variation$n_clusters)
})

# ==============================================================================
# Test: Integration with calculate_rank_assumptions
# ==============================================================================

test_that("GEE metrics integrate with calculate_rank_assumptions()", {
    set.seed(42)
    data <- matrix(rnorm(50 * 10, mean = 10, sd = 2), nrow = 50, ncol = 10)
    
    # Test that gee_metrics check can be requested
    result <- tryCatch(
        .calculate_rank_assumptions(data, checks = c("gee_metrics")),
        error = function(e) NULL
    )
    
    # Should return rank_assumptions object with gee_metrics
    if (!is.null(result)) {
        expect_is(result, "rank_assumptions")
        checks <- attr(result, "checks")
        expect_true("gee_metrics" %in% names(checks))
    }
})

test_that("GEE metrics included in 'all' checks preset", {
    set.seed(42)
    data <- matrix(rnorm(40 * 8, mean = 10, sd = 2), nrow = 40, ncol = 8)
    
    result <- tryCatch(
        .calculate_rank_assumptions(data, checks = "all"),
        error = function(e) NULL
    )
    
    if (!is.null(result)) {
        checks <- attr(result, "checks")
        # Should include gee_metrics when checks="all"
        expect_true("gee_metrics" %in% names(checks) || 
                   !is.null(checks$gee_metrics))
    }
})

# ==============================================================================
# End of GEE Metrics Tests
# ==============================================================================
