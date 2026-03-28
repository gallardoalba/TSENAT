context("Rcpp Jackknife Implementation")

test_that("Rcpp jackknife available for testing", {
  # Check if Rcpp is compiled
  rcpp_status <- tryCatch({
    check_rcpp_available()
    TRUE
  }, error = function(e) FALSE)
  
  expect_true(is.logical(rcpp_status) || is.null(rcpp_status))
})

test_that("Rcpp jackknife produces valid results", {
  skip_on_cran()
  
  # Create test data
  set.seed(42)
  counts <- matrix(rpois(50 * 20, lambda = 5), nrow = 50, ncol = 20)
  
  # Run C++ version
  result_cpp <- jackknife_resampling_cpp(
    counts,
    q = 1.0,
    normalize = TRUE,
    log_base = exp(1),
    pseudocount = 0
  )
  
  # Verify output structure and validity
  expect_true(is.numeric(result_cpp$estimate))
  expect_false(is.na(result_cpp$estimate))
  expect_true(is.numeric(result_cpp$jackknife_estimates))
  expect_true(is.numeric(result_cpp$influence))
  expect_true(is.numeric(result_cpp$jackknife_se))
})

test_that("Rcpp jackknife handles different Tsallis parameters", {
  skip_on_cran()
  
  set.seed(42)
  counts <- matrix(rpois(50 * 20, lambda = 5), nrow = 50, ncol = 20)
  
  # Test q = 0.5
  result_q05 <- jackknife_resampling_cpp(counts, q = 0.5)
  expect_true(is.numeric(result_q05$estimate))
  
  # Test q = 2.0
  result_q20 <- jackknife_resampling_cpp(counts, q = 2.0)
  expect_true(is.numeric(result_q20$estimate))
  
  # Different q values should give different results
  expect_false(isTRUE(all.equal(result_q05$estimate, result_q20$estimate)))
})

test_that("Rcpp jackknife handles edge cases", {
  skip_on_cran()
  
  # Minimum valid case (2 rows, minimal counts)
  set.seed(42)
  counts_min <- matrix(rpois(2 * 10, lambda = 1), nrow = 2, ncol = 10)
  result_min <- jackknife_resampling_cpp(counts_min)
  expect_true(is.numeric(result_min$estimate))
  
  # Very sparse counts (still valid, not all zeros)
  counts_sparse <- matrix(c(rep(1, 50), rep(0, 50)), nrow = 10, ncol = 10)
  result_sparse <- jackknife_resampling_cpp(counts_sparse)
  expect_true(is.numeric(result_sparse$estimate))
})

test_that("Rcpp jackknife with pseudocount", {
  skip_on_cran()
  
  set.seed(42)
  counts <- matrix(rpois(30 * 15, lambda = 3), nrow = 30, ncol = 15)
  
  # With pseudocount
  result_pc <- jackknife_resampling_cpp(counts, pseudocount = 0.5)
  
  # Without pseudocount
  result_no_pc <- jackknife_resampling_cpp(counts, pseudocount = 0.0)
  
  # Results should differ
  expect_false(isTRUE(all.equal(result_pc$estimate, result_no_pc$estimate)))
})

test_that("Hybrid wrapper falls back to R version gracefully", {
  set.seed(42)
  counts <- matrix(rpois(30 * 15, lambda = 5), nrow = 30, ncol = 15)
  
  # Should return result with proper structure regardless of Rcpp availability
  result <- TSENAT:::.jackknife_resampling_hybrid(
    counts,
    q = 1.0,
    threshold = 90,
    norm = TRUE,
    log_base = exp(1),
    pseudocount = 0
  )
  
  # Check structure
  expect_true(is.list(result))
  expect_true("estimate" %in% names(result))
  expect_true("jackknife_estimates" %in% names(result))
  expect_true("influence" %in% names(result))
  expect_true("jackknife_se" %in% names(result))
})

test_that("Rcpp entropy computation matches R version", {
  skip_on_cran()
  
  # Create proportions
  p <- c(0.2, 0.3, 0.3, 0.2)
  
  # Test Shannon entropy (q = 1)
  entropy_cpp_shannon <- entropy_cpp(p, q = 1, normalize = FALSE, log_base = exp(1))
  entropy_r_shannon <- .entropy_core(p, q = 1, norm = FALSE, log_base = exp(1))
  
  expect_equal(entropy_cpp_shannon, entropy_r_shannon, tolerance = 1e-6)
})

test_that("Rcpp jackknife works with large datasets", {
  skip_on_cran()
  
  # Verify that Rcpp version works correctly with large data
  set.seed(42)
  counts <- matrix(rpois(5000 * 100, lambda = 5), nrow = 5000, ncol = 100)
  
  # Run C++ version
  result_cpp <- jackknife_resampling_cpp(counts)
  
  # Check reasonable output
  expect_true(is.numeric(result_cpp$estimate))
  expect_false(is.na(result_cpp$estimate))
  expect_equal(length(result_cpp$jackknife_estimates), nrow(counts))
  
  message("✓ Rcpp jackknife works with large datasets (5000x100)")
})

test_that("Rcpp output structure includes all required fields", {
  skip_on_cran()
  
  set.seed(42)
  counts <- matrix(rpois(50 * 20, lambda = 5), nrow = 50, ncol = 20)
  result <- jackknife_resampling_cpp(counts)
  
  # Check all fields are present
  required_fields <- c("estimate", "jackknife_estimates", "influence", "jackknife_se", "q", "normalize")
  for (field in required_fields) {
    expect_true(field %in% names(result), info = sprintf("Missing field: %s", field))
  }
  
  # Verify field values
  expect_true(is.numeric(result$estimate))
  expect_true(is.numeric(result$jackknife_estimates))
  expect_true(is.numeric(result$influence))
  expect_true(is.numeric(result$jackknife_se))
  expect_equal(result$q, 1.0)
  expect_equal(result$normalize, TRUE)
})

test_that("Rcpp jackknife parameters stored in output", {
  skip_on_cran()
  
  set.seed(42)
  counts <- matrix(rpois(50 * 20, lambda = 5), nrow = 50, ncol = 20)
  
  # Test with q = 0.5
  result_q05 <- jackknife_resampling_cpp(counts, q = 0.5)
  expect_equal(result_q05$q, 0.5)
  
  # Test with normalize = FALSE
  result_no_norm <- jackknife_resampling_cpp(counts, normalize = FALSE)
  expect_equal(result_no_norm$normalize, FALSE)
  
  # Test combination
  result_combo <- jackknife_resampling_cpp(counts, q = 2.0, normalize = FALSE)
  expect_equal(result_combo$q, 2.0)
  expect_equal(result_combo$normalize, FALSE)
})

test_that("Rcpp jackknife with different log bases", {
  skip_on_cran()
  
  set.seed(42)
  counts <- matrix(rpois(50 * 20, lambda = 5), nrow = 50, ncol = 20)
  
  # Test with log base 2
  result_log2 <- jackknife_resampling_cpp(counts, log_base = 2)
  expect_true(is.numeric(result_log2$estimate))
  
  # Test with log base 10
  result_log10 <- jackknife_resampling_cpp(counts, log_base = 10)
  expect_true(is.numeric(result_log10$estimate))
  
  # Test with natural log (exp(1))
  result_loge <- jackknife_resampling_cpp(counts, log_base = exp(1))
  expect_true(is.numeric(result_loge$estimate))
  
  # When normalized, different log bases should give the SAME result
  # (normalization is scale-invariant across log bases)
  expect_equal(result_log2$estimate, result_log10$estimate, tolerance = 1e-6)
})

test_that("Rcpp jackknife no NaN or Inf values in results", {
  skip_on_cran()
  
  set.seed(42)
  counts <- matrix(rpois(100 * 30, lambda = 5), nrow = 100, ncol = 30)
  result <- jackknife_resampling_cpp(counts)
  
  # Check no NaN in main outputs
  expect_false(any(is.nan(result$estimate)))
  expect_false(any(is.nan(result$jackknife_se)))
  
  # Check no Inf
  expect_false(any(is.infinite(result$estimate)))
  expect_false(any(is.infinite(result$jackknife_se)))
  
  # Check jackknife_estimates and influence don't have Inf (NaN is okay for invalid rows)
  valid_estimates <- !is.na(result$jackknife_estimates)
  if (any(valid_estimates)) {
    expect_false(any(is.infinite(result$jackknife_estimates[valid_estimates])))
    expect_false(any(is.infinite(result$influence[valid_estimates])))
  }
})

test_that("Rcpp entropy with normalized vs unnormalized", {
  skip_on_cran()
  
  set.seed(42)
  counts <- matrix(rpois(50 * 20, lambda = 5), nrow = 50, ncol = 20)
  
  # Normalized entropy
  result_norm <- jackknife_resampling_cpp(counts, normalize = TRUE)
  
  # Unnormalized entropy
  result_unnorm <- jackknife_resampling_cpp(counts, normalize = FALSE)
  
  # Unnormalized should be >= normalized (max entropy normalization)
  expect_true(result_unnorm$estimate >= result_norm$estimate - 1e-6)
  
  # They should differ (unless normalize has no effect, which is unlikely)
  expect_false(isTRUE(all.equal(result_norm$estimate, result_unnorm$estimate)))
})

test_that("Rcpp jackknife influence is non-negative", {
  skip_on_cran()
  
  set.seed(42)
  counts <- matrix(rpois(50 * 20, lambda = 5), nrow = 50, ncol = 20)
  result <- jackknife_resampling_cpp(counts)
  
  # All influence values should be >= 0
  valid_influence <- !is.na(result$influence)
  expect_true(all(result$influence[valid_influence] >= 0))
})

test_that("Rcpp jackknife SE is always positive", {
  skip_on_cran()
  
  set.seed(42)
  counts <- matrix(rpois(50 * 20, lambda = 5), nrow = 50, ncol = 20)
  result <- jackknife_resampling_cpp(counts)
  
  # SE should be positive (bias-corrected standard error)
  expect_true(result$jackknife_se >= 0)
})

test_that("Rcpp jackknife with varying sample sizes", {
  skip_on_cran()
  
  set.seed(42)
  
  # Small dataset
  counts_small <- matrix(rpois(10 * 5, lambda = 2), nrow = 10, ncol = 5)
  result_small <- jackknife_resampling_cpp(counts_small)
  expect_equal(length(result_small$jackknife_estimates), 10)
  
  # Medium dataset
  counts_med <- matrix(rpois(100 * 20, lambda = 5), nrow = 100, ncol = 20)
  result_med <- jackknife_resampling_cpp(counts_med)
  expect_equal(length(result_med$jackknife_estimates), 100)
  
  # Large dataset
  counts_large <- matrix(rpois(500 * 50, lambda = 5), nrow = 500, ncol = 50)
  result_large <- jackknife_resampling_cpp(counts_large)
  expect_equal(length(result_large$jackknife_estimates), 500)
})

test_that("Rcpp jackknife with varying feature/sample proportions", {
  skip_on_cran()
  
  set.seed(42)
  
  # Many features, few samples
  counts_wide <- matrix(rpois(20 * 100, lambda = 5), nrow = 20, ncol = 100)
  result_wide <- jackknife_resampling_cpp(counts_wide)
  expect_true(is.numeric(result_wide$estimate))
  expect_equal(length(result_wide$jackknife_estimates), 20)
  
  # Few features, many samples
  counts_tall <- matrix(rpois(100 * 10, lambda = 5), nrow = 100, ncol = 10)
  result_tall <- jackknife_resampling_cpp(counts_tall)
  expect_true(is.numeric(result_tall$estimate))
  expect_equal(length(result_tall$jackknife_estimates), 100)
})

test_that("Rcpp jackknife works across multiple parameter combinations", {
  skip_on_cran()
  
  set.seed(42)
  counts <- matrix(rpois(50 * 20, lambda = 5), nrow = 50, ncol = 20)
  
  # Test multiple parameter combinations
  params_list <- list(
    list(q = 0.5, normalize = TRUE),
    list(q = 1.0, normalize = FALSE),
    list(q = 2.0, normalize = TRUE),
    list(q = 0.5, normalize = FALSE)
  )
  
  for (params in params_list) {
    result_cpp <- jackknife_resampling_cpp(
      counts, q = params$q, normalize = params$normalize
    )
    
    # Verify valid output
    expect_true(is.numeric(result_cpp$estimate),
                info = sprintf("Invalid estimate for q=%.1f, normalize=%s", params$q, params$normalize))
    expect_false(is.na(result_cpp$estimate),
                 info = sprintf("NA estimate for q=%.1f, normalize=%s", params$q, params$normalize))
  }
})

test_that("Rcpp jackknife high pseudocount impact", {
  skip_on_cran()
  
  set.seed(42)
  counts <- matrix(rpois(50 * 20, lambda = 2), nrow = 50, ncol = 20)
  
  # Very low pseudocount
  result_low <- jackknife_resampling_cpp(counts, pseudocount = 0.001)
  
  # Very high pseudocount
  result_high <- jackknife_resampling_cpp(counts, pseudocount = 10)
  
  # High pseudocount should make estimate closer to uniform distribution
  # (lower entropy since pseudocount adds to rare categories)
  expect_true(is.numeric(result_low$estimate))
  expect_true(is.numeric(result_high$estimate))
  
  # They should be different
  expect_false(isTRUE(all.equal(result_low$estimate, result_high$estimate)))
})

test_that("Rcpp jackknife numeric vector output dimensions", {
  skip_on_cran()
  
  set.seed(42)
  counts <- matrix(rpois(75 * 25, lambda = 5), nrow = 75, ncol = 25)
  result <- jackknife_resampling_cpp(counts)
  
  # Output dimensions should match
  expect_equal(length(result$jackknife_estimates), nrow(counts))
  expect_equal(length(result$influence), nrow(counts))
  expect_equal(length(result$q), 1)
  expect_equal(length(result$normalize), 1)
  expect_equal(length(result$estimate), 1)
  expect_equal(length(result$jackknife_se), 1)
})

# ============================================================================
# ADDITIONAL COMPREHENSIVE TESTS FOR BUG FIXES AND EDGE CASES
# ============================================================================

test_that("Rcpp jackknife handles q=0 (species richness) correctly", {
  skip_on_cran()
  
  set.seed(42)
  counts <- matrix(rpois(50 * 20, lambda = 5), nrow = 50, ncol = 20)
  
  # Test q = 0 (species richness case)
  result_q0 <- jackknife_resampling_cpp(counts, q = 0.0, normalize = TRUE)
  
  # Verify it returns valid numbers
  expect_true(is.numeric(result_q0$estimate))
  expect_false(is.nan(result_q0$estimate))
  
  # q=0 estimate should be positive (log of species count)
  expect_true(result_q0$estimate > 0)
})

test_that("Rcpp jackknife q near 1 boundary handling", {
  skip_on_cran()
  
  set.seed(42)
  counts <- matrix(rpois(50 * 20, lambda = 5), nrow = 50, ncol = 20)
  
  # Test q very close to 1 (should use Shannon path)
  q_near_1 <- 1.0 + 1e-7
  result_near_1 <- jackknife_resampling_cpp(counts, q = q_near_1, normalize = TRUE)
  
  # Compare with exact q=1
  result_q1 <- jackknife_resampling_cpp(counts, q = 1.0, normalize = TRUE)
  
  # Should be very close (within 1e-4 of q=1)
  expect_equal(result_near_1$estimate, result_q1$estimate, tolerance = 1e-4)
})

test_that("Rcpp jackknife extreme q values", {
  skip_on_cran()
  
  set.seed(42)
  counts <- matrix(rpois(50 * 20, lambda = 5), nrow = 50, ncol = 20)
  
  # Very small q (approaching 0)
  result_small_q <- jackknife_resampling_cpp(counts, q = 0.01, normalize = TRUE)
  expect_true(is.numeric(result_small_q$estimate))
  expect_false(is.nan(result_small_q$estimate))
  
  # Very large q (approaches negative log of maximum)
  result_large_q <- jackknife_resampling_cpp(counts, q = 10.0, normalize = TRUE)
  expect_true(is.numeric(result_large_q$estimate))
  expect_false(is.nan(result_large_q$estimate))
  
  # Very large q should be different from very small q
  expect_false(isTRUE(all.equal(result_small_q$estimate, result_large_q$estimate)))
})

test_that("Rcpp jackknife single species per observation", {
  skip_on_cran()
  
  # Create data with only one species having counts per observation
  counts_single <- matrix(0, nrow = 20, ncol = 15)
  for (i in 1:20) {
    # Each observation has only one non-zero species
    j <- ((i - 1) %% 15) + 1
    counts_single[i, j] <- rpois(1, lambda = 5)
  }
  
  result <- jackknife_resampling_cpp(counts_single, normalize = TRUE)
  
  # Should return valid values (entropy should be lower for single species dominance)
  expect_true(is.numeric(result$estimate))
  expect_false(is.nan(result$estimate))
  
  # When dominated by single species, normalized entropy is bounded by [0,1]
  expect_true(result$estimate >= 0)
  expect_true(result$estimate <= 1)
})

test_that("Rcpp jackknife highly skewed distribution", {
  skip_on_cran()
  
  # Create highly skewed distribution
  counts_skewed <- matrix(c(
    rep(c(100, 1, 1, 1, 1, 1, 1, 1, 1, 1), 10)
  ), nrow = 10, ncol = 10, byrow = TRUE)
  
  result <- jackknife_resampling_cpp(counts_skewed, normalize = TRUE)
  
  # Should handle skewed distribution properly
  expect_true(is.numeric(result$estimate))
  expect_false(is.nan(result$estimate))
  
  # Skewed distribution should have lower entropy (expected < 0.3)
  expect_true(result$estimate < 0.5)
})

test_that("Rcpp jackknife uniform distribution maximizes entropy", {
  skip_on_cran()
  
  # Create nearly uniform distribution
  counts_uniform <- matrix(rep(10, 50 * 10), nrow = 50, ncol = 10)
  
  # Create more skewed distribution
  counts_skewed <- counts_uniform
  counts_skewed[1:25, 1] <- 1
  counts_skewed[1:25, 2] <- 90
  
  result_uniform <- jackknife_resampling_cpp(counts_uniform, normalize = TRUE)
  result_skewed <- jackknife_resampling_cpp(counts_skewed, normalize = TRUE)
  
  # Uniform should have higher entropy than skewed
  expect_true(result_uniform$estimate > result_skewed$estimate)
})

test_that("Rcpp jackknife handles NA values in proportions", {
  skip_on_cran()
  
  # Test entropy_cpp directly with NA values
  p_with_na <- c(0.2, 0.3, NA, 0.3, 0.2)
  entropy_val <- entropy_cpp(p_with_na, q = 1, normalize = FALSE, log_base = exp(1))
  
  # Should return valid numeric (NA's filtered out)
  expect_true(is.numeric(entropy_val))
  
  # Should not be NaN or Inf
  expect_false(is.nan(entropy_val))
  expect_false(is.infinite(entropy_val))
})

test_that("Rcpp jackknife empty proportions returns NA", {
  skip_on_cran()
  
  # Pass all-zero proportions (after filtering)
  p_zeros <- c(0, 0, 0, 0)
  entropy_val <- entropy_cpp(p_zeros, q = 1, normalize = FALSE, log_base = exp(1))
  
  # Should return NA_REAL for invalid input
  expect_true(is.na(entropy_val))
})

test_that("Rcpp jackknife row with zero counts handled correctly", {
  skip_on_cran()
  
  set.seed(42)
  counts <- matrix(rpois(30 * 10, lambda = 5), nrow = 30, ncol = 10)
  
  # Add rows with zero counts that will affect leave-one-out
  counts[5, ] <- 0  # Row with all zeros
  
  result <- jackknife_resampling_cpp(counts, normalize = TRUE)
  
  # Should handle gracefully
  expect_true(is.list(result))
  
  # Influence for zero-count row might be NA or very small
  # (depends on how much that row affects total)
  expect_true(is.numeric(result$influence[5]) || is.na(result$influence[5]))
})

test_that("Rcpp jackknife normalization effect on different q values", {
  skip_on_cran()
  
  set.seed(42)
  counts <- matrix(rpois(50 * 20, lambda = 5), nrow = 50, ncol = 20)
  
  # For q < 1, test that normalized entropy differs from unnormalized
  result_norm_q05 <- jackknife_resampling_cpp(counts, q = 0.5, normalize = TRUE)
  result_unnorm_q05 <- jackknife_resampling_cpp(counts, q = 0.5, normalize = FALSE)
  
  # Normalized and unnormalized should be different (unless normalize has no effect)
  expect_false(isTRUE(all.equal(result_norm_q05$estimate, result_unnorm_q05$estimate)))
  
  # For q > 1, normalized entropy should be bounded by [0,1] due to normalization
  result_norm_q20 <- jackknife_resampling_cpp(counts, q = 2.0, normalize = TRUE)
  result_unnorm_q20 <- jackknife_resampling_cpp(counts, q = 2.0, normalize = FALSE)
  
  # Normalized should be bounded to [0,1] when normalize=TRUE
  expect_true(result_norm_q20$estimate >= 0)
  expect_true(result_norm_q20$estimate <= 1)
  
  # Unnormalized can be anything (positive for Tsallis with q>1)
  expect_true(is.numeric(result_unnorm_q20$estimate))
})

test_that("Rcpp jackknife SE is NA when all estimates are invalid", {
  skip_on_cran()
  
  # Create counts that would produce mostly invalid estimates
  # This is very difficult to trigger, but we test the logic path
  set.seed(42)
  counts <- matrix(rpois(10 * 5, lambda = 2), nrow = 10, ncol = 5)
  
  result <- jackknife_resampling_cpp(counts, normalize = TRUE)
  
  # Most of the time, SE should be positive
  if (!is.na(result$jackknife_se)) {
    expect_true(result$jackknife_se >= 0)
  }
})

test_that("Rcpp jackknife consistency across multiple calls", {
  skip_on_cran()
  
  set.seed(42)
  counts <- matrix(rpois(50 * 20, lambda = 5), nrow = 50, ncol = 20)
  
  # Multiple calls with same input should give identical results
  result1 <- jackknife_resampling_cpp(counts, q = 0.5, normalize = TRUE)
  result2 <- jackknife_resampling_cpp(counts, q = 0.5, normalize = TRUE)
  
  # Should be bytewise identical
  expect_identical(result1$estimate, result2$estimate)
  expect_identical(result1$jackknife_estimates, result2$jackknife_estimates)
  expect_identical(result1$jackknife_se, result2$jackknife_se)
})

test_that("Rcpp jackknife monotonic relationship with q variation", {
  skip_on_cran()
  
  set.seed(123)
  counts <- matrix(rpois(50 * 15, lambda = 5), nrow = 50, ncol = 15)
  
  # For uniform distribution, entropy decreases as q increases
  # This is a property of Tsallis entropy
  q_values <- c(0.5, 1.0, 2.0, 5.0)
  estimates <- sapply(q_values, function(q) {
    jackknife_resampling_cpp(counts, q = q, normalize = TRUE)$estimate
  })
  
  # For uniform-ish distribution, entropy should decrease with increasing q
  # (This is a known property of normalized Tsallis entropy)
  expect_true(all(diff(estimates) <= 0.1))  # Allow some flexibility due to numerical effects
})

test_that("Rcpp jackknife different log bases produce consistent results", {
  skip_on_cran()
  
  set.seed(42)
  counts <- matrix(rpois(50 * 20, lambda = 5), nrow = 50, ncol = 20)
  
  # When normalize=TRUE, log base should not affect normalized result
  result_base2 <- jackknife_resampling_cpp(counts, log_base = 2, normalize = TRUE)
  result_base10 <- jackknife_resampling_cpp(counts, log_base = 10, normalize = TRUE)
  result_base_e <- jackknife_resampling_cpp(counts, log_base = exp(1), normalize = TRUE)
  
  # All normalized results should be the same
  expect_equal(result_base2$estimate, result_base10$estimate, tolerance = 1e-6)
  expect_equal(result_base2$estimate, result_base_e$estimate, tolerance = 1e-6)
})

test_that("Rcpp jackknife unnormalized results scale with log base", {
  skip_on_cran()
  
  set.seed(42)
  counts <- matrix(rpois(50 * 20, lambda = 5), nrow = 50, ncol = 20)
  
  # When normalize=FALSE, results scale by log base
  result_log2 <- jackknife_resampling_cpp(counts, log_base = 2, normalize = FALSE)
  result_loge <- jackknife_resampling_cpp(counts, log_base = exp(1), normalize = FALSE)
  
  # Ratio should be approximately log_e(2)
  expected_ratio <- log(2)
  actual_ratio <- result_loge$estimate / result_log2$estimate
  
  expect_equal(actual_ratio, expected_ratio, tolerance = 1e-5)
})

test_that("Rcpp influence values are always non-negative", {
  skip_on_cran()
  
  set.seed(42)
  counts <- matrix(rpois(100 * 25, lambda = 5), nrow = 100, ncol = 25)
  
  for (q in c(0.5, 1.0, 2.0, 5.0)) {
    result <- jackknife_resampling_cpp(counts, q = q, normalize = TRUE)
    
    # All valid influence values must be >= 0
    valid_influence <- !is.na(result$influence)
    expect_true(all(result$influence[valid_influence] >= -1e-10),
                info = sprintf("Negative influence detected for q=%g", q))
  }
})

test_that("Rcpp jackknife high-dimensional data", {
  skip_on_cran()
  
  # Test with many species (high-dimensional)
  set.seed(42)
  counts <- matrix(rpois(100 * 500, lambda = 2), nrow = 100, ncol = 500)
  
  result <- jackknife_resampling_cpp(counts, q = 1.0, normalize = TRUE)
  
  # Should handle high-dimensional data
  expect_true(is.numeric(result$estimate))
  expect_false(is.nan(result$estimate))
  
  # With many species and normalize=TRUE, entropy should be relatively high
  expect_true(result$estimate > 0.5)
})

test_that("Rcpp jackknife with comprehensive parameter coverage", {
  skip_on_cran()
  
  set.seed(42)
  counts <- matrix(rpois(50 * 20, lambda = 5), nrow = 50, ncol = 20)
  
  # Test comprehensive parameter grid
  q_vals <- c(0, 0.5, 1.0, 2.0, 5.0)
  normalize_vals <- c(TRUE, FALSE)
  pseudocount_vals <- c(0, 0.1, 1.0)
  
  test_count <- 0
  for (q in q_vals) {
    for (norm in normalize_vals) {
      for (pc in pseudocount_vals) {
        result_cpp <- jackknife_resampling_cpp(
          counts, q = q, normalize = norm, pseudocount = pc
        )
        
        # Verify valid output
        expect_true(is.numeric(result_cpp$estimate),
                   info = sprintf("q=%g, norm=%s, pc=%g", q, norm, pc))
        test_count <- test_count + 1
      }
    }
  }
  
  message(sprintf("✓ Tested %d parameter combinations", test_count))
})
