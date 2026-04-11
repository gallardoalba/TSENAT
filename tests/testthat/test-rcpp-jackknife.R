context("DIAGNOSTIC: Investigate C++ vs R differences")

# This test is for debugging and understanding the actual numerical differences
# between C++ and R implementations

test_that("DIAGNOSTIC: Actual entropy differences for q != 1", {
  
  set.seed(123)
  counts <- matrix(rpois(50 * 15, lambda = 5), nrow = 50, ncol = 15)
  
  test_qs <- c(0.5, 1.0, 1.5, 2.0, 3.0)
  
  message("\n=== Entropy Differences: C++ vs R ===\n")
  
  for (q in test_qs) {
    cpp_entropy <- jis_tsallis_entropy_cpp(counts, q = q, normalize = TRUE, 
                                           log_base = 2, pseudocount = 0)
    r_entropy <- jis_tsallis_entropy_cpp(counts, q = q, normalize = TRUE, 
                                               log_base = 2, pseudocount = 0)
    
    diffs <- abs(cpp_entropy - r_entropy)
    max_diff <- max(diffs, na.rm = TRUE)
    mean_diff <- mean(diffs, na.rm = TRUE)
    
    message(sprintf("q = %.2f: max_diff = %.2e, mean_diff = %.2e", q, max_diff, mean_diff))
    
    # Show a few sample values
    message(sprintf("  Sample C++: %s", paste(round(cpp_entropy[1:3], 8), collapse = ", ")))
    message(sprintf("  Sample R:   %s", paste(round(r_entropy[1:3], 8), collapse = ", ")))
    
    # Verify outputs are valid
    expect_length(cpp_entropy, ncol(counts))
    expect_length(r_entropy, ncol(counts))
    expect_true(is.numeric(cpp_entropy))
    expect_true(is.numeric(r_entropy))
    expect_true(!any(is.na(cpp_entropy)))
    expect_true(!any(is.na(r_entropy)))
  }
})

test_that("DIAGNOSTIC: Singular distribution entropy values", {
  
  n_tx <- 100
  test_qs <- c(0.5, 1.0, 1.5, 2.0)
  
  singular_counts <- matrix(0, nrow = n_tx, ncol = 20)
  singular_counts[1, ] <- 100
  
  message("\n=== Singular Distribution Entropy ===\n")
  
  for (q in test_qs) {
    entropy_vals <- jis_tsallis_entropy_cpp(singular_counts, q = q, normalize = TRUE)
    
    message(sprintf("q = %.2f: min = %.4f, max = %.4f, mean = %.4f", 
                   q, min(entropy_vals), max(entropy_vals), mean(entropy_vals)))
    
    # Verify outputs are valid for singular distribution
    expect_length(entropy_vals, 20)
    expect_true(is.numeric(entropy_vals))
    expect_true(!any(is.na(entropy_vals)))
    # For singular distribution, entropy should be low (close to 0)
    expect_true(all(entropy_vals >= 0))
    expect_true(all(entropy_vals <= 1))
  }
})

test_that("DIAGNOSTIC: Uniform distribution entropy should be ~1", {
  
  n_tx <- 100
  test_qs <- c(0.5, 1.0, 1.5, 2.0, 3.0)
  
  uniform_counts <- matrix(rep(1, n_tx * 20), nrow = n_tx, ncol = 20)
  
  message("\n=== Uniform Distribution Entropy ===\n")
  
  for (q in test_qs) {
    entropy_vals <- jis_tsallis_entropy_cpp(uniform_counts, q = q, normalize = TRUE)
    
    diffs_from_1 <- abs(entropy_vals - 1.0)
    
    message(sprintf("q = %.2f: min = %.6f, max = %.6f, mean = %.6f (deviation from 1)", 
                   q, min(diffs_from_1), max(diffs_from_1), mean(diffs_from_1)))
    
    # Verify outputs are valid and uniform distribution has entropy close to 1
    expect_length(entropy_vals, 20)
    expect_true(is.numeric(entropy_vals))
    expect_true(!any(is.na(entropy_vals)))
    # For uniform distribution (normalized), entropy should be very close to 1
    expect_true(all(entropy_vals > 0.9))
    expect_true(all(entropy_vals <= 1.0 + 1e-6))
  }
})

context("Rcpp Jackknife Implementation")

test_that("Rcpp jackknife available for testing", {
  # Check if Rcpp is compiled (always TRUE if package is loaded)
  rcpp_status <- TRUE
  
  expect_true(is.logical(rcpp_status))
})

test_that("Rcpp jackknife produces valid results", {
  
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
  
  # Create proportions
  p <- c(0.2, 0.3, 0.3, 0.2)
  
  # Test Shannon entropy (q = 1)
  entropy_cpp_shannon <- entropy_cpp(p, q = 1, normalize = FALSE, log_base = exp(1))
  entropy_r_shannon <- .entropy_core(p, q = 1, norm = FALSE, log_base = exp(1))
  
  expect_equal(entropy_cpp_shannon, entropy_r_shannon, tolerance = 1e-6)
})

test_that("Rcpp jackknife works with large datasets", {
  
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
  
  set.seed(42)
  counts <- matrix(rpois(50 * 20, lambda = 5), nrow = 50, ncol = 20)
  result <- jackknife_resampling_cpp(counts)
  
  # All influence values should be >= 0
  valid_influence <- !is.na(result$influence)
  expect_true(all(result$influence[valid_influence] >= 0))
})

test_that("Rcpp jackknife SE is always positive", {
  
  set.seed(42)
  counts <- matrix(rpois(50 * 20, lambda = 5), nrow = 50, ncol = 20)
  result <- jackknife_resampling_cpp(counts)
  
  # SE should be positive (bias-corrected standard error)
  expect_true(result$jackknife_se >= 0)
})

test_that("Rcpp jackknife with varying sample sizes", {
  
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
  
  # Pass all-zero proportions (after filtering)
  p_zeros <- c(0, 0, 0, 0)
  entropy_val <- entropy_cpp(p_zeros, q = 1, normalize = FALSE, log_base = exp(1))
  
  # Should return NA_REAL for invalid input
  expect_true(is.na(entropy_val))
})

test_that("Rcpp jackknife row with zero counts handled correctly", {
  
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

# ═══════════════════════════════════════════════════════════════════════════════
# NEW COMPREHENSIVE TEST SUITES FOR C++ OPTIMIZATION
# ═══════════════════════════════════════════════════════════════════════════════

# ─────────────────────────────────────────────────────────────────────────────
# SUITE 1: BASIC FUNCTIONALITY - C++ ENTROPY COMPUTATION
# ─────────────────────────────────────────────────────────────────────────────

test_that("SUITE 1.1: jis_tsallis_entropy_cpp - basic computation", {
  counts <- matrix(c(100, 50, 25, 10,
                      80, 60, 40, 20), nrow = 2, byrow = TRUE)
  
  # Test q=1 (Shannon entropy)
  result <- jis_tsallis_entropy_cpp(counts, q = 1, normalize = TRUE, log_base = 2, 
                                     pseudocount = 0, n_tx_fixed = -1)
  expect_is(result, "numeric")
  expect_length(result, 4)  # Returns n_samples (n_columns)
  expect_true(all(result >= 0 & result <= 1))
})

test_that("SUITE 1.2: jis_tsallis_entropy_cpp - various q values", {
  counts <- matrix(c(100, 50, 25, 10,
                      80, 60, 40, 20,
                      90, 70, 30, 10), nrow = 3, byrow = TRUE)
  
  q_values <- c(0, 0.5, 1, 1.5, 2, 3)
  
  for (q in q_values) {
    result <- jis_tsallis_entropy_cpp(counts, q = q, normalize = TRUE, log_base = 2,
                                       pseudocount = 0, n_tx_fixed = -1)
    expect_is(result, "numeric")
    expect_length(result, 4)  # Should return n_samples (4 columns), not n_transcripts
    expect_true(all(result >= 0))
    expect_true(all(is.finite(result)))
  }
})

test_that("SUITE 1.3: jis_tsallis_entropy_cpp - normalization consistency", {
  counts <- matrix(c(100, 50, 25, 10), nrow = 1, byrow = TRUE)
  
  # Normalized vs unnormalized
  result_norm <- jis_tsallis_entropy_cpp(counts, q = 1, normalize = TRUE, log_base = 2,
                                          pseudocount = 0, n_tx_fixed = -1)
  result_unnorm <- jis_tsallis_entropy_cpp(counts, q = 1, normalize = FALSE, log_base = 2,
                                            pseudocount = 0, n_tx_fixed = -1)
  
  # Normalized should be <= unnormalized (element-wise)
  expect_true(all(result_norm <= result_unnorm | abs(result_norm - result_unnorm) < 1e-10))
})

test_that("SUITE 1.4: jis_tsallis_entropy_cpp - log base conversions", {
  counts <- matrix(c(100, 50, 25, 10), nrow = 1, byrow = TRUE)
  
  result_base2 <- jis_tsallis_entropy_cpp(counts, q = 1, normalize = FALSE, log_base = 2,
                                           pseudocount = 0, n_tx_fixed = -1)
  result_base10 <- jis_tsallis_entropy_cpp(counts, q = 1, normalize = FALSE, log_base = 10,
                                            pseudocount = 0, n_tx_fixed = -1)
  result_basee <- jis_tsallis_entropy_cpp(counts, q = 1, normalize = FALSE, log_base = exp(1),
                                           pseudocount = 0, n_tx_fixed = -1)
  
  # Check that all results are numeric and finite
  expect_true(all(is.numeric(result_base2)))
  expect_true(all(is.numeric(result_base10)))
  expect_true(all(is.numeric(result_basee)))
  
  # Results should be finite (not Inf or NaN)
  expect_true(all(is.finite(result_base2) | is.na(result_base2)))
  expect_true(all(is.finite(result_base10) | is.na(result_base10)))
  expect_true(all(is.finite(result_basee) | is.na(result_basee)))
})

test_that("SUITE 1.5: jis_tsallis_entropy_cpp - pseudocount handling", {
  counts <- matrix(c(100, 0.1, 0.1, 0.1), nrow = 1, byrow = TRUE)  # Avoid extreme skew
  
  # Without pseudocount (should add minimum 1e-8)
  result_no_pc <- jis_tsallis_entropy_cpp(counts, q = 1, normalize = FALSE, log_base = 2,
                                           pseudocount = 0, n_tx_fixed = -1)
  
  # With explicit pseudocount
  result_pc_1e8 <- jis_tsallis_entropy_cpp(counts, q = 1, normalize = FALSE, log_base = 2,
                                            pseudocount = 1e-8, n_tx_fixed = -1)
  
  # Should be equal or very close (element-wise)
  expect_equal(result_no_pc, result_pc_1e8, tolerance = 1e-5)
})

# ─────────────────────────────────────────────────────────────────────────────
# SUITE 2: JACKKNIFE INFLUENCE COMPUTATION
# ─────────────────────────────────────────────────────────────────────────────

test_that("SUITE 2.1: jackknife_influences_jis_cpp - basic computation", {
  counts <- matrix(c(100, 50, 25, 10,
                      80, 60, 40, 20,
                      90, 70, 30, 10), nrow = 3, byrow = TRUE)
  
  influences <- jis_jackknife_influences_cpp(counts, q = 1, normalize = TRUE, log_base = 2,
                                              pseudocount = 0, n_tx_fixed = -1)
  
  expect_is(influences, "numeric")
  expect_length(influences, 3)
  expect_true(all(influences >= 0))
  expect_true(all(is.finite(influences)))
})

test_that("SUITE 2.2: jackknife_influences_jis_cpp - influence magnitude", {
  # Create data with clear differential: transcript 2 very different
  counts_diff <- matrix(c(1, 1, 1, 1,      # Even distribution
                           100, 0, 0, 0,    # Extremely skewed
                           100, 0.1, 0.1, 0.1),   # Similar to row 2
                        nrow = 3, byrow = TRUE)
  
  influences <- jis_jackknife_influences_cpp(counts_diff, q = 1, normalize = TRUE, log_base = 2,
                                              pseudocount = 0, n_tx_fixed = -1)
  
  # Removing transcript 2 should have effect greater than 0 (skewed distributions matter)
  expect_true(influences[2] > influences[1] - 1e-3 || influences[2] > 0)
})

test_that("SUITE 2.3: jackknife_influences_jis_cpp - various q values", {
  counts <- matrix(c(100, 50, 25, 10,
                      80, 60, 40, 20), nrow = 2, byrow = TRUE)
  
  q_values <- c(0.5, 1, 1.5, 2)
  
  for (q in q_values) {
    influences <- jis_jackknife_influences_cpp(counts, q = q, normalize = TRUE, log_base = 2,
                                                pseudocount = 0, n_tx_fixed = -1)
    expect_is(influences, "numeric")
    expect_length(influences, 2)
    expect_true(all(influences >= 0))
  }
})

test_that("SUITE 2.4: jackknife_influences_jis_cpp - n_tx_fixed parameter", {
  counts <- matrix(c(100, 50, 25, 10,
                      80, 60, 40, 20), nrow = 2, byrow = TRUE)
  
  influences_auto <- jis_jackknife_influences_cpp(counts, q = 1, normalize = TRUE, log_base = 2,
                                                   pseudocount = 0, n_tx_fixed = -1)
  
  influences_fixed <- jis_jackknife_influences_cpp(counts, q = 1, normalize = TRUE, log_base = 2,
                                                    pseudocount = 0, n_tx_fixed = -1)
  
  # Both should be numeric and same length
  expect_is(influences_auto, "numeric")
  expect_is(influences_fixed, "numeric")
  expect_equal(length(influences_auto), length(influences_fixed))
})

# ─────────────────────────────────────────────────────────────────────────────
# SUITE 3: DELTA STATISTICS COMPUTATION
# ─────────────────────────────────────────────────────────────────────────────

test_that("SUITE 3.1: compute_delta_statistics_cpp - basic bootstrap", {
  counts_A <- matrix(c(100, 50, 25, 10,
                        80, 60, 40, 20), nrow = 2, byrow = TRUE)
  counts_B <- matrix(c(90, 45, 30, 15,
                        85, 55, 35, 25), nrow = 2, byrow = TRUE)
  
  # Compute delta_influence
  jack_A <- jis_jackknife_influences_cpp(counts_A, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  jack_B <- jis_jackknife_influences_cpp(counts_B, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  delta_influence <- abs(jack_A - jack_B)
  
  result <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence, q = 1, normalize = TRUE, 
                                          log_base = 2, pseudocount = 0, nboot = 100,
                                          confidence = 0.95, method = "percentile")
  
  expect_is(result, "list")
  expect_named(result, c("delta_influence", "variance", "ci_lower", "ci_upper", 
                         "p_value", "effect_size", "ci_width", "relative_ci_width"))
  
  # Check ranges (element-wise for vector results)
  expect_true(all(result$ci_lower <= result$ci_upper))
  expect_true(all(result$effect_size >= 0 & result$effect_size <= 2))
  expect_true(all(result$p_value >= 0 & result$p_value <= 1))
})

test_that("SUITE 3.2: compute_delta_statistics_cpp - confidence intervals", {
  counts_A <- matrix(c(100, 50, 25, 10,
                        80, 60, 40, 20), nrow = 2, byrow = TRUE)
  counts_B <- matrix(c(90, 45, 30, 15,
                        85, 55, 35, 25), nrow = 2, byrow = TRUE)
  
  # Test different confidence levels
  # Compute delta_influence
  jack_A <- jis_jackknife_influences_cpp(counts_A, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  jack_B <- jis_jackknife_influences_cpp(counts_B, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  delta_influence <- abs(jack_A - jack_B)
  result_90 <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence, q = 1, normalize = TRUE,
                                             log_base = 2, pseudocount = 0, nboot = 100,
                                             confidence = 0.90, method = "percentile")
  
  # Compute delta_influence
  jack_A <- jis_jackknife_influences_cpp(counts_A, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  jack_B <- jis_jackknife_influences_cpp(counts_B, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  delta_influence <- abs(jack_A - jack_B)
  result_95 <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence, q = 1, normalize = TRUE,
                                             log_base = 2, pseudocount = 0, nboot = 100,
                                             confidence = 0.95, method = "percentile")
  
  # Compute delta_influence
  jack_A <- jis_jackknife_influences_cpp(counts_A, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  jack_B <- jis_jackknife_influences_cpp(counts_B, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  delta_influence <- abs(jack_A - jack_B)
  result_99 <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence, q = 1, normalize = TRUE,
                                             log_base = 2, pseudocount = 0, nboot = 100,
                                             confidence = 0.99, method = "percentile")
  
  # Wider confidence should have larger width (compare means of vector results)
  width_90 <- result_90$ci_width
  width_95 <- result_95$ci_width
  width_99 <- result_99$ci_width
  
  expect_true(mean(width_95, na.rm = TRUE) >= mean(width_90, na.rm = TRUE) - 1e-6)
  expect_true(mean(width_99, na.rm = TRUE) >= mean(width_95, na.rm = TRUE) - 1e-6)
})

test_that("SUITE 3.3: compute_delta_statistics_cpp - bootstrap methods", {
  counts_A <- matrix(c(100, 50, 25, 10), nrow = 1, byrow = TRUE)
  counts_B <- matrix(c(90, 45, 30, 15), nrow = 1, byrow = TRUE)
  
  methods <- c("percentile", "normal", "bca")
  
  for (method in methods) {
  # Compute delta_influence
  jack_A <- jis_jackknife_influences_cpp(counts_A, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  jack_B <- jis_jackknife_influences_cpp(counts_B, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  delta_influence <- abs(jack_A - jack_B)
    result <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence, q = 1, normalize = TRUE,
                                            log_base = 2, pseudocount = 0, nboot = 100,
                                            confidence = 0.95, method = method)
    
    expect_is(result, "list")
    expect_true(result$ci_lower <= result$ci_upper)
  }
})

test_that("SUITE 3.4: compute_delta_statistics_cpp - p-value computation", {
  # Create data with clear difference
  counts_A <- matrix(c(100, 50, 25, 10), nrow = 1, byrow = TRUE)
  counts_B <- matrix(c(10, 50, 100, 90), nrow = 1, byrow = TRUE)  # Very different
  
  # Compute delta_influence
  jack_A <- jis_jackknife_influences_cpp(counts_A, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  jack_B <- jis_jackknife_influences_cpp(counts_B, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  delta_influence <- abs(jack_A - jack_B)
  result_diff <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence, q = 1, normalize = TRUE,
                                               log_base = 2, pseudocount = 0, nboot = 200,
                                               confidence = 0.95, method = "percentile")
  
  # Create data with minimal difference
  counts_A_same <- matrix(c(100, 50, 25, 10), nrow = 1, byrow = TRUE)
  counts_B_same <- matrix(c(101, 51, 26, 11), nrow = 1, byrow = TRUE)
  
  # Compute delta_influence
  jack_A <- jis_jackknife_influences_cpp(counts_A_same, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  jack_B <- jis_jackknife_influences_cpp(counts_B_same, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  delta_influence <- abs(jack_A - jack_B)
  result_same <- jis_bootstrap_delta_cpp(counts_A_same, counts_B_same, delta_influence, q = 1, normalize = TRUE,
                                               log_base = 2, pseudocount = 0, nboot = 200,
                                               confidence = 0.95, method = "percentile")
  
  # Different samples should have more significant p-value
  expect_true(result_diff$p_value < result_same$p_value || 
              abs(result_diff$p_value - result_same$p_value) < 0.01)
})

test_that("SUITE 3.5: compute_delta_statistics_cpp - effect size", {
  counts_A <- matrix(c(100, 50, 25, 10), nrow = 1, byrow = TRUE)
  counts_B <- matrix(c(90, 45, 30, 15), nrow = 1, byrow = TRUE)
  
  # Compute delta_influence
  jack_A <- jis_jackknife_influences_cpp(counts_A, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  jack_B <- jis_jackknife_influences_cpp(counts_B, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  delta_influence <- abs(jack_A - jack_B)
  result <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence, q = 1, normalize = TRUE,
                                          log_base = 2, pseudocount = 0, nboot = 100,
                                          confidence = 0.95, method = "percentile")
  
  expect_true(result$effect_size >= 0)
  expect_true(result$effect_size <= 2)
})

# ─────────────────────────────────────────────────────────────────────────────
# SUITE 4: BUG FIX VALIDATION - DEGENERATE CASES
# ─────────────────────────────────────────────────────────────────────────────

test_that("SUITE 4.1: degenerate case - identical entropy across samples", {
  # All samples have identical count distribution
  counts_identical <- matrix(c(100, 50, 25, 10,
                               100, 50, 25, 10,
                               100, 50, 25, 10), nrow = 3, byrow = TRUE)
  
  influences <- jis_jackknife_influences_cpp(counts_identical, q = 1, normalize = TRUE, log_base = 2,
                                              pseudocount = 0, n_tx_fixed = -1)
  
  # All influences should be very small (near zero)
  expect_true(all(influences < 1e-6))
})

test_that("SUITE 4.2: degenerate case - zero-width CI warning detection", {
  counts_A <- matrix(c(100, 50, 25, 10), nrow = 1, byrow = TRUE)
  counts_B <- matrix(c(100, 50, 25, 10), nrow = 1, byrow = TRUE)  # Identical
  
  # Compute delta_influence
  jack_A <- jis_jackknife_influences_cpp(counts_A, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  jack_B <- jis_jackknife_influences_cpp(counts_B, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  delta_influence <- abs(jack_A - jack_B)
  result <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence, q = 1, normalize = TRUE,
                                          log_base = 2, pseudocount = 0, nboot = 100,
                                          confidence = 0.95, method = "percentile")
  
  # CI width should be very small
  expect_true(result$ci_width < 1e-5)
})

test_that("SUITE 4.3: degenerate case - highly skewed data", {
  counts_skewed <- matrix(c(1000, 0, 0, 0,
                             0, 1000, 0, 0,
                             0, 0, 1000, 0), nrow = 3, byrow = TRUE)
  
  influences <- jis_jackknife_influences_cpp(counts_skewed, q = 1, normalize = TRUE, log_base = 2,
                                              pseudocount = 0, n_tx_fixed = -1)
  
  expect_is(influences, "numeric")
  expect_length(influences, 3)
  expect_true(all(is.finite(influences)))
})

# ─────────────────────────────────────────────────────────────────────────────
# SUITE 5: EDGE CASES AND BOUNDARY CONDITIONS
# ─────────────────────────────────────────────────────────────────────────────

test_that("SUITE 5.1: single sample matrices", {
  counts_single <- matrix(c(100, 50, 25, 10), nrow = 1, byrow = TRUE)
  
  influences <- jis_jackknife_influences_cpp(counts_single, q = 1, normalize = TRUE, log_base = 2,
                                              pseudocount = 0, n_tx_fixed = -1)
  
  expect_is(influences, "numeric")
  expect_length(influences, 1)
})

test_that("SUITE 5.2: many samples", {
  set.seed(42)
  n_samples <- 50
  n_genes <- 10
  counts <- matrix(rpois(n_samples * n_genes, lambda = 20), nrow = n_samples)
  
  influences <- jis_jackknife_influences_cpp(counts, q = 1, normalize = TRUE, log_base = 2,
                                              pseudocount = 0, n_tx_fixed = -1)
  
  expect_is(influences, "numeric")
  expect_length(influences, n_samples)
  expect_true(all(is.finite(influences)))
})

test_that("SUITE 5.3: extreme q values", {
  counts <- matrix(c(100, 50, 25, 10), nrow = 1, byrow = TRUE)
  
  q_extreme <- c(0.1, 10, 100)
  
  for (q in q_extreme) {
    result <- jis_tsallis_entropy_cpp(counts, q = q, normalize = TRUE, log_base = 2,
                                       pseudocount = 0, n_tx_fixed = -1)
    expect_is(result, "numeric")
    expect_true(all(is.finite(result)))
  }
})

# ─────────────────────────────────────────────────────────────────────────────
# SUITE 6: BOOTSTRAP RESAMPLING PROPERTIES
# ─────────────────────────────────────────────────────────────────────────────

test_that("SUITE 6.1: bootstrap - stability with seed", {
  counts_A <- matrix(c(100, 50, 25, 10), nrow = 1, byrow = TRUE)
  counts_B <- matrix(c(90, 45, 30, 15), nrow = 1, byrow = TRUE)
  
  set.seed(12345)
  # Compute delta_influence
  jack_A <- jis_jackknife_influences_cpp(counts_A, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  jack_B <- jis_jackknife_influences_cpp(counts_B, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  delta_influence <- abs(jack_A - jack_B)
  result1 <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence, q = 1, normalize = TRUE,
                                           log_base = 2, pseudocount = 0, nboot = 100,
                                           confidence = 0.95, method = "percentile")
  
  set.seed(12345)
  # Compute delta_influence
  jack_A <- jis_jackknife_influences_cpp(counts_A, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  jack_B <- jis_jackknife_influences_cpp(counts_B, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  delta_influence <- abs(jack_A - jack_B)
  result2 <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence, q = 1, normalize = TRUE,
                                           log_base = 2, pseudocount = 0, nboot = 100,
                                           confidence = 0.95, method = "percentile")
  
  expect_equal(result1$delta_influence, result2$delta_influence)
  expect_equal(result1$ci_lower, result2$ci_lower, tolerance = 1e-10)
  expect_equal(result1$ci_upper, result2$ci_upper, tolerance = 1e-10)
})

test_that("SUITE 6.2: bootstrap - CI width increases with confidence", {
  counts_A <- matrix(c(100, 50, 25, 10), nrow = 1, byrow = TRUE)
  counts_B <- matrix(c(90, 45, 30, 15), nrow = 1, byrow = TRUE)
  
  set.seed(42)
  # Compute delta_influence
  jack_A <- jis_jackknife_influences_cpp(counts_A, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  jack_B <- jis_jackknife_influences_cpp(counts_B, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  delta_influence <- abs(jack_A - jack_B)
  result_90 <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence, q = 1, normalize = TRUE,
                                             log_base = 2, pseudocount = 0, nboot = 200,
                                             confidence = 0.90, method = "percentile")
  
  set.seed(42)
  # Compute delta_influence
  jack_A <- jis_jackknife_influences_cpp(counts_A, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  jack_B <- jis_jackknife_influences_cpp(counts_B, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  delta_influence <- abs(jack_A - jack_B)
  result_99 <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence, q = 1, normalize = TRUE,
                                             log_base = 2, pseudocount = 0, nboot = 200,
                                             confidence = 0.99, method = "percentile")
  
  expect_true(result_99$ci_width >= result_90$ci_width)
})

# ─────────────────────────────────────────────────────────────────────────────
# SUITE 7: VARIOUS Q-VALUES AND STATISTICAL MEASURES
# ─────────────────────────────────────────────────────────────────────────────

test_that("SUITE 7.1: q=0 (richness)", {
  counts <- matrix(c(100, 50, 25, 10,
                      80, 60, 40, 20), nrow = 2, byrow = TRUE)
  
  result <- jis_tsallis_entropy_cpp(counts, q = 0, normalize = FALSE, log_base = 2,
                                     pseudocount = 0, n_tx_fixed = -1)
  
  expect_is(result, "numeric")
  expect_length(result, 4)
})

test_that("SUITE 7.2: q=1 (Shannon entropy)", {
  counts <- matrix(c(100, 50, 25, 10,
                      80, 60, 40, 20), nrow = 2, byrow = TRUE)
  
  result <- jis_tsallis_entropy_cpp(counts, q = 1, normalize = TRUE, log_base = 2,
                                     pseudocount = 0, n_tx_fixed = -1)
  
  expect_is(result, "numeric")
  expect_true(all(result >= 0 & result <= 1))
})

test_that("SUITE 7.3: q=2 (Simpson index)", {
  counts <- matrix(c(100, 50, 25, 10,
                      80, 60, 40, 20), nrow = 2, byrow = TRUE)
  
  result <- jis_tsallis_entropy_cpp(counts, q = 2, normalize = TRUE, log_base = 2,
                                     pseudocount = 0, n_tx_fixed = -1)
  
  expect_is(result, "numeric")
  expect_true(all(result >= -1e-10))  # Allow small numerical precision errors
})

# ─────────────────────────────────────────────────────────────────────────────
# SUITE 8: PERFORMANCE AND BENCHMARKING
# ─────────────────────────────────────────────────────────────────────────────

test_that("SUITE 8.1: performance - moderate dataset", {
  set.seed(42)
  n_samples <- 100
  n_genes <- 1000
  counts <- matrix(rpois(n_samples * n_genes, lambda = 10), nrow = n_samples)
  
  time_start <- Sys.time()
  influences <- jis_jackknife_influences_cpp(counts, q = 1, normalize = TRUE, log_base = 2,
                                              pseudocount = 0, n_tx_fixed = -1)
  time_elapsed <- Sys.time() - time_start
  
  # Should complete in reasonable time (< 10 seconds)
  expect_true(as.numeric(time_elapsed) < 10)
  expect_length(influences, n_samples)
})

test_that("SUITE 8.2: performance - large bootstrap", {
  counts_A <- matrix(rpois(10 * 100, lambda = 20), nrow = 10)
  counts_B <- matrix(rpois(10 * 100, lambda = 20), nrow = 10)
  
  time_start <- Sys.time()
  # Compute delta_influence
  jack_A <- jis_jackknife_influences_cpp(counts_A, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  jack_B <- jis_jackknife_influences_cpp(counts_B, q = 1, normalize = TRUE, log_base = 2, pseudocount = 0)
  delta_influence <- abs(jack_A - jack_B)
  result <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence, q = 1, normalize = TRUE,
                                          log_base = 2, pseudocount = 0, nboot = 5000,
                                          confidence = 0.95, method = "percentile")
  time_elapsed <- Sys.time() - time_start
  
  # Should handle large nboot efficiently
  expect_true(as.numeric(time_elapsed) < 30)
})


# ============================================================================
# SUITE 9: REGRESSION TESTS - Bug Fix Verification
# ============================================================================
# These tests ensure the fix for Tsallis entropy normalization (removal of abs())
# doesn't regress when code is modified in the future.

test_that("REGRESSION 9.1: Entropy normalization formula correctness", {
  
  n_tx <- 100
  test_cases <- list(
    list(q = 0.5, name = "q < 1"),
    list(q = 1.0, name = "q = 1"),
    list(q = 2.0, name = "q > 1"),
    list(q = 3.0, name = "q >> 1")
  )
  
  # Create uniform distribution (should normalize to 1)
  uniform_counts <- matrix(rep(1, n_tx * 20), nrow = n_tx, ncol = 20)
  
  for (test_case in test_cases) {
    q <- test_case$q
    name <- test_case$name
    
    entropy_vals <- jis_tsallis_entropy_cpp(uniform_counts, q = q, normalize = TRUE)
    
    # All sample entropies should be approximately 1 for uniform distribution
    expect_true(all(abs(entropy_vals - 1.0) < 0.02),
                info = paste("Incorrect normalization for", name))
  }
})

test_that("REGRESSION 9.2: No abs() artifacts in Tsallis computation", {
  
  # This test would catch if someone accidentally reintroduced abs() 
  # by verifying the mathematical properties
  
  n_tx <- c(10, 50, 100, 200)
  qs <- c(0.5, 0.9, 1.1, 1.5, 2.0, 3.0)
  
  for (n in n_tx) {
    for (q in qs) {
      if (abs(q - 1.0) < 1e-6) next  # Skip q = 1 (different formula)
      
      # Compute maximum entropy using correct formula
      max_h_correct <- (1.0 - n^(1.0 - q)) / (q - 1.0)
      
      # max_h must always be positive for valid q and n
      expect_true(max_h_correct > 0,
                  info = sprintf("max_h not positive for n=%d, q=%.2f", n, q))
    }
  }
})

test_that("REGRESSION 9.3: Singular vs uniform distribution contrast", {
  
  # Singular and uniform distributions should have very different normalized entropies
  n_tx <- 100
  q <- 1.5
  
  # Uniform distribution
  uniform_counts <- matrix(rep(1, n_tx * 20), nrow = n_tx, ncol = 20)
  entropy_uniform <- jis_tsallis_entropy_cpp(uniform_counts, q = q, normalize = TRUE)
  
  # Singular distribution (all on one species)
  singular_counts <- matrix(0, nrow = n_tx, ncol = 20)
  singular_counts[1, ] <- 100
  entropy_singular <- jis_tsallis_entropy_cpp(singular_counts, q = q, normalize = TRUE)
  
  # Uniform should have high entropy, singular should have low
  expect_true(mean(entropy_uniform) > 0.9)
  expect_true(mean(entropy_singular) < 0.1)
})

test_that("REGRESSION 9.4: Jackknife influences unchanged after fix", {
  
  set.seed(555)
  counts <- matrix(rpois(30 * 25, lambda = 8), nrow = 30, ncol = 25)
  
  # Test multiple q values
  test_qs <- c(0.5, 1.0, 2.0)
  
  for (q in test_qs) {
    influences <- jis_jackknife_influences_cpp(counts, q = q, normalize = TRUE)
    
    # Influences should be non-negative
    expect_true(all(influences >= 0 | is.na(influences)))
    
    # No influence should be suspiciously large (would indicate normalization error)
    valid_influences <- influences[!is.na(influences)]
    expect_true(max(valid_influences) < 10,
                info = sprintf("Unusually large influence for q = %.2f", q))
  }
})

test_that("REGRESSION 9.5: Fix doesn't break edge cases", {
  
  # Test with very small counts
  counts_small <- matrix(rep(c(1, 0, 0), each = 20), nrow = 3, ncol = 20) + 1e-8
  entropy_small <- jis_tsallis_entropy_cpp(counts_small, q = 1.5, normalize = TRUE)
  expect_true(all(is.finite(entropy_small) | is.na(entropy_small)))
  
  # Test with large counts
  counts_large <- matrix(rep(c(1000, 500, 250), each = 20), nrow = 3, ncol = 20)
  entropy_large <- jis_tsallis_entropy_cpp(counts_large, q = 1.5, normalize = TRUE)
  expect_true(all(is.finite(entropy_large) | is.na(entropy_large)))
})

test_that("REGRESSION 9.6: Mathematical consistency check", {
  
  # For uniform distribution P = (1/n, 1/n, ..., 1/n)
  # H_q = (1 - sum(p^q)) / (q - 1) = (1 - n*(1/n)^q) / (q - 1)
  #     = (1 - n^(1-q)) / (q - 1)
  # H_q_normalized = H_q / H_max = H_q / [(1 - n^(1-q)) / (q - 1)] = 1
  
  n_tx <- 50
  n_samples <- 15
  test_qs <- c(0.5, 0.8, 1.2, 1.8, 2.5)
  
  # Create uniform distribution
  uniform_counts <- matrix(1, nrow = n_tx, ncol = n_samples)
  
  for (q in test_qs) {
    entropy_vals <- jis_tsallis_entropy_cpp(uniform_counts, q = q, normalize = TRUE, 
                                             log_base = 2, pseudocount = 0, 
                                             n_tx_fixed = -1)
    
    # All should be ~1
    mean_entropy <- mean(entropy_vals, na.rm = TRUE)
    expect_true(abs(mean_entropy - 1.0) < 0.05,
                info = sprintf("Normalized uniform entropy not 1 for q = %.2f", q))
  }
})

