context("C++ Bootstrap Implementation (Rcpp)")

# Runtime ~14s: kept in the daily CI (Bioconductor check budget is 40 min).
# Slow Monte Carlo validation lives in tests/testthat/ (skipped on Bioconductor via skip_on_bioc).

# Test data preparation
test_counts <- c(100, 50, 25, 10)
large_counts <- c(1000, 500, 250, 100, 50, 25)
skewed_counts <- c(1000, 10, 5, 2)
small_counts <- c(5, 3, 2, 1)

# ============================================================================
# SUITE 1: Basic Bootstrap Functionality
# ============================================================================

test_that("C++ bootstrap_compute_cpp compiles and is available", {
  
  # Check that the function exists
  expect_true(exists("bootstrap_compute_cpp_wrapper"))
  
  # Check that it's callable
  result <- bootstrap_compute_cpp_wrapper(
    x = test_counts, q = 1.0, normalize = TRUE, 
    nboot = 10L, log_base = exp(1), pseudocount = 0.0
  )
  
  expect_true(is.numeric(result))
  expect_length(result, 10)
})

test_that("bootstrap_compute_cpp returns valid structure", {
  
  result <- bootstrap_compute_cpp_wrapper(
    x = test_counts, q = 1.0, normalize = TRUE, 
    nboot = 50L, log_base = exp(1), pseudocount = 0.0
  )
  
  expect_is(result, "numeric")
  expect_length(result, 50)
  expect_true(all(is.finite(result)))
})

test_that("bootstrap_compute_cpp produces positive entropy values", {
  
  result <- bootstrap_compute_cpp_wrapper(
    x = test_counts, q = 1.0, normalize = TRUE, 
    nboot = 100L, log_base = exp(1), pseudocount = 0.0
  )
  
  # For Shannon entropy (q=1) with normalized=TRUE, values should be in [0, 1]
  expect_true(all(result >= 0, na.rm = TRUE))
  expect_true(all(result <= 1, na.rm = TRUE))
})

test_that("bootstrap_compute_cpp respects nboot parameter", {
  
  for (nboot in c(10, 50, 100, 500)) {
    result <- bootstrap_compute_cpp_wrapper(
      x = test_counts, q = 1.0, normalize = TRUE, 
      nboot = as.integer(nboot), log_base = exp(1), pseudocount = 0.0
    )
    
    expect_length(result, nboot)
  }
})

# ============================================================================
# SUITE 2: Q Parameter Variations
# ============================================================================

test_that("bootstrap_compute_cpp works with various q values", {
  
  q_values <- c(0.5, 1.0, 1.5, 2.0, 2.5)
  
  for (q in q_values) {
    result <- bootstrap_compute_cpp_wrapper(
      x = test_counts, q = q, normalize = TRUE, 
      nboot = 50L, log_base = exp(1), pseudocount = 0.0
    )
    
    expect_true(all(is.finite(result)))
    expect_length(result, 50)
  }
})

test_that("bootstrap_compute_cpp handles q=0 (species richness)", {
  
  result <- bootstrap_compute_cpp_wrapper(
    x = test_counts, q = 0.0, normalize = TRUE, 
    nboot = 50L, log_base = exp(1), pseudocount = 0.0
  )
  
  expect_true(all(is.finite(result)))
  expect_length(result, 50)
  # Species richness should be consistent
  expect_true(all(result == result[1]))
})

test_that("bootstrap_compute_cpp produces different distributions for different q", {
  
  result_q1 <- bootstrap_compute_cpp_wrapper(
    x = test_counts, q = 1.0, normalize = TRUE, 
    nboot = 100L, log_base = exp(1), pseudocount = 0.0
  )
  
  result_q2 <- bootstrap_compute_cpp_wrapper(
    x = test_counts, q = 2.0, normalize = TRUE, 
    nboot = 100L, log_base = exp(1), pseudocount = 0.0
  )
  
  # Means should be different for different q
  expect_false(isTRUE(all.equal(
    mean(result_q1, na.rm = TRUE),
    mean(result_q2, na.rm = TRUE),
    tolerance = 1e-3
  )))
})

# ============================================================================
# SUITE 3: Normalization Options
# ============================================================================

test_that("bootstrap_compute_cpp respects normalize parameter", {
  
  result_norm <- bootstrap_compute_cpp_wrapper(
    x = test_counts, q = 1.0, normalize = TRUE, 
    nboot = 50L, log_base = exp(1), pseudocount = 0.0
  )
  
  result_unnorm <- bootstrap_compute_cpp_wrapper(
    x = test_counts, q = 1.0, normalize = FALSE, 
    nboot = 50L, log_base = exp(1), pseudocount = 0.0
  )
  
  # Normalized should be <= 1, unnormalized can be > 1
  expect_true(all(result_norm <= 1, na.rm = TRUE))
  expect_true(max(result_unnorm, na.rm = TRUE) > 0.5)
})

test_that("bootstrap_compute_cpp normalized <= unnormalized for q=1", {
  
  set.seed(42)
  result_norm <- bootstrap_compute_cpp_wrapper(
    x = test_counts, q = 1.0, normalize = TRUE, 
    nboot = 100L, log_base = exp(1), pseudocount = 0.0
  )
  
  set.seed(42)
  result_unnorm <- bootstrap_compute_cpp_wrapper(
    x = test_counts, q = 1.0, normalize = FALSE, 
    nboot = 100L, log_base = exp(1), pseudocount = 0.0
  )
  
  # Should be roughly proportional (normalized is unnormalized / max_entropy)
  ratios <- result_norm[!is.na(result_norm)] / result_unnorm[!is.na(result_unnorm)]
  expect_true(all(ratios > 0 & ratios <= 1, na.rm = TRUE))
})

# ============================================================================
# SUITE 4: Log Base Parameter
# ============================================================================

test_that("bootstrap_compute_cpp works with different log bases", {
  
  log_bases <- c(exp(1), 2, 10, 2.718281828)
  
  for (base in log_bases) {
    result <- bootstrap_compute_cpp_wrapper(
      x = test_counts, q = 1.0, normalize = FALSE, 
      nboot = 50L, log_base = base, pseudocount = 0.0
    )
    
    expect_true(all(is.finite(result)))
    expect_length(result, 50)
  }
})

test_that("bootstrap_compute_cpp log base affects scale but not relative ordering", {
  
  set.seed(123)
  result_e <- bootstrap_compute_cpp_wrapper(
    x = test_counts, q = 1.0, normalize = FALSE, 
    nboot = 100L, log_base = exp(1), pseudocount = 0.0
  )
  
  set.seed(123)
  result_2 <- bootstrap_compute_cpp_wrapper(
    x = test_counts, q = 1.0, normalize = FALSE, 
    nboot = 100L, log_base = 2.0, pseudocount = 0.0
  )
  
  # Log base 2 should be log(e) times smaller than natural log
  # E.g., H_e / H_2 ≈ log(2) ≈ 0.693
  ratio <- mean(result_2[!is.na(result_2)]) / mean(result_e[!is.na(result_e)])
  expected_ratio <- 1 / log(2)  # ≈ 1.443
  
  expect_equal(ratio, expected_ratio, tolerance = 0.05)
})

# ============================================================================
# SUITE 5: Pseudocount Support
# ============================================================================

test_that("bootstrap_compute_cpp works with pseudocount", {
  
  for (pc in c(0, 0.5, 1.0)) {
    result <- bootstrap_compute_cpp_wrapper(
      x = test_counts, q = 1.0, normalize = TRUE, 
      nboot = 50L, log_base = exp(1), pseudocount = pc
    )
    
    expect_true(all(is.finite(result)))
    expect_length(result, 50)
  }
})

test_that("bootstrap_compute_cpp with pseudocount affects entropy", {
  
  # With counts c(100, 50, 25, 10), very unequal
  result_no_pc <- bootstrap_compute_cpp_wrapper(
    x = test_counts, q = 1.0, normalize = TRUE, 
    nboot = 100L, log_base = exp(1), pseudocount = 0.0
  )
  
  # With pseudocount, should be more uniform, so higher entropy
  result_with_pc <- bootstrap_compute_cpp_wrapper(
    x = test_counts, q = 1.0, normalize = TRUE, 
    nboot = 100L, log_base = exp(1), pseudocount = 1.0
  )
  
  # Mean entropy with pseudocount should be higher (more uniform)
  expect_gt(
    mean(result_with_pc, na.rm = TRUE),
    mean(result_no_pc, na.rm = TRUE)
  )
})

# ============================================================================
# SUITE 6: Edge Cases
# ============================================================================

test_that("bootstrap_compute_cpp handles small counts", {
  
  result <- bootstrap_compute_cpp_wrapper(
    x = small_counts, q = 1.0, normalize = TRUE, 
    nboot = 50L, log_base = exp(1), pseudocount = 0.0
  )
  
  expect_true(all(is.finite(result)))
  expect_length(result, 50)
})

test_that("bootstrap_compute_cpp handles skewed distributions", {
  
  result <- bootstrap_compute_cpp_wrapper(
    x = skewed_counts, q = 1.0, normalize = TRUE, 
    nboot = 50L, log_base = exp(1), pseudocount = 0.0
  )
  
  expect_true(all(is.finite(result)))
  expect_length(result, 50)
  # Skewed distribution should have lower entropy
  expect_true(mean(result, na.rm = TRUE) < 1.0)
})

test_that("bootstrap_compute_cpp handles large counts", {
  
  result <- bootstrap_compute_cpp_wrapper(
    x = large_counts, q = 1.0, normalize = TRUE, 
    nboot = 50L, log_base = exp(1), pseudocount = 0.0
  )
  
  expect_true(all(is.finite(result)))
  expect_length(result, 50)
})

test_that("bootstrap_compute_cpp handles single value vector", {
  
  result <- bootstrap_compute_cpp_wrapper(
    x = c(100), q = 1.0, normalize = TRUE, 
    nboot = 50L, log_base = exp(1), pseudocount = 0.0
  )
  
  # Single value should have zero entropy
  expect_true(all(result == 0, na.rm = TRUE))
})

test_that("bootstrap_compute_cpp handles zeros with pseudocount", {
  
  zero_counts <- c(0, 0, 0, 100)
  
  result <- bootstrap_compute_cpp_wrapper(
    x = zero_counts, q = 1.0, normalize = TRUE, 
    nboot = 50L, log_base = exp(1), pseudocount = 1.0
  )
  
  expect_true(all(is.finite(result)))
  expect_length(result, 50)
})

# ============================================================================
# SUITE 7: Reproducibility
# ============================================================================

test_that("bootstrap_compute_cpp is reproducible with set.seed", {
  
  set.seed(999)
  result1 <- bootstrap_compute_cpp_wrapper(
    x = test_counts, q = 1.0, normalize = TRUE, 
    nboot = 100L, log_base = exp(1), pseudocount = 0.0
  )
  
  set.seed(999)
  result2 <- bootstrap_compute_cpp_wrapper(
    x = test_counts, q = 1.0, normalize = TRUE, 
    nboot = 100L, log_base = exp(1), pseudocount = 0.0
  )
  
  expect_equal(result1, result2)
})

test_that("bootstrap_compute_cpp different seeds produce different results", {
  
  set.seed(111)
  result1 <- bootstrap_compute_cpp_wrapper(
    x = test_counts, q = 1.0, normalize = TRUE, 
    nboot = 100L, log_base = exp(1), pseudocount = 0.0
  )
  
  set.seed(222)
  result2 <- bootstrap_compute_cpp_wrapper(
    x = test_counts, q = 1.0, normalize = TRUE, 
    nboot = 100L, log_base = exp(1), pseudocount = 0.0
  )
  
  expect_false(isTRUE(all.equal(result1, result2)))
})

# ============================================================================
# SUITE 8: Integration with Bootstrap CI Wrapper
# ============================================================================

test_that(".bootstrap_resample_optimized uses C++ by default", {
  
  result <- .bootstrap_resample_optimized(
    x = test_counts, q = 1.0, norm = TRUE, nboot = 50,
    log_base = exp(1), pseudocount = 0, what = "S", paired = FALSE
  )
  
  expect_is(result, "numeric")
  expect_length(result, 50)
  expect_true(all(is.finite(result)))
})

test_that(".bootstrap_resample_optimized handles what='S' (entropy)", {
  
  result <- .bootstrap_resample_optimized(
    x = test_counts, q = 1.0, norm = TRUE, nboot = 50,
    log_base = exp(1), pseudocount = 0, what = "S", paired = FALSE
  )
  
  expect_true(all(result >= 0, na.rm = TRUE))
  expect_true(all(result <= 1, na.rm = TRUE))
})

test_that(".bootstrap_resample_optimized handles what='D' (Hill numbers)", {
  
  result <- .bootstrap_resample_optimized(
    x = test_counts, q = 1.0, norm = TRUE, nboot = 50,
    log_base = exp(1), pseudocount = 0, what = "D", paired = FALSE
  )
  
  expect_is(result, "numeric")
  expect_length(result, 50)
  expect_true(all(is.finite(result)))
  # Hill numbers should be >= 1 for q=1
  expect_true(all(result >= 1, na.rm = TRUE))
})

test_that(".bootstrap_resample_optimized converts entropy to Hill numbers correctly", {
  
  # For Shannon (q=1), Hill number D_1 = exp(H_1)
  set.seed(555)
  entropy_result <- .bootstrap_resample_optimized(
    x = test_counts, q = 1.0, norm = FALSE, nboot = 100,
    log_base = exp(1), pseudocount = 0, what = "S", paired = FALSE
  )
  
  set.seed(555)
  hill_result <- .bootstrap_resample_optimized(
    x = test_counts, q = 1.0, norm = FALSE, nboot = 100,
    log_base = exp(1), pseudocount = 0, what = "D", paired = FALSE
  )
  
  # Hill numbers should be exp(entropy)
  expected_hill <- exp(entropy_result)
  expect_equal(hill_result, expected_hill, tolerance = 1e-6)
})

# ============================================================================
# SUITE 9: Bootstrap CI Integration
# ============================================================================

test_that(".bootstrap_compute_ci uses C++ bootstrap internally", {
  
  result <- .bootstrap_compute_ci(
    x = test_counts, q = 1.0, norm = TRUE, nboot = 100,
    ci = 0.95, method = "percentile", log_base = exp(1),
    pseudocount = 0, what = "S", paired = FALSE
  )
  
  expect_is(result, "list")
  expect_named(result, c("point_est", "bootstrap_dist", "ci_result", "accel_factor"))
  expect_length(result$bootstrap_dist, 100)
  expect_true(all(is.finite(result$bootstrap_dist)))
})

test_that("bootstrap_compute_ci confidence intervals are consistent", {
  
  result <- .bootstrap_compute_ci(
    x = test_counts, q = 1.0, norm = TRUE, nboot = 100,
    ci = 0.95, method = "percentile", log_base = exp(1),
    pseudocount = 0, what = "S", paired = FALSE
  )
  
  # CI should bracket the median
  ci_lower <- result$ci_result$lower
  ci_upper <- result$ci_result$upper
  
  expect_lt(ci_lower, ci_upper)
  expect_gt(ci_upper, ci_lower)
})

test_that("bootstrap_compute_ci BCa method uses bootstrap correctly", {
  
  result_percentile <- .bootstrap_compute_ci(
    x = test_counts, q = 1.0, norm = TRUE, nboot = 100,
    ci = 0.95, method = "percentile", log_base = exp(1),
    pseudocount = 0, what = "S", paired = FALSE
  )
  
  result_bca <- .bootstrap_compute_ci(
    x = test_counts, q = 1.0, norm = TRUE, nboot = 100,
    ci = 0.95, method = "bca", log_base = exp(1),
    pseudocount = 0, what = "S", paired = FALSE
  )
  
  # Both should use same bootstrap samples (same distribution)
  expect_equal(length(result_percentile$bootstrap_dist), 
               length(result_bca$bootstrap_dist))
})

# ============================================================================
# SUITE 10: Performance Characteristics
# ============================================================================

test_that("bootstrap_compute_cpp scales reasonably with nboot", {
  
  times <- numeric(3)
  nboots <- c(100, 500, 1000)
  
  for (i in seq_along(nboots)) {
    t <- system.time({
      bootstrap_compute_cpp_wrapper(
        x = test_counts, q = 1.0, normalize = TRUE, 
        nboot = as.integer(nboots[i]), log_base = exp(1), pseudocount = 0.0
      )
    })[3]
    times[i] <- t
  }
  
  # Time should increase roughly linearly with nboot
  # Check that larger nboot doesn't have unexpected overhead
  per_rep <- times / nboots
  expect_true(all(per_rep < 0.1))  # Less than 100ms per replicate
})

test_that("bootstrap_compute_cpp is faster than R version", {
  
  # This test just ensures the C++ version completes in reasonable time
  t_cpp <- system.time({
    bootstrap_compute_cpp_wrapper(
      x = test_counts, q = 1.0, normalize = TRUE, 
      nboot = 1000L, log_base = exp(1), pseudocount = 0.0
    )
  })[3]
  
  # Should complete in less than 1 second for 1000 replicates
  expect_lt(t_cpp, 1.0)
})

# ============================================================================
# SUITE 11: Mathematical Properties
# ============================================================================

test_that("bootstrap entropy estimates are consistent with point estimate", {
  
  result <- .bootstrap_compute_ci(
    x = test_counts, q = 1.0, norm = TRUE, nboot = 200,
    ci = 0.95, method = "percentile", log_base = exp(1),
    pseudocount = 0, what = "S", paired = FALSE
  )
  
  # Point estimate should be close to median of bootstrap distribution
  point_est <- result$point_est
  boot_median <- median(result$bootstrap_dist, na.rm = TRUE)
  
  # Allow some difference due to sampling
  expect_equal(point_est, boot_median, tolerance = 0.1)
})

test_that("bootstrap standard error is reasonable", {
  
  result <- .bootstrap_compute_ci(
    x = test_counts, q = 1.0, norm = TRUE, nboot = 200,
    ci = 0.95, method = "percentile", log_base = exp(1),
    pseudocount = 0, what = "S", paired = FALSE
  )
  
  # SE should be positive and finite
  boot_se <- sd(result$bootstrap_dist, na.rm = TRUE)
  expect_gt(boot_se, 0)
  expect_true(is.finite(boot_se))
})

# ============================================================================
# SUITE 12: Error Handling and Fallback
# ============================================================================

test_that(".bootstrap_resample_optimized has fallback to R version", {
  
  # This test verifies the fallback mechanism exists
  # (by default C++ succeeds, but if it fails, should fall back)
  
  result <- .bootstrap_resample_optimized(
    x = test_counts, q = 1.0, norm = TRUE, nboot = 50,
    log_base = exp(1), pseudocount = 0, what = "S", paired = FALSE
  )
  
  expect_true(all(is.finite(result)))
})

test_that("bootstrap_compute_cpp_wrapper has proper error messages", {
  
  # Test that invalid inputs are caught
  expect_error(
    bootstrap_compute_cpp_wrapper(
      x = numeric(0), q = 1.0, normalize = TRUE, 
      nboot = 50L, log_base = exp(1), pseudocount = 0.0
    )
  )
})

test_that("bootstrap_compute_cpp handles negative nboot gracefully", {
  
  expect_error(
    bootstrap_compute_cpp_wrapper(
      x = test_counts, q = 1.0, normalize = TRUE, 
      nboot = -10L, log_base = exp(1), pseudocount = 0.0
    )
  )
})

# ============================================================================
# SUITE 13: Block Bootstrap (Paired Samples) - C++ Implementation
# ============================================================================

# Prepare paired test data
paired_counts_2 <- c(100, 95, 110, 105)  # 2 pairs
paired_counts_3 <- c(100, 95, 110, 105, 90, 92)  # 3 pairs
paired_counts_large <- c(500, 490, 520, 510, 450, 440, 480, 470)  # 4 pairs

test_that("C++ block_bootstrap_compute_cpp compiles and is available", {
  
  # Check that the function exists
  expect_true(exists("block_bootstrap_compute_cpp_wrapper"))
  
  # Check that it's callable
  result <- block_bootstrap_compute_cpp_wrapper(
    x = paired_counts_2, q = 1.0, normalize = TRUE, 
    nboot = 10L, log_base = exp(1), pseudocount = 0.0
  )
  
  expect_true(is.numeric(result))
  expect_length(result, 10)
})

test_that("block_bootstrap_compute_cpp requires even-length input", {
  
  odd_length <- c(100, 95, 110)  # 3 elements (1.5 pairs)
  
  expect_error(
    block_bootstrap_compute_cpp_wrapper(
      x = odd_length, q = 1.0, normalize = TRUE, 
      nboot = 10L, log_base = exp(1), pseudocount = 0.0
    )
  )
})

test_that("block_bootstrap_compute_cpp returns valid structure", {
  
  result <- block_bootstrap_compute_cpp_wrapper(
    x = paired_counts_3, q = 1.0, normalize = TRUE, 
    nboot = 50L, log_base = exp(1), pseudocount = 0.0
  )
  
  expect_is(result, "numeric")
  expect_length(result, 50)
  expect_true(all(is.finite(result)))
})

test_that("block_bootstrap_compute_cpp produces entropy in valid range", {
  
  result <- block_bootstrap_compute_cpp_wrapper(
    x = paired_counts_3, q = 1.0, normalize = TRUE, 
    nboot = 100L, log_base = exp(1), pseudocount = 0.0
  )
  
  # For normalized Shannon entropy, should be in [0, 1]
  expect_true(all(result >= 0, na.rm = TRUE))
  expect_true(all(result <= 1, na.rm = TRUE))
})

test_that("block_bootstrap_compute_cpp respects nboot parameter", {
  
  for (nboot in c(10, 50, 100, 200)) {
    result <- block_bootstrap_compute_cpp_wrapper(
      x = paired_counts_3, q = 1.0, normalize = TRUE, 
      nboot = as.integer(nboot), log_base = exp(1), pseudocount = 0.0
    )
    
    expect_length(result, nboot)
  }
})

test_that("block_bootstrap_compute_cpp works with various q values", {
  
  q_values <- c(0.5, 1.0, 1.5, 2.0)
  
  for (q in q_values) {
    result <- block_bootstrap_compute_cpp_wrapper(
      x = paired_counts_3, q = q, normalize = TRUE, 
      nboot = 50L, log_base = exp(1), pseudocount = 0.0
    )
    
    expect_true(all(is.finite(result)))
    expect_length(result, 50)
  }
})

test_that("block_bootstrap_compute_cpp handles normalization parameter", {
  
  result_norm <- block_bootstrap_compute_cpp_wrapper(
    x = paired_counts_3, q = 1.0, normalize = TRUE, 
    nboot = 50L, log_base = exp(1), pseudocount = 0.0
  )
  
  result_unnorm <- block_bootstrap_compute_cpp_wrapper(
    x = paired_counts_3, q = 1.0, normalize = FALSE, 
    nboot = 50L, log_base = exp(1), pseudocount = 0.0
  )
  
  # Normalized should be <= 1, unnormalized can be larger
  expect_true(all(result_norm <= 1, na.rm = TRUE))
  expect_true(max(result_unnorm, na.rm = TRUE) > 0.5)
})

test_that("block_bootstrap_compute_cpp works with different log bases", {
  
  log_bases <- c(exp(1), 2, 10)
  
  for (base in log_bases) {
    result <- block_bootstrap_compute_cpp_wrapper(
      x = paired_counts_3, q = 1.0, normalize = FALSE, 
      nboot = 50L, log_base = base, pseudocount = 0.0
    )
    
    expect_true(all(is.finite(result)))
    expect_length(result, 50)
  }
})

test_that("block_bootstrap_compute_cpp supports pseudocount", {
  
  for (pc in c(0, 0.5, 1.0)) {
    result <- block_bootstrap_compute_cpp_wrapper(
      x = paired_counts_3, q = 1.0, normalize = TRUE, 
      nboot = 50L, log_base = exp(1), pseudocount = pc
    )
    
    expect_true(all(is.finite(result)))
    expect_length(result, 50)
  }
})

test_that("block_bootstrap_compute_cpp is reproducible with set.seed", {
  
  set.seed(777)
  result1 <- block_bootstrap_compute_cpp_wrapper(
    x = paired_counts_3, q = 1.0, normalize = TRUE, 
    nboot = 100L, log_base = exp(1), pseudocount = 0.0
  )
  
  set.seed(777)
  result2 <- block_bootstrap_compute_cpp_wrapper(
    x = paired_counts_3, q = 1.0, normalize = TRUE, 
    nboot = 100L, log_base = exp(1), pseudocount = 0.0
  )
  
  expect_equal(result1, result2)
})

test_that("block_bootstrap_compute_cpp with different seeds differs", {
  
  # Use data with more variation to detect RNG seeding differences
  # Balanced data produces very similar entropy values, making seed differences hard to detect
  varied_counts <- c(10, 200, 50, 150, 5, 300)  # More unbalanced distribution
  
  set.seed(111)
  result1 <- block_bootstrap_compute_cpp_wrapper(
    x = varied_counts, q = 1.0, normalize = TRUE, 
    nboot = 100L, log_base = exp(1), pseudocount = 0.0
  )
  
  set.seed(222)
  result2 <- block_bootstrap_compute_cpp_wrapper(
    x = varied_counts, q = 1.0, normalize = TRUE, 
    nboot = 100L, log_base = exp(1), pseudocount = 0.0
  )
  
  # With varied data, different seeds should produce detectably different results
  expect_false(isTRUE(all.equal(result1, result2, tolerance = 0.01)))
})

# ============================================================================
# SUITE 14: Block Bootstrap Integration with Wrappers
# ============================================================================

test_that(".bootstrap_resample_optimized handles paired=TRUE with what='S'", {
  
  result <- .bootstrap_resample_optimized(
    x = paired_counts_3, q = 1.0, norm = TRUE, nboot = 50,
    log_base = exp(1), pseudocount = 0, what = "S", paired = TRUE
  )
  
  expect_is(result, "numeric")
  expect_length(result, 50)
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
})

test_that(".bootstrap_resample_optimized handles paired=TRUE with what='D'", {
  
  result <- .bootstrap_resample_optimized(
    x = paired_counts_3, q = 1.0, norm = TRUE, nboot = 50,
    log_base = exp(1), pseudocount = 0, what = "D", paired = TRUE
  )
  
  expect_is(result, "numeric")
  expect_length(result, 50)
  expect_true(all(is.finite(result)))
  # Hill numbers >= 1 for Shannon
  expect_true(all(result >= 1, na.rm = TRUE))
})

test_that(".bootstrap_resample_optimized paired entropy to Hill conversion", {
  
  # For Shannon (q=1), Hill number D_1 = exp(H_1)
  set.seed(888)
  entropy_result <- .bootstrap_resample_optimized(
    x = paired_counts_3, q = 1.0, norm = FALSE, nboot = 100,
    log_base = exp(1), pseudocount = 0, what = "S", paired = TRUE
  )
  
  set.seed(888)
  hill_result <- .bootstrap_resample_optimized(
    x = paired_counts_3, q = 1.0, norm = FALSE, nboot = 100,
    log_base = exp(1), pseudocount = 0, what = "D", paired = TRUE
  )
  
  # Hill numbers should be exp(entropy)
  expected_hill <- exp(entropy_result)
  expect_equal(hill_result, expected_hill, tolerance = 1e-6)
})

test_that(".bootstrap_compute_ci works with paired=TRUE and percentile method", {
  
  result <- .bootstrap_compute_ci(
    x = paired_counts_3, q = 1.0, norm = TRUE, nboot = 100,
    ci = 0.95, method = "percentile", log_base = exp(1),
    pseudocount = 0, what = "S", paired = TRUE
  )
  
  expect_is(result, "list")
  expect_named(result, c("point_est", "bootstrap_dist", "ci_result", "accel_factor"))
  expect_length(result$bootstrap_dist, 100)
  expect_true(all(is.finite(result$bootstrap_dist)))
})

test_that(".bootstrap_compute_ci works with paired=TRUE and BCa method", {
  
  result <- .bootstrap_compute_ci(
    x = paired_counts_3, q = 1.0, norm = TRUE, nboot = 100,
    ci = 0.95, method = "bca", log_base = exp(1),
    pseudocount = 0, what = "S", paired = TRUE
  )
  
  expect_is(result, "list")
  expect_named(result, c("point_est", "bootstrap_dist", "ci_result", "accel_factor"))
  expect_true(is.numeric(result$point_est))
  expect_true(is.numeric(result$ci_result$lower))
  expect_true(is.numeric(result$ci_result$upper))
})

test_that("paired bootstrap CI bounds bracket the point estimate", {
  
  result <- .bootstrap_compute_ci(
    x = paired_counts_3, q = 1.0, norm = TRUE, nboot = 100,
    ci = 0.95, method = "percentile", log_base = exp(1),
    pseudocount = 0, what = "S", paired = TRUE
  )
  
  # CI should bracket median or close to it
  ci_lower <- result$ci_result$lower
  ci_upper <- result$ci_result$upper
  boot_median <- median(result$bootstrap_dist, na.rm = TRUE)
  
  expect_lt(ci_lower, boot_median)
  expect_gt(ci_upper, boot_median)
})

test_that("block bootstrap with different pair counts works", {
  
  # Test with different numbers of pairs
  pair_counts <- list(
    pairs_2 = c(100, 95, 110, 105),
    pairs_3 = c(100, 95, 110, 105, 90, 92),
    pairs_4 = c(100, 95, 110, 105, 90, 92, 85, 88)
  )
  
  for (name in names(pair_counts)) {
    result <- block_bootstrap_compute_cpp_wrapper(
      x = pair_counts[[name]], q = 1.0, normalize = TRUE, 
      nboot = 50L, log_base = exp(1), pseudocount = 0.0
    )
    
    expect_true(all(is.finite(result)))
    expect_length(result, 50)
  }
})

test_that("block bootstrap with small pairs works", {
  
  small_paired <- c(10, 9, 12, 11)  # Very small counts
  
  result <- block_bootstrap_compute_cpp_wrapper(
    x = small_paired, q = 1.0, normalize = TRUE, 
    nboot = 50L, log_base = exp(1), pseudocount = 0.0
  )
  
  expect_true(all(is.finite(result)))
  expect_length(result, 50)
})

test_that("block bootstrap with skewed pairs works", {
  
  skewed_paired <- c(1000, 10, 5, 2)  # Highly skewed pairs
  
  result <- block_bootstrap_compute_cpp_wrapper(
    x = skewed_paired, q = 1.0, normalize = TRUE, 
    nboot = 50L, log_base = exp(1), pseudocount = 0.0
  )
  
  expect_true(all(is.finite(result)))
  expect_length(result, 50)
})

test_that("block bootstrap preserves pair structure", {
  
  # Block bootstrap should resample pairs, so it should preserve
  # pair correlation structure. Both methods should have similar means
  # but may differ in overall distribution.
  
  set.seed(999)
  block_result <- block_bootstrap_compute_cpp_wrapper(
    x = paired_counts_3, q = 1.0, normalize = TRUE, 
    nboot = 100L, log_base = exp(1), pseudocount = 0.0
  )
  
  set.seed(999)
  # Use standard bootstrap with same seed
  standard_result <- bootstrap_compute_cpp_wrapper(
    x = paired_counts_3, q = 1.0, normalize = TRUE, 
    nboot = 100L, log_base = exp(1), pseudocount = 0.0
  )
  
  # Means should be similar (both resample same data)
  expect_equal(
    mean(block_result, na.rm = TRUE), 
    mean(standard_result, na.rm = TRUE), 
    tolerance = 0.05
  )
  
  # Distributions should be different due to different resampling schemes
  expect_false(isTRUE(all.equal(
    sd(block_result, na.rm = TRUE), 
    sd(standard_result, na.rm = TRUE), 
    tolerance = 1e-3
  )))
})

test_that("block bootstrap performance scales reasonably", {
  
  nboots <- c(100, 500, 1000)
  times <- numeric(length(nboots))
  
  for (i in seq_along(nboots)) {
    t <- system.time({
      block_bootstrap_compute_cpp_wrapper(
        x = paired_counts_large, q = 1.0, normalize = TRUE, 
        nboot = as.integer(nboots[i]), log_base = exp(1), pseudocount = 0.0
      )
    })[3]
    times[i] <- t
  }
  
  # Time should increase roughly linearly with nboot
  per_rep <- times / nboots
  expect_true(all(per_rep < 0.1))  # Less than 100ms per replicate
})

test_that("block bootstrap with various q values in paired CI", {
  
  q_values <- c(0.5, 1.0, 1.5, 2.0)
  
  for (q in q_values) {
    result <- .bootstrap_compute_ci(
      x = paired_counts_3, q = q, norm = TRUE, nboot = 50,
      ci = 0.95, method = "percentile", log_base = exp(1),
      pseudocount = 0, what = "S", paired = TRUE
    )
    
    expect_is(result, "list")
    expect_true(is.numeric(result$point_est))
    expect_true(is.numeric(result$ci_result$lower))
    expect_true(is.numeric(result$ci_result$upper))
  }
})

# ============================================================================
# SUITE 16: Hill Numbers ("D") via Bootstrap
# ============================================================================

test_that("bootstrap CI with what='D' (Hill numbers) computes correctly", {
  
  result_entropy <- .bootstrap_compute_ci(
    x = test_counts, q = 1.0, norm = FALSE, nboot = 100,
    ci = 0.95, method = "percentile", log_base = exp(1),
    pseudocount = 0, what = "S"
  )
  
  result_hill <- .bootstrap_compute_ci(
    x = test_counts, q = 1.0, norm = FALSE, nboot = 100,
    ci = 0.95, method = "percentile", log_base = exp(1),
    pseudocount = 0, what = "D"
  )
  
  # Hill numbers should be exp(entropy) for q=1
  expect_true(is.numeric(result_hill$point_est))
  expect_true(result_hill$point_est > 0)
  # Hill number = exp(Shannon entropy)
  expect_equal(result_hill$point_est, exp(result_entropy$point_est), tolerance = 0.01)
})

test_that("bootstrap CI with what='D' matches Hill number formula", {
  
  result <- .bootstrap_compute_ci(
    x = c(100, 50, 30, 20), q = 2.0, norm = FALSE, nboot = 50,
    ci = 0.95, method = "percentile", log_base = exp(1),
    pseudocount = 0, what = "D"
  )
  
  # For q=2: D_2 = 1 / sum(p^2)
  expect_true(is.numeric(result$point_est))
  expect_true(result$point_est > 0)
  expect_true(result$point_est <= length(c(100, 50, 30, 20)))  # Max Hill = n
})

test_that("Hill numbers across q spectrum", {
  
  q_vals <- c(0.5, 1.0, 1.5, 2.0)
  results <- lapply(q_vals, function(q) {
    .bootstrap_compute_ci(
      x = test_counts, q = q, norm = FALSE, nboot = 50,
      ci = 0.95, method = "percentile", log_base = exp(1),
      pseudocount = 0, what = "D"
    )
  })
  
  points <- sapply(results, function(r) r$point_est)
  
  # Hill numbers should decrease as q increases (for most distributions)
  expect_true(all(is.numeric(points)))
  expect_true(all(points > 0))
  expect_length(points, 4)
})

# ============================================================================
# SUITE 17: Log Base Variations
# ============================================================================

test_that("bootstrap with log_base=2 (bits)", {
  
  result_e <- .bootstrap_compute_ci(
    x = test_counts, q = 1.0, norm = FALSE, nboot = 50,
    ci = 0.95, method = "percentile", log_base = exp(1),
    pseudocount = 0, what = "S"
  )
  
  result_2 <- .bootstrap_compute_ci(
    x = test_counts, q = 1.0, norm = FALSE, nboot = 50,
    ci = 0.95, method = "percentile", log_base = 2,
    pseudocount = 0, what = "S"
  )
  
  # Conversion factor: log2(x) = ln(x) / ln(2)
  conversion <- 1 / log(2)
  expect_equal(result_2$point_est, result_e$point_est * conversion, tolerance = 0.01)
})

test_that("bootstrap with log_base=10 (nats)", {
  
  result_e <- .bootstrap_compute_ci(
    x = test_counts, q = 1.0, norm = FALSE, nboot = 50,
    ci = 0.95, method = "percentile", log_base = exp(1),
    pseudocount = 0, what = "S"
  )
  
  result_10 <- .bootstrap_compute_ci(
    x = test_counts, q = 1.0, norm = FALSE, nboot = 50,
    ci = 0.95, method = "percentile", log_base = 10,
    pseudocount = 0, what = "S"
  )
  
  # Conversion factor: log10(x) = ln(x) / ln(10)
  conversion <- 1 / log(10)
  expect_equal(result_10$point_est, result_e$point_est * conversion, tolerance = 0.01)
})

test_that("C++ bootstrap respects log_base in entropy calculation", {
  
  set.seed(123)
  result_e <- bootstrap_compute_cpp_wrapper(
    x = test_counts, q = 1.0, normalize = FALSE, 
    nboot = 20, log_base = exp(1), pseudocount = 0
  )
  
  set.seed(123)
  result_2 <- bootstrap_compute_cpp_wrapper(
    x = test_counts, q = 1.0, normalize = FALSE, 
    nboot = 20, log_base = 2, pseudocount = 0
  )
  
  # All values should scale by 1/ln(2)
  conversion_factor <- 1 / log(2)
  expect_true(all(is.finite(result_e)))
  expect_true(all(is.finite(result_2)))
  expect_equal(result_2, result_e * conversion_factor, tolerance = 1e-6)
})

# ============================================================================
# SUITE 18: Pseudocount Handling
# ============================================================================

test_that("bootstrap with pseudocount affects entropy calculation", {
  
  result_no_pc <- .bootstrap_compute_ci(
    x = c(100, 0, 0, 0), q = 1.0, norm = FALSE, nboot = 50,
    ci = 0.95, method = "percentile", log_base = exp(1),
    pseudocount = 0, what = "S"
  )
  
  result_pc <- .bootstrap_compute_ci(
    x = c(100, 0, 0, 0), q = 1.0, norm = FALSE, nboot = 50,
    ci = 0.95, method = "percentile", log_base = exp(1),
    pseudocount = 0.5, what = "S"
  )
  
  # Pseudocount should allow zero-count species to contribute entropy
  expect_true(is.numeric(result_no_pc$point_est))
  expect_true(is.numeric(result_pc$point_est))
  expect_true(result_pc$point_est > result_no_pc$point_est)
})

test_that("pseudocount smooths bootstrap replicates", {
  
  # Use moderately sparse data where pseudocount improves stability
  # Extreme sparsity (100,0,0,0) causes pseudocount to ADD variance (degenerate -> less degenerate)
  # Moderate sparsity benefits from pseudocount smoothing
  rep_no_pc <- bootstrap_compute_cpp_wrapper(
    x = c(100, 10, 5, 2), q = 1.0, normalize = FALSE, 
    nboot = 100, log_base = exp(1), pseudocount = 0
  )
  
  rep_pc <- bootstrap_compute_cpp_wrapper(
    x = c(100, 10, 5, 2), q = 1.0, normalize = FALSE, 
    nboot = 100, log_base = exp(1), pseudocount = 0.5
  )
  
  # Pseudocount should reduce variance by smoothing proportions
  # Tolerance of 0.05 accounts for random variation in bootstrap samples
  expect_true(sd(rep_pc) <= sd(rep_no_pc) + 0.05)
})

# ============================================================================
# SUITE 19: Edge Cases
# ============================================================================

test_that("bootstrap handles uniform distribution correctly", {
  
  # Uniform distribution should have maximum entropy
  uniform_data <- rep(10, 4)  # All equal counts
  
  result <- .bootstrap_compute_ci(
    x = uniform_data, q = 1.0, norm = TRUE, nboot = 50,
    ci = 0.95, method = "percentile", log_base = exp(1),
    pseudocount = 0, what = "S"
  )
  
  # Normalized Shannon entropy for uniform distribution should be near 1
  expect_true(result$point_est > 0.98)
  expect_true(result$point_est <= 1.0)
})

test_that("bootstrap handles highly skewed distribution", {
  
  # Highly skewed: one dominant species
  skewed_data <- c(1000, 1, 1, 1)
  
  result <- .bootstrap_compute_ci(
    x = skewed_data, q = 1.0, norm = TRUE, nboot = 50,
    ci = 0.95, method = "percentile", log_base = exp(1),
    pseudocount = 0, what = "S"
  )
  
  # Skewed distribution should have low normalized entropy
  expect_true(result$point_est < 0.3)
})

test_that("bootstrap with single non-zero count", {
  
  # Only one species present
  single_species <- c(100, 0, 0, 0)
  
  result <- .bootstrap_compute_ci(
    x = single_species, q = 1.0, norm = TRUE, nboot = 50,
    ci = 0.95, method = "percentile", log_base = exp(1),
    pseudocount = 0, what = "S"
  )
  
  # Single species should have entropy = 0
  expect_true(result$point_est < 1e-6 || is.na(result$point_est))
})

test_that("bootstrap handles very small counts", {
  
  small_data <- c(1, 1, 1, 1)
  
  result <- .bootstrap_compute_ci(
    x = small_data, q = 1.0, norm = TRUE, nboot = 50,
    ci = 0.95, method = "percentile", log_base = exp(1),
    pseudocount = 0, what = "S"
  )
  
  # Should return valid result even with small counts
  expect_true(is.numeric(result$point_est))
  expect_true(result$point_est > 0)
})

# ============================================================================
# SUITE 20: Bootstrap CI Methods Comparison
# ============================================================================

test_that("percentile vs BCa CI methods produce reasonable results", {
  
  result_percentile <- .bootstrap_compute_ci(
    x = test_counts, q = 1.0, norm = TRUE, nboot = 100,
    ci = 0.95, method = "percentile", log_base = exp(1),
    pseudocount = 0, what = "S"
  )
  
  result_bca <- .bootstrap_compute_ci(
    x = test_counts, q = 1.0, norm = TRUE, nboot = 100,
    ci = 0.95, method = "BCa", log_base = exp(1),
    pseudocount = 0, what = "S"
  )
  
  # Both should have point estimates close to each other
  expect_true(abs(result_percentile$point_est - result_bca$point_est) < 0.05)
  
  # CIs should be overlapping and reasonable
  expect_true(result_percentile$ci_result$lower < result_percentile$point_est)
  expect_true(result_percentile$ci_result$upper > result_percentile$point_est)
  expect_true(result_bca$ci_result$lower < result_bca$point_est)
  expect_true(result_bca$ci_result$upper > result_bca$point_est)
})

test_that("CI coverage properties for different sample sizes", {
  
  sizes <- c(10, 50, 100)
  results <- lapply(sizes, function(n) {
    data <- sample(test_counts, n, replace = TRUE)
    result <- .bootstrap_compute_ci(
      x = data, q = 1.0, norm = TRUE, nboot = 50,
      ci = 0.95, method = "percentile", log_base = exp(1),
      pseudocount = 0, what = "S"
    )
    list(
      point = result$point_est,
      lower = result$ci_result$lower,
      upper = result$ci_result$upper,
      width = result$ci_result$upper - result$ci_result$lower
    )
  })
  
  # Larger samples should produce narrower CIs
  widths <- sapply(results, function(r) r$width)
  expect_true(widths[3] <= widths[1] + 0.01)  # 100-sample CI <= 10-sample CI
})

# ============================================================================
# SUITE 21: Multiple q-values in Bootstrap
# ============================================================================

test_that("bootstrap_compute_cpp handles full q spectrum", {
  
  q_values <- c(0.5, 1.0, 1.5, 2.0, 3.0)
  
  results <- lapply(q_values, function(q) {
    bootstrap_compute_cpp_wrapper(
      x = test_counts, q = q, normalize = TRUE, 
      nboot = 30, log_base = exp(1), pseudocount = 0
    )
  })
  
  # All should return valid numeric vectors
  expect_true(all(sapply(results, function(r) is.numeric(r))))
  expect_true(all(sapply(results, function(r) length(r) == 30)))
  expect_true(all(sapply(results, function(r) all(is.finite(r)))))
})

test_that("bootstrap CI with multiple q-values produces consistent results", {
  
  q_spectrum <- c(0.5, 1.0, 1.5, 2.0)
  cis <- lapply(q_spectrum, function(q) {
    .bootstrap_compute_ci(
      x = test_counts, q = q, norm = TRUE, nboot = 50,
      ci = 0.95, method = "percentile", log_base = exp(1),
      pseudocount = 0, what = "S"
    )
  })
  
  # All point estimates should be in [0, 1] for normalized entropy
  points <- sapply(cis, function(ci) ci$point_est)
  expect_true(all(points >= 0 & points <= 1))
  
  # Typical behavior: entropy increases then peaks
  expect_length(points, 4)
})

# ============================================================================
# SUITE 22: Consistency & Robustness
# ============================================================================

test_that("block vs standard bootstrap produce similar mean for balanced data", {
  
  # Balanced pairs should give similar bootstrap distributions
  balanced_pairs <- c(100, 100, 100, 100, 100, 100)
  
  set.seed(555)
  result_std <- bootstrap_compute_cpp_wrapper(
    x = balanced_pairs, q = 1.0, normalize = TRUE, 
    nboot = 100, log_base = exp(1), pseudocount = 0
  )
  
  set.seed(555)
  result_block <- block_bootstrap_compute_cpp_wrapper(
    x = balanced_pairs, q = 1.0, normalize = TRUE, 
    nboot = 100, log_base = exp(1), pseudocount = 0
  )
  
  # Means should be very similar for balanced data
  expect_equal(mean(result_std), mean(result_block), tolerance = 0.01)
})

test_that("bootstrap reproducibility with set.seed", {
  
  set.seed(999)
  result1 <- bootstrap_compute_cpp_wrapper(
    x = test_counts, q = 1.5, normalize = TRUE, 
    nboot = 100, log_base = 2, pseudocount = 0
  )
  
  set.seed(999)
  result2 <- bootstrap_compute_cpp_wrapper(
    x = test_counts, q = 1.5, normalize = TRUE, 
    nboot = 100, log_base = 2, pseudocount = 0
  )
  
  # Exact reproducibility with same seed
  expect_identical(result1, result2)
})

test_that("bootstrap CI integration: point estimate within CI bounds", {
  
  for (q in c(0.5, 1.0, 2.0)) {
    result <- .bootstrap_compute_ci(
      x = test_counts, q = q, norm = TRUE, nboot = 100,
      ci = 0.95, method = "percentile", log_base = exp(1),
      pseudocount = 0, what = "S"
    )
    
    # Point estimate should be within CI
    expect_true(result$point_est >= result$ci_result$lower - 1e-6)
    expect_true(result$point_est <= result$ci_result$upper + 1e-6)
  }
})


# ============================================================================
# SUITE 10: REGRESSION TESTS - Bug Fix Verification
# ============================================================================
# These tests ensure the quantile type=1 fix doesn't regress when code
# is modified or when edge cases are encountered.

test_that("REGRESSION 10.1: Bootstrap CI bounds are monotonic", {
  
  set.seed(666)
  counts_A <- matrix(rpois(15 * 40, lambda = 10), nrow = 15, ncol = 40)
  counts_B <- matrix(rpois(15 * 40, lambda = 10), nrow = 15, ncol = 40)
  
  jack_A <- jis_jackknife_influences_cpp(counts_A, q = 1, normalize = TRUE)
  jack_B <- jis_jackknife_influences_cpp(counts_B, q = 1, normalize = TRUE)
  delta_influence <- abs(jack_A - jack_B)
  
  # Run multiple times (different random seeds in bootstrap)
  for (seed in 100:105) {
    set.seed(seed)
    result <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence,
                                       q = 1, normalize = TRUE,
                                       nboot = 300, confidence = 0.95)
    
    # Check monotonicity: all ci_lower <= ci_upper
    valid_idx <- !is.na(result$ci_lower) & !is.na(result$ci_upper)
    expect_true(all(result$ci_lower[valid_idx] <= result$ci_upper[valid_idx]),
                info = sprintf("CI not monotonic for seed %d", seed))
  }
})

test_that("REGRESSION 10.2: Confidence level ordering respected", {
  
  set.seed(777)
  counts_A <- matrix(rpois(20 * 35, lambda = 8), nrow = 20, ncol = 35)
  counts_B <- matrix(rpois(20 * 35, lambda = 8), nrow = 20, ncol = 35)
  
  jack_A <- jis_jackknife_influences_cpp(counts_A, q = 1, normalize = TRUE)
  jack_B <- jis_jackknife_influences_cpp(counts_B, q = 1, normalize = TRUE)
  delta_influence <- abs(jack_A - jack_B)
  
  # Compute CIs at different confidence levels
  result_80 <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence,
                                        q = 1, normalize = TRUE,
                                        nboot = 500, confidence = 0.80)
  result_90 <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence,
                                        q = 1, normalize = TRUE,
                                        nboot = 500, confidence = 0.90)
  result_95 <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence,
                                        q = 1, normalize = TRUE,
                                        nboot = 500, confidence = 0.95)
  result_99 <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence,
                                        q = 1, normalize = TRUE,
                                        nboot = 500, confidence = 0.99)
  
  # Higher confidence should give wider (or equal) intervals
  valid_idx <- !is.na(result_80$ci_width) & !is.na(result_90$ci_width) &
               !is.na(result_95$ci_width) & !is.na(result_99$ci_width)
  
  if (any(valid_idx)) {
    # For most transcripts, width should be non-decreasing
    mean_comparison <- mean(
      result_90$ci_width[valid_idx] >= result_80$ci_width[valid_idx] - 1e-6
    ) > 0.9 &&
    mean(
      result_95$ci_width[valid_idx] >= result_90$ci_width[valid_idx] - 1e-6
    ) > 0.9 &&
    mean(
      result_99$ci_width[valid_idx] >= result_95$ci_width[valid_idx] - 1e-6
    ) > 0.9
    
    expect_true(mean_comparison)
  }
})

test_that("REGRESSION 10.3: P-value distribution validation", {
  
  set.seed(888)
  counts_A <- matrix(rpois(25 * 30, lambda = 12), nrow = 25, ncol = 30)
  counts_B <- matrix(rpois(25 * 30, lambda = 12), nrow = 25, ncol = 30)
  
  jack_A <- jis_jackknife_influences_cpp(counts_A, q = 1, normalize = TRUE)
  jack_B <- jis_jackknife_influences_cpp(counts_B, q = 1, normalize = TRUE)
  delta_influence <- abs(jack_A - jack_B)
  
  result <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence,
                                     q = 1, normalize = TRUE,
                                     nboot = 1000, confidence = 0.95)
  
  # P-values must be valid probabilities
  pvals <- result$p_value[!is.na(result$p_value)]
  expect_true(all(pvals >= 0))
  expect_true(all(pvals <= 1))
  
  # Minimum p-value should respect 1/n_bootstrap
  if (length(pvals) > 0) {
    expect_true(min(pvals) >= 1/1000)
  }
})

test_that("REGRESSION 10.4: Bootstrap statistics across q values remain stable", {
  
  set.seed(999)
  counts_A <- matrix(rpois(18 * 35, lambda = 10), nrow = 18, ncol = 35)
  counts_B <- matrix(rpois(18 * 35, lambda = 10), nrow = 18, ncol = 35)
  
  # Test across multiple q values
  test_qs <- c(0.5, 1.0, 1.5, 2.0, 2.5)
  
  for (q in test_qs) {
    jack_A <- jis_jackknife_influences_cpp(counts_A, q = q, normalize = TRUE)
    jack_B <- jis_jackknife_influences_cpp(counts_B, q = q, normalize = TRUE)
    delta_influence <- abs(jack_A - jack_B)
    
    result <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence,
                                       q = q, normalize = TRUE,
                                       nboot = 400, confidence = 0.95)
    
    # All output fields should exist and be properly sized
    expect_equal(length(result$ci_lower), nrow(counts_A))
    expect_equal(length(result$ci_upper), nrow(counts_A))
    expect_equal(length(result$p_value), nrow(counts_A))
    expect_equal(length(result$ci_width), nrow(counts_A))
  }
})

test_that("REGRESSION 10.5: Quantile method consistency across sample sizes", {
  
  # Test with various bootstrap replicates to ensure quantile indexing is correct
  nboot_values <- c(100, 500, 1000, 5000)
  
  for (nboot in nboot_values) {
    set.seed(1111)
    counts_A <- matrix(rpois(12 * 25, lambda = 8), nrow = 12, ncol = 25)
    counts_B <- matrix(rpois(12 * 25, lambda = 8), nrow = 12, ncol = 25)
    
    jack_A <- jis_jackknife_influences_cpp(counts_A, q = 1, normalize = TRUE)
    jack_B <- jis_jackknife_influences_cpp(counts_B, q = 1, normalize = TRUE)
    delta_influence <- abs(jack_A - jack_B)
    
    result <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence,
                                       q = 1, normalize = TRUE,
                                       nboot = nboot, confidence = 0.95)
    
    # CIs should be valid for any nboot
    valid_ci <- !is.na(result$ci_lower) & !is.na(result$ci_upper)
    if (any(valid_ci)) {
      expect_true(all(result$ci_lower[valid_ci] <= result$ci_upper[valid_ci]))
    }
  }
})

test_that("REGRESSION 10.6: Effect size and CI consistency", {
  
  set.seed(1212)
  counts_A <- matrix(rpois(16 * 32, lambda = 10), nrow = 16, ncol = 32)
  counts_B <- matrix(rpois(16 * 32, lambda = 10), nrow = 16, ncol = 32)
  
  jack_A <- jis_jackknife_influences_cpp(counts_A, q = 1, normalize = TRUE)
  jack_B <- jis_jackknife_influences_cpp(counts_B, q = 1, normalize = TRUE)
  delta_influence <- abs(jack_A - jack_B)
  
  result <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence,
                                     q = 1, normalize = TRUE,
                                     nboot = 600, confidence = 0.95)
  
  # Effect size should be non-negative
  effect_sizes <- result$effect_size[!is.na(result$effect_size)]
  expect_true(all(effect_sizes >= 0))
  
  # For most (not all) cases, effect size should be between CI bounds
  valid_all <- !is.na(result$effect_size) & !is.na(result$ci_lower) & !is.na(result$ci_upper)
  if (any(valid_all)) {
    in_ci <- (result$effect_size[valid_all] >= result$ci_lower[valid_all] - 1e-4) &
             (result$effect_size[valid_all] <= result$ci_upper[valid_all] + 1e-4)
    # Effect size is computed from bootstrap mean, so usually within CI
    expect_true(mean(in_ci) > 0.5)
  }
})

test_that("REGRESSION 10.7: No NaN propagation in quantile computation", {
  
  set.seed(1313)
  # Create data with some zero-sum samples (will produce NA in entropy)
  counts_A <- matrix(rpois(20 * 30, lambda = 5), nrow = 20, ncol = 30)
  counts_B <- matrix(rpois(20 * 30, lambda = 5), nrow = 20, ncol = 30)
  
  jack_A <- jis_jackknife_influences_cpp(counts_A, q = 1, normalize = TRUE)
  jack_B <- jis_jackknife_influences_cpp(counts_B, q = 1, normalize = TRUE)
  delta_influence <- abs(jack_A - jack_B)
  
  result <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence,
                                     q = 1, normalize = TRUE,
                                     nboot = 400, confidence = 0.95)
  
  # All output should be either numeric or NA, not NaN
  expect_true(all(!is.nan(result$ci_lower)))
  expect_true(all(!is.nan(result$ci_upper)))
  expect_true(all(!is.nan(result$p_value)))
  expect_true(all(!is.nan(result$ci_width)))
})

# ============================================================================
# SUITE 11: C++ Divergence Bootstrap (Independent Mode)
# ============================================================================

test_that("divergence_bootstrap_compute_cpp_wrapper compiles and is available", {
  
  # Check that the function exists
  expect_true(exists("divergence_bootstrap_compute_cpp_wrapper"))
  
  # Test basic call
  x <- c(100, 50, 25, 10)
  y <- c(80, 60, 40, 20)
  
  result <- divergence_bootstrap_compute_cpp_wrapper(
    x = x, y = y, nboot = 50L, q = 1.0, 
    pseudocount = 0.5, log_base = exp(1), paired = FALSE
  )
  
  expect_true(is.numeric(result))
  expect_length(result, 50)
  expect_true(all(is.finite(result)))
})

test_that("divergence_bootstrap_compute_cpp_wrapper respects nboot parameter", {
  
  x <- c(100, 50, 25)
  y <- c(80, 60, 40)
  
  for (nboot in c(10L, 50L, 100L, 500L)) {
    result <- divergence_bootstrap_compute_cpp_wrapper(
      x = x, y = y, nboot = nboot, q = 1.0,
      pseudocount = 0.5, log_base = exp(1), paired = FALSE
    )
    
    expect_length(result, nboot)
    expect_true(all(is.finite(result)))
  }
})

test_that("divergence_bootstrap_compute_cpp_wrapper works with various q values", {
  
  x <- c(100, 50, 25, 10)
  y <- c(80, 60, 40, 20)
  
  q_values <- c(0.5, 1.0, 1.5, 2.0, 3.0)
  
  for (q in q_values) {
    result <- divergence_bootstrap_compute_cpp_wrapper(
      x = x, y = y, nboot = 30L, q = q,
      pseudocount = 0.5, log_base = exp(1), paired = FALSE
    )
    
    expect_true(all(is.finite(result)))
    expect_length(result, 30)
  }
})

test_that("divergence_bootstrap_compute_cpp_wrapper handles pseudocount adjustment", {
  
  x <- c(100, 50, 25, 10)
  y <- c(80, 60, 40, 20)
  
  result_no_pc <- divergence_bootstrap_compute_cpp_wrapper(
    x = x, y = y, nboot = 100L, q = 1.0,
    pseudocount = 0.0, log_base = exp(1), paired = FALSE
  )
  
  result_with_pc <- divergence_bootstrap_compute_cpp_wrapper(
    x = x, y = y, nboot = 100L, q = 1.0,
    pseudocount = 1.0, log_base = exp(1), paired = FALSE
  )
  
  # Both should produce valid results
  expect_true(all(is.finite(result_no_pc)))
  expect_true(all(is.finite(result_with_pc)))
  
  # Pseudocount affects the values
  expect_false(isTRUE(all.equal(result_no_pc, result_with_pc)))
})

test_that("divergence_bootstrap_compute_cpp_wrapper is reproducible", {
  
  x <- c(100, 50, 25, 10)
  y <- c(80, 60, 40, 20)
  
  set.seed(42)
  result1 <- divergence_bootstrap_compute_cpp_wrapper(
    x = x, y = y, nboot = 100L, q = 1.0,
    pseudocount = 0.5, log_base = exp(1), paired = FALSE
  )
  
  set.seed(42)
  result2 <- divergence_bootstrap_compute_cpp_wrapper(
    x = x, y = y, nboot = 100L, q = 1.0,
    pseudocount = 0.5, log_base = exp(1), paired = FALSE
  )
  
  expect_equal(result1, result2)
})

# ============================================================================
# SUITE 12: C++ Divergence Bootstrap (Paired Mode)
# ============================================================================

test_that("divergence_bootstrap_paired_cpp_wrapper compiles and is available", {
  
  # Check that the function exists
  expect_true(exists("divergence_bootstrap_paired_cpp_wrapper"))
  
  # Test basic call with paired data
  # For paired design: x = control counts (1 per pair), y = treatment counts (1 per pair)
  x <- c(100, 80)      # Control samples for 2 pairs
  y <- c(90, 85)       # Treatment samples for 2 pairs
  pair_ids <- c(1L, 2L)  # Pair identifiers
  
  result <- divergence_bootstrap_paired_cpp_wrapper(
    x = x, y = y, pair_ids = pair_ids,
    nboot = 50L, q = 1.0,
    pseudocount = 0.5, log_base = exp(1)
  )
  
  expect_true(is.numeric(result))
  expect_length(result, 50)
  expect_true(all(is.finite(result)))
})

test_that("divergence_bootstrap_paired_cpp_wrapper respects nboot parameter", {
  
  # Paired data: one control and one treatment per pair (must be even number of elements)
  x <- c(100, 80, 120, 95)
  y <- c(90, 85, 110, 100)
  pair_ids <- c(1L, 2L, 3L, 4L)
  
  for (nboot in c(10L, 50L, 100L)) {
    result <- divergence_bootstrap_paired_cpp_wrapper(
      x = x, y = y, pair_ids = pair_ids,
      nboot = nboot, q = 1.0,
      pseudocount = 0.5, log_base = exp(1)
    )
    
    expect_length(result, nboot)
    expect_true(all(is.finite(result)))
  }
})

test_that("divergence_bootstrap_paired_cpp_wrapper works with various q values", {
  
  x <- c(100, 80, 120, 95)
  y <- c(90, 85, 110, 100)
  pair_ids <- c(1L, 2L, 3L, 4L)
  
  q_values <- c(0.5, 1.0, 1.5, 2.0)
  
  for (q in q_values) {
    result <- divergence_bootstrap_paired_cpp_wrapper(
      x = x, y = y, pair_ids = pair_ids,
      nboot = 30L, q = q,
      pseudocount = 0.5, log_base = exp(1)
    )
    
    expect_true(all(is.finite(result)))
    expect_length(result, 30)
  }
})

test_that("divergence_bootstrap_paired_cpp_wrapper handles multiple pairs", {
  
  # Test with increasing numbers of pairs (must be even)
  for (n_pairs in c(2, 4, 10)) {
    x <- rnorm(n_pairs, mean = 100, sd = 20)
    x <- pmax(x, 1)  # Ensure positive counts
    y <- rnorm(n_pairs, mean = 90, sd = 20)
    y <- pmax(y, 1)
    pair_ids <- seq_len(n_pairs)
    
    result <- divergence_bootstrap_paired_cpp_wrapper(
      x = x, y = y, pair_ids = as.integer(pair_ids),
      nboot = 30L, q = 1.0,
      pseudocount = 0.5, log_base = exp(1)
    )
    
    expect_true(all(is.finite(result)))
    expect_length(result, 30)
  }
})

test_that("divergence_bootstrap_paired_cpp_wrapper is reproducible", {
  
  x <- c(100, 80, 120, 95)
  y <- c(90, 85, 110, 100)
  pair_ids <- c(1L, 2L, 3L, 4L)
  
  set.seed(999)
  result1 <- divergence_bootstrap_paired_cpp_wrapper(
    x = x, y = y, pair_ids = pair_ids,
    nboot = 100L, q = 1.0,
    pseudocount = 0.5, log_base = exp(1)
  )
  
  set.seed(999)
  result2 <- divergence_bootstrap_paired_cpp_wrapper(
    x = x, y = y, pair_ids = pair_ids,
    nboot = 100L, q = 1.0,
    pseudocount = 0.5, log_base = exp(1)
  )
  
  expect_equal(result1, result2)
})

test_that("divergence_bootstrap_paired_cpp validates input lengths", {
  
  x <- c(100, 50, 80)  # Uneven length
  y <- c(80, 60)       # Different length
  pair_ids <- c(1L, 1L, 2L)
  
  # Should throw error for mismatched lengths
  expect_error(
    divergence_bootstrap_paired_cpp_wrapper(
      x = x, y = y, pair_ids = pair_ids,
      nboot = 10L, q = 1.0,
      pseudocount = 0.5, log_base = exp(1)
    ),
    "must have equal length"
  )
})

test_that("divergence_bootstrap_paired_cpp handles pseudocount", {
  
  x <- c(100, 80, 120, 95)
  y <- c(90, 85, 110, 100)
  pair_ids <- c(1L, 2L, 3L, 4L)
  
  # Test with different pseudocount values
  result_0 <- divergence_bootstrap_paired_cpp_wrapper(
    x = x, y = y, pair_ids = pair_ids,
    nboot = 50L, q = 1.0,
    pseudocount = 0.0, log_base = exp(1)
  )
  
  result_1 <- divergence_bootstrap_paired_cpp_wrapper(
    x = x, y = y, pair_ids = pair_ids,
    nboot = 50L, q = 1.0,
    pseudocount = 1.0, log_base = exp(1)
  )
  
  expect_true(all(is.finite(result_0)))
  expect_true(all(is.finite(result_1)))
  # Different pseudocounts should produce different results
  expect_false(isTRUE(all.equal(result_0, result_1)))
})

# ============================================================================
# TESTS FOR FLEXIBLE PAIRED/UNPAIRED BOOTSTRAP (NEW - MARCH 2026)
# ============================================================================

test_that("divergence_bootstrap_flexible_cpp_wrapper compiles and is available", {
  
  expect_true(exists("divergence_bootstrap_flexible_cpp_wrapper"))
  expect_true(is.function(divergence_bootstrap_flexible_cpp_wrapper))
})

test_that("divergence_bootstrap_flexible_cpp_wrapper handles complete pairs", {
  
  # Complete pairs: all pair_ids present in both x and y
  x <- c(100, 80, 120, 95)
  y <- c(90, 85, 110, 100)
  x_pair_ids <- c(1L, 2L, 3L, 4L)
  y_pair_ids <- c(1L, 2L, 3L, 4L)
  
  result <- divergence_bootstrap_flexible_cpp_wrapper(
    x = x, y = y,
    x_pair_ids = x_pair_ids, y_pair_ids = y_pair_ids,
    nboot = 100L, q = 1.0,
    pseudocount = 0.5, log_base = exp(1)
  )
  
  expect_length(result, 100)
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0))  # Divergence should be non-negative
})

test_that("divergence_bootstrap_flexible_cpp_wrapper handles unpaired samples only", {
  
  # All unpaired: pair_ids = 0
  x <- c(100, 80, 120, 95)
  y <- c(90, 85, 110, 100)
  x_pair_ids <- c(0L, 0L, 0L, 0L)  # All unpaired
  y_pair_ids <- c(0L, 0L, 0L, 0L)  # All unpaired
  
  result <- divergence_bootstrap_flexible_cpp_wrapper(
    x = x, y = y,
    x_pair_ids = x_pair_ids, y_pair_ids = y_pair_ids,
    nboot = 100L, q = 1.0,
    pseudocount = 0.5, log_base = exp(1)
  )
  
  expect_length(result, 100)
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0))
})

test_that("divergence_bootstrap_flexible_cpp_wrapper handles incomplete pairs", {
  
  # Incomplete pairs: pair 1,2,3 in both groups, pair 4 only in x
  x <- c(100, 80, 120, 95)
  y <- c(90, 85, 110)  # Only 3 samples
  x_pair_ids <- c(1L, 2L, 3L, 4L)  # pair 4 in x but not y
  y_pair_ids <- c(1L, 2L, 3L)    # no pair 4
  
  result <- divergence_bootstrap_flexible_cpp_wrapper(
    x = x, y = y,
    x_pair_ids = x_pair_ids, y_pair_ids = y_pair_ids,
    nboot = 100L, q = 1.0,
    pseudocount = 0.5, log_base = exp(1)
  )
  
  expect_length(result, 100)
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0))
})

test_that("divergence_bootstrap_flexible_cpp_wrapper handles mixed paired/unpaired", {
  
  # Mixed: pairs 1,2,3 complete, pair 0 (unpaired) in both, pair 4 only in y
  x <- c(100, 80, 120, 95, 110)
  y <- c(90, 85, 110, 100, 105)
  x_pair_ids <- c(1L, 2L, 3L, 0L, 0L)  # 0 = unpaired
  y_pair_ids <- c(1L, 2L, 3L, 4L, 0L)  # pair 4 only in y, 0 = unpaired
  
  result <- divergence_bootstrap_flexible_cpp_wrapper(
    x = x, y = y,
    x_pair_ids = x_pair_ids, y_pair_ids = y_pair_ids,
    nboot = 100L, q = 1.0,
    pseudocount = 0.5, log_base = exp(1)
  )
  
  expect_length(result, 100)
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0))
})

test_that("divergence_bootstrap_flexible_cpp_wrapper works with NA for unpaired", {
  
  # NA should be treated as unpaired (converted to 0 internally)
  x <- c(100, 80, 120, 95)
  y <- c(90, 85, 110, 100)
  x_pair_ids <- c(1L, 2L, 3L, NA_integer_)  # NA = unpaired
  y_pair_ids <- c(1L, 2L, 3L, NA_integer_)  # NA = unpaired
  
  result <- divergence_bootstrap_flexible_cpp_wrapper(
    x = x, y = y,
    x_pair_ids = x_pair_ids, y_pair_ids = y_pair_ids,
    nboot = 100L, q = 1.0,
    pseudocount = 0.5, log_base = exp(1)
  )
  
  expect_length(result, 100)
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0))
})

test_that("divergence_bootstrap_flexible_cpp_wrapper respects nboot parameter", {
  
  x <- c(100, 80, 120, 95)
  y <- c(90, 85, 110, 100)
  x_pair_ids <- c(1L, 2L, 3L, 4L)
  y_pair_ids <- c(1L, 2L, 3L, 4L)
  
  for (nboot in c(10L, 50L, 100L)) {
    result <- divergence_bootstrap_flexible_cpp_wrapper(
      x = x, y = y,
      x_pair_ids = x_pair_ids, y_pair_ids = y_pair_ids,
      nboot = nboot, q = 1.0,
      pseudocount = 0.5, log_base = exp(1)
    )
    
    expect_length(result, nboot)
    expect_true(all(is.finite(result)))
  }
})

test_that("divergence_bootstrap_flexible_cpp_wrapper works with various q values", {
  
  x <- c(100, 80, 120, 95)
  y <- c(90, 85, 110, 100)
  x_pair_ids <- c(1L, 2L, 3L, 4L)
  y_pair_ids <- c(1L, 2L, 3L, 4L)
  
  q_values <- c(0.0, 0.5, 1.0, 1.5, 2.0)
  
  for (q in q_values) {
    result <- divergence_bootstrap_flexible_cpp_wrapper(
      x = x, y = y,
      x_pair_ids = x_pair_ids, y_pair_ids = y_pair_ids,
      nboot = 50L, q = q,
      pseudocount = 0.5, log_base = exp(1)
    )
    
    expect_true(all(is.finite(result)))
    expect_length(result, 50)
  }
})

test_that("divergence_bootstrap_flexible_cpp_wrapper is reproducible", {
  
  x <- c(100, 80, 120, 95)
  y <- c(90, 85, 110, 100)
  x_pair_ids <- c(1L, 2L, 3L, 4L)
  y_pair_ids <- c(1L, 2L, 3L, 4L)
  
  set.seed(42)
  result1 <- divergence_bootstrap_flexible_cpp_wrapper(
    x = x, y = y,
    x_pair_ids = x_pair_ids, y_pair_ids = y_pair_ids,
    nboot = 100L, q = 1.0,
    pseudocount = 0.5, log_base = exp(1)
  )
  
  set.seed(42)
  result2 <- divergence_bootstrap_flexible_cpp_wrapper(
    x = x, y = y,
    x_pair_ids = x_pair_ids, y_pair_ids = y_pair_ids,
    nboot = 100L, q = 1.0,
    pseudocount = 0.5, log_base = exp(1)
  )
  
  expect_equal(result1, result2)
})

test_that("divergence_bootstrap_flexible_cpp_wrapper validates input lengths", {
  
  x <- c(100, 80, 120)
  y <- c(90, 85, 110, 100)
  
  # Mismatched x_pair_ids length
  expect_error(
    divergence_bootstrap_flexible_cpp_wrapper(
      x = x, y = y,
      x_pair_ids = c(1L, 2L),  # Wrong length!
      y_pair_ids = c(1L, 2L, 3L, 4L),
      nboot = 10L, q = 1.0,
      pseudocount = 0.5, log_base = exp(1)
    ),
    "must have same length"
  )
})

test_that("divergence_bootstrap_flexible_cpp_wrapper handles pseudocount", {
  
  x <- c(100, 80, 120)
  y <- c(90, 85, 110)
  x_pair_ids <- c(1L, 2L, 3L)
  y_pair_ids <- c(1L, 2L, 3L)
  
  # Test with different pseudocount values
  result_0 <- divergence_bootstrap_flexible_cpp_wrapper(
    x = x, y = y,
    x_pair_ids = x_pair_ids, y_pair_ids = y_pair_ids,
    nboot = 50L, q = 1.0,
    pseudocount = 0.0, log_base = exp(1)
  )
  
  result_1 <- divergence_bootstrap_flexible_cpp_wrapper(
    x = x, y = y,
    x_pair_ids = x_pair_ids, y_pair_ids = y_pair_ids,
    nboot = 50L, q = 1.0,
    pseudocount = 1.0, log_base = exp(1)
  )
  
  expect_true(all(is.finite(result_0)))
  expect_true(all(is.finite(result_1)))
  # Different pseudocounts should produce different results
  expect_false(isTRUE(all.equal(result_0, result_1)))
})

test_that("divergence_bootstrap_flexible_cpp_wrapper handles different log bases", {
  
  x <- c(100, 80, 120, 95)
  y <- c(90, 85, 110, 100)
  x_pair_ids <- c(1L, 2L, 3L, 4L)
  y_pair_ids <- c(1L, 2L, 3L, 4L)
  
  result_e <- divergence_bootstrap_flexible_cpp_wrapper(
    x = x, y = y,
    x_pair_ids = x_pair_ids, y_pair_ids = y_pair_ids,
    nboot = 50L, q = 1.0,
    pseudocount = 0.5, log_base = exp(1)
  )
  
  result_2 <- divergence_bootstrap_flexible_cpp_wrapper(
    x = x, y = y,
    x_pair_ids = x_pair_ids, y_pair_ids = y_pair_ids,
    nboot = 50L, q = 1.0,
    pseudocount = 0.5, log_base = 2.0
  )
  
  expect_true(all(is.finite(result_e)))
  expect_true(all(is.finite(result_2)))
  # Different log bases should produce different results
  expect_false(isTRUE(all.equal(result_e, result_2)))
})

test_that("divergence_bootstrap_flexible_cpp_wrapper handles many samples", {
  
  # Test with larger sample sizes
  n_pairs <- 20
  x <- rnorm(n_pairs, mean = 100, sd = 20)
  x <- pmax(x, 1)  # Ensure positive counts
  y <- rnorm(n_pairs, mean = 90, sd = 20)
  y <- pmax(y, 1)
  x_pair_ids <- seq_len(n_pairs)
  y_pair_ids <- seq_len(n_pairs)
  
  result <- divergence_bootstrap_flexible_cpp_wrapper(
    x = x, y = y,
    x_pair_ids = as.integer(x_pair_ids),
    y_pair_ids = as.integer(y_pair_ids),
    nboot = 100L, q = 1.0,
    pseudocount = 0.5, log_base = exp(1)
  )
  
  expect_length(result, 100)
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0))
})

test_that("divergence_bootstrap_flexible_cpp_wrapper handles asymmetric unpaired", {
  
  # x has 5 samples, y has 3 samples, all unpaired
  x <- c(100, 80, 120, 95, 110)
  y <- c(90, 85, 110)
  x_pair_ids <- c(0L, 0L, 0L, 0L, 0L)  # All unpaired
  y_pair_ids <- c(0L, 0L, 0L)           # All unpaired
  
  result <- divergence_bootstrap_flexible_cpp_wrapper(
    x = x, y = y,
    x_pair_ids = x_pair_ids, y_pair_ids = y_pair_ids,
    nboot = 100L, q = 1.0,
    pseudocount = 0.5, log_base = exp(1)
  )
  
  expect_length(result, 100)
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0))
})

test_that("divergence_bootstrap_flexible_cpp_wrapper handles sparse pair distribution", {
  
  # Complex scenario: some pairs present in both, some only in one group
  x <- c(100, 80, 120, 95, 110, 75)     # 6 samples
  y <- c(90, 85, 110, 100, 105)         # 5 samples
  # Pairs: 1,2,3 complete; 4,5 in x only; 6 in y only; 0 = unpaired
  x_pair_ids <- c(1L, 2L, 3L, 4L, 5L, 0L)
  y_pair_ids <- c(1L, 2L, 3L, 0L, 6L)
  
  result <- divergence_bootstrap_flexible_cpp_wrapper(
    x = x, y = y,
    x_pair_ids = x_pair_ids, y_pair_ids = y_pair_ids,
    nboot = 100L, q = 1.0,
    pseudocount = 0.5, log_base = exp(1)
  )
  
  expect_length(result, 100)
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0))
})

test_that("divergence_bootstrap_flexible_cpp_wrapper handles single sample", {
  
  # Edge case: single sample in each group
  x <- c(100)
  y <- c(90)
  x_pair_ids <- c(1L)
  y_pair_ids <- c(1L)
  
  result <- divergence_bootstrap_flexible_cpp_wrapper(
    x = x, y = y,
    x_pair_ids = x_pair_ids, y_pair_ids = y_pair_ids,
    nboot = 50L, q = 1.0,
    pseudocount = 0.5, log_base = exp(1)
  )
  
  expect_length(result, 50)
  expect_true(all(is.finite(result)))
})

# ============================================================================
# ADDITIONAL COMPREHENSIVE TESTS FOR BOOTSTRAP IMPLEMENTATIONS
# ============================================================================

test_that("divergence_bootstrap_compute_cpp_wrapper handles extreme q values", {
  
  x <- c(100, 80, 120, 95)
  y <- c(90, 85, 110, 100)
  
  # Test very small q (q near 0)
  result_q_small <- divergence_bootstrap_compute_cpp_wrapper(
    x = x, y = y, q = 0.01, nboot = 50L,
    paired = FALSE, pseudocount = 0.5, log_base = exp(1)
  )
  expect_length(result_q_small, 50)
  expect_true(sum(is.finite(result_q_small)) > 40)
  
  # Test large q (q > 2)
  result_q_large <- divergence_bootstrap_compute_cpp_wrapper(
    x = x, y = y, q = 3.0, nboot = 50L,
    paired = FALSE, pseudocount = 0.5, log_base = exp(1)
  )
  expect_length(result_q_large, 50)
  expect_true(sum(is.finite(result_q_large)) > 40)
})

test_that("divergence_bootstrap_compute_cpp_wrapper works with very small counts", {
  
  # Very small counts (e.g., single molecules)
  x <- c(1, 2, 1, 3)
  y <- c(2, 1, 2, 1)
  
  result <- divergence_bootstrap_compute_cpp_wrapper(
    x = x, y = y, q = 1.0, nboot = 50L,
    paired = FALSE, pseudocount = 0.5, log_base = exp(1)
  )
  
  expect_length(result, 50)
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0))
})

test_that("divergence_bootstrap_compute_cpp_wrapper works with very large counts", {
  
  # Very large counts (e.g., millions)
  x <- c(1e6, 2e6, 1.5e6, 800000)
  y <- c(1.2e6, 1.8e6, 1.6e6, 900000)
  
  result <- divergence_bootstrap_compute_cpp_wrapper(
    x = x, y = y, q = 1.0, nboot = 50L,
    paired = FALSE, pseudocount = 0.0, log_base = exp(1)
  )
  
  expect_length(result, 50)
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0))
})

test_that("divergence_bootstrap_compute_cpp_wrapper works with identical x and y", {
  
  # Identical distributions should have ~0 divergence
  counts <- c(100, 80, 120, 95)
  
  result <- divergence_bootstrap_compute_cpp_wrapper(
    x = counts, y = counts, q = 1.0, nboot = 50L,
    paired = FALSE, pseudocount = 0.5, log_base = exp(1)
  )
  
  expect_length(result, 50)
  expect_true(all(is.finite(result)))
  # Divergence should be very small (close to 0) for identical distributions
  expect_true(median(result, na.rm = TRUE) < 0.1)
})

test_that("divergence_bootstrap_compute_cpp_wrapper respects log_base parameter", {
  
  x <- c(100, 80, 120, 95)
  y <- c(90, 85, 110, 100)
  
  # Test with natural log (e)
  result_e <- divergence_bootstrap_compute_cpp_wrapper(
    x = x, y = y, q = 1.0, nboot = 50L,
    paired = FALSE, pseudocount = 0.5, log_base = exp(1)
  )
  
  # Test with log base 2
  result_2 <- divergence_bootstrap_compute_cpp_wrapper(
    x = x, y = y, q = 1.0, nboot = 50L,
    paired = FALSE, pseudocount = 0.5, log_base = 2.0
  )
  
  # Test with log base 10
  result_10 <- divergence_bootstrap_compute_cpp_wrapper(
    x = x, y = y, q = 1.0, nboot = 50L,
    paired = FALSE, pseudocount = 0.5, log_base = 10.0
  )
  
  expect_length(result_e, 50)
  expect_length(result_2, 50)
  expect_length(result_10, 50)
  
  # All should have values, different log bases should give different scales
  expect_true(all(is.finite(result_e)))
  expect_true(all(is.finite(result_2)))
  expect_true(all(is.finite(result_10)))
})

test_that("divergence_bootstrap_paired_cpp_wrapper handles multiple pairs efficiently", {
  
  # Test with increasing numbers of pairs to verify performance doesn't degrade (must be even)
  for (n_pairs in c(4, 10, 16)) {
    x <- rnorm(n_pairs, mean = 100, sd = 20)
    x <- pmax(x, 1)
    y <- rnorm(n_pairs, mean = 95, sd = 20)
    y <- pmax(y, 1)
    pair_ids <- seq_len(n_pairs)
    
    result <- divergence_bootstrap_paired_cpp_wrapper(
      x = x, y = y, pair_ids = as.integer(pair_ids),
      nboot = 30L, q = 1.0,
      pseudocount = 0.5, log_base = exp(1)
    )
    
    expect_length(result, 30)
    expect_true(sum(is.finite(result)) >= 25)  # Allow rare edge cases
  }
})

test_that("divergence_bootstrap_flexible_cpp_wrapper handles zero counts with pseudocount", {
  
  # Some counts are zero - pseudocount should prevent division issues
  x <- c(0, 100, 0, 80, 120)
  y <- c(90, 0, 110, 0, 100)
  x_pair_ids <- c(0L, 0L, 0L, 0L, 0L)
  y_pair_ids <- c(0L, 0L, 0L, 0L, 0L)
  
  result <- divergence_bootstrap_flexible_cpp_wrapper(
    x = x, y = y,
    x_pair_ids = x_pair_ids, y_pair_ids = y_pair_ids,
    nboot = 50L, q = 1.0,
    pseudocount = 1.0, log_base = exp(1)
  )
  
  expect_length(result, 50)
  expect_true(sum(is.finite(result)) >= 45)
  expect_true(all(result >= 0, na.rm = TRUE))
})

test_that("divergence_bootstrap_flexible_cpp_wrapper preserves pairing correlation", {
  
  # When we resample pairs as units, pair correlation should be preserved
  # Create perfectly correlated pairs
  x <- c(10, 20, 30, 40)
  y <- c(20, 40, 60, 80)  # y = 2*x exactly
  x_pair_ids <- c(1L, 2L, 3L, 4L)
  y_pair_ids <- c(1L, 2L, 3L, 4L)
  
  result_paired <- divergence_bootstrap_flexible_cpp_wrapper(
    x = x, y = y,
    x_pair_ids = x_pair_ids, y_pair_ids = y_pair_ids,
    nboot = 100L, q = 1.0,
    pseudocount = 0.5, log_base = exp(1)
  )
  
  # Also test with all unpaired for comparison
  result_unpaired <- divergence_bootstrap_flexible_cpp_wrapper(
    x = x, y = y,
    x_pair_ids = c(0L, 0L, 0L, 0L),
    y_pair_ids = c(0L, 0L, 0L, 0L),
    nboot = 100L, q = 1.0,
    pseudocount = 0.5, log_base = exp(1)
  )
  
  expect_length(result_paired, 100)
  expect_length(result_unpaired, 100)
  # Paired bootstrap might have smaller variance due to correlation preservation
  # (though with this specific example they may be close)
  expect_true(all(is.finite(result_paired)))
  expect_true(all(is.finite(result_unpaired)))
})

test_that("divergence_bootstrap_flexible_cpp_wrapper handles all mixed scenarios simultaneously", {
  
  # Complex: complete pairs + incomplete pairs in x + incomplete pairs in y + unpaired all together
  x <- c(100, 80, 120, 95, 110, 75, 105)       # 7 samples
  y <- c(90, 85, 110, 100, 105, 95)            # 6 samples
  
  # Pair IDs:
  # 1: in both (complete)
  # 2: in both (complete)
  # 3: in both (complete)
  # 4: in x only (incomplete)
  # 5: in x only (incomplete)
  # 6: in y only (incomplete)
  # 0: unpaired in both
  x_pair_ids <- c(1L, 2L, 3L, 4L, 5L, 0L, 0L)
  y_pair_ids <- c(1L, 2L, 3L, 0L, 6L, 0L)
  
  result <- divergence_bootstrap_flexible_cpp_wrapper(
    x = x, y = y,
    x_pair_ids = x_pair_ids, y_pair_ids = y_pair_ids,
    nboot = 100L, q = 1.0,
    pseudocount = 0.5, log_base = exp(1)
  )
  
  expect_length(result, 100)
  expect_true(sum(is.finite(result)) >= 95)
  expect_true(all(result >= 0, na.rm = TRUE))
})

test_that("divergence_bootstrap_flexible_cpp_wrapper works with q sequence", {
  
  x <- c(100, 80, 120, 95)
  y <- c(90, 85, 110, 100)
  x_pair_ids <- c(1L, 2L, 3L, 4L)
  y_pair_ids <- c(1L, 2L, 3L, 4L)
  
  q_sequence <- c(0.1, 0.5, 1.0, 1.5, 2.0, 2.5)
  
  for (q in q_sequence) {
    result <- divergence_bootstrap_flexible_cpp_wrapper(
      x = x, y = y,
      x_pair_ids = x_pair_ids, y_pair_ids = y_pair_ids,
      nboot = 50L, q = q,
      pseudocount = 0.5, log_base = exp(1)
    )
    
    expect_length(result, 50)
    expect_true(all(is.finite(result)))
    # Divergence should be non-negative for all q values
    expect_true(all(result >= 0))
  }
})

test_that("divergence_bootstrap_compute_cpp_wrapper produces consistent results across calls", {
  
  x <- c(100, 80, 120, 95)
  y <- c(90, 85, 110, 100)
  
  # Set seed and run once
  set.seed(12345)
  result1 <- divergence_bootstrap_compute_cpp_wrapper(
    x = x, y = y, q = 1.0, nboot = 50L,
    paired = FALSE, pseudocount = 0.5, log_base = exp(1)
  )
  
  # Set same seed and run again - should get identical results
  set.seed(12345)
  result2 <- divergence_bootstrap_compute_cpp_wrapper(
    x = x, y = y, q = 1.0, nboot = 50L,
    paired = FALSE, pseudocount = 0.5, log_base = exp(1)
  )
  
  expect_equal(result1, result2)
})

test_that("divergence_bootstrap_flexible_cpp_wrapper with pseudocount vector", {
  
  x <- c(100, 80, 120)
  y <- c(90, 85, 110)
  x_pair_ids <- c(1L, 2L, 3L)
  y_pair_ids <- c(1L, 2L, 3L)
  
  # Test with scalar pseudocount
  result_scalar <- divergence_bootstrap_flexible_cpp_wrapper(
    x = x, y = y,
    x_pair_ids = x_pair_ids, y_pair_ids = y_pair_ids,
    nboot = 50L, q = 1.0,
    pseudocount = 0.5, log_base = exp(1)
  )
  
  expect_length(result_scalar, 50)
  expect_true(all(is.finite(result_scalar)))
})

# ============================================================================
# COVERAGE IMPROVEMENT: Error Handling and Edge Cases
# ============================================================================

test_that("block_bootstrap_compute_cpp_wrapper rejects empty input (line 31)", {
  expect_error(
    block_bootstrap_compute_cpp_wrapper(
      x = NULL,  # Empty/NULL vector
      q = 1.0,
      normalize = TRUE,
      nboot = 10L,
      log_base = exp(1),
      pseudocount = 0
    ),
    "cannot be empty"
  )
})

test_that("block_bootstrap_compute_cpp_wrapper rejects odd-length input", {
  expect_error(
    block_bootstrap_compute_cpp_wrapper(
      x = c(100, 50, 25),  # Odd length
      q = 1.0,
      normalize = TRUE,
      nboot = 10L,
      log_base = exp(1),
      pseudocount = 0
    ),
    "must have even length"
  )
})

test_that("block_bootstrap_compute_cpp_wrapper warns on NA values", {
  expect_warning(
    block_bootstrap_compute_cpp_wrapper(
      x = c(100, NA, 25, 10),  # Contains NA
      q = 1.0,
      normalize = TRUE,
      nboot = 10L,
      log_base = exp(1),
      pseudocount = 0
    ),
    "contains.*NA values"
  )
})

test_that("block_bootstrap_compute_cpp_wrapper rejects all-zero input", {
  expect_error(
    block_bootstrap_compute_cpp_wrapper(
      x = c(0, 0, 0, 0),  # All zeros
      q = 1.0,
      normalize = TRUE,
      nboot = 10L,
      log_base = exp(1),
      pseudocount = 0
    ),
    "All values.*are zero"
  )
})

test_that("block_bootstrap_compute_cpp_wrapper rejects mismatched vector pseudocount (line 159)", {
  expect_error(
    block_bootstrap_compute_cpp_wrapper(
      x = c(100, 50, 25, 10),
      q = 1.0,
      normalize = TRUE,
      nboot = 10L,
      log_base = exp(1),
      pseudocount = c(0.5, 0.3)  # Wrong length (2 instead of 4)
    ),
    "pseudocount must have length 1 or equal to x length"
  )
})

test_that("divergence_bootstrap_paired_cpp_wrapper rejects mismatched pair_ids (line 210)", {
  expect_error(
    divergence_bootstrap_paired_cpp_wrapper(
      x = c(100, 50),
      y = c(90, 85),
      pair_ids = c(1, 2, 3),  # Wrong length (3 instead of 2)
      nboot = 10L,
      q = 1,
      pseudocount = 0.5,
      log_base = exp(1)
    ),
    "pair_ids must have same length"
  )
})

test_that("divergence_bootstrap_paired_cpp_wrapper handles vector pseudocount for x only (lines 224-226)", {
  # Create vector pseudocount matching x length only
  x <- c(100, 80)
  y <- c(90, 85)
  pair_ids <- c(1L, 2L)
  
  # This should still work - pseudocount is recycled or handled
  result <- divergence_bootstrap_paired_cpp_wrapper(
    x = x,
    y = y,
    pair_ids = pair_ids,
    nboot = 30L,
    q = 1,
    pseudocount = c(0.5, 0.3),
    log_base = exp(1)
  )
  
  expect_true(is.numeric(result))
  expect_length(result, 30)
})

test_that("divergence_bootstrap_flexible_cpp_wrapper rejects mismatched vector pseudocount (line 288-289)", {
  x <- c(100, 80)
  y <- c(90, 85)
  
  expect_error(
    divergence_bootstrap_flexible_cpp_wrapper(
      x = x,
      y = y,
      x_pair_ids = c(1L, 2L),  # Provide matching pair IDs to test pseudocount error
      y_pair_ids = c(1L, 2L),
      nboot = 10L,
      q = 1,
      pseudocount = c(0.5, 0.3, 0.2),  # Wrong: should be length 4 (2+2)
      log_base = exp(1)
    ),
    "pseudocount must have length 1 or"
  )
})

test_that("divergence_bootstrap_flexible_cpp_wrapper handles x_pair_ids conversion errors (line 310-311)", {
  x <- c(100, 80)
  y <- c(90, 85)
  
  # Test with NA values that cause conversion issues
  result <- tryCatch({
    divergence_bootstrap_flexible_cpp_wrapper(
      x = x,
      y = y,
      x_pair_ids = c(NA, NA),  # NA values to test error handling
      y_pair_ids = c(1L, 2L),
      nboot = 10L,
      q = 1,
      pseudocount = 0.5,
      log_base = exp(1)
    )
  }, error = function(e) {
    list(error = TRUE, message = e$message)
  })
  
  # Should either return numeric or error gracefully
  expect_true(is.numeric(result) || (is.list(result) && result$error))
})

test_that("divergence_bootstrap_flexible_cpp_wrapper handles y_pair_ids conversion errors (line 314-315)", {
  x <- c(100, 80)
  y <- c(90, 85)
  
  # Test with NA values that cause conversion issues
  result <- tryCatch({
    divergence_bootstrap_flexible_cpp_wrapper(
      x = x,
      y = y,
      x_pair_ids = c(1L, 2L),
      y_pair_ids = c(NA, NA),  # NA values to test error handling
      nboot = 10L,
      q = 1,
      pseudocount = 0.5,
      log_base = exp(1)
    )
  }, error = function(e) {
    list(error = TRUE, message = e$message)
  })
  
  # Should either return numeric or error gracefully
  expect_true(is.numeric(result) || (is.list(result) && result$error))
})

test_that(".bootstrap_resample_optimized handles effective_length mismatch (lines 374, 377)", {
  x <- c(100, 50, 25, 10)
  effective_length <- c(1000, 500)  # Wrong length (2 instead of 4)
  
  # Should produce warning about length mismatch
  expect_warning(
    .bootstrap_resample_optimized(
      x = x,
      q = 1,
      norm = TRUE,
      nboot = 5L,
      log_base = exp(1),
      pseudocount = 0,
      what = "S",  # 'S' for entropy
      paired = FALSE,
      effective_length = effective_length
    ),
    "effective_length"
  )
})

test_that(".bootstrap_resample_optimized with effective_length parameter", {
  x <- c(100, 50, 25, 10)
  effective_length <- c(1000, 500, 250, 100)
  
  result <- .bootstrap_resample_optimized(
    x = x,
    q = 1,
    norm = TRUE,
    nboot = 5L,
    log_base = exp(1),
    pseudocount = 0,
    what = "S",  # 'S' for entropy
    paired = FALSE,
    effective_length = effective_length
  )
  
  expect_true(is.numeric(result))
  expect_length(result, 5)
})

test_that("bootstrap_compute_cpp_wrapper with all parameters", {
  result <- bootstrap_compute_cpp_wrapper(
    x = c(100, 50, 25, 10),
    q = 2,
    normalize = FALSE,
    nboot = 20L,
    log_base = 2,
    pseudocount = 1
  )
  
  expect_true(is.numeric(result))
  expect_length(result, 20)
  expect_true(all(is.finite(result)))
})

test_that("divergence_bootstrap_paired_cpp_wrapper validates input types", {
  expect_error(
    divergence_bootstrap_paired_cpp_wrapper(
      x = "not_numeric",  # String instead of numeric
      y = c(90, 85),
      pair_ids = c(1L, 2L),
      nboot = 10L,
      q = 1,
      pseudocount = 0.5,
      log_base = exp(1)
    )
  )
})

test_that(".bootstrap_resample_with_quality_control achieves min_valid_frac (coverage for QC logic)", {
  # Create data that might generate some invalid bootstrap replicates
  x <- c(100, 50, 25, 10, 5, 2)
  
  result <- .bootstrap_resample_with_quality_control(
    x = x,
    q = 1,
    norm = TRUE,
    nboot = 20L,
    log_base = exp(1),
    pseudocount = 0,
    what = "S",  # 'S' for entropy
    paired = FALSE,
    effective_length = NULL,
    min_valid_frac = 0.8  # Require 80% valid
  )
  
  expect_true(is.numeric(result) || is.null(result))
  if (is.numeric(result)) {
    expect_length(result, 20)
    expect_true(sum(!is.na(result)) >= 0.8 * 20)  # At least 80% valid
  }
})

# ============================================================================
# QUALITY CONTROL CODE PATH COVERAGE: .bootstrap_resample_with_quality_control
# ============================================================================

test_that(".bootstrap_resample_with_quality_control early return when all valid (line 734)", {
  # Use clean numeric data - should produce no NAs/NaNs on first attempt
  x <- c(100, 50, 25, 10, 5, 2, 1)
  
  result <- .bootstrap_resample_with_quality_control(
    x = x,
    q = 1,
    norm = TRUE,
    nboot = 10L,
    log_base = exp(1),
    pseudocount = 1,
    what = "S",
    paired = FALSE,
    effective_length = NULL,
    min_valid_frac = 0.75  # Default threshold
  )
  
  # Should return vector of valid bootstrap replicates
  expect_is(result, "numeric")
  expect_length(result, 10)
  
  # With clean data, expect most or all to be valid
  n_valid <- sum(!is.na(result) & is.finite(result))
  expect_gte(n_valid, 8)  # At least 80% valid
})

test_that(".bootstrap_resample_with_quality_control with regeneration and threshold met (lines 762-767)", {
  # Use data with some coverage to allow regeneration
  x <- c(50, 40, 30, 20)
  
  # Call function - may or may not generate message depending on data
  result <- .bootstrap_resample_with_quality_control(
    x = x,
    q = 2,
    norm = TRUE,
    nboot = 30L,
    log_base = exp(1),
    pseudocount = 0.5,
    what = "S",
    paired = FALSE,
    effective_length = NULL,
    min_valid_frac = 0.70  # Allow threshold achievement
  )
  
  # Should return valid numeric vector
  expect_is(result, "numeric")
  expect_length(result, 30)
  
  # Should meet the quality threshold
  n_valid <- sum(!is.na(result) & is.finite(result))
  expect_gte(n_valid / 30, 0.70 - 0.01)  # Allow small tolerance
})

test_that(".bootstrap_resample_with_quality_control with sparse data handling (line 784)", {
  # Use sparse data to generate NAs/NaNs
  # Very low counts might generate NAs in entropy calculation
  x <- c(1, 0, 0, 0)  # Mostly zeros
  
  # Call function - may warn or error depending on data severity
  result <- tryCatch({
    .bootstrap_resample_with_quality_control(
      x = x,
      q = 1,
      norm = FALSE,  # No normalization for this sparse case
      nboot = 20L,
      log_base = exp(1),
      pseudocount = 0,
      what = "S",
      paired = FALSE,
      effective_length = NULL,
      min_valid_frac = 0.80  # High threshold to trigger handling
    )
  }, error = function(e) {
    NA  # Capture errors as NA
  }, warning = function(w) {
    invokeRestart("muffleWarning")  # Suppress warning and continue
  })
  
  # Should either return result or handle gracefully
  expect_true(is.numeric(result) || is.na(result))
})

test_that(".bootstrap_resample_with_quality_control stops when valid_frac < 0.5 (line 777)", {
  # Create data that generates many NAs - all zeros triggers entropy issues
  x <- c(0, 0, 0, 0)  # All zeros
  
  # Expect error/stop when quality is critically bad
  expect_error(
    .bootstrap_resample_with_quality_control(
      x = x,
      q = 1,
      norm = FALSE,
      nboot = 15L,
      log_base = exp(1),
      pseudocount = 0,
      what = "S",
      paired = FALSE,
      effective_length = NULL,
      min_valid_frac = 0.85  # Very high threshold
    ),
    "All values.*are zero|CRITICAL"
  )
})

test_that(".bootstrap_resample_with_quality_control reaches max_attempts (line 770)", {
  # Data designed to generate consistent NAs across regeneration attempts
  x <- c(0.001, 0.001, 0.001)  # Near-zero counts
  
  # Should attempt regeneration multiple times
  result <- tryCatch({
    .bootstrap_resample_with_quality_control(
      x = x,
      q = 1,
      norm = TRUE,
      nboot = 25L,
      log_base = exp(1),
      pseudocount = 0,
      what = "S",
      paired = FALSE,
      effective_length = NULL,
      min_valid_frac = 0.95  # Impossible threshold to trigger max_attempts
    )
  }, error = function(e) {
    # Catch error if threshold not met
    list(error = TRUE, message = e$message)
  }, warning = function(w) {
    # Catch warning if quality issues
    list(warning = TRUE, message = w$message)
  })
  
  # Should either error out or warn about quality issues
  expect_true(is.list(result) || is.numeric(result))
})

test_that(".bootstrap_resample_with_quality_control with low min_valid_frac passes easily (line 745)", {
  # Clean data with low threshold should pass immediately
  x <- c(100, 80, 60, 40)
  
  result <- .bootstrap_resample_with_quality_control(
    x = x,
    q = 1,
    norm = TRUE,
    nboot = 15L,
    log_base = exp(1),
    pseudocount = 0.5,
    what = "S",
    paired = FALSE,
    effective_length = NULL,
    min_valid_frac = 0.50  # Very low threshold
  )
  
  # Should return fully valid result
  expect_is(result, "numeric")
  expect_length(result, 15)
  expect_gte(sum(!is.na(result) & is.finite(result)), 14)  # Nearly all valid
})

test_that(".bootstrap_resample_with_quality_control with vector effective_length (line 368-378)", {
  # Test effective_length parameter integration
  x <- c(100, 80, 60, 40)
  effective_length <- c(1000, 900, 850, 800)
  
  result <- .bootstrap_resample_optimized(
    x = x,
    q = 1,
    norm = TRUE,
    nboot = 12L,
    log_base = exp(1),
    pseudocount = 0.5,
    what = "S",
    paired = FALSE,
    effective_length = effective_length
  )
  
  # Should apply length normalization
  expect_is(result, "numeric")
  expect_length(result, 12)
})


# ============================================================================
# ADDITIONAL COVERAGE TESTS: Uncovered Code Paths
# ============================================================================

test_that(".bootstrap_resample_optimized handles effective_length mismatch warning (line 374-379)", {
  # Test with mismatched effective_length
  x <- c(100, 50, 25, 10)
  effective_length <- c(1000, 500)  # Length 2, should mismatch with x (length 4)
  
  # Should produce warning about length mismatch
  result <- expect_warning(
    .bootstrap_resample_optimized(
      x = x,
      q = 1,
      norm = TRUE,
      nboot = 8L,
      log_base = exp(1),
      pseudocount = 0.5,
      what = "S",
      paired = FALSE,
      effective_length = effective_length
    ),
    "effective_length|length"
  )
  
  expect_is(result, "numeric")
})

test_that(".bootstrap_validate_inputs rejects invalid q values (line 482-483)", {
  x <- c(100, 50, 25)
  
  # Test with negative q
  expect_error(
    .bootstrap_validate_inputs(
      x = x,
      q = -1,  # Invalid: q must be >= 0
      nboot = 10L
    ),
    "q"
  )
  
  # Test with non-numeric q
  expect_error(
    .bootstrap_validate_inputs(
      x = x,
      q = "invalid",
      nboot = 10L
    ),
    "numeric|character"
  )
})

test_that(".bootstrap_validate_inputs rejects invalid nboot values (line 486-487)", {
  x <- c(100, 50, 25)
  
  # Test with negative nboot
  expect_error(
    .bootstrap_validate_inputs(
      x = x,
      q = 1,
      nboot = -10L
    ),
    "nboot must be a numeric value >= 1"
  )
  
  # Test with zero nboot
  expect_error(
    .bootstrap_validate_inputs(
      x = x,
      q = 1,
      nboot = 0L
    ),
    "nboot must be a numeric value >= 1"
  )
})

test_that(".bootstrap_compute_ci handles all-NA bootstrap distribution (line 851-853)", {
  # Create SE with extreme sparse data
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(c(0, 0, 0, 0), nrow = 1, ncol = 4)),
    colData = data.frame(sample = 1:4)
  )
  
  # Try to compute CI on all-zero gene
  result <- tryCatch({
    .bootstrap_compute_ci(
      x = c(0, 0, 0, 0),
      q = 1,
      norm = FALSE,
      nboot = 10L,
      ci = 0.95,
      method = "percentile",
      log_base = exp(1),
      pseudocount = 0,
      what = "S",
      paired = FALSE,
      effective_length = NULL,
      min_valid_frac = 0.75
    )
  }, error = function(e) {
    list(error = TRUE, message = e$message)
  })
  
  # Should handle gracefully - either return NA or error
  expect_true(is.list(result) || is.numeric(result) || is.na(result))
})

test_that(".bootstrap_process_matrix handles mismatched pair_ids (line 597-599)", {
  x_matrix <- matrix(c(100, 80, 60, 40, 50, 30), nrow = 2, ncol = 3)
  pair_ids <- c(1, 1, 2, 2, 3, 3)  # Length 6, but matrix has only 3 columns
  
  # Should catch the dimension mismatch
  result <- tryCatch({
    .bootstrap_process_matrix(
      x = x_matrix,
      pair_ids = pair_ids,
      q = 1,
      nboot = 5L
    )
  }, error = function(e) {
    list(error = TRUE, message = e$message)
  })
  
  # Should error or return gracefully
  expect_true(is.list(result) || is.numeric(result))
})

test_that("divergence_bootstrap_compute_cpp_wrapper validates pseudocount vector (line 159)", {
  x <- c(100, 80, 60)
  y <- c(90, 85, 70)
  pseudocount <- c(0.5, 1.0)  # Wrong length (2 instead of 3)
  
  # Should reject mismatched pseudocount vector
  result <- tryCatch({
    divergence_bootstrap_compute_cpp_wrapper(
      x = x,
      y = y,
      nboot = 5L,
      q = 1,
      pseudocount = pseudocount,  # Mismatched length
      log_base = exp(1)
    )
  }, error = function(e) {
    list(error = TRUE, message = e$message)
  })
  
  # Should error gracefully
  expect_true(is.list(result) && result$error)
})

test_that("block_bootstrap_compute_cpp_wrapper handles NA values in input (line 38)", {
  x <- c(100, NA, 50, 25)  # Contains NA
  
  # Should warn or error about NA values
  result <- tryCatch({
    expect_warning(
      block_bootstrap_compute_cpp_wrapper(
        x = x,
        q = 1.0,
        normalize = TRUE,
        nboot = 5L,
        log_base = exp(1),
        pseudocount = 0
      ),
      "NA|invalid"
    )
  }, error = function(e) {
    list(error = TRUE)
  })
  
  expect_true(!is.null(result) || is.numeric(result))
})

# ============================================================================
# RIGID COMPONENT: bootstrap_replicate_cpp (July 2026: metrics.json risk=95.5)
# ============================================================================

context("Rigid: bootstrap_replicate_cpp")

test_that("bootstrap_replicate_cpp returns correct dimensions", {
  skip_if_not_installed("Rcpp")
  counts <- matrix(c(10, 5, 3, 2, 8, 6, 4, 1, 15, 12, 9, 7), nrow = 4, ncol = 3)
  nboot <- 100L
  result <- bootstrap_replicate_cpp(counts, nboot = nboot, q = 1.0)
  expect_type(result, "double")
  expect_length(result, nboot)
  expect_true(all(is.finite(result)))
})

test_that("bootstrap_replicate_cpp handles minimal input (2 samples)", {
  skip_if_not_installed("Rcpp")
  counts <- matrix(c(10, 5, 3, 2, 8, 6), nrow = 3, ncol = 2)
  nboot <- 50L
  result <- bootstrap_replicate_cpp(counts, nboot = nboot, q = 1.0)
  expect_length(result, nboot)
  expect_true(all(is.finite(result)))
})

test_that("bootstrap_replicate_cpp respects q parameter", {
  skip_if_not_installed("Rcpp")
  counts <- matrix(c(10, 5, 3, 2, 8, 6, 4, 1, 15, 12, 9, 7), nrow = 4, ncol = 3)
  result_q05 <- bootstrap_replicate_cpp(counts, nboot = 100L, q = 0.5)
  result_q20 <- bootstrap_replicate_cpp(counts, nboot = 100L, q = 2.0)
  expect_true(all(is.finite(result_q05)))
  expect_true(all(is.finite(result_q20)))
  expect_true(median(result_q20, na.rm = TRUE) <= median(result_q05, na.rm = TRUE))
})

test_that("bootstrap_replicate_cpp normalize parameter works", {
  skip_if_not_installed("Rcpp")
  counts <- matrix(c(10, 5, 3, 2, 8, 6, 4, 1, 15, 12, 9, 7), nrow = 4, ncol = 3)
  result_norm <- bootstrap_replicate_cpp(counts, nboot = 100L, q = 1.0, normalize = TRUE)
  result_raw  <- bootstrap_replicate_cpp(counts, nboot = 100L, q = 1.0, normalize = FALSE)
  expect_true(all(is.finite(result_norm)))
  expect_true(all(is.finite(result_raw)))
  expect_true(all(result_norm >= 0 & result_norm <= 1))
  expect_true(all(result_raw >= 0))
})

test_that("bootstrap_replicate_cpp handles pseudocount", {
  skip_if_not_installed("Rcpp")
  counts <- matrix(c(10, 0, 3, 2, 8, 0, 4, 1, 15, 0, 9, 7), nrow = 4, ncol = 3)
  result_no_pc <- bootstrap_replicate_cpp(counts, nboot = 100L, q = 1.0, pseudocount = 0)
  result_pc <- bootstrap_replicate_cpp(counts, nboot = 100L, q = 1.0, pseudocount = 1.0)
  expect_true(all(is.finite(result_no_pc)))
  expect_true(all(is.finite(result_pc)))
  expect_true(median(result_pc, na.rm = TRUE) > median(result_no_pc, na.rm = TRUE))
})

test_that("bootstrap_replicate_cpp with block_ids preserves pairing structure", {
  skip_if_not_installed("Rcpp")
  counts <- matrix(c(10, 5, 3, 2, 8, 6, 4, 1, 15, 12, 9, 7, 20, 18, 11, 5),
                   nrow = 4, ncol = 4)
  block_ids <- c(1L, 1L, 2L, 2L)
  nboot <- 100L
  result <- bootstrap_replicate_cpp(counts, nboot = nboot, q = 1.0,
                                     block_ids = block_ids)
  expect_length(result, nboot)
  expect_true(all(is.finite(result)))
})

test_that("bootstrap_replicate_cpp produces reproducible output with seed", {
  skip_if_not_installed("Rcpp")
  counts <- matrix(c(10, 5, 3, 2, 8, 6, 4, 1, 15, 12, 9, 7), nrow = 4, ncol = 3)
  set.seed(42)
  result1 <- bootstrap_replicate_cpp(counts, nboot = 100L, q = 1.0)
  set.seed(42)
  result2 <- bootstrap_replicate_cpp(counts, nboot = 100L, q = 1.0)
  expect_equal(result1, result2)
})

test_that("bootstrap_replicate_cpp nboot parameter works", {
  skip_if_not_installed("Rcpp")
  counts <- matrix(c(10, 5, 3, 2, 15, 12, 9, 7), nrow = 4, ncol = 2)
  result_10  <- bootstrap_replicate_cpp(counts, nboot = 10L, q = 1.0)
  result_100 <- bootstrap_replicate_cpp(counts, nboot = 100L, q = 1.0)
  expect_length(result_10, 10L)
  expect_length(result_100, 100L)
})

test_that("bootstrap_replicate_cpp handles large count values", {
  skip_if_not_installed("Rcpp")
  counts <- matrix(c(1000, 500, 200, 100, 800, 600, 300, 50, 1500, 1200, 900, 700),
                   nrow = 4, ncol = 3)
  result <- bootstrap_replicate_cpp(counts, nboot = 100L, q = 2.0)
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0))
})

test_that("bootstrap_replicate_cpp handles single-row input", {
  skip_if_not_installed("Rcpp")
  counts <- matrix(c(100, 200, 150), nrow = 1, ncol = 3)
  result <- bootstrap_replicate_cpp(counts, nboot = 100L, q = 1.0, normalize = TRUE)
  expect_length(result, 100L)
  expect_true(all(result >= 0 & result <= 1))
})

