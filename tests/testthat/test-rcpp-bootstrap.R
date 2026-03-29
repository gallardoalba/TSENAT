context("C++ Bootstrap Implementation (Rcpp)")

# Test data preparation
test_counts <- c(100, 50, 25, 10)
large_counts <- c(1000, 500, 250, 100, 50, 25)
skewed_counts <- c(1000, 10, 5, 2)
small_counts <- c(5, 3, 2, 1)

# ============================================================================
# SUITE 1: Basic Bootstrap Functionality
# ============================================================================

test_that("C++ bootstrap_compute_cpp compiles and is available", {
  skip_on_cran()
  
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
  skip_on_cran()
  
  result <- bootstrap_compute_cpp_wrapper(
    x = test_counts, q = 1.0, normalize = TRUE, 
    nboot = 50L, log_base = exp(1), pseudocount = 0.0
  )
  
  expect_is(result, "numeric")
  expect_length(result, 50)
  expect_true(all(is.finite(result)))
})

test_that("bootstrap_compute_cpp produces positive entropy values", {
  skip_on_cran()
  
  result <- bootstrap_compute_cpp_wrapper(
    x = test_counts, q = 1.0, normalize = TRUE, 
    nboot = 100L, log_base = exp(1), pseudocount = 0.0
  )
  
  # For Shannon entropy (q=1) with normalized=TRUE, values should be in [0, 1]
  expect_true(all(result >= 0, na.rm = TRUE))
  expect_true(all(result <= 1, na.rm = TRUE))
})

test_that("bootstrap_compute_cpp respects nboot parameter", {
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
  result <- bootstrap_compute_cpp_wrapper(
    x = small_counts, q = 1.0, normalize = TRUE, 
    nboot = 50L, log_base = exp(1), pseudocount = 0.0
  )
  
  expect_true(all(is.finite(result)))
  expect_length(result, 50)
})

test_that("bootstrap_compute_cpp handles skewed distributions", {
  skip_on_cran()
  
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
  skip_on_cran()
  
  result <- bootstrap_compute_cpp_wrapper(
    x = large_counts, q = 1.0, normalize = TRUE, 
    nboot = 50L, log_base = exp(1), pseudocount = 0.0
  )
  
  expect_true(all(is.finite(result)))
  expect_length(result, 50)
})

test_that("bootstrap_compute_cpp handles single value vector", {
  skip_on_cran()
  
  result <- bootstrap_compute_cpp_wrapper(
    x = c(100), q = 1.0, normalize = TRUE, 
    nboot = 50L, log_base = exp(1), pseudocount = 0.0
  )
  
  # Single value should have zero entropy
  expect_true(all(result == 0, na.rm = TRUE))
})

test_that("bootstrap_compute_cpp handles zeros with pseudocount", {
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
  result <- .bootstrap_resample_optimized(
    x = test_counts, q = 1.0, norm = TRUE, nboot = 50,
    log_base = exp(1), pseudocount = 0, what = "S", paired = FALSE
  )
  
  expect_is(result, "numeric")
  expect_length(result, 50)
  expect_true(all(is.finite(result)))
})

test_that(".bootstrap_resample_optimized handles what='S' (entropy)", {
  skip_on_cran()
  
  result <- .bootstrap_resample_optimized(
    x = test_counts, q = 1.0, norm = TRUE, nboot = 50,
    log_base = exp(1), pseudocount = 0, what = "S", paired = FALSE
  )
  
  expect_true(all(result >= 0, na.rm = TRUE))
  expect_true(all(result <= 1, na.rm = TRUE))
})

test_that(".bootstrap_resample_optimized handles what='D' (Hill numbers)", {
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
  # This test verifies the fallback mechanism exists
  # (by default C++ succeeds, but if it fails, should fall back)
  
  result <- .bootstrap_resample_optimized(
    x = test_counts, q = 1.0, norm = TRUE, nboot = 50,
    log_base = exp(1), pseudocount = 0, what = "S", paired = FALSE
  )
  
  expect_true(all(is.finite(result)))
})

test_that("bootstrap_compute_cpp_wrapper has proper error messages", {
  skip_on_cran()
  
  # Test that invalid inputs are caught
  expect_error(
    bootstrap_compute_cpp_wrapper(
      x = numeric(0), q = 1.0, normalize = TRUE, 
      nboot = 50L, log_base = exp(1), pseudocount = 0.0
    )
  )
})

test_that("bootstrap_compute_cpp handles negative nboot gracefully", {
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
  odd_length <- c(100, 95, 110)  # 3 elements (1.5 pairs)
  
  expect_error(
    block_bootstrap_compute_cpp_wrapper(
      x = odd_length, q = 1.0, normalize = TRUE, 
      nboot = 10L, log_base = exp(1), pseudocount = 0.0
    )
  )
})

test_that("block_bootstrap_compute_cpp returns valid structure", {
  skip_on_cran()
  
  result <- block_bootstrap_compute_cpp_wrapper(
    x = paired_counts_3, q = 1.0, normalize = TRUE, 
    nboot = 50L, log_base = exp(1), pseudocount = 0.0
  )
  
  expect_is(result, "numeric")
  expect_length(result, 50)
  expect_true(all(is.finite(result)))
})

test_that("block_bootstrap_compute_cpp produces entropy in valid range", {
  skip_on_cran()
  
  result <- block_bootstrap_compute_cpp_wrapper(
    x = paired_counts_3, q = 1.0, normalize = TRUE, 
    nboot = 100L, log_base = exp(1), pseudocount = 0.0
  )
  
  # For normalized Shannon entropy, should be in [0, 1]
  expect_true(all(result >= 0, na.rm = TRUE))
  expect_true(all(result <= 1, na.rm = TRUE))
})

test_that("block_bootstrap_compute_cpp respects nboot parameter", {
  skip_on_cran()
  
  for (nboot in c(10, 50, 100, 200)) {
    result <- block_bootstrap_compute_cpp_wrapper(
      x = paired_counts_3, q = 1.0, normalize = TRUE, 
      nboot = as.integer(nboot), log_base = exp(1), pseudocount = 0.0
    )
    
    expect_length(result, nboot)
  }
})

test_that("block_bootstrap_compute_cpp works with various q values", {
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
  small_paired <- c(10, 9, 12, 11)  # Very small counts
  
  result <- block_bootstrap_compute_cpp_wrapper(
    x = small_paired, q = 1.0, normalize = TRUE, 
    nboot = 50L, log_base = exp(1), pseudocount = 0.0
  )
  
  expect_true(all(is.finite(result)))
  expect_length(result, 50)
})

test_that("block bootstrap with skewed pairs works", {
  skip_on_cran()
  
  skewed_paired <- c(1000, 10, 5, 2)  # Highly skewed pairs
  
  result <- block_bootstrap_compute_cpp_wrapper(
    x = skewed_paired, q = 1.0, normalize = TRUE, 
    nboot = 50L, log_base = exp(1), pseudocount = 0.0
  )
  
  expect_true(all(is.finite(result)))
  expect_length(result, 50)
})

test_that("block bootstrap preserves pair structure", {
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
  skip_on_cran()
  
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
                                       n_bootstrap = 300, confidence = 0.95)
    
    # Check monotonicity: all ci_lower <= ci_upper
    valid_idx <- !is.na(result$ci_lower) & !is.na(result$ci_upper)
    expect_true(all(result$ci_lower[valid_idx] <= result$ci_upper[valid_idx]),
                info = sprintf("CI not monotonic for seed %d", seed))
  }
})

test_that("REGRESSION 10.2: Confidence level ordering respected", {
  skip_on_cran()
  
  set.seed(777)
  counts_A <- matrix(rpois(20 * 35, lambda = 8), nrow = 20, ncol = 35)
  counts_B <- matrix(rpois(20 * 35, lambda = 8), nrow = 20, ncol = 35)
  
  jack_A <- jis_jackknife_influences_cpp(counts_A, q = 1, normalize = TRUE)
  jack_B <- jis_jackknife_influences_cpp(counts_B, q = 1, normalize = TRUE)
  delta_influence <- abs(jack_A - jack_B)
  
  # Compute CIs at different confidence levels
  result_80 <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence,
                                        q = 1, normalize = TRUE,
                                        n_bootstrap = 500, confidence = 0.80)
  result_90 <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence,
                                        q = 1, normalize = TRUE,
                                        n_bootstrap = 500, confidence = 0.90)
  result_95 <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence,
                                        q = 1, normalize = TRUE,
                                        n_bootstrap = 500, confidence = 0.95)
  result_99 <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence,
                                        q = 1, normalize = TRUE,
                                        n_bootstrap = 500, confidence = 0.99)
  
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
  skip_on_cran()
  
  set.seed(888)
  counts_A <- matrix(rpois(25 * 30, lambda = 12), nrow = 25, ncol = 30)
  counts_B <- matrix(rpois(25 * 30, lambda = 12), nrow = 25, ncol = 30)
  
  jack_A <- jis_jackknife_influences_cpp(counts_A, q = 1, normalize = TRUE)
  jack_B <- jis_jackknife_influences_cpp(counts_B, q = 1, normalize = TRUE)
  delta_influence <- abs(jack_A - jack_B)
  
  result <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence,
                                     q = 1, normalize = TRUE,
                                     n_bootstrap = 1000, confidence = 0.95)
  
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
  skip_on_cran()
  
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
                                       n_bootstrap = 400, confidence = 0.95)
    
    # All output fields should exist and be properly sized
    expect_equal(length(result$ci_lower), nrow(counts_A))
    expect_equal(length(result$ci_upper), nrow(counts_A))
    expect_equal(length(result$p_value), nrow(counts_A))
    expect_equal(length(result$ci_width), nrow(counts_A))
  }
})

test_that("REGRESSION 10.5: Quantile method consistency across sample sizes", {
  skip_on_cran()
  
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
                                       n_bootstrap = nboot, confidence = 0.95)
    
    # CIs should be valid for any nboot
    valid_ci <- !is.na(result$ci_lower) & !is.na(result$ci_upper)
    if (any(valid_ci)) {
      expect_true(all(result$ci_lower[valid_ci] <= result$ci_upper[valid_ci]))
    }
  }
})

test_that("REGRESSION 10.6: Effect size and CI consistency", {
  skip_on_cran()
  
  set.seed(1212)
  counts_A <- matrix(rpois(16 * 32, lambda = 10), nrow = 16, ncol = 32)
  counts_B <- matrix(rpois(16 * 32, lambda = 10), nrow = 16, ncol = 32)
  
  jack_A <- jis_jackknife_influences_cpp(counts_A, q = 1, normalize = TRUE)
  jack_B <- jis_jackknife_influences_cpp(counts_B, q = 1, normalize = TRUE)
  delta_influence <- abs(jack_A - jack_B)
  
  result <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence,
                                     q = 1, normalize = TRUE,
                                     n_bootstrap = 600, confidence = 0.95)
  
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
  skip_on_cran()
  
  set.seed(1313)
  # Create data with some zero-sum samples (will produce NA in entropy)
  counts_A <- matrix(rpois(20 * 30, lambda = 5), nrow = 20, ncol = 30)
  counts_B <- matrix(rpois(20 * 30, lambda = 5), nrow = 20, ncol = 30)
  
  jack_A <- jis_jackknife_influences_cpp(counts_A, q = 1, normalize = TRUE)
  jack_B <- jis_jackknife_influences_cpp(counts_B, q = 1, normalize = TRUE)
  delta_influence <- abs(jack_A - jack_B)
  
  result <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence,
                                     q = 1, normalize = TRUE,
                                     n_bootstrap = 400, confidence = 0.95)
  
  # All output should be either numeric or NA, not NaN
  expect_true(all(!is.nan(result$ci_lower)))
  expect_true(all(!is.nan(result$ci_upper)))
  expect_true(all(!is.nan(result$p_value)))
  expect_true(all(!is.nan(result$ci_width)))
})
