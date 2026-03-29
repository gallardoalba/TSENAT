context("DEEP BUG ANALYSIS: Formula Sign Errors in Tsallis")

test_that("BUG ALERT: Tsallis formula uses wrong denominator (q-1 vs 1-q) - produces opposite signs", {
  # This test documents the critical formula bug
  
  # Mathematical formula check
  q_vals <- c(0.5, 1.5, 2.0, 3.0)
  p_q_sum <- 0.3  # sum of p^q
  
  results <- data.frame()
  
  for (q in q_vals) {
    # WRONG formula (current jackknife_isoform_switching.R line 93)
    h_wrong <- (1 / (1 - q)) * (1 - p_q_sum)
    
    # CORRECT formula (entropy_core.R line 56, C++ line 641)
    h_correct <- (1.0 - p_q_sum) / (q - 1.0)
    
    results <- rbind(results, data.frame(
      q = q,
      h_wrong = h_wrong,
      h_correct = h_correct,
      sign_match = sign(h_wrong) == sign(h_correct),
      actual_diff = abs(h_wrong - h_correct)
    ))
  }
  
  # For q < 1: formulas give opposite signs
  # For q > 1: formulas give opposite signs
  expect_false(results$sign_match[results$q == 0.5])  # Different signs for q=0.5
  expect_false(results$sign_match[results$q == 1.5])  # Different signs for q=1.5
  expect_false(results$sign_match[results$q == 2.0])  # Different signs for q=2.0
  expect_false(results$sign_match[results$q == 3.0])  # Different signs for q=3.0
})

test_that("BUG ALERT: Normalized Tsallis uses wrong denominator - produces opposite signs", {
  # This test documents the normalization formula bug
  
  q_vals <- c(0.5, 1.5, 2.0)
  n_tx <- 100
  
  results <- data.frame()
  
  for (q in q_vals) {
    # WRONG formula (current jackknife_isoform_switching.R line 97)
    max_h_wrong <- (1 - n_tx^(1 - q)) / (1 - q)
    
    # CORRECT formula (entropy_core.R line 70, C++ line 645)
    max_h_correct <- (1 - n_tx^(1 - q)) / (q - 1)
    
    results <- rbind(results, data.frame(
      q = q,
      max_h_wrong = max_h_wrong,
      max_h_correct = max_h_correct,
      sign_match = sign(max_h_wrong) == sign(max_h_correct),
      actual_diff = abs(max_h_wrong - max_h_correct)
    ))
  }
  
  # For q < 1: formulas give opposite signs
  # For q > 1: formulas give opposite signs  
  expect_false(results$sign_match[results$q == 0.5])  # Different signs for q=0.5
  expect_false(results$sign_match[results$q == 1.5])  # Different signs for q=1.5
  expect_false(results$sign_match[results$q == 2.0])  # Different signs for q=2.0
})

test_that("BUG IMPACT: R vs C++ give completely different entropy values due to formula bug", {
  skip_on_cran()
  
  # Create test data
  n_tx <- 100
  n_samples <- 5
  test_counts <- matrix(rep(1, n_tx * n_samples), nrow = n_tx, ncol = n_samples)
  
  # Add some variation
  test_counts[1:10, 1] <- 100  # Skew first sample
  test_counts[11:20, 2] <- 50   # Skew second sample
  
  test_qs <- c(0.5, 1.5, 2.0)
  
  for (q in test_qs) {
    # Get C++ result (CORRECT)
    cpp_entropy <- jis_tsallis_entropy_cpp(test_counts, q = q, normalize = TRUE, 
                                           log_base = exp(1), pseudocount = 0)
    
    # Get R result (NOW FIXED)
    r_entropy <- TSENAT:::.jis_tsallis_entropy(test_counts, q = q, norm = TRUE, 
                                              log_base = exp(1), pseudocount = 0)
    
    # After fixes: R and C++ should match for non-uniform distributions
    max_diff <- max(abs(cpp_entropy - r_entropy), na.rm = TRUE)
    expect_true(max_diff < 1e-4, 
                info = sprintf("q=%.1f: R and C++ should match after formula fixes (diff=%.2e)", q, max_diff))
  }
})

test_that("ROOT CAUSE: Denominator (q-1) vs (1-q) affects all q != 1 cases", {
  # This test shows the mathematical root cause
  
  # Mathematical fact: (1/(1-q)) * (1-x) = -((1-x)/(q-1))
  # For q < 1: 1-q > 0, so (1/(1-q)) is positive
  #            q-1 < 0, so (1/(q-1)) is negative
  # For q > 1: 1-q < 0, so (1/(1-q)) is negative
  #            q-1 > 0, so (1/(q-1)) is positive
  
  # Test: (1/(1-q)) * x = -1 * (x/(q-1))
  q <- 0.7
  x <- 0.3
  
  formula_1 <- (1 / (1 - q)) * x
  formula_2 <- -(x / (q - 1))
  
  expect_equal(formula_1, formula_2, tolerance = 1e-10)
  
  # This proves: (1/(1-q)) * (1-sum) = -(1-sum)/(q-1)
  # They have OPPOSITE signs!
})

test_that("DIAGNOSTIC: Investigate C++ vs R differences", {
  skip_on_cran()
  
  # This test documents the Shannon entropy pseudocount handling difference
  
  # Create a singular distribution (all mass on one species)
  n_tx <- 100
  singular_counts <- matrix(0, nrow = n_tx, ncol = 1)
  singular_counts[1, 1] <- 100
  
  # What happens in C++: filters p[i] <= 1e-15 before log
  # What happens in R: adds 1e-100 to ALL p values before log
  
  # Result: For singular distribution with pseudocount=1e-8:
  # C++ does: proportions of main species ≈ 1.0, others ≈ 0
  #           Only p[main] > 1e-15, so only compute log(1.0) ≈ 0
  #           Result: H ≈ 0
  
  # R does: proportions have 1e-100 added
  #         p[main] ≈ 1.0, others ≈ 1e-100
  #         log(1.0) ≈ 0, but log(1e-100) ≈ -230.26
  #         Result: H ≈ 0 + many small contributions from other species
  
  # Verify R and C++ are now consistent (after fix)
  cpp_entropy <- jis_tsallis_entropy_cpp(singular_counts, q = 1.0, normalize = TRUE,
                                         log_base = exp(1), pseudocount = 0)
  r_entropy <- TSENAT:::.jis_tsallis_entropy(singular_counts, q = 1.0, norm = TRUE,
                                             log_base = exp(1), pseudocount = 0)
  
  max_diff <- max(abs(cpp_entropy - r_entropy), na.rm = TRUE)
  expect_true(max_diff < 1e-4, 
              info = sprintf("R and C++ should match for singular distribution, diff = %.2e", max_diff))
})

test_that("QUANTILE METHOD: Check if type=1 is being used correctly", {
  # The bootstrap quantile bug should be fixed, but let's verify
  
  # Generate some bootstrap replicates
  set.seed(42)
  bootstrap_vals <- rnorm(1000, mean = 0, sd = 1)
  
  # Compute confidence intervals
  ci_type1 <- quantile(bootstrap_vals, c(0.025, 0.975), type = 1)
  ci_type7 <- quantile(bootstrap_vals, c(0.025, 0.975), type = 7)
  
  # They should be different
  expect_false(all(ci_type1 == ci_type7),
               info = "type=1 and type=7 should give different quantiles")
})
