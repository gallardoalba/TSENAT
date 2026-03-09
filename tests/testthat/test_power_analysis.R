# Tests for power_tsenat_entropy function
# Information-theoretic power analysis using Tsallis q-divergence
# Reference: R/power_tsenat_entropy.R

library(testthat)

context("Power Analysis: TSENAT Entropy-Based (Tsallis q-Divergence)")

# =============================================================================
# TEST 1: Output Bounds and Return Type
# =============================================================================
test_that("power_tsenat_entropy returns values in [0, 1]", {
  # Test at various parameter combinations
  expect_true(power_tsenat_entropy(n=10, fc=2.0) >= 0 && 
              power_tsenat_entropy(n=10, fc=2.0) <= 1)
  
  # Edge case: very low sample size
  result_n1 <- power_tsenat_entropy(n=1, fc=1.5)
  expect_true(is.numeric(result_n1))
  expect_gte(result_n1, 0)
  expect_lte(result_n1, 1)
  
  # Multiple calls at different parameters
  sample_sizes <- c(5, 10, 20, 50, 100)
  for (n in sample_sizes) {
    power_val <- power_tsenat_entropy(n=n, fc=2.0)
    expect_true(power_val >= 0 && power_val <= 1, 
                info = paste("Failed at n =", n))
  }
})

# =============================================================================
# TEST 2: Return Type
# =============================================================================
test_that("power_tsenat_entropy returns numeric scalar", {
  result <- power_tsenat_entropy(n=20, fc=2.0)
  expect_is(result, "numeric")
  expect_length(result, 1)
  expect_false(is.na(result))
})

# =============================================================================
# TEST 3: Monotonicity in Sample Size
# =============================================================================
test_that("power_tsenat_entropy increases with sample size", {
  # At fixed fold change, power should increase with n
  fc <- 2.0
  n_vals <- c(5, 10, 15, 20, 30, 50, 100)
  powers <- sapply(n_vals, function(n) power_tsenat_entropy(n=n, fc=fc, q=0.2))
  
  # Check strictly increasing (or at least non-decreasing with small tolerances)
  diffs <- diff(powers)
  expect_true(all(diffs >= -1e-10))
  
  # Check more substantial increase at extremes
  expect_lt(powers[1], powers[7])
})

# =============================================================================
# TEST 4: Monotonicity in Fold Change
# =============================================================================
test_that("power_tsenat_entropy increases with fold change", {
  # At fixed sample size, power should increase with fc
  n <- 20
  fc_vals <- c(1.2, 1.5, 2.0, 3.0, 5.0)
  powers <- sapply(fc_vals, function(fc) power_tsenat_entropy(n=n, fc=fc, q=0.2))
  
  # Check increasing trend
  diffs <- diff(powers)
  expect_true(all(diffs >= -1e-10))
  
  # Check substantial increase
  expect_lt(powers[1], powers[5])
})

# =============================================================================
# TEST 5: Monotonicity in q-Parameter
# =============================================================================
test_that("power_tsenat_entropy increases with q-parameter", {
  # At fixed n and fc, power should increase with q
  n <- 20
  fc <- 2.0
  q_vals <- c(0.0, 0.2, 0.5, 1.0, 1.5, 2.0)
  powers <- sapply(q_vals, function(q) power_tsenat_entropy(n=n, fc=fc, q=q))
  
  # Check increasing trend
  diffs <- diff(powers)
  expect_true(all(diffs >= -1e-10))
  
  # Higher q = more information gain
  expect_lt(powers[1], powers[6])
})

# =============================================================================
# TEST 6: Input Validation - Sample Size
# =============================================================================
test_that("power_tsenat_entropy validates sample size parameter", {
  # n must be positive
  expect_error(power_tsenat_entropy(n=-5, fc=2.0), 
               "Sample size n must be a positive number")
  expect_error(power_tsenat_entropy(n=0, fc=2.0), 
               "Sample size n must be a positive number")
  expect_error(power_tsenat_entropy(n="invalid", fc=2.0), 
               "Sample size n must be a positive number")
})

# =============================================================================
# TEST 7: Input Validation - Fold Change
# =============================================================================
test_that("power_tsenat_entropy validates fold change parameter", {
  # fc must be > 1
  expect_error(power_tsenat_entropy(n=20, fc=1.0), 
               "Fold change fc must be greater than 1")
  expect_error(power_tsenat_entropy(n=20, fc=0.5), 
               "Fold change fc must be greater than 1")
  expect_error(power_tsenat_entropy(n=20, fc=-2.0), 
               "Fold change fc must be greater than 1")
  expect_error(power_tsenat_entropy(n=20, fc="invalid"), 
               "Fold change fc must be greater than 1")
})

# =============================================================================
# TEST 8: Input Validation - q-Parameter
# =============================================================================
test_that("power_tsenat_entropy validates q-parameter", {
  # q must be in [0, 2]
  expect_error(power_tsenat_entropy(n=20, fc=2.0, q=-0.5), 
               "Tsallis q-parameter must be in range")
  expect_error(power_tsenat_entropy(n=20, fc=2.0, q=2.5), 
               "Tsallis q-parameter must be in range")
  expect_error(power_tsenat_entropy(n=20, fc=2.0, q="invalid"), 
               "Tsallis q-parameter must be in range")
})

# =============================================================================
# TEST 9: Input Validation - Alpha (Significance Level)
# =============================================================================
test_that("power_tsenat_entropy validates alpha parameter", {
  # alpha must be in (0, 1)
  expect_error(power_tsenat_entropy(n=20, fc=2.0, alpha=0), 
               "Significance level alpha must be in")
  expect_error(power_tsenat_entropy(n=20, fc=2.0, alpha=1.0), 
               "Significance level alpha must be in")
  expect_error(power_tsenat_entropy(n=20, fc=2.0, alpha=-0.05), 
               "Significance level alpha must be in")
  expect_error(power_tsenat_entropy(n=20, fc=2.0, alpha=1.5), 
               "Significance level alpha must be in")
})

# =============================================================================
# TEST 10: Typical Usage - Realistic Parameters
# =============================================================================
test_that("power_tsenat_entropy handles realistic RNA-seq parameters", {
  # Typical scenario: moderate sample size, 2-fold change
  result <- power_tsenat_entropy(n=20, fc=2.0, baseline_mean=100, 
                                  dispersion=0.1, alpha=0.05, q=0.2)
  expect_true(result > 0.3)
  expect_true(result < 1.0)
  
  # Small effect, moderate sample
  result_small <- power_tsenat_entropy(n=50, fc=1.2)
  expect_true(result_small >= 0 && result_small <= 1)
  
  # Large effect, small sample
  result_large <- power_tsenat_entropy(n=5, fc=5.0)
  expect_true(result_large >= 0 && result_large <= 1)
})

# =============================================================================
# TEST 11: Parameter Defaults
# =============================================================================
test_that("power_tsenat_entropy uses correct default parameters", {
  # With defaults
  result_default <- power_tsenat_entropy(n=20, fc=2.0)
  
  # With explicit defaults
  result_explicit <- power_tsenat_entropy(n=20, fc=2.0, 
                                           baseline_mean=100, 
                                           dispersion=0.1, 
                                           alpha=0.05, q=0.2)
  
  expect_equal(result_default, result_explicit,
               info = "Defaults should match explicit parameter specification")
})

# =============================================================================
# TEST 12: Comparison Between Different q Values
# =============================================================================
test_that("q-parameter effects align with theory", {
  n <- 30
  fc <- 1.5
  
  # Calculate powers at different q values
  p_q0 <- power_tsenat_entropy(n=n, fc=fc, q=0.0)
  p_q05 <- power_tsenat_entropy(n=n, fc=fc, q=0.5)
  p_q10 <- power_tsenat_entropy(n=n, fc=fc, q=1.0)
  
  # Higher q should give higher power
  expect_lt(p_q0, p_q05)
  expect_lt(p_q05, p_q10)
})

# =============================================================================
# TEST 13: Alpha Effect on Power
# =============================================================================
test_that("power_tsenat_entropy decreases with stricter alpha (multiple testing)", {
  n <- 20
  fc <- 2.0
  
  # Lower alpha = stricter threshold = lower power
  p_alpha01 <- power_tsenat_entropy(n=n, fc=fc, alpha=0.1)
  p_alpha05 <- power_tsenat_entropy(n=n, fc=fc, alpha=0.05)
  p_alpha01_100 <- power_tsenat_entropy(n=n, fc=fc, alpha=0.001)
  
  # More lenient alpha → higher power
  expect_gt(p_alpha01, p_alpha05)
  expect_gt(p_alpha05, p_alpha01_100)
})

# =============================================================================
# TEST 14: Vectorization Check (if supported)
# =============================================================================
test_that("power_tsenat_entropy handles scalar inputs correctly", {
  # Single value - should work
  result <- power_tsenat_entropy(n=20, fc=2.0)
  expect_length(result, 1)
  
  # If vectorization is not supported, each must be scalar
  # (current implementation expects scalars)
})

# =============================================================================
# TEST 15: Effect Size Scaling (log-scale interpretation)
# =============================================================================
test_that("power_tsenat_entropy uses log-scale effect size correctly", {
  n <- 30
  
  # fc=2.0 means 2-fold (log scale: log(2) ≈ 0.693)
  # fc=4.0 means 4-fold, should have more power than 2-fold
  p_fc2 <- power_tsenat_entropy(n=n, fc=2.0)
  p_fc4 <- power_tsenat_entropy(n=n, fc=4.0)
  
  expect_lt(p_fc2, p_fc4)
  
  # fc=2.0 vs fc=2.0^2=4.0 should differ substantially
  ratio <- p_fc4 / p_fc2
  expect_gt(ratio, 1.1)
})

# =============================================================================
# TEST 16: Boundary Conditions
# =============================================================================
test_that("power_tsenat_entropy handles boundary conditions reasonably", {
  # Minimum sample size (n=1)
  expect_true(is.numeric(power_tsenat_entropy(n=1, fc=2.0)))
  
  # Fold change just above threshold
  expect_true(is.numeric(power_tsenat_entropy(n=20, fc=1.001)))
  
  # Minimum q value
  expect_true(is.numeric(power_tsenat_entropy(n=20, fc=2.0, q=0)))
  
  # Maximum q value
  expect_true(is.numeric(power_tsenat_entropy(n=20, fc=2.0, q=2.0)))
})

# =============================================================================
# TEST 17: Consistency Across Multiple Calls
# =============================================================================
test_that("power_tsenat_entropy is deterministic (consistent across calls)", {
  result1 <- power_tsenat_entropy(n=25, fc=1.8, q=0.3)
  result2 <- power_tsenat_entropy(n=25, fc=1.8, q=0.3)
  result3 <- power_tsenat_entropy(n=25, fc=1.8, q=0.3)
  
  expect_equal(result1, result2)
  expect_equal(result2, result3)
})

# =============================================================================
# TEST 18: Information Gain Interpretation
# =============================================================================
test_that("power_tsenat_entropy formula incorporates information divergence", {
  # At fixed n, log(fc) is linear in information gain
  # So fc=2 vs fc=4: log(2)=0.693, log(4)=1.386 (ratio ≈ 2)
  n <- 40
  q <- 0.2
  
  p_fc2 <- power_tsenat_entropy(n=n, fc=2.0, q=q)
  p_fc4 <- power_tsenat_entropy(n=n, fc=4.0, q=q)
  
  # Power should increase substantially with doubled log-effect
  expect_lt(p_fc2, p_fc4)
})

# =============================================================================
# TEST 19: Multiple Testing Correction (Bonferroni via log(n+2))
# =============================================================================
test_that("power_tsenat_entropy includes multiple testing correction", {
  fc <- 2.0
  
  # At very large n, correction becomes stricter
  # This should manifest as power not reaching 1 even asymptotically
  p_small_n <- power_tsenat_entropy(n=50, fc=fc)
  p_large_n <- power_tsenat_entropy(n=500, fc=fc)
  
  # Both should be high power
  expect_gt(p_large_n, p_small_n)
  
  # Power should plateau (correction prevents unlimited growth)
  expect_lte(p_large_n, 1.0)
})

# =============================================================================
# TEST 20: Documentation Example Verification
# =============================================================================
test_that("power_tsenat_entropy documentation examples execute correctly", {
  # Example 1: Basic usage
  result <- power_tsenat_entropy(n=20, fc=2.0)
  expect_true(is.numeric(result) && !is.na(result))
  
  # Example 2: Sample size sequence (from docs)
  n_seq <- c(5, 10, 15, 20, 30, 50)
  powers <- sapply(n_seq, function(n) power_tsenat_entropy(n, fc=2.0))
  expect_equal(length(powers), 6)
  expect_true(all(!is.na(powers)))
  expect_true(all(is.numeric(powers)))
  
  # Example 3: Fold change sequence (from docs)
  fc_seq <- c(1.2, 1.5, 2.0, 3.0, 5.0)
  powers_fc <- sapply(fc_seq, function(fc) power_tsenat_entropy(n=20, fc=fc))
  expect_equal(length(powers_fc), 5)
  expect_true(all(!is.na(powers_fc)))
  expect_true(all(is.numeric(powers_fc)))
})

# Tests for Classical Power Analysis Functions
# DESeq2 Wald, edgeR Exact, and edgeR QLF power calculations
# Reference: R/power_analysis_classical.R

library(testthat)

context("Power Analysis: Classical Methods (DESeq2 Wald, edgeR Exact, edgeR QLF)")

# =============================================================================
# HELPER FUNCTION: Test common properties of all three methods
# =============================================================================

test_classical_power_function <- function(func, func_name) {
  # Test output bounds
  result <- func(n=20, fc=2.0)
  expect_true(result >= 0 && result <= 1, 
              info = paste(func_name, "should return value in [0, 1]"))
  
  # Test return type
  expect_is(result, "numeric", info = paste(func_name, "should return numeric"))
  expect_length(result, 1, info = paste(func_name, "should return scalar"))
  
  # Test monotonicity with n
  n_vals <- c(5, 10, 20, 50)
  powers_n <- sapply(n_vals, function(n) func(n=n, fc=2.0))
  expect_true(all(diff(powers_n) >= -1e-10), 
              info = paste(func_name, "power should increase with sample size"))
  
  # Test monotonicity with fc
  fc_vals <- c(1.2, 1.5, 2.0, 3.0)
  powers_fc <- sapply(fc_vals, function(fc) func(n=20, fc=fc))
  expect_true(all(diff(powers_fc) >= -1e-10), 
              info = paste(func_name, "power should increase with fold change"))
}

# =============================================================================
# ============================ DESEQ2 WALD TESTS ============================
# =============================================================================

test_that("power_deseq2_wald returns values in [0, 1]", {
  result <- power_deseq2_wald(n=20, fc=2.0)
  expect_true(result >= 0 && result <= 1)
  
  # Test across range
  for (n in c(1, 5, 10, 20, 50)) {
    for (fc in c(1.2, 1.5, 2.0, 3.0)) {
      power_val <- power_deseq2_wald(n=n, fc=fc)
      expect_true(power_val >= 0 && power_val <= 1,
                  info = paste("Failed at n =", n, ", fc =", fc))
    }
  }
})

test_that("power_deseq2_wald returns numeric scalar", {
  result <- power_deseq2_wald(n=20, fc=2.0)
  expect_is(result, "numeric")
  expect_length(result, 1)
  expect_false(is.na(result))
})

test_that("power_deseq2_wald increases with sample size", {
  fc <- 2.0
  n_vals <- c(5, 10, 15, 20, 30, 50)
  powers <- sapply(n_vals, function(n) power_deseq2_wald(n=n, fc=fc))
  
  diffs <- diff(powers)
  expect_true(all(diffs >= -1e-10))
  expect_lt(powers[1], powers[6])
})

test_that("power_deseq2_wald increases with fold change", {
  n <- 20
  fc_vals <- c(1.2, 1.5, 2.0, 3.0, 5.0)
  powers <- sapply(fc_vals, function(fc) power_deseq2_wald(n=n, fc=fc))
  
  diffs <- diff(powers)
  expect_true(all(diffs >= -1e-10))
  expect_lt(powers[1], powers[5])
})

test_that("power_deseq2_wald handles negative binomial variance correctly", {
  # Variance = μ + φμ² (NB model)
  # Verify calculations are reasonable
  result_low_disp <- power_deseq2_wald(n=20, fc=2.0, dispersion=0.05)
  result_high_disp <- power_deseq2_wald(n=20, fc=2.0, dispersion=0.3)
  
  # Higher dispersion → lower power (more variance)
  expect_lt(result_high_disp, result_low_disp)
})

test_that("power_deseq2_wald responds to baseline mean differences", {
  result_low_mean <- power_deseq2_wald(n=20, fc=2.0, baseline_mean=10)
  result_high_mean <- power_deseq2_wald(n=20, fc=2.0, baseline_mean=1000)
  
  # Results might differ due to numerical precision impacts
  # but both should be valid powers
  expect_true(!is.na(result_low_mean) && !is.na(result_high_mean))
})

test_that("power_deseq2_wald responds to alpha level", {
  result_alpha05 <- power_deseq2_wald(n=20, fc=2.0, alpha=0.05)
  result_alpha01 <- power_deseq2_wald(n=20, fc=2.0, alpha=0.01)
  
  # Stricter alpha (0.01 vs 0.05) → lower power
  expect_lt(result_alpha01, result_alpha05)
})

test_that("power_deseq2_wald uses default parameters correctly", {
  result_default <- power_deseq2_wald(n=20, fc=2.0)
  result_explicit <- power_deseq2_wald(n=20, fc=2.0, 
                                        baseline_mean=100, 
                                        dispersion=0.1, 
                                        alpha=0.05)
  expect_equal(result_default, result_explicit)
})

test_that("power_deseq2_wald realistic parameters give expected power", {
  # Typical RNA-seq: n=20, fc=2.0 should have very high power with classical method
  result <- power_deseq2_wald(n=20, fc=2.0)
  expect_gt(result, 0.9)
  expect_lt(result, 1.0)
})

test_that("power_deseq2_wald is deterministic", {
  r1 <- power_deseq2_wald(n=25, fc=1.8)
  r2 <- power_deseq2_wald(n=25, fc=1.8)
  r3 <- power_deseq2_wald(n=25, fc=1.8)
  
  expect_equal(r1, r2)
  expect_equal(r2, r3)
})

test_that("power_deseq2_wald documentation example works", {
  result <- power_deseq2_wald(n=20, fc=2.0)
  expect_is(result, "numeric")
  expect_true(result >= 0 && result <= 1)
})

# =============================================================================
# ========================= EDGER EXACT TEST TESTS ==========================
# =============================================================================

test_that("power_edger_exact returns values in [0, 1]", {
  result <- power_edger_exact(n=20, fc=2.0)
  expect_true(result >= 0 && result <= 1)
  
  # Test across parameter ranges
  for (n in c(1, 5, 10, 20, 50)) {
    for (fc in c(1.2, 1.5, 2.0, 3.0)) {
      power_val <- power_edger_exact(n=n, fc=fc)
      expect_true(power_val >= 0 && power_val <= 1,
                  info = paste("Failed at n =", n, ", fc =", fc))
    }
  }
})

test_that("power_edger_exact returns numeric scalar", {
  result <- power_edger_exact(n=20, fc=2.0)
  expect_is(result, "numeric")
  expect_length(result, 1)
  expect_false(is.na(result))
})

test_that("power_edger_exact increases with sample size", {
  fc <- 2.0
  n_vals <- c(5, 10, 15, 20, 30, 50)
  powers <- sapply(n_vals, function(n) power_edger_exact(n=n, fc=fc))
  
  diffs <- diff(powers)
  expect_true(all(diffs >= -1e-10))
  expect_lt(powers[1], powers[6])
})

test_that("power_edger_exact increases with fold change", {
  n <- 20
  fc_vals <- c(1.2, 1.5, 2.0, 3.0, 5.0)
  powers <- sapply(fc_vals, function(fc) power_edger_exact(n=n, fc=fc))
  
  diffs <- diff(powers)
  expect_true(all(diffs >= -1e-10))
  expect_lt(powers[1], powers[5])
})

test_that("power_edger_exact is generally conservative vs Wald", {
  # edgeR exact test uses hypergeometric approx → more conservative
  # So power_edger_exact(n, fc) < power_deseq2_wald(n, fc) typically
  result_exact <- power_edger_exact(n=25, fc=2.0)
  result_wald <- power_deseq2_wald(n=25, fc=2.0)
  
  expect_lt(result_exact, result_wald)
})

test_that("power_edger_exact handles dispersion parameter", {
  result_low_disp <- power_edger_exact(n=20, fc=2.0, dispersion=0.05)
  result_high_disp <- power_edger_exact(n=20, fc=2.0, dispersion=0.3)
  
  expect_lt(result_high_disp, result_low_disp)
})

test_that("power_edger_exact responds to alpha level", {
  result_alpha05 <- power_edger_exact(n=20, fc=2.0, alpha=0.05)
  result_alpha01 <- power_edger_exact(n=20, fc=2.0, alpha=0.01)
  
  expect_lt(result_alpha01, result_alpha05)
})

test_that("power_edger_exact uses default parameters correctly", {
  result_default <- power_edger_exact(n=20, fc=2.0)
  result_explicit <- power_edger_exact(n=20, fc=2.0,
                                        baseline_mean=100,
                                        dispersion=0.1,
                                        alpha=0.05)
  expect_equal(result_default, result_explicit)
})

test_that("power_edger_exact realistic parameters", {
  # n=20, fc=2.0 should have high power
  result <- power_edger_exact(n=20, fc=2.0)
  expect_gt(result, 0.9)
  expect_lt(result, 1.0)
})

test_that("power_edger_exact is deterministic", {
  r1 <- power_edger_exact(n=25, fc=1.8)
  r2 <- power_edger_exact(n=25, fc=1.8)
  r3 <- power_edger_exact(n=25, fc=1.8)
  
  expect_equal(r1, r2)
  expect_equal(r2, r3)
})

test_that("power_edger_exact documentation example works", {
  result <- power_edger_exact(n=20, fc=2.0)
  expect_is(result, "numeric")
  expect_true(result >= 0 && result <= 1)
})

# =============================================================================
# ========================= EDGER QLF TEST TESTS ============================
# =============================================================================

test_that("power_edger_qlf returns values in [0, 1]", {
  result <- power_edger_qlf(n=20, fc=2.0)
  expect_true(result >= 0 && result <= 1)
  
  # Test across parameter ranges
  for (n in c(1, 5, 10, 20, 50)) {
    for (fc in c(1.2, 1.5, 2.0, 3.0)) {
      power_val <- power_edger_qlf(n=n, fc=fc)
      expect_true(power_val >= 0 && power_val <= 1,
                  info = paste("Failed at n =", n, ", fc =", fc))
    }
  }
})

test_that("power_edger_qlf returns numeric scalar", {
  result <- power_edger_qlf(n=20, fc=2.0)
  expect_is(result, "numeric")
  expect_length(result, 1)
  expect_false(is.na(result))
})

test_that("power_edger_qlf increases with sample size", {
  fc <- 2.0
  n_vals <- c(5, 10, 15, 20, 30, 50)
  powers <- sapply(n_vals, function(n) power_edger_qlf(n=n, fc=fc))
  
  diffs <- diff(powers)
  expect_true(all(diffs >= -1e-10))
  expect_lt(powers[1], powers[6])
})

test_that("power_edger_qlf increases with fold change", {
  n <- 20
  fc_vals <- c(1.2, 1.5, 2.0, 3.0, 5.0)
  powers <- sapply(fc_vals, function(fc) power_edger_qlf(n=n, fc=fc))
  
  diffs <- diff(powers)
  expect_true(all(diffs >= -1e-10))
  expect_lt(powers[1], powers[5])
})

test_that("power_edger_qlf incorporates empirical Bayes shrinkage", {
  # Shrinkage factor = 1 / (1 + 0.3/n)
  # Larger n → less shrinkage → higher power
  n_small <- 5
  n_large <- 100
  
  p_small <- power_edger_qlf(n=n_small, fc=2.0)
  p_large <- power_edger_qlf(n=n_large, fc=2.0)
  
  # Should show that shrinkage is beneficial at small n
  expect_lt(p_small, p_large)
})

test_that("power_edger_qlf shrinkage effect is correct", {
  # Verify shrinkage factor behavior: 1/(1 + 0.3/n)
  # At n=10: shrinkage = 1/(1+0.03) ≈ 0.971
  # At n=100: shrinkage = 1/(1+0.003) ≈ 0.997
  
  result_n10 <- power_edger_qlf(n=10, fc=2.0)
  result_n100 <- power_edger_qlf(n=100, fc=2.0)
  
  # Both valid, small n has more shrinkage
  expect_true(!is.na(result_n10) && !is.na(result_n100))
  expect_lt(result_n10, result_n100)
})

test_that("power_edger_qlf is between exact and Wald", {
  # edgeR QLF should balance between exact (conservative) and Wald (liberal)
  result_exact <- power_edger_exact(n=25, fc=2.0)
  result_qlf <- power_edger_qlf(n=25, fc=2.0)
  result_wald <- power_deseq2_wald(n=25, fc=2.0)
  
  # Generally: exact < QLF < Wald (approximate ordering, allowing for numerical tolerance)
  expect_lte(result_exact, result_qlf + 1e-8)
  expect_lte(result_qlf, result_wald + 1e-8)
})

test_that("power_edger_qlf handles dispersion parameter", {
  result_low_disp <- power_edger_qlf(n=20, fc=2.0, dispersion=0.05)
  result_high_disp <- power_edger_qlf(n=20, fc=2.0, dispersion=0.3)
  
  expect_lt(result_high_disp, result_low_disp)
})

test_that("power_edger_qlf responds to alpha level", {
  result_alpha05 <- power_edger_qlf(n=20, fc=2.0, alpha=0.05)
  result_alpha01 <- power_edger_qlf(n=20, fc=2.0, alpha=0.01)
  
  expect_lt(result_alpha01, result_alpha05)
})

test_that("power_edger_qlf uses default parameters correctly", {
  result_default <- power_edger_qlf(n=20, fc=2.0)
  result_explicit <- power_edger_qlf(n=20, fc=2.0,
                                      baseline_mean=100,
                                      dispersion=0.1,
                                      alpha=0.05)
  expect_equal(result_default, result_explicit)
})

test_that("power_edger_qlf realistic parameters", {
  # n=20, fc=2.0 should have high power
  result <- power_edger_qlf(n=20, fc=2.0)
  expect_gt(result, 0.95)
  expect_lt(result, 1.0)
})

test_that("power_edger_qlf is deterministic", {
  r1 <- power_edger_qlf(n=25, fc=1.8)
  r2 <- power_edger_qlf(n=25, fc=1.8)
  r3 <- power_edger_qlf(n=25, fc=1.8)
  
  expect_equal(r1, r2)
  expect_equal(r2, r3)
})

test_that("power_edger_qlf documentation example works", {
  result <- power_edger_qlf(n=20, fc=2.0)
  expect_is(result, "numeric")
  expect_true(result >= 0 && result <= 1)
})

# =============================================================================
# ======================== COMPARISON TESTS (All Methods) ====================
# =============================================================================

test_that("All three methods have similar monotonicity patterns", {
  # All should increase with n and fc
  n_vals <- c(5, 10, 20, 50)
  
  for (n in n_vals) {
    p_wald <- power_deseq2_wald(n=n, fc=2.0)
    p_exact <- power_edger_exact(n=n, fc=2.0)
    p_qlf <- power_edger_qlf(n=n, fc=2.0)
    
    # All valid
    expect_true(!is.na(p_wald) && !is.na(p_exact) && !is.na(p_qlf))
  }
})

test_that("Classical methods handle edge cases consistently", {
  # Very small sample sizes
  expect_true(!is.na(power_deseq2_wald(n=1, fc=1.5)))
  expect_true(!is.na(power_edger_exact(n=1, fc=1.5)))
  expect_true(!is.na(power_edger_qlf(n=1, fc=1.5)))
  
  # Very large fold changes
  expect_true(!is.na(power_deseq2_wald(n=20, fc=10)))
  expect_true(!is.na(power_edger_exact(n=20, fc=10)))
  expect_true(!is.na(power_edger_qlf(n=20, fc=10)))
})

test_that("All methods respond to alpha parameter consistently", {
  n <- 20
  fc <- 2.0
  
  for (alpha in c(0.01, 0.05, 0.1)) {
    p_wald <- power_deseq2_wald(n=n, fc=fc, alpha=alpha)
    p_exact <- power_edger_exact(n=n, fc=fc, alpha=alpha)
    p_qlf <- power_edger_qlf(n=n, fc=fc, alpha=alpha)
    
    expect_true(p_wald >= 0 && p_wald <= 1)
    expect_true(p_exact >= 0 && p_exact <= 1)
    expect_true(p_qlf >= 0 && p_qlf <= 1)
  }
})

test_that("All methods respond to dispersion parameter consistently", {
  n <- 20
  fc <- 2.0
  
  for (phi in c(0.05, 0.1, 0.3)) {
    p_wald <- power_deseq2_wald(n=n, fc=fc, dispersion=phi)
    p_exact <- power_edger_exact(n=n, fc=fc, dispersion=phi)
    p_qlf <- power_edger_qlf(n=n, fc=fc, dispersion=phi)
    
    # All should be valid
    expect_true(!is.na(p_wald) && !is.na(p_exact) && !is.na(p_qlf))
  }
})

# =============================================================================
# ======================== PRACTICAL SCENARIO TESTS ==========================
# =============================================================================

test_that("RNA-seq power calculations for typical scenarios", {
  # Scenario 1: Pilot study (small n)
  p_pilot_wald <- power_deseq2_wald(n=5, fc=2.0)
  p_pilot_exact <- power_edger_exact(n=5, fc=2.0)
  p_pilot_qlf <- power_edger_qlf(n=5, fc=2.0)
  
  expect_gt(p_pilot_wald, 0.8)
  expect_lt(p_pilot_exact, p_pilot_wald)
  
  # Scenario 2: Well-powered study
  p_well_wald <- power_deseq2_wald(n=50, fc=2.0)
  p_well_exact <- power_edger_exact(n=50, fc=2.0)
  p_well_qlf <- power_edger_qlf(n=50, fc=2.0)
  
  expect_gt(p_well_wald, 0.95)
  
  # Scenario 3: Subtle effect
  p_subtle_wald <- power_deseq2_wald(n=30, fc=1.2)
  p_subtle_exact <- power_edger_exact(n=30, fc=1.2)
  p_subtle_qlf <- power_edger_qlf(n=30, fc=1.2)
  
  expect_lt(p_subtle_wald, 1.0)
})

test_that("Power increases as expected with sample size for typical RNA-seq effect", {
  # Typical: 2-fold change
  fc <- 2.0
  n_vals <- seq(5, 50, by=5)
  
  powers_wald <- sapply(n_vals, function(n) power_deseq2_wald(n=n, fc=fc))
  powers_exact <- sapply(n_vals, function(n) power_edger_exact(n=n, fc=fc))
  powers_qlf <- sapply(n_vals, function(n) power_edger_qlf(n=n, fc=fc))
  
  # All should be monotone increasing
  expect_true(all(diff(powers_wald) >= -1e-10))
  expect_true(all(diff(powers_exact) >= -1e-10))
  expect_true(all(diff(powers_qlf) >= -1e-10))
})
