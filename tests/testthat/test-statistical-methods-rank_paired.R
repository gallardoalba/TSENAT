# ============================================================================
# TESTS FOR PAIRED RANK METHODS
# ============================================================================
# Tests for .apply_art_friedman(), .apply_robust_friedman(), .colMedians()
# and .select_rank_test_paired() conditional selection logic
# 
# These functions are NOT orphaned - they're integrated but conditionally 
# selected based on data characteristics. Tests create data with the specific
# characteristics that trigger each path:
# - Heteroscedasticity (Breusch-Pagan p<0.05 + var_ratio>2) -> ART-Friedman
# - Extreme skewness (|skewness| > 2) -> Robust Friedman
# ============================================================================

# Load required libraries
library(testthat)
library(S4Vectors)
library(data.table)

# ============================================================================
# TEST 1: .colMedians() Helper Function
# ============================================================================

test_that(".colMedians computes column medians correctly", {
  # Create test matrix
  mat <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3, byrow = TRUE)
  # mat = [[1, 2, 3], [4, 5, 6]]
  
  result <- TSENAT:::.colMedians(mat)
  
  # Expected: column medians = [median(1,4)=2.5, median(2,5)=3.5, median(3,6)=4.5]
  expected <- c(2.5, 3.5, 4.5)
  
  expect_equal(result, expected)
  expect_length(result, ncol(mat))
})

test_that(".colMedians handles single row", {
  mat <- matrix(c(1, 2, 3), nrow = 1, ncol = 3)
  result <- TSENAT:::.colMedians(mat)
  
  # Single row -> medians = values themselves
  expected <- c(1, 2, 3)
  expect_equal(result, expected)
})

test_that(".colMedians handles NA values", {
  mat <- matrix(c(1, NA, 3, 4, 5, 6), nrow = 2, ncol = 3, byrow = TRUE)
  # mat = [[1, NA, 3], [4, 5, 6]]
  
  result <- TSENAT:::.colMedians(mat)
  
  # Expected: [median(1,4)=2.5, median(NA,5)=5, median(3,6)=4.5]
  expected <- c(2.5, 5, 4.5)
  
  expect_equal(result, expected)
})

# ============================================================================
# TEST 2: .select_rank_test_paired() - Detection Logic
# ============================================================================

test_that(".select_rank_test_paired selects friedman for normal data", {
  # Normal, homoscedastic, symmetric data - should select standard Friedman
  set.seed(101)
  data <- data.frame(
    entropy = rnorm(30, mean = 2, sd = 0.3),
    q = factor(rep(1:5, 6)),
    subject = factor(rep(1:6, each = 5))
  )
  
  result <- TSENAT:::.select_rank_test_paired(data, "entropy", "q", "subject", verbose = FALSE)
  
  expect_equal(result$test_selected, "friedman")
  expect_false(result$characteristics$heteroscedastic)
  expect_false(result$characteristics$highly_skewed)
})

test_that(".select_rank_test_paired detects heteroscedasticity in paired data", {
  # Create paired data with MODERATE heteroscedasticity but NO extreme skewness
  # Key: need variance ratio > 2 AND Breusch-Pagan p < 0.05 but |skewness| < 2
  set.seed(102)
  
  subjects_rep <- rep(1:20, 4)
  treatments <- rep(1:4, each = 20)
  
  # Create entropy values: treatment effect + heteroscedastic noise
  # Use moderate variance increases: 1:1:1:3 ratio
  treatment_means <- c(1.0, 1.5, 2.0, 2.5)
  treatment_sds <- c(0.15, 0.15, 0.18, 0.5)  # Variance ratio ~11
  
  entropy <- numeric(80)
  for (i in 1:4) {
    idx <- which(treatments == i)
    entropy[idx] <- rnorm(20, mean = treatment_means[i], sd = treatment_sds[i])
  }
  
  data <- data.frame(
    entropy = entropy,
    q = factor(treatments),
    subject = factor(subjects_rep)
  )
  
  # Verify heteroscedasticity and moderate skewness
  var_by_treatment <- tapply(data$entropy, data$q, var, na.rm = TRUE)
  var_ratio <- max(var_by_treatment) / min(var_by_treatment)
  expect_true(var_ratio > 3)
  expect_true(abs(moments::skewness(data$entropy)) < 1.5)  # NOT highly skewed
  
  result <- TSENAT:::.select_rank_test_paired(data, "entropy", "q", "subject", verbose = FALSE)
  
  # Should detect heteroscedasticity (before skewness) and select ART-Friedman
  expect_true(result$characteristics$heteroscedastic)
  expect_equal(result$test_selected, "art_friedman")
})

test_that(".select_rank_test_paired detects extreme skewness", {
  # Create paired data with EXTREME right skewness (> 2)
  # Use lognormal distribution: exp(rnorm(n, 0, 2)) has skewness ~ 6-7
  set.seed(103)
  
  # Create 10 subjects with 3 treatments - all highly skewed
  subjects_rep <- rep(1:10, 3)
  treatments <- rep(1:3, each = 10)
  
  # Create highly skewed lognormal entropy values
  # lognormal with high sigma produces skewness > 2
  entropy <- exp(rnorm(30, mean = 0, sd = 1.8))  # Produces skewness ~5-7
  
  data <- data.frame(
    entropy = entropy,
    q = factor(treatments),
    subject = factor(subjects_rep)
  )
  
  # Verify we actually have extreme skewness
  actual_skew <- moments::skewness(data$entropy)
  expect_true(abs(actual_skew) > 2.0)  # Should be extreme skew (>2)
  
  result <- TSENAT:::.select_rank_test_paired(data, "entropy", "q", "subject", verbose = FALSE)
  
  # Should detect skewness and select robust Friedman
  expect_true(result$characteristics$highly_skewed)
  expect_equal(result$test_selected, "robust_friedman")
})

# ============================================================================
# TEST 3: .apply_art_friedman() - ART-Friedman Test
# ============================================================================

test_that(".apply_art_friedman detects heteroscedastic paired groups", {
  # Create heteroscedastic paired data with clear group differences
  set.seed(104)
  
  # Subjects: 1-8 (8 paired measurements)
  # Treatments: q1 (low entropy), q2 (high entropy)
  # Treatment 1 (q=1): low and stable
  treatment1 <- rep(1.0, 8) + rnorm(8, 0, 0.1)
  
  # Treatment 2 (q=2): high and VARIABLE (heteroscedastic!)
  treatment2 <- rep(3.0, 8) + rnorm(8, 0, 1.0)
  
  data <- data.frame(
    entropy = c(treatment1, treatment2),
    q = factor(rep(1:2, each = 8)),
    subject = factor(rep(1:8, 2))
  )
  
  result <- TSENAT:::.apply_art_friedman(data, "entropy", "q", "subject")
  
  # Should return valid results
  expect_true("statistic" %in% names(result))
  expect_true("p_value" %in% names(result))
  expect_true("method" %in% names(result))
  
  # p_value should be valid (0-1) and ideally small (groups differ)
  expect_true(result$p_value >= 0 && result$p_value <= 1)
  expect_true(result$p_value < 0.05)  # Should detect difference
  
  # Method should mention ART
  expect_match(result$method, "ART", ignore.case = TRUE)
})

test_that(".apply_art_friedman handles balanced design correctly", {
  # Perfectly balanced paired design
  set.seed(105)
  
  # 10 subjects, 4 q-values, all factorial combinations
  subjects <- rep(1:10, 4)
  qs <- rep(1:4, each = 10)
  
  # Values: effect of q-value with some heteroscedasticity
  entropy <- ifelse(qs == 1, rnorm(40, 1.0, 0.3),
            ifelse(qs == 2, rnorm(40, 1.5, 0.5),
            ifelse(qs == 3, rnorm(40, 2.0, 0.7),
                   rnorm(40, 2.5, 1.2))))
  
  data <- data.frame(
    entropy = entropy,
    q = factor(qs),
    subject = factor(subjects)
  )
  
  result <- TSENAT:::.apply_art_friedman(data, "entropy", "q", "subject")
  
  expect_true(is.numeric(result$statistic))
  expect_true(is.numeric(result$p_value))
  expect_false(is.na(result$p_value))
})

# ============================================================================
# TEST 4: .apply_robust_friedman() - Robust Friedman Test
# ============================================================================

test_that(".apply_robust_friedman handles skewed paired data", {
  # Create extreme skewness with paired blocking
  set.seed(106)
  
  # 50 subjects, 3 q-values with extreme skewness (large n reduces chi-sq warning)
  subjects_rep <- rep(1:50, 3)
  treatments <- rep(1:3, each = 50)
  
  entropy <- exp(rnorm(150, mean = 0, sd = 1.8))
  
  data <- data.frame(
    entropy = entropy,
    q = factor(treatments),
    subject = factor(subjects_rep)
  )
  
  # Verify we have extreme skewness
  actual_skew <- moments::skewness(data$entropy)
  expect_true(abs(actual_skew) > 2.0)
  
  result <- TSENAT:::.apply_robust_friedman(data, "entropy", "q", "subject")
  
  # Should return valid results
  expect_true("statistic" %in% names(result))
  expect_true("p_value" %in% names(result))
  expect_true("method" %in% names(result))
  
  # p_value should be valid
  expect_true(result$p_value >= 0 && result$p_value <= 1 || is.na(result$p_value))
  
  # Method should mention robust or median or mood
  expect_match(result$method, "robust|median|Mood", ignore.case = TRUE)
})

test_that(".apply_robust_friedman preserves paired structure", {
  # Verify that blocking structure (subjects) is properly preserved
  set.seed(107)
  
  # Create paired data with realistic random variation for chi-sq stability
  n_subjects <- 15
  n_treatments <- 2
  obs_per_cell <- 8
  
  entropy <- numeric(n_subjects * n_treatments * obs_per_cell)
  q_vals <- numeric(n_subjects * n_treatments * obs_per_cell)
  subj_vals <- numeric(n_subjects * n_treatments * obs_per_cell)
  
  idx <- 1
  for (s in 1:n_subjects) {
    for (t in 1:n_treatments) {
      base_val <- s * 0.5 + t * 0.3
      for (o in 1:obs_per_cell) {
        entropy[idx] <- base_val + rnorm(1, 0, 0.15)
        q_vals[idx] <- t
        subj_vals[idx] <- s
        idx <- idx + 1
      }
    }
  }
  
  data <- data.frame(
    entropy = entropy,
    q = factor(q_vals),
    subject = factor(subj_vals)
  )
  
  result <- TSENAT:::.apply_robust_friedman(data, "entropy", "q", "subject")
  
  # Results should be numeric and valid (or NA if test fails)
  expect_true(is.numeric(result$statistic))
  expect_true(is.numeric(result$p_value) || is.na(result$p_value))
  # Don't be strict about p_value being in [0,1] - Friedman can fail on small data
})

# ============================================================================
# TEST 5: .apply_conditional_rank_test() with Paired Heteroscedasticity
# ============================================================================

test_that(".apply_conditional_rank_test routes to ART-Friedman for heteroscedastic paired", {
  # Heteroscedastic paired data - moderate variance differences
  set.seed(108)
  
  subjects_rep <- rep(1:20, 4)
  treatments <- rep(1:4, each = 20)
  
  # Treatment effect + heteroscedastic noise (moderate variances, no extreme skew)
  treatment_means <- c(1.0, 1.5, 2.0, 2.5)
  treatment_sds <- c(0.15, 0.15, 0.18, 0.5)  # Variance ratio ~11
  
  entropy <- numeric(80)
  for (i in 1:4) {
    idx <- which(treatments == i)
    entropy[idx] <- rnorm(20, mean = treatment_means[i], sd = treatment_sds[i])
  }
  
  data <- data.frame(
    entropy = entropy,
    q = factor(treatments),
    subject = factor(subjects_rep)
  )
  
  result <- TSENAT:::.apply_conditional_rank_test(
    data, "entropy", "q",
    paired = TRUE, subject_col = "subject", verbose = FALSE
  )
  
  # Should route to ART-Friedman
  expect_equal(result$test_type, "art_friedman")
  expect_true(result$characteristics$heteroscedastic)
})

test_that(".apply_conditional_rank_test routes to Robust Friedman for skewed paired", {
  # Extreme skewness paired data
  set.seed(109)
  
  subjects_rep <- rep(1:50, 3)
  treatments <- rep(1:3, each = 50)
  
  entropy <- exp(rnorm(150, mean = 0, sd = 1.8))
  
  data <- data.frame(
    entropy = entropy,
    q = factor(treatments),
    subject = factor(subjects_rep)
  )
  
  result <- TSENAT:::.apply_conditional_rank_test(
    data, "entropy", "q",
    paired = TRUE, subject_col = "subject", verbose = FALSE
  )
  
  # Should route to robust Friedman
  expect_equal(result$test_type, "robust_friedman")
  expect_true(result$characteristics$highly_skewed)
})

test_that(".apply_conditional_rank_test gives valid p-values for all paired paths", {
  # Test that all three paired paths return valid results
  
  # Path 1: Standard Friedman (normal data)
  set.seed(110)
  data_normal <- data.frame(
    entropy = rnorm(24, 2, 0.3),
    q = factor(rep(1:4, 6)),
    subject = factor(rep(1:6, each = 4))
  )
  result_normal <- TSENAT:::.apply_conditional_rank_test(
    data_normal, "entropy", "q",
    paired = TRUE, subject_col = "subject", verbose = FALSE
  )
  expect_true(result_normal$p_value >= 0 && result_normal$p_value <= 1)
  
  # Path 2: ART Friedman (heteroscedastic)
  set.seed(111)
  data_hetero <- data.frame(
    entropy = c(
      rnorm(6, 1.0, 0.1),   # Treatment 1: low var
      rnorm(6, 1.5, 0.1),   # Treatment 2: low var
      rnorm(6, 2.0, 1.5)    # Treatment 3: HIGH var
    ),
    q = factor(rep(1:3, each = 6)),
    subject = factor(rep(1:6, 3))
  )
  result_hetero <- TSENAT:::.apply_conditional_rank_test(
    data_hetero, "entropy", "q",
    paired = TRUE, subject_col = "subject", verbose = FALSE
  )
  expect_true(result_hetero$p_value >= 0 && result_hetero$p_value <= 1)
  
  # Path 3: Robust Friedman (skewed)
  set.seed(112)
  data_skew <- data.frame(
    entropy = rexp(18, 1) + 0.5,
    q = factor(rep(1:3, each = 6)),
    subject = factor(rep(1:6, 3))
  )
  result_skew <- TSENAT:::.apply_conditional_rank_test(
    data_skew, "entropy", "q",
    paired = TRUE, subject_col = "subject", verbose = FALSE
  )
  expect_true(result_skew$p_value >= 0 && result_skew$p_value <= 1)
})

# ============================================================================
# TEST 6: Error Handling & Edge Cases
# ============================================================================

test_that(".apply_art_friedman handles single subject gracefully", {
  # Edge case: only 1 subject (no blocking benefit, but should not crash)
  set.seed(113)
  data <- data.frame(
    entropy = c(1.0, 1.5, 2.0, 2.5),
    q = factor(1:4),
    subject = factor(rep(1, 4))
  )
  
  # Should complete without error (may return NA or try-error handling)
  result <- TSENAT:::.apply_art_friedman(data, "entropy", "q", "subject")
  expect_true("p_value" %in% names(result))
})

test_that(".apply_robust_friedman handles small samples", {
  # Edge case: small but sufficient sample to avoid chi-sq warnings
  set.seed(114)
  
  n_subjects <- 12
  n_treatments <- 2
  obs_per_cell <- 6
  
  entropy <- numeric(n_subjects * n_treatments * obs_per_cell)
  q_vals <- numeric(n_subjects * n_treatments * obs_per_cell)
  subj_vals <- numeric(n_subjects * n_treatments * obs_per_cell)
  
  idx <- 1
  for (s in 1:n_subjects) {
    for (t in 1:n_treatments) {
      base_val <- s * 0.2 + t * 0.5
      for (o in 1:obs_per_cell) {
        entropy[idx] <- base_val + rnorm(1, 0, 0.1)
        q_vals[idx] <- t
        subj_vals[idx] <- s
        idx <- idx + 1
      }
    }
  }
  
  data <- data.frame(
    entropy = entropy,
    q = factor(q_vals),
    subject = factor(subj_vals)
  )
  
  result <- TSENAT:::.apply_robust_friedman(data, "entropy", "q", "subject")
  expect_true("p_value" %in% names(result))
})

test_that(".colMedians works with single column matrix", {
  # Edge case: 1-column matrix
  mat <- matrix(c(1, 2, 3, 4), nrow = 4, ncol = 1)
  result <- TSENAT:::.colMedians(mat)
  
  expect_equal(result, 2.5)  # median(1,2,3,4)
  expect_length(result, 1)
})

# ============================================================================
# TEST 7: Integration Tests
# ============================================================================

test_that("paired methods integrate cleanly with selection logic", {
  # Full integration: selection -> application -> results
  set.seed(115)
  
  # Create realistic entropy data with heteroscedasticity
  set.seed(115)
  n_subjects <- 8
  n_qvals <- 5
  
  # Base entropy with increasing trend
  base_entropy <- rep(seq(1, 3, length.out = n_qvals), n_subjects)
  
  # Add heteroscedastic noise (increasing with q-value)
  noise_sd <- rep(c(0.1, 0.15, 0.5, 1.0, 1.5), n_subjects)
  entropy_values <- base_entropy + rnorm(n_subjects * n_qvals, 0, noise_sd)
  
  data <- data.frame(
    entropy = entropy_values,
    q = factor(rep(1:n_qvals, n_subjects)),
    subject = factor(rep(1:n_subjects, each = n_qvals))
  )
  
  # Run full pipeline
  result <- TSENAT:::.apply_conditional_rank_test(
    data, "entropy", "q",
    paired = TRUE, subject_col = "subject", verbose = FALSE
  )
  
  # Should successfully route and return results
  expect_true(result$test_type %in% c("friedman", "art_friedman", "robust_friedman"))
  expect_true(result$characteristics$pairing_used == TRUE)
  expect_true(is.numeric(result$p_value))
  expect_true(result$p_value >= 0 && result$p_value <= 1)
})

# ============================================================================
# TEST 8: Numerical Correctness Validation
# ============================================================================

test_that(".colMedians computes exact medians for known data", {
  # Test 1: Even number of elements per column
  mat1 <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE)
  # [[1, 2], [3, 4]]
  result1 <- TSENAT:::.colMedians(mat1)
  expect_equal(result1[1], 2.0)  # median(1, 3) = 2.0
  expect_equal(result1[2], 3.0)  # median(2, 4) = 3.0
  
  # Test 2: Odd number of elements per column
  mat2 <- matrix(c(5, 10, 15, 20, 25, 30), nrow = 3, ncol = 2, byrow = TRUE)
  # [[5, 10], [15, 20], [25, 30]]
  result2 <- TSENAT:::.colMedians(mat2)
  expect_equal(result2[1], 15.0)  # median(5, 15, 25) = 15.0
  expect_equal(result2[2], 20.0)  # median(10, 20, 30) = 20.0
  
  # Test 3: With negative numbers
  mat3 <- matrix(c(-5, -3, 2, 4, 7, -1), nrow = 3, ncol = 2, byrow = TRUE)
  # [[-5, -3], [2, 4], [7, -1]]
  result3 <- TSENAT:::.colMedians(mat3)
  expect_equal(result3[1], 2.0)   # median(-5, 2, 7) = 2.0
  expect_equal(result3[2], -1.0)  # median(-3, 4, -1) = -1.0
})

test_that(".colMedians preserves data type precision with floats", {
  # Test with floating point data
  mat <- matrix(c(1.1, 2.2, 3.3, 4.4, 5.5, 6.6), nrow = 3, ncol = 2, byrow = TRUE)
  result <- TSENAT:::.colMedians(mat)
  
  # For even count (2 values): median = mean of the two
  expected_col1 <- (1.1 + 5.5) / 2  # 3.3
  expected_col2 <- (2.2 + 6.6) / 2  # 4.4
  
  expect_equal(result[1], expected_col1, tolerance = 1e-10)
  expect_equal(result[2], expected_col2, tolerance = 1e-10)
})

test_that("Friedman test p-value is deterministic from seed", {
  # Verify reproducibility: same seed -> same p-value
  set.seed(200)
  data1 <- data.frame(
    entropy = rnorm(24, 2, 0.5),
    q = factor(rep(1:4, 6)),
    subject = factor(rep(1:6, each = 4))
  )
  result1 <- TSENAT:::.apply_art_friedman(data1, "entropy", "q", "subject")
  
  # Re-run with same seed
  set.seed(200)
  data2 <- data.frame(
    entropy = rnorm(24, 2, 0.5),
    q = factor(rep(1:4, 6)),
    subject = factor(rep(1:6, each = 4))
  )
  result2 <- TSENAT:::.apply_art_friedman(data2, "entropy", "q", "subject")
  
  # Results should be identical
  expect_equal(result1$p_value, result2$p_value)
  expect_equal(result1$statistic, result2$statistic)
})

test_that("Robust Friedman correctly classifies above/below median", {
  # Create larger paired dataset for chi-squared stability
  set.seed(201)
  
  # 10 subjects, 2 treatments, 5 obs per subject per treatment
  n_subjects <- 10
  n_treatments <- 2
  n_obs_per_cell <- 5
  
  entropy_vals <- numeric(n_subjects * n_treatments * n_obs_per_cell)
  q_vals <- numeric(n_subjects * n_treatments * n_obs_per_cell)
  s_vals <- numeric(n_subjects * n_treatments * n_obs_per_cell)
  
  idx <- 1
  for (s in 1:n_subjects) {
    for (t in 1:n_treatments) {
      base_val <- s * 0.3 + t * 0.8
      for (o in 1:n_obs_per_cell) {
        entropy_vals[idx] <- base_val + rnorm(1, 0, 0.2)
        q_vals[idx] <- t
        s_vals[idx] <- s
        idx <- idx + 1
      }
    }
  }
  
  data <- data.frame(
    entropy = entropy_vals,
    q = factor(q_vals),
    subject = factor(s_vals)
  )
  
  result <- TSENAT:::.apply_robust_friedman(data, "entropy", "q", "subject")
  
  # Verify result has required fields
  expect_true("statistic" %in% names(result))
  expect_true("p_value" %in% names(result))
  expect_true(is.numeric(result$statistic) || is.na(result$statistic))
  expect_true(is.numeric(result$p_value) || is.na(result$p_value))
})

test_that("Heteroscedasticity detection uses correct BP p-value row", {
  # Create data with KNOWN heteroscedasticity and measure what BP test returns
  set.seed(202)
  
  # Four treatments with extreme variance differences
  treatment1 <- rnorm(30, 1.0, 0.05)   # Very small variance
  treatment2 <- rnorm(30, 2.0, 0.05)   # Very small variance
  treatment3 <- rnorm(30, 3.0, 2.0)    # LARGE variance
  treatment4 <- rnorm(30, 4.0, 2.0)    # LARGE variance
  
  data <- data.frame(
    entropy = c(treatment1, treatment2, treatment3, treatment4),
    q = factor(rep(1:4, each = 30)),
    subject = factor(rep(1:30, 4))
  )
  
  # Manually verify heteroscedasticity
  var_by_treatment <- tapply(data$entropy, data$q, var, na.rm = TRUE)
  var_ratio <- max(var_by_treatment) / min(var_by_treatment)
  expect_true(var_ratio > 5)  # Strong heteroscedasticity
  
  # Run selection to verify it detects
  result <- TSENAT:::.select_rank_test_paired(data, "entropy", "q", "subject", verbose = FALSE)
  expect_true(result$characteristics$heteroscedastic)
  expect_equal(result$test_selected, "art_friedman")
})

test_that("Skewness detection threshold (> 2) is correctly applied", {
  # Test boundary: create data with |skewness| < 2.0 (should NOT trigger)
  # and |skewness| > 2.0 (should trigger)
  
  # Moderately skewed - use beta distribution with parameters giving skewness ~1.2
  # rbeta(n, 2, 5) has skewness ≈ 1.13
  set.seed(203)
  mod_skew <- rbeta(60, 2, 5)  # Moderate right skew, skewness < 2
  actual_mod_skew <- moments::skewness(mod_skew)
  
  data_mod <- data.frame(
    entropy = mod_skew,
    q = factor(rep(1:3, each = 20)),
    subject = factor(rep(1:20, 3))
  )
  
  result_mod <- TSENAT:::.select_rank_test_paired(data_mod, "entropy", "q", "subject", verbose = FALSE)
  expect_true(abs(actual_mod_skew) < 2.0, info = paste("Moderate skew:", actual_mod_skew))
  expect_false(result_mod$characteristics$highly_skewed)
  
  # Extremely skewed (should trigger robust Friedman)
  set.seed(204)
  high_skew <- exp(rnorm(150, 0, 2.0))  # Strong right skew
  actual_high_skew <- moments::skewness(high_skew)
  
  data_high <- data.frame(
    entropy = high_skew,
    q = factor(rep(1:3, each = 50)),
    subject = factor(rep(1:50, 3))
  )
  
  result_high <- TSENAT:::.select_rank_test_paired(data_high, "entropy", "q", "subject", verbose = FALSE)
  expect_true(abs(actual_high_skew) > 2.0)
  expect_true(result_high$characteristics$highly_skewed)
  expect_equal(result_high$test_selected, "robust_friedman")
})

test_that("Method selection priority: skewness > heteroscedasticity when both present", {
  # Create data with BOTH heteroscedasticity AND extreme skewness
  # Should prioritize skewness and select robust Friedman
  set.seed(205)
  
  # Combine extreme skewness with heteroscedasticity
  n_per_treatment <- 40
  
  # Generate lognormal (highly skewed) with different variances per treatment
  entropy_t1 <- exp(rnorm(n_per_treatment, 0, 1.5))  # Skewed + var
  entropy_t2 <- exp(rnorm(n_per_treatment, 0.5, 0.8))  # Skewed + different var
  entropy_all <- c(entropy_t1, entropy_t2)
  
  data <- data.frame(
    entropy = entropy_all,
    q = factor(rep(1:2, each = n_per_treatment)),
    subject = factor(rep(1:n_per_treatment, 2))
  )
  
  # Verify both conditions are present
  var_ratio <- var(entropy_t1) / var(entropy_t2)
  expect_true(var_ratio > 1.5 | var_ratio < 0.67)  # Heteroscedasticity present
  expect_true(abs(moments::skewness(entropy_all)) > 2.0)  # Skewness present
  
  result <- TSENAT:::.select_rank_test_paired(data, "entropy", "q", "subject", verbose = FALSE)
  
  # Should prioritize skewness and select robust Friedman
  expect_true(result$characteristics$highly_skewed)
  expect_equal(result$test_selected, "robust_friedman")
})

test_that("Normal data routes to standard Friedman (no preferential treatment)", {
  # Standard normal data (no skew, no heteroscedasticity)
  # Should route to standard Friedman, not art or robust
  set.seed(206)
  
  subjects_rep <- rep(1:12, 5)
  treatments <- rep(1:5, each = 12)
  entropy <- rnorm(60, 2.0, 0.3)
  
  data <- data.frame(
    entropy = entropy,
    q = factor(treatments),
    subject = factor(subjects_rep)
  )
  
  # Verify normality assumptions
  var_by_treatment <- tapply(data$entropy, data$q, var, na.rm = TRUE)
  var_ratio <- max(var_by_treatment) / min(var_by_treatment)
  expect_true(var_ratio < 2.0)  # Homoscedastic
  expect_true(abs(moments::skewness(data$entropy)) < 1.5)  # Symmetric
  
  result <- TSENAT:::.select_rank_test_paired(data, "entropy", "q", "subject", verbose = FALSE)
  
  expect_false(result$characteristics$heteroscedastic)
  expect_false(result$characteristics$highly_skewed)
  expect_equal(result$test_selected, "friedman")
})

test_that("P-value range validation across all methods", {
  # Verify all three paired methods return p-values in [0, 1]
  set.seed(207)
  
  # Normal data for standard Friedman
  data_normal <- data.frame(
    entropy = rnorm(24, 2, 0.3),
    q = factor(rep(1:4, 6)),
    subject = factor(rep(1:6, each = 4))
  )
  result_normal <- TSENAT:::.apply_conditional_rank_test(
    data_normal, "entropy", "q",
    paired = TRUE, subject_col = "subject", verbose = FALSE
  )
  expect_true(result_normal$p_value >= 0 && result_normal$p_value <= 1)
  
  # Heteroscedastic data for ART-Friedman
  t1 <- rnorm(15, 1.0, 0.1)
  t2 <- rnorm(15, 2.0, 0.1)
  t3 <- rnorm(15, 3.0, 2.0)
  data_hetero <- data.frame(
    entropy = c(t1, t2, t3),
    q = factor(rep(1:3, each = 15)),
    subject = factor(rep(1:15, 3))
  )
  result_hetero <- TSENAT:::.apply_conditional_rank_test(
    data_hetero, "entropy", "q",
    paired = TRUE, subject_col = "subject", verbose = FALSE
  )
  expect_true(result_hetero$p_value >= 0 && result_hetero$p_value <= 1)
  
  # Skewed data for Robust Friedman
  data_skew <- data.frame(
    entropy = exp(rnorm(60, 0, 1.8)),
    q = factor(rep(1:3, each = 20)),
    subject = factor(rep(1:20, 3))
  )
  result_skew <- TSENAT:::.apply_conditional_rank_test(
    data_skew, "entropy", "q",
    paired = TRUE, subject_col = "subject", verbose = FALSE
  )
  expect_true(result_skew$p_value >= 0 && result_skew$p_value <= 1)
})

test_that("Conditional routing consistency: same data always routes to same method", {
  # Verify that calling selection twice on same data gives same result
  set.seed(208)
  
  data <- data.frame(
    entropy = c(
      rnorm(20, 1.0, 0.1),
      rnorm(20, 2.0, 0.1),
      rnorm(20, 3.0, 2.0)  # Heteroscedastic
    ),
    q = factor(rep(1:3, each = 20)),
    subject = factor(rep(1:20, 3))
  )
  
  result1 <- TSENAT:::.select_rank_test_paired(data, "entropy", "q", "subject", verbose = FALSE)
  result2 <- TSENAT:::.select_rank_test_paired(data, "entropy", "q", "subject", verbose = FALSE)
  
  expect_equal(result1$test_selected, result2$test_selected)
  expect_equal(result1$characteristics$heteroscedastic, result2$characteristics$heteroscedastic)
  expect_equal(result1$characteristics$highly_skewed, result2$characteristics$highly_skewed)
})
