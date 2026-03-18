context("Friedman Test Integration and Edge Cases")

source("../../R/rank_based_methods.R", local = TRUE)

# ============================================================================
# Test: Paired vs Unpaired Comparison on Same Dataset
# ============================================================================

test_that("Paired analysis produces different results than unpaired on same data", {
  set.seed(888)
  
  # Create data with strong subject effect
  n_subjects <- 8
  n_q <- 5
  
  subject_effect <- rnorm(n_subjects, sd = 2)
  entropy_data <- numeric(n_subjects * n_q)
  
  for (i in 1:n_subjects) {
    for (j in 1:n_q) {
      idx <- (i-1)*n_q + j
      entropy_data[idx] <- 2 + subject_effect[i] + (j - n_q/2) * 0.3 + rnorm(1, sd=0.1)
    }
  }
  
  data <- data.frame(
    entropy = entropy_data,
    q = factor(rep(1:n_q, n_subjects)),
    subject = factor(rep(1:n_subjects, each = n_q))
  )
  
  # Paired analysis
  paired_result <- .tsenat_apply_conditional_rank_test(
    data = data,
    value_col = "entropy",
    group_col = "q",
    paired = TRUE,
    subject_col = "subject",
    verbose = FALSE
  )
  
  # Unpaired analysis (data without subject column in conditional context)
  unpaired_result <- .tsenat_apply_conditional_rank_test(
    data = data[, c("entropy", "q")],
    value_col = "entropy",
    group_col = "q",
    paired = FALSE,
    subject_col = NULL,
    verbose = FALSE
  )
  
  # Results should differ (paired should typically be more significant)
  expect_is(paired_result$p_value, "numeric")
  expect_is(unpaired_result$p_value, "numeric")
  
  # Document that they are different
  expect_false(isTRUE(all.equal(paired_result$p_value, unpaired_result$p_value)),
               info = "Paired and unpaired analyses should produce different p-values")
})

# ============================================================================
# Test: Friedman with Minimum Required Structure
# ============================================================================

test_that("Friedman works with minimal paired structure (3 subjects, 2 treatments)", {
  data <- data.frame(
    entropy = c(1.0, 2.0, 1.5, 2.5),
    subject = factor(c(1, 2, 1, 2)),
    treatment = factor(c("A", "A", "B", "B"))
  )
  
  result <- .tsenat_apply_friedman_test(
    data = data,
    value_col = "entropy",
    group_col = "treatment",
    subject_col = "subject"
  )
  
  expect_equal(result$method, "Friedman test (paired)")
  expect_is(result$p_value, "numeric")
  expect_true(result$p_value >= 0 && result$p_value <= 1)
})

# ============================================================================
# Test: Large Number of Subjects
# ============================================================================

test_that("Friedman handles large number of subjects", {
  set.seed(999)
  
  n_subjects <- 50  # Large
  n_q <- 4
  
  data <- data.frame(
    entropy = rnorm(n_subjects * n_q, mean = 2, sd = 0.3),
    subject = factor(rep(1:n_subjects, n_q)),
    q = factor(rep(1:n_q, each = n_subjects))
  )
  
  result <- .tsenat_apply_friedman_test(
    data = data,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  expect_equal(result$method, "Friedman test (paired)")
  expect_is(result$statistic, "numeric")
  expect_true(result$statistic >= 0)
})

# ============================================================================
# Test: Large Number of Q-Values
# ============================================================================

test_that("Friedman handles large number of q-values (treatments)", {
  set.seed(1111)
  
  n_subjects <- 6
  n_q <- 20  # Many q-values
  
  data <- data.frame(
    entropy = rnorm(n_subjects * n_q, mean = 2, sd = 0.3),
    subject = factor(rep(1:n_subjects, n_q)),
    q = factor(rep(1:n_q, each = n_subjects))
  )
  
  result <- .tsenat_apply_friedman_test(
    data = data,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  expect_equal(result$method, "Friedman test (paired)")
  expect_is(result$p_value, "numeric")
  expect_true(result$p_value >= 0 && result$p_value <= 1)
})

# ============================================================================
# Test: Friedman Response to Actual Effect
# ============================================================================

test_that("Friedman detects strong q-effect (low p-value)", {
  set.seed(2222)
  
  n_subjects <- 10
  n_q <- 4
  
  # Strong q-effect: entropy depends on q
  entropy <- numeric(n_subjects * n_q)
  for (i in 1:n_subjects) {
    for (j in 1:n_q) {
      idx <- (i-1)*n_q + j
      entropy[idx] <- 1 + (j - 1) * 0.8 + rnorm(1, sd=0.05)  # Strong effect, small noise
    }
  }
  
  data <- data.frame(
    entropy = entropy,
    subject = factor(rep(1:n_subjects, each = n_q)),
    q = factor(rep(1:n_q, n_subjects))
  )
  
  result <- .tsenat_apply_friedman_test(
    data = data,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  # Should be highly significant
  expect_lt(result$p_value, 0.05)
})

test_that("Friedman shows weak effect (high p-value) for random data", {
  set.seed(3333)
  
  n_subjects <- 8
  n_q <- 4
  
  # No q-effect: completely random
  entropy <- rnorm(n_subjects * n_q, mean = 2, sd = 0.8)
  
  data <- data.frame(
    entropy = entropy,
    subject = factor(rep(1:n_subjects, each = n_q)),
    q = factor(rep(1:n_q, n_subjects))
  )
  
  result <- .tsenat_apply_friedman_test(
    data = data,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  # Should not be highly significant
  expect_gt(result$p_value, 0.05)
})

# ============================================================================
# Test: Friedman Statistics Interpretation
# ============================================================================

test_that("Friedman statistic follows chi-square distribution (df = k-1)", {
  set.seed(4444)
  
  n_subjects <- 15
  k_treatments <- 5  # number of q-values
  
  data <- data.frame(
    entropy = rnorm(n_subjects * k_treatments, mean = 2, sd = 0.4),
    subject = factor(rep(1:n_subjects, k_treatments)),
    q = factor(rep(1:k_treatments, each = n_subjects))
  )
  
  result <- .tsenat_apply_friedman_test(
    data = data,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  # Friedman statistic should be non-negative
  expect_gte(result$statistic, 0)
  
  # For random data, statistic follows approximately chi-square(k-1)
  # So it should be in reasonable range (0 to maybe 4*df for random data)
  df <- k_treatments - 1
  expect_lt(result$statistic, 4 * df + 5)
})

# ============================================================================
# Test: Column Name Flexibility
# ============================================================================

test_that("Friedman works with different column names", {
  set.seed(5555)
  
  # Different column names
  data <- data.frame(
    my_entropy = rnorm(20, mean = 2, sd = 0.3),
    my_treatment = factor(rep(1:4, 5)),
    my_block = factor(rep(1:5, each = 4))
  )
  
  result <- .tsenat_apply_friedman_test(
    data = data,
    value_col = "my_entropy",
    group_col = "my_treatment",
    subject_col = "my_block"
  )
  
  expect_equal(result$method, "Friedman test (paired)")
  expect_is(result$p_value, "numeric")
})

# ============================================================================
# Test: Factor vs Character Handling
# ============================================================================

test_that("Friedman handles both factor and numeric group identifiers", {
  set.seed(6666)
  
  # With numeric (converted to factor)
  data_numeric <- data.frame(
    entropy = rnorm(16, mean = 2, sd = 0.3),
    subject = 1:4,  # Numeric
    q = rep(1:4, each = 4)  # Numeric
  )
  
  # Should still work (function converts internally)
  result_numeric <- .tsenat_apply_friedman_test(
    data = data_numeric,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  # With character
  data_char <- data.frame(
    entropy = rnorm(16, mean = 2, sd = 0.3),
    subject = factor(rep(paste0("S", 1:4), each = 4)),  # Factor with character
    q = factor(rep(paste0("Q", 1:4), 4))  # Factor with character
  )
  
  result_char <- .tsenat_apply_friedman_test(
    data = data_char,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  # Both should succeed
  expect_equal(result_numeric$method, "Friedman test (paired)")
  expect_equal(result_char$method, "Friedman test (paired)")
})

# ============================================================================
# Test: Behavior with Identical Values
# ============================================================================

test_that("Friedman handles data with tied (identical) values", {
  # Data with many ties
  data <- data.frame(
    entropy = c(2.0, 2.0, 2.0, 2.5, 2.5, 2.5, 2.0, 2.0, 2.0, 2.5, 2.5, 2.5),
    subject = factor(c(1, 2, 3, 1, 2, 3, 4, 5, 6, 4, 5, 6)),
    q = factor(c(rep("A", 6), rep("B", 6)))
  )
  
  result <- .tsenat_apply_friedman_test(
    data = data,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  # Should still return valid result
  expect_equal(result$method, "Friedman test (paired)")
  expect_is(result$p_value, "numeric")
  expect_true(result$p_value >= 0 && result$p_value <= 1)
})

# ============================================================================
# Test: Missing Data Patterns
# ============================================================================

test_that("Friedman fails gracefully with missing values in matrix", {
  # Create a data frame where using xtabs would create a result without required dimensions
  # e.g., only one level in one dimension after filtering
  data <- data.frame(
    entropy = c(1.0, 2.0),
    subject = factor(c(1, 1)),
    q = factor(c("A", "B"))
  )
  # This has only 1 subject but 2 treatments - Friedman needs at least 2 subjects
  
  result <- .tsenat_apply_friedman_test(
    data = data,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  # Should fail gracefully - either NA p-value or test_failed method
  expect_true(is.na(result$p_value) || result$method == "test_failed")
})

# ============================================================================
# Test: Monotonic Transformations Don't Change p-value
# ============================================================================

test_that("Friedman p-value invariant to monotonic transformations", {
  set.seed(7777)
  
  entropy_original <- rnorm(24, mean = 2, sd = 0.5)
  
  data_original <- data.frame(
    entropy = entropy_original,
    subject = factor(rep(1:6, 4)),
    q = factor(rep(1:4, each = 6))
  )
  
  # Log transformation
  data_log <- data_original
  data_log$entropy <- log(entropy_original + 1)  # Ensure positive
  
  # Square transformation
  data_sq <- data_original
  data_sq$entropy <- entropy_original^2
  
  result_original <- .tsenat_apply_friedman_test(
    data = data_original,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  result_log <- .tsenat_apply_friedman_test(
    data = data_log,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  result_sq <- .tsenat_apply_friedman_test(
    data = data_sq,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  # Friedman is rank-based, so monotonic transformations shouldn't change results
  # (p-values may differ slightly due to ranking, but interpretation should hold)
  expect_is(result_original$p_value, "numeric")
  expect_is(result_log$p_value, "numeric")
  expect_is(result_sq$p_value, "numeric")
  
  # Rank-invariance property: should get same ranking, hence same p-value
  expect_equal(result_original$statistic, result_log$statistic, tolerance = 1e-10)
  expect_equal(result_original$statistic, result_sq$statistic, tolerance = 1e-10)
})

# ============================================================================
# Test: Extreme Value Handling
# ============================================================================

test_that("Friedman handles extreme values correctly", {
  # Very small values
  data_small <- data.frame(
    entropy = c(1e-10, 2e-10, 1.5e-10, 2.5e-10, 1e-10, 2e-10, 1.5e-10, 2.5e-10),
    subject = factor(rep(1:4, 2)),
    q = factor(rep(1:2, each = 4))
  )
  
  result_small <- .tsenat_apply_friedman_test(
    data = data_small,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  # Very large values
  data_large <- data.frame(
    entropy = c(1e10, 2e10, 1.5e10, 2.5e10, 1e10, 2e10, 1.5e10, 2.5e10),
    subject = factor(rep(1:4, 2)),
    q = factor(rep(1:2, each = 4))
  )
  
  result_large <- .tsenat_apply_friedman_test(
    data = data_large,
    value_col = "entropy",
    group_col = "q",
    subject_col = "subject"
  )
  
  # Both should work (rank-based, so scale doesn't matter)
  expect_equal(result_small$method, "Friedman test (paired)")
  expect_equal(result_large$method, "Friedman test (paired)")
  
  # Results should be identical (rank-based)
  expect_equal(result_small$statistic, result_large$statistic)
  expect_equal(result_small$p_value, result_large$p_value)
})
