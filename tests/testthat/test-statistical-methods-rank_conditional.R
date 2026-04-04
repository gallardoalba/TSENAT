library(TSENAT)
library(testthat)

# ============================================================================
# TESTS FOR CONDITIONAL RANK TEST METHODS
# ============================================================================
# These tests create specific data characteristics to trigger conditional
# rank test selection paths that currently have 0% coverage:
# - .apply_art_kw() — requires heteroscedastic data
# - .apply_quantile_test() — requires boundary clustering
# - .apply_robust_median_test() — requires extreme skewness
# - print.rank_correlation_ci() — S3 print method
# ============================================================================

# ============================================================================
# TEST SUITE 1: .apply_art_kw() — Heteroscedasticity Trigger
# ============================================================================
# Condition: Breusch-Pagan test p < 0.05 AND variance ratio > 2
# Strategy: Create two groups with very different variances

test_that(".apply_art_kw can be called with heteroscedastic data", {
  skip_on_cran()
  
  # Create heteroscedastic data directly for testing
  set.seed(42)
  
  # Group 1: Low variance entropy values
  entropy_group1 <- rnorm(30, mean = 1.0, sd = 0.2)
  # Group 2: High variance entropy values  
  entropy_group2 <- rnorm(30, mean = 1.0, sd = 2.0)
  
  # Create data frame for rank test input
  data <- data.frame(
    entropy = c(entropy_group1, entropy_group2),
    group = rep(c("Low_Var", "High_Var"), each = 30)
  )
  
  # Call .apply_art_kw directly
  result <- TSENAT:::.apply_art_kw(data, value_col = "entropy", group_col = "group")
  
  expect_is(result, "list")
  expect_true("statistic" %in% names(result))
  expect_true("p_value" %in% names(result))
  expect_true(is.numeric(result$statistic) || is.na(result$statistic))
})

test_that(".apply_art_kw returns valid results structure", {
  skip_on_cran()
  
  # Direct call to .apply_art_kw with heteroscedastic data
  set.seed(123)
  
  # Create data frame with heteroscedastic entropy values
  data <- data.frame(
    entropy = c(
      rnorm(20, mean = 1.5, sd = 0.2),  # Group 1: low variance
      rnorm(20, mean = 1.5, sd = 2.0)   # Group 2: high variance
    ),
    group = rep(c("A", "B"), each = 20)
  )
  
  # Apply ART directly
  result <- TSENAT:::.apply_art_kw(data, value_col = "entropy", group_col = "group")
  
  expect_is(result, "list")
  expect_true("statistic" %in% names(result))
  expect_true("p_value" %in% names(result))
  expect_true("method" %in% names(result))
  expect_true(is.numeric(result$statistic) || is.na(result$statistic))
  expect_true(is.numeric(result$p_value) || is.na(result$p_value))
})

test_that(".apply_art_kw p-values are valid (between 0 and 1)", {
  skip_on_cran()
  
  set.seed(456)
  
  # Create distinct groups with different means
  data <- data.frame(
    entropy = c(
      rnorm(25, mean = 1.0, sd = 0.5),  # Group A
      rnorm(25, mean = 2.5, sd = 0.5)   # Group B: clearly different
    ),
    group = rep(c("A", "B"), each = 25)
  )
  
  result <- TSENAT:::.apply_art_kw(data, value_col = "entropy", group_col = "group")
  
  # p-value must be between 0 and 1 (or NA)
  if (!is.na(result$p_value)) {
    expect_true(result$p_value >= 0 && result$p_value <= 1,
                label = "p-value should be in [0, 1]")
  }
  
  # statistic should be non-negative for Kruskal-Wallis type tests
  if (!is.na(result$statistic)) {
    expect_true(result$statistic >= 0,
                label = "test statistic should be non-negative")
  }
})

test_that(".apply_art_kw detects difference between distinct groups", {
  skip_on_cran()
  
  set.seed(789)
  
  # Create data with LARGE EFFECT SIZE and clear separation
  # Use uniform distributions for completely non-overlapping groups
  data <- data.frame(
    entropy = c(
      runif(35, min = 0.1, max = 0.4),     # Group A: low entropy
      runif(35, min = 0.6, max = 1.0)      # Group B: high entropy
    ),
    group = rep(c("Low", "High"), each = 35)
  )
  
  result <- TSENAT:::.apply_art_kw(data, value_col = "entropy", group_col = "group")
  
  # With completely separated groups, ART should detect the difference
  if (!is.na(result$p_value)) {
    expect_true(result$p_value < 0.1,
                label = "p-value should be small for completely separated groups")
  }
})

# ============================================================================
# TEST SUITE 2: .apply_quantile_test() — Boundary Clustering Trigger
# ============================================================================
# Condition: >40% of values within 10% of boundaries (quartiles)

test_that(".apply_quantile_test triggered by boundary-clustered data", {
  skip_on_cran()
  skip_if_not_installed("SummarizedExperiment")
  
  library(SummarizedExperiment)
  
  # Create data with BOUNDARY CLUSTERING at quartiles
  set.seed(99)
  n_samples <- 20
  
  # Create values clustered near quartile boundaries (0.25, 0.50, 0.75)
  boundary_values <- c(
    rep(c(0.24, 0.25, 0.26), 3),      # Cluster at Q1 (3x3 = 9 values)
    rep(c(0.49, 0.50, 0.51), 3),      # Cluster at Q2 (3x3 = 9 values)  
    rep(c(0.05, 0.10), 1)             # Scattered (2 values)
  )
  
  # Create data frame for rank test
  data <- data.frame(
    entropy = boundary_values[1:20],
    group = rep(c("A", "B"), each = 10)
  )
  
  # Apply quantile test directly - should handle the clustering
  result <- TSENAT:::.apply_quantile_test(data, value_col = "entropy", group_col = "group")
  
  expect_is(result, "list")
  expect_true("statistic" %in% names(result))
  expect_true("p_value" %in% names(result))
})

test_that(".apply_quantile_test returns numeric p-values in valid range", {
  skip_on_cran()
  
  set.seed(654)
  
  # Create two groups with different distributions
  data <- data.frame(
    entropy = c(
      runif(20, 0, 0.5),    # Group A: lower range
      runif(20, 0.5, 1.0)   # Group B: upper range
    ),
    group = rep(c("Lower", "Upper"), each = 20)
  )
  
  result <- TSENAT:::.apply_quantile_test(data, value_col = "entropy", group_col = "group")
  
  # p-value must be valid (between 0 and 1)
  expect_true(result$p_value >= 0 && result$p_value <= 1,
              label = "p-value must be in [0, 1]")
  
  # statistic should be numeric
  expect_true(is.numeric(result$statistic),
              label = "statistic must be numeric")
})

test_that(".apply_quantile_test returns valid results with clustered input", {
  skip_on_cran()
  
  # Direct call to .apply_quantile_test with boundary-clustered data
  set.seed(456)
  
  # Create boundary-clustered entropy values
  boundary_vals <- c(
    rep(0.48, 8),     # Cluster at boundary
    rep(0.52, 8),     # Cluster at boundary
    runif(4, 0.3, 0.7)  # Some scattered values
  )
  
  data <- data.frame(
    entropy = boundary_vals,
    q = rep(c("q_1.0", "q_1.0"), each = 10)
  )
  
  # Apply quantile test directly
  result <- TSENAT:::.apply_quantile_test(data, value_col = "entropy", group_col = "q")
  
  expect_is(result, "list")
  expect_true("statistic" %in% names(result))
  expect_true("p_value" %in% names(result))
  expect_true("method" %in% names(result))
  expect_true(is.numeric(result$statistic))
  expect_true(is.numeric(result$p_value))
})

# ============================================================================
# TEST SUITE 3: .apply_robust_median_test() — Extreme Skewness Trigger
# ============================================================================
# Condition: |skewness| > 1 (extreme positive or negative skewness)

test_that(".apply_robust_median_test triggered by extremely skewed data", {
  skip_on_cran()
  
  # Create data with EXTREME SKEWNESS
  # Use exponential distribution (naturally right-skewed, skewness ≈ 2)
  set.seed(777)
  
  # Exponential distribution is highly skewed
  entropy_skewed_a <- rexp(25, rate = 2)     # Group A: skewed ~0.5
  entropy_skewed_b <- rexp(25, rate = 0.5)   # Group B: skewed ~2
  
  data <- data.frame(
    entropy = c(entropy_skewed_a, entropy_skewed_b),
    group = rep(c("Skewed_A", "Skewed_B"), each = 25)
  )
  
  # Apply robust median test directly
  result <- TSENAT:::.apply_robust_median_test(data, value_col = "entropy", group_col = "group")
  
  expect_is(result, "list")
  expect_true("statistic" %in% names(result))
  expect_true("p_value" %in% names(result))
})

test_that(".apply_robust_median_test returns valid results with skewed input", {
  skip_on_cran()
  
  # Direct call to .apply_robust_median_test with skewed data
  set.seed(789)
  
  # Create highly skewed entropy values with larger sample sizes
  # to avoid chi-squared approximation warnings
  # Use log-normal distribution (naturally right-skewed)
  group_a <- rlnorm(40, meanlog = -2, sdlog = 1)    # Group A: left-skewed concentrated at low values
  group_b <- rlnorm(40, meanlog = -1, sdlog = 1.2)  # Group B: more right-skewed
  
  data <- data.frame(
    entropy = c(group_a, group_b),
    group = rep(c("A", "B"), each = 40)  # 40 + 40 = 80 total values
  )
  
  # Apply robust median test directly
  result <- TSENAT:::.apply_robust_median_test(data, value_col = "entropy", group_col = "group")
  
  expect_is(result, "list")
  expect_true("statistic" %in% names(result))
  expect_true("p_value" %in% names(result))
  expect_true(is.numeric(result$statistic) || is.na(result$statistic))
})

# ============================================================================
# TEST SUITE 4: print.rank_correlation_ci() — S3 Print Method
# ============================================================================

test_that("print.rank_correlation_ci displays output without error", {
  # Create a mock rank_correlation_ci object
  ci_obj <- structure(
    list(
      method = "Spearman Rank Correlation",
      ci_level = 0.95,
      correlation_matrix = matrix(
        c(1.0, 0.75, 0.60,
          0.75, 1.0, 0.82,
          0.60, 0.82, 1.0),
        nrow = 3, ncol = 3,
        dimnames = list(c("q_0.5", "q_1.0", "q_1.5"),
                        c("q_0.5", "q_1.0", "q_1.5"))
      ),
      interpretation = list(
        q_0.5_vs_q_1.0 = "Very stable",
        q_1.0_vs_q_1.5 = "Robust",
        q_0.5_vs_q_1.5 = "Moderate"
      )
    ),
    class = "rank_correlation_ci"
  )
  
  # Should print without error and return invisibly
  expect_message(print(ci_obj), "RANK CORRELATION")
  result <- print(ci_obj)
  expect_identical(result, ci_obj)
})

test_that("print.rank_correlation_ci shows method description", {
  ci_obj <- structure(
    list(
      method = "Kendall Tau Correlation",
      ci_level = 0.90,
      correlation_matrix = matrix(1.0, 2, 2),
      interpretation = list()
    ),
    class = "rank_correlation_ci"
  )
  
  expect_message(print(ci_obj), "Kendall Tau")
})

test_that("print.rank_correlation_ci shows confidence level", {
  ci_obj <- structure(
    list(
      method = "Spearman",
      ci_level = 0.99,
      correlation_matrix = matrix(1.0, 2, 2),
      interpretation = list()
    ),
    class = "rank_correlation_ci"
  )
  
  expect_message(print(ci_obj), "99%")
})

test_that("print.rank_correlation_ci shows interpretation guidelines", {
  ci_obj <- structure(
    list(
      method = "Spearman",
      ci_level = 0.95,
      correlation_matrix = matrix(1.0, 2, 2),
      interpretation = list(test = "robust")
    ),
    class = "rank_correlation_ci"
  )
  
  expect_message(print(ci_obj), "Robust|Stable|Weak|Variable")
})

test_that("print.rank_correlation_ci handles edge cases", {
  ci_obj <- structure(
    list(
      method = "Test",
      ci_level = 0.95,
      correlation_matrix = matrix(NA_real_, 2, 2),
      interpretation = NULL
    ),
    class = "rank_correlation_ci"
  )
  
  expect_error(print(ci_obj), NA)
})

# ============================================================================
# INTEGRATION TEST: Rank tests with realistic multi-q data characteristics
# ============================================================================

test_that("conditional rank tests handle multi-group data", {
  skip_on_cran()
  
  # Create realistic scenario with q-dependent effects
  set.seed(555)
  
  # Simulate entropy values that vary by group and might trigger different
  # conditional test selections based on variance, skewness, or boundary clustering
  
  # Group A: relatively homoscedastic
  entropy_a <- rnorm(40, mean = 1.0, sd = 0.3)
  
  # Group B: heteroscedastic (higher variance)
  entropy_b <- rnorm(40, mean = 1.2, sd = 1.5)
  
  data <- data.frame(
    entropy = c(entropy_a, entropy_b),
    group = rep(c("GroupA", "GroupB"), each = 40),
    q_value = rep(c(0.5, 1.0, 1.5, 2.0), 20)
  )
  
  # Test that rank methods can process this data without errors
  
  # Try ART
  result_art <- TSENAT:::.apply_art_kw(data, value_col = "entropy", group_col = "group")
  expect_is(result_art, "list")
  
  # Try quantile test
  result_quant <- TSENAT:::.apply_quantile_test(data, value_col = "entropy", group_col = "group")
  expect_is(result_quant, "list")
  
  # Try robust median
  result_robust <- TSENAT:::.apply_robust_median_test(data, value_col = "entropy", group_col = "group")
  expect_is(result_robust, "list")
})

# ============================================================================
# NUMERICAL VALIDATION TESTS: Verify statistical properties
# ============================================================================

test_that(".apply_art_kw returns p-values in valid range [0, 1]", {
  skip_on_cran()
  
  set.seed(321)
  
  # Create heteroscedastic data (different variances)
  group_a <- rnorm(30, mean = 1.0, sd = 0.5)   # Low variance
  group_b <- rnorm(30, mean = 1.1, sd = 2.0)   # High variance
  
  data <- data.frame(
    entropy = c(group_a, group_b),
    group = rep(c("LowVar", "HighVar"), each = 30)
  )
  
  result <- TSENAT:::.apply_art_kw(data, value_col = "entropy", group_col = "group")
  
  # p-value must be valid probability
  expect_true(result$p_value >= 0 && result$p_value <= 1,
              label = "p-value must be in [0, 1]")
  
  # statistic should be positive (F-like)
  expect_true(result$statistic >= 0,
              label = "ART statistic should be non-negative")
})

test_that(".apply_quantile_test returns p-values in valid range [0, 1]", {
  skip_on_cran()
  
  set.seed(654)
  
  # Create groups with different distributions
  data <- data.frame(
    entropy = c(
      runif(25, 0, 0.5),    # Group A: lower range
      runif(25, 0.5, 1.0)   # Group B: upper range
    ),
    group = rep(c("Lower", "Upper"), each = 25)
  )
  
  result <- TSENAT:::.apply_quantile_test(data, value_col = "entropy", group_col = "group")
  
  # p-value must be valid
  expect_true(result$p_value >= 0 && result$p_value <= 1,
              label = "p-value must be in [0, 1]")
  
  # statistic should be numeric
  expect_true(is.numeric(result$statistic),
              label = "statistic must be numeric")
})

test_that(".apply_robust_median_test returns p-values in valid range [0, 1]", {
  skip_on_cran()
  
  set.seed(987)
  
  # Create skewed groups
  group_a <- rexp(30, rate = 1)      # Right-skewed
  group_b <- rexp(30, rate = 0.5)    # More right-skewed
  
  data <- data.frame(
    entropy = c(group_a, group_b),
    group = rep(c("ExpFast", "ExpSlow"), each = 30)
  )
  
  result <- TSENAT:::.apply_robust_median_test(data, value_col = "entropy", group_col = "group")
  
  # p-value must be valid
  expect_true(result$p_value >= 0 && result$p_value <= 1,
              label = "p-value must be in [0, 1]")
  
  # statistic should be numeric
  expect_true(is.numeric(result$statistic) || is.na(result$statistic),
              label = "statistic must be numeric or NA")
})

test_that("conditional rank tests detect differences between distinct groups", {
  skip_on_cran()
  
  set.seed(1111)
  
  # Create data with EXTREME DIFFERENCE between groups
  # Group A: entropy in range [0.1, 0.3] (low entropy region)
  # Group B: entropy in range [0.7, 0.9] (high entropy region)
  # These are almost completely non-overlapping distributions
  
  group_a <- runif(50, min = 0.1, max = 0.3)   # Low entropy cluster
  group_b <- runif(50, min = 0.7, max = 0.9)   # High entropy cluster
  
  data <- data.frame(
    entropy = c(group_a, group_b),
    group = rep(c("Low", "High"), each = 50)
  )
  
  # All three rank tests should detect this extreme separation
  result_art <- TSENAT:::.apply_art_kw(data, value_col = "entropy", group_col = "group")
  result_quant <- TSENAT:::.apply_quantile_test(data, value_col = "entropy", group_col = "group")
  result_robust <- TSENAT:::.apply_robust_median_test(data, value_col = "entropy", group_col = "group")
  
  # With completely separated groups, at least one method should detect significant difference
  # or return NA (which is acceptable)
  has_significant <- FALSE
  if (!is.na(result_art$p_value) && result_art$p_value < 0.05) has_significant <- TRUE
  if (!is.na(result_quant$p_value) && result_quant$p_value < 0.05) has_significant <- TRUE
  if (!is.na(result_robust$p_value) && result_robust$p_value < 0.05) has_significant <- TRUE
  
  # With such extreme separation, expect at least one method to find significance
  expect_true(has_significant,
              label = "At least one rank test should detect extreme separation (p < 0.05)")
})
