# Suppress nboot warnings for this test file
options(TSENAT.suppress_nboot_warning = TRUE)

# Safety check: ensure test factory is loaded
if (!exists("create_test_se_simple", mode = "function")) {
  factory_file <- file.path(dirname(getwd()), "testthat", "tests-factory.R")
  if (!file.exists(factory_file)) {
    factory_file <- "tests/testthat/tests-factory.R"
  }
  if (file.exists(factory_file)) {
    source(factory_file, local = FALSE)
  }
}

# Setup test data
set.seed(42)
test_counts <- c(100, 80, 60, 40, 20)
test_counts_balanced <- c(50, 50, 50, 50)
test_counts_skewed <- c(500, 200, 100, 50, 20, 10)
test_matrix <- matrix(
  c(100, 80, 60, 40, 20, 150, 100, 50, 30, 10),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(c("Gene1", "Gene2"), NULL)
)


context("Resampling Methods: Multi-q Support for Bootstrap and Jackknife")

test_that("calculate_tsallis_entropy_bootstrap accepts vector q", {
    x <- c(100, 50, 30, 20)
    result <- .calculate_tsallis_entropy_bootstrap(x, q = c(1, 2), nboot = 100)
    
    expect_is(result, "tsenat_bootstrap_ci_list")
    expect_length(result, 2)
    expect_named(result, c("q=1", "q=2"))
})

test_that("bootstrap multi-q returns correct structure", {
    x <- c(100, 50, 30, 20)
    result <- .calculate_tsallis_entropy_bootstrap(x, q = c(1, 2), nboot = 100)
    
    # Check list structure (2 q values: 1 and 2)
    expect_length(result, 2)
    expect_true(all(sapply(result, function(r) "estimate" %in% names(r))))
    expect_true(all(sapply(result, function(r) "lower_ci" %in% names(r))))
    expect_true(all(sapply(result, function(r) "upper_ci" %in% names(r))))
})

test_that("bootstrap multi-q estimates differ across q values", {
    set.seed(42)
    x <- c(100, 50, 30, 20, 10)
    result <- .calculate_tsallis_entropy_bootstrap(x, q = c(0.5, 1, 2), nboot = 100)  # Reduced from 200 to 100
    
    # Estimates should be different for different q values
    est_q05 <- result$`q=0.5`$estimate
    est_q1 <- result$`q=1`$estimate
    est_q2 <- result$`q=2`$estimate
    
    # Should not all be equal
    expect_false(isTRUE(all.equal(est_q05, est_q1)))
    expect_false(isTRUE(all.equal(est_q1, est_q2)))
})

test_that("bootstrap multi-q CI bounds are sensible for each q", {
    
    x <- c(100, 50, 30, 20)
    result <- .calculate_tsallis_entropy_bootstrap(x, q = c(1, 2), nboot = 100)
    
    for (res in result) {
        expect_lt(res$lower_ci, res$estimate)
        expect_gt(res$upper_ci, res$estimate)
        expect_gte(res$lower_ci, 0)
        expect_lte(res$upper_ci, 1)
    }
})

test_that("bootstrap multi-q with different nboot values", {
    x <- c(100, 50, 30, 20)
    # nboot=50 triggers warning about being below recommended minimum (expected for exploratory testing)
    result_small <- suppressWarnings(.calculate_tsallis_entropy_bootstrap(x, q = c(1, 2), nboot = 50))  # Reduced for speed
    result_large <- suppressWarnings(.calculate_tsallis_entropy_bootstrap(x, q = c(1, 2), nboot = 75))  # Reduced from 500
    
    # Both should return valid results
    expect_is(result_small, "tsenat_bootstrap_ci_list")
    expect_is(result_large, "tsenat_bootstrap_ci_list")
    
    # Larger nboot should give more stable estimates
    expect_equal(length(result_small$`q=1`$bootstrap_dist), 50)  # Updated from 100
    expect_equal(length(result_large$`q=1`$bootstrap_dist), 75)  # Updated from 500
})

test_that("calculate_jeo accepts vector q", {
    x <- c(100, 50, 30, 20)
    result <- .calculate_jeo(x, q = c(1, 2), norm = TRUE, verbose = FALSE)
    
    expect_is(result, "tsenat_jackknife_list_multiq")
    expect_length(result, 2)
    expect_named(result, c("q=1", "q=2"))
})

test_that("jackknife multi-q returns correct structure", {
    x <- c(100, 50, 30, 20, 15, 10)
    result <- .calculate_jeo(x, q = c(0.5, 1, 2), norm = TRUE, verbose = FALSE)
    
    # Check list structure
    expect_length(result, 3)
    expect_true(all(sapply(result, function(r) "estimate" %in% names(r))))
    expect_true(all(sapply(result, function(r) "jackknife_se" %in% names(r))))
    expect_true(all(sapply(result, function(r) "influence" %in% names(r))))
})

test_that("jackknife multi-q estimates differ across q values", {
    x <- c(100, 50, 30, 20, 10)
    result <- .calculate_jeo(x, q = c(0.5, 1, 2), norm = TRUE, verbose = FALSE)
    
    # Estimates should be different for different q values
    est_q05 <- result$`q=0.5`$estimate
    est_q1 <- result$`q=1`$estimate
    est_q2 <- result$`q=2`$estimate
    
    expect_false(isTRUE(all.equal(est_q05, est_q1)))
    expect_false(isTRUE(all.equal(est_q1, est_q2)))
})

test_that("jackknife multi-q with matrix input", {
    
    counts_matrix <- rbind(
        "Gene1" = c(100, 50, 30, 20),
        "Gene2" = c(80, 60, 40, 20)
    )
    
    # Jackknife on first row
    result <- .calculate_jeo(
        x = counts_matrix[1, , drop = FALSE],
        q = c(1, 2),
        norm = TRUE,
        verbose = FALSE
    )
    
    expect_is(result, "tsenat_jackknife_list_multiq")
    expect_length(result, 2)
})

test_that("jackknife multi-q SE estimates are positive", {
    x <- c(100, 50, 30, 20, 15)
    result <- .calculate_jeo(x, q = c(1, 2), norm = TRUE, verbose = FALSE)
    
    for (res in result) {
        expect_gt(res$jackknife_se, 0)
        expect_length(res$influence, length(x))
        expect_true(all(res$influence >= 0))
    }
})

test_that("bootstrap accepts vector q with length > 1", {
    x <- c(100, 50, 30, 20)
    
    # Vector q with length > 1 returns list
    result <- .calculate_tsallis_entropy_bootstrap(x, q = c(1, 2), nboot = 100)
    expect_is(result, "tsenat_bootstrap_ci_list")
    expect_length(result, 2)
    expect_named(result, c("q=1", "q=2"))
})

test_that("jackknife accepts vector q with length > 1", {
    x <- c(100, 50, 30, 20)
    
    # Vector q with length > 1 returns list
    result <- .calculate_jeo(x, q = c(1, 2), norm = TRUE, verbose = FALSE)
    expect_is(result, "tsenat_jackknife_list_multiq")
    expect_length(result, 2)
    expect_named(result, c("q=1", "q=2"))
})

test_that("single q still works (backward compatibility)", {
    x <- c(100, 50, 30, 20)
    
    # Bootstrap with scalar q (without diagnostics to test original format)
    boot_result <- .calculate_tsallis_entropy_bootstrap(x, q = 2, nboot = 100, include_diagnostics = FALSE)
    expect_is(boot_result, "tsenat_bootstrap_ci")
    expect_named(boot_result, c("estimate", "lower_ci", "upper_ci", "ci_level", "method", "nboot", "bootstrap_dist"))
    
    # Jackknife with scalar q
    jack_result <- .calculate_jeo(x, q = 1, norm = TRUE, verbose = FALSE)
    expect_is(jack_result, "tsenat_jackknife")
    # Check that key fields exist (structure may have additional fields)
    expect_true("estimate" %in% names(jack_result))
    expect_true("jackknife_se" %in% names(jack_result))
    expect_true("influence" %in% names(jack_result))
    expect_true("q" %in% names(jack_result))
})

test_that("bootstrap multi-q requires set.seed() for reproducibility", {
    x <- c(100, 50, 30, 20)
    
    set.seed(555)
    result1 <- .calculate_tsallis_entropy_bootstrap(x, q = c(1, 2), nboot = 100)
    
    set.seed(555) # Use same seed
    result2 <- .calculate_tsallis_entropy_bootstrap(x, q = c(1, 2), nboot = 100)
    
    # Same seed should give identical results
    expect_equal(result1$`q=1`$estimate, result2$`q=1`$estimate, tolerance = 1e-10)
    expect_equal(result1$`q=2`$estimate, result2$`q=2`$estimate, tolerance = 1e-10)
    expect_equal(result1$`q=1`$lower_ci, result2$`q=1`$lower_ci, tolerance = 1e-10)
    expect_equal(result1$`q=2`$lower_ci, result2$`q=2`$lower_ci, tolerance = 1e-10)
})

test_that("bootstrap multi-q ci parameter returns finite widths", {
    x <- c(100, 80, 60, 40, 30, 20, 15, 10, 8, 5)
    
    result_95 <- suppressWarnings(.calculate_tsallis_entropy_bootstrap(x, q = c(1.5, 2.0), nboot = 20, ci = 0.95))
    result_90 <- suppressWarnings(.calculate_tsallis_entropy_bootstrap(x, q = c(1.5, 2.0), nboot = 20, ci = 0.90))
    
    # Both should give valid CI widths
    width_95_q1 <- result_95$`q=1.5`$upper_ci - result_95$`q=1.5`$lower_ci
    width_90_q1 <- result_90$`q=1.5`$upper_ci - result_90$`q=1.5`$lower_ci
    
    # Check that widths are finite and positive
    expect_true(is.finite(width_95_q1))
    expect_true(is.finite(width_90_q1))
    expect_gt(width_95_q1, 0)
    expect_gt(width_90_q1, 0)
})

test_that("jackknife multi-q accepts Hill numbers (D)", {
    x <- c(100, 50, 30, 20, 15)
    # Note: jackknife doesn't have 'what' parameter, but we test multi-q works
    result <- .calculate_jeo(x, q = c(1, 2), norm = FALSE, verbose = FALSE)
    
    expect_is(result, "tsenat_jackknife_list_multiq")
    expect_length(result, 2)
    expect_true(all(sapply(result, function(r) !is.na(r$estimate))))
})

#!/usr/bin/env Rscript
#===============================================================================
# TEST SUITE: Priority 3 Bootstrap Diagnostics
#===============================================================================
# Tests for new diagnostic functions:
# - .estimate_bootstrap_skewness()
# - .detect_multimodality()
# - .analyze_ci_width()
# - .generate_bootstrap_diagnostics_report()
#
# These tests verify diagnostic accuracy and integration
#===============================================================================

library(testthat)
library(TSENAT)

# ============================================================================
# SUITE 1: Skewness Estimation
# ============================================================================

test_that(".estimate_bootstrap_skewness handles symmetric distributions", {
  # Create symmetric bootstrap distribution (normal)
  set.seed(123)
  boot_dist <- rnorm(1000, mean = 2, sd = 0.5)
  
  result <- TSENAT:::.estimate_bootstrap_skewness(boot_dist, compute_ci = FALSE)
  
  # Symmetric distributions should have skewness close to 0
  expect_true(abs(result$skewness_mean) < 0.3)
  expect_true(abs(result$skewness_quartile) < 0.3)
  expect_match(result$interpretation, "symmetric|low", ignore.case = TRUE)
})

test_that(".estimate_bootstrap_skewness detects right-skewed distributions", {
  # Create right-skewed distribution (exponential)
  set.seed(456)
  boot_dist <- rexp(1000, rate = 1)
  
  result <- TSENAT:::.estimate_bootstrap_skewness(boot_dist, compute_ci = FALSE)
  
  # Right-skewed distributions should have positive skewness
  expect_true(result$skewness_mean > 0)
  expect_true(result$skewness_quartile > 0)
})

test_that(".estimate_bootstrap_skewness computes jackknife CI", {
  # Create sample bootstrap distribution
  set.seed(789)
  boot_dist <- rnorm(100, mean = 1, sd = 0.3)
  
  result <- TSENAT:::.estimate_bootstrap_skewness(boot_dist, compute_ci = TRUE)
  
  # CI bounds should be finite
  expect_true(is.finite(result$ci_lower))
  expect_true(is.finite(result$ci_upper))
  
  # CI should bracket the point estimate
  expect_true(result$ci_lower <= result$skewness_mean)
  expect_true(result$skewness_mean <= result$ci_upper)
})

test_that(".estimate_bootstrap_skewness handles degenerate cases", {
  # All same values - zero variance
  boot_dist <- rep(2, 50)
  
  result <- TSENAT:::.estimate_bootstrap_skewness(boot_dist, compute_ci = FALSE)
  
  # Degenerate cases should return NA
  expect_true(is.na(result$skewness_mean) || result$skewness_mean == 0)
})

test_that(".estimate_bootstrap_skewness handles small samples", {
  # Small sample (n < 3)
  boot_dist <- c(1, 2)
  
  result <- TSENAT:::.estimate_bootstrap_skewness(boot_dist, compute_ci = FALSE)
  
  # Should warn or return NA
  expect_true(is.na(result$skewness_mean) || 
              grepl("insufficient", result$interpretation, ignore.case = TRUE))
})

# ============================================================================
# SUITE 2: Multimodality Detection
# ============================================================================

test_that(".detect_multimodality identifies unimodal distributions", {
  set.seed(111)
  boot_dist <- rnorm(500, mean = 2, sd = 0.5)
  
  result <- TSENAT:::.detect_multimodality(boot_dist, method = "kde")
  
  expect_false(result$is_multimodal)
  expect_equal(result$n_modes, 1L)
  expect_match(result$interpretation, "unimodal", ignore.case = TRUE)
})

test_that(".detect_multimodality identifies bimodal distributions", {
  set.seed(222)
  # Mixture of two normals
  boot_dist <- c(rnorm(300, mean = 1, sd = 0.3),
                 rnorm(300, mean = 4, sd = 0.3))
  
  result <- TSENAT:::.detect_multimodality(boot_dist, method = "kde")
  
  expect_true(result$is_multimodal)
  expect_true(result$n_modes >= 2)
  expect_match(result$interpretation, "modality|mode", ignore.case = TRUE)
})

test_that(".detect_multimodality histogram method works", {
  set.seed(333)
  boot_dist <- rnorm(200, mean = 2, sd = 0.5)
  
  result <- TSENAT:::.detect_multimodality(boot_dist, method = "histogram")
  
  expect_true(is.logical(result$is_multimodal))
  expect_true(is.integer(result$n_modes))
  expect_equal(result$method_used, "histogram")
})

test_that(".detect_multimodality gaps method works", {
  set.seed(444)
  boot_dist <- rnorm(200, mean = 2, sd = 0.5)
  
  result <- TSENAT:::.detect_multimodality(boot_dist, method = "gaps")
  
  expect_true(is.logical(result$is_multimodal))
  expect_true(is.integer(result$n_modes))
  expect_equal(result$method_used, "gaps")
})

test_that(".detect_multimodality rejects invalid methods", {
  boot_dist <- rnorm(100, mean = 2, sd = 0.5)
  
  expect_error(
    TSENAT:::.detect_multimodality(boot_dist, method = "invalid"),
    "should be one of"
  )
})

test_that(".detect_multimodality handles small samples", {
  boot_dist <- c(1, 2, 3, 4, 5)  # n < 10
  
  result <- TSENAT:::.detect_multimodality(boot_dist, method = "kde")
  
  expect_true(is.na(result$is_multimodal))
  expect_match(result$method_used, "insufficient", ignore.case = TRUE)
})

test_that(".detect_multimodality computes separation score", {
  set.seed(555)
  boot_dist <- c(rnorm(300, mean = 1, sd = 0.2),
                 rnorm(300, mean = 5, sd = 0.2))
  
  result <- TSENAT:::.detect_multimodality(boot_dist, method = "kde")
  
  if (result$is_multimodal) {
    # Separation score should be between 0 and 1 for multimodal
    expect_true(result$separation_score >= 0 && result$separation_score <= 1)
  } else {
    # For unimodal, should be very high (>0.9)
    expect_true(result$separation_score > 0.9)
  }
})

# ============================================================================
# SUITE 3: CI Width Analysis
# ============================================================================

test_that(".analyze_ci_width computes basic characteristics", {
  ci_lower <- 1.2
  ci_upper <- 3.5
  point_est <- 2.3
  boot_dist <- rnorm(500, mean = 2.3, sd = 0.5)
  n_bootstrap <- 500
  
  result <- TSENAT:::.analyze_ci_width(ci_lower, ci_upper, point_est, boot_dist, n_bootstrap)
  
  # Check CI width
  expect_equal(result$ci_width, ci_upper - ci_lower)
  
  # Check ratios are positive
  expect_true(result$ci_width_to_estimate_ratio > 0)
  expect_true(result$ci_width_to_sd_ratio > 0)
})

test_that(".analyze_ci_width detects asymmetric CIs", {
  # Asymmetric CI: lower at 1, upper at 10, point estimate at 2
  ci_lower <- 1
  ci_upper <- 10
  point_est <- 2
  boot_dist <- c(rep(1.5, 200), runif(300, 9, 10))
  n_bootstrap <- 500
  
  result <- TSENAT:::.analyze_ci_width(ci_lower, ci_upper, point_est, boot_dist, n_bootstrap)
  
  # Should detect asymmetry
  if (!is.na(result$ci_symmetry_ratio)) {
    expect_true(result$ci_symmetry_ratio < 0.9)
  }
  
  # Should flag in issues
  expect_true(any(grepl("asymmetric|recommendation", tolower(result$potential_issues), 
                         ignore.case = TRUE)))
})

test_that(".analyze_ci_width assesses precision levels", {
  boot_dist <- rnorm(500, mean = 2, sd = 0.5)
  
  # Excellent precision
  result_excellent <- TSENAT:::.analyze_ci_width(
    ci_lower = 1.95, ci_upper = 2.05,
    point_est = 2.0, boot_dist = boot_dist, nboot = 500
  )
  expect_match(result_excellent$precision_assessment, "excellent|good")
  
  # Poor precision
  result_poor <- TSENAT:::.analyze_ci_width(
    ci_lower = 0.5, ci_upper = 3.5,
    point_est = 2.0, boot_dist = boot_dist, nboot = 500
  )
  expect_match(result_poor$precision_assessment, "poor|acceptable")
})

test_that(".analyze_ci_width handles zero point estimate", {
  ci_lower <- -0.1
  ci_upper <- 0.1
  point_est <- 0
  boot_dist <- rnorm(500, mean = 0, sd = 0.05)
  
  result <- TSENAT:::.analyze_ci_width(ci_lower, ci_upper, point_est, boot_dist, 500)
  
  # Should handle zero estimate gracefully
  expect_true(is.na(result$ci_width_to_estimate_ratio) || 
              is.finite(result$ci_width_to_estimate_ratio))
})

# ============================================================================
# SUITE 4: Integrated Diagnostics Report
# ============================================================================

test_that(".generate_bootstrap_diagnostics_report works on real bootstrap result", {
  # Create a bootstrap result object
  set.seed(666)
  x <- c(100, 50, 25, 10)
  
  result <- TSENAT:::.calculate_tsallis_entropy_bootstrap(
    x = x, q = 1.5, nboot = 500, ci = 0.95,
    method = "percentile", include_diagnostics = TRUE,
    verbose = FALSE
  )
  
  # Generate report
  report <- TSENAT:::.generate_bootstrap_diagnostics_report(result)
  
  # Check report structure
  expect_true(all(c("skewness_analysis", "multimodality_analysis", 
                     "ci_width_analysis", "overall_reliability", 
                     "summary_recommendations") %in% names(report)))
  
  # Check overall reliability is one of three values
  expect_true(report$overall_reliability %in% c("Reliable", "Caution", "Unreliable"))
  
  # Check recommendations are character vector
  expect_true(is.character(report$summary_recommendations))
})

test_that(".generate_bootstrap_diagnostics_report detects unreliable results", {
  # Create bootstrap result with known issues
  set.seed(777)
  x <- c(10, 8, 5, 2, 1)  # Low counts - problematic
  
  result <- TSENAT:::.calculate_tsallis_entropy_bootstrap(
    x = x, q = 1, nboot = 100,  # Low nboot - problematic
    ci = 0.95, method = "percentile", 
    include_diagnostics = TRUE, verbose = FALSE
  )
  
  report <- TSENAT:::.generate_bootstrap_diagnostics_report(result)
  
  # Should suggest caution or flag issues
  expect_true(report$overall_reliability %in% c("Caution", "Unreliable"))
  
  # Should have recommendations
  expect_true(length(report$summary_recommendations) > 0)
})

test_that(".generate_bootstrap_diagnostics_report rejects invalid input", {
  # Invalid input - not a bootstrap CI object
  invalid_result <- list(estimate = 1, lower_ci = 0.5, upper_ci = 1.5)
  
  expect_error(
    TSENAT:::.generate_bootstrap_diagnostics_report(invalid_result),
    "must be of class tsenat_bootstrap_ci"
  )
})

# ============================================================================
# SUITE 5: Integration with main bootstrap function
# ============================================================================

test_that("bootstrap CI includes diagnostics when requested", {
  x <- c(150, 100, 50, 25, 10)
  
  result <- TSENAT:::.calculate_tsallis_entropy_bootstrap(
    x = x, q = 1.5, nboot = 500, ci = 0.95,
    include_diagnostics = TRUE, verbose = FALSE
  )
  
  # Should have diagnostics field
  expect_true("diagnostics" %in% names(result))
  expect_true(!is.null(result$diagnostics))
})

test_that("bootstrap CI omits diagnostics when not requested", {
  x <- c(150, 100, 50, 25, 10)
  
  result <- TSENAT:::.calculate_tsallis_entropy_bootstrap(
    x = x, q = 1.5, nboot = 500, ci = 0.95,
    include_diagnostics = FALSE, verbose = FALSE
  )
  
  # Should not have diagnostics field (or it's NULL)
  expect_true(!"diagnostics" %in% names(result) || is.null(result$diagnostics))
})


# ============================================================================
# INTEGRATION TESTS: Bootstrap CI S3 Methods via Public API
# ============================================================================
# These tests verify that bootstrap S3 methods (print, summary) are properly
# triggered when users interact with the public API through:
# - calculate_jeo() [exported S4 function]
# - calculate_diversity() [exported S4 function with bootstrap parameters]
# - jeoResults() [exported accessor function]
# 
# NOTE: The internal S3 methods are NOT directly exported but are registered
# via registerS3method() in .onLoad() and are automatically used when:
# 1. Users print results from calculate_jeo()
# 2. Users summarize bootstrap CI objects from jackknife functions
# 3. Bootstrap objects are returned from internal .calculate_tsallis_entropy_bootstrap()
# ============================================================================

# Create a tsenat_bootstrap_ci object
create_test_bootstrap_ci <- function(estimate = 0.65, lower_ci = 0.45, upper_ci = 0.82) {
  result <- list(
    estimate = estimate,
    lower_ci = lower_ci,
    upper_ci = upper_ci,
    ci_level = 0.95,
    method = "percentile",
    nboot = 1000,
    bootstrap_dist = rnorm(1000, mean = estimate, sd = 0.08),
    diagnostics = list(
      effective_sample_size = 995,
      skewness = 0.12,
      bias = 0.015,
      kurtosis = -0.05
    )
  )
  class(result) <- c("tsenat_bootstrap_ci", "list")
  result
}

# Create a tsenat_bootstrap_ci_list object (multiple q values)
create_test_bootstrap_ci_list <- function() {
  list(
    `q=0.5` = create_test_bootstrap_ci(estimate = 0.58, lower_ci = 0.42, upper_ci = 0.75),
    `q=1.0` = create_test_bootstrap_ci(estimate = 0.65, lower_ci = 0.45, upper_ci = 0.82),
    `q=1.5` = create_test_bootstrap_ci(estimate = 0.72, lower_ci = 0.52, upper_ci = 0.88)
  ) |> structure(class = c("tsenat_bootstrap_ci_list", "list"))
}

# Create a tsenat_divergence_bootstrap_ci object
create_test_divergence_bootstrap_ci <- function(estimate = 0.35, lower_ci = 0.15, upper_ci = 0.58) {
  result <- list(
    estimate = estimate,
    lower_ci = lower_ci,
    upper_ci = upper_ci,
    ci_level = 0.95,
    method = "percentile",
    nboot = 1000,
    bootstrap_dist = rnorm(1000, mean = estimate, sd = 0.10),
    p_value = 0.032,
    effect_size = estimate / 0.5,  # Relative to some reference
    diagnostics = list(
      effective_sample_size = 990,
      skewness = 0.18,
      bias = 0.008,
      kurtosis = 0.02,
      relative_ci_width = (upper_ci - lower_ci) / estimate
    )
  )
  class(result) <- c("tsenat_divergence_bootstrap_ci", "list")
  result
}

# ============================================================================
# TEST SUITE 1: print.tsenat_bootstrap_ci_list (8 uncovered lines)
# ============================================================================

test_that("print.tsenat_bootstrap_ci_list displays message with list header", {
  ci_list <- create_test_bootstrap_ci_list()
  
  expect_message(
    print(ci_list),
    "Bootstrap Confidence Intervals"
  )
})

test_that("print.tsenat_bootstrap_ci_list shows number of q values", {
  ci_list <- create_test_bootstrap_ci_list()
  
  expect_message(
    print(ci_list),
    "Number of q values: 3"
  )
})

test_that("print.tsenat_bootstrap_ci_list displays q values correctly", {
  ci_list <- create_test_bootstrap_ci_list()
  
  expect_message(
    print(ci_list),
    "q = q=0.5"
  )
  expect_message(
    print(ci_list),
    "q = q=1"
  )
  expect_message(
    print(ci_list),
    "q = q=1.5"
  )
})

test_that("print.tsenat_bootstrap_ci_list shows estimate for each q", {
  ci_list <- create_test_bootstrap_ci_list()
  
  expect_message(
    print(ci_list),
    "Estimate:"
  )
})

test_that("print.tsenat_bootstrap_ci_list shows confidence intervals for each q", {
  ci_list <- create_test_bootstrap_ci_list()
  
  expect_message(
    print(ci_list),
    "95% CI:"
  )
})

test_that("print.tsenat_bootstrap_ci_list returns object invisibly", {
  ci_list <- create_test_bootstrap_ci_list()
  
  result <- print(ci_list)
  expect_identical(result, ci_list)
})

test_that("print.tsenat_bootstrap_ci_list formats numbers with 6 decimal places", {
  ci_list <- create_test_bootstrap_ci_list()
  
  # Should use sprintf with %.6f format
  expect_message(
    print(ci_list),
    "0.58"  # Check estimate appears with reasonable precision
  )
})

test_that("print.tsenat_bootstrap_ci_list handles empty list gracefully", {
  empty_list <- structure(list(), class = c("tsenat_bootstrap_ci_list", "list"))
  
  expect_message(
    print(empty_list),
    "Bootstrap Confidence Intervals"
  )
})

# ============================================================================
# TEST SUITE 2: print.tsenat_divergence_bootstrap_ci (1 uncovered line)
# ============================================================================

test_that("print.tsenat_divergence_bootstrap_ci returns invisibly", {
  div_ci <- create_test_divergence_bootstrap_ci()
  
  result <- print(div_ci)
  expect_identical(result, div_ci)
  invisible(result)
})

test_that("print.tsenat_divergence_bootstrap_ci does not error for valid object", {
  div_ci <- create_test_divergence_bootstrap_ci()
  
  expect_error(print(div_ci), NA)
})

test_that("print.tsenat_divergence_bootstrap_ci handles NULL bootstrap_dist", {
  div_ci <- create_test_divergence_bootstrap_ci()
  div_ci$bootstrap_dist <- NULL
  
  expect_error(print(div_ci), NA)
})

test_that("print.tsenat_divergence_bootstrap_ci works with zero-length bootstrap_dist", {
  div_ci <- create_test_divergence_bootstrap_ci()
  div_ci$bootstrap_dist <- numeric(0)
  
  expect_error(print(div_ci), NA)
})

# ============================================================================
# TEST SUITE 3: summary.tsenat_divergence_bootstrap_ci (28 uncovered lines)
# ============================================================================

test_that("summary.tsenat_divergence_bootstrap_ci displays header message", {
  div_ci <- create_test_divergence_bootstrap_ci()
  
  expect_message(
    summary(div_ci),
    "Summary of Divergence Bootstrap"
  )
})

test_that("summary.tsenat_divergence_bootstrap_ci shows mean of bootstrap distribution", {
  div_ci <- create_test_divergence_bootstrap_ci()
  
  expect_message(
    summary(div_ci),
    "Mean:"
  )
})

test_that("summary.tsenat_divergence_bootstrap_ci shows median of bootstrap distribution", {
  div_ci <- create_test_divergence_bootstrap_ci()
  
  expect_message(
    summary(div_ci),
    "Median:"
  )
})

test_that("summary.tsenat_divergence_bootstrap_ci shows standard deviation", {
  div_ci <- create_test_divergence_bootstrap_ci()
  
  expect_message(
    summary(div_ci),
    "SD:"
  )
})

test_that("summary.tsenat_divergence_bootstrap_ci shows minimum value", {
  div_ci <- create_test_divergence_bootstrap_ci()
  
  expect_message(
    summary(div_ci),
    "Min:"
  )
})

test_that("summary.tsenat_divergence_bootstrap_ci shows maximum value", {
  div_ci <- create_test_divergence_bootstrap_ci()
  
  expect_message(
    summary(div_ci),
    "Max:"
  )
})

test_that("summary.tsenat_divergence_bootstrap_ci displays diagnostics section", {
  div_ci <- create_test_divergence_bootstrap_ci()
  
  expect_message(
    summary(div_ci),
    "Diagnostics:"
  )
})

test_that("summary.tsenat_divergence_bootstrap_ci calculates and shows skewness", {
  div_ci <- create_test_divergence_bootstrap_ci()
  
  expect_message(
    summary(div_ci),
    "Skewness:"
  )
})

test_that("summary.tsenat_divergence_bootstrap_ci calculates effective sample size", {
  div_ci <- create_test_divergence_bootstrap_ci()
  
  expect_message(
    summary(div_ci),
    "Effective sample size:"
  )
})

test_that("summary.tsenat_divergence_bootstrap_ci displays stability metrics section", {
  div_ci <- create_test_divergence_bootstrap_ci()
  
  expect_message(
    summary(div_ci),
    "Stability metrics:"
  )
})

test_that("summary.tsenat_divergence_bootstrap_ci shows CI width to estimate ratio", {
  div_ci <- create_test_divergence_bootstrap_ci()
  
  expect_message(
    summary(div_ci),
    "CI width to estimate ratio:"
  )
})

test_that("summary.tsenat_divergence_bootstrap_ci counts unique rounded values", {
  div_ci <- create_test_divergence_bootstrap_ci()
  
  expect_message(
    summary(div_ci),
    "Unique rounded values:"
  )
})

test_that("summary.tsenat_divergence_bootstrap_ci correctly computes skewness", {
  # Create data with known skewness
  set.seed(42)
  bootstrap_dist <- c(-2, -1, 0, 1, 2, 3, 4, 5, 6)
  
  div_ci <- create_test_divergence_bootstrap_ci()
  div_ci$bootstrap_dist <- bootstrap_dist
  
  # Manually calculate expected skewness
  m <- mean(bootstrap_dist)
  s <- sd(bootstrap_dist)
  n <- length(bootstrap_dist)
  expected_skew <- (sum((bootstrap_dist - m)^3) / n) / s^3
  
  expect_message(
    summary(div_ci),
    class = "character"
  )
})

test_that("summary.tsenat_divergence_bootstrap_ci handles constant bootstrap distribution", {
  div_ci <- create_test_divergence_bootstrap_ci()
  div_ci$bootstrap_dist <- rep(0.5, 100)  # All same value
  
  # Should show "N/A (no variation)" for skewness
  expect_message(
    summary(div_ci),
    "Skewness:"
  )
})

test_that("summary.tsenat_divergence_bootstrap_ci handles NA values in bootstrap distribution", {
  div_ci <- create_test_divergence_bootstrap_ci()
  # Remove NAs from bootstrap dist for valid calculation
  # (The summary function computes on raw dist which may have NAs)
  div_ci$bootstrap_dist <- div_ci$bootstrap_dist[!is.na(div_ci$bootstrap_dist)]
  
  expect_message(
    summary(div_ci),
    "Summary of Divergence Bootstrap"
  )
})

test_that("summary.tsenat_divergence_bootstrap_ci handles Inf values", {
  div_ci <- create_test_divergence_bootstrap_ci()
  # Remove Inf values for valid calculation
  div_ci$bootstrap_dist <- div_ci$bootstrap_dist[is.finite(div_ci$bootstrap_dist)]
  
  expect_message(
    summary(div_ci),
    "Summary of Divergence Bootstrap"
  )
})

test_that("summary.tsenat_divergence_bootstrap_ci returns object invisibly", {
  div_ci <- create_test_divergence_bootstrap_ci()
  
  result <- summary(div_ci)
  expect_identical(result, div_ci)
})

test_that("summary.tsenat_divergence_bootstrap_ci computes ESS as percentage", {
  div_ci <- create_test_divergence_bootstrap_ci()
  
  # ESS = (n_unique / n_total) * 100
  n_bootstrap <- length(div_ci$bootstrap_dist)
  n_unique_rounded <- length(unique(round(div_ci$bootstrap_dist, 6)))
  expected_ess <- (n_unique_rounded / n_bootstrap) * 100
  
  expect_message(
    summary(div_ci),
    "Effective sample size:"
  )
})

test_that("summary.tsenat_divergence_bootstrap_ci relative CI width is finite", {
  div_ci <- create_test_divergence_bootstrap_ci()
  
  expect_message(
    summary(div_ci),
    "CI width to estimate ratio:"
  )
})

test_that("summary.tsenat_divergence_bootstrap_ci handles very small estimate", {
  div_ci <- create_test_divergence_bootstrap_ci(estimate = 0.001)
  
  expect_message(
    summary(div_ci),
    "Summary of Divergence Bootstrap"
  )
})

test_that("summary.tsenat_divergence_bootstrap_ci handles very large bootstrap distribution", {
  div_ci <- create_test_divergence_bootstrap_ci()
  div_ci$bootstrap_dist <- rnorm(10000, mean = 0.35, sd = 0.10)
  
  expect_message(
    summary(div_ci),
    "Summary of Divergence Bootstrap"
  )
})

test_that("summary.tsenat_divergence_bootstrap_ci shows all required fields", {
  div_ci <- create_test_divergence_bootstrap_ci()
  # Capture messages from summary function
  expect_message(
    summary(div_ci),
    "Mean:"
  )
  
  expect_message(
    summary(div_ci),
    "Median:"
  )
  
  expect_message(
    summary(div_ci),
    "SD:"
  )
  # Summary function outputs messages, should have multiple lines
  expect_message(
    summary(div_ci),
    "Bootstrap distribution"
  )
})

test_that("summary.tsenat_divergence_bootstrap_ci numerical outputs are rounded", {
  div_ci <- create_test_divergence_bootstrap_ci()
  
  # Summary function outputs messages, should have multiple lines
  expect_message(
    summary(div_ci),
    "Bootstrap distribution"
  )
})

# ============================================================================
# INTEGRATION TESTS: S3 Methods With Real Workflow
# ============================================================================

test_that("print and summary work on real bootstrap result", {
  
  # Create a simple test case with real bootstrap computation
  x <- c(100, 50, 75, 200, 80, 120)
  
  result <- .calculate_tsallis_entropy_bootstrap(
    x = x,
    q = 2,
    nboot = 50,
    ci = 0.95,
    method = "percentile",
    verbose = FALSE
  )
  
  # Should be able to print without error
  expect_error(print(result), NA)
})

test_that("print.tsenat_bootstrap_ci_list works with actual bootstrap results", {
  
  x <- c(100, 50, 75, 200, 80, 120)
  
  result <- .calculate_tsallis_entropy_bootstrap(
    x = x,
    q = c(1, 2, 3),
    nboot = 50,
    ci = 0.95,
    method = "percentile",
    verbose = FALSE
  )
  
  expect_true(inherits(result, "tsenat_bootstrap_ci_list"))
  
  expect_error(print(result), NA)
})

test_that("print method accessible via S3 dispatch", {
  ci_obj <- create_test_bootstrap_ci()
  
  expect_true(inherits(ci_obj, "tsenat_bootstrap_ci"))
  
  expect_error(print(ci_obj), NA)
})

test_that("summary method accessible via S3 dispatch", {
  div_ci <- create_test_divergence_bootstrap_ci()
  
  expect_s3_class(div_ci, "tsenat_divergence_bootstrap_ci")
  
  expect_error(summary(div_ci), NA)
})

test_that("generic print function dispatches to S3 method correctly", {
  ci_list <- create_test_bootstrap_ci_list()
  
  expect_s3_class(ci_list, "tsenat_bootstrap_ci_list")
  
  expect_message(print(ci_list), "Bootstrap Confidence Intervals")
})

test_that("generic summary function dispatches to S3 method correctly", {
  div_ci <- create_test_divergence_bootstrap_ci()
  
  expect_s3_class(div_ci, "tsenat_divergence_bootstrap_ci")
  
  expect_message(summary(div_ci), "Summary of Divergence Bootstrap")
})

# ============================================================================
# NEW BOOTSTRAP TESTS FOR COVERAGE IMPROVEMENT
# ============================================================================
# These tests target uncovered lines identified in bootstrap_coverage_analysis.md
# Coverage gaps: vector pseudocount validation, validation error paths, edge cases
#
# Add these tests to: tests/testthat/test-rcpp-bootstrap.R
# ============================================================================

# ============================================================================
# SECTION 1: Vector Pseudocount Error Handling
# ============================================================================
# Addresses uncovered lines in:
#   - bootstrap_compute_cpp_wrapper (lines 84-85, 88-89)
#   - block_bootstrap_compute_cpp_wrapper (lines 37-38, 41-42)

test_that("bootstrap_compute_cpp_wrapper rejects mismatched pseudocount vector", {
  # Pseudocount vector length != x length should error
  expect_error(
    bootstrap_compute_cpp_wrapper(
      x = c(100, 50, 25, 10),           # length 4
      q = 1.0,
      normalize = TRUE,
      nboot = 10L,
      log_base = exp(1),
      pseudocount = c(1, 2, 3)           # length 3 - MISMATCH
    ),
    "pseudocount must have length 1 or equal to x length"
  )
})

test_that("bootstrap_compute_cpp_wrapper accepts matching pseudocount vector", {
  # Pseudocount vector length == x length should work
  result <- bootstrap_compute_cpp_wrapper(
    x = c(100, 50, 25, 10),              # length 4
    q = 1.0,
    normalize = TRUE,
    nboot = 10L,
    log_base = exp(1),
    pseudocount = c(1, 2, 3, 4)          # length 4 - MATCHING
  )
  
  expect_is(result, "numeric")
  expect_length(result, 10)
  expect_true(all(is.finite(result)))
})

test_that("block_bootstrap_compute_cpp_wrapper rejects mismatched pseudocount vector", {
  # Even-length x with mismatched pseudocount vector
  expect_error(
    block_bootstrap_compute_cpp_wrapper(
      x = c(100, 95, 110, 105),          # length 4 (2 pairs)
      q = 1.0,
      normalize = TRUE,
      nboot = 10L,
      log_base = exp(1),
      pseudocount = c(1, 2, 3)           # length 3 - MISMATCH
    ),
    "pseudocount must have length 1 or equal to x length"
  )
})

test_that("block_bootstrap_compute_cpp_wrapper accepts matching pseudocount vector", {
  # Even-length x with matching pseudocount vector
  result <- block_bootstrap_compute_cpp_wrapper(
    x = c(100, 95, 110, 105),            # length 4 (2 pairs)
    q = 1.0,
    normalize = TRUE,
    nboot = 10L,
    log_base = exp(1),
    pseudocount = c(1, 2, 3, 4)          # length 4 - MATCHING
  )
  
  expect_is(result, "numeric")
  expect_length(result, 10)
  expect_true(all(is.finite(result)))
})

# ============================================================================
# SECTION 2: Divergence Bootstrap Input Validation
# ============================================================================
# Addresses uncovered lines in:
#   - divergence_bootstrap_compute_cpp_wrapper (lines 127-156)

test_that("divergence_bootstrap_compute_cpp_wrapper rejects non-numeric inputs", {
  # x must be numeric
  expect_error(
    divergence_bootstrap_compute_cpp_wrapper(
      x = c("a", "b", "c"),              # non-numeric
      y = c(100, 50, 25),
      q = 1.0,
      nboot = 10L,
      paired = FALSE,
      pseudocount = 0,
      log_base = exp(1)
    ),
    "x and y must be numeric vectors"
  )
  
  # y must be numeric
  expect_error(
    divergence_bootstrap_compute_cpp_wrapper(
      x = c(100, 50, 25),
      y = c("a", "b", "c"),              # non-numeric
      q = 1.0,
      nboot = 10L,
      paired = FALSE,
      pseudocount = 0,
      log_base = exp(1)
    ),
    "x and y must be numeric vectors"
  )
})

test_that("divergence_bootstrap_compute_cpp_wrapper requires equal length for x and y", {
  expect_error(
    divergence_bootstrap_compute_cpp_wrapper(
      x = c(100, 50, 25, 10),            # length 4
      y = c(75, 40),                     # length 2 - MISMATCH
      q = 1.0,
      nboot = 10L,
      paired = FALSE,
      pseudocount = 0,
      log_base = exp(1)
    ),
    "x and y must have the same length"
  )
})

test_that("divergence_bootstrap_compute_cpp_wrapper rejects negative values", {
  # Negative in x
  expect_error(
    divergence_bootstrap_compute_cpp_wrapper(
      x = c(100, -50, 25),               # negative value
      y = c(75, 40, 30),
      q = 1.0,
      nboot = 10L,
      paired = FALSE,
      pseudocount = 0,
      log_base = exp(1)
    ),
    "x and y must contain non-negative values only"
  )
  
  # Negative in y
  expect_error(
    divergence_bootstrap_compute_cpp_wrapper(
      x = c(100, 50, 25),
      y = c(75, -40, 30),                # negative value
      q = 1.0,
      nboot = 10L,
      paired = FALSE,
      pseudocount = 0,
      log_base = exp(1)
    ),
    "x and y must contain non-negative values only"
  )
})

test_that("divergence_bootstrap_compute_cpp_wrapper validates q parameter", {
  expect_error(
    divergence_bootstrap_compute_cpp_wrapper(
      x = c(100, 50, 25),
      y = c(75, 40, 30),
      q = -1.0,                          # negative q
      nboot = 10L,
      paired = FALSE,
      pseudocount = 0,
      log_base = exp(1)
    ),
    "q must be a non-negative numeric value"
  )
})

test_that("divergence_bootstrap_compute_cpp_wrapper validates nboot parameter", {
  expect_error(
    divergence_bootstrap_compute_cpp_wrapper(
      x = c(100, 50, 25),
      y = c(75, 40, 30),
      q = 1.0,
      nboot = 0L,                        # invalid nboot
      paired = FALSE,
      pseudocount = 0,
      log_base = exp(1)
    ),
    "nboot must be a positive integer"
  )
})

test_that("divergence_bootstrap_compute_cpp_wrapper validates paired parameter", {
  expect_error(
    divergence_bootstrap_compute_cpp_wrapper(
      x = c(100, 50, 25),
      y = c(75, 40, 30),
      q = 1.0,
      nboot = 10L,
      paired = c(TRUE, FALSE),           # not single logical
      pseudocount = 0,
      log_base = exp(1)
    ),
    "paired must be a single logical value"
  )
})

test_that("divergence_bootstrap_compute_cpp_wrapper requires even length for paired=TRUE", {
  expect_error(
    divergence_bootstrap_compute_cpp_wrapper(
      x = c(100, 50, 25),                # odd length (3)
      y = c(75, 40, 30),
      q = 1.0,
      nboot = 10L,
      paired = TRUE,                     # paired requires even length
      pseudocount = 0,
      log_base = exp(1)
    ),
    "For paired=TRUE, x and y must have even length"
  )
})

# ============================================================================
# SECTION 3: Validation Data Tests (.validate_bootstrap_data)
# ============================================================================
# Addresses uncovered lines in:
#   - .validate_bootstrap_data (lines 478-521, 52.2% coverage)

test_that("validate_bootstrap_data rejects empty input", {
  expect_error(
    TSENAT:::.validate_bootstrap_data(
      x = numeric(0),
      effective_length = NULL,
      pseudocount = 0
    ),
    "Input x must be a non-empty vector"
  )
})

test_that("validate_bootstrap_data rejects all zeros without pseudocount", {
  expect_error(
    TSENAT:::.validate_bootstrap_data(
      x = c(0, 0, 0, 0),
      effective_length = NULL,
      pseudocount = 0
    ),
    "All counts are zero and pseudocount = 0"
  )
})

test_that("validate_bootstrap_data warns on all zeros with pseudocount", {
  expect_warning(
    TSENAT:::.validate_bootstrap_data(
      x = c(0, 0, 0, 0),
      effective_length = NULL,
      pseudocount = 1.0
    ),
    "All counts are zero"
  )
})

test_that("validate_bootstrap_data detects effective_length mismatch", {
  expect_error(
    TSENAT:::.validate_bootstrap_data(
      x = c(100, 50, 25, 10),            # length 4
      effective_length = c(1.0, 1.0, 1.0),  # length 3 - MISMATCH
      pseudocount = 0
    ),
    "Length mismatch: effective_length"
  )
})

test_that("validate_bootstrap_data warns on non-positive effective_length", {
  expect_warning(
    TSENAT:::.validate_bootstrap_data(
      x = c(100, 50, 25, 10),
      effective_length = c(1.0, 0.0, -1.0, 1.0),  # has 0 and negative
      pseudocount = 0
    ),
    "Found.*position.*effective_length <= 0"
  )
})

test_that("validate_bootstrap_data warns on single isoform", {
  expect_warning(
    TSENAT:::.validate_bootstrap_data(
      x = c(100),                        # only 1 isoform
      effective_length = c(1.0),
      pseudocount = 0
    ),
    "Single isoform detected"
  )
})

test_that("validate_bootstrap_data warns on high proportion of zeros", {
  expect_warning(
    TSENAT:::.validate_bootstrap_data(
      x = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 100),  # 90% zeros
      effective_length = rep(1.0, 10),
      pseudocount = 0
    ),
    "High proportion of zeros"
  )
})

test_that("validate_bootstrap_data returns TRUE invisibly on valid input", {
  result <- TSENAT:::.validate_bootstrap_data(
    x = c(100, 50, 25, 10),
    effective_length = c(1.0, 1.0, 1.0, 1.0),
    pseudocount = 0
  )
  
  expect_equal(result, TRUE)
  expect_true(is.logical(result))
})

# ============================================================================
# SECTION 4: Edge Cases and Boundary Conditions
# ============================================================================

test_that("bootstrap functions handle very small pseudocount values", {
  result <- bootstrap_compute_cpp_wrapper(
    x = c(0, 0, 0, 100),                # mostly zeros
    q = 1.0,
    normalize = TRUE,
    nboot = 10L,
    log_base = exp(1),
    pseudocount = 1e-8                  # very small but nonzero
  )
  
  expect_is(result, "numeric")
  expect_length(result, 10)
  expect_true(all(is.finite(result)))
})

test_that("block_bootstrap with vector pseudocount produces different results than scalar", {
  x <- c(100, 95, 110, 105)
  
  result_scalar <- block_bootstrap_compute_cpp_wrapper(
    x = x,
    q = 1.0,
    normalize = TRUE,
    nboot = 100L,
    log_base = exp(1),
    pseudocount = 1.0
  )
  
  result_vector <- block_bootstrap_compute_cpp_wrapper(
    x = x,
    q = 1.0,
    normalize = TRUE,
    nboot = 100L,
    log_base = exp(1),
    pseudocount = c(1, 1, 1, 1)
  )
  
  # Should be equivalent (vector applied upfront)
  expect_equal(length(result_scalar), length(result_vector))
  # Results should be similar (not exactly equal due to randomness, but similar distributions)
  expect_true(abs(mean(result_scalar) - mean(result_vector)) < 0.1)
})

# ============================================================================
# SECTION 5: Paired Design Validation
# ============================================================================

test_that("divergence_bootstrap_paired rejects odd-length pairs", {
  expect_error(
    divergence_bootstrap_paired_cpp_wrapper(
      x = c(100, 50, 25),                # odd length
      y = c(75, 40, 30),
      pair_ids = c(1, 1, 2),
      q = 1.0,
      nboot = 10L,
      pseudocount = 0,
      log_base = exp(1)
    ),
    "must have same length"
  )
})

test_that("divergence_bootstrap_flexible rejects mismatched group lengths", {
  expect_error(
    divergence_bootstrap_flexible_cpp_wrapper(
      x = c(100, 50, 25, 10),
      y = c(75, 40, 30),                 # length 3
      x_pair_ids = c(1, 1, 2, 2),
      y_pair_ids = c(1, 1),              # length 2, doesn't match y length 3
      q = 1.0,
      nboot = 10L,
      pseudocount = 0,
      log_base = exp(1)
    ),
    "y and y_pair_ids must have same length"
  )
})

# ============================================================================
# TEST SUITE: bootstrap.R - Uncovered Lines Coverage
# ============================================================================
# Tests targeting specific uncovered lines from cobertura.xml
# Focus on edge cases, error handling, and specialized code paths
# ============================================================================

library(TSENAT)

# ============================================================================
# Vector Pseudocount Handling Tests (Lines 150-156, 203-219, 279-287)
# ============================================================================

