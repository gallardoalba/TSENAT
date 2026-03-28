library(TSENAT)

context("entropy_core: Centralized Entropy Calculation")

# ============================================================================
# Test .entropy_core - Basic functionality
# ============================================================================

test_that(".entropy_core computes Shannon entropy correctly", {
  # Uniform distribution: H = log(n)
  proportions <- c(0.25, 0.25, 0.25, 0.25)
  H <- TSENAT:::.entropy_core(proportions, q = 1, norm = FALSE, log_base = exp(1))
  expected <- log(4)  # Shannon entropy of uniform distribution
  expect_equal(H, expected, tolerance = 1e-10)
})

test_that(".entropy_core computes Tsallis entropy correctly", {
  # Test with q = 2
  proportions <- c(0.5, 0.5)
  H <- TSENAT:::.entropy_core(proportions, q = 2, norm = FALSE, log_base = exp(1))
  # Tsallis(q=2) for [0.5, 0.5]: (1 - 0.5) / (2-1) = 0.5
  expected <- (1 - (0.5^2 + 0.5^2)) / 1
  expect_equal(H, expected, tolerance = 1e-10)
})

test_that(".entropy_core handles zero proportions", {
  # Zero values should be filtered out
  proportions <- c(0.5, 0.5, 0, 0)
  H <- TSENAT:::.entropy_core(proportions, q = 1, norm = FALSE, log_base = exp(1))
  expected <- log(2)  # Should be same as c(0.5, 0.5)
  expect_equal(H, expected, tolerance = 1e-10)
})

test_that(".entropy_core normalizes correctly", {
  proportions <- c(0.25, 0.25, 0.25, 0.25)
  H <- TSENAT:::.entropy_core(proportions, q = 1, norm = TRUE, log_base = exp(1))
  # Normalized Shannon: log(4) / log(4) = 1
  expect_equal(H, 1.0, tolerance = 1e-10)
})

test_that(".entropy_core returns NA for invalid input", {
  expect_true(is.na(TSENAT:::.entropy_core(numeric(0), q = 1)))
  expect_true(is.na(TSENAT:::.entropy_core(c(NA, NA), q = 1)))
})

test_that(".entropy_core respects q tolerance", {
  proportions <- c(0.5, 0.5)
  H1 <- TSENAT:::.entropy_core(proportions, q = 1.0, norm = FALSE)
  H2 <- TSENAT:::.entropy_core(proportions, q = 1.0000001, norm = FALSE, q_tol = 1e-6)
  # With q_tol = 1e-6, q=1.0000001 should be treated as Shannon entropy
  expect_equal(H1, H2, tolerance = 1e-5)
})

# ============================================================================
# Test .entropy_vectorized
# ============================================================================

test_that(".entropy_vectorized computes entropy for matrix rows", {
  counts <- matrix(c(10, 5, 5, 10), nrow = 2, byrow = TRUE)
  entropies <- TSENAT:::.entropy_vectorized(counts, q = 1, norm = FALSE)
  
  expect_length(entropies, 2)
  expect_true(all(!is.na(entropies)))
  # Both rows have same proportions (0.5, 0.5), so same entropy
  expect_equal(entropies[1], entropies[2], tolerance = 1e-10)
})

test_that(".entropy_vectorized handles pseudocount", {
  counts <- matrix(c(10, 0, 0, 10), nrow = 2, byrow = TRUE)
  H_no_pseudo <- TSENAT:::.entropy_vectorized(counts, pseudocount = 0)
  H_with_pseudo <- TSENAT:::.entropy_vectorized(counts, pseudocount = 1)
  
  # With pseudocount, second row should have higher entropy (adds species)
  expect_true(H_with_pseudo[2] > H_no_pseudo[2])
})

# ============================================================================
# Test .entropy_max
# ============================================================================

test_that(".entropy_max computes Shannon maximum correctly", {
  # For n species with q=1: H_max = log(n)
  H_max <- TSENAT:::.entropy_max(n_species = 4, q = 1, log_base = exp(1))
  expected <- log(4)
  expect_equal(H_max, expected, tolerance = 1e-10)
})

test_that(".entropy_max computes Tsallis maximum correctly", {
  # For n=2, q=2: H_max = (1 - 2^(1-2)) / (2-1) = (1 - 0.5) / 1 = 0.5
  H_max <- TSENAT:::.entropy_max(n_species = 2, q = 2, log_base = exp(1))
  expected <- (1 - 2^-1) / 1
  expect_equal(H_max, expected, tolerance = 1e-10)
})

# ============================================================================
# NOTE: .jackknife_resampling tests removed
# This function has been deprecated in favor of the C++ implementation
# The R version has been completely replaced with the optimized C++/Rcpp version
# See test-rcpp-jackknife.R for comprehensive C++ jackknife tests
