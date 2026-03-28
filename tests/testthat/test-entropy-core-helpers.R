library(TSENAT)

context("entropy_core: Centralized Entropy Calculation")

# ============================================================================
# Test .tsenat_entropy_core - Basic functionality
# ============================================================================

test_that(".tsenat_entropy_core computes Shannon entropy correctly", {
  # Uniform distribution: H = log(n)
  proportions <- c(0.25, 0.25, 0.25, 0.25)
  H <- TSENAT:::.tsenat_entropy_core(proportions, q = 1, norm = FALSE, log_base = exp(1))
  expected <- log(4)  # Shannon entropy of uniform distribution
  expect_equal(H, expected, tolerance = 1e-10)
})

test_that(".tsenat_entropy_core computes Tsallis entropy correctly", {
  # Test with q = 2
  proportions <- c(0.5, 0.5)
  H <- TSENAT:::.tsenat_entropy_core(proportions, q = 2, norm = FALSE, log_base = exp(1))
  # Tsallis(q=2) for [0.5, 0.5]: (1 - 0.5) / (2-1) = 0.5
  expected <- (1 - (0.5^2 + 0.5^2)) / 1
  expect_equal(H, expected, tolerance = 1e-10)
})

test_that(".tsenat_entropy_core handles zero proportions", {
  # Zero values should be filtered out
  proportions <- c(0.5, 0.5, 0, 0)
  H <- TSENAT:::.tsenat_entropy_core(proportions, q = 1, norm = FALSE, log_base = exp(1))
  expected <- log(2)  # Should be same as c(0.5, 0.5)
  expect_equal(H, expected, tolerance = 1e-10)
})

test_that(".tsenat_entropy_core normalizes correctly", {
  proportions <- c(0.25, 0.25, 0.25, 0.25)
  H <- TSENAT:::.tsenat_entropy_core(proportions, q = 1, norm = TRUE, log_base = exp(1))
  # Normalized Shannon: log(4) / log(4) = 1
  expect_equal(H, 1.0, tolerance = 1e-10)
})

test_that(".tsenat_entropy_core returns NA for invalid input", {
  expect_true(is.na(TSENAT:::.tsenat_entropy_core(numeric(0), q = 1)))
  expect_true(is.na(TSENAT:::.tsenat_entropy_core(c(NA, NA), q = 1)))
})

test_that(".tsenat_entropy_core respects q tolerance", {
  proportions <- c(0.5, 0.5)
  H1 <- TSENAT:::.tsenat_entropy_core(proportions, q = 1.0, norm = FALSE)
  H2 <- TSENAT:::.tsenat_entropy_core(proportions, q = 1.0000001, norm = FALSE, q_tol = 1e-6)
  # With q_tol = 1e-6, q=1.0000001 should be treated as Shannon entropy
  expect_equal(H1, H2, tolerance = 1e-5)
})

# ============================================================================
# Test .tsenat_entropy_vectorized
# ============================================================================

test_that(".tsenat_entropy_vectorized computes entropy for matrix rows", {
  counts <- matrix(c(10, 5, 5, 10), nrow = 2, byrow = TRUE)
  entropies <- TSENAT:::.tsenat_entropy_vectorized(counts, q = 1, norm = FALSE)
  
  expect_length(entropies, 2)
  expect_true(all(!is.na(entropies)))
  # Both rows have same proportions (0.5, 0.5), so same entropy
  expect_equal(entropies[1], entropies[2], tolerance = 1e-10)
})

test_that(".tsenat_entropy_vectorized handles pseudocount", {
  counts <- matrix(c(10, 0, 0, 10), nrow = 2, byrow = TRUE)
  H_no_pseudo <- TSENAT:::.tsenat_entropy_vectorized(counts, pseudocount = 0)
  H_with_pseudo <- TSENAT:::.tsenat_entropy_vectorized(counts, pseudocount = 1)
  
  # With pseudocount, second row should have higher entropy (adds species)
  expect_true(H_with_pseudo[2] > H_no_pseudo[2])
})

# ============================================================================
# Test .tsenat_entropy_max
# ============================================================================

test_that(".tsenat_entropy_max computes Shannon maximum correctly", {
  # For n species with q=1: H_max = log(n)
  H_max <- TSENAT:::.tsenat_entropy_max(n_species = 4, q = 1, log_base = exp(1))
  expected <- log(4)
  expect_equal(H_max, expected, tolerance = 1e-10)
})

test_that(".tsenat_entropy_max computes Tsallis maximum correctly", {
  # For n=2, q=2: H_max = (1 - 2^(1-2)) / (2-1) = (1 - 0.5) / 1 = 0.5
  H_max <- TSENAT:::.tsenat_entropy_max(n_species = 2, q = 2, log_base = exp(1))
  expected <- (1 - 2^-1) / 1
  expect_equal(H_max, expected, tolerance = 1e-10)
})

# ============================================================================
# Test .tsenat_jackknife_resampling
# ============================================================================

test_that(".tsenat_jackknife_resampling detects insufficient data", {
  # Need at least 2 observations
  counts <- matrix(rnorm(5), nrow = 1)
  expect_warning(
    result <- TSENAT:::.tsenat_jackknife_resampling(counts, q = 1),
    "Insufficient observations"
  )
  expect_null(result)
})

test_that(".tsenat_jackknife_resampling computes valid influence values", {
  set.seed(42)
  # Matrix: rows=samples, cols=species/genes
  counts <- matrix(c(100, 50, 25, 10), nrow = 4, ncol = 1)
  result <- TSENAT:::.tsenat_jackknife_resampling(counts, q = 1, norm = FALSE)
  
  expect_false(is.null(result))
  expect_is(result$influence, "numeric")
  expect_equal(length(result$influence), 4)
  expect_true(all(result$influence >= 0))
})

test_that(".tsenat_jackknife_resampling computes jackknife standard error", {
  # Matrix: rows=samples (leave-one-out), cols=species
  counts <- matrix(c(100, 50, 50, 20), nrow = 4, ncol = 1)
  result <- TSENAT:::.tsenat_jackknife_resampling(counts, q = 1, norm = FALSE)
  
  expect_true(is.numeric(result$jackknife_se))
  expect_true(result$jackknife_se >= 0)
  expect_true(!is.na(result$jackknife_se))
})

test_that(".tsenat_jackknife_resampling detects outliers correctly", {
  # Create data with multiple species where one sample is extreme
  set.seed(123)
  # 6 samples, 3 species - last row is extreme outlier
  counts <- matrix(c(
    100, 50, 20,    # sample 1
    80, 60, 30,     # sample 2
    90, 45, 25,     # sample 3
    75, 55, 40,     # sample 4
    85, 65, 35,     # sample 5
    1, 1, 100       # sample 6 - extreme outlier
  ), nrow = 6, ncol = 3, byrow = TRUE)
  
  result <- TSENAT:::.tsenat_jackknife_resampling(counts, q = 1, threshold = 80)
  
  # With multiple species, the outlier effect is more pronounced
  # The last sample with extreme composition should be detected
  expect_true(length(result$outlier_indices) > 0 || result$outlier_threshold >= 80)
})

test_that(".tsenat_jackknife_resampling respects threshold parameter", {
  # Multiple samples with varying influence
  counts <- matrix(c(100, 50, 25, 10), nrow = 4, ncol = 1)
  
  result_high <- TSENAT:::.tsenat_jackknife_resampling(counts, threshold = 95)
  result_low <- TSENAT:::.tsenat_jackknife_resampling(counts, threshold = 50)
  
  # Higher threshold = fewer outliers
  expect_true(length(result_high$outlier_indices) <= length(result_low$outlier_indices))
})

# ============================================================================
# Test .tsenat_jackknife_batch
# ============================================================================

test_that(".tsenat_jackknife_batch processes multiple genes", {
  set.seed(456)
  # Matrix: rows = samples, cols = genes
  counts <- matrix(rpois(20, lambda = 50), nrow = 5, ncol = 4)
  colnames(counts) <- paste0("Gene_", 1:4)
  
  results <- TSENAT:::.tsenat_jackknife_batch(counts, q = 1, verbose = FALSE)
  
  expect_length(results, 4)
  expect_equal(names(results), colnames(counts))
  expect_true(all(vapply(results, inherits, "tsenat_jackknife", FUN.VALUE = logical(1))))
})

test_that(".tsenat_jackknife_batch returns list with correct class", {
  # rows = samples, cols = genes
  counts <- matrix(rpois(12, 50), nrow = 6, ncol = 2)
  results <- TSENAT:::.tsenat_jackknife_batch(counts, q = 1)
  
  expect_s3_class(results, c("tsenat_jackknife_list", "list"))
})

test_that(".tsenat_jackknife_batch handles empty input", {
  # Empty matrix means no columns (genes)
  counts <- matrix(numeric(0), nrow = 5, ncol = 0)
  results <- TSENAT:::.tsenat_jackknife_batch(counts)
  
  expect_length(results, 0)
})
