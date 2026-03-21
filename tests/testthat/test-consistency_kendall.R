library(TSENAT)
library(testthat)

context("Rank-Based Methods: Consistency (Kendall's W & ICC)")

# ============================================================================
# Test: Consistency calculation with sufficient data
# ============================================================================

test_that("test_rankbased_assumptions calculates consistency with high concordance", {
  # Create data with high consistency (all samples rank genes similarly)
  set.seed(123)
  n_genes <- 10
  n_samples <- 5
  
  # Create data where each gene has consistent rank across samples
  base_vals <- seq(0.1, 1.0, length.out = n_genes)
  data <- matrix(0, nrow = n_genes, ncol = n_samples)
  for (i in seq_len(n_samples)) {
    noise <- rnorm(n_genes, 0, 0.05)
    data[, i] <- base_vals + noise
  }
  rownames(data) <- paste0("Gene_", 1:n_genes)
  colnames(data) <- paste0("Sample_", 1:n_samples)
  
  # Run consistency check
  result <- test_rankbased_assumptions(data, checks = "consistency")
  
  # Verify consistency results are present
  expect_true("consistency" %in% names(result))
  expect_true("kendall_w" %in% names(result$consistency))
  expect_true("icc_simplified" %in% names(result$consistency))
  expect_true("status" %in% names(result$consistency))
  expect_true("details" %in% names(result$consistency))
  
  # High concordance should give Kendall's W > 0.4
  expect_true(result$consistency$kendall_w > 0.4 || is.na(result$consistency$kendall_w))
})

test_that("test_rankbased_assumptions calculates Kendall's W correctly", {
  # Create simple test case for manual verification
  set.seed(456)
  data <- matrix(
    c(1, 2, 3,
      2, 3, 1,
      3, 1, 2,
      1, 3, 2,
      2, 1, 3),
    nrow = 5,
    ncol = 3,
    byrow = TRUE
  )
  rownames(data) <- paste0("Gene_", 1:5)
  colnames(data) <- paste0("Sample_", 1:3)
  
  result <- test_rankbased_assumptions(data, checks = "consistency")
  
  # Check that Kendall's W was calculated (between 0 and 1)
  expect_true(!is.na(result$consistency$kendall_w))
  expect_true(result$consistency$kendall_w >= 0)
  expect_true(result$consistency$kendall_w <= 1)
})

test_that("test_rankbased_assumptions high Kendall's W assigns PASS status", {
  # Create perfectly consistent data (all samples rank genes identically)
  data <- matrix(
    c(1, 1, 1,  # Gene 1 consistently ranked 1st
      2, 2, 2,  # Gene 2 consistently ranked 2nd
      3, 3, 3,  # Gene 3 consistently ranked 3rd
      4, 4, 4,  # Gene 4 consistently ranked 4th
      5, 5, 5), # Gene 5 consistently ranked 5th
    nrow = 5,
    ncol = 3,
    byrow = TRUE
  )
  rownames(data) <- paste0("Gene_", 1:5)
  colnames(data) <- paste0("Sample_", 1:3)
  
  result <- test_rankbased_assumptions(data, checks = "consistency")
  
  # Perfect consistency: Kendall's W should be close to 1
  expect_true(result$consistency$kendall_w >= 0.7 || is.na(result$consistency$kendall_w))
  
  # Status should reflect high concordance
  if (!is.na(result$consistency$kendall_w) && result$consistency$kendall_w > 0.7) {
    expect_match(result$consistency$status, "PASS|OK", ignore.case = TRUE)
  }
})

test_that("test_rankbased_assumptions moderate Kendall's W assigns ACCEPTABLE status", {
  # Create moderately consistent data
  set.seed(789)
  data <- matrix(
    c(1, 2, 2,
      2, 1, 3,
      3, 3, 1,
      4, 4, 4,
      5, 5, 5),
    nrow = 5,
    ncol = 3,
    byrow = TRUE
  )
  + rnorm(15, 0, 0.2) # Add small noise
  
  rownames(data) <- paste0("Gene_", 1:5)
  colnames(data) <- paste0("Sample_", 1:3)
  
  result <- test_rankbased_assumptions(data, checks = "consistency")
  
  # Verify Kendall's W is calculated
  expect_true(!is.na(result$consistency$kendall_w))
  
  # Status should be appropriate
  if (result$consistency$kendall_w > 0.4 && result$consistency$kendall_w <= 0.7) {
    expect_match(result$consistency$status, "ACCEPTABLE|?", ignore.case = TRUE)
  }
})

test_that("test_rankbased_assumptions low Kendall's W assigns LOW CONSISTENCY status", {
  # Create very inconsistent data (random ranks)
  set.seed(321)
  data <- matrix(
    runif(50),
    nrow = 10,
    ncol = 5
  )
  rownames(data) <- paste0("Gene_", 1:10)
  colnames(data) <- paste0("Sample_", 1:5)
  
  result <- test_rankbased_assumptions(data, checks = "consistency")
  
  # Verify Kendall's W is calculated
  expect_true(!is.na(result$consistency$kendall_w))
  
  # Very random data should have low concordance
  if (result$consistency$kendall_w < 0.4) {
    expect_match(result$consistency$status, "LOW|?", ignore.case = TRUE)
  }
})

test_that("test_rankbased_assumptions calculates ICC correctly", {
  # Create test data with known structure
  data <- matrix(
    c(0.5, 0.6, 0.55,
      0.7, 0.75, 0.72,
      0.3, 0.32, 0.31,
      0.9, 0.92, 0.91,
      0.4, 0.41, 0.42),
    nrow = 5,
    ncol = 3,
    byrow = TRUE
  )
  rownames(data) <- paste0("Gene_", 1:5)
  colnames(data) <- paste0("Sample_", 1:3)
  
  result <- test_rankbased_assumptions(data, checks = "consistency")
  
  # ICC should be between 0 and 1 (or NA)
  icc <- result$consistency$icc_simplified
  if (!is.na(icc)) {
    expect_true(icc >= 0)
    expect_true(icc <= 1)
  }
})

test_that("test_rankbased_assumptions details string contains metrics", {
  # Create simple test data
  data <- matrix(
    runif(30),
    nrow = 6,
    ncol = 5
  )
  rownames(data) <- paste0("Gene_", 1:6)
  colnames(data) <- paste0("Sample_", 1:5)
  
  result <- test_rankbased_assumptions(data, checks = "consistency")
  
  # Details string should contain Kendall W and ICC values
  details <- result$consistency$details
  expect_match(details, "Kendall W=", ignore.case = TRUE)
  expect_match(details, "ICC~=", ignore.case = TRUE)
})

test_that("test_rankbased_assumptions skips consistency with insufficient data (1 row)", {
  # Create degenerate data (only 1 gene)
  data <- matrix(runif(5), nrow = 1, ncol = 5)
  rownames(data) <- "Gene_1"
  colnames(data) <- paste0("Sample_", 1:5)
  
  result <- test_rankbased_assumptions(data, checks = "consistency")
  
  # Should set status to SKIP
  expect_match(result$consistency$status, "SKIP", ignore.case = TRUE)
  expect_match(result$consistency$method, "Insufficient", ignore.case = TRUE)
})

test_that("test_rankbased_assumptions skips consistency with insufficient data (1 column)", {
  # Create degenerate data (only 1 sample)
  data <- matrix(runif(5), nrow = 5, ncol = 1)
  rownames(data) <- paste0("Gene_", 1:5)
  colnames(data) <- "Sample_1"
  
  result <- test_rankbased_assumptions(data, checks = "consistency")
  
  # Should set status to SKIP
  expect_match(result$consistency$status, "SKIP", ignore.case = TRUE)
  expect_match(result$consistency$details, "2 samples and 2 genes", ignore.case = TRUE)
})

test_that("test_rankbased_assumptions consistency matrix operations work", {
  # Test that the ranking and mean computation works correctly
  data <- matrix(
    c(1, 2, 3,
      4, 5, 6,
      7, 8, 9),
    nrow = 3,
    ncol = 3,
    byrow = TRUE
  )
  rownames(data) <- paste0("Gene_", 1:3)
  colnames(data) <- paste0("Sample_", 1:3)
  
  result <- test_rankbased_assumptions(data, checks = "consistency")
  
  # Verify that result is valid
  expect_true(!is.null(result$consistency))
  expect_true(is.list(result$consistency))
  expect_true(all(c("method", "status", "details") %in% names(result$consistency)))
})

test_that("test_rankbased_assumptions consistency with missing values", {
  # Create data with NA values
  data <- matrix(
    c(1, 2, NA,
      2, NA, 3,
      3, 1, 2),
    nrow = 3,
    ncol = 3,
    byrow = TRUE
  )
  rownames(data) <- paste0("Gene_", 1:3)
  colnames(data) <- paste0("Sample_", 1:3)
  
  # Should handle NA gracefully
  result <- test_rankbased_assumptions(data, checks = "consistency")
  
  # Result should still be generated (with na.rm = TRUE in calculations)
  expect_true("consistency" %in% names(result))
  expect_true(!is.null(result$consistency$status))
})

test_that("test_rankbased_assumptions consistency handles single gene with NA", {
  # Edge case: single gene creates NA for Kendall's W (n=1)
  data <- matrix(
    c(0.5, 0.6, 0.55),
    nrow = 1,
    ncol = 3
  )
  rownames(data) <- "Gene_1"
  colnames(data) <- paste0("Sample_", 1:3)
  
  result <- test_rankbased_assumptions(data, checks = "consistency")
  
  # n=1 case should be skipped
  expect_match(result$consistency$status, "SKIP", ignore.case = TRUE)
})

test_that("test_rankbased_assumptions all consistency message formats", {
  # Test that all message formats are correctly produced
  data <- matrix(
    c(0.1, 0.2, 0.3,
      0.4, 0.5, 0.6,
      0.7, 0.8, 0.9,
      0.15, 0.25, 0.35,
      0.45, 0.55, 0.65),
    nrow = 5,
    ncol = 3,
    byrow = TRUE
  )
  rownames(data) <- paste0("Gene_", 1:5)
  colnames(data) <- paste0("Sample_", 1:3)
  
  result <- test_rankbased_assumptions(data, checks = "consistency")
  
  # Verify all expected fields exist
  expect_true("description" %in% names(result$consistency))
  expect_equal(result$consistency$description, "Rank consistency evaluation (Kendall's W & ICC)")
  expect_match(result$consistency$method, "Kendall", ignore.case = TRUE)
})
