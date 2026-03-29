# Tests for .estimate_pseudocount() function
# Library size normalization for pseudocount estimation

library(testthat)
library(TSENAT)

context("Pseudocount Estimation: Library Size Normalization")

test_that("estimate_pseudocount works with matrix input", {
  counts <- matrix(c(10, 5, 1, 20, 8, 3, 15, 10, 5), nrow = 3, ncol = 3)
  rownames(counts) <- c("Gene1", "Gene2", "Gene3")
  
  result <- .estimate_pseudocount(counts, verbose = FALSE)
  
  # Check return structure
  expect_is(result, "list")
  expect_true("scalar_pseudocount" %in% names(result))
  expect_true("size_factors" %in% names(result))
  expect_true("diagnostics" %in% names(result))
  
  # Check scalar pseudocount
  expect_is(result$scalar_pseudocount, "numeric")
  expect_true(result$scalar_pseudocount > 0)
  expect_true(is.finite(result$scalar_pseudocount))
  
  # Check size factors
  expect_is(result$size_factors, "numeric")
  expect_equal(length(result$size_factors), 3)
})

test_that("estimate_pseudocount works with SummarizedExperiment input", {
  skip_if_not_installed("SummarizedExperiment")
  
  counts <- matrix(c(10, 5, 1, 20, 8, 3, 15, 10, 5), nrow = 3, ncol = 3)
  rownames(counts) <- c("Gene1", "Gene2", "Gene3")
  se <- suppressWarnings(SummarizedExperiment::SummarizedExperiment(assay = list(data = counts)))
  
  result <- .estimate_pseudocount(se, verbose = FALSE)
  
  # Check return structure
  expect_is(result, "list")
  expect_equal(length(result$size_factors), 3)
  expect_true(result$scalar_pseudocount > 0)
})

test_that("estimate_pseudocount diagnostics are valid", {
  counts <- matrix(c(10, 5, 1, 20, 8, 3, 15, 10, 5), nrow = 3, ncol = 3)
  
  result <- .estimate_pseudocount(counts, verbose = FALSE)
  
  # Check diagnostics structure
  expect_is(result$diagnostics, "list")
  expect_true("n_genes" %in% names(result$diagnostics))
  expect_true("n_samples" %in% names(result$diagnostics))
  expect_true("mean_lib_size" %in% names(result$diagnostics))
  
  # Check values
  expect_equal(result$diagnostics$n_genes, 3)
  expect_equal(result$diagnostics$n_samples, 3)
  expect_true(result$diagnostics$mean_lib_size > 0)
})

test_that("estimate_pseudocount size factors are normalized", {
  counts <- matrix(c(10, 5, 1, 20, 8, 3, 15, 10, 5), nrow = 3, ncol = 3)
  rownames(counts) <- c("Gene1", "Gene2", "Gene3")
  
  result <- .estimate_pseudocount(counts, verbose = FALSE)
  
  # Size factors should normalize around 1
  expect_equal(mean(result$size_factors), 1, tolerance = 1e-6)
  expect_true(all(result$size_factors > 0))
  expect_true(all(is.finite(result$size_factors)))
})

test_that("estimate_pseudocount returns reasonable pseudocount values", {
  counts <- matrix(c(10, 5, 1, 20, 8, 3, 15, 10, 5), nrow = 3, ncol = 3)
  
  result <- .estimate_pseudocount(counts, verbose = FALSE)
  
  # Pseudocount should be reasonably small and positive
  expect_true(result$scalar_pseudocount > 0)
  expect_true(result$scalar_pseudocount < 20)
  
  # Should be finite
  expect_true(is.finite(result$scalar_pseudocount))
})

test_that("estimate_pseudocount handles sparse counts", {
  # Matrix with mostly zeros, but all samples have at least some counts
  counts <- matrix(0, nrow = 5, ncol = 10)
  counts[1, 1:3] <- c(100, 50, 20)
  counts[2, 4:6] <- c(80, 40, 15)
  counts[3, 7:10] <- c(30, 25, 20, 10)  # Ensure all samples have at least some counts
  
  result <- .estimate_pseudocount(counts, verbose = FALSE)
  
  expect_equal(length(result$size_factors), 10)
  expect_true(all(is.finite(result$size_factors)))
  # All size factors should be positive (since all samples now have counts)
  expect_true(all(result$size_factors > 0))
})

test_that("estimate_pseudocount handles all-zero matrix", {
  counts <- matrix(0, nrow = 3, ncol = 3)
  
  result <- .estimate_pseudocount(counts, verbose = FALSE)
  
  # Should still return valid structure
  expect_equal(length(result$size_factors), 3)
  expect_is(result$scalar_pseudocount, "numeric")
})

test_that("estimate_pseudocount rejects invalid input", {
  # Non-matrix, non-SummarizedExperiment input
  expect_error(
    .estimate_pseudocount(c(1, 2, 3), verbose = FALSE),
    "must be a SummarizedExperiment or matrix"
  )
})

test_that("estimate_pseudocount is consistent across multiple calls", {
  counts <- matrix(c(10, 5, 1, 20, 8, 3, 15, 10, 5), nrow = 3, ncol = 3)
  
  result1 <- .estimate_pseudocount(counts, verbose = FALSE)
  result2 <- .estimate_pseudocount(counts, verbose = FALSE)
  
  # Same input should give identical results
  expect_equal(result1$scalar_pseudocount, result2$scalar_pseudocount)
  expect_equal(result1$size_factors, result2$size_factors)
})

test_that("estimate_pseudocount scales appropriately with sequencing depth", {
  # Small library sizes
  counts_small <- matrix(c(10, 5, 1, 20, 8, 3), nrow = 2, ncol = 3)
  result_small <- .estimate_pseudocount(counts_small, verbose = FALSE)
  
  # Large library sizes (same proportions, scaled up by 10x)
  counts_large <- matrix(c(100, 50, 10, 200, 80, 30), nrow = 2, ncol = 3)
  result_large <- .estimate_pseudocount(counts_large, verbose = FALSE)
  
  # Both should return valid results
  expect_equal(length(result_small$size_factors), 3)
  expect_equal(length(result_large$size_factors), 3)
  
  # Larger library size should result in larger pseudocount
  expect_true(result_large$scalar_pseudocount > result_small$scalar_pseudocount)
})
