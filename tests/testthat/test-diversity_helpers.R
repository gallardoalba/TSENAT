# Comprehensive tests for diversity_helpers.R Bayesian functions
# Tests for empirical Bayes prior fitting, WLFC pseudocounts, and posterior inference

library(testthat)
library(TSENAT)

context("Bayesian Methods: Empirical Bayes Prior Estimation")

test_that("fit_empirical_beta_prior returns valid Beta parameters", {
  # Create simple count matrix
  counts <- matrix(
    c(10, 5, 1, 20, 8, 3, 15, 10, 5),
    nrow = 3, ncol = 3
  )
  
  prior <- fit_empirical_beta_prior(counts)
  
  # Check structure
  expect_is(prior, "list")
  expect_true("alpha" %in% names(prior))
  expect_true("beta" %in% names(prior))
  
  # Check positivity
  expect_true(prior$alpha > 0)
  expect_true(prior$beta > 0)
  
  # Check reasonable range (not extreme values)
  expect_true(prior$alpha < 100)
  expect_true(prior$beta < 100)
})

test_that("fit_empirical_beta_prior handles matrix with all zeros", {
  # All-zero matrix
  counts <- matrix(0, nrow = 3, ncol = 3)
  
  # Expect warning about zero variance when fitting prior on all-zeros
  prior <- expect_warning(
    fit_empirical_beta_prior(counts),
    "Variance near zero"
  )
  
  # Should fall back to uniform prior (alpha=1, beta=1) for zero variance
  expect_equal(prior$alpha, 1)
  expect_equal(prior$beta, 1)
})

test_that("fit_empirical_beta_prior handles single-row matrix", {
  counts <- matrix(c(1, 2, 3, 4, 5), nrow = 1, ncol = 5)
  
  # Single-row matrix can have valid variance when normalized across samples
  prior <- fit_empirical_beta_prior(counts)
  
  expect_true(prior$alpha > 0)
  expect_true(prior$beta > 0)
})

test_that("fit_empirical_beta_prior handles data.frame input", {
  counts_df <- data.frame(
    S1 = c(10, 5, 1),
    S2 = c(20, 8, 3),
    S3 = c(15, 10, 5)
  )
  
  prior <- fit_empirical_beta_prior(counts_df)
  
  expect_true(prior$alpha > 0)
  expect_true(prior$beta > 0)
})

test_that("fit_empirical_beta_prior rejects invalid input", {
  expect_error(fit_empirical_beta_prior(c(1, 2, 3)), 
               "must be a matrix or data.frame")
})

context("Bayesian Methods: Weighted Likelihood Fold Change Pseudocount")

test_that("compute_wlfc_pseudocounts returns scalar pseudocount", {
  counts <- c(10, 5, 2, 8, 3)
  alpha <- 0.5
  beta <- 4.0
  
  pc <- compute_wlfc_pseudocounts(counts, alpha, beta)
  
  # Check type and range
  expect_is(pc, "numeric")
  expect_length(pc, 1)
  expect_true(pc >= 0)
  expect_true(pc <= 1)  # Typically in [0, 1] range
})

test_that("compute_wlfc_pseudocounts handles zero counts", {
  counts <- c(0, 0, 0, 0, 0)
  alpha <- 0.5
  beta <- 4.0
  
  pc <- compute_wlfc_pseudocounts(counts, alpha, beta)
  
  # Should return default Jeffreys prior
  expect_equal(pc, 0.5)
})

test_that("compute_wlfc_pseudocounts varies with count magnitude", {
  alpha <- 0.5
  beta <- 4.0
  
  pc_low <- compute_wlfc_pseudocounts(c(1, 1, 1), alpha, beta)
  pc_high <- compute_wlfc_pseudocounts(c(100, 100, 100), alpha, beta)
  
  # Higher counts should give different (typically higher) pseudocounts
  expect_false(pc_low == pc_high)
})

test_that("compute_wlfc_pseudocounts rejects invalid parameters", {
  counts <- c(10, 5, 2)
  
  # Invalid alpha
  expect_error(compute_wlfc_pseudocounts(counts, -0.5, 4.0),
               "alpha must be a positive numeric value")
  
  # Invalid beta
  expect_error(compute_wlfc_pseudocounts(counts, 0.5, -4.0),
               "beta must be a positive numeric value")
})

context("Bayesian Methods: Posterior Distribution Extraction")

test_that("get_posterior_distribution returns valid posterior parameters", {
  counts <- c(10, 5, 2, 8, 3)
  alpha <- 0.5
  beta <- 4.0
  
  posterior <- get_posterior_distribution(counts, alpha, beta, ci = 0.95)
  
  # Check structure
  expect_is(posterior, "list")
  expect_true("posterior_alpha" %in% names(posterior))
  expect_true("posterior_beta" %in% names(posterior))
  expect_true("posterior_mean" %in% names(posterior))
  expect_true("posterior_variance" %in% names(posterior))
  expect_true("posterior_sd" %in% names(posterior))
  expect_true("ci_lower" %in% names(posterior))
  expect_true("ci_upper" %in% names(posterior))
  
  # Check positivity
  expect_true(posterior$posterior_alpha > 0)
  expect_true(posterior$posterior_beta > 0)
  expect_true(posterior$posterior_mean > 0)
  expect_true(posterior$posterior_mean < 1)
  expect_true(posterior$posterior_variance > 0)
  expect_true(posterior$posterior_sd > 0)
})

test_that("get_posterior_distribution computes valid credible interval", {
  counts <- c(10, 5, 2, 8, 3)
  alpha <- 0.5
  beta <- 4.0
  
  posterior <- get_posterior_distribution(counts, alpha, beta, ci = 0.95)
  
  # CI bounds should be in valid range
  expect_true(posterior$ci_lower >= 0)
  expect_true(posterior$ci_upper <= 1)
  expect_true(posterior$ci_lower < posterior$ci_upper)
  
  # Posterior mean should be within CI
  expect_true(posterior$posterior_mean >= posterior$ci_lower - 1e-6)
  expect_true(posterior$posterior_mean <= posterior$ci_upper + 1e-6)
})

test_that("get_posterior_distribution skips CI when ci=NULL", {
  counts <- c(10, 5, 2, 8, 3)
  alpha <- 0.5
  beta <- 4.0
  
  posterior <- get_posterior_distribution(counts, alpha, beta, ci = NULL)
  
  # CI components should not exist
  expect_false("ci_lower" %in% names(posterior))
  expect_false("ci_upper" %in% names(posterior))
  expect_false("ci_level" %in% names(posterior))
  
  # But posterior statistics should exist
  expect_true("posterior_mean" %in% names(posterior))
  expect_true("posterior_sd" %in% names(posterior))
})

test_that("get_posterior_distribution CI width decreases with more data", {
  alpha <- 0.5
  beta <- 4.0
  
  # Small data
  post_small <- get_posterior_distribution(c(5, 5), alpha, beta, ci = 0.95)
  ci_width_small <- post_small$ci_upper - post_small$ci_lower
  
  # Large data (same proportions)
  post_large <- get_posterior_distribution(c(500, 500), alpha, beta, ci = 0.95)
  ci_width_large <- post_large$ci_upper - post_large$ci_lower
  
  # More data should give narrower CI
  expect_true(ci_width_large < ci_width_small)
})

test_that("get_posterior_distribution rejects invalid input", {
  alpha <- 0.5
  beta <- 4.0
  
  # Empty counts
  expect_error(get_posterior_distribution(c(), alpha, beta),
               "counts must be a non-empty numeric vector")
  
  # Invalid CI
  expect_error(get_posterior_distribution(c(10, 5), alpha, beta, ci = 1.5),
               "ci must be NULL or a numeric value between 0 and 1")
})

context("Bayesian Methods: Posterior Credible Intervals")

test_that("compute_posterior_credible_intervals returns data frame with expected structure", {
  counts <- matrix(c(10, 5, 1, 20, 8, 3, 15, 10, 5), nrow = 3, ncol = 3)
  rownames(counts) <- c("Gene1", "Gene2", "Gene3")
  
  alpha <- 0.5
  beta <- 4.0
  
  cis <- compute_posterior_credible_intervals(counts, alpha, beta, ci = 0.95)
  
  # Check structure
  expect_is(cis, "data.frame")
  expect_equal(nrow(cis), 3)
  
  # Check columns
  expect_true("gene" %in% colnames(cis))
  expect_true("posterior_mean" %in% colnames(cis))
  expect_true("ci_lower" %in% colnames(cis))
  expect_true("ci_upper" %in% colnames(cis))
  expect_true("posterior_sd" %in% colnames(cis))
  
  # Check values
  expect_equal(cis$gene, c("Gene1", "Gene2", "Gene3"))
})

test_that("compute_posterior_credible_intervals CI bounds are valid", {
  counts <- matrix(c(10, 5, 1, 20, 8, 3, 15, 10, 5), nrow = 3, ncol = 3)
  
  cis <- compute_posterior_credible_intervals(counts, alpha = 0.5, beta = 4.0)
  
  # All CIs should be valid
  for (i in 1:nrow(cis)) {
    expect_true(cis$ci_lower[i] >= 0, info = paste("Gene", i, "CI lower"))
    expect_true(cis$ci_upper[i] <= 1, info = paste("Gene", i, "CI upper"))
    expect_true(cis$ci_lower[i] <= cis$ci_upper[i], 
                info = paste("Gene", i, "CI order"))
    expect_true(cis$posterior_mean[i] >= cis$ci_lower[i] - 1e-6,
                info = paste("Gene", i, "mean below lower"))
    expect_true(cis$posterior_mean[i] <= cis$ci_upper[i] + 1e-6,
                info = paste("Gene", i, "mean above upper"))
  }
})

test_that("compute_posterior_credible_intervals handles matrix without rownames", {
  counts <- matrix(c(10, 5, 1, 20, 8, 3, 15, 10, 5), nrow = 3, ncol = 3)
  
  cis <- compute_posterior_credible_intervals(counts, alpha = 0.5, beta = 4.0)
  
  # Should auto-generate gene names
  expect_equal(cis$gene, c("Gene_1", "Gene_2", "Gene_3"))
})

context("Bootstrap Methods: Confidence Intervals")

test_that(".tsenat_bootstrap_resample produces valid bootstrap distribution", {
  x <- c(10, 5, 2, 8, 3)
  
  boot_dist <- .tsenat_bootstrap_resample(
    x = x,
    q = 2,
    norm = TRUE,
    nboot = 100,
    log_base = exp(1),
    pseudocount = 0.1,
    what = "S"
  )
  
  # Check type and length
  expect_is(boot_dist, "numeric")
  expect_length(boot_dist, 100)
  
  # All bootstrap replicates should be finite
  expect_true(all(is.finite(boot_dist)))
})

test_that(".tsenat_ci_percentile computes valid percentile CI", {
  boot_dist <- rnorm(1000, mean = 0.5, sd = 0.1)
  
  ci <- .tsenat_ci_percentile(boot_dist, ci = 0.95)
  
  # Check structure
  expect_is(ci, "list")
  expect_true("lower" %in% names(ci))
  expect_true("upper" %in% names(ci))
  
  # Check ordering
  expect_true(ci$lower < ci$upper)
  
  # Check approximate coverage (should contain ~95% of data)
  coverage <- sum(boot_dist >= ci$lower & boot_dist <= ci$upper) / length(boot_dist)
  expect_true(coverage >= 0.90)  # Allow some tolerance
  expect_true(coverage <= 1.00)
})

context("Bayesian Methods: Shrinkage Parameter Estimation")

test_that(".tsenat_estimate_shrinkage_params returns valid structure", {
  # Create count matrix
  counts <- matrix(rpois(60, lambda = 5), nrow = 10, ncol = 6)
  genes <- rep(c("G1", "G2"), each = 5)
  
  # Create entropy matrix
  entropy_matrix <- matrix(runif(20, 0, 1), nrow = 10, ncol = 2)
  colnames(entropy_matrix) <- c("S1_q=2", "S2_q=2")
  rownames(entropy_matrix) <- paste0("Gene_", 1:10)
  
  params <- .tsenat_estimate_shrinkage_params(
    x = counts,
    genes = genes,
    entropy_matrix = entropy_matrix,
    q = 2,
    min_count = 1
  )
  
  # Check structure
  expect_is(params, "list")
  expect_true("global_mean" %in% names(params))
  expect_true("global_var" %in% names(params))
  expect_true("n_isoforms" %in% names(params))
  
  # Check n_isoforms output
  expect_equal(length(params$n_isoforms), 2)  # 2 unique genes
  expect_true(all(params$n_isoforms > 0))
})

test_that(".tsenat_apply_shrinkage produces output same size as input", {
  # Create test data
  entropy_matrix <- matrix(runif(50, 0, 1), nrow = 10, ncol = 5)
  rownames(entropy_matrix) <- paste0("Gene_", 1:10)
  
  params <- list(
    global_mean = c("q=2" = 0.5),
    global_var = c("q=2" = 0.05),
    n_isoforms = setNames(rep(3, 10), paste0("Gene_", 1:10))
  )
  
  colnames(entropy_matrix) <- c(
    "S1_q=2", "S2_q=2", "S3_q=2", "S4_q=2", "S5_q=2"
  )
  
  shrunk <- .tsenat_apply_shrinkage(entropy_matrix, params)
  
  # Check dimensions
  expect_equal(dim(shrunk), dim(entropy_matrix))
  
  # Check all values are finite
  expect_true(all(is.finite(shrunk)))
})

test_that(".tsenat_apply_shrinkage shrinks toward global mean", {
  # Create test data with known global mean
  global_mean <- 0.5
  entropy_matrix <- matrix(c(0.1, 0.9, 0.2, 0.8), nrow = 2, ncol = 2)
  rownames(entropy_matrix) <- c("Gene_1", "Gene_2")
  colnames(entropy_matrix) <- c("S1_q=2", "S2_q=2")
  
  params <- list(
    global_mean = c("q=2" = global_mean),
    global_var = c("q=2" = 0.01),
    n_isoforms = setNames(c(2, 2), c("Gene_1", "Gene_2"))
  )
  
  shrunk <- .tsenat_apply_shrinkage(entropy_matrix, params)
  
  # Shrunk values should be closer to global mean than originals
  for (i in 1:nrow(entropy_matrix)) {
    for (j in 1:ncol(entropy_matrix)) {
      orig_dist <- abs(entropy_matrix[i, j] - global_mean)
      shrunk_dist <- abs(shrunk[i, j] - global_mean)
      expect_true(shrunk_dist <= orig_dist + 1e-6,
                  info = paste("Cell", i, j, "should shrink toward mean"))
    }
  }
})

context("Bayesian Methods: Integration Tests")

test_that("Full Bayesian posterior workflow executes without error", {
  # Generate synthetic count data
  set.seed(123)
  counts <- matrix(rpois(300, lambda = 10), nrow = 30, ncol = 10)
  rownames(counts) <- paste0("Gene_", 1:30)
  
  # Fit empirical Bayes prior
  prior <- fit_empirical_beta_prior(counts)
  expect_true(prior$alpha > 0)
  expect_true(prior$beta > 0)
  
  # Compute WLFC pseudocounts for first gene
  pc <- compute_wlfc_pseudocounts(counts[1, ], prior$alpha, prior$beta)
  expect_true(is.finite(pc))
  expect_true(pc >= 0)
  
  # Get posterior for first gene
  posterior <- get_posterior_distribution(counts[1, ], prior$alpha, prior$beta)
  expect_true(is.finite(posterior$posterior_mean))
  expect_true(posterior$posterior_mean > 0 & posterior$posterior_mean < 1)
  
  # Compute posteriors for all genes
  all_posteriors <- compute_posterior_credible_intervals(
    counts, prior$alpha, prior$beta
  )
  expect_equal(nrow(all_posteriors), 30)
})

test_that("Posterior mean estimates are valid across different data sizes", {
  alpha <- 1.0
  beta <- 1.0
  
  # Test with different data sizes
  post_low <- get_posterior_distribution(c(2, 1), alpha, beta)
  post_mid <- get_posterior_distribution(c(20, 10), alpha, beta)
  post_high <- get_posterior_distribution(c(200, 100), alpha, beta)
  
  # All posterior means should be valid probabilities
  for (mean in c(post_low$posterior_mean, post_mid$posterior_mean, post_high$posterior_mean)) {
    expect_true(mean > 0)
    expect_true(mean < 1)
    expect_true(is.finite(mean))
  }
  
  # All posterior sds should be positive and decreasing with more data
  # (more data = lower posterior uncertainty)
  expect_true(post_low$posterior_sd > 0)
  expect_true(post_mid$posterior_sd > 0)
  expect_true(post_high$posterior_sd > 0)
  
  # Posterior SD should decrease with more data (narrower uncertainty)
  expect_true(post_low$posterior_sd >= post_mid$posterior_sd - 1e-6)
  expect_true(post_mid$posterior_sd >= post_high$posterior_sd - 1e-6)
})

context("Bayesian Methods: Integrated WLFC Workflow")

test_that("estimate_wlfc_pseudocounts works with matrix input", {
  counts <- matrix(c(10, 5, 1, 20, 8, 3, 15, 10, 5), nrow = 3, ncol = 3)
  rownames(counts) <- c("Gene1", "Gene2", "Gene3")
  
  result <- estimate_wlfc_pseudocounts(counts, verbose = FALSE)
  
  # Check return structure
  expect_is(result, "list")
  expect_true("pseudocounts" %in% names(result))
  expect_true("scalar_pseudocount" %in% names(result))
  expect_true("prior" %in% names(result))
  expect_true("diagnostics" %in% names(result))
  
  # Check pseudocounts
  expect_is(result$pseudocounts, "numeric")
  expect_equal(length(result$pseudocounts), 3)
  expect_equal(names(result$pseudocounts), c("Gene1", "Gene2", "Gene3"))
  
  # Check scalar pseudocount
  expect_is(result$scalar_pseudocount, "numeric")
  expect_equal(result$scalar_pseudocount, mean(result$pseudocounts))
  expect_true(result$scalar_pseudocount > 0)
})

test_that("estimate_wlfc_pseudocounts works with SummarizedExperiment input", {
  skip_if_not_installed("SummarizedExperiment")
  
  counts <- matrix(c(10, 5, 1, 20, 8, 3, 15, 10, 5), nrow = 3, ncol = 3)
  rownames(counts) <- c("Gene1", "Gene2", "Gene3")
  se <- suppressWarnings(SummarizedExperiment::SummarizedExperiment(assay = list(data = counts)))
  
  result <- estimate_wlfc_pseudocounts(se, verbose = FALSE)
  
  # Check return structure
  expect_is(result, "list")
  expect_equal(length(result$pseudocounts), 3)
  expect_equal(names(result$pseudocounts), c("Gene1", "Gene2", "Gene3"))
})

test_that("estimate_wlfc_pseudocounts prior parameters are positive", {
  counts <- matrix(c(10, 5, 1, 20, 8, 3, 15, 10, 5), nrow = 3, ncol = 3)
  
  result <- estimate_wlfc_pseudocounts(counts, verbose = FALSE)
  
  # Check prior structure
  expect_is(result$prior, "list")
  expect_true("alpha" %in% names(result$prior))
  expect_true("beta" %in% names(result$prior))
  
  # Check positivity
  expect_true(result$prior$alpha > 0)
  expect_true(result$prior$beta > 0)
  
  # Check reasonable range
  expect_true(result$prior$alpha < 100)
  expect_true(result$prior$beta < 100)
})

test_that("estimate_wlfc_pseudocounts returns diagnostic information", {
  counts <- matrix(c(10, 5, 1, 20, 8, 3, 15, 10, 5), nrow = 3, ncol = 3)
  rownames(counts) <- c("Gene1", "Gene2", "Gene3")
  
  result <- estimate_wlfc_pseudocounts(counts, verbose = FALSE)
  
  # Check diagnostics
  expect_is(result$diagnostics, "list")
  expect_equal(result$diagnostics$unique_rownames, 3)
  expect_equal(result$diagnostics$n_genes_filtered, 3)
  expect_equal(result$diagnostics$n_samples, 3)
  expect_null(result$diagnostics$duplicate_genes)
})

test_that("estimate_wlfc_pseudocounts detects duplicate gene rownames", {
  counts <- matrix(c(10, 5, 1, 20, 8, 3, 15, 10, 5), nrow = 3, ncol = 3)
  rownames(counts) <- c("Gene1", "Gene1", "Gene2")
  
  result <- estimate_wlfc_pseudocounts(counts, verbose = FALSE)
  
  # Check duplicate detection
  expect_equal(result$diagnostics$unique_rownames, 2)
  expect_true(!is.null(result$diagnostics$duplicate_genes))
  expect_true("Gene1" %in% result$diagnostics$duplicate_genes)
})

test_that("estimate_wlfc_pseudocounts pseudocount values are in valid range", {
  counts <- matrix(c(10, 5, 1, 20, 8, 3, 15, 10, 5), nrow = 3, ncol = 3)
  
  result <- estimate_wlfc_pseudocounts(counts, verbose = FALSE)
  
  # All pseudocounts should be non-negative
  expect_true(all(result$pseudocounts >= 0))
  
  # All should be less than 1 (typical for WLFC)
  expect_true(all(result$pseudocounts <= 1))
  
  # Scalar should be mean
  expect_equal(result$scalar_pseudocount, mean(result$pseudocounts))
})

test_that("estimate_wlfc_pseudocounts handles sparse counts", {
  # Matrix with mostly zeros
  counts <- matrix(0, nrow = 5, ncol = 10)
  counts[1, 1:3] <- c(100, 50, 20)
  counts[2, 4:6] <- c(80, 40, 15)
  
  result <- estimate_wlfc_pseudocounts(counts, verbose = FALSE)
  
  expect_equal(length(result$pseudocounts), 5)
  expect_true(all(is.finite(result$pseudocounts)))
  expect_true(all(result$pseudocounts >= 0))
})

test_that("estimate_wlfc_pseudocounts handles all-zero matrix", {
  counts <- matrix(0, nrow = 3, ncol = 3)
  
  result <- expect_warning(
    estimate_wlfc_pseudocounts(counts, verbose = FALSE),
    "Variance near zero|cannot compute"
  )
  
  # Should still return valid structure (with fallback values)
  expect_equal(length(result$pseudocounts), 3)
})

test_that("estimate_wlfc_pseudocounts rejects invalid input", {
  # Non-matrix, non-SummarizedExperiment input
  expect_error(
    estimate_wlfc_pseudocounts(c(1, 2, 3), verbose = FALSE),
    "must be a SummarizedExperiment or matrix"
  )
  
  # Data frame (not supported directly)
  expect_error(
    estimate_wlfc_pseudocounts(data.frame(x = 1:3), verbose = FALSE),
    "must be a SummarizedExperiment or matrix"
  )
})

test_that("estimate_wlfc_pseudocounts verbose output is informative", {
  counts <- matrix(c(10, 5, 1, 20, 8, 3, 15, 10, 5), nrow = 3, ncol = 3)
  rownames(counts) <- c("Gene1", "Gene2", "Gene3")
  
  # Capture output when verbose=TRUE
  captured_output <- capture.output({
    result <- estimate_wlfc_pseudocounts(counts, verbose = TRUE)
  })
  
  # Should contain diagnostic information
  expect_true(any(grepl("Diagnostic", captured_output)))
  expect_true(any(grepl("Rows", captured_output)))
  expect_true(any(grepl("Empirical Beta prior", captured_output)))
  expect_true(any(grepl("WLFC Pseudocount Distribution", captured_output)))
  expect_true(any(grepl("Mean", captured_output)))
})

test_that("estimate_wlfc_pseudocounts consistent across multiple calls", {
  counts <- matrix(c(10, 5, 1, 20, 8, 3, 15, 10, 5), nrow = 3, ncol = 3)
  
  result1 <- estimate_wlfc_pseudocounts(counts, verbose = FALSE)
  result2 <- estimate_wlfc_pseudocounts(counts, verbose = FALSE)
  
  # Same input should give identical results
  expect_equal(result1$pseudocounts, result2$pseudocounts)
  expect_equal(result1$scalar_pseudocount, result2$scalar_pseudocount)
  expect_equal(result1$prior$alpha, result2$prior$alpha)
  expect_equal(result1$prior$beta, result2$prior$beta)
})

test_that("estimate_wlfc_pseudocounts scales appropriately with sample size", {
  # Small sample
  counts_small <- matrix(c(10, 5, 1, 20, 8, 3), nrow = 2, ncol = 3)
  result_small <- estimate_wlfc_pseudocounts(counts_small, verbose = FALSE)
  
  # Large sample (same proportions, scaled up)
  counts_large <- matrix(c(100, 50, 10, 200, 80, 30), nrow = 2, ncol = 3)
  result_large <- estimate_wlfc_pseudocounts(counts_large, verbose = FALSE)
  
  # Both should return valid results with similar structure
  expect_equal(length(result_small$pseudocounts), 2)
  expect_equal(length(result_large$pseudocounts), 2)
  
  # Pseudocounts may differ (effect of scaling on posterior), but all valid
  expect_true(all(result_small$pseudocounts >= 0))
  expect_true(all(result_large$pseudocounts >= 0))
})
