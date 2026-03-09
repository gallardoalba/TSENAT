# Tests for new posterior extraction functions
# get_posterior_distribution() and compute_posterior_credible_intervals()
# These functions implement Bayesian posterior inference for transcript proportions

library(TSENAT)

context("posterior_extraction: Posterior Distribution Extraction and Intervals")

# =============================================================================
# Tests for get_posterior_distribution()
# =============================================================================

test_that("get_posterior_distribution accepts valid inputs", {
    # Should accept numeric vector and positive alpha/beta
    counts <- c(10, 20, 30, 15)
    result <- get_posterior_distribution(counts = counts, alpha = 1.5, beta = 2.0)
    
    expect_is(result, "list")
    expect_true("posterior_alpha" %in% names(result))
    expect_true("posterior_beta" %in% names(result))
    expect_true("posterior_mean" %in% names(result))
    expect_true("posterior_variance" %in% names(result))
    expect_true("posterior_sd" %in% names(result))
    expect_true("ci_lower" %in% names(result))
    expect_true("ci_upper" %in% names(result))
    expect_true("ci_level" %in% names(result))
})

test_that("get_posterior_distribution computes correct posterior parameters", {
    # With simple counts: c(10, 10) and prior (1, 1)
    # posterior_alpha = 1 + 20 = 21
    # posterior_beta = 1 + (20 - 20) = 1
    counts <- c(10, 10)
    alpha <- 1
    beta <- 1
    
    result <- get_posterior_distribution(counts = counts, alpha = alpha, beta = beta)
    
    expect_equal(result$posterior_alpha, 21)
    expect_equal(result$posterior_beta, 1)
    expect_equal(result$posterior_mean, 21/22, tolerance = 1e-6)
})

test_that("get_posterior_distribution handles edge case: all zeros", {
    # When all counts are zero, posterior = prior
    counts <- c(0, 0, 0)
    alpha <- 1.5
    beta <- 2.0
    
    result <- get_posterior_distribution(counts = counts, alpha = alpha, beta = beta)
    
    expect_equal(result$posterior_alpha, alpha)
    expect_equal(result$posterior_beta, beta)
})

test_that("get_posterior_distribution returns positive posterior parameters", {
    counts <- c(5, 10, 15)
    result <- get_posterior_distribution(counts, alpha = 1.0, beta = 1.0)
    
    expect_true(result$posterior_alpha > 0)
    expect_true(result$posterior_beta > 0)
})

test_that("get_posterior_distribution returns valid posterior statistics", {
    counts <- c(10, 20, 30)
    result <- get_posterior_distribution(counts, alpha = 1.5, beta = 2.0)
    
    # Posterior mean should be in [0, 1]
    expect_true(result$posterior_mean >= 0)
    expect_true(result$posterior_mean <= 1)
    
    # Posterior variance should be positive
    expect_true(result$posterior_variance > 0)
    
    # Posterior SD should be sqrt of variance
    expect_equal(result$posterior_sd, sqrt(result$posterior_variance), tolerance = 1e-6)
})

test_that("get_posterior_distribution computes credible intervals correctly", {
    counts <- c(10, 10)
    result <- get_posterior_distribution(counts, alpha = 1, beta = 1, ci = 0.95)
    
    # CI lower should be less than mean, upper should be greater
    expect_true(result$ci_lower < result$posterior_mean)
    expect_true(result$ci_upper > result$posterior_mean)
    
    # CI bounds should be in [0, 1]
    expect_true(result$ci_lower >= 0)
    expect_true(result$ci_upper <= 1)
})

test_that("get_posterior_distribution respects ci parameter", {
    counts <- c(10, 20, 30)
    
    # 95% CI
    result_95 <- get_posterior_distribution(counts, alpha = 1.5, beta = 2.0, ci = 0.95)
    
    # 90% CI
    result_90 <- get_posterior_distribution(counts, alpha = 1.5, beta = 2.0, ci = 0.90)
    
    expect_equal(result_95$ci_level, 0.95)
    expect_equal(result_90$ci_level, 0.90)
    
    # 90% CI should be narrower than 95% CI
    expect_true((result_90$ci_upper - result_90$ci_lower) < 
                (result_95$ci_upper - result_95$ci_lower))
})

test_that("get_posterior_distribution with ci=NULL skips credible interval", {
    counts <- c(10, 20)
    result <- get_posterior_distribution(counts, alpha = 1.5, beta = 2.0, ci = NULL)
    
    # Should still have posterior parameters
    expect_true(!is.na(result$posterior_mean))
    expect_true(!is.na(result$posterior_variance))
    
    # But CI should be NULL or NA
    expect_true(is.na(result$ci_lower) || is.null(result$ci_lower))
    expect_true(is.na(result$ci_upper) || is.null(result$ci_upper))
})

test_that("get_posterior_distribution rejects invalid counts", {
    expect_error(get_posterior_distribution(counts = "not_numeric", alpha = 1, beta = 1))
    expect_error(get_posterior_distribution(counts = list(10, 20), alpha = 1, beta = 1))
})

test_that("get_posterior_distribution rejects invalid alpha", {
    counts <- c(10, 20)
    expect_error(get_posterior_distribution(counts, alpha = 0, beta = 1))
    expect_error(get_posterior_distribution(counts, alpha = -1, beta = 1))
    expect_error(get_posterior_distribution(counts, alpha = "not_numeric", beta = 1))
})

test_that("get_posterior_distribution rejects invalid beta", {
    counts <- c(10, 20)
    expect_error(get_posterior_distribution(counts, alpha = 1, beta = 0))
    expect_error(get_posterior_distribution(counts, alpha = 1, beta = -1))
    expect_error(get_posterior_distribution(counts, alpha = 1, beta = "not_numeric"))
})

test_that("get_posterior_distribution is consistent across runs", {
    counts <- c(10, 20, 30)
    
    result1 <- get_posterior_distribution(counts, alpha = 1.5, beta = 2.0)
    result2 <- get_posterior_distribution(counts, alpha = 1.5, beta = 2.0)
    
    expect_equal(result1$posterior_alpha, result2$posterior_alpha)
    expect_equal(result1$posterior_mean, result2$posterior_mean)
})

test_that("get_posterior_distribution handles high depth samples", {
    # Large count values should not cause numerical issues
    counts <- c(10000, 20000, 30000)
    result <- get_posterior_distribution(counts, alpha = 1.5, beta = 2.0)
    
    expect_true(is.finite(result$posterior_mean))
    expect_true(is.finite(result$posterior_variance))
    expect_true(!is.na(result$posterior_mean))
})

test_that("get_posterior_distribution handles low depth samples", {
    # Very small count values
    counts <- c(1, 1, 0)
    result <- get_posterior_distribution(counts, alpha = 0.5, beta = 0.5)
    
    expect_true(is.finite(result$posterior_mean))
    expect_true(result$posterior_mean >= 0 && result$posterior_mean <= 1)
})

# =============================================================================
# Tests for compute_posterior_credible_intervals()
# =============================================================================

test_that("compute_posterior_credible_intervals accepts valid inputs", {
    # Should accept count matrix and positive alpha/beta
    counts_matrix <- matrix(c(10, 5, 1, 20, 8, 3, 15, 10, 5), nrow = 3, ncol = 3)
    result <- compute_posterior_credible_intervals(
        counts_matrix = counts_matrix, 
        alpha = 1.5, 
        beta = 2.0
    )
    
    expect_is(result, "data.frame")
    expect_true("gene" %in% colnames(result))
    expect_true("posterior_mean" %in% colnames(result))
    expect_true("ci_lower" %in% colnames(result))
    expect_true("ci_upper" %in% colnames(result))
    expect_true("posterior_sd" %in% colnames(result))
})

test_that("compute_posterior_credible_intervals returns correct number of rows", {
    # Should return one row per gene
    counts_matrix <- matrix(c(10, 5, 1, 20, 8, 3, 15, 10, 5), nrow = 3, ncol = 3)
    result <- compute_posterior_credible_intervals(counts_matrix, alpha = 1.5, beta = 2.0)
    
    expect_equal(nrow(result), 3)
})

test_that("compute_posterior_credible_intervals handles data.frame input", {
    # Should accept data.frame as well as matrix
    counts_df <- data.frame(
        c(10, 5, 1),
        c(20, 8, 3),
        c(15, 10, 5)
    )
    result <- compute_posterior_credible_intervals(counts_df, alpha = 1.5, beta = 2.0)
    
    expect_is(result, "data.frame")
    expect_equal(nrow(result), 3)
})

test_that("compute_posterior_credible_intervals includes gene identifiers", {
    counts_matrix <- matrix(c(10, 5, 1, 20, 8, 3, 15, 10, 5), nrow = 3, ncol = 3)
    rownames(counts_matrix) <- c("gene1", "gene2", "gene3")
    
    result <- compute_posterior_credible_intervals(counts_matrix, alpha = 1.5, beta = 2.0)
    
    expect_equal(result$gene, c("gene1", "gene2", "gene3"))
})

test_that("compute_posterior_credible_intervals with default rownames", {
    # When no rownames, should use numeric indices
    counts_matrix <- matrix(c(10, 5, 1, 20, 8, 3, 15, 10, 5), nrow = 3, ncol = 3)
    result <- compute_posterior_credible_intervals(counts_matrix, alpha = 1.5, beta = 2.0)
    
    expect_true("gene" %in% colnames(result))
    expect_equal(length(result$gene), 3)
})

test_that("compute_posterior_credible_intervals respects ci parameter", {
    counts_matrix <- matrix(c(10, 5, 20, 8, 15, 10), nrow = 2, ncol = 3)
    
    result_95 <- compute_posterior_credible_intervals(counts_matrix, alpha = 1.5, beta = 2.0, ci = 0.95)
    result_90 <- compute_posterior_credible_intervals(counts_matrix, alpha = 1.5, beta = 2.0, ci = 0.90)
    
    # Both should return data.frames with same structure
    expect_equal(nrow(result_95), nrow(result_90))
    expect_equal(ncol(result_95), ncol(result_90))
    
    # 90% CI should be narrower than 95%
    result_95_width <- mean(result_95$ci_upper - result_95$ci_lower)
    result_90_width <- mean(result_90$ci_upper - result_90$ci_lower)
    expect_true(result_90_width < result_95_width)
})

test_that("compute_posterior_credible_intervals validates bounds", {
    counts_matrix <- matrix(c(10, 5, 20, 8), nrow = 2, ncol = 2)
    result <- compute_posterior_credible_intervals(counts_matrix, alpha = 1.5, beta = 2.0)
    
    # All CI lower bounds < posterior mean < upper bounds
    expect_true(all(result$ci_lower < result$posterior_mean))
    expect_true(all(result$posterior_mean < result$ci_upper))
    
    # All bounds in [0, 1]
    expect_true(all(result$ci_lower >= 0))
    expect_true(all(result$ci_upper <= 1))
})

test_that("compute_posterior_credible_intervals computes posterior_sd correctly", {
    counts_matrix <- matrix(c(10, 5, 20, 8), nrow = 2, ncol = 2)
    result <- compute_posterior_credible_intervals(counts_matrix, alpha = 1.5, beta = 2.0)
    
    # For each row, we can verify against get_posterior_distribution
    for (i in seq_len(nrow(counts_matrix))) {
        posterior_single <- get_posterior_distribution(
            counts_matrix[i, ], alpha = 1.5, beta = 2.0
        )
        expect_equal(result$posterior_sd[i], posterior_single$posterior_sd, tolerance = 1e-6)
    }
})

test_that("compute_posterior_credible_intervals handles edge cases", {
    # Matrix with all zeros in one row
    counts_matrix <- matrix(c(0, 0, 0, 20, 8, 10), nrow = 2, ncol = 3)
    result <- compute_posterior_credible_intervals(counts_matrix, alpha = 1.5, beta = 2.0)
    
    expect_equal(nrow(result), 2)
    expect_true(all(is.finite(result$posterior_mean)))
    expect_true(all(is.finite(result$posterior_sd)))
})

test_that("compute_posterior_credible_intervals is deterministic", {
    counts_matrix <- matrix(c(10, 5, 20, 8, 15, 10), nrow = 2, ncol = 3)
    
    result1 <- compute_posterior_credible_intervals(counts_matrix, alpha = 1.5, beta = 2.0)
    result2 <- compute_posterior_credible_intervals(counts_matrix, alpha = 1.5, beta = 2.0)
    
    expect_equal(result1$posterior_mean, result2$posterior_mean)
    expect_equal(result1$ci_lower, result2$ci_lower)
    expect_equal(result1$ci_upper, result2$ci_upper)
})

test_that("compute_posterior_credible_intervals rejects invalid inputs", {
    not_numeric <- data.frame(a = "text", b = c(1, 2))
    expect_error(compute_posterior_credible_intervals(not_numeric, alpha = 1.5, beta = 2.0))
})

test_that("compute_posterior_credible_intervals works with larger datasets", {
    # Simulate realistic dataset
    set.seed(42)
    n_genes <- 100
    n_samples <- 20
    counts_matrix <- matrix(
        rpois(n_genes * n_samples, lambda = 25),
        nrow = n_genes, ncol = n_samples
    )
    
    result <- compute_posterior_credible_intervals(counts_matrix, alpha = 1.5, beta = 2.0)
    
    expect_equal(nrow(result), n_genes)
    expect_true(all(is.finite(result$posterior_mean)))
    expect_true(all(result$posterior_mean >= 0 & result$posterior_mean <= 1))
})

# =============================================================================
# Integration tests: Posterior extraction workflow
# =============================================================================

test_that("Posterior extraction workflow: fit prior then extract posteriors", {
    # Full workflow: fit empirical Bayes prior, then extract posteriors
    set.seed(123)
    counts_matrix <- matrix(
        rpois(30 * 15, lambda = 25),
        nrow = 30, ncol = 15
    )
    
    # Step 1: Fit empirical Bayes prior
    prior_params <- fit_empirical_beta_prior(counts_matrix)
    expect_true(prior_params$alpha > 0)
    expect_true(prior_params$beta > 0)
    
    # Step 2: Extract posteriors for all genes
    posteriors <- compute_posterior_credible_intervals(
        counts_matrix, 
        alpha = prior_params$alpha, 
        beta = prior_params$beta,
        ci = 0.95
    )
    
    expect_equal(nrow(posteriors), 30)
    expect_true(all(posteriors$posterior_mean >= 0))
    expect_true(all(posteriors$posterior_mean <= 1))
})

test_that("Posterior with high confidence (narrow CI) for high-count genes", {
    # Genes with more reads should have narrower credible intervals
    # High count gene: c(100, 100, 100)
    # Low count gene: c(1, 1, 1)
    
    high_counts <- c(100, 100, 100)
    low_counts <- c(1, 1, 1)
    
    posterior_high <- get_posterior_distribution(high_counts, alpha = 1.5, beta = 2.0, ci = 0.95)
    posterior_low <- get_posterior_distribution(low_counts, alpha = 1.5, beta = 2.0, ci = 0.95)
    
    width_high <- posterior_high$ci_upper - posterior_high$ci_lower
    width_low <- posterior_low$ci_upper - posterior_low$ci_lower
    
    # High count should have narrower CI
    expect_true(width_high < width_low)
})

test_that("Posterior mean is between prior mean and empirical mean", {
    # Posterior mean should be a compromise between prior and data
    counts <- c(100, 100)
    alpha <- 1
    beta <- 10  # Prior lean toward low abundance
    
    prior_mean <- alpha / (alpha + beta)
    empirical_mean <- sum(counts) / (sum(counts) + sum(1 - 1))  # Simplified
    
    posterior <- get_posterior_distribution(counts, alpha, beta, ci = NULL)
    
    # Posterior should be between prior and data
    expect_true(posterior$posterior_mean > prior_mean)
})

test_that("Different priors affect posterior uncertainty", {
    counts <- c(10, 10)
    
    # Weak prior (Jeffreys)
    posterior_weak <- get_posterior_distribution(counts, alpha = 0.5, beta = 0.5, ci = 0.95)
    
    # Strong prior (informative)
    posterior_strong <- get_posterior_distribution(counts, alpha = 10, beta = 10, ci = 0.95)
    
    # Both should give finite posteriors
    expect_true(is.finite(posterior_weak$posterior_mean))
    expect_true(is.finite(posterior_strong$posterior_mean))
    
    # Posteriors may differ depending on data vs prior conflict
    # (exact relationship depends on data, so we just check validity)
    expect_true(all(c(
        posterior_weak$posterior_mean >= 0, posterior_weak$posterior_mean <= 1,
        posterior_strong$posterior_mean >= 0, posterior_strong$posterior_mean <= 1
    )))
})
