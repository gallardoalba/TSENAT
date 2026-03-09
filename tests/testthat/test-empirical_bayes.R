# Tests for empirical Bayes functions: fit_empirical_beta_prior and compute_wlfc_pseudocounts
# These functions estimate data-adaptive pseudocounts using Beta priors

library(TSENAT)

context("empirical_bayes: Shrinkage and Posterior Estimation")

# =============================================================================
# Tests for fit_empirical_beta_prior()
# =============================================================================

test_that("fit_empirical_beta_prior accepts matrix input", {
    # Should accept numeric matrix
    counts <- matrix(c(10, 5, 1, 20, 8, 3, 15, 10, 5), nrow = 3, ncol = 3)
    result <- fit_empirical_beta_prior(counts)
    
    expect_is(result, "list")
    expect_true("alpha" %in% names(result))
    expect_true("beta" %in% names(result))
})

test_that("fit_empirical_beta_prior accepts data.frame input", {
    # Should accept data.frame as well
    counts <- data.frame(c(10, 5, 1), c(20, 8, 3), c(15, 10, 5))
    result <- fit_empirical_beta_prior(counts)
    
    expect_is(result, "list")
    expect_true(is.numeric(result$alpha))
    expect_true(is.numeric(result$beta))
})

test_that("fit_empirical_beta_prior rejects invalid input", {
    # Should reject non-matrix, non-dataframe inputs
    expect_error(fit_empirical_beta_prior(c(10, 20, 30)))
    expect_error(fit_empirical_beta_prior("not_a_matrix"))
    expect_error(fit_empirical_beta_prior(list(10, 20, 30)))
})

test_that("fit_empirical_beta_prior returns positive parameters for normal data", {
    # For typical count data, both alpha and beta should be positive
    counts <- matrix(c(
        c(100, 50, 25, 10),
        c(80, 60, 40, 20),
        c(120, 70, 30, 10)
    ), nrow = 3, ncol = 4, byrow = TRUE)
    
    result <- fit_empirical_beta_prior(counts)
    
    expect_true(result$alpha > 0, info = "alpha should be positive")
    expect_true(result$beta > 0, info = "beta should be positive")
})

test_that("fit_empirical_beta_prior handles uniform distribution", {
    # For uniform counts, should return positive parameters and warn about variance
    counts <- matrix(rep(50, 12), nrow = 3, ncol = 4)
    
    # Uniform distribution causes zero variance, which triggers warning
    expect_warning(
        result <- fit_empirical_beta_prior(counts),
        "Variance near zero"
    )
    
    expect_is(result, "list")
    expect_true(result$alpha > 0)
    expect_true(result$beta > 0)
})

test_that("fit_empirical_beta_prior handles zero variance gracefully", {
    # When variance is very low, should return fallback prior and warn
    counts <- matrix(rep(50, 12), nrow = 3, ncol = 4)  # All same values
    
    # This might trigger a variance warning
    expect_warning(
        result <- fit_empirical_beta_prior(counts),
        info = "Should warn when variance is near zero"
    )
    
    # Should return some reasonable prior
    expect_is(result, "list")
    expect_true(result$alpha > 0)
    expect_true(result$beta > 0)
})

test_that("fit_empirical_beta_prior works with realistic larger count matrix", {
    # Test with realistic count data (larger matrix simulating real dataset)
    set.seed(42)
    n_genes <- 100
    n_samples <- 20
    # Simulate realistic count data with some genes high/low abundance
    counts <- matrix(
        c(rpois(n_genes * n_samples * 0.7, lambda = 50),
          rpois(n_genes * n_samples * 0.3, lambda = 5)),
        nrow = n_genes, ncol = n_samples
    )
    
    result <- fit_empirical_beta_prior(counts)
    
    expect_is(result, "list")
    expect_true(is.numeric(result$alpha))
    expect_true(is.numeric(result$beta))
    expect_true(result$alpha > 0)
    expect_true(result$beta > 0)
})

test_that("fit_empirical_beta_prior returns reasonable prior parameters", {
    # Parameters should be in reasonable range for Beta distribution
    counts <- matrix(
        c(100, 50, 25, 10, 80, 60, 40, 20, 120, 70, 30, 10),
        nrow = 3, ncol = 4, byrow = TRUE
    )
    
    result <- fit_empirical_beta_prior(counts)
    
    # Alpha and beta are usually < 10 for typical count data
    expect_true(result$alpha > 0.1, info = "alpha should be substantially positive")
    expect_true(result$beta > 0.1, info = "beta should be substantially positive")
})

test_that("fit_empirical_beta_prior is consistent across runs", {
    # Function should be deterministic
    counts <- matrix(c(10, 5, 1, 20, 8, 3, 15, 10, 5), nrow = 3, ncol = 3)
    
    result1 <- fit_empirical_beta_prior(counts)
    result2 <- fit_empirical_beta_prior(counts)
    
    expect_equal(result1$alpha, result2$alpha)
    expect_equal(result1$beta, result2$beta)
})

# =============================================================================
# Tests for compute_wlfc_pseudocounts()
# =============================================================================

test_that("compute_wlfc_pseudocounts accepts valid inputs", {
    # Should accept numeric vector and positive alpha/beta
    counts <- c(10, 20, 30, 15)
    result <- compute_wlfc_pseudocounts(counts, alpha = 1.5, beta = 2.0)
    
    expect_is(result, "numeric")
    expect_length(result, 1)
})

test_that("compute_wlfc_pseudocounts rejects non-numeric counts", {
    # Should reject non-numeric counts
    expect_error(compute_wlfc_pseudocounts(c("a", "b", "c"), 1.5, 2.0))
    expect_error(compute_wlfc_pseudocounts(list(10, 20), 1.5, 2.0))
})

test_that("compute_wlfc_pseudocounts rejects invalid alpha", {
    # Should reject non-positive or non-numeric alpha
    expect_error(compute_wlfc_pseudocounts(c(10, 20), alpha = 0, beta = 2.0))
    expect_error(compute_wlfc_pseudocounts(c(10, 20), alpha = -1, beta = 2.0))
    expect_error(compute_wlfc_pseudocounts(c(10, 20), alpha = "not_numeric", beta = 2.0))
})

test_that("compute_wlfc_pseudocounts rejects invalid beta", {
    # Should reject non-positive or non-numeric beta
    expect_error(compute_wlfc_pseudocounts(c(10, 20), alpha = 1.5, beta = 0))
    expect_error(compute_wlfc_pseudocounts(c(10, 20), alpha = 1.5, beta = -1))
    expect_error(compute_wlfc_pseudocounts(c(10, 20), alpha = 1.5, beta = "not_numeric"))
})

test_that("compute_wlfc_pseudocounts handles all-zero counts", {
    # When counts are all zero, should return default ~0.5 (Jeffreys prior)
    counts <- c(0, 0, 0, 0)
    result <- compute_wlfc_pseudocounts(counts, alpha = 1.5, beta = 2.0)
    
    expect_is(result, "numeric")
    expect_equal(result, 0.5, info = "Should return Jeffreys prior for all-zero counts")
})

test_that("compute_wlfc_pseudocounts returns positive value", {
    # For positive counts and valid priors, result should be positive
    counts <- c(10, 20, 30, 15)
    result <- compute_wlfc_pseudocounts(counts, alpha = 1.5, beta = 2.0)
    
    expect_true(result > 0, info = "WLFC pseudocount should be positive")
})

test_that("compute_wlfc_pseudocounts returns bounded value", {
    # Result should typically be in [0, 1] range
    counts <- c(10, 20, 30, 15)
    result <- compute_wlfc_pseudocounts(counts, alpha = 1.5, beta = 2.0)
    
    expect_true(result >= 0 && result <= 1,
                info = "WLFC pseudocount should be in [0,1] range"
    )
})

test_that("compute_wlfc_pseudocounts works with different priors", {
    # Function should work with different prior strengths
    counts <- c(5, 10, 15, 20)
    
    # Weak prior (Jeffreys)
    result_weak <- compute_wlfc_pseudocounts(counts, alpha = 0.5, beta = 0.5)
    
    # Strong prior
    result_strong <- compute_wlfc_pseudocounts(counts, alpha = 10, beta = 10)
    
    # Both should return valid positive values
    expect_true(result_weak > 0)
    expect_true(result_strong > 0)
    expect_true(is.finite(result_weak))
    expect_true(is.finite(result_strong))
})

test_that("compute_wlfc_pseudocounts is consistent", {
    # Function should be deterministic
    counts <- c(10, 20, 30, 15)
    
    result1 <- compute_wlfc_pseudocounts(counts, alpha = 1.5, beta = 2.0)
    result2 <- compute_wlfc_pseudocounts(counts, alpha = 1.5, beta = 2.0)
    
    expect_equal(result1, result2)
})

test_that("compute_wlfc_pseudocounts with uniform vs non-uniform counts", {
    # Test with different abundance patterns
    counts_uniform <- c(50, 50, 50, 50)
    counts_skewed <- c(100, 20, 10, 20)
    
    result_uniform <- compute_wlfc_pseudocounts(counts_uniform, alpha = 1.5, beta = 2.0)
    result_skewed <- compute_wlfc_pseudocounts(counts_skewed, alpha = 1.5, beta = 2.0)
    
    # Both should be valid and comparable
    expect_true(result_uniform > 0)
    expect_true(result_skewed > 0)
    expect_true(is.finite(result_uniform))
    expect_true(is.finite(result_skewed))
})

test_that("compute_wlfc_pseudocounts with single sample", {
    # Should work with single count value
    counts <- c(25)
    result <- compute_wlfc_pseudocounts(counts, alpha = 1.5, beta = 2.0)
    
    expect_is(result, "numeric")
    expect_length(result, 1)
    expect_true(result > 0)
})

# =============================================================================
# Integration tests for empirical Bayes workflow
# =============================================================================

test_that("empirical Bayes workflow: fit prior then compute pseudocounts for realistic dataset", {
    # Full workflow with larger realistic dataset
    set.seed(123)
    n_genes <- 50
    n_samples <- 15
    counts_matrix <- matrix(
        rpois(n_genes * n_samples, lambda = 25),
        nrow = n_genes, ncol = n_samples
    )
    
    # Step 1: Fit prior
    prior_params <- fit_empirical_beta_prior(counts_matrix)
    expect_is(prior_params, "list")
    
    # Step 2: Compute pseudocounts for each gene
    pseudocounts <- apply(
        counts_matrix, 1,
        compute_wlfc_pseudocounts,
        prior_params$alpha,
        prior_params$beta
    )
    
    # Check results
    expect_is(pseudocounts, "numeric")
    expect_equal(length(pseudocounts), n_genes, info = "Should have one pseudocount per gene")
    expect_true(all(pseudocounts > 0), info = "All pseudocounts should be positive")
    expect_true(all(pseudocounts <= 1), info = "All pseudocounts should be <= 1")
})

test_that("empirical Bayes workflow with realistic dataset", {
    # Full workflow with synthetic but realistic data
    set.seed(456)
    n_genes <- 30
    n_samples <- 10
    readcounts_sim <- matrix(
        rpois(n_genes * n_samples, lambda = 30),
        nrow = n_genes, ncol = n_samples
    )
    
    # Fit prior
    prior_params <- fit_empirical_beta_prior(readcounts_sim)
    
    # Compute pseudocounts
    pseudocounts <- apply(
        readcounts_sim, 1,
        compute_wlfc_pseudocounts,
        prior_params$alpha,
        prior_params$beta
    )
    
    # Verify
    expect_equal(length(pseudocounts), nrow(readcounts_sim))
    expect_true(all(pseudocounts > 0))
    expect_true(all(is.finite(pseudocounts)))
})

test_that("empirical Bayes pseudocounts can be used in entropy calculation", {
    # Verify pseudocounts are compatible with calculate_tsallis_entropy
    set.seed(789)
    n_genes <- 20
    n_samples <- 8
    readcounts_test <- matrix(
        rpois(n_genes * n_samples, lambda = 25),
        nrow = n_genes, ncol = n_samples
    )
    
    # Get empirical Bayes pseudocounts
    suppressWarnings({
        prior_params <- fit_empirical_beta_prior(readcounts_test)
        pseudocounts <- apply(
            readcounts_test, 1,
            compute_wlfc_pseudocounts,
            prior_params$alpha,
            prior_params$beta
        )
        
        # Use in entropy calculation for first gene
        gene_idx <- 1
        entropy <- calculate_tsallis_entropy(
            readcounts_test[gene_idx, ],
            q = 1.0,
            pseudocount = pseudocounts[gene_idx],
            norm = TRUE
        )
        
        expect_is(entropy, "numeric")
        expect_true(is.finite(entropy))
        expect_true(entropy >= 0)
    })
})

test_that("empirical Bayes vs Jeffreys prior: WLFC computation", {
    # Test computation of WLFC vs default Jeffreys prior value
    suppressWarnings({
        counts <- c(5, 10, 15, 20)  # Normal counts
        
        prior_params <- fit_empirical_beta_prior(matrix(counts, nrow = 1))
        wlfc_pc <- compute_wlfc_pseudocounts(counts, 
                                             prior_params$alpha, 
                                             prior_params$beta)
        
        # WLFC should be a valid positive number
        expect_true(is.numeric(wlfc_pc))
        expect_true(wlfc_pc > 0)
        expect_true(is.finite(wlfc_pc))
    })
})

test_that("fit_empirical_beta_prior handles extremely small counts", {
    # Should handle edge case of very small count values
    counts <- matrix(c(1, 0, 0, 1, 0, 1), nrow = 2, ncol = 3)
    
    suppressWarnings({
        result <- fit_empirical_beta_prior(counts)
        
        expect_true(result$alpha > 0)
        expect_true(result$beta > 0)
    })
})

test_that("empirical Bayes parameters vary by count distribution", {
    # Different count distributions should yield valid priors
    suppressWarnings({
        counts_uniform <- matrix(rep(50, 12), nrow = 3, ncol = 4)
        counts_skewed <- matrix(c(100, 1, 1, 1, 95, 5, 0, 0, 90, 10, 0, 0),
                                nrow = 3, ncol = 4, byrow = TRUE)
        
        result_uniform <- fit_empirical_beta_prior(counts_uniform)
        result_skewed <- fit_empirical_beta_prior(counts_skewed)
        
        # Results should both return valid priors
        expect_true(result_uniform$alpha > 0)
        expect_true(result_uniform$beta > 0)
        expect_true(result_skewed$alpha > 0)
        expect_true(result_skewed$beta > 0)
    })
})
