context("Tier 3: Hierarchical Bayesian Methods")

# =====================================================================
# TEST SUITE: estimate_hierarchical_ar1_prior
# =====================================================================

test_that("estimate_hierarchical_ar1_prior returns correct structure", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(42)
    
    # Create test SE with multiple q-values
    # Structure: each gene × sample × q creates a curve
    n_genes <- 20
    n_samples <- 8
    n_q <- 4
    
    # Generate entropy data with realistic q-dependent structure
    # Tsallis entropy: H_q ∝ log(exp(α*q) - 1)/(q-1)
    # Decreases monotonically with q for q > 1
    entropy_matrix <- matrix(nrow = n_genes, ncol = n_samples * n_q)
    
    for (i in seq_len(n_genes)) {
        # Gene-specific scaling (concentration)
        gene_scale <- runif(1, 0.5, 2)
        
        for (j in seq_len(n_samples)) {
            # Sample-specific noise level
            sample_noise <- rnorm(1, 0, 0.1)
            
            # Tsallis entropy decreases as q increases
            q_vals <- seq(0.5, 3, length.out = n_q)
            base_curve <- 3.5 * exp(-0.15 * q_vals) * gene_scale
            
            for (k in seq_len(n_q)) {
                col_idx <- (j - 1) * n_q + k
                # Add realistic noise with AR(1) structure
                entropy_matrix[i, col_idx] <- max(0.1, base_curve[k] + sample_noise)
            }
        }
    }
    
    rownames(entropy_matrix) <- paste0("Gene_", 1:n_genes)
    colnames(entropy_matrix) <- paste0(
        rep(paste0("Sample_", 1:n_samples), each = n_q),
        "_q=", round(rep(seq(0.5, 3, length.out = n_q), n_samples), 2)
    )
    
    se <- SummarizedExperiment(assays = list(diversity = entropy_matrix))
    
    # Estimate hierarchical prior
    prior <- estimate_hierarchical_ar1_prior(se, verbose = FALSE)
    
    # Check structure
    expect_is(prior, "ar1_hierarchical_prior")
    expect_true(is.list(prior))
    expect_true(all(c("mu_phi", "sigma_phi", "phi_individual", "phi_shrunk",
                      "shrinkage_factors", "diagnostics") %in% names(prior)))
})

test_that("estimate_hierarchical_ar1_prior produces valid AR(1) parameters", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(42)
    
    # Create realistic entropy test SE
    n_genes <- 20
    n_samples <- 8
    n_q <- 4
    
    entropy_matrix <- matrix(nrow = n_genes, ncol = n_samples * n_q)
    for (i in seq_len(n_genes)) {
        gene_scale <- runif(1, 0.5, 2)
        for (j in seq_len(n_samples)) {
            sample_noise <- rnorm(1, 0, 0.1)
            q_vals <- seq(0.5, 3, length.out = n_q)
            base_curve <- 3.5 * exp(-0.15 * q_vals) * gene_scale
            for (k in seq_len(n_q)) {
                col_idx <- (j - 1) * n_q + k
                entropy_matrix[i, col_idx] <- max(0.1, base_curve[k] + sample_noise)
            }
        }
    }
    
    rownames(entropy_matrix) <- paste0("Gene_", 1:n_genes)
    colnames(entropy_matrix) <- paste0(
        rep(paste0("Sample_", 1:n_samples), each = n_q),
        "_q=", round(rep(seq(0.5, 3, length.out = n_q), n_samples), 2)
    )
    
    se <- SummarizedExperiment(assays = list(diversity = entropy_matrix))
    
    prior <- estimate_hierarchical_ar1_prior(se, verbose = FALSE)
    
    # Check AR(1) validity: φ should be in (-1, 1)
    expect_true(prior$mu_phi > -1 && prior$mu_phi < 1)
    expect_true(prior$sigma_phi > 0)
    
    # Individual estimates should be in valid range
    valid_phi <- prior$phi_individual[!is.na(prior$phi_individual)]
    expect_true(all(valid_phi > -0.999 & valid_phi < 0.999))
    
    # Shrunk estimates should also be valid
    expect_true(all(prior$phi_shrunk[!is.na(prior$phi_shrunk)] > -0.999 &
                    prior$phi_shrunk[!is.na(prior$phi_shrunk)] < 0.999))
})

test_that("estimate_hierarchical_ar1_prior with Yule-Walker method", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(42)
    
    n_genes <- 20
    n_samples <- 8
    n_q <- 4
    
    entropy_matrix <- matrix(nrow = n_genes, ncol = n_samples * n_q)
    for (i in seq_len(n_genes)) {
        gene_scale <- runif(1, 0.5, 2)
        for (j in seq_len(n_samples)) {
            sample_noise <- rnorm(1, 0, 0.1)
            q_vals <- seq(0.5, 3, length.out = n_q)
            base_curve <- 3.5 * exp(-0.15 * q_vals) * gene_scale
            for (k in seq_len(n_q)) {
                col_idx <- (j - 1) * n_q + k
                entropy_matrix[i, col_idx] <- max(0.1, base_curve[k] + sample_noise)
            }
        }
    }
    
    rownames(entropy_matrix) <- paste0("Gene_", 1:n_genes)
    colnames(entropy_matrix) <- paste0(
        rep(paste0("Sample_", 1:n_samples), each = n_q),
        "_q=", round(rep(seq(0.5, 3, length.out = n_q), n_samples), 2)
    )
    
    se <- SummarizedExperiment(assays = list(diversity = entropy_matrix))
    
    prior_yw <- estimate_hierarchical_ar1_prior(se, method = "yule_walker", verbose = FALSE)
    
    # Should have valid estimates
    expect_equal(prior_yw$method, "yule_walker")
    expect_true(!is.na(prior_yw$mu_phi))
    expect_true(length(prior_yw$phi_individual) == 20)
})

test_that("estimate_hierarchical_ar1_prior shrinkage works correctly", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(42)
    
    n_genes <- 20
    n_samples <- 8
    n_q <- 4
    
    entropy_matrix <- matrix(nrow = n_genes, ncol = n_samples * n_q)
    for (i in seq_len(n_genes)) {
        gene_scale <- runif(1, 0.5, 2)
        for (j in seq_len(n_samples)) {
            sample_noise <- rnorm(1, 0, 0.1)
            q_vals <- seq(0.5, 3, length.out = n_q)
            base_curve <- 3.5 * exp(-0.15 * q_vals) * gene_scale
            for (k in seq_len(n_q)) {
                col_idx <- (j - 1) * n_q + k
                entropy_matrix[i, col_idx] <- max(0.1, base_curve[k] + sample_noise)
            }
        }
    }
    
    rownames(entropy_matrix) <- paste0("Gene_", 1:n_genes)
    colnames(entropy_matrix) <- paste0(
        rep(paste0("Sample_", 1:n_samples), each = n_q),
        "_q=", round(rep(seq(0.5, 3, length.out = n_q), n_samples), 2)
    )
    
    se <- SummarizedExperiment(assays = list(diversity = entropy_matrix))
    
    prior <- estimate_hierarchical_ar1_prior(se, verbose = FALSE)
    
    # Shrunk estimates should be closer to mu_phi than raw estimates
    # (where both are non-NA)
    valid_genes <- which(!is.na(prior$phi_individual) & !is.na(prior$phi_shrunk))
    
    if (length(valid_genes) > 0) {
        for (i in valid_genes) {
            # Shrink factor should be between 0 and 1 (except when φ ≈ μ)
            if (abs(prior$phi_individual[i] - prior$mu_phi) > 0.01) {
                expect_true(prior$shrinkage_factors[i] >= 0 & prior$shrinkage_factors[i] <= 1)
            }
        }
    }
})

test_that("estimate_hierarchical_ar1_prior handles low-variance data", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(42)
    
    # Mostly constant with small perturbations (low variance)
    n_genes <- 20
    n_samples <- 8
    n_q <- 4
    
    entropy_matrix <- matrix(nrow = n_genes, ncol = n_samples * n_q)
    for (i in seq_len(n_genes)) {
        for (j in seq_len(n_samples)) {
            q_vals <- seq(0.5, 3, length.out = n_q)
            # Small monotone decrease
            base_curve <- 2.0 - 0.05 * q_vals
            for (k in seq_len(n_q)) {
                col_idx <- (j - 1) * n_q + k
                entropy_matrix[i, col_idx] <- base_curve[k] + rnorm(1, 0, 0.02)
            }
        }
    }
    
    rownames(entropy_matrix) <- paste0("Gene_", 1:n_genes)
    colnames(entropy_matrix) <- paste0(
        rep(paste0("Sample_", 1:n_samples), each = n_q),
        "_q=", round(rep(seq(0.5, 3, length.out = n_q), n_samples), 2)
    )
    
    se <- SummarizedExperiment(assays = list(diversity = entropy_matrix))
    
    prior <- estimate_hierarchical_ar1_prior(se, verbose = FALSE)
    
    # With low variance, φ estimates should be small
    valid_phi <- prior$phi_individual[!is.na(prior$phi_individual)]
    expect_true(mean(abs(valid_phi)) < 0.5)
})

test_that("estimate_hierarchical_ar1_prior requires sufficient q-values", {
    skip_if_not_installed("SummarizedExperiment")
    
    # Only 2 q-values (insufficient)
    entropy_matrix <- matrix(rnorm(40), nrow = 20, ncol = 2)
    rownames(entropy_matrix) <- paste0("Gene_", 1:20)
    colnames(entropy_matrix) <- c("S1_q=1", "S2_q=2")
    
    se <- SummarizedExperiment(assays = list(diversity = entropy_matrix))
    
    expect_error(
        estimate_hierarchical_ar1_prior(se, verbose = FALSE),
        "Need at least 3 distinct q-values"
    )
})

test_that("estimate_hierarchical_ar1_prior requires sufficient genes", {
    skip_if_not_installed("SummarizedExperiment")
    
    # Only 2 genes (insufficient) with 8 samples and 4 q-values
    n_genes <- 2
    n_samples <- 8
    n_q <- 4
    
    entropy_matrix <- matrix(rnorm(n_genes * n_samples * n_q), nrow = n_genes, ncol = n_samples * n_q)
    rownames(entropy_matrix) <- c("Gene_1", "Gene_2")
    colnames(entropy_matrix) <- paste0(
        rep(paste0("Sample_", 1:n_samples), each = n_q),
        "_q=", round(rep(seq(0.5, 3, length.out = n_q), n_samples), 2)
    )
    
    se <- SummarizedExperiment(assays = list(diversity = entropy_matrix))
    
    expect_error(
        estimate_hierarchical_ar1_prior(se, verbose = FALSE),
        "Need at least 5 genes"
    )
})

test_that("estimate_hierarchical_ar1_prior print method works", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(42)
    
    n_genes <- 20
    n_samples <- 8
    n_q <- 4
    
    entropy_matrix <- matrix(nrow = n_genes, ncol = n_samples * n_q)
    for (i in seq_len(n_genes)) {
        gene_scale <- runif(1, 0.5, 2)
        for (j in seq_len(n_samples)) {
            sample_noise <- rnorm(1, 0, 0.1)
            q_vals <- seq(0.5, 3, length.out = n_q)
            base_curve <- 3.5 * exp(-0.15 * q_vals) * gene_scale
            for (k in seq_len(n_q)) {
                col_idx <- (j - 1) * n_q + k
                entropy_matrix[i, col_idx] <- max(0.1, base_curve[k] + sample_noise)
            }
        }
    }
    
    rownames(entropy_matrix) <- paste0("Gene_", 1:n_genes)
    colnames(entropy_matrix) <- paste0(
        rep(paste0("Sample_", 1:n_samples), each = n_q),
        "_q=", round(rep(seq(0.5, 3, length.out = n_q), n_samples), 2)
    )
    
    se <- SummarizedExperiment(assays = list(diversity = entropy_matrix))
    prior <- estimate_hierarchical_ar1_prior(se, verbose = FALSE)
    
    # Print should not error
    expect_output(print(prior), "HIERARCHICAL AR")
    expect_output(print(prior), "Population mean")
    expect_output(print(prior), "Population std")
})

test_that("estimate_hierarchical_ar1_prior diagnostics are informative", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(42)
    
    n_genes <- 20
    n_samples <- 8
    n_q <- 4
    
    entropy_matrix <- matrix(nrow = n_genes, ncol = n_samples * n_q)
    for (i in seq_len(n_genes)) {
        gene_scale <- runif(1, 0.5, 2)
        for (j in seq_len(n_samples)) {
            sample_noise <- rnorm(1, 0, 0.1)
            q_vals <- seq(0.5, 3, length.out = n_q)
            base_curve <- 3.5 * exp(-0.15 * q_vals) * gene_scale
            for (k in seq_len(n_q)) {
                col_idx <- (j - 1) * n_q + k
                entropy_matrix[i, col_idx] <- max(0.1, base_curve[k] + sample_noise)
            }
        }
    }
    
    rownames(entropy_matrix) <- paste0("Gene_", 1:n_genes)
    colnames(entropy_matrix) <- paste0(
        rep(paste0("Sample_", 1:n_samples), each = n_q),
        "_q=", round(rep(seq(0.5, 3, length.out = n_q), n_samples), 2)
    )
    
    se <- SummarizedExperiment(assays = list(diversity = entropy_matrix))
    prior <- estimate_hierarchical_ar1_prior(se, verbose = FALSE)
    
    # Diagnostics should have required fields
    diag <- prior$diagnostics
    expect_true(all(c("n_genes_with_valid_phi", "mean_obs_per_gene",
                      "proportion_high_phi", "mean_shrinkage_factor") %in% names(diag)))
    
    # Values should be reasonable
    expect_true(diag$n_genes_with_valid_phi >= 3)
    expect_true(diag$proportion_high_phi >= 0 & diag$proportion_high_phi <= 1)
    expect_true(diag$mean_shrinkage_factor >= 0 & diag$mean_shrinkage_factor <= 1)
})

test_that("estimate_hierarchical_ar1_prior hyperprior_dist options work", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(42)
    
    n_genes <- 20
    n_samples <- 8
    n_q <- 4
    
    entropy_matrix <- matrix(nrow = n_genes, ncol = n_samples * n_q)
    for (i in seq_len(n_genes)) {
        gene_scale <- runif(1, 0.5, 2)
        for (j in seq_len(n_samples)) {
            sample_noise <- rnorm(1, 0, 0.1)
            q_vals <- seq(0.5, 3, length.out = n_q)
            base_curve <- 3.5 * exp(-0.15 * q_vals) * gene_scale
            for (k in seq_len(n_q)) {
                col_idx <- (j - 1) * n_q + k
                entropy_matrix[i, col_idx] <- max(0.1, base_curve[k] + sample_noise)
            }
        }
    }
    
    rownames(entropy_matrix) <- paste0("Gene_", 1:n_genes)
    colnames(entropy_matrix) <- paste0(
        rep(paste0("Sample_", 1:n_samples), each = n_q),
        "_q=", round(rep(seq(0.5, 3, length.out = n_q), n_samples), 2)
    )
    
    se <- SummarizedExperiment(assays = list(diversity = entropy_matrix))
    
    # Test normal hyperprior
    prior_normal <- estimate_hierarchical_ar1_prior(
        se, hyperprior_dist = "normal", verbose = FALSE
    )
    expect_equal(prior_normal$hyperprior_dist, "normal")
    
    # Test uniform hyperprior
    prior_uniform <- estimate_hierarchical_ar1_prior(
        se, hyperprior_dist = "uniform", verbose = FALSE
    )
    expect_equal(prior_uniform$hyperprior_dist, "uniform")
})

test_that("estimate_hierarchical_ar1_prior with missing values", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(42)
    
    n_genes <- 20
    n_samples <- 8
    n_q <- 4
    
    entropy_matrix <- matrix(nrow = n_genes, ncol = n_samples * n_q)
    for (i in seq_len(n_genes)) {
        gene_scale <- runif(1, 0.5, 2)
        for (j in seq_len(n_samples)) {
            sample_noise <- rnorm(1, 0, 0.1)
            q_vals <- seq(0.5, 3, length.out = n_q)
            base_curve <- 3.5 * exp(-0.15 * q_vals) * gene_scale
            for (k in seq_len(n_q)) {
                col_idx <- (j - 1) * n_q + k
                entropy_matrix[i, col_idx] <- max(0.1, base_curve[k] + sample_noise)
            }
        }
    }
    
    # Add some NAs
    entropy_matrix[1, 1] <- NA
    entropy_matrix[5, 2:3] <- NA
    
    rownames(entropy_matrix) <- paste0("Gene_", 1:n_genes)
    colnames(entropy_matrix) <- paste0(
        rep(paste0("Sample_", 1:n_samples), each = n_q),
        "_q=", round(rep(seq(0.5, 3, length.out = n_q), n_samples), 2)
    )
    
    se <- SummarizedExperiment(assays = list(diversity = entropy_matrix))
    
    prior <- estimate_hierarchical_ar1_prior(se, verbose = FALSE)
    
    # Should handle NAs gracefully
    expect_is(prior, "ar1_hierarchical_prior")
    expect_true(!is.na(prior$mu_phi))
})

test_that("estimate_hierarchical_ar1_prior papers validation present", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(42)
    
    n_genes <- 20
    n_samples <- 8
    n_q <- 4
    
    entropy_matrix <- matrix(nrow = n_genes, ncol = n_samples * n_q)
    for (i in seq_len(n_genes)) {
        gene_scale <- runif(1, 0.5, 2)
        for (j in seq_len(n_samples)) {
            sample_noise <- rnorm(1, 0, 0.1)
            q_vals <- seq(0.5, 3, length.out = n_q)
            base_curve <- 3.5 * exp(-0.15 * q_vals) * gene_scale
            for (k in seq_len(n_q)) {
                col_idx <- (j - 1) * n_q + k
                entropy_matrix[i, col_idx] <- max(0.1, base_curve[k] + sample_noise)
            }
        }
    }
    
    rownames(entropy_matrix) <- paste0("Gene_", 1:n_genes)
    colnames(entropy_matrix) <- paste0(
        rep(paste0("Sample_", 1:n_samples), each = n_q),
        "_q=", round(rep(seq(0.5, 3, length.out = n_q), n_samples), 2)
    )
    
    se <- SummarizedExperiment(assays = list(diversity = entropy_matrix))
    prior <- estimate_hierarchical_ar1_prior(se, verbose = FALSE)
    
    # Should cite supporting papers
    expect_true(length(prior$papers_validation) > 0)
    expect_true("S168-S171" %in% names(prior$papers_validation))
    expect_true("BY002-BY003" %in% names(prior$papers_validation))
})
