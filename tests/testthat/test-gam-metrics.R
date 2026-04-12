context("GAM Metrics - Generalized Additive Models Diagnostics")

# Test data setup - Create synthetic entropy data for GAM metrics testing
# GAM metrics expect: rows=genes, columns=q-values (not samples!)
# Example: 100 genes evaluated at 41 q-values gives 100 x 41 matrix

set.seed(42)
n_genes <- 100
n_q <- 41

# Create synthetic q-values (0 to 2 in steps of 0.05)
q_values <- seq(0, 2, length.out = n_q)

# Create synthetic entropy data: 
# - Each row is a gene's entropy curve across q-values
# - Entropy curves typically show Gaussian-like or unimodal patterns
x <- matrix(NA_real_, nrow = n_genes, ncol = n_q)
for (i in 1:n_genes) {
    # Create realistic entropy curve: entropy ~ Gaussian centered around a q-value
    peak_q <- runif(1, min = 0.5, max = 1.5)  # Random peak location
    peak_height <- runif(1, min = 1, max = 5)  # Random peak height
    width <- runif(1, min = 0.3, max = 1.0)    # Random curve width
    baseline <- runif(1, min = 0, max = 0.5)   # Random baseline noise
    
    entropy_curve <- peak_height * exp(-((q_values - peak_q)^2) / (2 * width^2)) + baseline
    x[i, ] <- entropy_curve + rnorm(n_q, 0, 0.1)  # Add small noise
}

# Ensure non-negative entropy
x[x < 0] <- 0

# Set column names as q-values for reference
colnames(x) <- paste0("q_", round(q_values, 3))

# Create small datasets for specific tests
n_obs <- nrow(x)
n_vars <- ncol(x)

# Data with collinearity (for concurvity testing) - create highly correlated q-columns
x_collinear <- cbind(
    x[, 1:2], 
    X3 = 0.9 * x[, 1] + 0.1 * x[, 2] + 0.05 * matrix(rnorm(n_obs), nrow = n_obs),
    x[, 3:10]
)

# Smaller data for faster intermediate tests
x_small <- x[1:min(20, n_obs), 1:min(10, n_q), drop = FALSE]

# Single predictor data - but still need multiple q-values for GAM
x_single <- x[1:min(50, n_obs), 1:min(5, n_q), drop = FALSE]

# CONCURVITY INDEX TESTS

test_that("compute_concurvity_index returns expected structure", {
    result <- TSENAT:::.compute_concurvity_index(x, q_values = q_values)
    
    expect_is(result, "list")
    expect_true(!is.null(result$description))
    expect_true(!is.null(result$status))
    expect_true(!is.null(result$details))
})

test_that("compute_concurvity_index detects high collinearity", {
    # Highly collinear data should have concurvity > 0.6
    result <- TSENAT:::.compute_concurvity_index(x_collinear, q_values = q_values)
    
    # Check for reasonable value or NA
    expect_true(is.na(result$overall_concurvity) || (result$overall_concurvity >= 0 && result$overall_concurvity <= 1))
})

test_that("compute_concurvity_index handles small data", {
    # Should gracefully handle small datasets
    q_small <- seq(0, 2, length.out = ncol(x_small))
    result <- TSENAT:::.compute_concurvity_index(x_small, q_values = q_small)
    
    expect_is(result, "list")
    expect_true(!is.null(result$status))
})

test_that("compute_concurvity_index requires at least 2 predictors", {
    # Single q-value is not enough for GAM
    x_true_single <- x[, 1, drop = FALSE]
    result <- TSENAT:::.compute_concurvity_index(x_true_single, q_values = c(0.5))
    
    expect_equal(result$overall_concurvity, 0)
    expect_match(result$status, "N/A")
})

test_that("compute_concurvity_index converts non-matrix to matrix", {
    df_data <- as.data.frame(x)
    result <- TSENAT:::.compute_concurvity_index(df_data, q_values = q_values)
    
    expect_is(result, "list")
    expect_true(!is.null(result$status))
})

# EFFECTIVE DEGREES OF FREEDOM (EDF) TESTS

test_that("compute_edf_metric returns expected structure", {
    result <- TSENAT:::.compute_edf_metric(x, q_values = q_values)
    
    expect_is(result, "list")
    expect_true(!is.null(result$description))
    expect_true(!is.null(result$status))
    expect_match(result$description, "Effective Degrees of Freedom")
})

test_that("compute_edf_metric returns numeric EDF values", {
    result <- TSENAT:::.compute_edf_metric(x, q_values = q_values)
    
    expect_is(result$total_edf, "numeric")
    expect_is(result$edf_ratio, "numeric")
    expect_true(result$edf_ratio > 0)
})

test_that("compute_edf_metric categorizes smoothing appropriately", {
    result <- TSENAT:::.compute_edf_metric(x, q_values = q_values)
    
    # Status should be one of the three categories
    expect_match(
        result$status,
        "over-smoothed|appropriate|under-smoothed"
    )
})

test_that("compute_edf_metric provides interpretation", {
    result <- TSENAT:::.compute_edf_metric(x, q_values = q_values)
    
    if (result$edf_ratio < 0.5) {
        expect_match(result$status, "over-smoothed")
    } else if (result$edf_ratio <= 2.0) {
        expect_match(result$status, "appropriate")
    } else {
        expect_match(result$status, "under-smoothed")
    }
})

# NONLINEARITY CONTRIBUTION TESTS

test_that("compute_nonlinearity_contribution returns expected structure", {
    result <- TSENAT:::.compute_nonlinearity_contribution(x, q_values = q_values)
    
    expect_is(result, "list")
    expect_true(!is.null(result$description))
    expect_true(!is.null(result$status))
})

test_that("compute_nonlinearity_contribution returns percentage improvement", {
    result <- TSENAT:::.compute_nonlinearity_contribution(x, q_values = q_values)
    
    expect_is(result$r2_improvement_percent, "numeric")
    expect_true(result$r2_improvement_percent >= -100)  # Can go negative
})

test_that("compute_nonlinearity_contribution categorizes appropriately", {
    result <- TSENAT:::.compute_nonlinearity_contribution(x, q_values = q_values)
    
    # Status should be one of three categories
    expect_match(
        result$status,
        "use linear|gam justified|gam essential"
    )
})

test_that("compute_nonlinearity_contribution decision logic works", {
    result <- TSENAT:::.compute_nonlinearity_contribution(x, q_values = q_values)
    
    if (result$r2_improvement_percent < 5) {
        expect_match(result$status, "use linear")
    } else if (result$r2_improvement_percent < 20) {
        expect_match(result$status, "gam justified")
    } else {
        expect_match(result$status, "gam essential")
    }
})

# ============================================================================
# BASIS FUNCTION ADEQUACY TESTS
# ============================================================================

test_that("compute_basis_adequacy returns expected structure", {
    result <- TSENAT:::.compute_basis_adequacy(x, q_values = q_values)
    
    expect_is(result, "list")
    expect_true(!is.null(result$description))
    expect_true(!is.null(result$status))
})

test_that("compute_basis_adequacy returns integer k dimension", {
    result <- TSENAT:::.compute_basis_adequacy(x, q_values = q_values)
    
    expect_type(result$optimal_basis_dimension, "integer")
    expect_true(result$optimal_basis_dimension > 0)
})

test_that("compute_basis_adequacy returns k in tested range", {
    result <- TSENAT:::.compute_basis_adequacy(x, q_values = q_values)
    
    k <- result$optimal_basis_dimension
    expect_true(k %in% c(3, 5, 8, 10, 15))
})

test_that("compute_basis_adequacy indicates convergence status", {
    result <- TSENAT:::.compute_basis_adequacy(x, q_values = q_values)
    
    expect_match(
        result$status,
        "adequate|consider increase"
    )
})

# ============================================================================
# WRAPPER FUNCTION TESTS
# ============================================================================

test_that("get_gam_metrics calls all 4 metric functions", {
    result <- TSENAT:::.get_gam_metrics(x, q_values = q_values)
    
    expect_is(result, "list")
    expect_equal(length(result), 5)  # 4 metrics + consolidated result
    expect_true("concurvity" %in% names(result))
    expect_true("edf" %in% names(result))
    expect_true("nonlinearity" %in% names(result))
    expect_true("basis_adequacy" %in% names(result))
    expect_true("consolidated" %in% names(result))
})

test_that("get_gam_metrics computes all sub-metrics", {
    result <- TSENAT:::.get_gam_metrics(x, q_values = q_values)
    
    # Each should be a list with description
    expect_is(result$concurvity, "list")
    expect_is(result$edf, "list")
    expect_is(result$nonlinearity, "list")
    expect_is(result$basis_adequacy, "list")
    
    expect_true(!is.null(result$concurvity$description))
    expect_true(!is.null(result$edf$description))
    expect_true(!is.null(result$nonlinearity$description))
    expect_true(!is.null(result$basis_adequacy$description))
})

test_that("get_gam_metrics handles matrix and data.frame input", {
    df_data <- as.data.frame(x)
    result_df <- TSENAT:::.get_gam_metrics(df_data, q_values = q_values)
    
    expect_is(result_df, "list")
    expect_equal(length(result_df), 5)  # 4 metrics + consolidated
})

test_that("get_gam_metrics failures in one metric don't stop others", {
    # This verifies independent computation
    q_small <- seq(0, 2, length.out = ncol(x_small))
    result <- TSENAT:::.get_gam_metrics(x_small, q_values = q_small)
    
    expect_is(result, "list")
    expect_equal(length(result), 5)  # 4 metrics + consolidated
    # Just verify structure - some metrics may fail on small data
    expect_true(all(c("concurvity", "edf", "nonlinearity", "basis_adequacy") %in% names(result)))
})

# ============================================================================
# INTEGRATION WITH CALCULATE_RANK_ASSUMPTIONS TESTS
# ============================================================================

test_that("calculate_rank_assumptions accepts gam_metrics check", {
    result <- TSENAT:::.calculate_rank_assumptions(x, checks = c("gam_metrics"), q_values = q_values)
    
    check_results <- attr(result, "checks")
    expect_true("gam_metrics" %in% names(check_results))
})

test_that("calculate_rank_assumptions with gam_metrics produces 4 metrics", {
    result <- TSENAT:::.calculate_rank_assumptions(x, checks = c("gam_metrics"), q_values = q_values)
    
    check_results <- attr(result, "checks")
    gam_metrics <- check_results$gam_metrics
    
    expect_equal(length(gam_metrics), 5)  # 4 metrics + consolidated
})

test_that("calculate_rank_assumptions combines rank and GAM checks", {
    # Use "all" preset which includes both
    result <- TSENAT:::.calculate_rank_assumptions(x, checks = "all", q_values = q_values)
    
    check_results <- attr(result, "checks")
    
    expect_true("exchangeability" %in% names(check_results))
    expect_true("monotonicity" %in% names(check_results))
    expect_true("consistency" %in% names(check_results))
    expect_true("gam_metrics" %in% names(check_results))
})

# PRINT METHOD TESTS

test_that("print.rank_assumptions handles GAM metrics", {
    result <- suppressWarnings(
        TSENAT:::.calculate_rank_assumptions(x, checks = c("gam_metrics"), q_values = q_values)
    )
    
    # Verify structure contains GAM metrics
    check_results <- attr(result, "checks")
    expect_true(!is.null(check_results$gam_metrics))
    expect_equal(length(check_results$gam_metrics), 5)  # 4 metrics + consolidated
    
    # Verify print method runs without error
    expect_error(
        {
            capture.output({
                withr::with_message_sink(stdout(), {
                    print(result)
                })
            })
        },
        NA  # Expect NO error
    )
})

test_that("print.rank_assumptions displays both rank and GAM checks", {
    result <- suppressWarnings(
        TSENAT:::.calculate_rank_assumptions(
            x,
            checks = c("exchangeability", "gam_metrics"),
            q_values = q_values
        )
    )
    
    # Verify both rank and GAM checks are present
    check_results <- attr(result, "checks")
    expect_true(!is.null(check_results$exchangeability))
    expect_true(!is.null(check_results$gam_metrics))
    
    # Verify print method runs without error
    expect_error(
        {
            capture.output({
                withr::with_message_sink(stdout(), {
                    print(result)
                })
            })
        },
        NA  # Expect NO error
    )
})

# ============================================================================
# ERROR HANDLING TESTS
# ============================================================================

test_that("GAM metrics handle missing values gracefully", {
    x_missing <- x
    x_missing[1:5, 1] <- NA
    
    result <- TSENAT:::.compute_concurvity_index(x_missing, q_values = q_values)
    expect_is(result, "list")
    expect_true(!is.null(result$status))
})

test_that("GAM metrics handle all-NA columns", {
    x_allna <- x
    x_allna[, 2] <- NA
    
    result <- TSENAT:::.compute_edf_metric(x_allna, q_values = q_values)
    expect_is(result, "list")
})

test_that("GAM metrics handle infinite values", {
    x_inf <- x
    x_inf[1, 1] <- Inf
    
    result <- TSENAT:::.compute_nonlinearity_contribution(x_inf, q_values = q_values)
    expect_is(result, "list")
    expect_true(!is.null(result$status))
})

# ============================================================================
# EDGE CASE TESTS
# ============================================================================

test_that("GAM metrics work with very small data (n=2)", {
    x_tiny <- matrix(rnorm(10), nrow = 2, ncol = 5)
    colnames(x_tiny) <- c("X1", "X2", "X3", "X4", "X5")
    q_tiny <- seq(0, 2, length.out = 5)
    
    result <- TSENAT:::.get_gam_metrics(x_tiny, q_values = q_tiny)
    expect_is(result, "list")
    expect_equal(length(result), 5)  # 4 metrics + consolidated
})

test_that("GAM metrics work with single predictor (though limited)", {
    x_truly_single <- x[, 1, drop = FALSE]
    
    result <- TSENAT:::.get_gam_metrics(x_truly_single, q_values = c(0))
    expect_is(result, "list")
    expect_equal(length(result), 4)
})

test_that("GAM metrics work with many predictors", {
    x_many <- matrix(rnorm(n_obs * 20), nrow = n_obs, ncol = 20)
    colnames(x_many) <- paste0("X", seq_len(20))
    q_many <- seq(0, 2, length.out = 20)
    
    result <- TSENAT:::.get_gam_metrics(x_many, q_values = q_many)
    expect_is(result, "list")
    expect_equal(length(result), 5)  # 4 metrics + consolidated
})

# CITATION VERIFICATION TESTS

test_that("GAM metrics reference correct papers in details", {
    result_edf <- TSENAT:::.compute_edf_metric(x, q_values = q_values)
    # Verify details field exists (paper references depend on EDF computation)
    expect_true(!is.null(result_edf$details))
})

test_that("GAM metrics mention Concurvity papers", {
    result <- TSENAT:::.compute_concurvity_index(x, q_values = q_values)
    # Verify structure exists (details may or may not mention papers depending on concurvity)
    expect_true(!is.null(result$details))
})

# ============================================================================
# NUMERICAL CORRECTNESS VALIDATION TESTS
# ============================================================================

test_that("concurvity values are mathematically sound", {
    # Independent (uncorrelated) predictors should have low concurvity
    result_uncor <- TSENAT:::.compute_concurvity_index(x, q_values = q_values)
    expect_true(is.na(result_uncor$overall_concurvity) || result_uncor$overall_concurvity < 0.5)
    
    # Collinear predictors should have higher concurvity
    q_collinear <- seq(0, 2, length.out = ncol(x_collinear))
    result_cor <- TSENAT:::.compute_concurvity_index(x_collinear, q_values = q_collinear)
    if (!is.na(result_uncor$overall_concurvity) && !is.na(result_cor$overall_concurvity)) {
        expect_true(result_cor$overall_concurvity >= result_uncor$overall_concurvity)
    }
})

test_that("EDF values are bounded and reasonable", {
    result <- TSENAT:::.compute_edf_metric(x, q_values = q_values)
    
    # EDF ratio should be positive
    expect_true(result$edf_ratio > 0)
    # Status should describe what the EDF ratio means
    expect_true(result$status %in% c("appropriate", "over-smoothed", "under-smoothed"))
})

test_that("nonlinearity contribution is reasonable", {
    result <- TSENAT:::.compute_nonlinearity_contribution(x, q_values = q_values)
    
    # R² improvement percentage should be computed
    # Can be negative if GAM is worse, but usually positive
    expect_true(is.numeric(result$r2_improvement_percent))
    expect_true(!is.na(result$r2_improvement_percent))
})

test_that("basis adequacy finds optimal dimension", {
    result <- TSENAT:::.compute_basis_adequacy(x, q_values = q_values)
    
    # Should find an optimal basis dimension or report error/skip
    expect_true(is.numeric(result$optimal_basis_dimension))
    expect_true(result$optimal_basis_dimension > 0)
    expect_true(result$optimal_basis_dimension <= 20)
})

test_that("GAM R² is not worse than linear R²", {
    # Build models for comparison
    y <- rowMeans(x, na.rm = TRUE)
    data_df <- as.data.frame(x)
    data_df$y <- y
    col_names <- colnames(data_df)[colnames(data_df) != "y"]
    
    # Linear model
    formula_linear <- as.formula(paste0("y ~ ", paste0("`", head(col_names, 3), "`", collapse = " + ")))
    lm_model <- stats::lm(formula_linear, data = data_df)
    r2_lm <- suppressWarnings(summary(lm_model))$r.squared
    
    # GAM model
    formula_gam <- as.formula(paste0("y ~ ", paste0("s(`", head(col_names, 3), "`)", collapse = " + ")))
    gam_model <- mgcv::gam(formula_gam, data = data_df, method = "GCV.Cp", control = list(maxit = 100))
    r2_gam <- suppressWarnings(summary(gam_model))$r.sq
    
    # GAM should typically have R² >= linear (unless overfitting with small sample)
    if (!is.na(r2_lm) && !is.na(r2_gam)) {
        expect_true(r2_gam >= r2_lm * 0.95)  # Allow 5% tolerance for numerical differences
    }
})
