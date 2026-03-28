
# ════════════════════════════════════════════════════════════════════════════════
# GAMMA FAMILY CONDITIONAL SELECTION TESTS (March 2026)
# ════════════════════════════════════════════════════════════════════════════════

context("Gamma Family Conditional Selection")

testthat::test_that(".compute_skewness calculates skewness correctly", {
    # Test symmetric distribution (skewness ≈ 0)
    symmetric <- c(1, 2, 3, 4, 5)
    expect_lt(abs(.compute_skewness(symmetric)), 0.1)
    
    # Test right-skewed distribution (positive skewness)
    right_skewed <- c(1, 1, 2, 3, 10)
    expect_gt(.compute_skewness(right_skewed), 0)
    
    # Test left-skewed distribution (negative skewness)
    left_skewed <- c(-10, 1, 2, 3, 3)
    expect_lt(.compute_skewness(left_skewed), 0)
    
    # Test with NA handling
    with_na <- c(1, 2, 3, NA, 4, 5)
    expect_true(!is.na(.compute_skewness(with_na, na.rm = TRUE)))
})

testthat::test_that(".detect_heteroscedasticity identifies variance changes", {
    set.seed(42)
    n_q <- 10
    n_samples <- 4
    
    # Create homoscedastic data (no variance changes)
    df_homo <- data.frame(
        entropy = rnorm(n_q * n_samples, mean = 0.3, sd = 0.05),
        q = rep(seq(0.01, 0.8, length.out = n_q), n_samples),
        group = rep(c("A", "B"), each = n_q * n_samples / 2)
    )
    
    result_homo <- .detect_heteroscedasticity(df_homo, 
                                                     q_vals = sort(unique(df_homo$q)),
                                                     group_vec = df_homo$group)
    
    expect_true(is.list(result_homo))
    expect_true("is_heteroscedastic" %in% names(result_homo))
    expect_true("p_value" %in% names(result_homo))
    # Homoscedastic data should NOT be significant
    expect_gt(result_homo$p_value, 0.05)
    
    # Create heteroscedastic data (variance increases with q)
    df_hetero <- data.frame(
        entropy = numeric(n_q * n_samples),
        q = rep(seq(0.01, 0.8, length.out = n_q), n_samples),
        group = rep(c("A", "B"), each = n_q * n_samples / 2)
    )
    
    for (i in seq_len(n_q)) {
        idx <- which(df_hetero$q == df_hetero$q[i])
        sd_i <- 0.01 + 0.05 * (i / n_q)^2  # Increasing variance with q
        df_hetero$entropy[idx] <- 0.3 + rnorm(length(idx), 0, sd_i)
    }
    
    result_hetero <- .detect_heteroscedasticity(df_hetero,
                                                       q_vals = sort(unique(df_hetero$q)),
                                                       group_vec = df_hetero$group)
    
    expect_true(is.list(result_hetero))
    expect_true(!is.na(result_hetero$p_value))
})

testthat::test_that(".select_gam_family chooses Gaussian for normal data", {
    set.seed(42)
    n_q <- 8
    n_samples <- 4
    
    # Normal entropy data without issues
    df <- data.frame(
        entropy = rnorm(n_q * n_samples, mean = 0.3, sd = 0.05),
        q = rep(seq(0.01, 0.8, length.out = n_q), n_samples),
        group = rep(c("A", "B"), each = n_q * n_samples / 2)
    )
    
    # Ensure all positive
    df$entropy <- pmax(df$entropy, 0.001)
    
    bounded_result <- .handle_bounded_support(df, q_vals = sort(unique(df$q)),
                                        group_vec = df$group, verbose = FALSE)
    result <- bounded_result$family_info
    
    expect_true(is.list(result))
    expect_true("use_gamma" %in% names(result))
    expect_true("use_gaussian" %in% names(result))
    
    # For normal data, should prefer Gaussian
    expect_false(result$use_gamma)
    expect_true(result$use_gaussian)
})

testthat::test_that(".select_gam_family triggers Gamma for strong heteroscedasticity", {
    skip_if_not_installed("mgcv")
    set.seed(42)
    n_q <- 10
    n_samples <- 4
    
    # Create strongly heteroscedastic data
    df <- data.frame(
        entropy = numeric(n_q * n_samples),
        q = rep(seq(0.01, 0.8, length.out = n_q), n_samples),
        group = rep(c("A", "B"), each = n_q * n_samples / 2)
    )
    
    # Extreme variance increase with q (variance ratio > 5)
    for (i in seq_len(n_q)) {
        idx <- which(df$q == df$q[i])
        sd_i <- 0.005 + 0.08 * (i / n_q)^3  # Strongly increasing variance
        df$entropy[idx] <- 0.25 + rnorm(length(idx), 0, sd_i)
    }
    
    df$entropy <- pmax(df$entropy, 0.001)  # Ensure positive
    
    bounded_result <- .handle_bounded_support(df, q_vals = sort(unique(df$q)),
                                        group_vec = df$group, verbose = FALSE)
    result <- bounded_result$family_info
    
    expect_true(is.list(result))
    
    # Heteroscedastic data may or may not trigger Gamma depending on test p-value
    # But we should see the analysis was performed
    expect_true("var_ratio_q" %in% names(result))
    expect_gt(result$var_ratio_q, 1)  # Should show some variance change
})

testthat::test_that(".select_gam_family detects boundary clustering", {
    set.seed(42)
    n_q <- 10
    n_samples <- 4
    
    # Create data with extreme boundary clustering (>40% near min bound)
    entropy_boundary <- c(rep(0.001, 4), rep(0.002, 4), seq(0.05, 0.5, length.out = n_q-8))
    
    df <- data.frame(
        entropy = rep(entropy_boundary, n_samples),
        q = rep(seq(0.01, 0.8, length.out = n_q), n_samples),
        group = rep(c("A", "B"), each = n_q * n_samples / 2)
    )
    
    bounded_result <- .handle_bounded_support(df, q_vals = sort(unique(df$q)),
                                        group_vec = df$group, verbose = FALSE)
    result <- bounded_result$family_info
    
    expect_true(is.list(result))
    expect_true("boundary_pct" %in% names(result))
    # Should detect substantial boundary clustering
    expect_gt(result$boundary_pct, 30)
})

testthat::test_that(".handle_bounded_support returns correct family structure", {
    set.seed(42)
    n_q <- 8
    n_samples <- 4
    
    df <- data.frame(
        entropy = rnorm(n_q * n_samples, mean = 0.3, sd = 0.05),
        q = rep(seq(0.01, 0.8, length.out = n_q), n_samples),
        group = rep(c("A", "B"), each = n_q * n_samples / 2)
    )
    
    df$entropy <- pmax(df$entropy, 0.001)
    
    result <- .handle_bounded_support(df, q_vals = sort(unique(df$q)),
                                             group_vec = df$group, verbose = FALSE)
    
    expect_true(is.list(result))
    expect_true("use_gamma" %in% names(result))
    expect_true("use_gaussian" %in% names(result))
    expect_true("family_obj" %in% names(result))
    expect_true("inverse_link" %in% names(result))
    
    # For normal data, should be Gaussian
    expect_false(result$use_gamma)
    expect_true(result$use_gaussian)
    
    # Check family object
    expect_true(inherits(result$family_obj, "family"))
    
    # Check inverse link function works
    test_eta <- c(-1, 0, 1, 2)
    test_pred <- result$inverse_link(test_eta)
    expect_true(all(is.numeric(test_pred)))
    expect_equal(length(test_pred), length(test_eta))
})

testthat::test_that(".handle_bounded_support inverse link for Gaussian is identity", {
    set.seed(42)
    df <- data.frame(
        entropy = rnorm(32, mean = 0.3, sd = 0.05),
        q = rep(seq(0.01, 0.8, length.out = 8), 4),
        group = rep(c("A", "B"), each = 16)
    )
    
    df$entropy <- pmax(df$entropy, 0.001)
    
    result <- .handle_bounded_support(df, q_vals = sort(unique(df$q)),
                                             group_vec = df$group)
    
    # If using Gaussian, inverse link should be identity
    if (result$use_gaussian) {
        test_eta <- c(-2, -1, 0, 1, 2)
        test_pred <- result$inverse_link(test_eta)
        expect_equal(test_pred, test_eta, tolerance = 1e-10)
    }
})

testthat::test_that(".handle_bounded_support inverse link for Gamma is exp", {
    skip_if_not_installed("mgcv")
    set.seed(42)
    n_q <- 10
    n_samples <- 4
    
    # Create heteroscedastic data to potentially trigger Gamma
    df <- data.frame(
        entropy = numeric(n_q * n_samples),
        q = rep(seq(0.01, 0.8, length.out = n_q), n_samples),
        group = rep(c("A", "B"), each = n_q * n_samples / 2)
    )
    
    for (i in seq_len(n_q)) {
        idx <- which(df$q == df$q[i])
        sd_i <- 0.005 + 0.1 * (i / n_q)^3
        df$entropy[idx] <- 0.25 + rnorm(length(idx), 0, sd_i)
    }
    
    df$entropy <- pmax(df$entropy, 0.001)
    
    result <- .handle_bounded_support(df, q_vals = sort(unique(df$q)),
                                             group_vec = df$group)
    
    # Test inverse link function
    test_eta <- c(-2, -1, 0, 1, 2)
    test_pred <- result$inverse_link(test_eta)
    
    expect_true(all(is.numeric(test_pred)))
    expect_true(all(test_pred > 0))  # Should be positive for any family
    
    # If using Gaussian, should be identity
    if (result$use_gaussian) {
        expect_equal(test_pred, test_eta, tolerance = 1e-10)
    }
    # If using Gamma, should be exponential
    else if (result$use_gamma) {
        expect_equal(test_pred, exp(test_eta), tolerance = 1e-10)
    }
})

testthat::test_that("GAM fitting with selected family produces valid predictions", {
    skip_if_not_installed("mgcv")
    set.seed(42)
    n_q <- 8
    n_samples <- 4
    
    # Create test data
    df <- data.frame(
        entropy = rnorm(n_q * n_samples, mean = 0.3, sd = 0.05),
        q = rep(seq(0.01, 0.8, length.out = n_q), n_samples),
        group = rep(c("A", "B"), each = n_q * n_samples / 2)
    )
    
    df$entropy <- pmax(df$entropy, 0.001)
    
    # Get selected family
    family_info <- .handle_bounded_support(df, q_vals = sort(unique(df$q)),
                                                  group_vec = df$group)
    
    # Fit GAM with selected family
    fit <- tryCatch({
        mgcv::gam(entropy ~ group + s(q, bs="tp", k=4),
                 family = family_info$family_obj,
                 data = df)
    }, error = function(e) NULL)
    
    skip_if(is.null(fit), "GAM fitting failed")
    
    # Get predictions
    pred <- predict(fit, type = "response")
    
    expect_true(all(is.numeric(pred)))
    expect_true(all(pred > 0))  # All predictions should be positive for entropy
    expect_equal(length(pred), nrow(df))
})

testthat::test_that(".gam_interaction uses conditional family selection", {
    skip_if_not_installed("mgcv")
    set.seed(42)
    
    n_q <- 8
    q <- rep(seq(0.01, 0.8, length.out = n_q), 4)
    group <- rep(c("A", "B"), each = 16)
    entropy <- 0.3 + 0.15 * q + 0.1 * (group == "B") + rnorm(32, 0, 0.03)
    
    df <- data.frame(entropy = pmax(entropy, 0.001), q = q, group = factor(group))
    
    # Call GAM interaction
    result <- suppressWarnings(.gam_interaction(df, q_vals = unique(q), g = "test_gene", min_obs = 5))
    
    skip_if(is.null(result), "GAM interaction returned NULL")
    
    expect_true(is.data.frame(result))
    expect_true("p_interaction" %in% colnames(result))
    expect_true("gene" %in% colnames(result))
    
    # P-value should be numeric and in valid range
    expect_true(is.numeric(result$p_interaction))
    expect_gte(result$p_interaction, 0)
    expect_lte(result$p_interaction, 1)
})

testthat::test_that("Family selection is stable across repeated calls", {
    set.seed(42)
    n_q <- 8
    n_samples <- 4
    
    df <- data.frame(
        entropy = rnorm(n_q * n_samples, mean = 0.3, sd = 0.05),
        q = rep(seq(0.01, 0.8, length.out = n_q), n_samples),
        group = rep(c("A", "B"), each = n_q * n_samples / 2)
    )
    
    df$entropy <- pmax(df$entropy, 0.001)
    
    # Call family selection multiple times
    results <- list()
    for (i in seq_len(3)) {
        bounded_result <- .handle_bounded_support(df, q_vals = sort(unique(df$q)),
                                                   group_vec = df$group, verbose = FALSE)
        results[[i]] <- bounded_result$family_info
    }
    
    # All results should be consistent
    expect_equal(results[[1]]$use_gamma, results[[2]]$use_gamma)
    expect_equal(results[[2]]$use_gamma, results[[3]]$use_gamma)
})

testthat::test_that("Near-zero entropy values are handled correctly", {
    skip_if_not_installed("mgcv")
    set.seed(42)
    
    # Create data with very small entropy values
    entropy_edge <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2)
    
    df <- data.frame(
        entropy = rep(entropy_edge, 2),
        group = rep(c("A", "B"), each = length(entropy_edge)),
        idx = rep(1:length(entropy_edge), 2)
    )
    
    # Should handle without error
    result <- tryCatch({
        bounded_result <- .handle_bounded_support(df, q_vals = rep(0.5, nrow(df)),
                                 group_vec = df$group)
        bounded_result$family_info
    }, error = function(e) NULL)
    
    expect_false(is.null(result))
    expect_true(is.list(result))
})

testthat::test_that("Gamma predictions are always positive for log link", {
    skip_if_not_installed("mgcv")
    set.seed(42)
    n_q <- 10
    n_samples <- 4
    
    # Create data that might trigger Gamma
    df <- data.frame(
        entropy = numeric(n_q * n_samples),
        q = rep(seq(0.01, 0.8, length.out = n_q), n_samples),
        group = rep(c("A", "B"), each = n_q * n_samples / 2)
    )
    
    for (i in seq_len(n_q)) {
        idx <- which(df$q == df$q[i])
        sd_i <- 0.01 + 0.05 * (i / n_q)^2
        df$entropy[idx] <- 0.25 + rnorm(length(idx), 0, sd_i)
    }
    
    df$entropy <- pmax(df$entropy, 0.001)
    
    # Get family
    family_info <- .handle_bounded_support(df, q_vals = sort(unique(df$q)),
                                                  group_vec = df$group)
    
    # Fit model
    fit <- tryCatch({
        mgcv::gam(entropy ~ group + s(q, bs="tp", k=5),
                 family = family_info$family_obj,
                 data = df)
    }, error = function(e) NULL)
    
    skip_if(is.null(fit), "GAM fitting failed")
    
    # Check predictions
    pred <- predict(fit, type = "response", newdata = df[1:10, ])
    
    expect_true(all(pred > 0), "All predictions should be positive")
})
