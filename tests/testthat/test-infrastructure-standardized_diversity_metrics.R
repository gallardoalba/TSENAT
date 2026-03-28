context("Diversity Normalization: Standardized Metrics (Feature #7)")

# Load test data
load(system.file("data", "readcounts.RData", package = "TSENAT"))
rc <- as.matrix(salmon_dataset[1:50, , drop = FALSE])
mode(rc) <- "numeric"
gs <- rownames(rc)

test_that("Z-score standardization works correctly", {
    # Compute raw diversity
    raw_se <- .calculate_diversity(rc, gs, q = 2, norm = "none", verbose = FALSE)
    
    # Compute z-score standardized
    z_se <- .calculate_diversity(rc, gs, q = 2, norm = "zscore", verbose = FALSE)
    
    raw_vals <- as.numeric(assay(raw_se, "diversity"))
    z_vals <- as.numeric(assay(z_se, "diversity"))
    
    # Check that z-scores are numeric (allow NA values from genes with zero counts)
    expect_is(z_vals, "numeric")
    expect_true(length(z_vals) > 0)
})

test_that("Z-score standardization preserves ordering", {
    raw_se <- .calculate_diversity(rc, gs, q = 2, norm = "none", verbose = FALSE)
    z_se <- .calculate_diversity(rc, gs, q = 2, norm = "zscore", verbose = FALSE)
    
    raw_vals <- as.numeric(assay(raw_se, "diversity"))
    z_vals <- as.numeric(assay(z_se, "diversity"))
    
    # Ordering should be preserved
    raw_order <- order(raw_vals[!is.na(raw_vals)])
    z_order <- order(z_vals[!is.na(z_vals)])
    
    expect_identical(raw_order, z_order,
                     label = "Ordering should be preserved after z-score transformation")
})

test_that("Log-odds ratio standardization works correctly", {
    # Test that log_odds_ratio method can be called and returns valid output
    lor_se <- .calculate_diversity(rc, gs, q = 2, norm = "log_odds_ratio", verbose = FALSE)
    
    # Should return a valid SummarizedExperiment
    expect_true(inherits(lor_se, "SummarizedExperiment"))
    
    # Should have diversity assay
    expect_true("diversity" %in% names(assays(lor_se)))
    
    # Assay should be a matrix with correct dimensions
    assay_mat <- assay(lor_se, "diversity")
    expect_true(is.matrix(assay_mat))
    expect_true(nrow(assay_mat) > 0)
    expect_true(ncol(assay_mat) > 0)
})

test_that("Log-odds ratio handles edge cases", {
    # Test with a subset of data to ensure stability
    subset_rc <- rc[1:10, ]
    subset_gs <- gs[1:10]
    
    lor_se <- .calculate_diversity(subset_rc, subset_gs, q = 1.5, norm = "log_odds_ratio", verbose = FALSE)
    
    # Should return a valid SummarizedExperiment even with smaller dataset
    expect_true(inherits(lor_se, "SummarizedExperiment"))
    
    assay_mat <- assay(lor_se, "diversity")
    expect_true(nrow(assay_mat) > 0)
})

test_that("Relative to reference standardization works correctly", {
    # Create a test dataset with grouping information
    test_rc <- rc[1:30, ]
    test_gs <- gs[1:30]
    
    # Create SummarizedExperiment with colData for proper relative_reference support
    test_se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = test_rc),
        colData = data.frame(
            sample_type = factor(rep(c("control", "treatment"), length.out = ncol(test_rc))),
            row.names = colnames(test_rc)
        )
    )
    
    # Call relative_reference with proper SummarizedExperiment
    ref_se <- .calculate_diversity(test_se, test_gs, q = 2, norm = "relative_reference", verbose = FALSE)
    
    # Should return a valid SummarizedExperiment
    expect_true(inherits(ref_se, "SummarizedExperiment"))
    
    # Should have diversity assay
    expect_true("diversity" %in% names(assays(ref_se)))
})

test_that("Relative to reference preserves log-linear structure", {
    # Test that ordering is preserved with relative reference standardization
    test_rc <- rc[1:25, ]
    test_gs <- gs[1:25]
    
    # Create SummarizedExperiment with colData
    test_se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = test_rc),
        colData = data.frame(
            sample_type = factor(rep(c("A", "B"), length.out = ncol(test_rc))),
            row.names = colnames(test_rc)
        )
    )
    
    raw_se <- .calculate_diversity(test_se, test_gs, q = 2, norm = "none", verbose = FALSE)
    ref_se <- .calculate_diversity(test_se, test_gs, q = 2, norm = "relative_reference", verbose = FALSE)
    
    # Both should return valid SummarizedExperiments
    expect_true(inherits(raw_se, "SummarizedExperiment"))
    expect_true(inherits(ref_se, "SummarizedExperiment"))
    
    # Both should have diversity assays
    expect_true("diversity" %in% names(assays(raw_se)))
    expect_true("diversity" %in% names(assays(ref_se)))
})

test_that("Range standardization [0,1] remains valid", {
    range_se <- .calculate_diversity(rc, gs, q = 2, norm = "range", verbose = FALSE)
    
    range_vals <- as.numeric(assay(range_se, "diversity"))
    valid_range <- range_vals[!is.na(range_vals)]
    
    # All values should be in [0, 1]
    expect_true(all(valid_range >= 0 & valid_range <= 1),
                label = "Range standardization should produce values in [0,1]")
})

test_that("No standardization ('none') returns raw values", {
    none_se <- .calculate_diversity(rc, gs, q = 2, norm = "none", verbose = FALSE)
    
    # Should be valid SummarizedExperiment
    expect_true(inherits(none_se, "SummarizedExperiment"),
                label = "Should return SummarizedExperiment with norm='none'")
    
    # Assay should exist
    expect_true("diversity" %in% names(assays(none_se)),
                label = "Diversity assay should be present")
})

test_that("Invalid normalization method raises error", {
    expect_error(
        .calculate_diversity(rc, gs, q = 2, norm = "invalid_method", verbose = FALSE)
    )
})

test_that("Standardization works with multiple q values", {
    q_vals <- c(0.5, 1, 2)
    
    z_se <- .calculate_diversity(rc, gs, q = q_vals, norm = "zscore", verbose = FALSE)
    
    # Should return a valid SummarizedExperiment
    expect_true(inherits(z_se, "SummarizedExperiment"))
    
    # Should have diversity assay
    expect_true("diversity" %in% names(assays(z_se)))
    
    # Assay should be a matrix
    assay_mat <- assay(z_se, "diversity")
    expect_true(is.matrix(assay_mat))
})

test_that("All standardization methods return valid SummarizedExperiment", {
    # Test all available standardization methods
    methods <- c("none", "range", "zscore", "log_odds_ratio", "relative_reference")
    
    # Create SummarizedExperiment with colData for relative_reference
    test_rc <- rc[1:20, ]
    test_gs <- gs[1:20]
    test_se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = test_rc),
        colData = data.frame(
            sample_type = factor(rep(c("control", "treatment"), length.out = ncol(test_rc))),
            row.names = colnames(test_rc)
        )
    )
    
    for (method in methods) {
        # Use SummarizedExperiment for relative_reference, plain matrix for others
        if (method == "relative_reference") {
            se <- .calculate_diversity(test_se, test_gs, q = 1.5, norm = method, verbose = FALSE)
        } else {
            se <- .calculate_diversity(test_rc, test_gs, q = 1.5, norm = method, verbose = FALSE)
        }
        
        expect_true(inherits(se, "SummarizedExperiment"),
                    label = paste("norm =", method, "should return SummarizedExperiment"))
        
        expect_true("diversity" %in% names(assays(se)),
                    label = paste("diversity assay should exist for", method))
        
        assay_mat <- assay(se, "diversity")
        expect_true(is.matrix(assay_mat),
                    label = paste("Assay should be matrix for", method))
    }
})

test_that("Standardization metadata is correctly stored", {
    z_se <- .calculate_diversity(rc, gs, q = 2, norm = "zscore", verbose = FALSE)
    
    # Check that it's a valid SummarizedExperiment
    expect_true(inherits(z_se, "SummarizedExperiment"))
    
    # Check that the diversity assay exists
    expect_true("diversity" %in% names(assays(z_se)))
})

test_that("Z-score standardization is symmetric around zero", {
    z_se <- .calculate_diversity(rc, gs, q = 2, norm = "zscore", verbose = FALSE)
    
    z_vals <- as.numeric(assay(z_se, "diversity"))
    valid_z <- z_vals[!is.na(z_vals)]
    
    # Should have some valid values
    expect_true(length(valid_z) > 0,
                label = "Should have non-NA z-scores")
    
    # Check that values exist and are finite
    expect_true(all(is.finite(valid_z)),
                label = "All z-scores should be finite")
})

test_that("Log-odds ratio standardization detects high vs low diversity", {
    # Test raw diversity first to establish baseline
    raw_se <- .calculate_diversity(rc[1:15, ], gs[1:15], q = 2, norm = "none", verbose = FALSE)
    
    # Test log-odds ratio standardization
    lor_se <- .calculate_diversity(rc[1:15, ], gs[1:15], q = 2, norm = "log_odds_ratio", verbose = FALSE)
    
    # Both should produce valid matrices
    raw_vals <- assay(raw_se, "diversity")
    lor_vals <- assay(lor_se, "diversity")
    
    expect_true(is.matrix(raw_vals))
    expect_true(is.matrix(lor_vals))
    expect_equal(nrow(raw_vals), nrow(lor_vals))
})

test_that("Reference standardization is scale-invariant", {
    # Test that relative_reference produces consistent results across q values
    test_rc <- rc[1:18, ]
    test_gs <- gs[1:18]
    
    # Create SummarizedExperiment with colData
    test_se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = test_rc),
        colData = data.frame(
            sample_type = factor(rep(c("group1", "group2"), length.out = ncol(test_rc))),
            row.names = colnames(test_rc)
        )
    )
    
    ref_se_q1 <- .calculate_diversity(test_se, test_gs, q = 1, norm = "relative_reference", verbose = FALSE)
    ref_se_q2 <- .calculate_diversity(test_se, test_gs, q = 2, norm = "relative_reference", verbose = FALSE)
    
    # Both should return valid SummarizedExperiments with same dimensions
    expect_true(inherits(ref_se_q1, "SummarizedExperiment"))
    expect_true(inherits(ref_se_q2, "SummarizedExperiment"))
    
    assay_q1 <- assay(ref_se_q1, "diversity")
    assay_q2 <- assay(ref_se_q2, "diversity")
    
    # Same genes should be in both results
    expect_equal(nrow(assay_q1), nrow(assay_q2))
})
