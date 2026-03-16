context("calculate_divergence: Bootstrap Divergence CI Computation")

library(TSENAT)
library(SummarizedExperiment)

test_that("calculate_divergence works with basic SE input", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(42)
    
    # Create a basic SummarizedExperiment
    counts <- matrix(rpois(20 * 8, lambda = 100), nrow = 20, ncol = 8)
    rownames(counts) <- paste0("Gene_", 1:20)
    colnames(counts) <- paste0("Sample_", 1:8)
    
    metadata <- data.frame(
        group = factor(c(rep("Control", 4), rep("Treatment", 4))),
        row.names = colnames(counts)
    )
    
    se <- SummarizedExperiment(
        assays = list(counts = counts),
        colData = metadata
    )
    
    # Test with bootstrap=TRUE (default)
    result_bootstrap <- calculate_divergence(
        se = se,
        bootstrap = TRUE,
        nboot = 100,
        control_group = "Control",
        progress = FALSE
    )
    
    # Check return type (should be SummarizedExperiment)
    expect_is(result_bootstrap, "SummarizedExperiment")
    
    # Check rowData has required columns
    rd <- rowData(result_bootstrap)
    expect_true("gene_name" %in% colnames(rd))
    expect_true("estimate" %in% colnames(rd))
    expect_true("lower_ci" %in% colnames(rd))
    expect_true("upper_ci" %in% colnames(rd))
    
    # Check number of genes processed (should be 20, all genes in SE)
    expect_equal(nrow(result_bootstrap), 20)
    
    # Check estimates are computed
    expect_false(all(is.na(rd$estimate)))
    expect_false(all(is.na(rd$lower_ci)))
})

test_that("calculate_divergence bootstrap parameter works correctly", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(43)
    
    counts <- matrix(rpois(20 * 8, lambda = 100), nrow = 20, ncol = 8)
    rownames(counts) <- paste0("Gene_", 1:20)
    colnames(counts) <- paste0("Sample_", 1:8)
    metadata <- data.frame(group = factor(c(rep("A", 4), rep("B", 4))),
                           row.names = colnames(counts))
    se <- SummarizedExperiment(assays = list(counts = counts), colData = metadata)
    
    # Point estimates only (bootstrap=FALSE)
    result_point <- calculate_divergence(
        se = se,
        bootstrap = FALSE,
        control_group = "A",
        progress = FALSE
    )
    
    # With bootstrap (bootstrap=TRUE)
    result_boot <- calculate_divergence(
        se = se,
        bootstrap = TRUE,
        nboot = 50,
        control_group = "A",
        progress = FALSE
    )
    
    # Verify both return SummarizedExperiment
    expect_is(result_point, "SummarizedExperiment")
    expect_is(result_boot, "SummarizedExperiment")
    
    # Both should return estimates in rowData
    rd_point <- rowData(result_point)
    rd_boot <- rowData(result_boot)
    
    expect_false(all(is.na(rd_point$estimate)))
    expect_false(all(is.na(rd_boot$estimate)))
    
    # Point estimates should have NA CIs
    expect_true(all(is.na(rd_point$lower_ci)))
    expect_true(all(is.na(rd_point$upper_ci)))
    
    # Bootstrap should have CIs
    expect_false(all(is.na(rd_boot$lower_ci)))
    expect_false(all(is.na(rd_boot$upper_ci)))
})

test_that("calculate_divergence input validation", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(44)
    
    counts <- matrix(rpois(20 * 8, lambda = 100), nrow = 20, ncol = 8)
    rownames(counts) <- paste0("Gene_", 1:20)
    colnames(counts) <- paste0("Sample_", 1:8)
    metadata <- data.frame(group = factor(c(rep("A", 4), rep("B", 4))),
                           row.names = colnames(counts))
    se <- SummarizedExperiment(assays = list(counts = counts), colData = metadata)
    
    # Invalid SE (not SummarizedExperiment)
    expect_error(
        calculate_divergence(se = list()),
        "must be a SummarizedExperiment"
    )
    
    # Invalid bootstrap parameter
    expect_error(
        calculate_divergence(se = se, bootstrap = "yes"),
        "bootstrap must be a logical"
    )
})

test_that("calculate_divergence handles parallel processing", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(45)
    
    counts <- matrix(rpois(20 * 8, lambda = 100), nrow = 20, ncol = 8)
    rownames(counts) <- paste0("Gene_", 1:20)
    colnames(counts) <- paste0("Sample_", 1:8)
    metadata <- data.frame(group = factor(c(rep("A", 4), rep("B", 4))),
                           row.names = colnames(counts))
    se <- SummarizedExperiment(assays = list(counts = counts), colData = metadata)
    
    # Sequential (nthreads=1) - uses default nboot=1000, so no warning expected
    result_seq <- calculate_divergence(
        se = se,
        nthreads = 1,
        control_group = "A",
        progress = FALSE
    )
    
    expect_is(result_seq, "SummarizedExperiment")
    expect_false(all(is.na(rowData(result_seq)$estimate)))
})
test_that("calculate_divergence auto-detects paired samples", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(46)
    
    # Create SE with paired_samples column (expected for TSENAT readcounts)
    counts <- matrix(rpois(20 * 8, lambda = 100), nrow = 20, ncol = 8)
    rownames(counts) <- paste0("Gene_", 1:20)
    colnames(counts) <- c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8")
    
    # Metadata with paired_samples column (A-D pairs, matched normal-tumor)
    metadata <- data.frame(
        condition = factor(c(rep("Normal", 4), rep("Tumor", 4))),
        paired_samples = c("A", "B", "C", "D",  # Paired group IDs
                          "A", "B", "C", "D"),   # Same IDs for paired samples
        row.names = colnames(counts)
    )
    
    se <- SummarizedExperiment(
        assays = list(counts = counts),
        colData = metadata
    )
    
    # Test with bootstrap=TRUE (triggers auto-detection)
    result <- calculate_divergence(
        se = se,
        bootstrap = TRUE,
        nboot = 100,
        group_col = "condition",
        control_group = "Normal",
        progress = FALSE
    )
    
    # Should return valid SummarizedExperiment
    expect_is(result, "SummarizedExperiment")
    
    # Should have estimates computed
    rd <- rowData(result)
    expect_false(all(is.na(rd$estimate)))
    
    # Should have CIs from bootstrap
    expect_false(all(is.na(rd$lower_ci)))
    expect_false(all(is.na(rd$upper_ci)))
})

test_that("calculate_divergence works without paired_samples column", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(47)
    
    # SE without paired_samples column (should use independent bootstrap)
    counts <- matrix(rpois(20 * 8, lambda = 100), nrow = 20, ncol = 8)
    rownames(counts) <- paste0("Gene_", 1:20)
    colnames(counts) <- paste0("Sample_", 1:8)
    
    metadata <- data.frame(
        group = factor(c(rep("Control", 4), rep("Treatment", 4))),
        row.names = colnames(counts)
    )
    
    se <- SummarizedExperiment(
        assays = list(counts = counts),
        colData = metadata
    )
    
    # Test with bootstrap=TRUE (no pairing detected)
    result <- calculate_divergence(
        se = se,
        bootstrap = TRUE,
        nboot = 100,
        control_group = "Control",
        progress = FALSE
    )
    
    # Should work fine with independent bootstrap
    expect_is(result, "SummarizedExperiment")
    rd <- rowData(result)
    expect_false(all(is.na(rd$estimate)))
})

test_that(".detect_pair_ids correctly identifies paired structures", {
    skip_if_not_installed("SummarizedExperiment")
    
    # Test 1: paired_samples column detected
    counts <- matrix(rpois(20 * 6, lambda = 100), nrow = 20, ncol = 6)
    colnames(counts) <- paste0("Sample_", 1:6)
    
    metadata <- data.frame(
        paired_samples = c("Pair_A", "Pair_B", "Pair_C", "Pair_A", "Pair_B", "Pair_C"),
        row.names = colnames(counts)
    )
    
    se <- SummarizedExperiment(assays = list(counts = counts), colData = metadata)
    
    detected <- TSENAT:::.detect_pair_ids(se)
    
    expect_equal(detected$num_pairs, 3)
    expect_equal(detected$column_name, "paired_samples")
    expect_is(detected$pair_ids, "character")
    expect_equal(length(detected$pair_ids), 6)
    
    # Test 2: No paired column (should return NULL)
    metadata_no_pairs <- data.frame(
        group = factor(c(rep("A", 3), rep("B", 3))),
        row.names = colnames(counts)
    )
    
    se_no_pairs <- SummarizedExperiment(assays = list(counts = counts), colData = metadata_no_pairs)
    detected_no_pairs <- TSENAT:::.detect_pair_ids(se_no_pairs)
    
    expect_equal(detected_no_pairs$num_pairs, 0)
    expect_null(detected_no_pairs$pair_ids)
})

test_that("calculate_divergence normalization modes work correctly", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(42)
    
    # Create a basic SummarizedExperiment
    counts <- matrix(rpois(20 * 8, lambda = 100), nrow = 20, ncol = 8)
    rownames(counts) <- paste0("Gene_", 1:20)
    colnames(counts) <- paste0("Sample_", 1:8)
    
    metadata <- data.frame(
        group = factor(c(rep("Control", 4), rep("Treatment", 4))),
        row.names = colnames(counts)
    )
    
    se <- SummarizedExperiment(
        assays = list(counts = counts),
        colData = metadata
    )
    
    # Test norm="none" (no normalization)
    result_none <- calculate_divergence(
        se = se,
        bootstrap = FALSE,
        norm = "none",
        control_group = "Control",
        progress = FALSE
    )
    
    expect_is(result_none, "SummarizedExperiment")
    expect_equal(metadata(result_none)$normalization, "none")
    
    # Test norm="range" (range normalization [0,1])
    result_range <- calculate_divergence(
        se = se,
        bootstrap = FALSE,
        norm = "range",
        control_group = "Control",
        progress = FALSE
    )
    
    expect_is(result_range, "SummarizedExperiment")
    expect_equal(metadata(result_range)$normalization, "range")
    
    # Test norm="zscore" (z-score standardization)
    result_zscore <- calculate_divergence(
        se = se,
        bootstrap = FALSE,
        norm = "zscore",
        control_group = "Control",
        progress = FALSE
    )
    
    expect_is(result_zscore, "SummarizedExperiment")
    expect_equal(metadata(result_zscore)$normalization, "zscore")
    
    # Test backward compatibility: norm=TRUE should equal norm="range"
    result_true <- calculate_divergence(
        se = se,
        bootstrap = FALSE,
        norm = TRUE,
        control_group = "Control",
        progress = FALSE
    )
    
    expect_equal(metadata(result_true)$normalization, "range")
    
    # Test backward compatibility: norm=FALSE should equal norm="none"
    result_false <- calculate_divergence(
        se = se,
        bootstrap = FALSE,
        norm = FALSE,
        control_group = "Control",
        progress = FALSE
    )
    
    expect_equal(metadata(result_false)$normalization, "none")
})

test_that("calculate_divergence all normalization modes are supported", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(42)
    
    # Create a basic SummarizedExperiment
    counts <- matrix(rpois(20 * 8, lambda = 100), nrow = 20, ncol = 8)
    rownames(counts) <- paste0("Gene_", 1:20)
    colnames(counts) <- paste0("Sample_", 1:8)
    
    metadata <- data.frame(
        group = factor(c(rep("Control", 4), rep("Treatment", 4))),
        row.names = colnames(counts)
    )
    
    se <- SummarizedExperiment(
        assays = list(counts = counts),
        colData = metadata
    )
    
    # Test all 5 normalization modes
    norm_modes <- c("none", "range", "zscore", "log_odds_ratio", "relative_reference")
    results_list <- list()
    
    for (norm_mode in norm_modes) {
        result <- calculate_divergence(
            se = se,
            bootstrap = FALSE,
            norm = norm_mode,
            control_group = "Control",
            progress = FALSE
        )
        
        results_list[[norm_mode]] <- result
        expect_equal(metadata(result)$normalization, norm_mode)
        expect_is(result, "SummarizedExperiment")
        expect_true(nrow(result) > 0)
    }
})

test_that("calculate_divergence normalization produces valid ranges", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(42)
    
    counts <- matrix(rpois(20 * 8, lambda = 100), nrow = 20, ncol = 8)
    rownames(counts) <- paste0("Gene_", 1:20)
    colnames(counts) <- paste0("Sample_", 1:8)
    
    metadata <- data.frame(
        group = factor(c(rep("Control", 4), rep("Treatment", 4))),
        row.names = colnames(counts)
    )
    
    se <- SummarizedExperiment(
        assays = list(counts = counts),
        colData = metadata
    )
    
    # Test range normalization produces values in [0,1]
    result_range <- calculate_divergence(
        se = se,
        bootstrap = FALSE,
        norm = "range",
        control_group = "Control",
        progress = FALSE
    )
    
    estimates_range <- rowData(result_range)$estimate
    valid_estimates <- estimates_range[!is.na(estimates_range)]
    
    if (length(valid_estimates) > 0) {
        expect_true(all(valid_estimates >= -1e-6), info = "Range: values >= 0")
        expect_true(all(valid_estimates <= 1 + 1e-6), info = "Range: values <= 1")
    }
})

test_that("calculate_divergence norm parameter validation", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(42)
    
    counts <- matrix(rpois(20 * 8, lambda = 100), nrow = 20, ncol = 8)
    rownames(counts) <- paste0("Gene_", 1:20)
    colnames(counts) <- paste0("Sample_", 1:8)
    
    metadata <- data.frame(
        group = factor(c(rep("Control", 4), rep("Treatment", 4))),
        row.names = colnames(counts)
    )
    
    se <- SummarizedExperiment(
        assays = list(counts = counts),
        colData = metadata
    )
    
    # Test invalid norm value raises error
    expect_error(
        calculate_divergence(
            se = se,
            bootstrap = FALSE,
            norm = "invalid_mode",
            control_group = "Control",
            progress = FALSE
        )
    )
    
    # Test valid character values all work
    for (valid_mode in c("none", "range", "zscore", "log_odds_ratio", "relative_reference")) {
        result <- calculate_divergence(
            se = se,
            bootstrap = FALSE,
            norm = valid_mode,
            control_group = "Control",
            progress = FALSE
        )
        expect_equal(metadata(result)$normalization, valid_mode)
    }
})

# =====================================================================
# Bayesian Credible Intervals (Tier 2 Integration)
# =====================================================================

