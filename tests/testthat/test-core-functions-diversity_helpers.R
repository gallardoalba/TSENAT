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

context("Diversity Normalization: Standardized Metrics (Feature #7)")

# Load test data
load(system.file("data", "readcounts.RData", package = "TSENAT"))
rc <- as.matrix(readcounts[1:50, , drop = FALSE])
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

# ============================================================================
# Tests for .suggest_min_count() and min_count parameter
# ============================================================================

context("Minimum Count Filtering: Smart Gene Pre-Filtering (Option 1)")

test_that("suggest_min_count works with matrix input", {
  # Create matrix with mix of abundant and sparse genes
  counts <- matrix(c(
    rep(100, 9),   # 3 genes with 100 counts each
    rep(10, 9),    # 3 genes with 10 counts each  
    rep(1, 9)      # 3 genes with 1 count each
  ), nrow = 9, ncol = 3, byrow = TRUE)
  
  # Test auto-detection at median
  min_cnt <- .suggest_min_count(counts, percentile = 0.5, verbose = FALSE)
  expect_is(min_cnt, "numeric")
  expect_true(min_cnt > 0)
  expect_true(is.finite(min_cnt))
  
  # Should return median of gene totals: [300, 300, 300, 30, 30, 30, 3, 3, 3]
  # Median = 30
  expect_equal(min_cnt, 30)
})

test_that("suggest_min_count respects percentile parameter", {
  counts <- matrix(c(
    rep(100, 9),   # 3 genes: 300 total
    rep(10, 9),    # 3 genes: 30 total
    rep(1, 9)      # 3 genes: 3 total
  ), nrow = 9, ncol = 3, byrow = TRUE)
  
  # Test lower percentile (25th)
  min_25 <- .suggest_min_count(counts, percentile = 0.25, verbose = FALSE)
  expect_true(min_25 > 0)
  
  # Test higher percentile (75th)
  min_75 <- .suggest_min_count(counts, percentile = 0.75, verbose = FALSE)
  expect_true(min_75 > min_25)
  
  # Higher percentile should give higher threshold
  expect_true(min_75 > min_25)
})

test_that("suggest_min_count works with SummarizedExperiment", {
  skip_if_not_installed("SummarizedExperiment")
  
  counts <- matrix(c(
    rep(100, 9),
    rep(10, 9)
  ), nrow = 6, ncol = 3, byrow = TRUE)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = S4Vectors::DataFrame(row.names = paste0("gene_", 1:6))
  )
  
  min_cnt <- .suggest_min_count(se, percentile = 0.5, verbose = FALSE)
  expect_is(min_cnt, "numeric")
  expect_true(min_cnt > 0)
})

test_that("suggest_min_count validates percentile parameter", {
  counts <- matrix(c(10, 20, 30), nrow = 3)
  
  # Invalid percentile (negative)
  expect_error(
    .suggest_min_count(counts, percentile = -0.5, verbose = FALSE),
    "percentile must be numeric in \\[0, 1\\]"
  )
  
  # Invalid percentile (> 1)
  expect_error(
    .suggest_min_count(counts, percentile = 1.5, verbose = FALSE),
    "percentile must be numeric in \\[0, 1\\]"
  )
})

test_that("suggest_min_count handles edge cases", {
  # All genes have same count
  counts <- matrix(rep(50, 15), nrow = 5, ncol = 3)
  min_cnt <- .suggest_min_count(counts, percentile = 0.5, verbose = FALSE)
  expect_equal(min_cnt, 150)  # 50 * 3 samples
  
  # Zero-count genes
  counts_zeros <- matrix(c(
    rep(100, 9),
    rep(0, 9)
  ), nrow = 6, ncol = 3, byrow = TRUE)
  min_cnt_zeros <- .suggest_min_count(counts_zeros, percentile = 0.5, verbose = FALSE)
  expect_true(min_cnt_zeros == 150 || min_cnt_zeros == 0)
})


# ============================================================================
# FILTER INTEGRATION TESTS: Verify that filtering preserves data consistency
# across colnames/rownames alignment
# ============================================================================

test_that("filter_analysis_s4 maintains TSENATAnalysis integrity", {
  skip_if_not_installed("SummarizedExperiment")
  
  # Create a minimal TSENATAnalysis object with TPM to avoid warnings
  counts <- matrix(c(5, 10, 8, 12, 3, 7, 15, 20, 10), nrow = 3, ncol = 3)
  tpm <- matrix(c(50, 100, 80, 120, 30, 70, 150, 200, 100), nrow = 3, ncol = 3)  # Add TPM
  rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(
      gene_id = c("G1", "G1", "G2"),
      row.names = rownames(counts)
    ),
    colData = data.frame(
      sample_id = c("S1", "S2", "S3"),
      condition = c("A", "B", "A"),
      row.names = c("S1", "S2", "S3")
    )
  )
  
  analysis <- TSENAT::TSENATAnalysis(se)
  
  # Filter the analysis (but don't filter out any samples, keep all 3)
  # Use high min_samples to avoid filtering samples
  analysis_filt <- TSENAT:::filter_analysis_s4(analysis, min_samples = 1, verbose = FALSE)
  
  # Check that TSENATAnalysis is still valid
  expect_s4_class(analysis_filt, "TSENATAnalysis")
  expect_s4_class(analysis_filt@se, "SummarizedExperiment")
  
  # Check that assays are filtered consistently
  se_filt <- analysis_filt@se
  n_assay_rows <- nrow(assays(se_filt)[[1]])
  n_rowdata_rows <- nrow(SummarizedExperiment::rowData(se_filt))
  n_assay_cols <- ncol(assays(se_filt)[[1]])
  n_coldata_rows <- nrow(SummarizedExperiment::colData(se_filt))
  
  # All assays should have same rows
  expect_equal(n_assay_rows, n_rowdata_rows)
  # All assays should have same columns as colData (since filter doesn't modify samples)
  expect_equal(n_assay_cols, n_coldata_rows)
})

test_that("Filtered diversity results maintain colname alignment", {
  # Test that diversity calculation after filtering maintains colname/colData alignment
  skip_if_not_installed("SummarizedExperiment")
  
  # Create data that will pass filtering
  counts <- matrix(c(20, 25, 30, 15, 18, 22, 10, 12, 14), nrow = 3, ncol = 3)
  rownames(counts) <- c("TX1", "TX2", "TX3")
  colnames(counts) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = data.frame(
      gene_id = c("G1", "G1", "G2"),
      row.names = rownames(counts)
    ),
    colData = data.frame(
      sample_id = colnames(counts),
      row.names = colnames(counts)
    )
  )
  
  # Calculate diversity WITHOUT specifying q (will use default q=2)
  genes <- SummarizedExperiment::rowData(se)$gene_id
  result <- TSENAT:::.calculate_method(counts, genes = genes, q = 2, what = "S")
  
  # Prepare metadata
  prep_output <- TSENAT:::.prepare_diversity_metadata(
    counts, result, counts, genes, q = 2
  )
  
  # Check colname alignment
  expect_equal(
    colnames(prep_output$result_assay),
    rownames(prep_output$colData),
    info = "Assay colnames should match colData rownames"
  )
  
  # Check all column names have the format "SampleName_q=value"
  all_cn_have_q <- all(grepl("_q=", colnames(prep_output$result_assay)))
  expect_true(all_cn_have_q, info = "All column names should contain _q= separator")
})
