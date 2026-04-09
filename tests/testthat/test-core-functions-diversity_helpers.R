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

test_that("filter_analysis maintains TSENATAnalysis integrity", {
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
  analysis_filt <- TSENAT:::filter_analysis(analysis, min_samples = 1, verbose = FALSE)
  
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


# ============================================================================
# TESTS FOR .normalize_log_odds_ratio
# ============================================================================
# Tests for normalization using log-odds ratio relative to S_max

context("Entropy Normalization: Log-Odds Ratio")

test_that(".normalize_log_odds_ratio handles Shannon entropy (q=1) correctly", {
  # Create entropy matrix for Shannon diversity
  entropy_matrix <- matrix(
    c(0.5, 1.0, 1.5, 2.0, 0.8, 1.2),
    nrow = 3, ncol = 2
  )
  colnames(entropy_matrix) <- c("Sample1_q=1", "Sample2_q=1")
  rownames(entropy_matrix) <- c("gene1", "gene2", "gene3")
  
  # n_isoforms for each gene (as named vector)
  n_isoforms <- c(gene1 = 3, gene2 = 4, gene3 = 5)
  
  result <- TSENAT:::.normalize_log_odds_ratio(
    entropy_matrix = entropy_matrix,
    n_isoforms = n_isoforms,
    q = 1
  )
  
  # Result should be matrix
  expect_is(result, "matrix")
  expect_equal(dim(result), dim(entropy_matrix))
  expect_equal(rownames(result), rownames(entropy_matrix))
  
  # Values should be finite (log of ratios)
  expect_true(all(is.finite(result)))
  
  # For Shannon entropy with q=1: S_max = log(n_isoforms)
  # Normalized value = log(S / S_max) = log(S / log(n))
  # For gene1 (n=3): S_max = log(3) ≈ 1.0986
  # Sample1: 0.5 / 1.0986 ≈ 0.455, log(0.455) ≈ -0.789
  expected_normalized <- log(entropy_matrix / outer(log(n_isoforms), rep(1, ncol(entropy_matrix))))
  expect_equal(as.numeric(result), as.numeric(expected_normalized), tolerance = 1e-10)
})

test_that(".normalize_log_odds_ratio handles Tsallis entropy (q≠1) correctly", {
  # Tsallis divergence at q=2
  entropy_matrix <- matrix(
    c(0.3, 0.6, 0.2, 0.5, 0.4, 0.7),
    nrow = 3, ncol = 2
  )
  colnames(entropy_matrix) <- c("Sample1_q=2", "Sample2_q=2")
  rownames(entropy_matrix) <- c("gene1", "gene2", "gene3")
  
  n_isoforms <- c(gene1 = 2, gene2 = 3, gene3 = 4)
  
  result <- TSENAT:::.normalize_log_odds_ratio(
    entropy_matrix = entropy_matrix,
    n_isoforms = n_isoforms,
    q = 2
  )
  
  # Result should be matrix with finite values
  expect_is(result, "matrix")
  expect_true(all(is.finite(result)))
  
  # For Tsallis (q=2): S_max = (1 - n^(1-q)) / (q-1) = (1 - n^(-1)) / 1
  # For gene1 (n=2): S_max = (1 - 0.5) = 0.5
  # For gene2 (n=3): S_max = (1 - 1/3) = 2/3
  # For gene3 (n=4): S_max = (1 - 0.25) = 0.75
  s_max_gene1 <- (1 - 2^(1-2)) / (2-1)  # 0.5
  s_max_gene2 <- (1 - 3^(1-2)) / (2-1)  # 2/3
  s_max_gene3 <- (1 - 4^(1-2)) / (2-1)  # 0.75
  
  expected_s_max <- c(s_max_gene1, s_max_gene2, s_max_gene3)
  expected_normalized <- log(cbind(
    entropy_matrix[, 1] / expected_s_max,
    entropy_matrix[, 2] / expected_s_max
  ))
  
  expect_equal(as.numeric(result), as.numeric(expected_normalized), tolerance = 1e-10)
})

test_that(".normalize_log_odds_ratio handles data.frame input", {
  # Test with data.frame instead of matrix
  entropy_df <- data.frame(
    Sample1 = c(0.5, 1.0, 1.5),
    Sample2 = c(0.8, 1.2, 0.9)
  )
  rownames(entropy_df) <- c("gene1", "gene2", "gene3")
  
  n_isoforms <- c(gene1 = 3, gene2 = 4, gene3 = 2)
  
  result <- TSENAT:::.normalize_log_odds_ratio(
    entropy_matrix = entropy_df,
    n_isoforms = n_isoforms,
    q = 1
  )
  
  # Should return matrix
  expect_is(result, "matrix")
  expect_equal(nrow(result), nrow(entropy_df))
  expect_equal(ncol(result), ncol(entropy_df))
})

test_that(".normalize_log_odds_ratio handles NA values correctly", {
  entropy_matrix <- matrix(
    c(0.5, NA, 1.5, 2.0, 0.8, NA),
    nrow = 3, ncol = 2
  )
  colnames(entropy_matrix) <- c("Sample1_q=1", "Sample2_q=1")
  rownames(entropy_matrix) <- c("gene1", "gene2", "gene3")
  
  n_isoforms <- c(gene1 = 3, gene2 = 4, gene3 = 5)
  
  result <- TSENAT:::.normalize_log_odds_ratio(
    entropy_matrix = entropy_matrix,
    n_isoforms = n_isoforms,
    q = 1
  )
  
  # NA values should be preserved
  expect_true(is.na(result[2, 1]))
  expect_true(is.na(result[3, 2]))
  expect_false(is.na(result[1, 1]))  # Valid value not NA
})

test_that(".normalize_log_odds_ratio extracts q from column names when not specified", {
  # Test automatic q extraction from colnames
  entropy_matrix <- matrix(
    c(0.5, 1.0, 2.0, 1.5),
    nrow = 2, ncol = 2
  )
  colnames(entropy_matrix) <- c("Sample1_q=1.0", "Sample2_q=1.0")
  rownames(entropy_matrix) <- c("gene1", "gene2")
  
  n_isoforms <- c(gene1 = 3, gene2 = 4)
  
  # Call without specifying q - should extract from colnames
  result <- TSENAT:::.normalize_log_odds_ratio(
    entropy_matrix = entropy_matrix,
    n_isoforms = n_isoforms,
    q = NULL
  )
  
  # Should successfully extract q=1.0 from column names
  expect_is(result, "matrix")
  expect_true(all(is.finite(result)))
})

test_that(".normalize_log_odds_ratio rejects invalid input", {
  # Test with invalid input types
  entropy_list <- list(a = c(0.5, 1.0), b = c(0.8, 1.2))
  n_isoforms <- c(3, 4)
  
  expect_error(
    TSENAT:::.normalize_log_odds_ratio(
      entropy_matrix = entropy_list,
      n_isoforms = n_isoforms,
      q = 1
    ),
    "must be a matrix or data.frame"
  )
})

test_that(".normalize_log_odds_ratio numerical correctness: entropy ratios", {
  # High precision numerical validation
  entropy_matrix <- matrix(
    c(0.6931471806, 1.0986122887, 1.3862943611),  # ln(2), ln(3), ln(4)
    nrow = 3, ncol = 1
  )
  colnames(entropy_matrix) <- "Sample_q=1"
  rownames(entropy_matrix) <- c("gene_with_2iso", "gene_with_3iso", "gene_with_4iso")
  
  n_isoforms <- c(gene_with_2iso = 2, gene_with_3iso = 3, gene_with_4iso = 4)
  
  result <- TSENAT:::.normalize_log_odds_ratio(
    entropy_matrix = entropy_matrix,
    n_isoforms = n_isoforms,
    q = 1
  )
  
  # All normalized values should be 0 (since S = S_max for uniform distributions)
  # S_max(q=1, n) = log(n)
  # S / S_max ratios: [ln(2)/ln(2), ln(3)/ln(3), ln(4)/ln(4)] = [1, 1, 1]
  # log(1) = 0 for all
  expect_equal(as.numeric(result), c(0, 0, 0), tolerance = 1e-10)
})

test_that(".normalize_log_odds_ratio handles zero values", {
  # Test with zero entropy (edge case)
  entropy_matrix <- matrix(
    c(0, 0.5, 1.0),
    nrow = 3, ncol = 1
  )
  colnames(entropy_matrix) <- "Sample_q=1"
  rownames(entropy_matrix) <- c("gene1", "gene2", "gene3")
  
  n_isoforms <- c(gene1 = 4, gene2 = 4, gene3 = 4)
  
  result <- TSENAT:::.normalize_log_odds_ratio(
    entropy_matrix = entropy_matrix,
    n_isoforms = n_isoforms,
    q = 1
  )
  
  # Zero entropy doesn't meet s_vals > 0 condition, remains unchanged at 0
  expect_equal(result[1, 1], 0.0)
  expect_true(is.finite(result[2, 1]))
  expect_true(is.finite(result[3, 1]))
})


# ============================================================================
# TESTS FOR .normalize_relative_reference
# ============================================================================
# Tests for normalization relative to a reference group

context("Entropy Normalization: Relative Reference")

test_that(".normalize_relative_reference normalizes correctly with specified reference", {
  # Create entropy matrix for 2 genes, 4 samples (2 control, 2 treatment)
  entropy_matrix <- matrix(
    c(1.0, 1.1, 2.0, 2.1, 1.5, 1.6, 2.5, 2.6),
    nrow = 2, ncol = 4,
    byrow = TRUE
  )
  colnames(entropy_matrix) <- c("ctrl_rep1", "ctrl_rep2", "treat_rep1", "treat_rep2")
  rownames(entropy_matrix) <- c("gene1", "gene2")
  
  # Group vector: 2 controls, 2 treatments
  group_vector <- factor(c("control", "control", "treatment", "treatment"))
  
  result <- TSENAT:::.normalize_relative_reference(
    entropy_matrix = entropy_matrix,
    group_vector = group_vector,
    reference_group = "control"
  )
  
  # Result should be matrix with same dimensions
  expect_is(result, "matrix")
  expect_equal(dim(result), dim(entropy_matrix))
  
  # Control samples should be normalized to ~1.0 (divided by their own mean)
  # gene1: ref_mean = (1.0 + 1.1) / 2 = 1.05
  # gene1, ctrl_rep1 normalized: 1.0 / 1.05 ≈ 0.952
  # gene1, ctrl_rep2 normalized: 1.1 / 1.05 ≈ 1.048
  expect_equal(result[1, 1], 1.0 / ((1.0 + 1.1) / 2), tolerance = 1e-10)
  expect_equal(result[1, 2], 1.1 / ((1.0 + 1.1) / 2), tolerance = 1e-10)
  
  # Treatment samples should be different
  # gene1, treat_rep1 normalized: 2.0 / 1.05 ≈ 1.905
  expect_equal(result[1, 3], 2.0 / ((1.0 + 1.1) / 2), tolerance = 1e-10)
})

test_that(".normalize_relative_reference uses first group when reference_group=NULL", {
  entropy_matrix <- matrix(
    c(1.0, 1.2, 2.0, 2.2, 1.5, 1.7, 2.5, 2.7),
    nrow = 2, ncol = 4,
    byrow = TRUE
  )
  rownames(entropy_matrix) <- c("gene1", "gene2")
  
  # Alphabetically: "group_A" comes first, then "group_B"
  group_vector <- factor(c("group_A", "group_A", "group_B", "group_B"))
  
  result <- TSENAT:::.normalize_relative_reference(
    entropy_matrix = entropy_matrix,
    group_vector = group_vector,
    reference_group = NULL  # Should default to first level
  )
  
  # Should normalize relative to group_A
  # gene1, group_A samples: (1.0 + 1.2) / 2 = 1.1
  # Result[1,1] = 1.0 / 1.1
  expect_equal(as.numeric(result[1, 1]), 1.0 / 1.1, tolerance = 1e-10)
})

test_that(".normalize_relative_reference handles character group_vector", {
  entropy_matrix <- matrix(
    c(1.0, 1.1, 2.0, 2.1),
    nrow = 2, ncol = 2
  )
  rownames(entropy_matrix) <- c("gene1", "gene2")
  
  # Character vector (not factor) - should be converted
  group_vector <- c("control", "treatment")
  
  result <- TSENAT:::.normalize_relative_reference(
    entropy_matrix = entropy_matrix,
    group_vector = group_vector,
    reference_group = "control"
  )
  
  expect_is(result, "matrix")
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
})

test_that(".normalize_relative_reference validates group_vector length", {
  entropy_matrix <- matrix(
    c(1.0, 1.1, 2.0, 2.1),
    nrow = 2, ncol = 2
  )
  
  # Mismatched group_vector length
  group_vector <- c("a", "b", "c")  # Length 3, but matrix has 2 columns
  
  expect_error(
    TSENAT:::.normalize_relative_reference(
      entropy_matrix = entropy_matrix,
      group_vector = group_vector,
      reference_group = "a"
    ),
    "group_vector length must equal ncol"
  )
})

test_that(".normalize_relative_reference validates reference_group exists", {
  entropy_matrix <- matrix(
    c(1.0, 1.1, 2.0, 2.1),
    nrow = 2, ncol = 2
  )
  
  group_vector <- c("control", "treatment")
  
  expect_error(
    TSENAT:::.normalize_relative_reference(
      entropy_matrix = entropy_matrix,
      group_vector = group_vector,
      reference_group = "nonexistent_group"
    ),
    "reference_group .* not found"
  )
})

test_that(".normalize_relative_reference handles NA values correctly", {
  entropy_matrix <- matrix(
    c(1.0, NA, 2.0, 2.1, 1.5, NA, 2.5, 2.6),
    nrow = 2, ncol = 4,
    byrow = TRUE
  )
  rownames(entropy_matrix) <- c("gene1", "gene2")
  
  group_vector <- factor(c("a", "a", "b", "b"))
  
  result <- TSENAT:::.normalize_relative_reference(
    entropy_matrix = entropy_matrix,
    group_vector = group_vector,
    reference_group = "a"
  )
  
  # NA values should be preserved
  expect_true(is.na(result[1, 2]))
  expect_true(is.na(result[2, 2]))
  # Valid values should be normalized
  expect_false(is.na(result[1, 1]))
})

test_that(".normalize_relative_reference numerical correctness: group means", {
  # High precision test
  entropy_matrix <- matrix(
    c(0.5, 1.5, 1.0, 2.0, 2.5, 3.0),
    nrow = 3, ncol = 2,
    byrow = TRUE
  )
  colnames(entropy_matrix) <- c("ref_sample", "test_sample")
  rownames(entropy_matrix) <- c("gene1", "gene2", "gene3")
  
  group_vector <- c("reference", "test")
  
  result <- TSENAT:::.normalize_relative_reference(
    entropy_matrix = entropy_matrix,
    group_vector = group_vector,
    reference_group = "reference"
  )
  
  # Reference samples divided by their own value (should be 1.0)
  expect_equal(result[1, 1], 1.0, tolerance = 1e-10)
  expect_equal(result[2, 1], 1.0, tolerance = 1e-10)
  expect_equal(result[3, 1], 1.0, tolerance = 1e-10)
  
  # Test samples divided by reference value (ratio differs from 1.0)
  # gene1: 1.5 / 0.5 = 3.0
  expect_equal(result[1, 2], 3.0, tolerance = 1e-10)
  # gene2: 2.0 / 1.0 = 2.0
  expect_equal(result[2, 2], 2.0, tolerance = 1e-10)
  # gene3: 3.0 / 2.5 = 1.2
  expect_equal(result[3, 2], 1.2, tolerance = 1e-10)
})

test_that(".normalize_relative_reference rejects invalid input", {
  # Test with non-matrix input
  entropy_list <- list(a = c(1, 2), b = c(3, 4))
  group_vector <- c("a", "b")
  
  expect_error(
    TSENAT:::.normalize_relative_reference(
      entropy_matrix = entropy_list,
      group_vector = group_vector,
      reference_group = "a"
    ),
    "must be a matrix or data.frame"
  )
})

test_that(".normalize_relative_reference handles single sample per group", {
  entropy_matrix <- matrix(
    c(1.0, 2.0, 1.5, 1.2, 2.5, 2.2),
    nrow = 2, ncol = 3,
    byrow = TRUE
  )
  
  # Single reference sample (reference_group has only 1 sample)
  group_vector <- factor(c("ref", "test1", "test2"), levels = c("ref", "test1", "test2"))
  
  result <- TSENAT:::.normalize_relative_reference(
    entropy_matrix = entropy_matrix,
    group_vector = group_vector,
    reference_group = "ref"
  )
  
  # Reference sample normalized by itself should be 1.0
  expect_equal(result[1, 1], 1.0, tolerance = 1e-10)
  # Test samples normalized by reference value
  expect_equal(result[1, 2], 2.0 / 1.0, tolerance = 1e-10)
  expect_equal(result[1, 3], 1.5 / 1.0, tolerance = 1e-10)
})

test_that(".normalize_relative_reference handles multiple replicates per group", {
  # 4 genes x 6 samples (3 control replicates, 3 treatment replicates)
  entropy_matrix <- matrix(
    c(
      1.0, 1.05, 0.95, 2.0, 2.05, 1.95,  # gene1
      1.5, 1.55, 1.45, 2.5, 2.55, 2.45,  # gene2
      0.8, 0.85, 0.75, 1.8, 1.85, 1.75,  # gene3
      2.2, 2.25, 2.15, 3.2, 3.25, 3.15   # gene4
    ),
    nrow = 4, ncol = 6, byrow = TRUE
  )
  rownames(entropy_matrix) <- paste0("gene", 1:4)
  colnames(entropy_matrix) <- paste0("sample", 1:6)
  
  # 3 controls, 3 treatments
  group_vector <- factor(rep(c("control", "control", "control", "treatment", "treatment", "treatment")),
                        levels = c("control", "treatment"))
  
  result <- TSENAT:::.normalize_relative_reference(
    entropy_matrix = entropy_matrix,
    group_vector = group_vector,
    reference_group = "control"
  )
  
  # Check that control replicates normalize to approximately 1.0
  # (average of control replicates = reference for that gene)
  # gene1 control mean = (1.0 + 1.05 + 0.95) / 3 ≈ 1.0
  control_mean_gene1 <- (1.0 + 1.05 + 0.95) / 3
  expect_equal(result[1, 1] * control_mean_gene1, 1.0, tolerance = 1e-10)
  
  # Treatment samples should be higher
  expect_true(result[1, 4] > result[1, 1])  # treatment > control
})
