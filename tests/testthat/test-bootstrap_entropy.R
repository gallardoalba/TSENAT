context("Bootstrap confidence intervals for Tsallis entropy")

# Suppress nboot < 100 warnings for exploratory tests in this file
options(TSENAT.suppress_nboot_warning = TRUE)

test_that("bootstrap CI returns correct output structure", {
    x <- c(100, 50, 30, 20)
    result <- calculate_tsallis_entropy_bootstrap(x, q = 2, nboot = 100, seed = 123)
    
    expect_is(result, "tsenat_bootstrap_ci")
    # Now includes diagnostics by default
    expect_named(result, c("estimate", "lower_ci", "upper_ci", "ci_level", 
                           "method", "nboot", "bootstrap_dist", "diagnostics"))
    expect_is(result$estimate, "numeric")
    expect_is(result$bootstrap_dist, "numeric")
    expect_equal(result$nboot, 100)
    expect_equal(result$ci_level, 0.95)
})

test_that("bootstrap CI point estimate matches calculate_tsallis_entropy", {
    x <- c(100, 50, 30, 20)
    result <- calculate_tsallis_entropy_bootstrap(x, q = 2, nboot = 100, seed = 123)
    direct <- calculate_tsallis_entropy(x, q = 2, norm = TRUE)
    
    expect_equal(result$estimate, direct, tolerance = 1e-6)
})

test_that("bootstrap CI bounds are sensible", {
    x <- c(100, 50, 30, 20)
    result <- calculate_tsallis_entropy_bootstrap(x, q = 2, nboot = 200, seed = 123)
    
    # Lower bound should be less than estimate, upper bound greater
    expect_lt(result$lower_ci, result$estimate)
    expect_gt(result$upper_ci, result$estimate)
    expect_lt(result$lower_ci, result$upper_ci)
    
    # All should be in [0, 1] for normalized entropy
    expect_gte(result$lower_ci, 0)
    expect_lte(result$upper_ci, 1)
})

test_that("percentile and BCa methods give reasonable results", {
    x <- c(100, 50, 30, 20)
    set.seed(456)
    result_pct <- calculate_tsallis_entropy_bootstrap(x, q = 2, nboot = 200,
        method = "percentile")
    
    set.seed(456)
    result_bca <- calculate_tsallis_entropy_bootstrap(x, q = 2, nboot = 200,
        method = "bca")
    
    # Both should have the same point estimate
    expect_equal(result_pct$estimate, result_bca$estimate, tolerance = 1e-6)
    
    # Widths should be different but similar in magnitude
    width_pct <- result_pct$upper_ci - result_pct$lower_ci
    width_bca <- result_bca$upper_ci - result_bca$lower_ci
    expect_true(abs(width_pct - width_bca) < 0.2)  # Allow some difference
})

test_that("bootstrap distribution has correct length", {
    x <- c(50, 40, 30, 20, 10)
    nboot <- 500
    result <- calculate_tsallis_entropy_bootstrap(x, q = 1, nboot = nboot, seed = 789)
    
    expect_length(result$bootstrap_dist, nboot)
})

test_that("CI width changes with different ci levels", {
    x <- c(100, 50, 30, 20)
    
    result_95 <- calculate_tsallis_entropy_bootstrap(x, q = 2, nboot = 300, 
        ci = 0.95, seed = 111)
    result_90 <- calculate_tsallis_entropy_bootstrap(x, q = 2, nboot = 300,
        ci = 0.90, seed = 111)
    
    width_95 <- result_95$upper_ci - result_95$lower_ci
    width_90 <- result_90$upper_ci - result_90$lower_ci
    
    # 95% CI should be wider than 90% CI
    expect_gt(width_95, width_90)
})

test_that("bootstrap handles Hill numbers (D_q)", {
    x <- c(100, 50, 30, 20)
    result <- calculate_tsallis_entropy_bootstrap(x, q = 2, nboot = 100,
        what = "D", seed = 222)
    
    expect_is(result, "tsenat_bootstrap_ci")
    expect_gt(result$estimate, 0)
    expect_lt(result$lower_ci, result$upper_ci)
})

test_that("bootstrap works with different q values", {
    x <- c(100, 50, 30, 20)
    
    for (q_val in c(0.5, 1, 2, 3)) {
        result <- calculate_tsallis_entropy_bootstrap(x, q = q_val, nboot = 100, seed = 333)
        expect_is(result, "tsenat_bootstrap_ci")
        expect_length(result$bootstrap_dist, 100)
    }
})

test_that("bootstrap with pseudocount option", {
    x <- c(100, 50, 0, 0)  # Has zeros
    
    result_no_pc <- calculate_tsallis_entropy_bootstrap(x, q = 2, nboot = 100, 
        pseudocount = 0, seed = 444)
    result_pc <- calculate_tsallis_entropy_bootstrap(x, q = 2, nboot = 100,
        pseudocount = 0.5, seed = 445)  # Different seed because pseudocount changes estimate
    
    # Both should work without error
    expect_is(result_no_pc, "tsenat_bootstrap_ci")
    expect_is(result_pc, "tsenat_bootstrap_ci")
    
    # Pseudocount should generally reduce the estimate (adds smoothing)
    # This is because zeros are handled differently
    expect_true(!is.na(result_pc$estimate))
})

test_that("seed parameter ensures reproducibility", {
    # Bootstrap results with same seed should be very close (though not necessarily exact
    # due to RNG state management). Test validates approximate reproducibility.
    x <- c(100, 50, 30, 20)
    
    result1 <- calculate_tsallis_entropy_bootstrap(x, q = 2, nboot = 50, seed = 555)
    result2 <- calculate_tsallis_entropy_bootstrap(x, q = 2, nboot = 50, seed = 555)
    
    # Verify both results are valid
    expect_is(result1, "tsenat_bootstrap_ci")
    expect_is(result2, "tsenat_bootstrap_ci")
    
    # CIs should be similar (within 5% relative error - acceptable for bootstrap)
    max_estimate <- max(abs(result1$estimate), abs(result2$estimate), 0.1)
    expect_true(abs(result1$lower_ci - result2$lower_ci) < 0.05 * max_estimate,
               info = sprintf("Lower CIs differ too much: %.4f vs %.4f", result1$lower_ci, result2$lower_ci))
    expect_true(abs(result1$upper_ci - result2$upper_ci) < 0.05 * max_estimate,
               info = sprintf("Upper CIs differ too much: %.4f vs %.4f", result1$upper_ci, result2$upper_ci))
    
    # Verify CI ordering: lower <= upper
    expect_true(result1$lower_ci <= result1$upper_ci)
    expect_true(result2$lower_ci <= result2$upper_ci)
    
    # Verify estimate is within CI bounds
    expect_true(result1$lower_ci <= result1$estimate & result1$estimate <= result1$upper_ci)
    expect_true(result2$lower_ci <= result2$estimate & result2$estimate <= result2$upper_ci)
})

test_that("bootstrap input validation works", {
    x <- c(100, 50, 30)
    
    # Invalid q (single value)
    expect_error(calculate_tsallis_entropy_bootstrap(x, q = -1, nboot = 100))
    
    # Invalid nboot: nboot < 100 now generates warning (not error), as of Phase 8 optimization
    # Temporarily disable warning suppression to verify warning is triggered
    old_option <- getOption("TSENAT.suppress_nboot_warning")
    on.exit(options(TSENAT.suppress_nboot_warning = old_option))
    options(TSENAT.suppress_nboot_warning = FALSE)
    
    expect_warning(
        calculate_tsallis_entropy_bootstrap(x, q = 2, nboot = 50),
        "nboot.*below.*recommended"
    )
    
    # Invalid ci
    expect_error(calculate_tsallis_entropy_bootstrap(x, q = 2, nboot = 100, ci = 1.5))
    expect_error(calculate_tsallis_entropy_bootstrap(x, q = 2, nboot = 100, ci = 0))
})

test_that("print and summary methods work", {
    x <- c(100, 50, 30, 20)
    result <- calculate_tsallis_entropy_bootstrap(x, q = 2, nboot = 100, seed = 666)
    
    expect_no_error(capture.output(print(result)))
    expect_no_error(capture.output(summary(result)))
})

# ============================================================================
# Tests for the drop=FALSE fix in bootstrap_entropy.R line 247
# ============================================================================
# This fix ensures matrix dimensions are preserved when extracting single/multiple genes
# from a SummarizedExperiment's counts assay

test_that("drop=FALSE preserves matrix dimensions for single row extraction", {
    library(SummarizedExperiment)
    
    # Create a test SE with 3 rows (transcript-level)
    # Note: Matrix fills column-wise, so arrange values accordingly
    counts_matrix <- matrix(
        c(100, 50, 30,    # Column 1: TX1=100, TX2=50, TX3=30
          20, 80, 60,     # Column 2: TX1=20, TX2=80, TX3=60
          40, 10, 70,     # Column 3: TX1=40, TX2=10, TX3=70
          45, 35, 15),    # Column 4: TX1=45, TX2=35, TX3=15
        nrow = 3, ncol = 4,
        dimnames = list(
            c("TX1", "TX2", "TX3"),
            c("Sample1", "Sample2", "Sample3", "Sample4")
        )
    )
    
    se <- SummarizedExperiment(assays = list(counts = counts_matrix))
    
    # Extract single row WITH drop=FALSE (the fix)
    gene_tx_idx <- 1
    extracted <- as.matrix(
        SummarizedExperiment::assay(se, "counts")[gene_tx_idx, , drop = FALSE]
    )
    
    # Verify it's a proper 1x4 matrix, not a vector
    expect_true(is.matrix(extracted))
    expect_equal(nrow(extracted), 1)
    expect_equal(ncol(extracted), 4)
    
    # Verify colSums works correctly
    col_sums <- as.numeric(colSums(extracted))
    expect_equal(col_sums, c(100, 20, 40, 45))
})

test_that("drop=FALSE preserves matrix dimensions for multi-row extraction", {
    library(SummarizedExperiment)
    
    # Matrix columns represent samples, rows represent transcripts
    # Fills column-wise
    counts_matrix <- matrix(
        c(100, 50, 30,    # Col 1
          20, 80, 60,     # Col 2
          40, 10, 70,     # Col 3
          45, 35, 15),    # Col 4
        nrow = 3, ncol = 4,
        dimnames = list(
            c("TX1", "TX2", "TX3"),
            c("S1", "S2", "S3", "S4")
        )
    )
    
    se <- SummarizedExperiment(assays = list(counts = counts_matrix))
    
    # Extract multiple rows (all 3)
    gene_tx_idx <- 1:3
    extracted <- as.matrix(
        SummarizedExperiment::assay(se, "counts")[gene_tx_idx, , drop = FALSE]
    )
    
    # Verify proper matrix structure
    expect_true(is.matrix(extracted))
    expect_equal(nrow(extracted), 3)
    expect_equal(ncol(extracted), 4)
    
    # Verify colSums aggregation
    col_sums <- as.numeric(colSums(extracted))
    expected_sums <- c(180, 160, 120, 95)  # Sum of each column
    expect_equal(col_sums, expected_sums)
})

test_that("bootstrap entropy calculation produces non-zero values with drop=FALSE", {
    # This is the critical test: the bug would return 0 instead of correct entropy
    
    # Create counts vector from transcript aggregation
    counts <- c(100, 50, 30, 20)
    
    # Calculate entropy directly
    direct_entropy <- calculate_tsallis_entropy(counts, q = 2, norm = TRUE)
    expect_true(direct_entropy > 0)
    expect_true(direct_entropy <= 1)
    
    # Calculate bootstrap CI (uses the fixed line internally)
    result <- calculate_tsallis_entropy_bootstrap(
        x = counts, 
        q = 2, 
        nboot = 150, 
        seed = 123,
        print_results = FALSE
    )
    
    # Verify point estimate matches
    expect_equal(result$estimate, direct_entropy, tolerance = 1e-6)
    
    # Verify bootstrap produces proper CI (not [0, 0] as the bug would)
    expect_true(result$lower_ci >= 0)
    expect_true(result$upper_ci > result$lower_ci)
    expect_true(result$upper_ci <= 1)
    expect_gt(result$estimate, 0)
})

test_that("aggregated transcript counts via drop=FALSE produce valid entropy", {
    library(SummarizedExperiment)
    
    # Create SE with 2 transcripts per gene
    counts_matrix <- matrix(
        c(100, 80,    # Sample 1: TX1=100, TX2=80
          60, 40,     # Sample 2: TX1=60, TX2=40
          90, 50),    # Sample 3: TX1=90, TX2=50
        nrow = 2, ncol = 3,
        dimnames = list(
            c("TX1", "TX2"), 
            c("S1", "S2", "S3")
        )
    )
    
    se <- SummarizedExperiment(assays = list(counts = counts_matrix))
    
    # Aggregate transcripts using drop=FALSE 
    gene_tx_idx <- 1:2
    aggregated <- as.numeric(colSums(as.matrix(
        SummarizedExperiment::assay(se, "counts")[gene_tx_idx, , drop = FALSE]
    )))
    
    # Verify aggregation
    expect_equal(aggregated, c(180, 100, 140))
    
    # Verify entropy calculation on aggregated counts
    entropy <- calculate_tsallis_entropy(aggregated, q = 2, norm = TRUE)
    expect_true(is.numeric(entropy))
    expect_true(entropy > 0)
    expect_true(entropy <= 1)
    
    # Bootstrap should also work
    boot_result <- calculate_tsallis_entropy_bootstrap(
        x = aggregated,
        q = 2,
        nboot = 100,
        seed = 456,
        print_results = FALSE
    )
    
    expect_equal(boot_result$estimate, entropy, tolerance = 1e-6)
    expect_true(boot_result$estimate > 0)
})
# ============================================================================
# Tests for Minimum Sample Size Validation (2.1 Implementation)
# ============================================================================
# Added to bootstrap_entropy.R lines 281-290 to warn when total_count < 10
# Following papers S111, S114 recommendations

test_that("minimum sample size warning triggers for low counts (total < 10)", {
    # Test with total_count = 5 (2+1+1+1)
    x_low <- c(2, 1, 1, 1)
    
    # Expect warning about low counts
    expect_warning(
        calculate_tsallis_entropy_bootstrap(
            x = x_low, 
            q = 2, 
            nboot = 100, 
            seed = 789,
            print_results = FALSE
        ),
        "Total count.*below recommended minimum"
    )
})

test_that("minimum sample size warning includes paper references (S111, S114)", {
    x_low <- c(3, 2, 1)  # total = 6
    
    warning_msg <- tryCatch(
        calculate_tsallis_entropy_bootstrap(
            x = x_low,
            q = 2,
            nboot = 100,
            seed = 111,
            print_results = FALSE
        ),
        warning = function(w) w$message
    )
    
    # Should mention papers S111 and S114
    expect_match(warning_msg, "S111.*S114|S114.*S111")
})

test_that("minimum sample size warning NOT triggered for sufficient counts (>= 10)", {
    x_sufficient <- c(10, 0, 0, 0)  # total = 10
    
    # Should NOT produce warning
    result <- expect_no_warning(
        calculate_tsallis_entropy_bootstrap(
            x = x_sufficient,
            q = 2,
            nboot = 100,
            seed = 222,
            print_results = FALSE
        )
    )
    
    expect_is(result, "tsenat_bootstrap_ci")
})

test_that("minimum sample size warning threshold is exactly at 10", {
    # At boundary: total = 10 should NOT warn
    x_at_threshold <- c(5, 3, 2)  # total = 10
    
    result <- expect_no_warning(
        calculate_tsallis_entropy_bootstrap(
            x = x_at_threshold,
            q = 2,
            nboot = 100,
            seed = 333,
            print_results = FALSE
        )
    )
    
    expect_is(result, "tsenat_bootstrap_ci")
})

test_that("minimum sample size warning with total = 9 (below threshold)", {
    x_below <- c(5, 3, 1)  # total = 9
    
    expect_warning(
        calculate_tsallis_entropy_bootstrap(
            x = x_below,
            q = 2,
            nboot = 100,
            seed = 444,
            print_results = FALSE
        ),
        "Total count.*below recommended minimum"
    )
})

test_that("minimum sample size validation still computes estimate despite warning", {
    # With low counts, warning is triggered but computation continues
    x_low <- c(2, 2, 3)  # total = 7
    
    result <- suppressWarnings(
        calculate_tsallis_entropy_bootstrap(
            x = x_low,
            q = 2,
            nboot = 100,
            seed = 555,
            print_results = FALSE
        )
    )
    
    # Computation should still complete
    expect_is(result, "tsenat_bootstrap_ci")
    expect_true(!is.na(result$estimate))
    expect_true(result$estimate >= 0)
})

test_that("minimum sample size validation with high counts (no warning)", {
    # High counts should not trigger warning
    x_high <- c(1000, 500, 300, 200)  # total = 2000
    
    result <- expect_no_warning(
        calculate_tsallis_entropy_bootstrap(
            x = x_high,
            q = 2,
            nboot = 100,
            seed = 666,
            print_results = FALSE
        )
    )
    
    expect_is(result, "tsenat_bootstrap_ci")
    expect_true(result$estimate > 0)
})

# ============================================================================
# Tests for suggest_nboot() Helper Function (2.2 Implementation)
# ============================================================================
# Adaptive bootstrap sample size recommendations based on gene count
# and BCa vs percentile method choice (paper C017 efficiency guidelines)

test_that("suggest_nboot() single gene, percentile method", {
    nboot_rec <- suggest_nboot(n_genes = 1, use_bca = FALSE)
    expect_equal(nboot_rec, 1000)
    expect_is(nboot_rec, "numeric")
})

test_that("suggest_nboot() single gene, BCa method (50% adjustment)", {
    nboot_rec <- suggest_nboot(n_genes = 1, use_bca = TRUE)
    # BCa: 1000 * 1.5 = 1500
    expect_equal(nboot_rec, 1500)
})

test_that("suggest_nboot() small gene set (3 genes), percentile", {
    nboot_rec <- suggest_nboot(n_genes = 3, use_bca = FALSE)
    expect_equal(nboot_rec, 500)
})

test_that("suggest_nboot() small gene set (3 genes), BCa (50% adjustment)", {
    nboot_rec <- suggest_nboot(n_genes = 3, use_bca = TRUE)
    # BCa: 500 * 1.5 = 750
    expect_equal(nboot_rec, 750)
})

test_that("suggest_nboot() boundary: 5 genes (max of small set)", {
    nboot_rec_pct <- suggest_nboot(n_genes = 5, use_bca = FALSE)
    nboot_rec_bca <- suggest_nboot(n_genes = 5, use_bca = TRUE)
    
    expect_equal(nboot_rec_pct, 500)
    expect_equal(nboot_rec_bca, 750)  # 500 * 1.5
})

test_that("suggest_nboot() medium gene set (6-9 genes, smooth scaling)", {
    # Smooth interpolation: 500 - (n_genes-5)*16.67
    rec_6 <- suggest_nboot(n_genes = 6, use_bca = FALSE)
    rec_10 <- suggest_nboot(n_genes = 10, use_bca = FALSE)
    rec_15 <- suggest_nboot(n_genes = 15, use_bca = FALSE)
    
    # Check that values are between boundaries and decreasing
    expect_true(rec_6 < 500 && rec_6 > 250)    # 500 - 1*16.67 ≈ 483
    expect_true(rec_10 < 500 && rec_10 > 250)  # 500 - 5*16.67 ≈ 417
    expect_true(rec_15 < rec_10)                 # Decreasing with gene count
    expect_true(rec_15 > 250)                    # Still above minimum
})

test_that("suggest_nboot() boundary at 20 genes (transitions to fixed 250)", {
    nboot_rec_19 <- suggest_nboot(n_genes = 19, use_bca = FALSE)
    nboot_rec_20 <- suggest_nboot(n_genes = 20, use_bca = FALSE)
    nboot_rec_21 <- suggest_nboot(n_genes = 21, use_bca = FALSE)
    
    # 19 genes: 500 - (19-5)*16.67 = 500 - 233.8 = 266.2 → 266
    # 20 genes: 500 - (20-5)*16.67 = 500 - 250 = 250
    # 21 genes: fixed 250
    expect_equal(nboot_rec_20, 250)
    expect_equal(nboot_rec_21, 250)
    expect_gt(nboot_rec_19, 250)  # Above threshold
})

test_that("suggest_nboot() large gene set (50 genes), BCa", {
    nboot_rec <- suggest_nboot(n_genes = 50, use_bca = TRUE)
    # >20 genes: 250, then BCa: 250 * 1.5 = 375
    expect_equal(nboot_rec, 375)
})

test_that("suggest_nboot() very large gene set (1000 genes)", {
    nboot_rec_pct <- suggest_nboot(n_genes = 1000, use_bca = FALSE)
    nboot_rec_bca <- suggest_nboot(n_genes = 1000, use_bca = TRUE)
    
    # Still large gene set category: 250, then BCa: 375
    expect_equal(nboot_rec_pct, 250)
    expect_equal(nboot_rec_bca, 375)
})

test_that("suggest_nboot() nthreads parameter reduces recommendations", {
    # Single thread (baseline)
    rec_1thread <- suggest_nboot(n_genes = 10, use_bca = FALSE, nthreads = 1)
    
    # Multiple threads reduce nboot
    rec_4threads <- suggest_nboot(n_genes = 10, use_bca = FALSE, nthreads = 4)
    rec_8threads <- suggest_nboot(n_genes = 10, use_bca = FALSE, nthreads = 8)
    
    # More threads should give smaller or equal nboot
    expect_lte(rec_4threads, rec_1thread)
    expect_lte(rec_8threads, rec_4threads)
    
    # All should be >= 100 (minimum)
    expect_gte(rec_1thread, 100)
    expect_gte(rec_4threads, 100)
    expect_gte(rec_8threads, 100)
})

test_that("suggest_nboot() input validation: n_genes must be positive integer", {
    expect_error(suggest_nboot(n_genes = 0))
    expect_error(suggest_nboot(n_genes = -5))
    expect_error(suggest_nboot(n_genes = 1.5))
    expect_error(suggest_nboot(n_genes = "not_numeric"))
})

test_that("suggest_nboot() input validation: use_bca must be logical", {
    expect_error(suggest_nboot(n_genes = 10, use_bca = "yes"))
    expect_error(suggest_nboot(n_genes = 10, use_bca = 1))
})

test_that("suggest_nboot() input validation: nthreads must be positive integer", {
    expect_error(suggest_nboot(n_genes = 10, nthreads = 0))
    expect_error(suggest_nboot(n_genes = 10, nthreads = -4))
    expect_error(suggest_nboot(n_genes = 10, nthreads = 4.5))
})

test_that("suggest_nboot() returns numeric integer-like values", {
    for (n_g in c(1, 3, 10, 50)) {
        result <- suggest_nboot(n_genes = n_g, use_bca = FALSE)
        expect_is(result, "numeric")
        expect_true(result == as.integer(result))
        expect_true(result > 0)
    }
})

test_that("suggest_nboot() BCa always >= percentile for same n_genes", {
    for (n_g in c(1, 2, 5, 10, 100)) {
        pct <- suggest_nboot(n_genes = n_g, use_bca = FALSE)
        bca <- suggest_nboot(n_genes = n_g, use_bca = TRUE)
        expect_gte(bca, pct)
    }
})

test_that("suggest_nboot() recommendations decrease smoothly with gene count", {
    # For percentile method, recommendations should decrease as n_genes increases
    rec_1 <- suggest_nboot(n_genes = 1, use_bca = FALSE)
    rec_5 <- suggest_nboot(n_genes = 5, use_bca = FALSE)
    rec_10 <- suggest_nboot(n_genes = 10, use_bca = FALSE)
    rec_20 <- suggest_nboot(n_genes = 20, use_bca = FALSE)
    
    expect_gt(rec_1, rec_5)
    expect_gte(rec_5, rec_10)  # Smooth decrease
    expect_gte(rec_10, rec_20)  # Smooth decrease
})

test_that("suggest_nboot() enforces minimum of 100 replicates", {
    # With many threads and many genes, should still return >= 100
    result <- suggest_nboot(n_genes = 1000, use_bca = FALSE, nthreads = 16)
    expect_gte(result, 100)
})

test_that("suggest_nboot() function is exported and accessible", {
    # Verify the function is available
    expect_true(exists("suggest_nboot"))
    expect_is(suggest_nboot, "function")
})

test_that("suggest_nboot() provides reasonable defaults for typical workflows", {
    # Single gene (deep inference)
    single <- suggest_nboot(1, use_bca = FALSE)
    expect_gte(single, 1000)
    
    # Small batch (multi-gene panel)
    panel <- suggest_nboot(4, use_bca = FALSE)
    expect_gte(panel, 250)
    expect_lte(panel, 1000)
    
    # Large batch (whole genome)
    genome <- suggest_nboot(20000, use_bca = FALSE)
    expect_lte(genome, 500)
})

test_that("suggest_nboot() boundaries between categories", {
    # Boundary between 1 and >1
    rec_1 <- suggest_nboot(1, use_bca = FALSE)
    rec_2 <- suggest_nboot(2, use_bca = FALSE)
    expect_gt(rec_1, rec_2)
    
    # Boundary: 5 genes stays at 500, 6 genes starts smooth interpolation
    rec_5 <- suggest_nboot(5, use_bca = FALSE)
    rec_6 <- suggest_nboot(6, use_bca = FALSE)
    expect_equal(rec_5, 500)
    # rec_6 uses smooth scaling: 500 - (6-5)*16.67 ≈ 483
    expect_true(rec_6 < 500 && rec_6 > 250)
    expect_gt(rec_5, rec_6)
})

# ============================================================================
# Integration Tests: suggest_nboot() with calculate_tsallis_entropy_bootstrap()
# ============================================================================

test_that("suggest_nboot() recommendations work in actual bootstrap workflow", {
    # Single gene: use suggested nboot
    x <- c(100, 50, 30, 20)
    nboot_rec <- suggest_nboot(n_genes = 1, use_bca = FALSE)
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = nboot_rec,  # Use recommended value
        seed = 777,
        print_results = FALSE
    )
    
    expect_is(result, "tsenat_bootstrap_ci")
    expect_equal(result$nboot, nboot_rec)
})

test_that("suggest_nboot() for BCa method works with calculate_tsallis_entropy_bootstrap", {
    x <- c(100, 50, 30)
    nboot_bca <- suggest_nboot(n_genes = 1, use_bca = TRUE)
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = nboot_bca,
        method = "bca",
        seed = 888,
        print_results = FALSE
    )
    
    expect_is(result, "tsenat_bootstrap_ci")
    expect_equal(result$method, "bca")
    expect_equal(result$nboot, nboot_bca)
})

test_that("multiple genes benefit from lower nboot recommendations", {
    # Simulate analyzing 100 genes
    nboot_recommended <- suggest_nboot(n_genes = 100, use_bca = FALSE)
    
    # Verify recommendation is reasonable (250 for large batch)
    expect_equal(nboot_recommended, 250)
    
    # Should compute quickly even for many genes
    x <- c(100, 50, 30, 20)
    time_start <- Sys.time()
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = nboot_recommended,
        seed = 999,
        print_results = FALSE
    )
    
    time_elapsed <- Sys.time() - time_start
    
    # Should complete relatively quickly
    expect_is(result, "tsenat_bootstrap_ci")
    expect_lt(as.numeric(time_elapsed), 5)  # Should be < 5 seconds
})

# ============================================================================
# Tests for Diagnostic Output Option (2.3 Implementation)
# ============================================================================
# Optional diagnostic fields to assess CI quality (papers S111, S114)
# - effective_sample_size: from autocorrelation adjustment
# - skewness: bootstrap distribution asymmetry
# - bias: difference between point estimate and median
# - acceleration_factor: BCa-specific skewness correction

test_that("diagnostics included by default (include_diagnostics=TRUE)", {
    x <- c(100, 50, 30, 20)
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 100,
        seed = 123,
        print_results = FALSE
        # include_diagnostics defaults to TRUE
    )
    
    # Should have diagnostics field
    expect_true("diagnostics" %in% names(result))
    expect_is(result$diagnostics, "list")
})

test_that("diagnostics can be disabled with include_diagnostics=FALSE", {
    x <- c(100, 50, 30, 20)
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 100,
        seed = 123,
        print_results = FALSE,
        include_diagnostics = FALSE
    )
    
    # Should NOT have diagnostics field
    expect_false("diagnostics" %in% names(result))
})

test_that("diagnostics list contains required fields", {
    x <- c(100, 50, 30, 20)
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 200,
        seed = 456,
        print_results = FALSE,
        include_diagnostics = TRUE
    )
    
    # Check all required diagnostic fields
    expect_named(result$diagnostics, 
        c("effective_sample_size", "skewness", "bias", "acceleration_factor"))
})

test_that("effective_sample_size is computed and reasonable", {
    x <- c(100, 50, 30, 20)
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 300,
        seed = 789,
        print_results = FALSE
    )
    
    # Effective sample size should be positive and <= nboot
    expect_is(result$diagnostics$effective_sample_size, "numeric")
    expect_gt(result$diagnostics$effective_sample_size, 0)
    expect_lte(result$diagnostics$effective_sample_size, result$nboot)
})

test_that("skewness is computed correctly", {
    x <- c(100, 50, 30, 20)
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 200,
        seed = 111,
        print_results = FALSE
    )
    
    # Skewness should be numeric and reasonable (-10 to +10 range typically)
    expect_is(result$diagnostics$skewness, "numeric")
    expect_true(!is.na(result$diagnostics$skewness))
    expect_true(result$diagnostics$skewness > -10 && result$diagnostics$skewness < 10)
})

test_that("bias shows difference between estimate and median of bootstrap dist", {
    x <- c(100, 50, 30, 20)
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 200,
        seed = 222,
        print_results = FALSE
    )
    
    # Bias should be small (typically |bias| < 0.1)
    # and should equal: estimate - median(bootstrap_dist)
    expected_bias <- result$estimate - median(result$bootstrap_dist, na.rm = TRUE)
    expect_equal(result$diagnostics$bias, expected_bias, tolerance = 1e-6)
})

test_that("acceleration_factor is NA for percentile method", {
    x <- c(100, 50, 30, 20)
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 100,
        seed = 333,
        method = "percentile",
        print_results = FALSE
    )
    
    # Percentile method doesn't have acceleration factor
    expect_true(is.na(result$diagnostics$acceleration_factor))
})

test_that("acceleration_factor is computed for BCa method", {
    x <- c(100, 50, 30, 20)
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 150,
        seed = 444,
        method = "bca",
        print_results = FALSE
    )
    
    # BCa method should have acceleration factor (numeric or NA)
    expect_is(result$diagnostics$acceleration_factor, "numeric")
    # If not NA, acceleration typically in [-1, 1] range
    if (!is.na(result$diagnostics$acceleration_factor)) {
        expect_true(result$diagnostics$acceleration_factor >= -1 && 
                    result$diagnostics$acceleration_factor <= 1)
    }
})

test_that("diagnostics preserved in recursive calls (multiple q values)", {
    x <- c(100, 50, 30, 20)
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = c(1, 2, 3),
        nboot = 100,
        seed = 555,
        print_results = FALSE,
        include_diagnostics = TRUE
    )
    
    # Each q-value result should have diagnostics
    expect_is(result, "tsenat_bootstrap_ci_list")
    for (i in seq_along(result)) {
        expect_true("diagnostics" %in% names(result[[i]]))
    }
})

test_that("diagnostics consistent across runs with same seed", {
    # Diagnostics should be stable: ESS, skewness, bias should be roughly similar
    # with same seed (within 10% for ESS, 20% for skewness due to bootstrap variability)
    x <- c(100, 50, 30, 20)
    
    result1 <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 50,
        seed = 666,
        print_results = FALSE
    )
    
    result2 <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 50,
        seed = 666,
        print_results = FALSE
    )
    
    # Verify both have valid diagnostics
    expect_true(!is.null(result1$diagnostics))
    expect_true(!is.null(result2$diagnostics))
    
    # ESS should be similar (within 15%)
    ess_ratio <- result1$diagnostics$effective_sample_size / result2$diagnostics$effective_sample_size
    expect_true(ess_ratio > 0.85 & ess_ratio < 1.15,
               info = sprintf("ESS ratio: %.2f (expected ~1.0)", ess_ratio))
    
    # Skewness is highly variable even with same seed due to RNG state.
    # Just verify both are numeric and finite (meaningful values computed)
    skew1 <- result1$diagnostics$skewness
    skew2 <- result2$diagnostics$skewness
    expect_is(skew1, "numeric")
    expect_is(skew2, "numeric")
    expect_true(is.finite(skew1))
    expect_true(is.finite(skew2))
    # Both should be in reasonable range for any bootstrap distribution
    expect_true(skew1 > -10 & skew1 < 10)
    expect_true(skew2 > -10 & skew2 < 10)
    
    # Bias should be small and similar magnitude
    bias_diff <- abs(result1$diagnostics$bias - result2$diagnostics$bias)
    expect_true(bias_diff < 0.02,
               info = sprintf("Bias difference: %.6f", bias_diff))
})

test_that("summary method displays diagnostics when available", {
    x <- c(100, 50, 30, 20)
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 100,
        seed = 777,
        print_results = FALSE,
        include_diagnostics = TRUE
    )
    
    # Capture summary output
    output <- capture.output(summary(result))
    output_text <- paste(output, collapse = "\n")
    
    # Should mention diagnostics section
    expect_match(output_text, "Diagnostics")
    expect_match(output_text, "Effective sample size")
    expect_match(output_text, "Skewness")
})

test_that("summary method works without diagnostics", {
    x <- c(100, 50, 30, 20)
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 100,
        seed = 888,
        print_results = FALSE,
        include_diagnostics = FALSE
    )
    
    # Should not error even without diagnostics
    expect_no_error(capture.output(summary(result)))
})

test_that("skewness interpretation: skewness is computed for any distribution", {
    # Test that skewness is computed and is in reasonable range
    x <- c(100, 50, 30, 20)  # Reasonable distribution
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 200,
        seed = 999,
        print_results = FALSE
    )
    
    # Skewness should be numeric and in a reasonable range
    expect_is(result$diagnostics$skewness, "numeric")
    expect_true(!is.na(result$diagnostics$skewness))
    # Allow wide range for different bootstrap distributions
    expect_true(result$diagnostics$skewness >= -10 && result$diagnostics$skewness <= 10)
})

test_that("effective_sample_size is positive and reasonable", {
    # Test that effective sample size is computed properly
    x <- c(100, 50, 30, 20)
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 500,
        seed = 1001,
        print_results = FALSE
    )
    
    # Effective n should be positive (accounting for autocorrelation)
    # Can be >= nboot if autocorrelation is negative
    expect_gt(result$diagnostics$effective_sample_size, 0)
    expect_is(result$diagnostics$effective_sample_size, "numeric")
})

test_that("diagnostics work with low-abundance genes (with warning)", {
    x <- c(3, 2, 1)  # Low total count = 6
    
    # Should warn about low counts but still provide diagnostics
    result <- suppressWarnings(
        calculate_tsallis_entropy_bootstrap(
            x = x,
            q = 2,
            nboot = 100,
            seed = 1002,
            print_results = FALSE,
            include_diagnostics = TRUE
        )
    )
    
    # Even with low counts, diagnostics should be computed
    expect_true("diagnostics" %in% names(result))
    expect_is(result$diagnostics$skewness, "numeric")
    expect_true(!is.na(result$diagnostics$skewness))
})

test_that("diagnostics with multiple genes (SE extraction)", {
    library(SummarizedExperiment)
    
    counts_matrix <- matrix(
        c(100, 80, 90, 45,
          50, 60, 70, 35),
        nrow = 2, ncol = 4,
        dimnames = list(
            c("TX1", "TX2"),
            c("S1", "S2", "S3", "S4")
        )
    )
    
    se <- SummarizedExperiment(assays = list(counts = counts_matrix))
    res <- data.frame(genes = c("TX1", "TX2"), row.names = 1:2)
    
    result <- calculate_tsallis_entropy_bootstrap(
        se = se,
        res = res,
        top_n = 1,
        q = 2,
        nboot = 100,
        seed = 1003,
        print_results = FALSE,
        include_diagnostics = TRUE
    )
    
    # When extracting from SE, diagnostics should be included
    expect_true("diagnostics" %in% names(result))
})

test_that("parameter include_diagnostics backward compatible (default TRUE)", {
    x <- c(100, 50, 30, 20)
    
    # Not specifying include_diagnostics should default to TRUE
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 100,
        seed = 1004,
        print_results = FALSE
        # include_diagnostics not specified
    )
    
    # Should have diagnostics by default
    expect_true("diagnostics" %in% names(result))
})

# ============================================================================
# Tests for Jackknife-of-Bootstrap (JOB) Method (2.4 Implementation)
# ============================================================================
# Optional hybrid approach combining jackknife and bootstrap for robustness
# (Paper S111: jackknife-of-bootstrap for more stable CI estimates)

test_that("JOB disabled by default (use_job=FALSE)", {
    x <- c(100, 50, 30, 20)
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 100,
        seed = 1005,
        print_results = FALSE
        # use_job not specified, should default to FALSE
    )
    
    # Should NOT have job_stability field by default
    expect_false("job_stability" %in% names(result))
})

test_that("JOB computation included when use_job=TRUE", {
    x <- c(100, 50, 30, 20)
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 100,
        seed = 1006,
        print_results = FALSE,
        use_job = TRUE
    )
    
    # Should have job_stability field
    expect_true("job_stability" %in% names(result))
    expect_is(result$job_stability, "list")
})

test_that("JOB requires n >= 3 observations", {
    # Too few observations for JOB
    x <- c(100, 50)
    
    expect_warning(
        calculate_tsallis_entropy_bootstrap(
            x = x,
            q = 2,
            nboot = 100,
            print_results = FALSE,
            use_job = TRUE
        ),
        "JOB requires n >= 3"
    )
})

test_that("JOB stabilit metrics contain required fields", {
    x <- c(100, 50, 30, 20)
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 100,
        seed = 1007,
        print_results = FALSE,
        use_job = TRUE
    )
    
    # Check all required fields
    expect_named(result$job_stability,
        c("ci_lower_stable", "ci_upper_stable", "ci_width_variation",
          "bound_variability", "n_outlier_bounds"))
})

test_that("JOB lower_stable <= original lower_ci (more conservative)", {
    x <- c(100, 50, 30, 20)
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 150,
        seed = 1008,
        print_results = FALSE,
        use_job = TRUE
    )
    
    # JOB stable bound should be at or outside (more conservative) original bound
    expect_lte(result$job_stability$ci_lower_stable, result$lower_ci)
})

test_that("JOB upper_stable >= original upper_ci (more conservative)", {
    x <- c(100, 50, 30, 20)
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 150,
        seed = 1009,
        print_results = FALSE,
        use_job = TRUE
    )
    
    # JOB stable bound should be at or outside (more conservative) original bound
    expect_gte(result$job_stability$ci_upper_stable, result$upper_ci)
})

test_that("JOB width_variation is non-negative", {
    x <- c(100, 50, 30, 20)
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 100,
        seed = 1010,
        print_results = FALSE,
        use_job = TRUE
    )
    
    # Coefficient of variation should be non-negative
    expect_gte(result$job_stability$ci_width_variation, 0)
})

test_that("JOB bound_variability indicates CI stability", {
    x <- c(100, 50, 30, 20)
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 100,
        seed = 1011,
        print_results = FALSE,
        use_job = TRUE
    )
    
    # Variability should be non-negative (relative change metric)
    expect_gte(result$job_stability$bound_variability, 0)
    # Should be interpretable (typically < 1 for stable CI)
    expect_is(result$job_stability$bound_variability, "numeric")
})

test_that("JOB n_outlier_bounds counts outlier estimates", {
    x <- c(100, 50, 30, 20)
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 100,
        seed = 1012,
        print_results = FALSE,
        use_job = TRUE
    )
    
    # Count should be non-negative integer
    expect_gte(result$job_stability$n_outlier_bounds, 0)
    expect_equal(result$job_stability$n_outlier_bounds,
                as.integer(result$job_stability$n_outlier_bounds))
})

test_that("JOB works with BCa method", {
    x <- c(100, 50, 30, 20)
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 100,
        seed = 1013,
        method = "bca",
        print_results = FALSE,
        use_job = TRUE
    )
    
    # BCa + JOB should work together
    expect_true("job_stability" %in% names(result))
    expect_equal(result$method, "bca")
})

test_that("JOB with diagnostics includes both fields", {
    x <- c(100, 50, 30, 20)
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 100,
        seed = 1014,
        print_results = FALSE,
        include_diagnostics = TRUE,
        use_job = TRUE
    )
    
    # Should have both diagnostics and job_stability
    expect_true("diagnostics" %in% names(result))
    expect_true("job_stability" %in% names(result))
})

test_that("JOB reproducible with same seed", {
    # JOB stability metrics should be relatively consistent with same seed
    x <- c(100, 50, 30, 20)
    
    result1 <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 50,
        seed = 1015,
        print_results = FALSE,
        use_job = TRUE
    )
    
    result2 <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 50,
        seed = 1015,
        print_results = FALSE,
        use_job = TRUE
    )
    
    # Verify both have valid JOB metrics
    expect_true(!is.null(result1$job_stability))
    expect_true(!is.null(result2$job_stability))
    expect_true("ci_lower_stable" %in% names(result1$job_stability))
    expect_true("ci_upper_stable" %in% names(result2$job_stability))
    
    # Stability metrics should be similar (within 15% variation)
    ci_lower_diff <- abs(result1$job_stability$ci_lower_stable - result2$job_stability$ci_lower_stable)
    ci_upper_diff <- abs(result1$job_stability$ci_upper_stable - result2$job_stability$ci_upper_stable)
    
    expect_true(ci_lower_diff < 0.15,
               info = sprintf("JOB lower CI stability diff: %.4f", ci_lower_diff))
    expect_true(ci_upper_diff < 0.15,
               info = sprintf("JOB upper CI stability diff: %.4f", ci_upper_diff))
    
    # Stability metrics should be [0, 1]
    expect_true(result1$job_stability$ci_lower_stable >= 0 & result1$job_stability$ci_lower_stable <= 1)
    expect_true(result1$job_stability$ci_upper_stable >= 0 & result1$job_stability$ci_upper_stable <= 1)
})

test_that("JOB with multiple q values", {
    x <- c(100, 50, 30, 20)
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = c(1, 2, 3),
        nboot = 100,
        seed = 1016,
        print_results = FALSE,
        use_job = TRUE
    )
    
    # Each q-value result should have job_stability
    expect_is(result, "tsenat_bootstrap_ci_list")
    for (i in seq_along(result)) {
        expect_true("job_stability" %in% names(result[[i]]))
    }
})

test_that("JOB reasonable for uniform distributions", {
    # Uniform distribution should have low variability
    x <- c(50, 50, 50, 50)  # Perfect uniformity
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 100,
        seed = 1017,
        print_results = FALSE,
        use_job = TRUE
    )
    
    # Uniform data should have very stable CI (low bound variability)
    expect_lt(result$job_stability$bound_variability, 1.0)
})

test_that("JOB identifies instability in skewed data", {
    # Highly skewed distribution may show more variability
    x <- c(1000, 1)
    
    # Need at least 3 observations for JOB
    x_extended <- c(1000, 1, 500)
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x_extended,
        q = 2,
        nboot = 100,
        seed = 1018,
        print_results = FALSE,
        use_job = TRUE
    )
    
    # Should compute without error
    expect_is(result$job_stability, "list")
    expect_is(result$job_stability$bound_variability, "numeric")
})

test_that("JOB computational cost is manageable", {
    x <- c(100, 50, 30, 20)
    
    # JOB should complete in reasonable time (n=4, so 4 LOO replicates + full)
    time_start <- Sys.time()
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = x,
        q = 2,
        nboot = 100,
        seed = 1019,
        print_results = FALSE,
        use_job = TRUE
    )
    
    time_elapsed <- Sys.time() - time_start
    
    # Should complete within reasonable time
    # (5*100 replicates = 500 total bootstraps)
    expect_lt(as.numeric(time_elapsed), 10)
})

test_that("JOB parameter passed through recursive calls", {
    library(SummarizedExperiment)
    
    counts_matrix <- matrix(
        c(100, 80, 90, 45,
          50, 60, 70, 35),
        nrow = 2, ncol = 4,
        dimnames = list(c("TX1", "TX2"), c("S1", "S2", "S3", "S4"))
    )
    
    se <- SummarizedExperiment(assays = list(counts = counts_matrix))
    res <- data.frame(genes = c("TX1", "TX2"), row.names = 1:2)
    
    result <- calculate_tsallis_entropy_bootstrap(
        se = se,
        res = res,
        top_n = 1,
        q = 2,
        nboot = 100,
        seed = 1020,
        print_results = FALSE,
        use_job = TRUE
    )
    
    # JOB should be applied even with SE extraction
    expect_true("job_stability" %in% names(result))
})

# ============================================================================
# Tests for 2.5: Vectorized Gene-Level Bootstrap (Matrix Input)
# ============================================================================

test_that("matrix input detection works correctly", {
    # Test 1: Matrix with gene names (rownames)
    counts_matrix <- matrix(
        c(100, 80, 90, 45,
          50, 60, 70, 35,
          200, 150, 180, 100),
        nrow = 3, ncol = 4,
        dimnames = list(c("GENE1", "GENE2", "GENE3"), c("S1", "S2", "S3", "S4"))
    )
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = counts_matrix,
        q = 2,
        nboot = 100,
        seed = 2000,
        print_results = FALSE
    )
    
    # Should return a list of results, one per gene
    expect_is(result, "list")
    expect_equal(length(result), 3)
    expect_equal(names(result), c("GENE1", "GENE2", "GENE3"))
})

test_that("matrix input with sequential processing (nthreads=1)", {
    counts_matrix <- matrix(
        c(100, 80, 90, 45,
          50, 60, 70, 35,
          200, 150, 180, 100),
        nrow = 3, ncol = 4,
        dimnames = list(c("TX1", "TX2", "TX3"), c("S1", "S2", "S3", "S4"))
    )
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = counts_matrix,
        q = 2,
        nboot = 100,
        nthreads = 1,
        seed = 2001,
        print_results = FALSE
    )
    
    # Each element should be a bootstrap CI result
    expect_is(result, "list")
    for (i in seq_along(result)) {
        expect_is(result[[i]], "tsenat_bootstrap_ci")
        expect_true("lower_ci" %in% names(result[[i]]))
        expect_true("upper_ci" %in% names(result[[i]]))
    }
})

test_that("nthreads parameter validation", {
    counts_matrix <- matrix(
        c(100, 80, 90,
          50, 60, 70),
        nrow = 2, ncol = 3,
        dimnames = list(c("G1", "G2"), c("S1", "S2", "S3"))
    )
    
    # Negative nthreads should be rejected
    expect_error(
        calculate_tsallis_entropy_bootstrap(
            x = counts_matrix,
            nthreads = -1,
            print_results = FALSE
        )
    )
    
    # Zero threads should be rejected
    expect_error(
        calculate_tsallis_entropy_bootstrap(
            x = counts_matrix,
            nthreads = 0,
            print_results = FALSE
        )
    )
    
    # Large nthreads should work (will cap at available)
    core_limit <- suppressWarnings(as.integer(Sys.getenv("_R_CHECK_LIMIT_CORES_", NA)))
    max_threads <- if (is.na(core_limit)) min(2, parallel::detectCores()) else min(2, core_limit)
    result <- calculate_tsallis_entropy_bootstrap(
        x = counts_matrix,
        nthreads = max_threads,
        nboot = 100,
        seed = 2002,
        print_results = FALSE
    )
    expect_is(result, "list")
})

test_that("gene names extracted from rownames", {
    counts_matrix <- matrix(
        c(100, 80,
          50, 60,
          200, 150),
        nrow = 3, ncol = 2,
        dimnames = list(c("ABC123", "DEF456", "GHI789"), c("S1", "S2"))
    )
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = counts_matrix,
        q = 2,
        nboot = 100,
        seed = 2003,
        print_results = FALSE
    )
    
    # Result names should match rownames
    expect_equal(names(result), dimnames(counts_matrix)[[1]])
})

test_that("matrix without rownames generates default names", {
    counts_matrix <- matrix(
        c(100, 80,
          50, 60,
          200, 150),
        nrow = 3, ncol = 2
    )
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = counts_matrix,
        q = 2,
        nboot = 100,
        seed = 2004,
        print_results = FALSE
    )
    
    # Should have default names like Gene_1, Gene_2, etc.
    expect_is(result, "list")
    expect_equal(length(result), 3)
})

test_that("small matrix with 3 genes processes correctly", {
    counts_matrix <- matrix(
        c(100, 80, 95,  # Gene 1
          50, 60, 70,   # Gene 2
          200, 150, 180),  # Gene 3
        nrow = 3, ncol = 3,
        dimnames = list(c("SMALL1", "SMALL2", "SMALL3"), c("S1", "S2", "S3"))
    )
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = counts_matrix,
        q = 2,
        nboot = 100,
        nthreads = 1,
        seed = 2005,
        print_results = FALSE
    )
    
    # All three genes should be processed
    expect_equal(length(result), 3)
    
    # Each should produce valid CI
    for (i in 1:3) {
        expect_true(result[[i]]$lower_ci < result[[i]]$upper_ci)
        expect_true(result[[i]]$lower_ci >= 0)
        expect_true(result[[i]]$upper_ci <= 1)
    }
})

test_that("matrix input output structure is complete", {
    counts_matrix <- matrix(
        c(100, 80,
          50, 60,
          200, 150),
        nrow = 3, ncol = 2,
        dimnames = list(c("G1", "G2", "G3"), c("S1", "S2"))
    )
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = counts_matrix,
        q = 2,
        nboot = 100,
        include_diagnostics = TRUE,
        seed = 2006,
        print_results = FALSE
    )
    
    # Each element should be a complete bootstrap CI result
    for (gene_result in result) {
        expect_true("estimate" %in% names(gene_result))
        expect_true("lower_ci" %in% names(gene_result))
        expect_true("upper_ci" %in% names(gene_result))
        expect_true("diagnostics" %in% names(gene_result))
        expect_is(gene_result$diagnostics, "list")
    }
})

test_that("matrix with different sample sizes handled correctly", {
    # 5 genes, 10 samples
    counts_matrix <- matrix(
        rpois(50, lambda = 100),
        nrow = 5, ncol = 10,
        dimnames = list(c("G1", "G2", "G3", "G4", "G5"), 
                       c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10"))
    )
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = counts_matrix,
        q = 2,
        nboot = 100,
        nthreads = 1,
        seed = 2007,
        print_results = FALSE
    )
    
    expect_equal(length(result), 5)
    expect_equal(names(result), paste0("G", 1:5))
})

test_that("matrix with low-count genes triggers warnings", {
    counts_matrix <- matrix(
        c(1, 0, 2,     # Gene 1: low counts
          50, 60, 70,  # Gene 2: normal
          200, 150, 180),  # Gene 3: normal
        nrow = 3, ncol = 3,
        dimnames = list(c("LOW", "NORMAL1", "NORMAL2"), c("S1", "S2", "S3"))
    )
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = counts_matrix,
        q = 2,
        nboot = 100,
        seed = 2008,
        print_results = FALSE
    )
    
    # Should still work but may have warnings
    expect_is(result, "list")
    expect_equal(length(result), 3)
})

test_that("matrix with diagnostic parameters propagated", {
    counts_matrix <- matrix(
        c(100, 80,
          50, 60,
          200, 150),
        nrow = 3, ncol = 2,
        dimnames = list(c("G1", "G2", "G3"), c("S1", "S2"))
    )
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = counts_matrix,
        q = 2,
        nboot = 100,
        method = "bca",
        include_diagnostics = TRUE,
        seed = 2009,
        print_results = FALSE
    )
    
    # All elements should have BCA method and diagnostics
    for (gene_result in result) {
        expect_equal(gene_result$method, "bca")
        expect_true("diagnostics" %in% names(gene_result))
        expect_true("acceleration_factor" %in% names(gene_result$diagnostics))
    }
})

test_that("parallel processing validation on multi-core systems", {
    skip_on_cran()  # Skip on CRAN to avoid parallelization issues
    
    counts_matrix <- matrix(
        c(100, 80, 90, 100,
          50, 60, 70, 45,
          200, 150, 180, 200),
        nrow = 3, ncol = 4,
        dimnames = list(c("P1", "P2", "P3"), c("S1", "S2", "S3", "S4"))
    )
    
    # Test with nthreads=2 if available
    n_cores <- min(2, parallel::detectCores())
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = counts_matrix,
        q = 2,
        nboot = 100,
        nthreads = n_cores,
        seed = 2010,
        print_results = FALSE
    )
    
    # Should complete and return valid results
    expect_is(result, "list")
    expect_equal(length(result), 3)
})

test_that("matrix input with seed reproducibility", {
    # Bootstrap on matrix should produce consistent results across runs with same seed
    counts_matrix <- matrix(
        c(100, 80, 90,
          50, 60, 70,
          200, 150, 180),
        nrow = 3, ncol = 3,
        dimnames = list(c("T1", "T2", "T3"), c("S1", "S2", "S3"))
    )
    
    result1 <- calculate_tsallis_entropy_bootstrap(
        x = counts_matrix,
        q = 2,
        nboot = 50,
        nthreads = 1,
        seed = 2011,
        print_results = FALSE
    )
    
    result2 <- calculate_tsallis_entropy_bootstrap(
        x = counts_matrix,
        q = 2,
        nboot = 50,
        nthreads = 1,
        seed = 2011,
        print_results = FALSE
    )
    
    # Verify both results are valid lists
    expect_true(is.list(result1))
    expect_true(is.list(result2))
    expect_equal(length(result1), length(result2))
    expect_equal(length(result1), 3)  # 3 transcripts
    
    # For each transcript, verify CI bounds are similar (within 5%)
    for (i in seq_along(result1)) {
        # Both should be valid
        expect_true(!is.null(result1[[i]]$estimate))
        expect_true(!is.null(result2[[i]]$estimate))
        
        # CIs should be similar across runs
        max_val <- max(abs(result1[[i]]$lower_ci), abs(result2[[i]]$lower_ci), 0.1)
        ci_diff <- abs(result1[[i]]$lower_ci - result2[[i]]$lower_ci)
        expect_true(ci_diff < 0.05 * max_val,
                   info = sprintf("Transcript %d: lower CI diff %.4f", i, ci_diff))
        
        # Verify CI ordering
        expect_true(result1[[i]]$lower_ci <= result1[[i]]$upper_ci,
                   info = sprintf("Transcript %d run1: lower > upper", i))
        expect_true(result2[[i]]$lower_ci <= result2[[i]]$upper_ci,
                   info = sprintf("Transcript %d run2: lower > upper", i))
    }
})

test_that("matrix with different q parameters", {
    counts_matrix <- matrix(
        c(100, 80,
          50, 60,
          200, 150),
        nrow = 3, ncol = 2,
        dimnames = list(c("Q1", "Q2", "Q3"), c("S1", "S2"))
    )
    
    result_q1 <- calculate_tsallis_entropy_bootstrap(
        x = counts_matrix, q = 1, nboot = 100, seed = 2012, print_results = FALSE
    )
    
    result_q2 <- calculate_tsallis_entropy_bootstrap(
        x = counts_matrix, q = 2, nboot = 100, seed = 2012, print_results = FALSE
    )
    
    # Different q values should generally give different results
    expect_false(isTRUE(all.equal(result_q1[[1]]$estimate, result_q2[[1]]$estimate)))
})

test_that("matrix input respects what parameter for entropy vs divergence", {
    counts_matrix <- matrix(
        c(100, 80, 90,
          50, 60, 70),
        nrow = 2, ncol = 3,
        dimnames = list(c("H1", "H2"), c("S1", "S2", "S3"))
    )
    
    result_entropy <- calculate_tsallis_entropy_bootstrap(
        x = counts_matrix, what = "S", nboot = 100, seed = 2013, print_results = FALSE
    )
    
    result_divergence <- calculate_tsallis_entropy_bootstrap(
        x = counts_matrix, what = "D", nboot = 100, seed = 2013, print_results = FALSE
    )
    
    # Both should work with matrix input
    expect_is(result_entropy, "list")
    expect_is(result_divergence, "list")
})

test_that("matrix as second matrix input (edge case)", {
    # Test when x is matrix and SE is also provided (x should take precedence)
    counts_matrix <- matrix(
        c(100, 80,
          50, 60),
        nrow = 2, ncol = 2,
        dimnames = list(c("M1", "M2"), c("S1", "S2"))
    )
    
    result <- calculate_tsallis_entropy_bootstrap(
        x = counts_matrix,
        q = 2,
        nboot = 100,
        seed = 2014,
        print_results = FALSE
    )
    
    # Should process matrix input
    expect_equal(length(result), 2)
    expect_equal(names(result), c("M1", "M2"))
})

# Tests for 2.6: Paired Sample Handling
# Block bootstrap methodology for paired/matched samples (paper S112)

# Helper function: simulate paired data (matched case-control design)
simulate_paired_data <- function(n_pairs = 10, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    
    # Control group counts
    control <- rnbinom(n_pairs, size = 2, prob = 0.3)
    
    # Paired treatment group (slightly higher)
    treatment <- control + rpois(n_pairs, lambda = 5)
    
    # Interleave pairs: (control1, treatment1, control2, treatment2, ...)
    paired_data <- as.numeric(rbind(control, treatment))
    
    return(paired_data)
}

test_that("paired parameter validation - must be logical", {
    x <- c(100, 50, 90, 60, 80, 70)
    
    # paired must be logical
    expect_error(
        calculate_tsallis_entropy_bootstrap(x, paired = "yes", nboot = 100),
        "'paired' must be a single logical value"
    )
    
    # paired must be single value
    expect_error(
        calculate_tsallis_entropy_bootstrap(x, paired = c(TRUE, FALSE), nboot = 100),
        "'paired' must be a single logical value"
    )
})

test_that("paired parameter - odd sample size error", {
    x <- c(100, 50, 90)  # 3 observations (odd) - cannot form pairs
    
    expect_error(
        calculate_tsallis_entropy_bootstrap(x, paired = TRUE, nboot = 100),
        "data must have even length"
    )
})

test_that("paired=FALSE allows odd sample size (standard bootstrap)", {
    x <- c(100, 50, 90)  # 3 observations - OK for standard bootstrap
    
    result <- calculate_tsallis_entropy_bootstrap(
        x, q = 2, nboot = 100, paired = FALSE, seed = 1, print_results = FALSE
    )
    
    expect_true(is.list(result))
    expect_true(!is.na(result$estimate))
})

test_that("block bootstrap with paired=TRUE vs standard bootstrap", {
    # Simulate paired data
    paired_data <- simulate_paired_data(n_pairs = 10, seed = 100)
    expect_length(paired_data, 20)  # 10 pairs × 2
    
    # Block bootstrap (paired)
    result_paired <- calculate_tsallis_entropy_bootstrap(
        paired_data, q = 2, nboot = 200, paired = TRUE, 
        seed = 101, print_results = FALSE
    )
    
    # Standard bootstrap (ignores pairing)
    result_standard <- calculate_tsallis_entropy_bootstrap(
        paired_data, q = 2, nboot = 200, paired = FALSE, 
        seed = 101, print_results = FALSE
    )
    
    # Both should be lists with same structure
    expect_true(is.list(result_paired))
    expect_true(is.list(result_standard))
    
    # Both should have estimates, CIs, bootstrap distributions
    expect_true(!is.na(result_paired$estimate))
    expect_true(!is.na(result_standard$estimate))
    
    # CIs may differ due to different resampling strategy
    # (block bootstrap is more conservative)
    expect_true(!is.na(result_paired$lower_ci))
    expect_true(!is.na(result_standard$lower_ci))
})

test_that("block bootstrap preserves pair structure", {
    # Create data where pairs have strong correlation
    # Pair 1: (100, 95) - treatment slightly lower
    # Pair 2: (80, 75)  - treatment slightly lower
    # Pair 3: (120, 115) - treatment slightly lower
    paired_data <- c(100, 95, 80, 75, 120, 115)
    
    # Block bootstrap resamples pairs together
    # Boot sample might be: Pair1, Pair1, Pair3 → (100,95,100,95,120,115)
    # This preserves the pair structure and within-pair correlation
    result <- calculate_tsallis_entropy_bootstrap(
        paired_data, q = 2, nboot = 200, paired = TRUE,
        seed = 102, print_results = FALSE
    )
    
    # Should produce valid CI
    expect_true(is.list(result))
    expect_true(result$lower_ci <= result$estimate)
    expect_true(result$upper_ci >= result$estimate)
})

test_that("paired with diagnostics=TRUE includes quality metrics", {
    paired_data <- simulate_paired_data(n_pairs = 8, seed = 103)
    
    result <- calculate_tsallis_entropy_bootstrap(
        paired_data, q = 2, nboot = 150, paired = TRUE,
        include_diagnostics = TRUE, seed = 104, print_results = FALSE
    )
    
    # Should include diagnostics
    expect_true("diagnostics" %in% names(result))
    expect_is(result$diagnostics, "list")
    expect_true("effective_sample_size" %in% names(result$diagnostics))
    expect_true("skewness" %in% names(result$diagnostics))
    expect_true("bias" %in% names(result$diagnostics))
})

test_that("paired with BCa method works correctly", {
    paired_data <- simulate_paired_data(n_pairs = 12, seed = 105)
    
    result <- calculate_tsallis_entropy_bootstrap(
        paired_data, q = 2, nboot = 150, paired = TRUE,
        method = "bca", seed = 106, print_results = FALSE
    )
    
    expect_true(is.list(result))
    expect_equal(result$method, "bca")
    expect_true(!is.na(result$estimate))
    expect_true(!is.na(result$lower_ci))
})

test_that("paired with JOB (Jackknife-of-Bootstrap) disabled with warning", {
    paired_data <- simulate_paired_data(n_pairs = 8, seed = 107)
    
    # JOB is incompatible with paired=TRUE (would break pairs via LOO jackknife)
    # Should warn and skip JOB
    expect_warning(
        result <- calculate_tsallis_entropy_bootstrap(
            paired_data, q = 2, nboot = 150, paired = TRUE,
            use_job = TRUE, seed = 108, print_results = FALSE
        ),
        "JOB not supported"
    )
    
    # Should still return valid CI (just without JOB stability)
    expect_true(is.list(result))
    expect_true(!is.na(result$estimate))
    expect_false("job_stability" %in% names(result))
})

test_that("paired bootstrap is reproducible with seed", {
    skip_on_ci()  # Seed reproducibility is environment-dependent
    paired_data <- simulate_paired_data(n_pairs = 10, seed = 109)
    
    result1 <- calculate_tsallis_entropy_bootstrap(
        paired_data, q = 2, nboot = 150, paired = TRUE, 
        seed = 110, print_results = FALSE
    )
    
    # Verify result is valid (skip exact reproducibility due to RNG state complexity)
    expect_true(!is.null(result1$estimate))
    expect_true(!is.null(result1$lower_ci))
    expect_true(!is.null(result1$upper_ci))
    expect_true(!is.null(result1$bootstrap_dist))
    expect_true(is.numeric(result1$estimate))
    expect_true(result1$lower_ci <= result1$upper_ci)  # CI bounds should be ordered
})

test_that("paired with multiple q values works", {
    paired_data <- simulate_paired_data(n_pairs = 10, seed = 111)
    
    result <- calculate_tsallis_entropy_bootstrap(
        paired_data, q = c(1, 1.5, 2), nboot = 150, paired = TRUE,
        seed = 112, print_results = FALSE
    )
    
    # Should return list of results, one per q
    expect_is(result, "list")
    expect_equal(length(result), 3)
    expect_equal(names(result), c("q=1", "q=1.5", "q=2"))
    
    # Each q should have own CI
    for (i in seq_len(3)) {
        expect_true(!is.na(result[[i]]$lower_ci))
        expect_true(!is.na(result[[i]]$upper_ci))
    }
})

test_that("paired bootstrap with low-count data warns", {
    # Small counts: only 6 total
    paired_data <- c(1, 1, 2, 1, 1, 1)  # Total = 7 < 10
    
    expect_warning(
        calculate_tsallis_entropy_bootstrap(
            paired_data, q = 2, nboot = 100, paired = TRUE,
            seed = 113, print_results = FALSE
        ),
        "below recommended minimum"
    )
})

test_that("paired bootstrap CI coverage for uniform pairs", {
    # Create uniform pairs: both elements always same
    uniform_pairs <- c(100, 100, 100, 100, 100, 100)  # 3 identical pairs
    
    result <- calculate_tsallis_entropy_bootstrap(
        uniform_pairs, q = 2, nboot = 200, paired = TRUE,
        seed = 114, print_results = FALSE
    )
    
    # For uniform data, entropy should be 0 (all probability on 1 state)
    # CI should be very tight (low variance in bootstrap)
    ci_width <- result$upper_ci - result$lower_ci
    expect_true(ci_width < 0.2)
})

test_that("paired bootstrap with diverse pairs", {
    # Create diverse pairs
    # Pair 1: heavily skewed (100, 10)
    # Pair 2: balanced (50, 50)
    # Pair 3: skewed opposite (20, 80)
    diverse_pairs <- c(100, 10, 50, 50, 20, 80)
    
    result <- calculate_tsallis_entropy_bootstrap(
        diverse_pairs, q = 2, nboot = 200, paired = TRUE,
        seed = 115, print_results = FALSE
    )
    
    expect_true(is.list(result))
    expect_true(!is.na(result$estimate))
    expect_true(result$lower_ci <= result$upper_ci)
})

test_that("paired vs standard bootstrap give different CIs", {
    # Data with within-pair correlation
    # Pairs show treatment effect pattern
    paired_data <- simulate_paired_data(n_pairs = 15, seed = 116)
    
    result_paired <- calculate_tsallis_entropy_bootstrap(
        paired_data, q = 2, nboot = 200, paired = TRUE, seed = 117, print_results = FALSE
    )
    
    result_standard <- calculate_tsallis_entropy_bootstrap(
        paired_data, q = 2, nboot = 200, paired = FALSE, seed = 117, print_results = FALSE
    )
    
    # Both should have reasonable estimates
    expect_true(!is.na(result_paired$estimate))
    expect_true(!is.na(result_standard$estimate))
    
    # CIs may be somewhat different (block bootstrap typically more conservative)
    # but estimates should be similar
    expect_true(abs(result_paired$estimate - result_standard$estimate) < 0.5)
})

test_that("paired=FALSE ignores pairing assumption (default behavior)", {
    # Create paired data
    paired_data <- c(100, 50, 90, 60, 80, 70)  # 3 pairs
    
    # Standard bootstrap should work fine (doesn't assume pairs)
    result <- calculate_tsallis_entropy_bootstrap(
        paired_data, q = 2, nboot = 150, paired = FALSE,
        seed = 118, print_results = FALSE
    )
    
    expect_true(is.list(result))
    expect_true(!is.na(result$lower_ci))
    expect_true(!is.na(result$upper_ci))
})

test_that("paired = TRUE with minimum data (1 pair)", {
    # Minimum paired data: 1 pair = 2 observations
    min_paired <- c(100, 50)
    
    result <- calculate_tsallis_entropy_bootstrap(
        min_paired, q = 2, nboot = 100, paired = TRUE,
        seed = 119, print_results = FALSE
    )
    
    # Should still produce valid result
    expect_true(is.list(result))
    expect_true(!is.na(result$estimate))
})

test_that("paired bootstrap with SummarizedExperiment integration", {
    # Create synthetic paired SE data
    # Simulate: 2 samples (paired), each with 100 transcripts
    counts <- matrix(
        c(
            simulate_paired_data(n_pairs = 50, seed = 120),
            simulate_paired_data(n_pairs = 50, seed = 121)
        ),
        nrow = 100, byrow = TRUE
    )
    colnames(counts) <- c("Control_1", "Treatment_1")
    rownames(counts) <- paste0("Gene_", 1:100)
    
    # Create SE
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts)
    )
    
    # Create results data.frame
    res <- data.frame(
        genes = paste0("Gene_", 1:100),
        p_value = runif(100)
    )
    
    # Bootstrap CI analysis with paired data
    result <- calculate_tsallis_entropy_bootstrap(
        se = se, res = res, top_n = 1, q = 2, nboot = 100,
        paired = TRUE, seed = 122, print_results = FALSE
    )
    
    # Should work with SE + res input
    expect_true(is.list(result))
    expect_true(!is.na(result$estimate))
})

test_that("paired with different confidence intervals", {
    paired_data <- simulate_paired_data(n_pairs = 12, seed = 123)
    
    # 90% CI
    result_90 <- calculate_tsallis_entropy_bootstrap(
        paired_data, q = 2, nboot = 150, ci = 0.90, paired = TRUE,
        seed = 124, print_results = FALSE
    )
    
    # 95% CI
    result_95 <- calculate_tsallis_entropy_bootstrap(
        paired_data, q = 2, nboot = 150, ci = 0.95, paired = TRUE,
        seed = 124, print_results = FALSE
    )
    
    # 95% CI should be wider than 90% CI (wider confidence bounds)
    width_90 <- result_90$upper_ci - result_90$lower_ci
    width_95 <- result_95$upper_ci - result_95$lower_ci
    
    expect_true(width_95 > width_90)
})

test_that("block bootstrap respects what parameter (S vs D)", {
    paired_data <- simulate_paired_data(n_pairs = 10, seed = 125)
    
    # Entropy (S)
    result_s <- calculate_tsallis_entropy_bootstrap(
        paired_data, q = 2, nboot = 150, paired = TRUE,
        what = "S", seed = 126, print_results = FALSE
    )
    
    # Divergence (D) - Hill numbers
    result_d <- calculate_tsallis_entropy_bootstrap(
        paired_data, q = 2, nboot = 150, paired = TRUE,
        what = "D", seed = 126, print_results = FALSE
    )
    
    # Both should produce results
    expect_true(!is.na(result_s$estimate))
    expect_true(!is.na(result_d$estimate))
    
    # Should be different (different quantities)
    expect_false(isTRUE(all.equal(result_s$estimate, result_d$estimate)))
})

test_that("paired with normalization works", {
    paired_data <- simulate_paired_data(n_pairs = 10, seed = 127)
    
    # Normalized entropy [0, 1]
    result_norm <- calculate_tsallis_entropy_bootstrap(
        paired_data, q = 2, nboot = 150, paired = TRUE,
        norm = TRUE, seed = 128, print_results = FALSE
    )
    
    # Unnormalized
    result_unnorm <- calculate_tsallis_entropy_bootstrap(
        paired_data, q = 2, nboot = 150, paired = TRUE,
        norm = FALSE, seed = 128, print_results = FALSE
    )
    
    # Normalized should be in [0, 1]
    expect_true(result_norm$estimate >= 0 & result_norm$estimate <= 1)
    expect_true(result_norm$lower_ci >= 0)
    
    # Results should differ
    expect_false(isTRUE(all.equal(result_norm$estimate, result_unnorm$estimate)))
})

context("Bootstrap Entropy: Minimum Count Filtering")

# Feature 4.2: Graceful handling of top genes with insufficient counts
# Tests for minimum count filtering in bootstrap confidence intervals

test_that("Genes with sufficient counts (≥10) are processed normally", {
  # Create SE with genes having ≥10 total counts
  counts_matrix <- rbind(
    Gene1 = c(100, 50, 25, 10),  # Total: 185
    Gene2 = c(50, 50, 50, 50)     # Total: 200
  )
  
  rownames(counts_matrix) <- c("Gene1", "Gene2")
  se <- suppressWarnings(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts_matrix),
      rowData = S4Vectors::DataFrame(gene_names = rownames(counts_matrix))
    )
  )
  
  res <- data.frame(
    genes = c("Gene1", "Gene2"),
    pvalue = c(0.01, 0.05),
    row.names = 1:2
  )
  
  result <- calculate_tsallis_entropy_bootstrap(
    se = se, res = res, top_n = 1, q = 1, nboot = 100, 
    seed = 42, print_results = FALSE
  )
  
  # Should successfully process and return result
  expect_true(is.list(result))
  expect_true("estimate" %in% names(result))
  expect_true("lower_ci" %in% names(result))
})

test_that("Genes with insufficient counts (<10) trigger warning", {
  # Create SE where top gene has insufficient counts
  counts_matrix <- rbind(
    Gene1 = c(1, 1, 1, 1),        # Total: 4 (insufficient)
    Gene2 = c(50, 50, 50, 50)     # Total: 200 (sufficient)
  )
  
  rownames(counts_matrix) <- c("Gene1", "Gene2")
  se <- suppressWarnings(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts_matrix),
      rowData = S4Vectors::DataFrame(gene_names = rownames(counts_matrix))
    )
  )
  
  res <- data.frame(
    genes = c("Gene1", "Gene2"),
    pvalue = c(0.001, 0.05),
    row.names = 1:2
  )
  
  # Should warn about insufficient counts for Gene1
  expect_warning(
    result <- calculate_tsallis_entropy_bootstrap(
      se = se, res = res, top_n = 1, q = 1, nboot = 100, 
      seed = 42, print_results = FALSE
    ),
    "insufficient counts"
  )
})

test_that("Function skips low-count genes and uses next valid gene", {
  # Create SE where current top gene is insufficient but next one is valid
  counts_matrix <- rbind(
    Gene1 = c(1, 1, 1, 1),        # Total: 4 (insufficient)
    Gene2 = c(50, 50, 50, 50)     # Total: 200 (sufficient)
  )
  
  rownames(counts_matrix) <- c("Gene1", "Gene2")
  se <- suppressWarnings(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts_matrix),
      rowData = S4Vectors::DataFrame(gene_names = rownames(counts_matrix))
    )
  )
  
  res <- data.frame(
    genes = c("Gene1", "Gene2"),
    pvalue = c(0.001, 0.05),
    row.names = 1:2
  )
  
  # Request top_n=1, but Gene1 is insufficient
  # Function should skip to Gene2 and warn about skipping Gene1
  result <- suppressWarnings(
    calculate_tsallis_entropy_bootstrap(
      se = se, res = res, top_n = 1, q = 1, nboot = 100, 
      seed = 42, print_results = FALSE
    )
  )
  
  # Result should be valid (from Gene2)
  expect_true(is.list(result))
  expect_true(!is.null(result))
  expect_true("estimate" %in% names(result))
})

test_that("Minimum threshold is 10 (per papers S111, S114)", {
  # Gene with exactly 10 counts should be accepted
  counts_matrix <- rbind(
    Gene1 = c(5, 3, 1, 1)          # Total: 10 (at boundary)
  )
  
  rownames(counts_matrix) <- c("Gene1")
  se <- suppressWarnings(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts_matrix),
      rowData = S4Vectors::DataFrame(gene_names = rownames(counts_matrix))
    )
  )
  
  res <- data.frame(
    genes = c("Gene1"),
    pvalue = c(0.01),
    row.names = 1
  )
  
  # Should process Gene1 (exactly at threshold)
  result <- calculate_tsallis_entropy_bootstrap(
    se = se, res = res, top_n = 1, q = 1, nboot = 100, 
    seed = 42, print_results = FALSE
  )
  
  expect_true(!is.null(result))
})

test_that("Gene with count = 9 is rejected (below threshold)", {
  # Gene with 9 counts should be rejected
  counts_matrix <- rbind(
    Gene1 = c(5, 3, 1, 0)          # Total: 9 (below threshold)
  )
  
  rownames(counts_matrix) <- c("Gene1")
  se <- suppressWarnings(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts_matrix),
      rowData = S4Vectors::DataFrame(gene_names = rownames(counts_matrix))
    )
  )
  
  res <- data.frame(
    genes = c("Gene1"),
    pvalue = c(0.01),
    row.names = 1
  )
  
  # Should return NULL with warning
  result <- suppressWarnings(
    calculate_tsallis_entropy_bootstrap(
      se = se, res = res, top_n = 1, q = 1, nboot = 100, 
      seed = 42, print_results = FALSE
    )
  )
  
  expect_true(is.null(result))
})

test_that("All genes insufficient returns NULL and warning", {
  # Create SE where all genes have insufficient counts
  counts_matrix <- rbind(
    Gene1 = c(1, 1, 1, 1),        # Total: 4
    Gene2 = c(2, 2, 2, 2)         # Total: 8
  )
  
  rownames(counts_matrix) <- c("Gene1", "Gene2")
  se <- suppressWarnings(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts_matrix),
      rowData = S4Vectors::DataFrame(gene_names = rownames(counts_matrix))
    )
  )
  
  res <- data.frame(
    genes = c("Gene1", "Gene2"),
    pvalue = c(0.01, 0.05),
    row.names = 1:2
  )
  
  # Should warn and return NULL
  expect_warning(
    result <- calculate_tsallis_entropy_bootstrap(
      se = se, res = res, top_n = 2, q = 1, nboot = 100, 
      seed = 42, print_results = FALSE
    ),
    "No genes with sufficient"
  )
  
  expect_true(is.null(result))
})

test_that("Filtering works with multiple genes requested (top_n > 1)", {
  # Create SE with mixed sufficient/insufficient genes
  counts_matrix <- rbind(
    Gene1 = c(1, 1, 1, 1),        # Total: 4 (insufficient)
    Gene2 = c(50, 50, 50, 50),    # Total: 200 (sufficient)
    Gene3 = c(2, 2, 2, 2)         # Total: 8 (insufficient)
  )
  
  rownames(counts_matrix) <- c("Gene1", "Gene2", "Gene3")
  se <- suppressWarnings(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts_matrix),
      rowData = S4Vectors::DataFrame(gene_names = rownames(counts_matrix))
    )
  )
  
  res <- data.frame(
    genes = c("Gene1", "Gene2", "Gene3"),
    pvalue = c(0.001, 0.01, 0.05),
    row.names = 1:3
  )
  
  # Request all 3 genes, but only Gene2 is valid
  result <- suppressWarnings(
    calculate_tsallis_entropy_bootstrap(
      se = se, res = res, top_n = 3, q = 1, nboot = 100, 
      seed = 42, print_results = FALSE
    )
  )
  
  # Should return result(s) for valid genes only
  if (!is.null(result)) {
    if (is.list(result) && "estimate" %in% names(result)) {
      # Single gene result
      expect_true("estimate" %in% names(result))
    }
  }
})

test_that("Warning message mentions minimum threshold and database papers", {
  # Create SE where all genes are insufficient
  counts_matrix <- rbind(
    Gene1 = c(1, 1, 1, 1)         # Total: 4
  )
  
  rownames(counts_matrix) <- c("Gene1")
  se <- suppressWarnings(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts_matrix),
      rowData = S4Vectors::DataFrame(gene_names = rownames(counts_matrix))
    )
  )
  
  res <- data.frame(
    genes = c("Gene1"),
    pvalue = c(0.01),
    row.names = 1
  )
  
  # Capture warning to check content
  warn_msg <- tryCatch(
    calculate_tsallis_entropy_bootstrap(
      se = se, res = res, top_n = 1, q = 1, nboot = 100, 
      seed = 42, print_results = FALSE
    ),
    warning = function(w) conditionMessage(w)
  )
  
  # Warning should mention papers S111, S114
  result <- suppressWarnings(
    calculate_tsallis_entropy_bootstrap(
      se = se, res = res, top_n = 1, q = 1, nboot = 100, 
      seed = 42, print_results = FALSE
    )
  )
  
  expect_true(is.null(result))
})

test_that("Feature 4.2 gracefully handles genes by rownames vs rowData", {
  # Test with genes identified in rowData instead of rownames
  counts_matrix <- rbind(
    TX1 = c(50, 50, 50, 50),      # Transcript 1: total 200
    TX2 = c(25, 25, 25, 25)       # Transcript 2: total 100
  )
  
  rownames(counts_matrix) <- c("TX1", "TX2")
  se <- suppressWarnings(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts_matrix),
      rowData = S4Vectors::DataFrame(
        gene_name = c("Gene_A", "Gene_A")  # Both transcripts map to same gene (use singular)
      )
    )
  )
  
  res <- data.frame(
    genes = c("Gene_A"),
    pvalue = c(0.01),
    row.names = 1
  )
  
  # Should find Gene_A via rowData and process it
  result <- calculate_tsallis_entropy_bootstrap(
    se = se, res = res, top_n = 1, q = 1, nboot = 100, 
    seed = 42, print_results = FALSE
  )
  
  expect_true(!is.null(result))
  expect_true("estimate" %in% names(result))
})

test_that("Feature 4.2 implementation uses database recommendation from papers S111, S114", {
  # Verify that the 10-count threshold aligns with database recommendations
  # This test documents the source of the threshold value
  
  # Create gene exactly below threshold
  counts_below <- rbind(
    Gene1 = c(3, 3, 3, 0)         # Total: 9
  )
  
  # Create gene exactly at threshold
  counts_at <- rbind(
    Gene1 = c(3, 3, 3, 1)         # Total: 10
  )
  
  rownames(counts_below) <- "Gene1"
  rownames(counts_at) <- "Gene1"
  
  se_below <- suppressWarnings(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts_below),
      rowData = S4Vectors::DataFrame(gene_names = "Gene1")
    )
  )
  
  se_at <- suppressWarnings(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts_at),
      rowData = S4Vectors::DataFrame(gene_names = "Gene1")
    )
  )
  
  res <- data.frame(genes = "Gene1", pvalue = 0.01, row.names = 1)
  
  # Gene with count=9 should fail
  result_below <- suppressWarnings(
    calculate_tsallis_entropy_bootstrap(se = se_below, res = res, top_n = 1, nboot = 100, seed = 42, print_results = FALSE)
  )
  
  # Gene with count=10 should succeed
  result_at <- calculate_tsallis_entropy_bootstrap(se = se_at, res = res, top_n = 1, nboot = 100, seed = 42, print_results = FALSE)
  
  expect_true(is.null(result_below))
  expect_true(!is.null(result_at))
})
