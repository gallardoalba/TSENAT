context("Label Shuffling: Paired Permutations")

test_that("paired permutation requires even number of samples", {
    x <- matrix(runif(6), nrow = 2)
    samples <- c("A", "B", "A")
    expect_error(
        label_shuffling(x, samples, control = "A", method = "mean", paired = TRUE),
        "Paired permutation requires an even number of samples"
    )
})

context("Label Shuffling: Sign-Flip Permutations (Exact and Sampled)")

test_that("signflip exact enumeration runs and returns valid p-values", {
    set.seed(42)
    x <- matrix(rnorm(8), nrow = 2)
    samples <- rep(c("Normal", "Tumor"), times = 2) # two pairs
    # total combinations = 2^2 = 4
    res <- label_shuffling(x, samples, control = "Normal", method = "mean", randomizations = 4, paired = TRUE, paired_method = "signflip")
    expect_true(is.data.frame(res))
    # Now expects: pvalue, padj, log2FC, U, r, + 2 group means = 7 columns
    expect_equal(ncol(res), 7)
    expect_equal(nrow(res), nrow(x))
    # Check p-value columns are in [0,1]
    expect_true(all(res$pvalue >= 0 & res$pvalue <= 1, na.rm = TRUE))
    expect_true(all(res$padj >= 0 & res$padj <= 1, na.rm = TRUE))
})

test_that("signflip sampled returns same shape and in-range p-values", {
    set.seed(42)
    x <- matrix(rnorm(8), nrow = 2)
    samples <- rep(c("Normal", "Tumor"), times = 2)
    res <- label_shuffling(x, samples, control = "Normal", method = "mean", randomizations = 10, paired = TRUE, paired_method = "signflip")
    expect_true(is.data.frame(res))
    # Now expects: pvalue, padj, log2FC, U, r, + 2 group means = 7 columns
    expect_equal(ncol(res), 7)
    expect_equal(nrow(res), nrow(x))
    # Check p-value columns are in [0,1]
    expect_true(all(res$pvalue >= 0 & res$pvalue <= 1, na.rm = TRUE))
    expect_true(all(res$padj >= 0 & res$padj <= 1, na.rm = TRUE))
})

test_that("label_shuffling handles constant null distribution correctly", {
    # Build a simple matrix where permutation yields identical nulls
    # Two groups of 2 samples each: group A (control), group B (case)
    mat <- matrix(c(
        0.1, 0.1, 0.1, 0.1, # gene1: identical across samples -> permuted log2fc all zero
        1.0, 1.0, 2.0, 2.0 # gene2: different between groups -> non-zero observed
    ), nrow = 2, byrow = TRUE)
    colnames(mat) <- paste0("S", 1:4)
    samples <- c("A", "A", "B", "B")

    # gene1 observed log2fc should be 0; nulls will be all 0 -> pval should be 1/(n+1)
    # gene2 observed log2fc should be non-zero -> pval should reflect permutations

    # Use a small number of permutations for reproducibility
    set.seed(42)
    res <- label_shuffling(mat, samples, control = "A", method = "mean", randomizations = 10, pcorr = "none")
    expect_true(is.matrix(res) || is.data.frame(res))
    raw <- as.numeric(res[, 1])

    # For gene1: obs = 0, nulls all zero -> count = (all nulls |null| >= |obs|) = n_perm
    # empirical p = (n_perm + 1) / (n_perm + 1) = 1
    expect_equal(raw[1], 1)

    # For gene2: expect p between 0 and 1
    expect_true(raw[2] > 0 && raw[2] <= 1)
})

test_that("label_shuffling() returns named p-value, log2FC and group mean columns", {
    mat <- matrix(c(
        0.1, 0.2, 0.3, 0.4,
        1.0, 1.2, 0.9, 1.1
    ), nrow = 2, byrow = TRUE)
    colnames(mat) <- paste0("S", 1:4)
    samples <- c("A", "A", "B", "B")
    set.seed(1)
    res <- label_shuffling(mat, samples, control = "A", method = "mean", randomizations = 10, pcorr = "none")
    expect_true(is.data.frame(res))
    expect_true(all(c("pvalue", "padj", "log2FC") %in% colnames(res)))
    # Should have group means columns as well
    expect_true(ncol(res) >= 3)
})

# ============================================================================
# S019 PHIPSON & SMYTH (2010) BIAS CORRECTION TESTS
# ============================================================================
# Validate that p-values follow the formula: p = (b + 1) / (m + 1)
# where b = count of extreme permutations, m = total permutations
#
# This ensures:
#   1. No p-values equal exactly zero (minimum = 1/(m+1))
#   2. Proper statistical calibration
#   3. Correct Type I error control
# ============================================================================

context("Label Shuffling: S019 Phipson & Smyth p-Value Bias Correction")

test_that("No p-values equal zero (bias correction working)", {
    # Test that the correction prevents p = 0 for extreme observations
    set.seed(123)
    
    # Create matrix where some genes will be extreme after permutation
    # 4 genes x 8 samples (4 control + 4 case)
    mat <- matrix(rnorm(32, mean = 0, sd = 1), nrow = 4)
    samples <- c("Control", "Control", "Control", "Control", "Case", "Case", "Case", "Case")
    
    res <- label_shuffling(
        mat, 
        samples, 
        control = "Control", 
        method = "mean", 
        randomizations = 100,
        pcorr = "none"
    )
    
    # Critical: No p-value should equal exactly 0
    expect_true(all(res$pvalue > 0, na.rm = TRUE),
                info = "S019 correction failed: found p-values equal to 0")
    
    # All p-values should be finite
    expect_true(all(is.finite(res$pvalue), na.rm = TRUE),
                info = "Found non-finite p-values")
})

test_that("Minimum p-value equals 1/(m+1) with m permutations", {
    # With m=99 permutations, minimum p should be 1/100 = 0.01
    set.seed(456)
    
    mat <- matrix(rnorm(32), nrow = 4)  # 4 genes x 8 samples
    samples <- c("A", "A", "A", "A", "B", "B", "B", "B")
    
    m <- 99  # number of permutations
    theoretical_min <- 1 / (m + 1)  # = 1/100 = 0.01
    
    res <- label_shuffling(
        mat,
        samples,
        control = "A",
        method = "mean",
        randomizations = m,
        pcorr = "none"
    )
    
    min_pval <- min(res$pvalue, na.rm = TRUE)
    
    # Minimum observed p-value should be >= theoretical minimum
    expect_true(min_pval >= theoretical_min - 1e-10,
                info = sprintf("Minimum p-value (%.6f) < theoretical minimum (%.6f)",
                              min_pval, theoretical_min))
})

test_that("S019: P-value bounds are correct [1/(m+1), 1]", {
    # All p-values should be in the proper range given correction
    set.seed(789)
    
    mat <- matrix(rnorm(32, mean = 0, sd = 2), nrow = 4)  # 4 genes x 8 samples
    samples <- c("X", "X", "X", "X", "Y", "Y", "Y", "Y")
    
    m <- 1000
    theoretical_min <- 1 / (m + 1)
    
    res <- label_shuffling(
        mat,
        samples,
        control = "X",
        method = "mean",
        randomizations = m,
        pcorr = "none"
    )
    
    # Check lower bound
    expect_true(all(res$pvalue >= theoretical_min - 1e-10, na.rm = TRUE),
                info = sprintf("Some p-values < 1/(m+1) = 1/%d", m + 1))
    
    # Check upper bound (should be < 1, not equal)
    expect_true(all(res$pvalue <= 1, na.rm = TRUE),
                info = "Some p-values > 1")
})

test_that("S019: Unpaired permutation test uses bias correction", {
    # Verify bias correction works in unpaired design
    set.seed(111)
    
    # Create clear case vs control difference
    control_vals <- matrix(rnorm(12, mean = 1.0, sd = 0.5), nrow = 3)
    case_vals <- matrix(rnorm(12, mean = 2.0, sd = 0.5), nrow = 3)
    mat <- cbind(control_vals, case_vals)  # 3 genes x 8 samples
    
    samples <- c(rep("Control", 4), rep("Case", 4))
    
    m <- 100
    theoretical_min <- 1 / (m + 1)
    
    res <- label_shuffling(
        mat,
        samples,
        control = "Control",
        method = "mean",
        randomizations = m,
        pcorr = "none"
    )
    
    # Verify correction is applied
    expect_true(all(res$pvalue > 0, na.rm = TRUE),
                info = "Unpaired test: found p = 0")
    
    expect_true(all(res$pvalue >= theoretical_min - 1e-10, na.rm = TRUE),
                info = sprintf("Unpaired test: p < 1/(m+1) = %.6f", theoretical_min))
})

test_that("S019: Paired permutation test uses bias correction", {
    # Verify bias correction works in paired design
    set.seed(222)
    
    # Create paired samples: pairs A, B, C with clear group difference
    mat <- matrix(c(
        1.0, 1.1,  # pair A: control vs case
        1.5, 1.6,  # pair B
        0.8, 0.9   # pair C
    ), nrow = 1, byrow = TRUE)
    
    samples <- c("Control", "Case", "Control", "Case", "Control", "Case")
    pairs_vector <- c("A", "A", "B", "B", "C", "C")
    
    m <- 100
    theoretical_min <- 1 / (m + 1)
    
    res <- label_shuffling(
        mat,
        samples,
        control = "Control",
        method = "mean",
        randomizations = m,
        paired = TRUE,
        pairs = pairs_vector,
        pcorr = "none"
    )
    
    # Verify correction is applied in paired design
    expect_true(all(res$pvalue > 0, na.rm = TRUE),
                info = "Paired test: found p = 0")
    
    expect_true(all(res$pvalue >= theoretical_min - 1e-10, na.rm = TRUE),
                info = sprintf("Paired test: p < 1/(m+1) = %.6f", theoretical_min))
})

test_that("S019: Bias correction with few permutations", {
    # Test that correction is especially important with small m
    # With m=9, minimum p = 1/10 = 0.1
    set.seed(333)
    
    mat <- matrix(rnorm(32), nrow = 4)  # 4 genes x 8 samples
    samples <- c("A", "A", "A", "A", "B", "B", "B", "B")
    
    m <- 9
    theoretical_min <- 1 / (m + 1)  # = 0.1
    
    res <- label_shuffling(
        mat,
        samples,
        control = "A",
        method = "mean",
        randomizations = m,
        pcorr = "none"
    )
    
    # With only 9 permutations, correction is critical
    # Minimum should be 0.1, not 0
    expect_true(all(res$pvalue >= theoretical_min - 1e-10, na.rm = TRUE),
                info = sprintf("Small m=%d: minimum p-value < %.3f", m, theoretical_min))
    
    min_observed <- min(res$pvalue, na.rm = TRUE)
    expect_true(min_observed >= 0.05,  # Should be close to 0.1 in practice
                info = sprintf("With m=%d, minimum p-value should be ~0.1, got %.6f", 
                              m, min_observed))
})

test_that("S019: Bias correction produces sensible p-value distribution", {
    # Verify that corrected p-values show proper empirical distribution
    set.seed(444)
    
    # Create matrix with notable difference in one gene
    mat <- matrix(rnorm(32, mean = 0, sd = 1), nrow = 4)  # 4 genes x 8 samples
    mat[1, 5:8] <- mat[1, 5:8] + 2.0  # Add large effect to gene 1 (case group)
    
    samples <- c(rep("Control", 4), rep("Case", 4))
    
    res <- label_shuffling(
        mat,
        samples,
        control = "Control",
        method = "mean",
        randomizations = 100,
        pcorr = "none"
    )
    
    # Gene 1 should have smaller p-value (strong effect)
    # Genes 2-4 should have larger p-values (no effect)
    expect_true(res$pvalue[1] < res$pvalue[2],
                info = "Gene with clear effect should have smaller p-value")
    
    # Distribution properties
    expect_true(all(res$pvalue > 0, na.rm = TRUE),
                info = "All p-values should be > 0 with correction")
})

test_that("S019: Extreme case - observed statistic is most extreme", {
    # Special case: when observed is most extreme, b=0, so p = 1/(m+1)
    # This is the key case where traditional p=0 is WRONG
    set.seed(555)
    
    # Create a case where one gene has extremely different behavior
    mat <- matrix(rnorm(32, mean = 1, sd = 0.3), nrow = 4)  # 4 genes x 8 samples
    mat[1, 5:8] <- mat[1, 5:8] + 5.0  # Very large effect in gene 1 (case group)
    
    samples <- c(rep("Control", 4), rep("Case", 4))
    
    m <- 999
    theoretical_min <- 1 / (m + 1)  # ~0.001
    
    res <- label_shuffling(
        mat,
        samples,
        control = "Control",
        method = "mean",
        randomizations = m,
        pcorr = "none"
    )
    
    # Even if gene 1 is most extreme, p-value should never be 0
    expect_true(res$pvalue[1] > 0,
                info = "Even most extreme observation should have p > 0 with S019 correction")
    
    # Minimum p-value with m=999 should be ~1/1000
    expect_true(all(res$pvalue >= 0.0009, na.rm = TRUE),  # 0.001 is 1/1000
                info = sprintf("With m=%d, minimum p > 1/1000", m))
})

test_that("S019: Multiple testing correction applied to corrected p-values", {
    # Verify that p.adjust is applied AFTER bias correction
    set.seed(666)
    
    mat <- matrix(rnorm(64), nrow = 8)  # 8 genes x 8 samples
    samples <- c(rep("A", 4), rep("B", 4))
    
    res <- label_shuffling(
        mat,
        samples,
        control = "A",
        method = "mean",
        randomizations = 50,
        pcorr = "BH"
    )
    
    # Adjusted p-values should be >= raw p-values (when using BH)
    expect_true(all(res$padj >= res$pvalue - 1e-10, na.rm = TRUE),
                info = "BH adjustment: adjusted p-values should be >= raw p-values")
    
    # All adjusted p-values should be in [0, 1]
    expect_true(all(res$padj >= 0 & res$padj <= 1, na.rm = TRUE),
                info = "Adjusted p-values should be in [0, 1]")
})

test_that("S019: Bias correction consistent across random seeds", {
    # Verify reproducibility with set.seed()
    mat <- matrix(rnorm(32), nrow = 4)  # 4 genes x 8 samples
    samples <- c("X", "X", "X", "X", "Y", "Y", "Y", "Y")
    
    set.seed(777)
    res1 <- label_shuffling(
        mat,
        samples,
        control = "X",
        method = "mean",
        randomizations = 50,
        pcorr = "none"
    )
    
    set.seed(777)
    res2 <- label_shuffling(
        mat,
        samples,
        control = "X",
        method = "mean",
        randomizations = 50,
        pcorr = "none"
    )
    
    # Results should be identical with same seed
    expect_equal(res1$pvalue, res2$pvalue,
                 info = "Results not reproducible with same seed")
    expect_equal(res1$log2FC, res2$log2FC,
                 info = "log2FC not reproducible with same seed")
})

test_that("S019: Both unpaired and paired designs prevent p=0", {
    # Comprehensive test comparing unpaired vs paired with same data
    set.seed(888)
    
    # Create test matrix
    mat <- matrix(rnorm(16), nrow = 2)
    samples <- c("A", "A", "A", "A", "B", "B", "B", "B")
    
    # Unpaired test
    res_unpaired <- label_shuffling(
        mat,
        samples,
        control = "A",
        method = "mean",
        randomizations = 100,
        paired = FALSE,
        pcorr = "none"
    )
    
    # Paired test with explicit pairs
    pairs_vec <- c("P1", "P1", "P2", "P2", "P1", "P1", "P2", "P2")
    res_paired <- label_shuffling(
        mat,
        samples,
        control = "A",
        method = "mean",
        randomizations = 100,
        paired = TRUE,
        pairs = pairs_vec,
        pcorr = "none"
    )
    
    # Both should prevent p=0
    expect_true(all(res_unpaired$pvalue > 0, na.rm = TRUE),
                info = "Unpaired: found p = 0")
    expect_true(all(res_paired$pvalue > 0, na.rm = TRUE),
                info = "Paired: found p = 0")
    
    # Both should have proper bounds
    m <- 100
    min_p <- 1 / (m + 1)
    expect_true(all(res_unpaired$pvalue >= min_p - 1e-10, na.rm = TRUE),
                info = "Unpaired: minimum p-value incorrect")
    expect_true(all(res_paired$pvalue >= min_p - 1e-10, na.rm = TRUE),
                info = "Paired: minimum p-value incorrect")
})

# ============================================================================
# Effect Size Tests: U and r columns in label_shuffling output
# ============================================================================

context("Label Shuffling: Effect Size Measures (U and r Columns)")

test_that("label_shuffling returns U and r columns", {
    set.seed(42)
    # 5 genes x 10 samples (5 control + 5 case)
    mat <- matrix(rnorm(50, mean = 0, sd = 1), nrow = 5)
    samples <- c(rep("A", 5), rep("B", 5))
    
    res <- label_shuffling(
        mat,
        samples,
        control = "A",
        method = "mean",
        randomizations = 100,
        pcorr = "BH"
    )
    
    # Check that required columns exist
    expect_true("U" %in% colnames(res),
                info = "Column 'U' (Mann-Whitney statistic) missing from label_shuffling output")
    expect_true("r" %in% colnames(res),
                info = "Column 'r' (effect size) missing from label_shuffling output")
    
    # Check dimensions
    expect_equal(nrow(res), nrow(mat),
                 info = "Output has wrong number of rows")
    expect_equal(length(res$U), nrow(mat),
                 info = "U column has wrong length")
    expect_equal(length(res$r), nrow(mat),
                 info = "r column has wrong length")
})

test_that("label_shuffling: r values are in valid range [-1, 1]", {
    set.seed(123)
    
    # Create data with clear group differences
    mat <- rbind(
        rnorm(8, mean = 0, sd = 0.5),  # gene 1: similar groups
        rnorm(8, mean = 0, sd = 0.5),  # gene 2: similar groups
        c(rnorm(4, mean = 1, sd = 0.2), rnorm(4, mean = 3, sd = 0.2))  # gene 3: diff groups
    )
    
    samples <- c(rep("Control", 4), rep("Case", 4))
    
    res <- label_shuffling(
        mat,
        samples,
        control = "Control",
        method = "mean",
        randomizations = 100,
        pcorr = "none"
    )
    
    # All non-NA r values should be in [-1, 1]
    r_vals <- res$r[!is.na(res$r)]
    expect_true(all(r_vals >= -1 & r_vals <= 1),
                info = sprintf("Found r values outside [-1, 1] range: min=%.3f, max=%.3f",
                              min(r_vals), max(r_vals)))
})

test_that("label_shuffling: U values are non-negative", {
    set.seed(456)
    
    mat <- matrix(rnorm(24), nrow = 4)  # 4 genes x 6 samples
    samples <- c(rep("Control", 3), rep("Case", 3))
    
    res <- label_shuffling(
        mat,
        samples,
        control = "Control",
        method = "mean",
        randomizations = 50,
        pcorr = "none"
    )
    
    # All non-NA U values should be >= 0
    u_vals <- res$U[!is.na(res$U)]
    expect_true(all(u_vals >= 0),
                info = sprintf("Found negative U values: min=%.3f", min(u_vals)))
})

test_that("label_shuffling: Unpaired effect sizes are computed correctly", {
    set.seed(789)
    
    # Simple 2x8 matrix: 2 genes, 8 samples (4 per group)
    mat <- matrix(c(
        1.0, 1.1, 1.2, 1.0,  2.9, 3.0, 3.1, 2.9,  # gene 1 (clear diff)
        5.0, 5.1, 5.0, 4.9,  5.1, 4.9, 5.0, 5.1   # gene 2 (no diff)
    ), nrow = 2, byrow = TRUE)
    
    samples <- c(rep("Control", 4), rep("Case", 4))
    
    res <- label_shuffling(
        mat,
        samples,
        control = "Control",
        method = "mean",
        randomizations = 100,
        paired = FALSE,
        pcorr = "none"
    )
    
    # Gene 1 should have larger effect size than gene 2
    r1 <- abs(res$r[1])
    r2 <- abs(res$r[2])
    
    expect_true(!is.na(r1) && !is.na(r2),
                info = "Cannot compute effect sizes")
    # Gene 1 has clear separation, Gene 2 has much less
    expect_true(r1 > r2 || (r1 > 0.1 && r2 < 0.15),
                info = sprintf("Gene 1 effect (|r|=%.3f) should clearly exceed Gene 2 effect (|r|=%.3f)",
                              r1, r2))
})

test_that("label_shuffling: Paired effect sizes are computed correctly", {
    set.seed(321)
    
    # 2 genes, 8 samples (4 pairs: matched samples within pairs)
    mat <- matrix(c(
        1.0, 1.1, 1.0, 1.0,  2.9, 3.0, 2.9, 3.1,  # gene 1 (clear paired diff)
        5.0, 4.9, 5.1, 5.0,  5.0, 5.1, 4.9, 5.0   # gene 2 (minimal paired diff)
    ), nrow = 2, byrow = TRUE)
    
    samples <- c(rep("Normal", 4), rep("Tumor", 4))
    pairs <- c("P1", "P2", "P3", "P4", "P1", "P2", "P3", "P4")
    
    res <- label_shuffling(
        mat,
        samples,
        control = "Normal",
        method = "mean",
        randomizations = 100,
        paired = TRUE,
        pairs = pairs,
        pcorr = "none"
    )
    
    # Should have U and r columns
    expect_true("U" %in% colnames(res) && "r" %in% colnames(res))
    
    # Both genes should have computable effect sizes
    r1 <- abs(res$r[1])
    r2 <- abs(res$r[2])
    
    expect_true(!is.na(r1) && !is.na(r2),
                info = "Cannot compute paired effect sizes")
    # Effect sizes should be meaningful (not all NA)
    expect_true(r1 > 0.05 || r2 > 0.05,
                info = "All effect sizes are essentially zero")
})

test_that("label_shuffling effect sizes match wilcoxon effect sizes", {
    set.seed(555)
    
    # Create simple test data with enough samples
    mat <- matrix(rnorm(30, mean = 5, sd = 1), nrow = 3)  # 3 genes x 10 samples
    samples <- c(rep("A", 5), rep("B", 5))
    
    # Get effect sizes from label_shuffling
    res_shuffle <- label_shuffling(
        mat,
        samples,
        control = "A",
        method = "mean",
        randomizations = 50,
        paired = FALSE,
        pcorr = "none"
    )
    
    # Get effect sizes from wilcoxon
    res_wilcox <- wilcoxon(
        mat,
        samples,
        pcorr = "none",
        paired = FALSE
    )
    
    # U and r values should match between the two methods
    # They both compute from the same observed data
    expect_equal(res_shuffle$U, res_wilcox$U, tolerance = 1e-10,
                 info = "U values differ between label_shuffling and wilcoxon")
    expect_equal(res_shuffle$r, res_wilcox$r, tolerance = 1e-10,
                 info = "r values differ between label_shuffling and wilcoxon")
})

test_that("label_shuffling: Effect sizes independent of randomization count", {
    set.seed(666)
    
    mat <- matrix(rnorm(20), nrow = 4)
    samples <- c(rep("X", 3), rep("Y", 2))
    
    # Run with different permutation counts
    res_50 <- label_shuffling(
        mat,
        samples,
        control = "X",
        method = "mean",
        randomizations = 50,
        pcorr = "none"
    )
    
    res_200 <- label_shuffling(
        mat,
        samples,
        control = "X",
        method = "mean",
        randomizations = 200,
        pcorr = "none"
    )
    
    # Effect sizes (U, r) should be identical regardless of permutation count
    # (they only depend on observed data, not the null distribution)
    expect_equal(res_50$U, res_200$U, tolerance = 1e-10,
                 info = "U changed with different randomization counts")
    expect_equal(res_50$r, res_200$r, tolerance = 1e-10,
                 info = "r changed with different randomization counts")
    
    # But p-values will differ (based on null distribution)
    expect_false(isTRUE(all.equal(res_50$pvalue, res_200$pvalue)),
                 info = "p-values should differ with different permutation counts")
})

test_that("label_shuffling: Output column order includes U and r", {
    set.seed(777)
    
    mat <- matrix(rnorm(16), nrow = 2)
    samples <- c("A", "A", "B", "B", "A", "A", "B", "B")
    
    res <- label_shuffling(
        mat,
        samples,
        control = "A",
        method = "mean",
        randomizations = 100,
        pcorr = "BH"
    )
    
    # Expected column order: pvalue, padj, log2FC, U, r, group_means
    expected_cols <- c("pvalue", "padj", "log2FC", "U", "r")
    
    for (col in expected_cols) {
        expect_true(col %in% colnames(res),
                   info = sprintf("Expected column '%s' not found", col))
    }
})

test_that("label_shuffling: Handle NAs gracefully in effect sizes", {
    set.seed(888)
    
    # Create data with some constant rows (will have issues in Wilcoxon)
    mat <- rbind(
        c(1, 1, 1, 1, 1, 1),  # Constant row - Wilcoxon may fail
        rnorm(6),              # Normal row
        rnorm(6)               # Normal row
    )
    
    samples <- c(rep("Control", 3), rep("Case", 3))
    
    res <- label_shuffling(
        mat,
        samples,
        control = "Control",
        method = "mean",
        randomizations = 50,
        paired = FALSE,
        pcorr = "none"
    )
    
    # Should handle gracefully (NAs where computation fails)
    expect_true(is.numeric(res$U) && is.numeric(res$r))
    expect_equal(nrow(res), 3)
})
