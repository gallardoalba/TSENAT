context("M-Estimation in calculate_difference")

# ============================================================================
# Test 1: M-estimation method validation
# ============================================================================

test_that("M-estimation method is recognized and validated", {
    genes <- paste0("g", seq_len(5))
    mat <- matrix(rnorm(5 * 6), nrow = 5)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- c("A", "A", "A", "B", "B", "B")
    
    # M-estimation with wilcoxon test should work
    result <- suppressWarnings(.calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "m_estimate",
        test = "wilcoxon"
    ))
    
    expect_true(is.data.frame(result))
    expect_true(nrow(result) == 5)
    expect_true("log2_fold_change" %in% colnames(result))
})

# ============================================================================
# Test 2: M-estimation with shuffle permutation test
# ============================================================================

test_that("M-estimation method works with shuffle permutation test", {
    genes <- paste0("g", seq_len(8))
    mat <- matrix(rnorm(8 * 8, mean = 5, sd = 1), nrow = 8)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("A", "B"), each = 4)
    
    # Test with shuffle
    result <- suppressWarnings(.calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "m_estimate",
        test = "shuffle",
        randomizations = 20
    ))
    
    expect_true(is.data.frame(result))
    expect_true(nrow(result) == 8)
    expect_true(all(c("gene_id", "log2_fold_change", "pvalue", "padj") %in% colnames(result)))
})

# ============================================================================
# Test 3: M-estimation robustness with outliers
# ============================================================================

test_that("M-estimation is robust to outliers", {
    genes <- c("g1", "g2")
    # Group B has an outlier at 100
    # Create 2x8 matrix: g1 has values 3-100, g2 has values 10-100
    mat <- matrix(c(
        3, 4, 5, 6, 10, 11, 12, 100,      # g1
        1, 2, 3, 4, 5, 6, 7, 8            # g2
    ), nrow = 2, byrow = TRUE)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- c("A", "A", "A", "A", "B", "B", "B", "B")
    
    result_mean <- suppressWarnings(.calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "mean",
        test = "wilcoxon"
    ))
    
    result_mest <- suppressWarnings(.calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "m_estimate",
        test = "wilcoxon"
    ))
    
    # M-estimate should be smaller than mean (less influenced by outlier)
    mean_log2fc <- result_mean$log2_fold_change[1]
    mest_log2fc <- result_mest$log2_fold_change[1]
    
    expect_true(abs(mest_log2fc) < abs(mean_log2fc))
})

# ============================================================================
# Test 4: Robust loss and scale function parameters
# ============================================================================

test_that("Robust loss and scale parameters are accepted", {
    genes <- paste0("g", seq_len(5))
    mat <- matrix(rnorm(5 * 8), nrow = 5)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("A", "B"), each = 4)
    
    # Test with huber loss and MAD scale
    result_huber <- suppressWarnings(.calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "m_estimate",
        robust_loss_type = "huber",
        robust_scale_method = "mad",
        test = "wilcoxon"
    ))
    expect_true(is.data.frame(result_huber))
    
    # Test with tukey loss and proposal2 scale
    result_tukey <- suppressWarnings(.calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "m_estimate",
        robust_loss_type = "tukey",
        robust_scale_method = "proposal2",
        test = "wilcoxon"
    ))
    expect_true(is.data.frame(result_tukey))
})

# ============================================================================
# Test 5: Consistency of M-estimation across methods
# ============================================================================

test_that("M-estimation produces consistent log2 fold changes", {
    genes <- paste0("g", seq_len(8))
    mat <- matrix(rnorm(8 * 8), nrow = 8)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("A", "B"), each = 4)
    
    # Same data should produce identical log2 fold changes
    result1 <- suppressWarnings(.calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "m_estimate",
        test = "shuffle",
        randomizations = 20,
        seed = 42
    ))
    
    result2 <- suppressWarnings(.calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "m_estimate",
        test = "shuffle",
        randomizations = 20,
        seed = 42
    ))
    
    # Log2 fold changes should be identical (deterministic)
    expect_equal(result1$log2_fold_change, result2$log2_fold_change)
})

# ============================================================================
# Test 6: M-estimation vs mean/median comparison
# ============================================================================

test_that("All three location estimators produce valid results", {
    genes <- paste0("g", seq_len(6))
    mat <- matrix(rnorm(6 * 8), nrow = 6)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("A", "B"), each = 4)
    
    for (method in c("mean", "median", "m_estimate")) {
        result <- suppressWarnings(.calculate_difference(
            df,
            condition_col = samples,
            control = "A",
            method = method,
            test = "wilcoxon"
        ))
        
        expect_true(is.data.frame(result))
        expect_true(nrow(result) == 6)
        expect_true(all(!is.na(result$log2_fold_change)))
    }
})

# ============================================================================
# Test 7: calculate_fc function with m_estimate
# ============================================================================

test_that("calculate_fc works with m_estimate method", {
    mat <- matrix(rnorm(5 * 6), nrow = 5)
    samples <- c("A", "A", "A", "B", "B", "B")
    
    result <- TSENAT:::.calculate_fc(
        mat,
        samples = samples,
        control = "A",
        method = "m_estimate",
        robust_loss_type = "huber",
        robust_scale_method = "mad"
    )
    
    expect_true(is.data.frame(result))
    expect_true(nrow(result) == 5)
    expect_true(all(c("m_estimate_difference", "log2_fold_change") %in% colnames(result)))
})

# ============================================================================
# Test 8: M-estimation with edge cases
# ============================================================================

test_that("M-estimation handles NA values gracefully", {
    genes <- c("g1", "g2")
    mat <- matrix(c(
        1, 2, 3, 4, 5, 6,         # g1 normal
        1, NA, 3, 4, 5, NA        # g2 with NAs
    ), nrow = 2, byrow = TRUE)
    
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- c("A", "A", "A", "B", "B", "B")
    
    result <- suppressWarnings(.calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "m_estimate",
        test = "wilcoxon"
    ))
    
    expect_true(is.data.frame(result))
    expect_true(nrow(result) == 2)
})

# ============================================================================
# Test 9: M-estimation with paired samples
# ============================================================================

test_that("M-estimation works with paired samples", {
    genes <- paste0("g", seq_len(5))
    mat <- matrix(rnorm(5 * 8), nrow = 5)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- c("A", "A", "A", "A", "B", "B", "B", "B")
    pairs <- c(1, 2, 3, 4, 1, 2, 3, 4)
    
    result <- suppressWarnings(.calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "m_estimate",
        test = "shuffle",
        paired = TRUE,
        pairs = pairs,
        randomizations = 20
    ))
    
    expect_true(is.data.frame(result))
    expect_true(nrow(result) == 5)
})

# ============================================================================
# Test 10: M-estimation effect size and p-value relationships
# ============================================================================

test_that("M-estimation works with multiple genes", {
    genes <- paste0("g", seq_len(10))
    set.seed(123)
    mat <- matrix(rnorm(10 * 8), nrow = 10)
    
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("A", "B"), each = 4)
    
    result <- suppressWarnings(.calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "m_estimate",
        test = "shuffle",
        randomizations = 20
    ))
    
    expect_true(is.data.frame(result))
    expect_true(nrow(result) == 10)
})

context("M-Estimation Integration Tests")

# ============================================================================
# Test 11: SummarizedExperiment support with m_estimate
# ============================================================================

test_that("M-estimation works with SummarizedExperiment input", {
    # Create a simple SummarizedExperiment
    genes <- paste0("g", seq_len(8))
    mat <- matrix(rnorm(8 * 8), nrow = 8)
    dimnames(mat) <- list(genes, paste0("s", 1:8))
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        colData = data.frame(
            sample_type = rep(c("A", "B"), each = 4),
            row.names = colnames(mat)
        )
    )
    
    result <- suppressWarnings(.calculate_difference(
        se,
        condition_col = "sample_type",
        control = "A",
        method = "m_estimate",
        test = "wilcoxon"
    ))
    
    expect_true(is.data.frame(result))
    expect_true(nrow(result) == 8)
})

# ============================================================================
# Test 12: M-estimation reproducibility with multiple seeds
# ============================================================================

test_that("Different seeds produce consistent log2 fold changes", {
    genes <- paste0("g", seq_len(6))
    mat <- matrix(rnorm(6 * 8), nrow = 6)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("A", "B"), each = 4)
    
    result1 <- suppressWarnings(.calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "m_estimate",
        test = "shuffle",
        randomizations = 20,
        seed = 42
    ))
    
    result2 <- suppressWarnings(.calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "m_estimate",
        test = "shuffle",
        randomizations = 20,
        seed = 99
    ))
    
    # Log2 fold changes should be identical (deterministic)
    expect_equal(result1$log2_fold_change, result2$log2_fold_change)
})
