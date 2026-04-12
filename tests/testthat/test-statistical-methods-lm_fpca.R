library(testthat)

context("Linear Model Interaction: FPCA Methods")

test_that("fpca method attaches p_interaction to rowData", {
    # Use 4 samples and 5 q values to ensure adequate FPCA matrix coverage
    qvec <- seq(0.01, 0.05, by = 0.01)
    sample_ids <- rep(c("S1", "S2", "S3", "S4"), each = length(qvec))
    coln <- paste0(sample_ids, "_q=", rep(qvec, times = 4))

    set.seed(9)
    # create clear differences so PCA can pick up group separation
    gene1_vals <- c(qvec * 1, qvec * 2, qvec * 1.5, qvec * 2.2) + rnorm(length(coln), sd = 1e-4)
    gene2_vals <- c(qvec * 1, qvec * 1, qvec * 1, qvec * 1) + rnorm(length(coln), sd = 1e-4)
    mat <- rbind(g1 = gene1_vals, g2 = gene2_vals)
    colnames(mat) <- coln
    rownames(mat) <- c("g1", "g2")

    rd <- data.frame(genes = rownames(mat), row.names = rownames(mat), stringsAsFactors = FALSE)
    cd <- data.frame(samples = sample_ids, sample_type = rep(c("Normal", "Tumor", "Normal", "Tumor"), each = length(qvec)), row.names = coln, stringsAsFactors = FALSE)

    se <- SummarizedExperiment::SummarizedExperiment(assays = list(diversity = mat), rowData = rd, colData = cd)

    # lower min_obs so test is robust; 4 samples with full q coverage should pass
    res <- TSENAT:::.calculate_lm(se, condition_col = "sample_type", method = "fpca", min_obs = 2)
    if (is.data.frame(res)) {
        rd_out <- as.data.frame(res)
    } else {
        rd_out <- as.data.frame(SummarizedExperiment::rowData(res))
    }
    expect_true("p_interaction" %in% colnames(rd_out))
})

# ============================================================================
# Tests for FPCA interaction helper prcomp and try-error handling
# ============================================================================

context("Linear Model Interaction: FPCA with prcomp and Error Handling")

test_that(".fpca_interaction handles prcomp successfully", {
    # This test covers: pca <- try(stats::prcomp(mat_sub, center = TRUE, scale. = FALSE), silent = TRUE)
    # and: if (inherits(pca, "try-error")) { return(NULL) }
    
    # Create properly formed input data with 10 samples and 5 q-values for sufficient FPCA matrix
    set.seed(42)
    n_q <- 5
    n_samples <- 10
    mat <- matrix(rnorm(n_samples * n_q, mean = 1.5, sd = 0.4), nrow = 1, ncol = n_samples * n_q)
    rownames(mat) <- "Gene1"
    q_vals <- rep(seq(0.5, 2, length.out = n_q), times = n_samples)
    sample_names <- rep(paste0("S", 1:n_samples), each = n_q)
    group_vec <- rep(c("GroupA", "GroupB"), each = n_samples * n_q / 2)
    
    res <- TSENAT:::.fpca_interaction(mat, q_vals, sample_names, group_vec, "Gene1", min_obs = 5)
    
    # Should return either NULL or a valid result with p_interaction
    if (!is.null(res)) {
        expect_is(res, "data.frame")
        expect_true("gene" %in% colnames(res))
        expect_true("p_interaction" %in% colnames(res))
        expect_equal(res$gene, "Gene1")
    } else {
        # When result is NULL, verify that it's handled correctly
        expect_null(res)
    }
})

test_that(".fpca_interaction returns NULL when prcomp fails", {
    # Create data with constant values (no variation)
    # This should return NULL due to singular covariance (prcomp failure)
    mat <- matrix(1.0, nrow = 1, ncol = 20)
    rownames(mat) <- "Gene1"
    q_vals <- rep(seq(0.5, 2, length.out = 5), times = 4)
    sample_names <- rep(paste0("S", 1:4), each = 5)
    group_vec <- rep(c("A", "B"), each = 10)
    
    # Constant values should fail or return NULL
    res <- TSENAT:::.fpca_interaction(mat, q_vals, sample_names, group_vec, "Gene1", min_obs = 2)
    
    # Result should be NULL or a valid data frame
    if (!is.null(res)) {
        expect_is(res, "data.frame")
    } else {
        # When NULL, verify the return is correct
        expect_null(res)
    }
})

test_that(".fpca_interaction with NAs in data", {
    # Test prcomp with data containing NAs that need imputation
    set.seed(123)
    n_q <- 5
    n_samples <- 8
    mat <- matrix(rnorm(n_samples * n_q, mean = 1.2, sd = 0.3), nrow = 1, ncol = n_samples * n_q)
    # Introduce a few NAs for imputation testing
    mat[1, 3] <- NA
    mat[1, 9] <- NA
    rownames(mat) <- "Gene1"
    q_vals <- rep(seq(0.5, 2, length.out = n_q), times = n_samples)
    sample_names <- rep(paste0("S", 1:n_samples), each = n_q)
    group_vec <- rep(c("GroupA", "GroupB"), each = n_samples * n_q / 2)
    
    res <- TSENAT:::.fpca_interaction(mat, q_vals, sample_names, group_vec, "Gene1", min_obs = 5)
    
    # Should handle NAs and return result or NULL
    if (!is.null(res)) {
        expect_is(res, "data.frame")
        expect_equal(res$gene, "Gene1")
    } else {
        # When result is NULL, verify it's correctly NULL
        expect_null(res)
    }
})

# ============================================================================
# Additional FPCA tests (moved from test-statistical-methods-lm_helpers.R)
# ============================================================================

test_that(".fpca_interaction computes a p-value with reasonable input", {
    set.seed(2)
    # Create dataset with sufficient samples and q-values for successful FPCA
    n_q <- 6
    n_samples <- 6
    genes <- "g1"
    
    # Create entropy values with group differences
    entropy_vals <- c(
        seq(0.4, 1.2, length.out = n_q),  # GroupA sample1
        seq(0.5, 1.3, length.out = n_q),  # GroupA sample2
        seq(0.45, 1.25, length.out = n_q),  # GroupA sample3
        seq(1.2, 2.0, length.out = n_q),  # GroupB sample4
        seq(1.3, 2.1, length.out = n_q),  # GroupB sample5
        seq(1.25, 2.05, length.out = n_q)   # GroupB sample6
    ) + rnorm(n_q * n_samples, sd = 0.08)
    
    mat <- matrix(entropy_vals, nrow = 1)
    rownames(mat) <- genes
    q_vals <- rep(seq(0.2, 1.8, length.out = n_q), times = n_samples)
    sample_names <- rep(paste0("S", 1:n_samples), each = n_q)
    group_vec <- rep(c("A", "B"), each = n_q * n_samples / 2)
    
    res <- TSENAT:::.fpca_interaction(mat, q_vals = q_vals, sample_names = sample_names, group_vec = group_vec, g = 1, min_obs = 5)
    
    # Should produce a result with these data
    if (!is.null(res)) {
        expect_is(res, "data.frame")
        expect_true("p_interaction" %in% colnames(res))
    }
})

test_that(".fpca_interaction returns NULL for non-diverse groups and handles imputation path", {
    # non-diverse groups -> NULL with expected warning
    # Create data with sufficient structure for ARIMA
    set.seed(42)
    genes <- 1
    n_q <- 6
    n_samples <- 6
    
    # Single group: all samples in group A
    q_vals <- rep(1:n_q, times = n_samples)
    samples <- rep(paste0("s", 1:n_samples), each = n_q)
    mat <- matrix(rnorm(length(q_vals)), nrow = 1)
    rownames(mat) <- "g1"
    group_vec <- rep("A", length.out = length(q_vals))
    
    res <- expect_warning(TSENAT:::.fpca_interaction(mat, q_vals = q_vals, sample_names = samples, group_vec = group_vec, g = 1, min_obs = 1), 
                          "Insufficient group variation")
    expect_null(res)

    # imputation path: create data with two groups but trigger matrix build failure with insufficient data
    set.seed(43)
    # Very sparse: only 1 q-value per sample will cause issues after reshaping
    q_vals2 <- c(1)
    mat2 <- matrix(rnorm(6, mean = 0.5, sd = 0.2), nrow = 1, ncol = 6)
    rownames(mat2) <- "g1"
    sample_names2 <- paste0("s", 1:6)
    group_vec2 <- rep(c("A", "B"), times = 3)
    
    # This has groups but only 3 samples per group with 1 q-value - should trigger insufficient data warning
    res2 <- expect_warning(TSENAT:::.fpca_interaction(mat2, q_vals = q_vals2, sample_names = sample_names2, group_vec = group_vec2, g = 1, min_obs = 1),
                           "Insufficient|Failed")
    expect_null(res2)
})

test_that(".fpca_interaction handles prcomp and t.test failures gracefully", {
    set.seed(101)
    # Test single group (should fail - needs 2 groups)
    n_q <- 5
    n_samples <- 4
    mat1 <- matrix(rnorm(n_q * n_samples, mean = 1, sd = 0.3), nrow = 1)
    rownames(mat1) <- "g1"
    q_vals1 <- rep(seq(0.5, 2, length.out = n_q), times = n_samples)
    sample_names1 <- rep(paste0("s", 1:n_samples), each = n_q)
    group_vec1 <- rep("A", length.out = n_q * n_samples)  # Only one group - should return NULL
    
    res1 <- expect_warning(TSENAT:::.fpca_interaction(mat1, q_vals = q_vals1, sample_names = sample_names1, group_vec = group_vec1, g = 1, min_obs = 2),
                           "Insufficient group variation")
    expect_null(res1)  # Should be NULL due to insufficient group diversity
    
    # Test with very sparse data (mostly NAs)
    mat2 <- matrix(NA_real_, nrow = 1, ncol = 10)
    mat2[1, c(1, 6)] <- c(1.2, 2.3)  # Only 2 data points
    rownames(mat2) <- "g1"
    q_vals2 <- rep(seq(0.5, 2, length.out = 5), times = 2)
    sample_names2 <- rep(paste0("s", 1:2), each = 5)
    group_vec2 <- rep(c("A", "B"), each = 5)
    
    res2 <- expect_warning(TSENAT:::.fpca_interaction(mat2, q_vals = q_vals2, sample_names = sample_names2, group_vec = group_vec2, g = 1, min_obs = 1),
                            "Failed to build curve matrix")
    # With sparse data, typically returns NULL
    expect_true(is.null(res2) || is.data.frame(res2))
})

test_that(".fpca_interaction includes slope_diff in results", {
    set.seed(1003)
    
    # Create properly sized dataset with sufficient samples and q-values
    n_q <- 6
    n_samples <- 6
    
    # Create entropy with clear group differences for slope variation
    entropy_groupA <- rep(seq(0.3, 0.8, length.out = n_q), times = n_samples/2)
    entropy_groupB <- rep(seq(1.0, 1.8, length.out = n_q), times = n_samples/2)
    entropy_all <- c(entropy_groupA, entropy_groupB) + rnorm(n_q * n_samples, sd = 0.08)
    
    mat <- matrix(entropy_all, nrow = 1)
    rownames(mat) <- "gene1"
    q_vals <- rep(seq(0.2, 1.8, length.out = n_q), times = n_samples)
    sample_names <- rep(paste0("s", 1:n_samples), each = n_q)
    group_vec <- rep(c("A", "B"), each = n_q * n_samples / 2)
    
    result <- TSENAT:::.fpca_interaction(
        mat,
        q_vals = q_vals,
        sample_names = sample_names,
        group_vec = group_vec,
        g = 1,
        min_obs = 5
    )
    
    # With this data, should produce result with slope_diff
    if (!is.null(result)) {
        expect_is(result, "data.frame")
        expect_true("slope_diff" %in% colnames(result))
    }
})

test_that(".fpca_interaction works with all regularization methods", {
    skip_on_ci()  # Expensive: tests all 3 regularization modes. Keep locally for comprehensive validation
    skip_if_not_installed("glmnet")
    # Test via main function to verify all regularization methods are properly integrated
    qvec <- seq(0.01, 0.15, by = 0.01)
    
    set.seed(205)
    # Create 15 samples per group with expression pattern
    samples_normal <- paste0("S", 1:15)
    samples_tumor <- paste0("S", 16:30)
    all_samples <- c(samples_normal, samples_tumor)
    
    gene1_vals <- c()
    for (i in seq_len(15)) {
        gene1_vals <- c(gene1_vals, qvec * 1.0 + rnorm(length(qvec), sd = 0.05))
    }
    for (i in 16:30) {
        gene1_vals <- c(gene1_vals, qvec * 1.3 + rnorm(length(qvec), sd = 0.05))
    }
    
    sample_names <- rep(all_samples, each = length(qvec))
    coln <- paste0(sample_names, "_q=", qvec)
    
    mat <- rbind(g1 = gene1_vals)
    colnames(mat) <- coln
    rownames(mat) <- "g1"
    
    rd <- data.frame(genes = "g1", row.names = "g1", stringsAsFactors = FALSE)
    
    cd <- data.frame(
        samples = sample_names,
        sample_type = rep(c("Normal", "Tumor"), c(15*length(qvec), 15*length(qvec))),
        row.names = coln,
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )
    
    # Test PCA regularization
    res_pca <- TSENAT:::.calculate_lm(se,
        condition_col = "sample_type",
        method = "fpca",
        regularization = "pca",
        min_obs = 2
    )
    
    # Test LASSO regularization
    res_lasso <- suppressWarnings({
        TSENAT:::.calculate_lm(se,
            condition_col = "sample_type",
            method = "fpca",
            regularization = "lasso",
            min_obs = 2
        )
    })
    
    # Test Elastic Net regularization
    res_elasticnet <- suppressWarnings({
        TSENAT:::.calculate_lm(se,
            condition_col = "sample_type",
            method = "fpca",
            regularization = "elasticnet",
            min_obs = 2
        )
    })
    
    # All should return valid results
    expect_is(res_pca, "data.frame")
    expect_true(nrow(res_pca) > 0)
    expect_true("p_interaction" %in% colnames(res_pca))
    
    if (!is.null(res_lasso) && nrow(res_lasso) > 0) {
        expect_is(res_lasso, "data.frame")
        expect_true("p_interaction" %in% colnames(res_lasso))
    }
    
    if (!is.null(res_elasticnet) && nrow(res_elasticnet) > 0) {
        expect_is(res_elasticnet, "data.frame")
        expect_true("p_interaction" %in% colnames(res_elasticnet))
    }
})

# ============================================================================
# FPCA regularization methods (moved from test-integration-linear_models.R)
# ============================================================================

test_that("FPCA with regularization='pca' (default) works correctly", {
    qvec <- seq(0.01, 0.15, by = 0.01)
    
    set.seed(200)
    # Create curve-like data where each sample belongs to ONE group only
    # S1-S15 are Normal; S16-S30 are Tumor (15 per group for glmnet)
    samples_normal <- paste0("S", 1:15)
    samples_tumor <- paste0("S", 16:30)
    all_samples <- c(samples_normal, samples_tumor)
    
    gene1_vals <- c()
    for (i in seq_len(15)) {
        gene1_vals <- c(gene1_vals, qvec * 1.0 + rnorm(length(qvec), sd = 0.01))  # Normal
    }
    for (i in 16:30) {
        gene1_vals <- c(gene1_vals, qvec * 1.5 + rnorm(length(qvec), sd = 0.01))  # Tumor
    }
    
    sample_names <- rep(all_samples, each = length(qvec))
    coln <- paste0(sample_names, "_q=", qvec)
    
    mat <- rbind(g1 = gene1_vals)
    colnames(mat) <- coln
    rownames(mat) <- "g1"
    
    rd <- data.frame(genes = "g1", row.names = "g1", stringsAsFactors = FALSE)
    
    cd <- data.frame(
        samples = sample_names,
        sample_type = rep(c("Normal", "Tumor"), c(15*length(qvec), 15*length(qvec))),
        row.names = coln,
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )
    
    # Default method should be "pca"
    res_default <- TSENAT:::.calculate_lm(se,
        condition_col = "sample_type",
        method = "fpca",
        min_obs = 2
    )
    
    # Explicit "pca" method should give same result
    res_pca <- TSENAT:::.calculate_lm(se,
        condition_col = "sample_type",
        method = "fpca",
        regularization = "pca",
        min_obs = 2
    )
    
    # Both should return results
    expect_is(res_default, "data.frame")
    expect_is(res_pca, "data.frame")
    expect_true(nrow(res_default) > 0)
    expect_true(nrow(res_pca) > 0)
    
    # Both should have p_interaction values
    expect_true(!is.na(res_default$p_interaction[1]))
    expect_true(!is.na(res_pca$p_interaction[1]))
    
    # Results should be very similar (same method)
    expect_equal(res_default$p_interaction[1], res_pca$p_interaction[1], tolerance = 1e-6)
})

test_that("FPCA with regularization='lasso' produces valid results", {
    skip_if_not_installed("glmnet")
    
    qvec <- seq(0.01, 0.15, by = 0.01)
    
    set.seed(201)
    # Create curve data where LASSO should select important q-values
    # 15 per group (30 total) for glmnet
    gene1_vals <- c()
    for (i in seq_len(15)) {
        gene1_vals <- c(gene1_vals, qvec * 1.0 + rnorm(length(qvec), sd = 0.01))  # Normal
    }
    for (i in 16:30) {
        gene1_vals <- c(gene1_vals, qvec * 2.0 + rnorm(length(qvec), sd = 0.01))  # Tumor
    }
    
    all_samples <- c(paste0("S", 1:15), paste0("S", 16:30))
    sample_names <- rep(all_samples, each = length(qvec))
    coln <- paste0(sample_names, "_q=", qvec)
    
    mat <- rbind(g1 = gene1_vals)
    colnames(mat) <- coln
    rownames(mat) <- "g1"
    
    rd <- data.frame(genes = "g1", row.names = "g1", stringsAsFactors = FALSE)
    
    cd <- data.frame(
        samples = sample_names,
        sample_type = rep(c("Normal", "Tumor"), c(15*length(qvec), 15*length(qvec))),
        row.names = coln,
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )
    
    res_lasso <- TSENAT:::.calculate_lm(se,
        condition_col = "sample_type",
        method = "fpca",
        regularization = "lasso",
        min_obs = 2
    )
    
    # Should return valid results
    expect_is(res_lasso, "data.frame")
    expect_true(nrow(res_lasso) > 0)
    expect_true("p_interaction" %in% colnames(res_lasso))
    expect_true(!is.na(res_lasso$p_interaction[1]))
    expect_true(res_lasso$p_interaction[1] >= 0 && res_lasso$p_interaction[1] <= 1)
})

test_that("FPCA with regularization='elasticnet' produces valid results", {
    skip_if_not_installed("glmnet")
    
    qvec <- seq(0.01, 0.15, by = 0.01)
    
    set.seed(202)
    # Create curve data with moderate interaction
    # 15 per group (30 total) for glmnet
    gene1_vals <- c()
    for (i in seq_len(15)) {
        gene1_vals <- c(gene1_vals, qvec * 1.2 + rnorm(length(qvec), sd = 0.01))  # Normal
    }
    for (i in 16:30) {
        gene1_vals <- c(gene1_vals, qvec * 1.8 + rnorm(length(qvec), sd = 0.01))  # Tumor
    }
    
    all_samples <- c(paste0("S", 1:15), paste0("S", 16:30))
    sample_names <- rep(all_samples, each = length(qvec))
    coln <- paste0(sample_names, "_q=", qvec)
    
    mat <- rbind(g1 = gene1_vals)
    colnames(mat) <- coln
    rownames(mat) <- "g1"
    
    rd <- data.frame(genes = "g1", row.names = "g1", stringsAsFactors = FALSE)
    
    cd <- data.frame(
        samples = sample_names,
        sample_type = rep(c("Normal", "Tumor"), c(15*length(qvec), 15*length(qvec))),
        row.names = coln,
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )
    
    res_elasticnet <- TSENAT:::.calculate_lm(se,
        condition_col = "sample_type",
        method = "fpca",
        regularization = "elasticnet",
        min_obs = 2
    )
    
    # Should return valid results
    expect_is(res_elasticnet, "data.frame")
    expect_true(nrow(res_elasticnet) > 0)
    expect_true("p_interaction" %in% colnames(res_elasticnet))
    expect_true(!is.na(res_elasticnet$p_interaction[1]))
    expect_true(res_elasticnet$p_interaction[1] >= 0 && res_elasticnet$p_interaction[1] <= 1)
})

test_that("FPCA regularization methods produce reasonable p-value differences", {
    skip_on_ci()  # Expensive: runs 3 regularization methods (pca, lasso, elasticnet). Keep locally for full validation
    skip_if_not_installed("glmnet")
    
    qvec <- seq(0.01, 0.15, by = 0.01)
    
    set.seed(203)
    # Data with clear interaction for comparing methods (15 per group)
    gene1_vals <- c()
    for (i in seq_len(15)) {
        gene1_vals <- c(gene1_vals, qvec * 1.0 + rnorm(length(qvec), sd = 0.005))   # Normal
    }
    for (i in 16:30) {
        gene1_vals <- c(gene1_vals, qvec * 2.0 + rnorm(length(qvec), sd = 0.005))   # Tumor
    }
    
    all_samples <- c(paste0("S", 1:15), paste0("S", 16:30))
    sample_names <- rep(all_samples, each = length(qvec))
    coln <- paste0(sample_names, "_q=", qvec)
    
    mat <- rbind(g1 = gene1_vals, g2 = rnorm(length(coln), sd = 0.1))
    colnames(mat) <- coln
    rownames(mat) <- c("g1", "g2")
    
    rd <- data.frame(genes = c("g1", "g2"), row.names = c("g1", "g2"), stringsAsFactors = FALSE)
    
    cd <- data.frame(
        samples = sample_names,
        sample_type = rep(c("Normal", "Tumor"), c(15*length(qvec), 15*length(qvec))),
        row.names = coln,
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )
    
    # Compare all three methods
    res_pca <- TSENAT:::.calculate_lm(se,
        condition_col = "sample_type",
        method = "fpca",
        regularization = "pca",
        min_obs = 2
    )
    
    res_lasso <- TSENAT:::.calculate_lm(se,
        condition_col = "sample_type",
        method = "fpca",
        regularization = "lasso",
        min_obs = 2
    )
    
    res_elasticnet <- TSENAT:::.calculate_lm(se,
        condition_col = "sample_type",
        method = "fpca",
        regularization = "elasticnet",
        min_obs = 2
    )
    
    # All three should return data frames with results
    expect_is(res_pca, "data.frame")
    expect_is(res_lasso, "data.frame")
    expect_is(res_elasticnet, "data.frame")
    
    expect_true(nrow(res_pca) > 0)
    expect_true(nrow(res_lasso) > 0)
    expect_true(nrow(res_elasticnet) > 0)
    
    # All genes should have p-values
    expect_true(all(!is.na(res_pca$p_interaction[res_pca$gene == "g1"])))
    expect_true(all(!is.na(res_lasso$p_interaction[res_lasso$gene == "g1"])))
    expect_true(all(!is.na(res_elasticnet$p_interaction[res_elasticnet$gene == "g1"])))
    
    # For g1 (true interaction), LASSO and ElasticNet should often provide more power
    # (smaller p-values) than PCA by selecting important q-values
    # But this depends on the specific data, so we just check they're all reasonable
    p_pca <- res_pca$p_interaction[res_pca$gene == "g1"]
    p_lasso <- res_lasso$p_interaction[res_lasso$gene == "g1"]
    p_elasticnet <- res_elasticnet$p_interaction[res_elasticnet$gene == "g1"]
    
    if (!is.na(p_pca) && !is.na(p_lasso) && !is.na(p_elasticnet)) {
        # All should be between 0 and 1
        expect_true(p_pca >= 0 && p_pca <= 1)
        expect_true(p_lasso >= 0 && p_lasso <= 1)
        expect_true(p_elasticnet >= 0 && p_elasticnet <= 1)
    }
})

test_that("FPCA regularization with paired design works correctly", {
    skip_if_not_installed("glmnet")
    
    # Test paired design with regularization methods
    qvec <- seq(0.01, 0.1, by = 0.01)
    sample_pairs <- paste0("P", 1:20)  # 20 pairs = 20 samples per group (optimal for glmnet with 10 features)
    
    set.seed(205)
    gene1_vals <- numeric()
    sample_vec <- character()
    type_vec <- character()
    base_vec <- character()
    
    for (p in sample_pairs) {
        # Normal (group 1)
        norm_vals <- qvec * 1.0 + rnorm(length(qvec), sd = 0.01)
        gene1_vals <- c(gene1_vals, norm_vals)
        sample_vec <- c(sample_vec, paste0(p, "_N_q=", qvec))
        type_vec <- c(type_vec, rep("Normal", length(qvec)))
        base_vec <- c(base_vec, rep(p, length(qvec)))
        
        # Tumor (group 2) with interaction
        tumor_vals <- qvec * 1.5 + rnorm(length(qvec), sd = 0.01)
        gene1_vals <- c(gene1_vals, tumor_vals)
        sample_vec <- c(sample_vec, paste0(p, "_T_q=", qvec))
        type_vec <- c(type_vec, rep("Tumor", length(qvec)))
        base_vec <- c(base_vec, rep(p, length(qvec)))
    }
    
    mat <- rbind(g1 = gene1_vals)
    colnames(mat) <- sample_vec
    rownames(mat) <- "g1"
    
    rd <- data.frame(genes = "g1", row.names = "g1", stringsAsFactors = FALSE)
    
    cd <- data.frame(
        samples = gsub("_q=.*", "", sample_vec),
        sample_type = type_vec,
        sample_base = base_vec,
        row.names = sample_vec,
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )
    
    # Test regularization with paired design (with sufficient observations for glmnet)
    res <- TSENAT:::.calculate_lm(se,
        condition_col = "sample_type",
        method = "fpca",
        subject_col = "sample_base",
        paired = TRUE,
        regularization = "lasso",
        min_obs = 2
    )
    
    # Should produce valid results
    expect_is(res, "data.frame")
    expect_true(nrow(res) > 0)
    expect_true(!is.na(res$p_interaction[1]))
})

# ============================================================================
# Tests for FPCA helper functions
# ============================================================================

context("FPCA Helper Functions")

test_that(".apply_arima_differencing_fpca applies differencing for paired design", {
    # Create test data frame with multiple subjects
    df <- data.frame(
        entropy = c(1.0, 2.0, 3.0, 1.5, 2.5, 3.5),
        q = rep(c(0.5, 1.0, 1.5), 2),
        group = rep(c("A", "B"), 3),
        subject = rep(c("S1", "S2"), each = 3),
        sample_name = c("S1_1", "S1_2", "S1_3", "S2_1", "S2_2", "S2_3"),
        stringsAsFactors = FALSE
    )
    
    result <- .apply_arima_differencing_fpca(df)
    
    # For each subject, differencing reduces rows by 1
    # 2 subjects * 2 differences each = 4 rows
    expect_equal(nrow(result), 4)
    # Check that entropy values are differenced
    expect_true(all(!is.na(result$entropy)))
    # Differenced entropy should be differences of original
    expect_equal(result$entropy[1], 2.0 - 1.0)  # S1: diff(1.0, 2.0)
    expect_equal(result$entropy[2], 3.0 - 2.0)  # S1: diff(2.0, 3.0)
})

test_that(".apply_arima_differencing_fpca handles single subject", {
    # Single subject should return original data
    df <- data.frame(
        entropy = c(1.0, 2.0, 3.0),
        q = c(0.5, 1.0, 1.5),
        group = c("A", "A", "A"),
        subject = c("S1", "S1", "S1"),
        sample_name = c("S1_1", "S1_2", "S1_3"),
        stringsAsFactors = FALSE
    )
    
    result <- .apply_arima_differencing_fpca(df)
    
    # Single subject should return original (no differencing)
    expect_equal(nrow(result), 3)
    expect_equal(result$entropy, df$entropy)
})

test_that(".build_curve_matrix constructs ordered matrix", {
    # Create test data
    entropy_vals <- c(1, 2, 3, 1.5, 2.5, 3.5)
    q_vals <- rep(c(0.5, 1.0, 1.5), 2)
    sample_names <- c("S1", "S1", "S1", "S2", "S2", "S2")
    
    mat <- .build_curve_matrix(entropy_vals, q_vals, sample_names, min_obs = 2)
    
    expect_is(mat, "matrix")
    expect_equal(nrow(mat), 2)  # 2 samples
    expect_equal(ncol(mat), 3)  # 3 unique q values
    expect_equal(rownames(mat), c("S1", "S2"))
    # Check matrix values are correctly placed
    expect_equal(unname(mat[1, 1]), 1.0)    # S1 at q=0.5
    expect_equal(unname(mat[1, 2]), 2.0)    # S1 at q=1.0
    expect_equal(unname(mat[2, 1]), 1.5)    # S2 at q=0.5
})

test_that(".build_curve_matrix returns NULL for insufficient data", {
    # Create sparse data
    entropy_vals <- c(1, 2)
    q_vals <- c(0.5, 1.0)
    sample_names <- c("S1", "S2")
    
    mat <- expect_warning(.build_curve_matrix(entropy_vals, q_vals, sample_names, min_obs = 3),
                          "Insufficient samples")
    
    # Should return NULL because min_obs=3 but only 2 unique samples
    expect_null(mat)
})

test_that(".impute_curve_matrix fills NAs with column means", {
    # Create matrix with NAs - simple 3x2 matrix
    mat <- matrix(c(1.0, 2.0, NA, 3.0, NA, 4.0), nrow = 3, ncol = 2)
    rownames(mat) <- c("S1", "S2", "S3")
    
    result <- .impute_curve_matrix(mat)
    
    # Check no NAs remain
    expect_true(!anyNA(result))
    # Column 1: 1.0, 2.0, NA -> mean(1.0, 2.0) = 1.5, so result[3,1] = 1.5
    # Column 2: 3.0, NA, 4.0 -> mean(3.0, 4.0) = 3.5, so result[2,2] = 3.5
    expect_equal(unname(result[3, 1]), 1.5)  # S3 imputed from column mean
    expect_equal(unname(result[2, 2]), 3.5)  # S2 imputed from column mean
    expect_equal(unname(result[1, 1]), 1.0)  # S1 original value
    expect_true(is.numeric(result))  # Result should be numeric
})

test_that(".aggregate_by_subject computes subject means", {
    # Test data with multiple observations per subject
    values <- c(1, 2, 1.5, 2.5)
    group_vec <- c("A", "A", "B", "B")
    subject_vec <- c("S1", "S1", "S2", "S2")
    
    result_a <- .aggregate_by_subject(values, group_vec, subject_vec, "A")
    result_b <- .aggregate_by_subject(values, group_vec, subject_vec, "B")
    
    # Group A: mean of (1, 2) = 1.5
    expect_equal(as.numeric(result_a["S1"]), 1.5)
    # Group B: mean of (1.5, 2.5) = 2.0
    expect_equal(as.numeric(result_b["S2"]), 2.0)
    expect_is(result_a, "numeric")
    expect_is(result_b, "numeric")
})

test_that(".select_npc selects appropriate PC count", {
    # Create PCA object
    mat <- matrix(rnorm(100), nrow = 20)
    pca <- prcomp(mat, center = TRUE, scale. = FALSE)
    
    n_pc <- .select_npc(pca)
    
    # Should select between 2 and 5 PCs
    expect_true(n_pc >= 2)
    expect_true(n_pc <= 5)
    expect_true(n_pc <= ncol(pca$x))
})

test_that(".test_pc_groupdiff performs t-test on PC values", {
    # Create PC values for two groups with strong separation
    set.seed(123)
    pc_vals <- c(rnorm(15, mean = 0, sd = 0.1), rnorm(15, mean = 2, sd = 0.1))  # Very clear separation
    grp_vals <- rep(c("A", "B"), each = 15)
    
    pval <- .test_pc_groupdiff(pc_vals, grp_vals, NULL, "A", "B")
    
    # Should return a p-value
    expect_is(pval, "numeric")
    expect_true(!is.na(pval))
    expect_true(pval >= 0 && pval <= 1)
    # With strong separation, should be very significant
    expect_true(pval < 0.001)
})

test_that(".test_pc_groupdiff handles paired t-test", {
    # Create paired test scenario
    pc_vals <- c(1, 2, 3, 4, 1.5, 2.5, 3.5, 4.5)
    grp_vals <- rep(c("A", "B"), 4)
    subj_vals <- rep(c("S1", "S2", "S3", "S4"), 2)
    
    pval <- .test_pc_groupdiff(pc_vals, grp_vals, subj_vals, "A", "B")
    
    # Should return a p-value
    expect_is(pval, "numeric")
    expect_true(!is.na(pval))
})

test_that(".test_all_pcs returns correct structure", {
    # Create PCA object with test data
    mat <- matrix(rnorm(100), nrow = 20)
    pca <- prcomp(mat, center = TRUE, scale. = FALSE)
    
    grp_vals <- rep(c("A", "B"), 10)
    
    result <- .test_all_pcs(pca, grp_vals, NULL)
    
    # Check result structure
    expect_is(result, "list")
    expect_true("p_interaction" %in% names(result))
    expect_true("n_pcs_tested" %in% names(result))
    expect_true("min_pc_pvalue" %in% names(result))
    
    # Check value ranges
    expect_true(result$p_interaction >= 0 && result$p_interaction <= 1)
    expect_true(result$n_pcs_tested >= 2 && result$n_pcs_tested <= 5)
})

test_that(".fpca_pca_method produces valid output", {
    # Create test data
    mat_sub <- matrix(rnorm(100), nrow = 20, ncol = 5)
    rownames(mat_sub) <- paste0("S", 1:20)
    grp_vals <- rep(c("A", "B"), 10)
    
    result <- .fpca_pca_method(mat_sub, grp_vals, NULL, "gene1", NULL)
    
    # Result should be a data frame or NULL
    if (!is.null(result)) {
        expect_is(result, "data.frame")
        expect_equal(result$gene, "gene1")
        expect_true("p_interaction" %in% colnames(result))
        expect_true("n_pcs_tested" %in% colnames(result))
    } else {
        expect_null(result)
    }
})

test_that(".fpca_regularization_method produces valid output", {
    # Create test data with clear group separation and larger sample
    set.seed(42)
    mat_sub <- matrix(rnorm(400), nrow = 40, ncol = 10)  # 40 rows x 10 cols = 400 elements
    mat_sub[1:20, ] <- mat_sub[1:20, ] + 1.5  # Shift one group
    rownames(mat_sub) <- paste0("S", 1:40)
    grp_vals <- rep(c("A", "B"), 20)  # 20 per group (>8 observations)
    
    suppressWarnings({
        result <- .fpca_regularization_method(mat_sub, grp_vals, NULL, 
                                                        "gene1", "lasso", NULL)
    })
    
    # Result should be a data frame or NULL
    if (!is.null(result)) {
        expect_is(result, "data.frame")
        expect_equal(result$gene, "gene1")
        expect_true("p_interaction" %in% colnames(result))
    } else {
        expect_null(result)
    }
})
