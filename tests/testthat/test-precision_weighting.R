# Tests for precision weighting functionality
# calculate_difference() with use_precision_weights parameter
# .apply_precision_weighting_to_test() helper function

library(TSENAT)

context("precision_weighting: Precision Weighting for Statistical Tests")

# =============================================================================
# Tests for calculate_difference() with precision weighting
# =============================================================================

test_that("calculate_difference with use_precision_weights=FALSE (default) works", {
    # Create sample data
    diversity_mat <- as.data.frame(matrix(runif(80), ncol = 8))
    rownames(diversity_mat) <- paste0("gene_", 1:nrow(diversity_mat))
    colnames(diversity_mat) <- paste0("S", 1:8)
    samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
    
    # Without precision weighting (default)
    result <- calculate_difference(
        x = cbind(genes = rownames(diversity_mat), diversity_mat),
        samples = samples,
        control = "Healthy",
        test = "wilcoxon",
        use_precision_weights = FALSE,
        verbose = FALSE
    )
    
    expect_is(result, "data.frame")
    expect_true("pvalue" %in% colnames(result))
    expect_true("padj" %in% colnames(result))
})

test_that("calculate_difference with precision weighting requires count data", {
    diversity_mat <- matrix(runif(80), ncol = 8)
    samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
    
    # Should error when use_precision_weights=TRUE but counts=NULL
    expect_error(
        calculate_difference(
            x = diversity_mat,
            samples = samples,
            control = "Healthy",
            test = "wilcoxon",
            use_precision_weights = TRUE,
            counts = NULL
        )
    )
})

test_that("calculate_difference with precision weighting requires alpha parameter", {
    diversity_mat <- matrix(runif(80), ncol = 8)
    counts_mat <- matrix(rpois(80, lambda = 20), ncol = 8)
    samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
    
    # Should error when alpha is missing
    expect_error(
        calculate_difference(
            x = diversity_mat,
            samples = samples,
            control = "Healthy",
            test = "wilcoxon",
            use_precision_weights = TRUE,
            counts = counts_mat,
            alpha = NULL,
            beta = 1.5
        )
    )
})

test_that("calculate_difference with precision weighting requires beta parameter", {
    diversity_mat <- matrix(runif(80), ncol = 8)
    counts_mat <- matrix(rpois(80, lambda = 20), ncol = 8)
    samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
    
    # Should error when beta is missing
    expect_error(
        calculate_difference(
            x = diversity_mat,
            samples = samples,
            control = "Healthy",
            test = "wilcoxon",
            use_precision_weights = TRUE,
            counts = counts_mat,
            alpha = 1.5,
            beta = NULL
        )
    )
})

test_that("calculate_difference with precision weighting validates alpha > 0", {
    diversity_mat <- matrix(runif(80), ncol = 8)
    counts_mat <- matrix(rpois(80, lambda = 20), ncol = 8)
    samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
    
    # Should error when alpha <= 0
    expect_error(
        calculate_difference(
            x = diversity_mat,
            samples = samples,
            control = "Healthy",
            test = "wilcoxon",
            use_precision_weights = TRUE,
            counts = counts_mat,
            alpha = 0,
            beta = 1.5
        )
    )
})

test_that("calculate_difference with precision weighting validates beta > 0", {
    diversity_mat <- matrix(runif(80), ncol = 8)
    counts_mat <- matrix(rpois(80, lambda = 20), ncol = 8)
    samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
    
    # Should error when beta <= 0
    expect_error(
        calculate_difference(
            x = diversity_mat,
            samples = samples,
            control = "Healthy",
            test = "wilcoxon",
            use_precision_weights = TRUE,
            counts = counts_mat,
            alpha = 1.5,
            beta = -1
        )
    )
})

test_that("calculate_difference with precision weighting validates counts dimensions", {
    diversity_mat <- matrix(runif(80), ncol = 8)
    counts_mat <- matrix(rpois(96, lambda = 20), ncol = 8)  # Wrong nrow! (12 rows instead of 10)
    samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
    
    # Should error when counts has different number of rows
    expect_error(
        calculate_difference(
            x = diversity_mat,
            samples = samples,
            control = "Healthy",
            test = "wilcoxon",
            use_precision_weights = TRUE,
            counts = counts_mat,
            alpha = 1.5,
            beta = 2.0
        )
    )
})

test_that("calculate_difference with precision weighting returns expected columns", {
    diversity_mat <- as.data.frame(matrix(runif(80), ncol = 8))
    rownames(diversity_mat) <- paste0("gene_", 1:nrow(diversity_mat))
    colnames(diversity_mat) <- paste0("S", 1:8)
    
    counts_mat <- matrix(rpois(80, lambda = 20), ncol = 8, nrow = nrow(diversity_mat))
    rownames(counts_mat) <- rownames(diversity_mat)
    colnames(counts_mat) <- colnames(diversity_mat)
    
    samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
    
    result <- calculate_difference(
        x = cbind(genes = rownames(diversity_mat), diversity_mat),
        samples = samples,
        control = "Healthy",
        test = "wilcoxon",
        use_precision_weights = TRUE,
        counts = counts_mat,
        alpha = 1.5,
        beta = 2.0
    )
    
    expect_is(result, "data.frame")
    expect_true("pvalue_original" %in% colnames(result))
    expect_true("pvalue" %in% colnames(result))
    expect_true("padj" %in% colnames(result))
    expect_true("precision" %in% colnames(result))
    # U statistic is the test statistic from wilcoxon test
    expect_true("U" %in% colnames(result))
})

test_that("calculate_difference precision-weighted results have valid p-values", {
    diversity_mat <- as.data.frame(matrix(runif(80), ncol = 8))
    rownames(diversity_mat) <- paste0("gene_", 1:nrow(diversity_mat))
    colnames(diversity_mat) <- paste0("S", 1:8)
    
    counts_mat <- matrix(rpois(80, lambda = 20), ncol = 8, nrow = nrow(diversity_mat))
    rownames(counts_mat) <- rownames(diversity_mat)
    colnames(counts_mat) <- colnames(diversity_mat)
    
    samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
    
    result <- calculate_difference(
        x = cbind(genes = rownames(diversity_mat), diversity_mat),
        samples = samples,
        control = "Healthy",
        test = "wilcoxon",
        use_precision_weights = TRUE,
        counts = counts_mat,
        alpha = 1.5,
        beta = 2.0
    )
    
    # P-values should be in [0, 1]
    expect_true(all(result$pvalue >= 0 & result$pvalue <= 1))
    expect_true(all(result$padj >= 0 & result$padj <= 1))
    expect_true(all(result$pvalue_original >= 0 & result$pvalue_original <= 1))
})

test_that("calculate_difference precision weights are normalized to [0,1]", {
    diversity_mat <- as.data.frame(matrix(runif(80), ncol = 8))
    rownames(diversity_mat) <- paste0("gene_", 1:nrow(diversity_mat))
    colnames(diversity_mat) <- paste0("S", 1:8)
    
    counts_mat <- matrix(rpois(80, lambda = 20), ncol = 8, nrow = nrow(diversity_mat))
    rownames(counts_mat) <- rownames(diversity_mat)
    colnames(counts_mat) <- colnames(diversity_mat)
    
    samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
    
    result <- calculate_difference(
        x = cbind(genes = rownames(diversity_mat), diversity_mat),
        samples = samples,
        control = "Healthy",
        test = "wilcoxon",
        use_precision_weights = TRUE,
        counts = counts_mat,
        alpha = 1.5,
        beta = 2.0
    )
    
    # Precision weights should be in [0, 1] (normalized)
    expect_true(all(result$precision >= 0 & result$precision <= 1))
    # Maximum precision should be close to 1
    expect_true(max(result$precision) > 0.9)
})

test_that("calculate_difference precision weighting with wilcoxon method", {
    # Create data as data.frame with gene IDs
    diversity_mat <- as.data.frame(matrix(runif(80), ncol = 8))
    rownames(diversity_mat) <- paste0("gene_", 1:nrow(diversity_mat))
    colnames(diversity_mat) <- paste0("S", 1:8)
    
    counts_mat <- matrix(rpois(80, lambda = 20), ncol = 8, nrow = nrow(diversity_mat))
    rownames(counts_mat) <- rownames(diversity_mat)
    colnames(counts_mat) <- colnames(diversity_mat)
    
    samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
    
    result <- calculate_difference(
        x = cbind(genes = rownames(diversity_mat), diversity_mat),
        samples = samples,
        control = "Healthy",
        test = "wilcoxon",
        use_precision_weights = TRUE,
        counts = counts_mat,
        alpha = 1.5,
        beta = 2.0
    )
    
    expect_equal(all(result$method == "wilcoxon"), TRUE)
    expect_equal(nrow(result), nrow(diversity_mat))
})

test_that("calculate_difference precision weighting with shuffle method", {
    # Create data as data.frame with gene IDs
    diversity_mat <- as.data.frame(matrix(runif(80), ncol = 8))
    rownames(diversity_mat) <- paste0("gene_", 1:nrow(diversity_mat))
    colnames(diversity_mat) <- paste0("S", 1:8)
    
    counts_mat <- matrix(rpois(80, lambda = 20), ncol = 8, nrow = nrow(diversity_mat))
    rownames(counts_mat) <- rownames(diversity_mat)
    colnames(counts_mat) <- colnames(diversity_mat)
    
    samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
    
    result <- suppressWarnings(
        calculate_difference(
            x = cbind(genes = rownames(diversity_mat), diversity_mat),
            samples = samples,
            control = "Healthy",
            test = "shuffle",
            use_precision_weights = TRUE,
            counts = counts_mat,
            alpha = 1.5,
            beta = 2.0,
            randomizations = 50  # Reduced for faster testing
        )
    )
    
    expect_equal(all(result$method == "shuffle"), TRUE)
    expect_equal(nrow(result), nrow(diversity_mat))
})

test_that("calculate_difference precision weighting is deterministic with seed", {
    # Create data as data.frame with gene IDs
    diversity_mat <- as.data.frame(matrix(runif(80), ncol = 8))
    rownames(diversity_mat) <- paste0("gene_", 1:nrow(diversity_mat))
    colnames(diversity_mat) <- paste0("S", 1:8)
    
    counts_mat <- matrix(rpois(80, lambda = 20), ncol = 8, nrow = nrow(diversity_mat))
    rownames(counts_mat) <- rownames(diversity_mat)
    colnames(counts_mat) <- colnames(diversity_mat)
    
    samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
    
    set.seed(123)
    result1 <- calculate_difference(
        x = cbind(genes = rownames(diversity_mat), diversity_mat),
        samples = samples,
        control = "Healthy",
        test = "wilcoxon",
        use_precision_weights = TRUE,
        counts = counts_mat,
        alpha = 1.5,
        beta = 2.0,
        seed = 123
    )
    
    set.seed(123)
    result2 <- calculate_difference(
        x = cbind(genes = rownames(diversity_mat), diversity_mat),
        samples = samples,
        control = "Healthy",
        test = "wilcoxon",
        use_precision_weights = TRUE,
        counts = counts_mat,
        alpha = 1.5,
        beta = 2.0,
        seed = 123
    )
    
    expect_equal(result1$pvalue, result2$pvalue)
    expect_equal(result1$precision, result2$precision)
})

test_that("calculate_difference precision weighting affects p-values (comparison test)", {
    # Create data with high variance in counts per gene
    set.seed(456)
    diversity_mat <- as.data.frame(matrix(runif(80), ncol = 8))
    rownames(diversity_mat) <- paste0("gene_", 1:nrow(diversity_mat))
    colnames(diversity_mat) <- paste0("S", 1:8)
    
    counts_mat <- matrix(rpois(80, lambda = 20), ncol = 8, nrow = nrow(diversity_mat))
    rownames(counts_mat) <- rownames(diversity_mat)
    colnames(counts_mat) <- colnames(diversity_mat)
    # Make some genes have very high or very low counts
    counts_mat[1, ] <- rpois(8, lambda = 100)  # High count gene
    counts_mat[2, ] <- rpois(8, lambda = 1)    # Low count gene
    
    samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
    
    # Without precision weighting
    result_no_weight <- calculate_difference(
        x = cbind(genes = rownames(diversity_mat), diversity_mat),
        samples = samples,
        control = "Healthy",
        test = "wilcoxon",
        use_precision_weights = FALSE
    )
    
    # With precision weighting
    result_weight <- calculate_difference(
        x = cbind(genes = rownames(diversity_mat), diversity_mat),
        samples = samples,
        control = "Healthy",
        test = "wilcoxon",
        use_precision_weights = TRUE,
        counts = counts_mat,
        alpha = 1.5,
        beta = 2.0
    )
    
    # Results should differ (precision weighting should change p-values)
    # At least one p-value should be different
    expect_true(any(result_no_weight$pvalue != result_weight$pvalue))
})

# =============================================================================
# Tests for .apply_precision_weighting_to_test() helper function
# =============================================================================

test_that(".apply_precision_weighting_to_test requires proper test_results format", {
    # Mock test_results (minimal required columns)
    test_results <- data.frame(
        pvalue = c(0.01, 0.05, 0.1),
        padj = c(0.03, 0.15, 0.3),
        statistic = c(100, 200, 150),
        method = c("wilcoxon", "wilcoxon", "wilcoxon")
    )
    
    counts <- matrix(rpois(24, lambda = 20), nrow = 3, ncol = 8)
    samples <- c(rep("A", 4), rep("B", 4))
    
    # Function should work with proper inputs
    result <- TSENAT:::.apply_precision_weighting_to_test(
        test_results = test_results,
        counts = counts,
        samples = samples,
        alpha = 1.5,
        beta = 2.0,
        pcorr = "BH"
    )
    
    expect_is(result, "data.frame")
})

test_that(".apply_precision_weighting_to_test returns correct column structure", {
    test_results <- data.frame(
        pvalue = c(0.01, 0.05),
        padj = c(0.03, 0.15),
        statistic = c(100, 200),
        method = c("wilcoxon", "wilcoxon")
    )
    
    counts <- matrix(rpois(16, lambda = 20), nrow = 2, ncol = 8)
    samples <- c(rep("A", 4), rep("B", 4))
    
    result <- TSENAT:::.apply_precision_weighting_to_test(
        test_results = test_results,
        counts = counts,
        samples = samples,
        alpha = 1.5,
        beta = 2.0,
        pcorr = "BH"
    )
    
    expect_equal(ncol(result), 6)  # 6 columns
    expect_true("pvalue_original" %in% colnames(result))
    expect_true("pvalue" %in% colnames(result))
    expect_true("padj" %in% colnames(result))
    expect_true("statistic" %in% colnames(result))
    expect_true("precision" %in% colnames(result))
    expect_true("method" %in% colnames(result))
})

test_that(".apply_precision_weighting_to_test preserves row count", {
    n_genes <- 10
    test_results <- data.frame(
        pvalue = runif(n_genes),
        padj = runif(n_genes),
        statistic = rnorm(n_genes),
        method = rep("wilcoxon", n_genes)
    )
    
    counts <- matrix(rpois(n_genes * 8, lambda = 20), nrow = n_genes, ncol = 8)
    samples <- c(rep("A", 4), rep("B", 4))
    
    result <- TSENAT:::.apply_precision_weighting_to_test(
        test_results = test_results,
        counts = counts,
        samples = samples,
        alpha = 1.5,
        beta = 2.0,
        pcorr = "BH"
    )
    
    expect_equal(nrow(result), n_genes)
})

test_that(".apply_precision_weighting_to_test computes valid precision weights", {
    test_results <- data.frame(
        pvalue = c(0.001, 0.01, 0.05),
        padj = c(0.003, 0.03, 0.15),
        statistic = c(150, 100, 50),
        method = c("wilcoxon", "wilcoxon", "wilcoxon")
    )
    
    counts <- matrix(rpois(24, lambda = 20), nrow = 3, ncol = 8)
    samples <- c(rep("A", 4), rep("B", 4))
    
    result <- TSENAT:::.apply_precision_weighting_to_test(
        test_results = test_results,
        counts = counts,
        samples = samples,
        alpha = 1.5,
        beta = 2.0,
        pcorr = "BH"
    )
    
    # Precision should be in [0, 1]
    expect_true(all(result$precision >= 0 & result$precision <= 1))
    # Max precision should be 1 (normalized)
    expect_equal(max(result$precision), 1, tolerance = 1e-6)
})

test_that(".apply_precision_weighting_to_test preserves original p-values", {
    original_pvalues <- c(0.001, 0.01, 0.05)
    test_results <- data.frame(
        pvalue = original_pvalues,
        padj = c(0.003, 0.03, 0.15),
        statistic = c(150, 100, 50),
        method = c("wilcoxon", "wilcoxon", "wilcoxon")
    )
    
    counts <- matrix(rpois(24, lambda = 20), nrow = 3, ncol = 8)
    samples <- c(rep("A", 4), rep("B", 4))
    
    result <- TSENAT:::.apply_precision_weighting_to_test(
        test_results = test_results,
        counts = counts,
        samples = samples,
        alpha = 1.5,
        beta = 2.0,
        pcorr = "BH"
    )
    
    # Original p-values should be preserved
    expect_equal(result$pvalue_original, original_pvalues)
})

test_that(".apply_precision_weighting_to_test modifies p-values", {
    test_results <- data.frame(
        pvalue = c(0.001, 0.01, 0.05),
        padj = c(0.003, 0.03, 0.15),
        statistic = c(150, 100, 50),
        method = c("wilcoxon", "wilcoxon", "wilcoxon")
    )
    
    counts <- matrix(rpois(24, lambda = 20), nrow = 3, ncol = 8)
    samples <- c(rep("A", 4), rep("B", 4))
    
    result <- TSENAT:::.apply_precision_weighting_to_test(
        test_results = test_results,
        counts = counts,
        samples = samples,
        alpha = 1.5,
        beta = 2.0,
        pcorr = "BH"
    )
    
    # Weighted p-values should be present
    expect_true(all(result$pvalue >= 0 & result$pvalue <= 1))
    # Adjusted p-values should be present
    expect_true(all(result$padj >= 0 & result$padj <= 1))
})

test_that(".apply_precision_weighting_to_test returns valid p-value ranges", {
    test_results <- data.frame(
        pvalue = c(0.001, 0.01, 0.5),
        padj = c(0.003, 0.03, 1.0),
        statistic = c(150, 100, 5),
        method = c("wilcoxon", "wilcoxon", "wilcoxon")
    )
    
    counts <- matrix(rpois(24, lambda = 20), nrow = 3, ncol = 8)
    samples <- c(rep("A", 4), rep("B", 4))
    
    result <- TSENAT:::.apply_precision_weighting_to_test(
        test_results = test_results,
        counts = counts,
        samples = samples,
        alpha = 1.5,
        beta = 2.0,
        pcorr = "BH"
    )
    
    # All p-values should be valid
    expect_true(all(is.finite(result$pvalue)))
    expect_true(all(is.finite(result$padj)))
    expect_true(all(is.finite(result$precision)))
})

# =============================================================================
# Integration tests: Full precision weighting workflow
# =============================================================================

test_that("Full workflow: fit prior, compute posteriors, apply precision weighting", {
    # Simulate realistic RNA-seq data
    set.seed(789)
    n_genes <- 50
    n_samples <- 16
    
    # Create count matrix
    counts_matrix <- matrix(rpois(n_genes * n_samples, lambda = 25), 
                           nrow = n_genes, ncol = n_samples)
    rownames(counts_matrix) <- paste0("gene_", 1:n_genes)
    colnames(counts_matrix) <- paste0("S", 1:n_samples)
    
    # Create diversity measures as data.frame with gene IDs
    diversity_mat <- as.data.frame(matrix(runif(n_genes * n_samples), nrow = n_genes, ncol = n_samples))
    rownames(diversity_mat) <- rownames(counts_matrix)
    colnames(diversity_mat) <- colnames(counts_matrix)
    
    # Sample labels
    samples <- c(rep("Control", 8), rep("Treatment", 8))
    
    # Step 1: Fit empirical Bayes prior
    prior_params <- fit_empirical_beta_prior(counts_matrix)
    
    # Step 2: Run differential analysis WITHOUT precision weighting
    result_no_weight <- calculate_difference(
        x = cbind(genes = rownames(diversity_mat), diversity_mat),
        samples = samples,
        control = "Control",
        test = "wilcoxon",
        use_precision_weights = FALSE
    )
    
    # Step 3: Run differential analysis WITH precision weighting
    result_weight <- calculate_difference(
        x = cbind(genes = rownames(diversity_mat), diversity_mat),
        samples = samples,
        control = "Control",
        test = "wilcoxon",
        use_precision_weights = TRUE,
        counts = counts_matrix,
        alpha = prior_params$alpha,
        beta = prior_params$beta
    )
    
    # Verify structure
    expect_equal(nrow(result_no_weight), nrow(result_weight))
    expect_equal(nrow(result_weight), n_genes)
    
    # Verify precision-weighted result has additional columns
    expect_true("precision" %in% colnames(result_weight))
    expect_true("pvalue_original" %in% colnames(result_weight))
})

test_that("Precision weighting affects genes with variable counts", {
    # Create data where precision weighting should have measurable effect
    set.seed(999)
    n_genes <- 20
    n_samples <- 10
    
    # Diverse gene expression levels
    gene_means <- c(rep(5, 5), rep(50, 5), rep(500, 5), rep(5000, 5))
    counts_matrix <- matrix(0, nrow = n_genes, ncol = n_samples)
    rownames(counts_matrix) <- paste0("gene_", 1:n_genes)
    colnames(counts_matrix) <- paste0("S", 1:n_samples)
    
    for (i in seq_len(n_genes)) {
        counts_matrix[i, ] <- rpois(n_samples, lambda = gene_means[i])
    }
    
    # Create corresponding diversity matrix as data.frame
    diversity_mat <- as.data.frame(matrix(runif(n_genes * n_samples, min = 0.1, max = 0.9), 
                           nrow = n_genes, ncol = n_samples))
    rownames(diversity_mat) <- rownames(counts_matrix)
    colnames(diversity_mat) <- colnames(counts_matrix)
    
    samples <- c(rep("A", 5), rep("B", 5))
    
    # Fit prior
    prior_params <- fit_empirical_beta_prior(counts_matrix)
    
    # Run with precision weighting
    result <- calculate_difference(
        x = cbind(genes = rownames(diversity_mat), diversity_mat),
        samples = samples,
        control = "A",
        test = "wilcoxon",
        use_precision_weights = TRUE,
        counts = counts_matrix,
        alpha = prior_params$alpha,
        beta = prior_params$beta
    )
    
    # Genes with higher expression (higher precision) should be identified
    high_expr_precision <- mean(result$precision[gene_means > 100])
    low_expr_precision <- mean(result$precision[gene_means < 50])
    
    expect_true(high_expr_precision > low_expr_precision)
})
