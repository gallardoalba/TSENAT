context("Permutation Tests: Westfall-Young Stepdown for Multi-q Analysis")

test_that("Basic Westfall-Young functionality with simple synthetic data", {
    # Create simple dataset with clear signal at q=1.0
    set.seed(123)
    n_genes <- 20
    n_samples <- 8  # 4 per group
    n_q <- 3  # q = 0.5, 1.0, 1.5
    q_vals <- c(0.5, 1.0, 1.5)
    
    # Set up column names: alternates (sample, q-value) for each column
    col_names <- rep(paste0("Sample", 1:n_samples), n_q)
    col_q <- rep(q_vals, times = n_samples)
    col_names <- paste0(col_names, "_q=", col_q)
    
    # Create diversity matrix
    diversity_data <- matrix(rnorm(n_genes * n_samples * n_q, mean = 0, sd = 1), 
                             nrow = n_genes)
    colnames(diversity_data) <- col_names
    rownames(diversity_data) <- paste0("Gene", 1:n_genes)
    
    # Add signal for Gene 1 at q=1.0
    signal_cols <- which(grepl("_q=1\\.0$", colnames(diversity_data)) & 
                         grepl("^Sample[1-4]", colnames(diversity_data)))
    if (length(signal_cols) > 0) {
        diversity_data[1, signal_cols] <- 
            diversity_data[1, signal_cols] + 2  # Mean difference of ~2
    }
    
    # Sample labels
    samples <- rep(c("Control", "Treatment"), each = n_samples/2)
    
    # Run Westfall-Young
    result <- label_shuffling_westfall_young(
        x = diversity_data,
        q_values = q_vals,
        samples = samples,
        control = "Control",
        randomizations = 100,  # Smaller for testing
        verbose = FALSE
    )
    
    # Check output structure
    expect_equal(nrow(result), n_genes * n_q)
    expect_gte(ncol(result), 6)
    expect_true(all(c("gene", "q_value", "log2FC", "pvalue_raw", "pvalue_wy", "padj_bh") %in% colnames(result)))
    expect_true(all(is.numeric(result$pvalue_raw)))
    expect_true(all(is.numeric(result$pvalue_wy)))
    expect_true(all(result$pvalue_raw >= 0 & result$pvalue_raw <= 1, na.rm = TRUE))
    expect_true(all(result$pvalue_wy >= 0 & result$pvalue_wy <= 1, na.rm = TRUE))
    
    # Check that WY-adjusted p-values are more conservative than raw
    expect_gte(mean(result$pvalue_wy, na.rm = TRUE), 
               mean(result$pvalue_raw, na.rm = TRUE))
    
    # Gene 1, q=1.0 should have smaller p-values than other tests
    gene1_q1 <- result[result$gene == "Gene1" & result$q_value == 1.0, ]
    other_tests <- result[!(result$gene == "Gene1" & result$q_value == 1.0), ]
    
    expect_lte(gene1_q1$pvalue_wy, median(other_tests$pvalue_wy, na.rm = TRUE))
})

test_that("Westfall-Young p-values are monotonically non-decreasing", {
    set.seed(456)
    n_genes <- 15
    n_samples <- 6
    n_q <- 2
    q_vals_test <- c(1.0, 2.0)
    
    diversity_data <- matrix(rnorm(n_genes * n_samples * n_q), nrow = n_genes)
    col_names <- rep(paste0("S", 1:n_samples), n_q)
    col_q <- rep(q_vals_test, times = n_samples)
    colnames(diversity_data) <- paste0(col_names, "_q=", col_q)
    rownames(diversity_data) <- paste0("G", 1:n_genes)
    
    samples <- rep(c("A", "B"), each = n_samples/2)
    
    result <- label_shuffling_westfall_young(
        x = diversity_data,
        q_values = q_vals_test,
        samples = samples,
        control = "A",
        randomizations = 50,
        verbose = FALSE
    )
    
    # Check monotonicity: as we go down sorted p-values, WY p-values should increase
    result_sorted <- result[order(result$pvalue_raw), ]
    
    wy_diffs <- diff(result_sorted$pvalue_wy)
    min_wy_diff <- min(wy_diffs[!is.na(wy_diffs)])
    
    expect_gte(min_wy_diff, -1e-10)  # Allow small numerical errors
})

test_that("Results are reproducible with same random seed", {
    setup_data <- function() {
        set.seed(789)
        diversity_data <- matrix(rnorm(50 * 8 * 2), nrow = 50)
        col_names <- rep(paste0("Sample", 1:8), 2)
        col_q <- rep(c(0.5, 1.5), times = 8)
        colnames(diversity_data) <- paste0(col_names, "_q=", col_q)
        rownames(diversity_data) <- paste0("Gene", 1:50)
        return(diversity_data)
    }
    
    samples <- rep(c("Ctrl", "Trt"), each = 4)
    
    # Run 1 with seed 999
    set.seed(999)
    result1 <- label_shuffling_westfall_young(
        x = setup_data(),
        q_values = c(0.5, 1.5),
        samples = samples,
        control = "Ctrl",
        randomizations = 50,
        verbose = FALSE
    )
    
    # Run 2 with same seed 999
    set.seed(999)
    result2 <- label_shuffling_westfall_young(
        x = setup_data(),
        q_values = c(0.5, 1.5),
        samples = samples,
        control = "Ctrl",
        randomizations = 50,
        verbose = FALSE
    )
    
    # Results should be identical
    expect_equal(nrow(result1), nrow(result2))
    expect_equal(result1$gene, result2$gene)
    expect_equal(result1$q_value, result2$q_value)
    expect_true(max(abs(result1$pvalue_wy - result2$pvalue_wy), na.rm = TRUE) < 1e-10)
})

test_that("FWER control under null hypothesis", {
    # Generate data with NO signal (null hypothesis)
    set.seed(111)
    n_genes <- 10
    n_samples <- 6
    n_q <- 2
    q_vals_test <- c(1.0, 2.0)
    
    # Pure noise, no signal
    diversity_data <- matrix(rnorm(n_genes * n_samples * n_q, mean = 0, sd = 1), 
                             nrow = n_genes)
    col_names <- rep(paste0("S", 1:n_samples), n_q)
    col_q <- rep(q_vals_test, times = n_samples)
    colnames(diversity_data) <- paste0(col_names, "_q=", col_q)
    rownames(diversity_data) <- paste0("G", 1:n_genes)
    
    samples <- rep(c("A", "B"), each = n_samples/2)
    
    result <- label_shuffling_westfall_young(
        x = diversity_data,
        q_values = q_vals_test,
        samples = samples,
        control = "A",
        randomizations = 200,
        verbose = FALSE
    )
    
    # Under null hypothesis with FWER ≤ 0.05, expect low number of false positives
    # With true FWER=0.05, P(≥1 false positives) = 0.05 for the family
    n_sig_wy <- sum(result$pvalue_wy < 0.05, na.rm = TRUE)
    
    # Basic sanity check: shouldn't have excessive false positives
    expect_lte(n_sig_wy, nrow(result) * 0.2)  # Allow up to 20% false positives under null
})

test_that("Single gene edge case is handled correctly", {
    set.seed(222)
    q_vals_test <- c(0.5, 1.0, 1.5)
    diversity_data <- matrix(rnorm(1 * 8 * 3), nrow = 1)
    col_names <- rep(paste0("S", 1:8), 3)
    col_q <- rep(q_vals_test, times = 8)
    colnames(diversity_data) <- paste0(col_names, "_q=", col_q)
    rownames(diversity_data) <- "Gene1"
    
    samples <- rep(c("A", "B"), each = 4)
    
    result <- label_shuffling_westfall_young(
        x = diversity_data,
        q_values = q_vals_test,
        samples = samples,
        control = "A",
        randomizations = 50,
        verbose = FALSE
    )
    
    expect_equal(nrow(result), 3)
    expect_equal(unique(result$gene), "Gene1")
    expect_equal(sort(unique(result$q_value)), c(0.5, 1.0, 1.5))
    expect_true(all(is.numeric(result$pvalue_wy)))
})

test_that("Single q-value edge case is handled correctly", {
    set.seed(333)
    n_genes <- 20
    n_samples <- 8
    
    diversity_data <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)
    colnames(diversity_data) <- paste0("Sample", 1:n_samples, "_q=1.0")
    rownames(diversity_data) <- paste0("Gene", 1:n_genes)
    
    samples <- rep(c("A", "B"), each = n_samples/2)
    
    result <- label_shuffling_westfall_young(
        x = diversity_data,
        q_values = 1.0,
        samples = samples,
        control = "A",
        randomizations = 50,
        verbose = FALSE
    )
    
    expect_equal(nrow(result), n_genes)
    expect_true(all(result$q_value == 1.0))
    expect_true(all(!is.na(result$pvalue_wy)))
})
