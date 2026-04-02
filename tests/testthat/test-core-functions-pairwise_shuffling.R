context("Pairwise Label Shuffling Tests")

# Create test data ============================================================

# Smaller test matrix for quick tests
create_test_matrix <- function(nrows = 5, ncols = 6, seed = 42) {
    set.seed(seed)
    matrix(rnorm(nrows * ncols, mean = 5, sd = 2), nrow = nrows, ncol = ncols)
}

# Create test sample labels (3 control, 3 case)
create_test_samples <- function() {
    c("Normal", "Normal", "Normal", "Tumor", "Tumor", "Tumor")
}

# Create test pairing information
create_test_pairs <- function() {
    c("Pair1", "Pair2", "Pair3", "Pair1", "Pair2", "Pair3")
}

# Create test data frame for testing
create_test_data_frame <- function() {
    df <- data.frame(
        feature_id = paste0("Gene", 1:5),
        pvalue = c(0.01, 0.05, 0.1, 0.2, 0.8),
        ctrl_mean = c(5.2, 4.8, 5.1, 5.0, 5.3),
        case_mean = c(7.1, 6.2, 5.5, 5.1, 5.2),
        stringsAsFactors = FALSE
    )
    df$log2FC <- log2(df$case_mean / df$ctrl_mean)
    df
}

# Tests for .validate_label_shuffling_inputs ================================

test_that(".validate_label_shuffling_inputs: accepts valid inputs", {
    samples <- create_test_samples()
    control <- "Normal"
    pairs <- create_test_pairs()
    n_samples <- 6
    
    # Should not raise an error
    expect_silent(.validate_label_shuffling_inputs(samples, control, pairs, n_samples, paired = FALSE))
})

test_that(".validate_label_shuffling_inputs: rejects invalid control group", {
    samples <- create_test_samples()
    control <- "Unknown"
    pairs <- NULL
    n_samples <- 6
    
    expect_error(
        .validate_label_shuffling_inputs(samples, control, pairs, n_samples, paired = FALSE),
        "Control group.*not found"
    )
})

test_that(".validate_label_shuffling_inputs: rejects non-binary sample groups", {
    samples <- c("A", "A", "B", "B", "C", "C")  # 3 groups instead of 2
    control <- "A"
    pairs <- NULL
    n_samples <- 6
    
    expect_error(
        .validate_label_shuffling_inputs(samples, control, pairs, n_samples, paired = FALSE),
        "requires exactly 2 sample groups"
    )
})

test_that(".validate_label_shuffling_inputs: rejects mismatched pairs length", {
    samples <- create_test_samples()
    control <- "Normal"
    pairs <- c("P1", "P2", "P3")  # Length 3, should be 6
    n_samples <- 6
    
    expect_error(
        .validate_label_shuffling_inputs(samples, control, pairs, n_samples, paired = TRUE),
        "pairs.*length equal to n_samples"
    )
})

# Tests for .compute_pseudocount =============================================

test_that(".compute_pseudocount: returns minimum positive value / 2", {
    fc_result <- data.frame(
        ctrl_mean = c(5.0, 3.0, 2.5, 1.0),
        case_mean = c(7.0, 5.0, 4.5, 3.0),
        log2FC = c(0.485, 0.737, 0.848, 1.585),
        pvalue = c(0.01, 0.05, 0.1, 0.2)
    )
    
    pseudocount <- .compute_pseudocount(fc_result)
    
    # Minimum positive value should be 1.0, so pseudocount = 0.5
    expect_equal(pseudocount, 0.5, tolerance = 1e-10)
})

test_that(".compute_pseudocount: handles NAs gracefully", {
    fc_result <- data.frame(
        ctrl_mean = c(5.0, NA, 2.5, 1.0),
        case_mean = c(7.0, 5.0, NA, 3.0),
        log2FC = c(0.485, 0.737, 0.848, 1.585),
        pvalue = c(0.01, 0.05, 0.1, 0.2)
    )
    
    pseudocount <- .compute_pseudocount(fc_result)
    expect_true(is.numeric(pseudocount))
    expect_true(pseudocount > 0)
})

test_that(".compute_pseudocount: returns default when all values are NA", {
    fc_result <- data.frame(
        ctrl_mean = c(NA, NA, NA),
        case_mean = c(NA, NA, NA),
        log2FC = c(NA, NA, NA),
        pvalue = c(NA, NA, NA)
    )
    
    pseudocount <- .compute_pseudocount(fc_result)
    expect_equal(pseudocount, 1e-06)
})

# Tests for .prepare_pair_indices =============================================

test_that(".prepare_pair_indices: correctly maps pair identifiers to indices", {
    pairs <- c("P1", "P2", "P3", "P1", "P2", "P3")
    unique_pairs <- c("P1", "P2", "P3")
    
    pair_indices <- .prepare_pair_indices(pairs, unique_pairs)
    
    expect_equal(pair_indices[[1]], c(1, 4))  # P1
    expect_equal(pair_indices[[2]], c(2, 5))  # P2
    expect_equal(pair_indices[[3]], c(3, 6))  # P3
})

test_that(".prepare_pair_indices: handles single pair per group", {
    pairs <- c("Pair1", "Pair1", "Pair2", "Pair2")
    unique_pairs <- c("Pair1", "Pair2")
    
    pair_indices <- .prepare_pair_indices(pairs, unique_pairs)
    
    expect_equal(length(pair_indices), 2)
    expect_equal(pair_indices[[1]], c(1, 2))
    expect_equal(pair_indices[[2]], c(3, 4))
})

# Tests for .wilcox_effect_sizes_unpaired ====================================

test_that(".wilcox_effect_sizes_unpaired: computes valid U and r statistics", {
    x <- create_test_matrix(nrows = 1, ncols = 6)
    samples <- create_test_samples()
    groups <- c("Normal", "Tumor")
    
    result <- .wilcox_effect_sizes_unpaired(x, feature_idx = 1, samples, groups)
    
    expect_equal(length(result), 2)
    expect_named(result, c("U", "r"))
    expect_true(is.numeric(result["U"]))
    expect_true(is.numeric(result["r"]))
    expect_true(!is.na(result["U"]))
    expect_true(!is.na(result["r"]))
    # r should be between -1 and 1
    expect_true(result["r"] >= -1 && result["r"] <= 1)
})

test_that(".wilcox_effect_sizes_unpaired: returns NA when groups empty", {
    x <- create_test_matrix(nrows = 1, ncols = 6)
    samples <- c("A", "A", "A", "B", "B", "B")
    groups <- c("C", "D")  # Groups not in samples
    
    result <- .wilcox_effect_sizes_unpaired(x, feature_idx = 1, samples, groups)
    
    expect_true(is.na(result["U"]))
    expect_true(is.na(result["r"]))
})

# Tests for .wilcox_effect_sizes_paired ========================================

test_that(".wilcox_effect_sizes_paired: computes valid U and r for paired data", {
    x <- create_test_matrix(nrows = 1, ncols = 6)
    samples <- create_test_samples()
    pairs <- create_test_pairs()
    groups <- c("Normal", "Tumor")
    
    result <- .wilcox_effect_sizes_paired(x, feature_idx = 1, pairs, samples, groups)
    
    expect_equal(length(result), 2)
    expect_named(result, c("U", "r"))
    expect_true(!is.na(result["U"]))
    expect_true(!is.na(result["r"]))
    # r should be between -1 and 1
    expect_true(result["r"] >= -1 && result["r"] <= 1)
})

test_that(".wilcox_effect_sizes_paired: returns NA when insufficient pairs", {
    x <- create_test_matrix(nrows = 1, ncols = 2)
    samples <- c("Normal", "Tumor")
    pairs <- c("P1", "P1")
    groups <- c("Normal", "Tumor")
    
    result <- .wilcox_effect_sizes_paired(x, feature_idx = 1, pairs, samples, groups)
    
    # With only 1 pair, Wilcoxon test needs at least 2 non-zero differences
    expect_true(is.na(result["U"]) || is.na(result["r"]))
})

# Tests for .compute_pvalues_from_permutations ================================

test_that(".compute_pvalues_from_permutations: returns valid p-values", {
    set.seed(123)
    log2_fc <- c(0.5, 1.0, 1.5, 2.0, 0.1)
    perm_mat <- matrix(rnorm(5 * 100, mean = 0, sd = 0.3), nrow = 5, ncol = 100)
    
    pvalues <- .compute_pvalues_from_permutations(log2_fc, perm_mat, nthreads = 1)
    
    expect_equal(length(pvalues), nrow(perm_mat))
    expect_true(all(pvalues >= 0))
    expect_true(all(pvalues <= 1))
})

test_that(".compute_pvalues_from_permutations: applies S019 pseudocount correction", {
    log2_fc <- c(10.0)  # Extreme observed value
    perm_mat <- matrix(c(0.5, 0.5, 0.5, 0.5, 0.5), nrow = 1, ncol = 5)
    
    pvalues <- .compute_pvalues_from_permutations(log2_fc, perm_mat, nthreads = 1)
    
    # With 5 permutations and observed value more extreme than all,
    # p = (0 + 1) / (5 + 1) = 1/6 ≈ 0.1667 (not 0)
    expect_equal(pvalues, 1/6, tolerance = 1e-10)
})

test_that(".compute_pvalues_from_permutations: handles all NA permutations", {
    log2_fc <- c(0.5, 1.0)
    perm_mat <- matrix(NA_real_, nrow = 2, ncol = 100)
    
    pvalues <- .compute_pvalues_from_permutations(log2_fc, perm_mat, nthreads = 1)
    
    expect_equal(pvalues, c(1, 1))  # Return 1 when all permutations are NA
})

# Tests for .format_pvalue_output =============================================

test_that(".format_pvalue_output: creates properly formatted data frame", {
    raw_p <- c(0.01, 0.05, 0.1, 0.2, 0.5)
    adj_p <- c(0.02, 0.1, 0.2, 0.4, 1.0)
    log2fc <- c(0.5, 1.0, 0.8, 0.3, -0.2)
    effect_stats <- list(
        U = c(10, 15, 8, 12, 13),
        r = c(0.6, 0.7, 0.5, 0.4, 0.3)
    )
    group_means <- data.frame(
        Normal_mean = c(5.0, 4.5, 5.2, 5.1, 4.9),
        Tumor_mean = c(5.7, 5.8, 5.7, 5.3, 4.8)
    )
    
    result <- .format_pvalue_output(raw_p, adj_p, log2fc, effect_stats, group_means)
    
    expect_equal(nrow(result), 5)
    expect_equal(ncol(result), 7)  # pvalue, padj, log2FC, U, r, group_means
    expect_named(result, c("pvalue", "padj", "log2FC", "U", "r", "Normal", "Tumor"))
    expect_equal(result$pvalue, raw_p)
    expect_equal(result$padj, adj_p)
})

test_that(".format_pvalue_output: handles named group means", {
    raw_p <- 0.01
    adj_p <- 0.02
    log2fc <- 0.5
    effect_stats <- list(U = 15, r = 0.7)
    group_means <- data.frame(
        Control = 5.0,
        Treatment = 5.7
    )
    
    result <- .format_pvalue_output(raw_p, adj_p, log2fc, effect_stats, group_means)
    
    expect_equal(colnames(result)[6:7], c("Control", "Treatment"))
})

# Tests for .compute_all_effect_sizes =========================================

test_that(".compute_all_effect_sizes: computes U and r for all features", {
    x <- create_test_matrix(nrows = 5, ncols = 6)
    samples <- create_test_samples()
    pairs <- NULL
    groups <- c("Normal", "Tumor")
    
    effect_stats <- .compute_all_effect_sizes(x, samples, pairs, groups, nthreads = 1)
    
    expect_equal(length(effect_stats$U), nrow(x))
    expect_equal(length(effect_stats$r), nrow(x))
    expect_true(all(!is.na(effect_stats$U)))
    expect_true(all(!is.na(effect_stats$r)))
})

test_that(".compute_all_effect_sizes: handles paired data", {
    x <- create_test_matrix(nrows = 3, ncols = 6)
    samples <- create_test_samples()
    pairs <- create_test_pairs()
    groups <- c("Normal", "Tumor")
    
    effect_stats <- .compute_all_effect_sizes(x, samples, pairs, groups, nthreads = 1)
    
    expect_equal(length(effect_stats$U), nrow(x))
    expect_equal(length(effect_stats$r), nrow(x))
})

# Tests for .generate_unpaired_permutations ===================================

test_that(".generate_unpaired_permutations: creates correct matrix dimensions", {
    x <- create_test_matrix(nrows = 5, ncols = 6)
    samples <- create_test_samples()
    control <- "Normal"
    case_group <- "Tumor"
    method <- "mean"
    randomizations <- 100
    pseudocount <- 0.5
    
    perm_mat <- .generate_unpaired_permutations(
        x, samples, control, case_group, method, randomizations,
        pseudocount, robust_loss_type = "huber", robust_scale_method = "mad"
    )
    
    expect_equal(nrow(perm_mat), nrow(x))
    expect_equal(ncol(perm_mat), randomizations)
})

test_that(".generate_unpaired_permutations: produces numeric results", {
    x <- create_test_matrix(nrows = 2, ncols = 6)
    samples <- create_test_samples()
    control <- "Normal"
    case_group <- "Tumor"
    method <- "mean"
    randomizations <- 50
    pseudocount <- 0.5
    
    perm_mat <- .generate_unpaired_permutations(
        x, samples, control, case_group, method, randomizations,
        pseudocount, robust_loss_type = "huber", robust_scale_method = "mad"
    )
    
    expect_true(all(!is.na(perm_mat)))
    expect_true(all(is.numeric(perm_mat)))
})

# Tests for .prepare_pair_indices (comprehensive) ==============================

test_that(".prepare_pair_indices: preserves pair ordering", {
    pairs <- c("P1", "P3", "P2", "P1", "P3", "P2")
    unique_pairs <- unique(pairs)
    
    pair_indices <- .prepare_pair_indices(pairs, unique_pairs)
    
    # Number of lists should match number of unique pairs
    expect_equal(length(pair_indices), length(unique_pairs))
})

test_that(".prepare_pair_indices: creates contiguous index vectors for pairs", {
    pairs <- c("A", "B", "A", "B", "A", "B")
    unique_pairs <- c("A", "B")
    
    pair_indices <- .prepare_pair_indices(pairs, unique_pairs)
    
    expect_equal(pair_indices[[1]], c(1, 3, 5))
    expect_equal(pair_indices[[2]], c(2, 4, 6))
})

# Integration Tests ===========================================================

test_that(".label_shuffling: main function returns valid results for unpaired data", {
    x <- create_test_matrix(nrows = 10, ncols = 6)
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 50, pcorr = "BH", paired = FALSE, nthreads = 1
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(colnames(result) %in% c("pvalue", "padj", "log2FC", "U", "r", "Normal", "Tumor")))
    expect_true(all(result$pvalue >= 0 & result$pvalue <= 1, na.rm = TRUE))
    expect_true(all(result$padj >= 0 & result$padj <= 1, na.rm = TRUE))
    # Adjusted p-values should be >= raw p-values (due to multiple testing correction)
    # Use tolerance for floating-point comparison and handle NAs
    non_na_idx <- !is.na(result$pvalue) & !is.na(result$padj)
    expect_true(all(result$padj[non_na_idx] >= result$pvalue[non_na_idx] - 1e-10))
})

test_that(".label_shuffling: main function returns valid results for paired data", {
    x <- create_test_matrix(nrows = 8, ncols = 6)
    samples <- create_test_samples()
    pairs <- create_test_pairs()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "median",
        randomizations = 50, pcorr = "bonferroni", paired = TRUE,
        paired_method = "swap", pairs = pairs, nthreads = 1
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(result$pvalue >= 0 & result$pvalue <= 1, na.rm = TRUE))
})

test_that(".label_shuffling: different methods produce different results", {
    x <- create_test_matrix(nrows = 5, ncols = 6, seed = 789)
    samples <- create_test_samples()
    
    result_mean <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 50, pcorr = "BH", paired = FALSE, nthreads = 1
    )
    
    result_median <- .label_shuffling(
        x, samples, control = "Normal", method = "median",
        randomizations = 50, pcorr = "BH", paired = FALSE, nthreads = 1
    )
    
    # Results should be different (not identical)
    expect_false(isTRUE(all.equal(result_mean$log2FC, result_median$log2FC)))
})

test_that(".label_shuffling: p-correction methods produce ordered adjusted p-values", {
    x <- create_test_matrix(nrows = 10, ncols = 6)
    samples <- create_test_samples()
    
    result_bh <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 50, pcorr = "BH", paired = FALSE, nthreads = 1
    )
    
    result_bonf <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 50, pcorr = "bonferroni", paired = FALSE, nthreads = 1
    )
    
    # Bonferroni should be more conservative (larger adjusted p-values)
    non_na_idx <- !is.na(result_bonf$padj) & !is.na(result_bh$padj)
    expect_true(all(result_bonf$padj[non_na_idx] >= result_bh$padj[non_na_idx] - 1e-10))
})

test_that(".label_shuffling: handles high randomization counts", {
    x <- create_test_matrix(nrows = 3, ncols = 6)
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 500, pcorr = "BH", paired = FALSE, nthreads = 1
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(result$pvalue >= 0 & result$pvalue <= 1, na.rm = TRUE))
})

# Comprehensive Argument Combination Tests ====================================

test_that(".label_shuffling: works with method='mean'", {
    x <- create_test_matrix(nrows = 5, ncols = 6)
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 30, pcorr = "BH"
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(result$pvalue >= 0 & result$pvalue <= 1, na.rm = TRUE))
})

test_that(".label_shuffling: works with method='median'", {
    x <- create_test_matrix(nrows = 5, ncols = 6)
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "median",
        randomizations = 30, pcorr = "BH"
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(result$pvalue >= 0 & result$pvalue <= 1, na.rm = TRUE))
})

test_that(".label_shuffling: works with pcorr='none'", {
    x <- create_test_matrix(nrows = 5, ncols = 6)
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 30, pcorr = "none"
    )
    
    expect_equal(nrow(result), nrow(x))
    # With pcorr='none', padj should equal pvalue
    expect_true(all(abs(result$padj - result$pvalue) < 1e-10, na.rm = TRUE))
})

test_that(".label_shuffling: works with pcorr='holm'", {
    x <- create_test_matrix(nrows = 5, ncols = 6)
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 30, pcorr = "holm"
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(result$padj >= result$pvalue - 1e-10, na.rm = TRUE))
})

test_that(".label_shuffling: works with pcorr='hochberg'", {
    x <- create_test_matrix(nrows = 5, ncols = 6)
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 30, pcorr = "hochberg"
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(result$padj >= result$pvalue - 1e-10, na.rm = TRUE))
})

test_that(".label_shuffling: works with small randomizations", {
    x <- create_test_matrix(nrows = 5, ncols = 6)
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 10, pcorr = "BH"
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(result$pvalue >= 0 & result$pvalue <= 1, na.rm = TRUE))
})

test_that(".label_shuffling: works with medium randomizations", {
    x <- create_test_matrix(nrows = 5, ncols = 6)
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 100, pcorr = "BH"
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(result$pvalue >= 0 & result$pvalue <= 1, na.rm = TRUE))
})

test_that(".label_shuffling: paired=FALSE with unpaired data", {
    x <- create_test_matrix(nrows = 5, ncols = 6)
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 30, pcorr = "BH", paired = FALSE
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(result$pvalue >= 0 & result$pvalue <= 1, na.rm = TRUE))
})

test_that(".label_shuffling: paired=TRUE with swap method", {
    x <- create_test_matrix(nrows = 5, ncols = 6)
    samples <- create_test_samples()
    pairs <- create_test_pairs()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 30, pcorr = "BH", paired = TRUE,
        paired_method = "swap", pairs = pairs
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(result$pvalue >= 0 & result$pvalue <= 1, na.rm = TRUE))
})

test_that(".label_shuffling: paired=TRUE with signflip method", {
    x <- create_test_matrix(nrows = 5, ncols = 6)
    samples <- create_test_samples()
    pairs <- create_test_pairs()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 30, pcorr = "BH", paired = TRUE,
        paired_method = "signflip", pairs = pairs
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(result$pvalue >= 0 & result$pvalue <= 1, na.rm = TRUE))
})

test_that(".label_shuffling: mean + BH + unpaired", {
    x <- create_test_matrix(nrows = 4, ncols = 6)
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 25, pcorr = "BH", paired = FALSE
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(is.finite(result$pvalue)))
})

test_that(".label_shuffling: median + bonferroni + unpaired", {
    x <- create_test_matrix(nrows = 4, ncols = 6)
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "median",
        randomizations = 25, pcorr = "bonferroni", paired = FALSE
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(is.finite(result$pvalue)))
})

test_that(".label_shuffling: mean + BH + paired swap", {
    x <- create_test_matrix(nrows = 4, ncols = 6)
    samples <- create_test_samples()
    pairs <- create_test_pairs()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 25, pcorr = "BH", paired = TRUE,
        paired_method = "swap", pairs = pairs
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(is.finite(result$pvalue)))
})

test_that(".label_shuffling: median + bonferroni + paired signflip", {
    x <- create_test_matrix(nrows = 4, ncols = 6)
    samples <- create_test_samples()
    pairs <- create_test_pairs()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "median",
        randomizations = 25, pcorr = "bonferroni", paired = TRUE,
        paired_method = "signflip", pairs = pairs
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(is.finite(result$pvalue)))
})

test_that(".label_shuffling: robust_loss_type='huber' (default)", {
    x <- create_test_matrix(nrows = 4, ncols = 6)
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 25, pcorr = "BH",
        robust_loss_type = "huber", robust_scale_method = "mad"
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(is.finite(result$pvalue)))
})

test_that(".label_shuffling: robust_loss_type='tukey'", {
    x <- create_test_matrix(nrows = 4, ncols = 6)
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 25, pcorr = "BH",
        robust_loss_type = "tukey", robust_scale_method = "mad"
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(is.finite(result$pvalue)))
})

test_that(".label_shuffling: robust_scale_method='mad'", {
    x <- create_test_matrix(nrows = 4, ncols = 6)
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 25, pcorr = "BH",
        robust_loss_type = "huber", robust_scale_method = "mad"
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(is.finite(result$pvalue)))
})

test_that(".label_shuffling: control group as second alphabetically", {
    x <- create_test_matrix(nrows = 5, ncols = 6)
    # Reverse the sample order so control is "Tumor" (comes after "Normal")
    samples <- c("Tumor", "Tumor", "Tumor", "Normal", "Normal", "Normal")
    
    result <- .label_shuffling(
        x, samples, control = "Tumor", method = "mean",
        randomizations = 25, pcorr = "BH"
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(result$pvalue >= 0 & result$pvalue <= 1, na.rm = TRUE))
})

test_that(".label_shuffling: different random seed produces different samples", {
    x <- create_test_matrix(nrows = 5, ncols = 6, seed = 100)
    samples <- create_test_samples()
    
    set.seed(111)
    result1 <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 30, pcorr = "BH"
    )
    
    set.seed(222)
    result2 <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 30, pcorr = "BH"
    )
    
    # Results should differ (different random permutations)
    expect_false(isTRUE(all.equal(result1$pvalue, result2$pvalue)))
})

test_that(".label_shuffling: reproducible with same seed", {
    x <- create_test_matrix(nrows = 5, ncols = 6, seed = 100)
    samples <- create_test_samples()
    
    set.seed(333)
    result1 <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 30, pcorr = "BH"
    )
    
    set.seed(333)
    result2 <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 30, pcorr = "BH"
    )
    
    # Results should be identical with same seed
    expect_true(isTRUE(all.equal(result1$pvalue, result2$pvalue)))
})

test_that(".label_shuffling: larger dataset with multiple features", {
    x <- create_test_matrix(nrows = 50, ncols = 10)
    samples <- c(rep("Control", 5), rep("Treatment", 5))
    
    result <- .label_shuffling(
        x, samples, control = "Control", method = "mean",
        randomizations = 40, pcorr = "BH"
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_equal(ncol(result), 7)  # pvalue, padj, log2FC, U, r, Control, Treatment
    expect_true(all(is.finite(result$pvalue)))
    expect_true(all(is.finite(result$padj)))
})

test_that(".label_shuffling: minimum sample size", {
    x <- matrix(rnorm(8), nrow = 2, ncol = 4)
    samples <- c("A", "A", "B", "B")
    
    result <- .label_shuffling(
        x, samples, control = "A", method = "mean",
        randomizations = 20, pcorr = "BH"
    )
    
    expect_equal(nrow(result), 2)
    expect_true(all(is.finite(result$pvalue)))
})

test_that(".label_shuffling: produces correct output structure across all combinations", {
    x <- create_test_matrix(nrows = 6, ncols = 6)
    samples <- create_test_samples()
    pairs <- create_test_pairs()
    
    # Test unpaired
    result_unpaired <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 20, pcorr = "BH", paired = FALSE
    )
    
    expect_true(all(c("pvalue", "padj", "log2FC", "U", "r") %in% colnames(result_unpaired)))
    expect_equal(nrow(result_unpaired), nrow(x))
    
    # Test paired
    result_paired <- .label_shuffling(
        x, samples, control = "Normal", method = "median",
        randomizations = 20, pcorr = "bonferroni", paired = TRUE,
        paired_method = "swap", pairs = pairs
    )
    
    expect_true(all(c("pvalue", "padj", "log2FC", "U", "r") %in% colnames(result_paired)))
    expect_equal(nrow(result_paired), nrow(x))
})

test_that(".label_shuffling: handles extreme p-value corrections properly", {
    x <- create_test_matrix(nrows = 3, ncols = 6)
    samples <- create_test_samples()
    
    # Test with very small p-values (with high randomizations)
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 1000, pcorr = "BH"
    )
    
    # Adjusted p-values should always be >= raw p-values
    non_na_idx <- !is.na(result$pvalue) & !is.na(result$padj)
    expect_true(all(result$padj[non_na_idx] >= result$pvalue[non_na_idx] - 1e-10))
})

# Edge Cases and Missing Scenarios ============================================

test_that(".label_shuffling: single feature (nrows=1)", {
    x <- matrix(rnorm(6), nrow = 1, ncol = 6)
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 25, pcorr = "BH"
    )
    
    expect_equal(nrow(result), 1)
    expect_true(is.finite(result$pvalue))
})

test_that(".label_shuffling: many features (nrows=100)", {
    x <- create_test_matrix(nrows = 100, ncols = 6)
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 15, pcorr = "BH"
    )
    
    expect_equal(nrow(result), 100)
    expect_true(all(is.finite(result$pvalue)))
})

test_that(".label_shuffling: unbalanced groups (2 control, 4 case)", {
    x <- create_test_matrix(nrows = 5, ncols = 6)
    samples <- c("Control", "Control", "Case", "Case", "Case", "Case")
    
    result <- .label_shuffling(
        x, samples, control = "Control", method = "mean",
        randomizations = 20, pcorr = "BH"
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(is.finite(result$pvalue)))
})

test_that(".label_shuffling: unbalanced groups (4 control, 2 case)", {
    x <- create_test_matrix(nrows = 5, ncols = 6)
    samples <- c("Control", "Control", "Control", "Control", "Case", "Case")
    
    result <- .label_shuffling(
        x, samples, control = "Control", method = "mean",
        randomizations = 20, pcorr = "BH"
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(is.finite(result$pvalue)))
})

test_that(".label_shuffling: very small randomizations (1 permutation)", {
    x <- create_test_matrix(nrows = 3, ncols = 6)
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 1, pcorr = "BH"
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(result$pvalue >= 0 & result$pvalue <= 1, na.rm = TRUE))
})

test_that(".label_shuffling: minimal randomizations (2 permutations)", {
    x <- create_test_matrix(nrows = 3, ncols = 6)
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 2, pcorr = "BH"
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(result$pvalue >= 0 & result$pvalue <= 1, na.rm = TRUE))
})

test_that(".label_shuffling: data with uniform values", {
    x <- matrix(5, nrow = 5, ncol = 6)  # All values = 5
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 20, pcorr = "BH"
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(length(result) > 0)
})

test_that(".label_shuffling: data with extreme values", {
    x <- create_test_matrix(nrows = 5, ncols = 6)
    x[1, 1] <- 1e6  # Add extreme value
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 20, pcorr = "BH"
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(is.finite(result$pvalue)))
})

test_that(".label_shuffling: data with negative values", {
    x <- create_test_matrix(nrows = 5, ncols = 6)
    x <- x - 10  # Make all values negative
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 20, pcorr = "BH"
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(length(result) > 0)
})

test_that(".label_shuffling: paired data with 2 pairs", {
    x <- create_test_matrix(nrows = 5, ncols = 4)
    samples <- c("Normal", "Normal", "Tumor", "Tumor")
    pairs <- c("Pair1", "Pair2", "Pair1", "Pair2")
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 20, pcorr = "BH", paired = TRUE,
        paired_method = "swap", pairs = pairs
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(is.finite(result$pvalue)))
})

test_that(".label_shuffling: paired data with many pairs", {
    x <- create_test_matrix(nrows = 5, ncols = 20)
    samples <- rep(c("A", "B"), 10)
    pairs <- rep(1:10, 2)
    
    result <- .label_shuffling(
        x, samples, control = "A", method = "mean",
        randomizations = 15, pcorr = "BH", paired = TRUE,
        paired_method = "swap", pairs = pairs
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(is.finite(result$pvalue)))
})

test_that(".label_shuffling: mean='mean' produces finite log2FC", {
    x <- create_test_matrix(nrows = 5, ncols = 6)
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 20, pcorr = "BH"
    )
    
    expect_true(all(is.finite(result$log2FC)))
})

test_that(".label_shuffling: method='median' produces finite log2FC", {
    x <- create_test_matrix(nrows = 5, ncols = 6)
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "median",
        randomizations = 20, pcorr = "BH"
    )
    
    expect_true(all(is.finite(result$log2FC)))
})

test_that(".label_shuffling: effect sizes U within reasonable bounds", {
    x <- create_test_matrix(nrows = 5, ncols = 6)
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 20, pcorr = "BH"
    )
    
    # U statistic should be non-negative for Wilcoxon test
    non_na_idx <- !is.na(result$U)
    expect_true(all(result$U[non_na_idx] >= 0))
})

test_that(".label_shuffling: effect size r within [-1, 1]", {
    x <- create_test_matrix(nrows = 5, ncols = 6)
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 20, pcorr = "BH"
    )
    
    # r should be between -1 and 1
    non_na_idx <- !is.na(result$r)
    expect_true(all(result$r[non_na_idx] >= -1 & result$r[non_na_idx] <= 1))
})

test_that(".label_shuffling: output column ordering consistent", {
    x <- create_test_matrix(nrows = 5, ncols = 6)
    samples <- create_test_samples()
    
    result1 <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 15, pcorr = "BH"
    )
    
    result2 <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 15, pcorr = "BH"
    )
    
    expect_equal(colnames(result1), colnames(result2))
})

test_that(".label_shuffling: different seeds produce different p-values", {
    x <- create_test_matrix(nrows = 4, ncols = 6)
    samples <- create_test_samples()
    
    set.seed(1000)
    result1 <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 40, pcorr = "BH"
    )
    
    set.seed(2000)
    result2 <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 40, pcorr = "BH"
    )
    
    # At least some p-values should differ
    expect_false(isTRUE(all.equal(result1$pvalue, result2$pvalue)))
})

test_that(".label_shuffling: median method differs from mean method", {
    set.seed(999)
    x <- create_test_matrix(nrows = 4, ncols = 6, seed = 999)
    samples <- create_test_samples()
    
    result_mean <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 20, pcorr = "BH"
    )
    
    result_median <- .label_shuffling(
        x, samples, control = "Normal", method = "median",
        randomizations = 20, pcorr = "BH"
    )
    
    # log2FC should generally differ between methods
    expect_false(isTRUE(all.equal(result_mean$log2FC, result_median$log2FC)))
})

test_that(".label_shuffling: swap and signflip methods produce different results", {
    set.seed(888)
    x <- create_test_matrix(nrows = 4, ncols = 6)
    samples <- create_test_samples()
    pairs <- create_test_pairs()
    
    result_swap <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 20, pcorr = "BH", paired = TRUE,
        paired_method = "swap", pairs = pairs
    )
    
    result_signflip <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 20, pcorr = "BH", paired = TRUE,
        paired_method = "signflip", pairs = pairs
    )
    
    # Results should differ (different permutation schemes)
    expect_false(isTRUE(all.equal(result_swap$pvalue, result_signflip$pvalue)))
})

test_that(".label_shuffling: control='Normal' vs control='Tumor'", {
    x <- create_test_matrix(nrows = 4, ncols = 6)
    samples <- create_test_samples()
    
    result_normal <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 20, pcorr = "BH"
    )
    
    result_tumor <- .label_shuffling(
        x, samples, control = "Tumor", method = "mean",
        randomizations = 20, pcorr = "BH"
    )
    
    # log2FC should flip sign when control group changes
    expect_true(all(sign(result_normal$log2FC) == -sign(result_tumor$log2FC)))
})

test_that(".label_shuffling: column names include both group names", {
    x <- create_test_matrix(nrows = 4, ncols = 6)
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 20, pcorr = "BH"
    )
    
    expect_true("Normal" %in% colnames(result))
    expect_true("Tumor" %in% colnames(result))
})

test_that(".label_shuffling: BH and FDR produce similar results", {
    x <- create_test_matrix(nrows = 4, ncols = 6)
    samples <- create_test_samples()
    
    result_bh <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 20, pcorr = "BH"
    )
    
    result_fdr <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 20, pcorr = "fdr"
    )
    
    # BH and FDR are the same algorithm but may have small numerical differences
    # They should be highly correlated
    expect_true(cor(result_bh$padj, result_fdr$padj, use = "complete.obs") > 0.95)
})

test_that(".label_shuffling: works with integer matrix input", {
    x <- matrix(as.integer(rnorm(30, mean = 100, sd = 20)), nrow = 5, ncol = 6)
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 20, pcorr = "BH"
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(is.finite(result$pvalue)))
})

test_that(".label_shuffling: works with numeric matrix input", {
    x <- matrix(as.numeric(rnorm(30)), nrow = 5, ncol = 6)
    samples <- create_test_samples()
    
    result <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 20, pcorr = "BH"
    )
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(is.finite(result$pvalue)))
})

test_that(".label_shuffling: default parameters work", {
    x <- create_test_matrix(nrows = 5, ncols = 6)
    samples <- create_test_samples()
    
    # Call with all defaults except required parameters
    result <- .label_shuffling(x, samples, control = "Normal", method = "mean")
    
    expect_equal(nrow(result), nrow(x))
    expect_true(all(is.finite(result$pvalue)))
})

test_that(".label_shuffling: explicit defaults match implicit defaults", {
    x <- create_test_matrix(nrows = 4, ncols = 6)
    samples <- create_test_samples()
    
    # Implicit defaults
    result_implicit <- .label_shuffling(
        x, samples, control = "Normal", method = "mean"
    )
    
    # Explicit defaults
    result_explicit <- .label_shuffling(
        x, samples, control = "Normal", method = "mean",
        randomizations = 100, pcorr = "BH", paired = FALSE,
        paired_method = "swap", nthreads = 1, pairs = NULL,
        robust_loss_type = "huber", robust_scale_method = "mad"
    )
    
    expect_equal(colnames(result_implicit), colnames(result_explicit))
    expect_equal(nrow(result_implicit), nrow(result_explicit))
})

test_that(".label_shuffling: all pcorr methods produce valid results", {
    x <- create_test_matrix(nrows = 3, ncols = 6)
    samples <- create_test_samples()
    pcorr_methods <- c("none", "bonferroni", "holm", "hochberg", "BH", "BY", "fdr")
    
    for (method in pcorr_methods) {
        result <- .label_shuffling(
            x, samples, control = "Normal", method = "mean",
            randomizations = 15, pcorr = method
        )
        
        expect_equal(nrow(result), nrow(x), info = paste("Failed for pcorr =", method))
        expect_true(all(is.finite(result$padj)), info = paste("NAs in padj for pcorr =", method))
    }
})

# ============================================================================
# REDISTRIBUTED TESTS FROM test-infrastructure-statistical_validation.R
# ============================================================================

context("Label Shuffling: Type I Error Control & Power")

test_that("label_shuffling produces p-values consistent with null (no effect)", {
    set.seed(100)
    
    # Generate 50 genes with NO effect (null hypothesis)
    n_genes <- 50
    n_samples <- 10
    
    # Random values with no difference between groups
    mat <- matrix(rnorm(n_genes * n_samples, mean = 0.5, sd = 0.1), 
                  nrow = n_genes)
    samples <- c(rep("A", 5), rep("B", 5))
    
    # Run permutation test
    result <- .label_shuffling(mat, samples, control = "A", 
                             method = "mean", randomizations = 100, 
                             pcorr = "none")
    
    p_values <- result[, "pvalue"]
    
    # Under null, p-values should be roughly uniformly distributed
    # With 50 genes and 100 permutations, expect reasonable spread
    # Check that we don't have obviously skewed distribution
    proportion_small <- mean(p_values < 0.3)
    proportion_large <- mean(p_values > 0.5)
    
    # Expect some p-values in different ranges (roughly uniform)
    expect_true(proportion_large > 0.1)  # Some p-values > 0.5
    expect_true(proportion_small > 0.1)  # Some p-values < 0.3
})

test_that("label_shuffling detects true differences (power test)", {
    set.seed(101)
    
    # Create 30 genes: 10 with effect, 20 with no effect
    n_effect <- 10
    n_null <- 20
    n_samples <- 8
    
    # Effect genes: A has mean 0.3, B has mean 0.7
    mat_effect <- matrix(NA, nrow = n_effect, ncol = n_samples)
    mat_effect[, 1:4] <- rnorm(n_effect * 4, mean = 0.3, sd = 0.05)
    mat_effect[, 5:8] <- rnorm(n_effect * 4, mean = 0.7, sd = 0.05)
    
    # Null genes: both groups mean 0.5
    mat_null <- matrix(rnorm(n_null * n_samples, mean = 0.5, sd = 0.1), 
                       nrow = n_null)
    
    mat <- rbind(mat_effect, mat_null)
    samples <- c(rep("A", 4), rep("B", 4))
    
    result <- .label_shuffling(mat, samples, control = "A", 
                             method = "mean", randomizations = 200, 
                             pcorr = "none")
    
    p_effect <- result[1:n_effect, "pvalue"]
    p_null <- result[(n_effect + 1):(n_effect + n_null), "pvalue"]
    
    # Effect genes should have significantly lower p-values
    expect_true(median(p_effect) < median(p_null))
    # At least 50% of effect genes should be detected at p < 0.1
    expect_true(mean(p_effect < 0.1) >= 0.5)
})

context("Label Shuffling: Pseudocount and FDR Correction")

test_that("pseudocount selection is data-driven and prevents negative log(0)", {
    # Data with zeros
    mat <- matrix(c(0, 0, 1, 5, 0, 0, 2, 3), nrow = 2, byrow = TRUE)
    samples <- c("A", "A", "B", "B")
    
    # Auto pseudocount (min/2 = 0.5)
    result <- TSENAT:::.calculate_fc(mat, samples, control = "A", pseudocount = 0)
    
    # Should not produce NaN or Inf
    expect_true(!any(is.nan(result$log2_fold_change)))
    expect_true(!any(is.infinite(result$log2_fold_change)))
})

test_that("BH FDR correction reduces false positives", {
    set.seed(105)
    
    # 100 null genes
    mat <- matrix(rnorm(100 * 8, mean = 0.5, sd = 0.1), nrow = 100)
    samples <- c(rep("A", 4), rep("B", 4))
    
    result <- .label_shuffling(mat, samples, control = "A", 
                             method = "mean", randomizations = 100,
                             pcorr = "BH")
    
    # Under BH correction, expect FDR < 0.05
    fdr <- mean(result[, "padj"] < 0.05)
    expect_true(fdr < 0.1)  # Lenient threshold for stochastic test
})

context("Scale Invariance: Fold Change")

test_that("fold change is scale invariant for log scale", {
    mat1 <- matrix(c(1, 2, 5, 10, 3, 6), nrow = 2, ncol = 3)
    mat2 <- mat1 * 1000  # Scale by 1000x
    
    samples <- c("Normal", "Tumor", "Tumor")
    
    result1 <- TSENAT:::.calculate_fc(mat1, samples, control = "Normal", pseudocount = 1e-6)
    result2 <- TSENAT:::.calculate_fc(mat2, samples, control = "Normal", pseudocount = 1e-6)
    
    # Log2 FC should be identical
    expect_equal(result1[, 4], result2[, 4], 
                 tolerance = 1e-10)
})

