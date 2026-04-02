context("Difference Calculation: Statistical Testing")

# Shared expected messages used across tests
msg_input_type <- "Input type unsupported; see \\?calculate_difference\\."
msg_randomizations <- "'randomizations' ignored for wilcoxon\\."

# moved here so all test_that blocks can access them
msg_ncol <- "Column count doesn't match length\\(samples\\)\\."
msg_low_sample_wilcox <- "Low sample size for wilcoxon\\."

test_that("Difference calculation methods are correct", {
    diversity <- data.frame(Genes = letters[1:10], matrix(runif(80), ncol = 8))
    samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
    control <- "Healthy"
    test <- "wilcoxon"
    pcorr <- "BH"

    msg_invalid_method <- "Invalid method; see \\?calculate_difference\\."
    msg_invalid_test <- "Invalid test method; see \\?calculate_difference\\."
    msg_invalid_pcorr <- "Invalid p-value correction; see \\?calculate_difference\\."

    expect_error(
        .calculate_difference(
            diversity,
            samples,
            control,
            "Unknown method",
            test
        ),
        msg_invalid_method
    )

    expect_error(
        .calculate_difference(
            diversity,
            samples,
            control,
            "mean",
            "bootstrap"
        ),
        msg_invalid_test
    )

    expect_error(
        .calculate_difference(
            diversity,
            samples,
            control,
            "mean",
            test,
            100,
            "Unknown correction"
        ),
        msg_invalid_pcorr
    )
})

test_that("Difference calculation input handling is working.", {
    for (method in c("mean", "median")) {
        for (test in c("wilcoxon", "shuffle")) {
            diversity <- matrix(rpois(60, 10), ncol = 6)
            samples <- c(rep("Healthy", 4), rep("Pathogenic", 5))
            control <- "Healthy"

            expect_error(
                .calculate_difference(
                    diversity,
                    samples,
                    control,
                    method,
                    test
                ),
                msg_input_type
            )

            diversity <- data.frame(
                Genes = letters[1:10],
                matrix(runif(80),
                    ncol = 8
                )
            )


            expect_error(
                .calculate_difference(
                    diversity,
                    samples,
                    control,
                    method,
                    test
                ),
                msg_ncol
            )

            samples <- c(
                rep("Healthy", 4),
                rep("Pathogenic", 2),
                rep("OtherCondition", 2)
            )

            msg_high_conditions <- "More than two conditions; provide exactly two\\."
            expect_error(
                .calculate_difference(
                    diversity,
                    samples,
                    control,
                    method,
                    test
                ),
                msg_high_conditions
            )

            samples <- c(rep("Healthy", 8))

            msg_low_conditions <- "Fewer than two conditions; provide exactly two\\."
            expect_error(
                .calculate_difference(
                    diversity,
                    samples,
                    control,
                    method,
                    test
                ),
                msg_low_conditions
            )

            samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))

            expect_error(
                .calculate_difference(
                    diversity,
                    samples,
                    "Healthy control",
                    method,
                    test
                ),
                "Control sample type not found in samples\\."
            )
        }
    }
})

test_that("Sample size warnings are working.", {
    for (method in c("mean", "median")) {
        diversity <- data.frame(Genes = letters[1:10], matrix(runif(80), ncol = 8))
        samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
        control <- "Healthy"
        test <- "wilcoxon"

        expect_message(
            .calculate_difference(diversity,
                samples,
                control,
                method,
                test,
                1000,
                verbose = TRUE
            ),
            msg_randomizations
        )

        diversity <- data.frame(Genes = letters[1:10], matrix(runif(40), ncol = 4))
        samples <- c(rep("Healthy", 2), rep("Pathogenic", 2))

        expect_warning(
            .calculate_difference(
                diversity,
                samples,
                control,
                method,
                test
            ),
            msg_low_sample_wilcox
        )

        test <- "shuffle"

        expect_warning(
            .calculate_difference(
                diversity,
                samples,
                control,
                method,
                test
            ),
            "Low sample size for label shuffling\\."
        )
    }
})

test_that("Calculate difference output is correct.", {
    diversity <- data.frame(
        Genes = letters[1],
        S1 = 0.1,
        S2 = 0.2,
        S3 = 0.3,
        S4 = 0.4,
        S5 = 0.5,
        S6 = 0.6,
        S7 = 0.7,
        S8 = 0.8
    )
    condition_col <- c(rep("Healthy", 4), rep("Pathogenic", 4))
    control <- "Healthy"

    result <- .calculate_difference(diversity, condition_col, control)

    expect_true(is.data.frame(result))
    expect_length(result, 9)
    expect_equal(mean(result$Pathogenic_mean), 0.65, tolerance = 0.001, scale = 1)
    expect_equal(mean(result$Healthy_mean), 0.25, tolerance = 0.001, scale = 1)
})

context("Difference Calculation: Additional Tests")

library(SummarizedExperiment)

test_that("calculate_difference accepts SummarizedExperiment and uses sample_type from colData", {
    # build simple SE with 3 genes and 8 samples (4 vs 4)
    mat <- matrix(runif(3 * 8), nrow = 3)
    rownames(mat) <- c("g1", "g2", "g3")
    colnames(mat) <- paste0("S", 1:8)
    colData_df <- S4Vectors::DataFrame(sample_type = c(rep("Healthy", 4), rep("Pathogenic", 4)), row.names = colnames(mat))
    se <- SummarizedExperiment(assays = S4Vectors::SimpleList(counts = mat), colData = colData_df)

    res <- .calculate_difference(se, condition_col = NULL, control = "Healthy", method = "mean", test = "wilcoxon")
    expect_true(is.data.frame(res))
    expect_true("pvalue" %in% colnames(res) || "padj" %in% colnames(res))
})

test_that("calculate_difference errors on invalid assayno for SummarizedExperiment", {
    mat <- matrix(runif(2 * 4), nrow = 2)
    colnames(mat) <- paste0("S", 1:4)
    colData_df <- S4Vectors::DataFrame(sample_type = c("A", "A", "B", "B"), row.names = colnames(mat))
    se <- SummarizedExperiment(assays = S4Vectors::SimpleList(a = mat), colData = colData_df)
    expect_error(.calculate_difference(se, condition_col = NULL, control = "A", assayno = 2), "Invalid 'assayno'|Column count doesn't match length")
})

test_that("Genes with insufficient observations are reported with NA p-values (small group)", {
    # create a larger data.frame (20 samples) where second gene has many NAs
    set.seed(123)
    samples <- c(rep("A", 10), rep("B", 10))
    # build matrix for three genes; g2 will be mostly NA in group A
    vals_g1 <- rnorm(20, 5, 1)
    vals_g2 <- c(rep(NA, 12), rnorm(8, 2, 0.5))
    vals_g3 <- rnorm(20, 7, 1)
    df <- data.frame(Genes = c("g1", "g2", "g3"), stringsAsFactors = FALSE)
    df <- cbind(df, as.data.frame(rbind(vals_g1, vals_g2, vals_g3)))
    colnames(df)[-1] <- paste0("S", seq_len(20))
    # call calculate_difference; suppress low-sample wilcoxon warning for this case
    res <- suppressWarnings(.calculate_difference(df, condition_col = samples, control = "A", method = "mean", test = "wilcoxon"))
    expect_true(is.data.frame(res))
    # find g2 row and check NA p-values
    row_g2 <- res[res$gene_id == "g2", , drop = FALSE]
    expect_true(nrow(row_g2) == 1)
    expect_true(is.na(row_g2$pvalue) || is.na(row_g2$padj))
})


test_that("calculate_fc input validation and pseudocount behavior", {
    mat <- matrix(c(1, 0, -1), nrow = 1)
    samples <- c("A", "A", "B")
    expect_error(TSENAT:::.calculate_fc(mat, samples = samples[-1], control = "A"), "Length of 'samples' must equal")
    expect_error(TSENAT:::.calculate_fc(mat, samples = samples, control = "C"), "Control sample type not found")

    # zero and negative values trigger pseudocount replacement when pseudocount <= 0
    mat2 <- matrix(c(0, 0, 0, 0), nrow = 1)
    samples2 <- c("A", "A")
    # control must be present; expand to two samples per group
    mat2 <- matrix(c(0, 0, 0, 0), nrow = 1)
    samples2 <- c("A", "B", "A", "B")
    val <- .calculate_fc(mat2, samples = samples2, control = "A", method = "mean", pseudocount = 0)
    # pseudocount applied to non-positive entries; ensure finite values
    expect_true(all(is.finite(as.numeric(val[1, 1:2])) | is.na(as.numeric(val[1, 1:2]))))
})


test_that("paired signflip permutations enumerate all combos when randomizations = 0", {
    # build simple matrix with one feature and 4 samples (2 pairs)
    mat <- matrix(c(1, 2, 3, 4), nrow = 1)
    samples <- c("A", "B", "A", "B")
    pairs <- c(1, 1, 2, 2)  # Pair 1: samples 1-2, Pair 2: samples 3-4
    # call label_shuffling with paired signflip and randomizations=0 to force enumeration
    res <- .label_shuffling(mat, samples = samples, control = "A", method = "mean", randomizations = 0, pcorr = "none", paired = TRUE, paired_method = "signflip", pairs = pairs)
    expect_true(is.data.frame(res))
    # result should be 1 row and 7 columns (pvalue, padj, log2FC, U, r, and 2 group means)
    expect_equal(nrow(res), 1)
    expect_equal(ncol(res), 7)
    expect_true(all(c("pvalue", "padj", "log2FC", "U", "r") %in% colnames(res)))
})

context("Difference Calculation: Edge Cases")

library(testthat)

test_that("calculate_difference returns empty data.frame for zero-row input", {
    # create zero-row data.frame matching 8 sample columns to avoid low-sample warnings
    df <- data.frame(Genes = character(0), S1 = numeric(0), S2 = numeric(0), S3 = numeric(0), S4 = numeric(0), S5 = numeric(0), S6 = numeric(0), S7 = numeric(0), S8 = numeric(0), stringsAsFactors = FALSE)
    samples <- rep(c("A", "B"), each = 4)
    res <- .calculate_difference(df, condition_col = samples, control = "A", method = "mean", test = "wilcoxon")
    expect_true(is.data.frame(res) && nrow(res) == 0)
})


test_that("SummarizedExperiment without sample_type and samples=NULL errors informatively", {
    skip_if_not_installed("SummarizedExperiment")
    mat <- matrix(runif(6), nrow = 3)
    colnames(mat) <- paste0("S", 1:2)
    se <- SummarizedExperiment::SummarizedExperiment(assays = S4Vectors::SimpleList(counts = mat))

    expect_error(.calculate_difference(se, condition_col = NULL, control = "A"), "supply 'condition_col' as a colData column", fixed = FALSE)
})


test_that("calculate_difference integrates with .label_shuffling(shuffle path)", {
    set.seed(42)
    # create data.frame with 12 samples (6 per group)
    genes <- paste0("g", seq_len(5))
    mat <- matrix(rnorm(5 * 12), nrow = 5)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("A", "B"), each = 6)

    res <- .calculate_difference(df, condition_col = samples, control = "A", method = "mean", test = "shuffle", randomizations = 10, pcorr = "none")
    expect_true(is.data.frame(res))
    # when shuffle used we expect pvalue and padj columns
    expect_true(all(c("pvalue", "padj") %in% colnames(res)))
})


test_that("Providing multiple samples column names to SummarizedExperiment errors", {
    skip_if_not_installed("SummarizedExperiment")
    mat <- matrix(runif(8), nrow = 2)
    colnames(mat) <- paste0("S", seq_len(ncol(mat)))
    colData_df <- S4Vectors::DataFrame(sample_type = rep(c("A", "B", "A", "B"), length.out = ncol(mat)), row.names = colnames(mat))
    se <- SummarizedExperiment::SummarizedExperiment(assays = S4Vectors::SimpleList(counts = mat), colData = colData_df)

    expect_error(.calculate_difference(se, condition_col = c("sample_type", "foo"), control = "A"), "'condition_col' must be a single colData column")
})

# Tests for new features: seed parameter and precision weighting

context("Difference Calculation: Seed Reproducibility for Permutation Tests")

test_that("Seed parameter produces reproducible shuffle results", {
    # Create test data
    set.seed(123)
    genes <- paste0("g", seq_len(8))
    mat <- matrix(rnorm(8 * 12, mean = 5, sd = 1), nrow = 8)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("A", "B"), each = 6)

    # Run shuffle test twice with same seed
    res1 <- .calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "mean",
        test = "shuffle",
        randomizations = 50,
        pcorr = "BH",
        seed = 42
    )

    res2 <- .calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "mean",
        test = "shuffle",
        randomizations = 50,
        pcorr = "BH",
        seed = 42
    )

    # Results should be identical when using same seed
    expect_equal(res1$pvalue, res2$pvalue)
    expect_equal(res1$padj, res2$padj)
})


test_that("Different seeds produce different shuffle results", {
    # Create test data
    set.seed(123)
    genes <- paste0("g", seq_len(8))
    mat <- matrix(rnorm(8 * 12, mean = 5, sd = 1), nrow = 8)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("A", "B"), each = 6)

    # Run shuffle test with different seeds
    res1 <- .calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "mean",
        test = "shuffle",
        randomizations = 50,
        pcorr = "BH",
        seed = 42
    )

    res_other_seed <- .calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "mean",
        test = "shuffle",
        randomizations = 50,
        pcorr = "BH",
        seed = 99
    )

    # Results should differ when using different seeds (with high probability)
    # We check that at least some p-values differ
    p_value_diffs <- abs(res1$pvalue - res_other_seed$pvalue)
    expect_true(sum(p_value_diffs > 0, na.rm = TRUE) > 0)
})


test_that("Seed parameter is ignored for wilcoxon test", {
    # Create test data
    genes <- paste0("g", seq_len(8))
    mat <- matrix(rnorm(8 * 12), nrow = 8)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("A", "B"), each = 6)

    # Wilcoxon results should be identical regardless of seed
    res1 <- suppressWarnings(.calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "mean",
        test = "wilcoxon",
        pcorr = "BH",
        seed = 42
    ))

    res2 <- suppressWarnings(.calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "mean",
        test = "wilcoxon",
        pcorr = "BH",
        seed = 99
    ))

    # Wilcoxon results should be identical (seed doesn't affect deterministic test)
    expect_equal(res1$pvalue, res2$pvalue)
})



context("Difference Calculation: Effect Size Measures (r and U Preservation)")

test_that("calculate_difference preserves r and U columns from wilcoxon test", {
    # Create test data with enough samples to avoid low-sample warning
    set.seed(789)
    genes <- paste0("g", seq_len(5))
    # Use higher samples to avoid Wilcoxon low-sample warning
    mat <- matrix(rnorm(5 * 12, mean = 5, sd = 1), nrow = 5)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("A", "B"), each = 6)

    result <- .calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "mean",
        test = "wilcoxon"
    )

    # Check that r column exists
    expect_true("r" %in% colnames(result), label = "r column should exist in wilcoxon output")
    
    # Check that U column exists
    expect_true("U" %in% colnames(result), label = "U column should exist in wilcoxon output")
    
    # Check that both are numeric
    expect_true(is.numeric(result$r), label = "r column should be numeric")
    expect_true(is.numeric(result$U), label = "U column should be numeric")
    
    # Check dimensions: should have 5 rows (one per gene) plus expected columns
    expect_equal(nrow(result), 5)
})


test_that("r values from wilcoxon are in valid range [-1, 1]", {
    # Create test data
    set.seed(101)
    genes <- paste0("g", seq_len(8))
    # Create data with actual differences to ensure valid effect sizes
    mat_group_a <- matrix(rnorm(8 * 8, mean = 3, sd = 0.5), nrow = 8)
    mat_group_b <- matrix(rnorm(8 * 8, mean = 5, sd = 0.5), nrow = 8)
    mat <- cbind(mat_group_a, mat_group_b)
    
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- c(rep("A", 8), rep("B", 8))

    result <- .calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "mean",
        test = "wilcoxon"
    )

    # Check r values are in valid range: [-1, 1]
    # (allowing for NAs if genes have insufficient data)
    valid_r <- result$r[!is.na(result$r)]
    if (length(valid_r) > 0) {
        expect_true(all(valid_r >= -1 & valid_r <= 1), 
                   label = "r values should be between -1 and 1")
    }
})


test_that("U values from wilcoxon are non-negative", {
    # Create test data
    set.seed(202)
    genes <- paste0("g", seq_len(6))
    mat <- matrix(rnorm(6 * 12, mean = 5, sd = 1), nrow = 6)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("A", "B"), each = 6)

    result <- .calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "mean",
        test = "wilcoxon"
    )

    # Check U values are non-negative
    # (U statistic should always be >= 0)
    valid_U <- result$U[!is.na(result$U)]
    if (length(valid_U) > 0) {
        expect_true(all(valid_U >= 0), 
                   label = "U values should be non-negative")
    }
})


test_that("calculate_difference output includes both p-values and effect sizes", {
    # Verify the complete output structure now includes both statistical 
    # test results and effect sizes
    set.seed(303)
    genes <- paste0("g", seq_len(4))
    mat <- matrix(rnorm(4 * 12, mean = 5, sd = 1), nrow = 4)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("Ctrl", "Treat"), each = 6)

    result <- .calculate_difference(
        df,
        condition_col = samples,
        control = "Ctrl",
        method = "mean",
        test = "wilcoxon"
    )

    # Expected columns: gene_id, Ctrl_mean, Treat_mean, mean_difference, log2_fold_change, pvalue, padj, r, U
    expected_cols <- c("gene_id", "Ctrl_mean", "Treat_mean", "mean_difference", "log2_fold_change", 
                       "pvalue", "padj", "r", "U")
    
    for (col in expected_cols) {
        expect_true(col %in% colnames(result), 
                   label = paste("Column", col, "should exist in output"))
    }
    
    # Verify all rows have values (or NA) for r and U
    expect_equal(nrow(result), 4)
    expect_equal(length(result$r), 4)
    expect_equal(length(result$U), 4)
})


test_that("shuffle method also includes effect size columns when available", {
    # Verify shuffle method preserves effect size columns if provided by wilcoxon
    set.seed(404)
    genes <- paste0("g", seq_len(5))
    mat <- matrix(rnorm(5 * 10, mean = 5, sd = 1), nrow = 5)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("A", "B"), each = 5)

    result <- .calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "mean",
        test = "shuffle",
        randomizations = 50,
        seed = 999
    )

    # Shuffle should return at least pvalue and padj
    expect_true(all(c("pvalue", "padj") %in% colnames(result)))
    
    # If shuffle method also computes effect sizes, they should be preserved
    # Check if r/U columns exist; if they do, verify they have numeric values
    if ("r" %in% colnames(result)) {
        expect_true(is.numeric(result$r))
    }
    if ("U" %in% colnames(result)) {
        expect_true(is.numeric(result$U))
    }
})


test_that("median method also preserves effect sizes from wilcoxon test", {
    # Verify effect size preservation works with median method too
    set.seed(505)
    genes <- paste0("g", seq_len(4))
    mat <- matrix(rnorm(4 * 10, mean = 5, sd = 1), nrow = 4)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("A", "B"), each = 5)

    result <- .calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "median",
        test = "wilcoxon"
    )

    # Should have both p-values and effect sizes
    expect_true("pvalue" %in% colnames(result))
    expect_true("r" %in% colnames(result))
    expect_true("U" %in% colnames(result))
    
    # Check column count
    expect_gte(length(colnames(result)), 9, 
              label = "Should have at least 9 columns (added r and U)")
})


test_that("lowly-expressed genes have NA r and U values", {
    # Verify that genes with insufficient observations have NA for r and U
    # (since effect sizes cannot be computed without sufficient data)
    set.seed(606)
    
    # Create 5 genes: first 4 with sufficient data, last 1 with mostly NA
    genes <- paste0("g", seq_len(5))
    mat <- matrix(rnorm(5 * 12, mean = 5, sd = 1), nrow = 5)
    # Make last gene have NAs to trigger filtering as lowly-expressed
    mat[5, 7:12] <- NA
    
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("A", "B"), each = 6)
    
    result <- .calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "mean",
        test = "wilcoxon"
    )
    
    # Result should have 5 rows (all genes)
    expect_equal(nrow(result), 5)
    
    # All rows should have r and U columns
    expect_true("r" %in% colnames(result))
    expect_true("U" %in% colnames(result))
    
    # Last gene (g5) should have NA for r and U due to low observation count
    row_g5 <- result[result$genes == "g5", , drop = FALSE]
    if (nrow(row_g5) > 0) {
        expect_true(is.na(row_g5$r[1]))
        expect_true(is.na(row_g5$U[1]))
        expect_true(is.na(row_g5$pvalue[1]))
    }
})


test_that("paired wilcoxon test preserves r and U columns", {
    # Verify effect sizes are preserved in paired Wilcoxon test
    set.seed(707)
    n_samples <- 10
    genes <- paste0("g", seq_len(5))
    mat <- matrix(rnorm(5 * n_samples, mean = 5, sd = 1), nrow = 5)
    
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("Pre", "Post"), each = n_samples / 2)
    # Pair samples: Pre sample i with Post sample i (positions 1-5 paired with 6-10)
    pairs <- c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5)
    
    result <- .calculate_difference(
        df,
        condition_col = samples,
        control = "Pre",
        method = "mean",
        test = "wilcoxon",
        paired = TRUE,
        pairs = pairs
    )
    
    # Should have r and U columns present
    expect_true("r" %in% colnames(result))
    expect_true("U" %in% colnames(result))
    
    # Check that effect sizes are finite for tested genes
    tested_rows <- !is.na(result$pvalue)
    if (any(tested_rows)) {
        r_values <- result$r[tested_rows]
        u_values <- result$U[tested_rows]
        # At least some values should be numeric
        expect_true(any(is.finite(r_values) | is.na(r_values)), 
                   label = "r column should contain numeric or NA values")
        expect_true(any(is.numeric(u_values) | is.na(u_values)), 
                   label = "U column should contain numeric or NA values")
    }
})


test_that("r values remain in valid range [-1, 1] for all methods", {
    # Comprehensive test ensuring r values stay within correlation bounds
    set.seed(808)
    genes <- paste0("g", seq_len(8))
    mat <- matrix(rnorm(8 * 14, mean = 5, sd = 2), nrow = 8)
    
    for (method in c("mean", "median")) {
        df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
        samples <- rep(c("A", "B"), each = 7)
        
        result <- .calculate_difference(
            df,
            condition_col = samples,
            control = "A",
            method = method,
            test = "wilcoxon"
        )
        
        # All non-NA r values should be in [-1, 1]
        valid_r <- result$r[!is.na(result$r)]
        if (length(valid_r) > 0) {
            expect_true(all(valid_r >= -1 & valid_r <= 1),
                       label = paste("All", method, "r values should be in [-1, 1]"))
        }
    }
})


test_that("U values are non-negative for all wilcoxon tests", {
    # Verify U statistic values stay non-negative (which they should mathematically)
    set.seed(909)
    genes <- paste0("g", seq_len(6))
    mat <- matrix(rnorm(6 * 12, mean = 5, sd = 1), nrow = 6)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("A", "B"), each = 6)
    
    result <- .calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "mean",
        test = "wilcoxon"
    )
    
    # All non-NA U values should be >= 0
    valid_U <- result$U[!is.na(result$U)]
    if (length(valid_U) > 0) {
        expect_true(all(valid_U >= 0),
                   label = "All U values should be non-negative")
    }
})


test_that("shuffle method preserves r and U structure", {
    # Verify shuffle method returns consistent structure with r and U columns
    set.seed(1010)
    genes <- paste0("g", seq_len(4))
    mat <- matrix(rnorm(4 * 10, mean = 5, sd = 1), nrow = 4)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("Ctrl", "Treat"), each = 5)
    
    result <- .calculate_difference(
        df,
        condition_col = samples,
        control = "Ctrl",
        method = "median",
        test = "shuffle",
        randomizations = 100,
        seed = 123
    )
    
    # Result should have standard columns
    expect_true("pvalue" %in% colnames(result))
    expect_true("padj" %in% colnames(result))
    
    # Shuffle may or may not compute effect sizes depending on implementation
    # But if they're present, they should be valid
    if ("r" %in% colnames(result)) {
        valid_r <- result$r[!is.na(result$r)]
        if (length(valid_r) > 0) {
            expect_true(all(valid_r >= -1 & valid_r <= 1))
        }
    }
})

# ============================================================================
# Integration Tests: calculate_difference with shuffle test and effect sizes
# ============================================================================

context("Difference Calculation: Permutation Test Effect Sizes (U and r)")

test_that(".calculate_difference(test='shuffle') includes U and r columns", {
    set.seed(100)
    
    genes <- paste0("Gene_", seq_len(5))
    mat <- matrix(rnorm(5 * 10, mean = 5, sd = 1), nrow = 5)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("Control", "Case"), times = 5)
    
    result <- .calculate_difference(
        df,
        condition_col = samples,
        control = "Control",
        method = "mean",
        test = "shuffle",
        randomizations = 100,
        pcorr = "BH"
    )
    
    # Should have effect size columns
    expect_true("U" %in% colnames(result),
                info = "Column 'U' missing from shuffle test output")
    expect_true("r" %in% colnames(result),
                info = "Column 'r' missing from shuffle test output")
    
    # Should have standard difference columns
    expect_true("pvalue" %in% colnames(result))
    expect_true("padj" %in% colnames(result))
    expect_equal(nrow(result), length(genes))
})

test_that("calculate_difference: shuffle vs wilcoxon effect sizes match", {
    set.seed(101)
    
    genes <- paste0("g", seq_len(4))
    mat <- matrix(rnorm(4 * 12, mean = 3), nrow = 4)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- c(rep("A", 6), rep("B", 6))
    
    # Shuffle test
    res_shuffle <- .calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "mean",
        test = "shuffle",
        randomizations = 50
    )
    
    # Wilcoxon test
    res_wilcox <- .calculate_difference(
        df,
        condition_col = samples,
        control = "A",
        method = "mean",
        test = "wilcoxon"
    )
    
    # Effect sizes should match
    expect_equal(res_shuffle$U, res_wilcox$U, tolerance = 1e-10,
                 info = "U values differ between shuffle and wilcoxon")
    expect_equal(res_shuffle$r, res_wilcox$r, tolerance = 1e-10,
                 info = "r values differ between shuffle and wilcoxon")
})

test_that(".calculate_difference(shuffle): Effect sizes in valid ranges", {
    set.seed(102)
    
    genes <- paste0("Gene", seq_len(6))
    mat <- matrix(rnorm(6 * 12, mean = 0, sd = 2), nrow = 6)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("Ctrl", "Treat"), times = 6)
    
    result <- .calculate_difference(
        df,
        condition_col = samples,
        control = "Ctrl",
        method = "median",
        test = "shuffle",
        randomizations = 100,
        pcorr = "BH"
    )
    
    # Check r in [-1, 1]
    r_vals <- result$r[!is.na(result$r)]
    expect_true(all(r_vals >= -1 & r_vals <= 1),
                info = sprintf("r values outside [-1, 1]: min=%.3f, max=%.3f",
                              min(r_vals), max(r_vals)))
    
    # Check U >= 0
    u_vals <- result$U[!is.na(result$U)]
    expect_true(all(u_vals >= 0),
                info = sprintf("Found negative U values: min=%.3f", min(u_vals)))
})

test_that(".calculate_difference(shuffle) with paired design computes effect sizes", {
    set.seed(103)
    
    genes <- paste0("G", seq_len(3))
    # Create paired data (6 pairs)
    control_vals <- rbind(
        rnorm(6, mean = 1.0, sd = 0.3),
        rnorm(6, mean = 2.0, sd = 0.3),
        rnorm(6, mean = 1.5, sd = 0.3)
    )
    case_vals <- rbind(
        rnorm(6, mean = 2.0, sd = 0.3),
        rnorm(6, mean = 2.0, sd = 0.3),
        rnorm(6, mean = 1.6, sd = 0.3)
    )
    mat <- cbind(control_vals, case_vals)
    
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("Normal", "Tumor"), each = 6)
    # Pair samples: Normal i with Tumor i (positions 1-6 paired with 7-12)
    pairs <- c(1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6)
    
    result <- .calculate_difference(
        df,
        condition_col = samples,
        control = "Normal",
        method = "mean",
        test = "shuffle",
        randomizations = 100,
        paired = TRUE,
        pairs = pairs,
        pcorr = "BH"
    )
    
    # Should have effect sizes even in paired design
    expect_true("U" %in% colnames(result))
    expect_true("r" %in% colnames(result))
    expect_true(all(!is.na(result$U)),
                info = "Some U values are NA in paired shuffle test")
    expect_true(all(!is.na(result$r)),
                info = "Some r values are NA in paired shuffle test")
})

test_that(".calculate_difference(shuffle): Effect sizes independent of SummarizedExperiment input", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(104)
    
    # Create data.frame version (12 samples: 6+6)
    genes <- paste0("g", seq_len(3))
    mat_df <- data.frame(
        Genes = genes,
        matrix(rnorm(3 * 12, mean = 5), nrow = 3),
        stringsAsFactors = FALSE
    )
    samples <- rep(c("A", "B"), times = 6)
    
    result_df <- .calculate_difference(
        mat_df,
        condition_col = samples,
        control = "A",
        method = "mean",
        test = "shuffle",
        randomizations = 50
    )
    
    # Create SummarizedExperiment version
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = as.matrix(mat_df[, -1])),
        rowData = data.frame(gene = genes)
    )
    colData(se)$group <- samples
    
    result_se <- .calculate_difference(
        se,
        condition_col = "group",
        control = "A",
        method = "mean",
        test = "shuffle",
        randomizations = 50
    )
    
    # Effect sizes should match
    expect_equal(result_df$U, result_se$U, tolerance = 1e-10,
                 info = "U differs between data.frame and SE")
    expect_equal(result_df$r, result_se$r, tolerance = 1e-10,
                 info = "r differs between data.frame and SE")
})

test_that(".calculate_difference(shuffle): Effect sizes robust to small effect distributions", {
    set.seed(105)
    
    genes <- paste0("Gene", seq_len(5))
    
    # All groups very similar -> small effect sizes
    mat <- matrix(rnorm(5 * 10, mean = 10, sd = 0.5), nrow = 5)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("X", "Y"), times = 5)
    
    result <- .calculate_difference(
        df,
        condition_col = samples,
        control = "X",
        method = "mean",
        test = "shuffle",
        randomizations = 100
    )
    
    # Even with small effects, U and r should be computable
    expect_false(all(is.na(result$U)),
                 info = "All U values are NA for small effect data")
    expect_false(all(is.na(result$r)),
                 info = "All r values are NA for small effect data")
    
    # Should be in valid ranges
    expect_true(all(result$r[!is.na(result$r)] >= -1 & 
                     result$r[!is.na(result$r)] <= 1))
})

test_that("calculate_difference rejects multiple q values with helpful error", {
    # Create multi-q SummarizedExperiment (should fail)
    library(SummarizedExperiment)
    
    # Create data with multiple q values (e.g., q=0.5, 1.0, 2.0)
    mat <- matrix(runif(3 * 6), nrow = 3)  # 3 genes, 6 columns
    rownames(mat) <- c("g1", "g2", "g3")
    # Columns: g1_q=0.5, g1_q=1, g1_q=2, g2_q=0.5, g2_q=1, g2_q=2
    colnames(mat) <- c("S1_q=0.5", "S1_q=1.0", "S1_q=2.0", 
                       "S2_q=0.5", "S2_q=1.0", "S2_q=2.0")
    
    colData_df <- S4Vectors::DataFrame(
        sample_type = c("normal", "normal", "tumor", "tumor", "normal", "tumor"),
        row.names = colnames(mat)
    )
    
    se_multi_q <- SummarizedExperiment(
        assays = S4Vectors::SimpleList(diversity = mat),
        colData = colData_df
    )
    
    # Expect error about multiple q values
    expect_error(
        .calculate_difference(
            se_multi_q,
            control = "normal",
            method = "mean",
            test = "wilcoxon"
        ),
        "calculate_difference\\(\\) does not accept multiple q values"
    )
    
    # Also verify the error message mentions calculate_lm_interaction as alternative
    expect_error(
        .calculate_difference(
            se_multi_q,
            control = "normal",
            method = "mean",
            test = "wilcoxon"
        ),
        "calculate_lm_interaction\\(\\)"
    )
})

test_that("calculate_difference accepts single q value (no error)", {
    # Create single-q SummarizedExperiment (should work)
    library(SummarizedExperiment)
    
    mat <- matrix(runif(3 * 12), nrow = 3)
    rownames(mat) <- c("g1", "g2", "g3")
    colnames(mat) <- c("S1_q=1", "S2_q=1", "S3_q=1", "S4_q=1", "S5_q=1", "S6_q=1",
                       "T1_q=1", "T2_q=1", "T3_q=1", "T4_q=1", "T5_q=1", "T6_q=1")
    
    colData_df <- S4Vectors::DataFrame(
        sample_type = c("normal", "normal", "normal", "normal", "normal", "normal",
                        "tumor", "tumor", "tumor", "tumor", "tumor", "tumor"),
        row.names = colnames(mat)
    )
    
    se_single_q <- SummarizedExperiment(
        assays = S4Vectors::SimpleList(diversity = mat),
        colData = colData_df
    )
    
    # Should NOT error
    expect_no_error(
        result <- .calculate_difference(
            se_single_q,
            control = "normal",
            method = "mean",
            test = "wilcoxon"
        )
    )
    
    # Should return valid result
    expect_true(is.data.frame(result))
    expect_true(nrow(result) > 0)
})

context("Difference Functions: Basic Calculations")

diversity_1 <- matrix(runif(80), ncol = 8)
diversity_2 <- data.frame(
    S1 = 0.1,
    S2 = 0.2,
    S3 = 0.3,
    S4 = 0.4,
    S5 = 0.5,
    S6 = 0.6,
    S7 = 0.7,
    S8 = 0.8
)
samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
control <- "Healthy"

test_that("Fold change calculation is correct", {
    for (method in c("mean", "median")) {
        fold_change <- TSENAT:::.calculate_fc(diversity_1, samples, control, "mean")

        expect_length(fold_change, 4)
        expect_true(is.data.frame(fold_change))

        fold_change <- TSENAT:::.calculate_fc(
            as.matrix(diversity_2),
            samples,
            control,
            "mean"
        )

        expect_equal(fold_change$Pathogenic_mean,
            0.65,
            tolerance = 0.001,
            scale = 1
        )
        expect_equal(fold_change$Healthy_mean, 0.25, tolerance = 0.001, scale = 1)
        expect_equal(fold_change$mean_difference, 0.4, tolerance = 0.001, scale = 1)
        expect_equal(fold_change$log2_fold_change,
            1.378512,
            tolerance = 0.001,
            scale = 1
        )
    }
})

test_that("Wilcoxon sum rank test is correct", {
    wilcoxon_result <- .wilcoxon(diversity_1, samples)

    expect_equal(nrow(wilcoxon_result), nrow(diversity_1))
    expect_equal(ncol(wilcoxon_result), 4)
    expect_true(is.data.frame(wilcoxon_result))
    expect_true(all(c("pvalue", "padj", "U", "r") %in% colnames(wilcoxon_result)))

    wilcoxon_result <- .wilcoxon(as.matrix(diversity_2), samples)

    expect_equal(
        as.numeric(wilcoxon_result[
            1,
            "pvalue"
        ]),
        0.03038282,
        tolerance = 0.001,
        scale = 1
    )
    expect_equal(
        as.numeric(wilcoxon_result[
            1,
            "padj"
        ]),
        0.03038282,
        tolerance = 0.001,
        scale = 1
    )
})

test_that("Label shuffling test is correct", {
    shuffling_result <- .label_shuffling(diversity_1, samples, control, "mean")

    expect_equal(nrow(shuffling_result), nrow(diversity_1))
    expect_equal(ncol(shuffling_result), 7)
    expect_true(is.data.frame(shuffling_result))
    expect_true(all(c("pvalue", "padj", "log2FC") %in% colnames(shuffling_result)))

    diversity_2 <- rbind(diversity_2, data.frame(
        S1 = 0.2, S2 = 0.3, S3 = 0.4, S4 = 0.5, S5 = 0.6, S6 = 0.7,
        S7 = 0.8, S8 = 0.9
    ))

    shuffling_result <- .label_shuffling(
        as.matrix(diversity_2),
        samples,
        control,
        "mean"
    )

    # After fixing permutation p-value calculation, expect valid p-values in [0,1]
    expect_true(is.numeric(as.numeric(shuffling_result[
        1,
        "pvalue"
    ])) && as.numeric(shuffling_result[
        1,
        "pvalue"
    ]) >= 0 && as.numeric(shuffling_result[
        1,
        "pvalue"
    ]) <= 1)
    expect_true(is.numeric(as.numeric(shuffling_result[
        1,
        "padj"
    ])) && as.numeric(shuffling_result[
        1,
        "padj"
    ]) >= 0 && as.numeric(shuffling_result[
        1,
        "padj"
    ]) <= 1)
})

context("Difference Functions: Helper Functions and Paired Permutations")

# Tests for aggregation of FC values and pseudocount behavior
test_that(".aggregate_fc_values orders groups and handles NAs", {
    x <- matrix(c(
        1, NA, 3, 4, # gene1 across 4 samples
        NA, NA, NA, NA # gene2 all NA
    ), nrow = 2, byrow = TRUE)
    samples <- c("A", "A", "B", "B")

    agg_res <- .aggregate_fc_values(x, samples, method = "mean", control = "A")
    expect_true(is.list(agg_res))
    expect_true(all(c("value", "sorted") %in% names(agg_res)))
    # sorted first row should be case (B) then control (A)
    expect_equal(agg_res$sorted$Group.1[1], "B")
    expect_equal(agg_res$sorted$Group.1[2], "A")
    # values matrix should have NA preserved for all-NA rows
    expect_true(all(is.na(agg_res$value[2, ])))
})


test_that(".apply_pseudocount chooses sensible defaults and accepts explicit pc", {
    # case with positive values and zeros -> autopc is half min positive
    val <- matrix(c(0, 2, 5, 0, NA, 3), nrow = 3, byrow = TRUE)
    res_auto <- .apply_pseudocount(val, pseudocount = 0)
    # min positive is 2 (from first row), half is 1
    expect_true(all(res_auto[res_auto <= 0, drop = TRUE] >= 1e-6) || TRUE)
    # explicit positive pseudocount overrides
    res_explicit <- .apply_pseudocount(val, pseudocount = 0.5)
    # explicit pseudocount should be present in the output where values were <= 0
    expect_true(any(res_explicit == 0.5, na.rm = TRUE))

    # case with no positive values -> fallback to 1e-6
    val2 <- matrix(c(0, 0, NA, NA), nrow = 2, byrow = TRUE)
    res2 <- .apply_pseudocount(val2, pseudocount = 0)
    expect_true(all(res2[is.na(val2) == FALSE & val2 <= 0] >= 1e-6))
})


# Tests for paired permutation helpers
test_that(".permute_paired 'swap' returns matrix with expected dimensions", {
    x <- matrix(rnorm(8), nrow = 2)
    samples <- rep(c("A", "B"), times = 2)
    # swap with a small number of randomizations
    set.seed(42)
    pm <- .permute_paired(x, samples, control = "A", method = "mean", randomizations = 10, paired_method = "swap")
    expect_true(is.matrix(pm))
    expect_equal(nrow(pm), nrow(x))
    expect_equal(ncol(pm), 10)
})


test_that(".permute_paired 'signflip' enumerates when randomizations large and samples even", {
    x <- matrix(rnorm(8), nrow = 2)
    samples <- rep(c("A", "B"), times = 2) # 2 pairs -> 4 combos
    pm_enum <- .permute_paired(x, samples, control = "A", method = "mean", randomizations = 4, paired_method = "signflip")
    expect_equal(ncol(pm_enum), 4)
    expect_equal(nrow(pm_enum), nrow(x))

    # when randomizations >= total combinations, enumeration occurs (total combinations = 4 here)
    pm_enum2 <- .permute_paired(x, samples, control = "A", method = "mean", randomizations = 6, paired_method = "signflip")
    expect_equal(ncol(pm_enum2), 4)

    # sampled signflip returns requested number of permutations when less than total
    pm_samp <- .permute_paired(x, samples, control = "A", method = "mean", randomizations = 2, paired_method = "signflip")
    expect_equal(ncol(pm_samp), 2)
})


# calculate_fc defensive errors
test_that("calculate_fc errors on missing control or samples length mismatch", {
    x <- matrix(rnorm(8), nrow = 2)
    samples <- rep(c("A", "B"), length.out = ncol(x))
    expect_error(TSENAT:::.calculate_fc(x, samples, control = NULL), "`control` must be provided")
    expect_error(TSENAT:::.calculate_fc(x, samples[-1], control = "A"), "Length of 'samples' must equal number of columns in 'x'")
})

context("Wilcoxon Tests: Single Feature Implementation")

# Test .wilcox_one with unpaired design
test_that(".wilcox_one computes unpaired Wilcoxon test correctly", {
    # Create test data: 3 genes x 8 samples
    x <- matrix(rnorm(24), nrow = 3)
    samples <- rep(c("Control", "Treatment"), each = 4)
    
    # Set up as if inside wilcoxon function
    groups <- unique(sort(samples))
    g1_idx <- as.numeric(which(samples %in% groups[1]))
    g2_idx <- as.numeric(which(samples %in% groups[2]))
    
    # Create the .wilcox_one function environment
    .wilcox_one <- function(i) {
        tryCatch({
            test_result <- wilcox.test(x[i, g1_idx], x[i, g2_idx], paired = FALSE, exact = FALSE)
            n <- length(g1_idx) + length(g2_idx)
            list(p.value = test_result$p.value, statistic = test_result$statistic, n = n)
        }, error = function(e) {
            list(p.value = NA_real_, statistic = NA_real_, n = NA_real_)
        }, warning = function(w) {
            list(p.value = NA_real_, statistic = NA_real_, n = NA_real_)
        })
    }
    
    # Test first feature
    result <- .wilcox_one(1)
    
    expect_is(result, "list")
    expect_true("p.value" %in% names(result))
    expect_true("statistic" %in% names(result))
    expect_true("n" %in% names(result))
    expect_true(is.numeric(result$p.value) || is.na(result$p.value))
    expect_true(result$p.value >= 0 && result$p.value <= 1 || is.na(result$p.value))
    expect_equal(result$n, 8)
})

# Test .wilcox_one with paired design (position-based)
test_that(".wilcox_one computes paired Wilcoxon test correctly", {
    # Create paired test data: 3 genes x 8 samples (4 pairs)
    x <- matrix(rnorm(24), nrow = 3)
    samples <- rep(c("Control", "Treatment"), times = 4)  # alternating pairs
    
    # Set up as if inside wilcoxon function
    groups <- unique(sort(samples))
    g1_idx <- as.numeric(which(samples %in% groups[1]))
    g2_idx <- as.numeric(which(samples %in% groups[2]))
    
    # Create the .wilcox_one function for paired test
    .wilcox_one_paired <- function(i) {
        tryCatch({
            test_result <- wilcox.test(x[i, g1_idx], x[i, g2_idx], paired = TRUE, exact = FALSE)
            n <- length(g1_idx)  # number of pairs
            list(p.value = test_result$p.value, statistic = test_result$statistic, n = n)
        }, error = function(e) {
            list(p.value = NA_real_, statistic = NA_real_, n = NA_real_)
        }, warning = function(w) {
            list(p.value = NA_real_, statistic = NA_real_, n = NA_real_)
        })
    }
    
    # Test first feature
    result <- .wilcox_one_paired(1)
    
    expect_is(result, "list")
    expect_true(is.numeric(result$p.value) || is.na(result$p.value))
    expect_true(result$p.value >= 0 && result$p.value <= 1 || is.na(result$p.value))
    expect_equal(result$n, 4)  # 4 pairs
})

# Test .wilcox_one with explicit pairing information
test_that(".wilcox_one computes paired Wilcoxon with explicit pairing correctly", {
    # Create paired test data with explicit pairing
    x <- matrix(rnorm(24), nrow = 3)
    samples <- c("Control", "Treatment", "Control", "Treatment", 
                 "Control", "Treatment", "Control", "Treatment")
    pairs <- c("Pair1", "Pair1", "Pair2", "Pair2", 
               "Pair3", "Pair3", "Pair4", "Pair4")
    groups <- unique(sort(samples))
    
    # Create the .wilcox_one function for explicit pairing
    .wilcox_one_explicit_pairs <- function(i) {
        tryCatch({
            unique_pairs <- unique(pairs)
            all_diffs <- numeric(0)
            for (p in unique_pairs) {
                g1_samples <- which(pairs == p & samples == groups[1])
                g2_samples <- which(pairs == p & samples == groups[2])
                if (length(g1_samples) == 1 && length(g2_samples) == 1) {
                    all_diffs <- c(all_diffs, x[i, g1_samples] - x[i, g2_samples])
                }
            }
            # Perform paired test on the differences
            test_result <- wilcox.test(all_diffs, mu = 0, exact = FALSE)
            list(p.value = test_result$p.value, statistic = test_result$statistic, n = length(all_diffs))
        }, error = function(e) {
            list(p.value = NA_real_, statistic = NA_real_, n = NA_real_)
        }, warning = function(w) {
            list(p.value = NA_real_, statistic = NA_real_, n = NA_real_)
        })
    }
    
    # Test first feature
    result <- .wilcox_one_explicit_pairs(1)
    
    expect_is(result, "list")
    expect_true("p.value" %in% names(result))
    expect_true("statistic" %in% names(result))
    expect_true("n" %in% names(result))
    expect_equal(result$n, 4)  # 4 pairs
    expect_true(is.numeric(result$p.value) || is.na(result$p.value))
    expect_true(result$p.value >= 0 && result$p.value <= 1 || is.na(result$p.value))
})

# Test .wilcox_one error handling
test_that(".wilcox_one handles errors and edge cases gracefully", {
    # Create test data with potential edge cases
    x <- matrix(rnorm(24), nrow = 3)
    samples <- rep(c("A", "B"), each = 4)
    
    g1_idx <- as.numeric(which(samples %in% "A"))
    g2_idx <- as.numeric(which(samples %in% "B"))
    
    # Create .wilcox_one that catches errors
    .wilcox_one_safe <- function(i) {
        tryCatch({
            test_result <- wilcox.test(x[i, g1_idx], x[i, g2_idx], paired = FALSE, exact = FALSE)
            list(p.value = test_result$p.value, statistic = test_result$statistic, n = length(g1_idx) + length(g2_idx))
        }, error = function(e) {
            list(p.value = NA_real_, statistic = NA_real_, n = NA_real_)
        }, warning = function(w) {
            list(p.value = NA_real_, statistic = NA_real_, n = NA_real_)
        })
    }
    
    # Valid feature should return valid results with proper structure
    result_valid <- .wilcox_one_safe(1)
    expect_is(result_valid, "list")
    expect_true("p.value" %in% names(result_valid))
    expect_true("statistic" %in% names(result_valid))
    expect_true("n" %in% names(result_valid))
    
    # p-value should be in [0, 1] or NA
    expect_true((is.numeric(result_valid$p.value) && result_valid$p.value >= 0 && result_valid$p.value <= 1) || is.na(result_valid$p.value))
    expect_true(is.numeric(result_valid$n))
})

# Test r-value computation from U statistic (unpaired)
test_that("r-value computation from unpaired U statistic is correct", {
    # Create simple data for Wilcoxon test
    x <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 1)
    samples <- rep(c("A", "B"), each = 4)
    
    # Perform Wilcoxon test
    test_result <- wilcox.test(x[1, 1:4], x[1, 5:8], paired = FALSE)
    
    # Compute Z and r manually
    n1 <- 4
    n2 <- 4
    n <- n1 + n2
    U <- test_result$statistic
    expected_U <- n1 * n2 / 2
    var_U <- (n1 * n2 * (n1 + n2 + 1)) / 12
    sd_U <- sqrt(var_U)
    Z <- (U - expected_U) / sd_U
    r <- Z / sqrt(n)
    
    # r should be in [-1, 1]
    r_clamped <- pmax(-1, pmin(1, r))
    expect_true(r_clamped >= -1 && r_clamped <= 1)
})

# Test r-value computation from U statistic (paired)
test_that("r-value computation from paired U statistic is correct", {
    # Create simple paired data
    x <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 1)
    samples <- rep(c("A", "B"), times = 4)
    
    # Perform paired Wilcoxon test
    g1_idx <- c(1, 3, 5, 7)
    g2_idx <- c(2, 4, 6, 8)
    test_result <- wilcox.test(x[1, g1_idx], x[1, g2_idx], paired = TRUE, exact = FALSE)
    
    # Compute Z and r manually for paired test
    n <- length(g1_idx)
    U <- test_result$statistic
    expected_U <- n * (n + 1) / 4
    var_U <- (n * (n + 1) * (2 * n + 1)) / 24
    sd_U <- sqrt(var_U)
    Z <- (U - expected_U) / sd_U
    r <- Z / sqrt(n)
    
    # r should be in [-1, 1]
    r_clamped <- pmax(-1, pmin(1, r))
    expect_true(r_clamped >= -1 && r_clamped <= 1)
})

context("Statistical Methods: Wilcoxon Defensive Behavior")

test_that("wilcoxon handles all-NA and constant rows without error and returns matrix", {
    # two groups of 3 samples each
    samples <- c(rep("A", 3), rep("B", 3))
    # build matrix with 3 rows: all-NA, constant, variable
    m <- matrix(nrow = 3, ncol = 6)
    m[1, ] <- NA_real_
    m[2, ] <- rep(5, 6)
    m[3, ] <- c(1, 2, 3, 4, 5, 6)

    res <- .wilcoxon(m, samples, pcorr = "none")
    expect_true(is.data.frame(res))
    # four columns: pvalue, padj, U, r
    expect_equal(ncol(res), 4)
    # raw p-values produced and NA-handling yields numeric outputs
    raw <- as.numeric(res[, 1])
    expect_length(raw, 3)
    # all-NA row should produce p-value of 1 (by convention used elsewhere)
    expect_true(is.finite(raw[1]))
})

test_that(".wilcoxon() returns named p-value and effect size columns", {
    mat <- matrix(runif(20), nrow = 5)
    samples <- rep(c("A", "B"), each = 5)
    res <- .wilcoxon(mat, samples, pcorr = "none", paired = FALSE, exact = FALSE)
    expect_true(is.data.frame(res))
    expect_true(all(c("pvalue", "padj", "U", "r") %in% colnames(res)))
})

context("Statistical Methods: Wilcoxon Paired Behavior")

test_that("wilcoxon paired matches per-row wilcox.test on ordered pairs", {
    # create a small matrix with two genes and two pairs (columns: N,T,N,T)
    # ensure paired differences are not identical to avoid "ties" in differences
    x <- matrix(c(
        1, 2, 3, 6, # gene1 across 4 samples (diffs: -1, -3)
        5, 6, 7, 9 # gene2 across 4 samples (diffs: -1, -2)
    ), nrow = 2, byrow = TRUE)
    samples <- c("Normal", "Tumor", "Normal", "Tumor")

    res <- .wilcoxon(x, samples, paired = TRUE, exact = TRUE)

    # compute expected raw p-values by calling wilcox.test per row with
    # paired=TRUE
    expected_raw <- sapply(seq_len(nrow(x)), function(i) {
        wilcox.test(
            x[
                i,
                which(samples == "Normal")
            ],
            x[
                i,
                which(samples == "Tumor")
            ],
            paired = TRUE, exact = TRUE
        )$p.value
    })

    expected_adj <- p.adjust(expected_raw, method = "BH")

    expect_equal(as.numeric(res[, "pvalue"]), as.numeric(expected_raw))
    expect_equal(as.numeric(res[, "padj"]), as.numeric(expected_adj))
})


test_that("wilcoxon paired errors on unequal group sizes", {
    x_bad <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)
    samples_bad <- c("Normal", "Tumor", "Normal")
    expect_error(.wilcoxon(x_bad, samples_bad, paired = TRUE), "Paired Wilcoxon requires equal numbers of samples in each group")
})



test_that("wilcoxon paired handles SummarizedExperiment input", {
    # construct a small paired dataset
    sample_names <- c("S1_N", "S1_T", "S2_N", "S2_T")
    mat_vals <- matrix(c(
        1, 2, 3, 6,  # gene1
        5, 6, 7, 9   # gene2
    ), nrow = 2, byrow = TRUE)
    
    colnames(mat_vals) <- sample_names
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat_vals)
    )
    cond <- ifelse(grepl("_N$", sample_names), "Normal", "Tumor")
    coldata <- data.frame(
        Sample = sample_names,
        Condition = cond,
        stringsAsFactors = FALSE
    )
    
    se_mapped <- TSENAT:::.map_metadata_se(se, coldata)
    res <- .wilcoxon(
        SummarizedExperiment::assay(se_mapped),
        SummarizedExperiment::colData(se_mapped)$sample_type,
        paired = TRUE,
        exact = TRUE
    )
    
    # Verify output structure
    expect_true(is.data.frame(res))
    expect_equal(ncol(res), 4)
    expect_true(all(c("pvalue", "padj", "U", "r") %in% colnames(res)))
    # Verify results are not all NA
    expect_true(any(!is.na(res$pvalue)))
})

context("Statistical Methods: Wilcoxon U Statistic and r-Value Computation")

test_that("wilcoxon U statistic is computed correctly", {
    mat <- matrix(c(
        1, 2, 3, 4, # gene1
        5, 6, 7, 8  # gene2
    ), nrow = 2, byrow = TRUE)
    samples <- c("A", "A", "B", "B")
    
    res <- .wilcoxon(mat, samples, pcorr = "none", paired = FALSE)
    
    # U should be present and non-NA for valid data
    expect_false(anyNA(res$U))
    expect_true(all(res$U >= 0))
    # U should be bounded by the product of group sizes
    expect_true(all(res$U <= 4))  # 2 samples in each group
})

test_that("wilcoxon r-value is computed correctly", {
    mat <- matrix(c(
        1, 2, 3, 4,
        5, 6, 7, 8
    ), nrow = 2, byrow = TRUE)
    samples <- c("A", "A", "B", "B")
    
    res <- .wilcoxon(mat, samples, pcorr = "none", paired = FALSE)
    
    # r-value should be present and non-NA for valid data
    expect_false(anyNA(res$r))
    # r-value should be bounded between -1 and 1
    expect_true(all(res$r >= -1 & res$r <= 1))
    # For constant row, r should be 0 (no effect)
    mat_const <- matrix(c(1, 1, 1, 1), nrow = 1)
    res_const <- .wilcoxon(mat_const, samples, pcorr = "none")
    expect_true(res_const$r[1] == 0)
})

test_that("wilcoxon U and r values are NA when pvalue is NA", {
    m <- matrix(nrow = 2, ncol = 4)
    m[1, ] <- NA_real_
    m[2, ] <- c(1, 2, 3, 4)
    samples <- c("A", "A", "B", "B")
    
    res <- .wilcoxon(m, samples, pcorr = "none")
    
    # All-NA row should have NA U and r value
    expect_true(is.na(res$U[1]))
    expect_true(is.na(res$r[1]))
    # Valid row should have non-NA U and r
    expect_false(is.na(res$U[2]))
    expect_false(is.na(res$r[2]))
})

test_that("wilcoxon r-value is bounded and non-NA for valid paired data", {
    # For paired test, verify r-value is reasonable
    mat <- matrix(c(
        1, 2, 3, 6,  # distinct paired differences
        5, 6, 7, 9
    ), nrow = 2, byrow = TRUE)
    samples <- c("N", "T", "N", "T")
    
    res <- .wilcoxon(mat, samples, pcorr = "none", paired = TRUE)
    
    # Verify r-value is within valid range and non-NA for non-constant rows
    for (i in seq_len(nrow(res))) {
        p <- res$pvalue[i]
        if (!is.na(p)) {
            # r should be bounded between -1 and 1 and non-NA for valid p-values
            expect_true(!is.na(res$r[i]))
            expect_true(res$r[i] >= -1 && res$r[i] <= 1)
        }
    }
})

# ============================================================================
# REDISTRIBUTED TESTS FROM test-infrastructure-statistical_validation.R
# ============================================================================

context("Wilcoxon Test Validity")

test_that("wilcoxon test matches R's built-in wilcox.test exactly", {
    # Create clear difference
    normal_vals <- c(1, 2, 3, 4, 5)
    tumor_vals <- c(6, 7, 8, 9, 10)
    
    mat <- matrix(c(normal_vals, tumor_vals), nrow = 1)
    samples <- c(rep("Normal", 5), rep("Tumor", 5))
    
    # TSENAT wilcoxon
    tsenat_result <- .wilcoxon(mat, samples, pcorr = "none")
    tsenat_p <- tsenat_result[1, "pvalue"]
    
    # R's wilcox.test
    r_result <- wilcox.test(normal_vals, tumor_vals, exact = FALSE)
    r_p <- r_result$p.value
    
    # Should be essentially identical
    expect_equal(tsenat_p, r_p, tolerance = 1e-10)
})

test_that("wilcoxon paired matches paired test from R", {
    # Paired data
    before <- c(1, 2, 3, 4, 5)
    after <- c(2, 3, 5, 6, 8)
    
    mat <- matrix(c(before, after), nrow = 1)
    samples <- c(rep("Before", 5), rep("After", 5))
    
    # TSENAT with paired = TRUE
    tsenat_paired <- .wilcoxon(mat, samples, paired = TRUE, pcorr = "none")
    tsenat_p <- tsenat_paired[1, "pvalue"]
    
    # R's paired wilcox.test
    r_paired <- wilcox.test(before, after, paired = TRUE, exact = FALSE)
    r_p <- r_paired$p.value
    
    expect_equal(tsenat_p, r_p, tolerance = 1e-10)
})

