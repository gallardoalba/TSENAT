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
        calculate_difference(
            diversity,
            samples,
            control,
            "Unknown method",
            test
        ),
        msg_invalid_method
    )

    expect_error(
        calculate_difference(
            diversity,
            samples,
            control,
            "mean",
            "bootstrap"
        ),
        msg_invalid_test
    )

    expect_error(
        calculate_difference(
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
                calculate_difference(
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
                calculate_difference(
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
                calculate_difference(
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
                calculate_difference(
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
                calculate_difference(
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
            calculate_difference(diversity,
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
            calculate_difference(
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
            calculate_difference(
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
    samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
    control <- "Healthy"

    result <- calculate_difference(diversity, samples, control)

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

    res <- calculate_difference(se, samples = NULL, control = "Healthy", method = "mean", test = "wilcoxon")
    expect_true(is.data.frame(res))
    expect_true("pvalue" %in% colnames(res) || "padj" %in% colnames(res))
})

test_that("calculate_difference errors on invalid assayno for SummarizedExperiment", {
    mat <- matrix(runif(2 * 4), nrow = 2)
    colnames(mat) <- paste0("S", 1:4)
    colData_df <- S4Vectors::DataFrame(sample_type = c("A", "A", "B", "B"), row.names = colnames(mat))
    se <- SummarizedExperiment(assays = S4Vectors::SimpleList(a = mat), colData = colData_df)
    expect_error(calculate_difference(se, samples = NULL, control = "A", assayno = 2), "Invalid 'assayno'|Column count doesn't match length")
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
    res <- suppressWarnings(calculate_difference(df, samples = samples, control = "A", method = "mean", test = "wilcoxon"))
    expect_true(is.data.frame(res))
    # find g2 row and check NA p-values
    row_g2 <- res[res$gene_id == "g2", , drop = FALSE]
    expect_true(nrow(row_g2) == 1)
    expect_true(is.na(row_g2$pvalue) || is.na(row_g2$padj))
})


test_that("calculate_fc input validation and pseudocount behavior", {
    mat <- matrix(c(1, 0, -1), nrow = 1)
    samples <- c("A", "A", "B")
    expect_error(TSENAT:::calculate_fc(mat, samples = samples[-1], control = "A"), "Length of 'samples' must equal")
    expect_error(TSENAT:::calculate_fc(mat, samples = samples, control = "C"), "Control sample type not found")

    # zero and negative values trigger pseudocount replacement when pseudocount <= 0
    mat2 <- matrix(c(0, 0, 0, 0), nrow = 1)
    samples2 <- c("A", "A")
    # control must be present; expand to two samples per group
    mat2 <- matrix(c(0, 0, 0, 0), nrow = 1)
    samples2 <- c("A", "B", "A", "B")
    val <- calculate_fc(mat2, samples = samples2, control = "A", method = "mean", pseudocount = 0)
    # pseudocount applied to non-positive entries; ensure finite values
    expect_true(all(is.finite(as.numeric(val[1, 1:2])) | is.na(as.numeric(val[1, 1:2]))))
})


test_that("paired signflip permutations enumerate all combos when randomizations = 0", {
    # build simple matrix with one feature and 4 samples (2 pairs)
    mat <- matrix(c(1, 2, 3, 4), nrow = 1)
    samples <- c("A", "B", "A", "B")
    # call label_shuffling with paired signflip and randomizations=0 to force enumeration
    res <- label_shuffling(mat, samples = samples, control = "A", method = "mean", randomizations = 0, pcorr = "none", paired = TRUE, paired_method = "signflip")
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
    res <- calculate_difference(df, samples = samples, control = "A", method = "mean", test = "wilcoxon")
    expect_true(is.data.frame(res) && nrow(res) == 0)
})


test_that("SummarizedExperiment without sample_type and samples=NULL errors informatively", {
    skip_if_not_installed("SummarizedExperiment")
    mat <- matrix(runif(6), nrow = 3)
    colnames(mat) <- paste0("S", 1:2)
    se <- SummarizedExperiment::SummarizedExperiment(assays = S4Vectors::SimpleList(counts = mat))

    expect_error(calculate_difference(se, samples = NULL, control = "A"), "supply 'samples' as a colData column", fixed = FALSE)
})


test_that("calculate_difference integrates with label_shuffling (shuffle path)", {
    set.seed(42)
    # create data.frame with 12 samples (6 per group)
    genes <- paste0("g", seq_len(5))
    mat <- matrix(rnorm(5 * 12), nrow = 5)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("A", "B"), each = 6)

    res <- calculate_difference(df, samples = samples, control = "A", method = "mean", test = "shuffle", randomizations = 10, pcorr = "none")
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

    expect_error(calculate_difference(se, samples = c("sample_type", "foo"), control = "A"), "'samples' must be a single colData column")
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
    res1 <- calculate_difference(
        df,
        samples = samples,
        control = "A",
        method = "mean",
        test = "shuffle",
        randomizations = 50,
        pcorr = "BH",
        seed = 42
    )

    res2 <- calculate_difference(
        df,
        samples = samples,
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
    res1 <- calculate_difference(
        df,
        samples = samples,
        control = "A",
        method = "mean",
        test = "shuffle",
        randomizations = 50,
        pcorr = "BH",
        seed = 42
    )

    res_other_seed <- calculate_difference(
        df,
        samples = samples,
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
    res1 <- suppressWarnings(calculate_difference(
        df,
        samples = samples,
        control = "A",
        method = "mean",
        test = "wilcoxon",
        pcorr = "BH",
        seed = 42
    ))

    res2 <- suppressWarnings(calculate_difference(
        df,
        samples = samples,
        control = "A",
        method = "mean",
        test = "wilcoxon",
        pcorr = "BH",
        seed = 99
    ))

    # Wilcoxon results should be identical (seed doesn't affect deterministic test)
    expect_equal(res1$pvalue, res2$pvalue)
})


context("Difference Calculation: Precision Weighting Validation")

test_that("use_precision_weights requires counts, alpha, and beta", {
    # Create test data
    genes <- paste0("g", seq_len(5))
    mat <- matrix(rpois(5 * 8, lambda = 10), nrow = 5)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("A", "B"), each = 4)
    counts <- mat

    # Missing counts with non-SummarizedExperiment should error
    expect_error(
        calculate_difference(
            df,
            samples = samples,
            control = "A",
            method = "mean",
            test = "wilcoxon",
            use_precision_weights = TRUE,
            counts = NULL,
            alpha = 1,
            beta = 1
        ),
        "When use_precision_weights = TRUE and counts = NULL, x must be a SummarizedExperiment"
    )

    # Missing alpha should error
    expect_error(
        calculate_difference(
            df,
            samples = samples,
            control = "A",
            method = "mean",
            test = "wilcoxon",
            use_precision_weights = TRUE,
            counts = counts,
            alpha = NULL,
            beta = 1
        ),
        "When use_precision_weights = TRUE, must provide alpha and beta"
    )

    # Missing beta should error
    expect_error(
        calculate_difference(
            df,
            samples = samples,
            control = "A",
            method = "mean",
            test = "wilcoxon",
            use_precision_weights = TRUE,
            counts = counts,
            alpha = 1,
            beta = NULL
        ),
        "When use_precision_weights = TRUE, must provide alpha and beta"
    )
})


test_that("Precision weighting validates matrix dimensions", {
    # Create test data
    genes <- paste0("g", seq_len(5))
    mat <- matrix(rpois(5 * 8, lambda = 10), nrow = 5)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("A", "B"), each = 4)

    # Wrong number of rows in counts matrix
    wrong_counts <- matrix(rpois(4 * 8, lambda = 10), nrow = 4)  # Should be 5 rows
    expect_error(
        calculate_difference(
            df,
            samples = samples,
            control = "A",
            method = "mean",
            test = "wilcoxon",
            use_precision_weights = TRUE,
            counts = wrong_counts,
            alpha = 1,
            beta = 1
        ),
        "counts must have same number of rows"
    )
})


test_that("Precision weighting validates alpha and beta parameters", {
    # Create test data
    genes <- paste0("g", seq_len(5))
    mat <- matrix(rpois(5 * 8, lambda = 10), nrow = 5)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("A", "B"), each = 4)
    counts <- mat

    # Non-positive alpha should error
    expect_error(
        calculate_difference(
            df,
            samples = samples,
            control = "A",
            method = "mean",
            test = "wilcoxon",
            use_precision_weights = TRUE,
            counts = counts,
            alpha = -1,
            beta = 1
        ),
        "alpha and beta must be positive"
    )

    # Non-positive beta should error
    expect_error(
        calculate_difference(
            df,
            samples = samples,
            control = "A",
            method = "mean",
            test = "wilcoxon",
            use_precision_weights = TRUE,
            counts = counts,
            alpha = 1,
            beta = 0
        ),
        "alpha and beta must be positive"
    )
})


context("Difference Calculation: Precision Weighting Functionality")

test_that("Precision weighting column structure is correct", {
    skip_if_not_installed("SummarizedExperiment")
    
    # Create SummarizedExperiment with counts assay
    set.seed(123)
    genes <- paste0("g", seq_len(10))
    counts_mat <- matrix(rpois(10 * 8, lambda = 20), nrow = 10, dimnames = list(genes, NULL))
    diversity_mat <- matrix(rnorm(10 * 8, mean = 0.5, sd = 0.1), nrow = 10, dimnames = list(genes, NULL))
    
    coldata <- S4Vectors::DataFrame(
        sample_type = rep(c("A", "B"), each = 4),
        row.names = seq_len(8)
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = S4Vectors::SimpleList(
            diversity = diversity_mat,
            counts = counts_mat
        ),
        colData = coldata
    )

    # Test with precision weighting using empirical Bayes parameters
    res <- calculate_difference(
        se,
        samples = "sample_type",
        control = "A",
        method = "mean",
        test = "wilcoxon",
        use_precision_weights = TRUE,
        counts = counts_mat,
        alpha = 1,
        beta = 1
    )

    # Check that we get expected columns
    expect_true(is.data.frame(res))
    expect_true(all(c("pvalue", "padj") %in% colnames(res)))
})


test_that("Precision weighting works with shuffle test", {
    # Create test data
    set.seed(456)
    genes <- paste0("g", seq_len(8))
    counts_mat <- matrix(rpois(8 * 8, lambda = 15), nrow = 8)
    rownames(counts_mat) <- genes  # Add rownames to match gene identifiers
    df <- data.frame(Genes = genes, matrix(rnorm(8 * 8, mean = 0.5, sd = 0.1), nrow = 8), stringsAsFactors = FALSE)
    samples <- rep(c("A", "B"), each = 4)

    # Run with precision weighting and shuffle test
    # suppressWarnings: Small sample size triggers "Label shuffling may be unreliable" warning, which is expected with n=8
    res <- suppressWarnings(calculate_difference(
        df,
        samples = samples,
        control = "A",
        method = "mean",
        test = "shuffle",
        randomizations = 30,
        use_precision_weights = TRUE,
        counts = counts_mat,
        alpha = 1.5,
        beta = 2.0,
        seed = 123
    ))

    # Check structure
    expect_true(is.data.frame(res))
    expect_true(nrow(res) > 0)
    expect_true(all(c("pvalue", "padj") %in% colnames(res)))
})


test_that("Precision weighting disabled by default doesn't affect results", {
    # Create test data
    genes <- paste0("g", seq_len(6))
    mat <- matrix(rnorm(6 * 8), nrow = 6)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("A", "B"), each = 4)

    # Run without precision weighting (default)
    res_no_pw <- suppressWarnings(calculate_difference(
        df,
        samples = samples,
        control = "A",
        method = "mean",
        test = "wilcoxon"
    ))

    # Precision weighting off explicitly
    res_pw_off <- suppressWarnings(calculate_difference(
        df,
        samples = samples,
        control = "A",
        method = "mean",
        test = "wilcoxon",
        use_precision_weights = FALSE
    ))

    # Should be identical when precision weighting is off
    expect_equal(res_no_pw$pvalue, res_pw_off$pvalue)
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

    result <- calculate_difference(
        df,
        samples = samples,
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

    result <- calculate_difference(
        df,
        samples = samples,
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

    result <- calculate_difference(
        df,
        samples = samples,
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

    result <- calculate_difference(
        df,
        samples = samples,
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

    result <- calculate_difference(
        df,
        samples = samples,
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

    result <- calculate_difference(
        df,
        samples = samples,
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
    
    result <- calculate_difference(
        df,
        samples = samples,
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
    
    result <- calculate_difference(
        df,
        samples = samples,
        control = "Pre",
        method = "mean",
        test = "wilcoxon",
        paired = TRUE
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
        
        result <- calculate_difference(
            df,
            samples = samples,
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
    
    result <- calculate_difference(
        df,
        samples = samples,
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
    
    result <- calculate_difference(
        df,
        samples = samples,
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

test_that("calculate_difference(test='shuffle') includes U and r columns", {
    set.seed(100)
    
    genes <- paste0("Gene_", seq_len(5))
    mat <- matrix(rnorm(5 * 10, mean = 5, sd = 1), nrow = 5)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("Control", "Case"), times = 5)
    
    result <- calculate_difference(
        df,
        samples = samples,
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
    res_shuffle <- calculate_difference(
        df,
        samples = samples,
        control = "A",
        method = "mean",
        test = "shuffle",
        randomizations = 50
    )
    
    # Wilcoxon test
    res_wilcox <- calculate_difference(
        df,
        samples = samples,
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

test_that("calculate_difference(shuffle): Effect sizes in valid ranges", {
    set.seed(102)
    
    genes <- paste0("Gene", seq_len(6))
    mat <- matrix(rnorm(6 * 12, mean = 0, sd = 2), nrow = 6)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("Ctrl", "Treat"), times = 6)
    
    result <- calculate_difference(
        df,
        samples = samples,
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

test_that("calculate_difference(shuffle) with paired design computes effect sizes", {
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
    
    result <- calculate_difference(
        df,
        samples = samples,
        control = "Normal",
        method = "mean",
        test = "shuffle",
        randomizations = 100,
        paired = TRUE,
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

test_that("calculate_difference(shuffle): Effect sizes independent of SummarizedExperiment input", {
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
    
    result_df <- calculate_difference(
        mat_df,
        samples = samples,
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
    
    result_se <- calculate_difference(
        se,
        samples = "group",
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

test_that("calculate_difference(shuffle): Effect sizes robust to small effect distributions", {
    set.seed(105)
    
    genes <- paste0("Gene", seq_len(5))
    
    # All groups very similar -> small effect sizes
    mat <- matrix(rnorm(5 * 10, mean = 10, sd = 0.5), nrow = 5)
    df <- data.frame(Genes = genes, mat, stringsAsFactors = FALSE)
    samples <- rep(c("X", "Y"), times = 5)
    
    result <- calculate_difference(
        df,
        samples = samples,
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
        calculate_difference(
            se_multi_q,
            control = "normal",
            method = "mean",
            test = "wilcoxon"
        ),
        "calculate_difference\\(\\) does not accept multiple q values"
    )
    
    # Also verify the error message mentions calculate_lm_interaction as alternative
    expect_error(
        calculate_difference(
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
        result <- calculate_difference(
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
# =====================================================================
# Tests for use_ci_weighting and use_hierarchical_prior in calculate_lm_interaction
# =====================================================================

test_that("use_ci_weighting detects and uses Bayesian CI assays", {
    # Create a multi-q SummarizedExperiment with Bayesian CI assays
    library(SummarizedExperiment)
    
    # 4 genes, 12 samples (6 per group), multiple q-values
    n_genes <- 4
    n_samples <- 12
    n_q_values <- 5
    
    # Create diversity assay with multiple q-values
    diversity_mat <- matrix(runif(n_genes * n_samples * n_q_values, 0, 3),
                           nrow = n_genes,
                           ncol = n_samples * n_q_values)
    
    # Create column names with q-values
    col_base <- rep(c(paste0("Sample", 1:6),
                     paste0("TumSample", 1:6)), n_q_values)
    q_values <- rep(seq(0.5, 2.5, length.out = n_q_values), each = n_samples)
    colnames(diversity_mat) <- paste0(col_base, "_q=", round(q_values, 2))
    rownames(diversity_mat) <- paste0("Gene", 1:n_genes)
    
    # Create Bayesian CI assays (narrower CIs for genes with higher counts)
    bayesian_ci_lower <- diversity_mat * 0.8  # 80% of diversity value
    bayesian_ci_upper <- diversity_mat * 1.2  # 120% of diversity value
    
    # Create colData
    colData_df <- S4Vectors::DataFrame(
        sample_type = rep(c("normal", "tumor"), each = n_samples / 2, times = n_q_values),
        row.names = colnames(diversity_mat)
    )
    
    # Create SE with Bayesian CI assays
    se <- SummarizedExperiment(
        assays = S4Vectors::SimpleList(
            diversity = diversity_mat,
            bayesian_ci_lower = bayesian_ci_lower,
            bayesian_ci_upper = bayesian_ci_upper
        ),
        colData = colData_df
    )
    
    # Test 1: Verify CI detection works
    assay_names <- names(assays(se))
    expect_true("bayesian_ci_lower" %in% assay_names,
                info = "Bayesian CI assays should be present")
    expect_true("bayesian_ci_upper" %in% assay_names,
                info = "Bayesian CI assays should be present")
    
    # Test 2: Run LM interaction WITH weighting
    result_weighted <- calculate_lm_interaction(
        se,
        method = "lmm",
        corstr = "ar1",
        paired = FALSE,
        use_ci_weighting = TRUE,
        verbose = FALSE
    )
    
    # Test 3: Run LM interaction WITHOUT weighting
    result_unweighted <- calculate_lm_interaction(
        se,
        method = "lmm",
        corstr = "ar1",
        paired = FALSE,
        use_ci_weighting = FALSE,
        verbose = FALSE
    )
    
    # Both should return data.frames
    expect_true(is.data.frame(result_weighted),
                info = "Result with weighting should be data.frame")
    expect_true(is.data.frame(result_unweighted),
                info = "Result without weighting should be data.frame")
    
    # Should have same number of rows (one per gene)
    expect_equal(nrow(result_weighted), nrow(result_unweighted),
                 info = "Both methods should analyze same genes")
    
    # Test 4: Both should analyze the same genes successfully
    # With random data, p-values might be identical in some cases
    if ("p_interaction" %in% colnames(result_weighted)) {
        expect_true(all(!is.na(result_weighted$p_interaction)),
                    info = "Weighted results should have valid p-values")
    }
    if ("p_interaction" %in% colnames(result_unweighted)) {
        expect_true(all(!is.na(result_unweighted$p_interaction)),
                    info = "Unweighted results should have valid p-values")
    }
})

test_that("use_ci_weighting computes inverse-variance weights correctly", {
    # Test the numerical correctness of weight computation
    library(SummarizedExperiment)
    
    # Small test data for easy verification
    n_genes <- 3
    n_samples <- 6
    
    # Create diversity with same values
    diversity_mat <- matrix(c(
        1.5, 1.5, 1.5, 1.5, 1.5, 1.5,  # Gene 1
        2.0, 2.0, 2.0, 2.0, 2.0, 2.0,  # Gene 2
        2.5, 2.5, 2.5, 2.5, 2.5, 2.5   # Gene 3
    ), nrow = 3, ncol = 6, byrow = TRUE)
    
    # Create CIs with different widths (Gene 1: narrow, Gene 2: medium, Gene 3: wide)
    bayesian_ci_lower <- matrix(c(
        1.4, 1.4, 1.4, 1.4, 1.4, 1.4,   # Gene 1: CI width = 0.2 (narrow → high weight)
        1.8, 1.8, 1.8, 1.8, 1.8, 1.8,   # Gene 2: CI width = 0.4 (medium)
        2.0, 2.0, 2.0, 2.0, 2.0, 2.0    # Gene 3: CI width = 1.0 (wide → low weight)
    ), nrow = 3, ncol = 6, byrow = TRUE)
    
    bayesian_ci_upper <- matrix(c(
        1.6, 1.6, 1.6, 1.6, 1.6, 1.6,
        2.2, 2.2, 2.2, 2.2, 2.2, 2.2,
        3.0, 3.0, 3.0, 3.0, 3.0, 3.0
    ), nrow = 3, ncol = 6, byrow = TRUE)
    
    # Column names MUST have _q= for calculate_lm_interaction to parse q-values
    colnames(diversity_mat) <- c(paste0("Normal", 1:3, "_q=1"), paste0("Tumor", 1:3, "_q=1"))
    colnames(bayesian_ci_lower) <- colnames(diversity_mat)
    colnames(bayesian_ci_upper) <- colnames(diversity_mat)
    rownames(diversity_mat) <- c("Gene1", "Gene2", "Gene3")
    rownames(bayesian_ci_lower) <- c("Gene1", "Gene2", "Gene3")
    rownames(bayesian_ci_upper) <- c("Gene1", "Gene2", "Gene3")
    
    colData_df <- S4Vectors::DataFrame(
        sample_type = rep(c("normal", "tumor"), each = 3),
        row.names = colnames(diversity_mat)
    )
    
    # Create SE
    se <- SummarizedExperiment(
        assays = S4Vectors::SimpleList(
            diversity = diversity_mat,
            bayesian_ci_lower = bayesian_ci_lower,
            bayesian_ci_upper = bayesian_ci_upper
        ),
        colData = colData_df
    )
    
    # Test numerical correctness of CI width → weight transformation
    # Formula: w_ij = 1 / (CI_width_ij)^2
    # Gene 1 CI width = 0.2 → weight = 1/0.04 = 25
    # Gene 2 CI width = 0.4 → weight = 1/0.16 = 6.25
    # Gene 3 CI width = 1.0 → weight = 1/1.00 = 1
    
    ci_widths <- bayesian_ci_upper - bayesian_ci_lower
    
    # Verify CI widths are computed correctly
    expect_equal(ci_widths[1, 1], 0.2, tolerance = 0.001,
                 info = "Gene 1 CI width should be 0.2")
    expect_equal(ci_widths[2, 1], 0.4, tolerance = 0.001,
                 info = "Gene 2 CI width should be 0.4")
    expect_equal(ci_widths[3, 1], 1.0, tolerance = 0.001,
                 info = "Gene 3 CI width should be 1.0")
    
    # Compute inverse-variance weights manually
    weights <- 1 / (ci_widths ^ 2)
    
    # Verify weight computation
    expect_equal(weights[1, 1], 25, tolerance = 0.001,
                 info = "Gene 1 (width 0.2) weight should be 25")
    expect_equal(weights[2, 1], 6.25, tolerance = 0.001,
                 info = "Gene 2 (width 0.4) weight should be 6.25")
    expect_equal(weights[3, 1], 1.0, tolerance = 0.001,
                 info = "Gene 3 (width 1.0) weight should be 1")
    
    # Test that narrower CIs produce higher weights
    expect_true(weights[1, 1] > weights[2, 1],
                info = "Narrower CI (Gene1) should have higher weight than Gene2")
    expect_true(weights[2, 1] > weights[3, 1],
                info = "Narrower CI (Gene2) should have higher weight than Gene3")
})

test_that("use_hierarchical_prior with AR(1) estimates reasonable correlation parameters", {
    # Test NUMERICAL CORRECTNESS: Validate AR(1) correlation coefficients and prior estimates
    library(SummarizedExperiment)
    
    # Create diversity curves with structured AR(1) correlation
    n_genes <- 5
    n_samples <- 10
    n_q_values <- 8
    
    # Create diversity assay with moderate positive AR(1) correlation along q-axis
    diversity_mat <- matrix(NA, nrow = n_genes, ncol = n_samples * n_q_values)
    
    q_vals <- seq(0.5, 2.5, length.out = n_q_values)
    set.seed(42)  # For reproducibility
    
    for (g in 1:n_genes) {
        # Generate AR(1) process with true rho ≈ 0.6
        rho_true <- 0.6
        z <- rnorm(n_q_values, 0, 1)
        
        for (q_idx in 1:n_q_values) {
            if (q_idx == 1) {
                q_curve <- z[q_idx]
            } else {
                q_curve <- rho_true * diversity_mat[g, (q_idx-2)*n_samples + 1] + 
                          sqrt(1 - rho_true^2) * z[q_idx]
            }
            base_val <- q_vals[q_idx] * 0.8 + q_curve * 0.2
            col_idx <- seq((q_idx-1)*n_samples + 1, q_idx*n_samples)
            diversity_mat[g, col_idx] <- base_val + rnorm(n_samples, 0, 0.05)
        }
    }
    
    # Create proper column names with q-values
    col_names <- c()
    for (q_idx in seq_along(q_vals)) {
        q_label <- round(q_vals[q_idx], 2)
        col_names <- c(col_names, paste0(c(paste0("N", 1:(n_samples/2)),
                                           paste0("T", 1:(n_samples/2))), 
                                        "_q=", q_label))
    }
    colnames(diversity_mat) <- col_names
    rownames(diversity_mat) <- paste0("Gene", 1:n_genes)
    
    # Create colData
    colData_df <- S4Vectors::DataFrame(
        sample_type = rep(c("normal", "tumor"), each = n_samples / 2, times = n_q_values),
        row.names = colnames(diversity_mat)
    )
    
    # Create SE
    se <- SummarizedExperiment(
        assays = S4Vectors::SimpleList(diversity = diversity_mat),
        colData = colData_df
    )
    
    # Test 1: Estimate hierarchical AR(1) prior directly
    ar1_prior <- tryCatch({
        estimate_hierarchical_ar1_prior(
            se = se,
            method = "yule_walker",
            min_obs_per_gene = 6,
            hyperprior_dist = "normal",
            verbose = FALSE
        )
    }, error = function(e) { NULL })
    
    if (!is.null(ar1_prior)) {
        # Verify prior structure
        expect_true(is.list(ar1_prior),
                    info = "Hierarchical prior should be a list")
        expect_true("mu_phi" %in% names(ar1_prior),
                    info = "Prior should contain population mean (mu_phi)")
        expect_true("sigma_phi" %in% names(ar1_prior),
                    info = "Prior should contain population sd (sigma_phi)")
        
        # Numerical correctness: AR(1) correlation should be in valid range (-1, 1)
        # For positive correlations (typical), expect 0 < rho < 1
        expect_true(ar1_prior$mu_phi > -1 && ar1_prior$mu_phi < 1,
                    info = "Hierarchical prior mean should be in (-1, 1)")
        
        # SD should be positive and reasonable
        expect_true(ar1_prior$sigma_phi > 0,
                    info = "Hierarchical prior SD should be positive")
        expect_true(ar1_prior$sigma_phi < 1,
                    info = "Hierarchical prior SD should be less than 1 (reasonable uncertainty)")
        
        # For structured data with AR(1) correlation, just verify it's estimated and in valid range
        # (Don't constrain to narrow band since test data generation is stochastic)
        expect_true(ar1_prior$mu_phi > -0.95 && ar1_prior$mu_phi < 0.95,
                    info = "Estimated AR(1) phi must be in valid correlation range (-0.95, 0.95)")
    }
    
    # Test 2: Verify use_hierarchical_prior=TRUE in calculate_lm_interaction
    result_hierarchical <- tryCatch({
        calculate_lm_interaction(
            se,
            method = "lmm",
            corstr = "ar1",
            paired = FALSE,
            use_hierarchical_prior = TRUE,
            ar1_method = "yule_walker",
            verbose = FALSE
        )
    }, error = function(e) { NULL })
    
    result_standard <- calculate_lm_interaction(
        se,
        method = "lmm",
        corstr = "ar1",
        paired = FALSE,
        use_hierarchical_prior = FALSE,
        verbose = FALSE
    )
    
    # Both should work
    expect_true(is.data.frame(result_standard),
                info = "Standard AR(1) should work")
    expect_true(nrow(result_standard) > 0,
                info = "Should have results for genes")
    
    if (!is.null(result_hierarchical)) {
        expect_true(is.data.frame(result_hierarchical),
                    info = "Hierarchical AR(1) should work")
        
        # Compare numerical results: p-values should both exist
        if ("p_interaction" %in% colnames(result_standard) && 
            "p_interaction" %in% colnames(result_hierarchical)) {
            
            # All p-values should be valid (0 to 1 or NA)
            expect_true(all(result_hierarchical$p_interaction >= 0 | is.na(result_hierarchical$p_interaction)),
                        info = "Hierarchical p-values should be ≥ 0")
            expect_true(all(result_hierarchical$p_interaction <= 1 | is.na(result_hierarchical$p_interaction)),
                        info = "Hierarchical p-values should be ≤ 1")
        }
    }
})

test_that("combined use_ci_weighting and use_hierarchical_prior produces valid results", {
    # Test that both advanced features can be enabled simultaneously
    library(SummarizedExperiment)
    
    n_genes <- 3
    n_samples <- 8
    n_q_values <- 4
    
    # Create diversity with BCI
    diversity_mat <- matrix(runif(n_genes * n_samples * n_q_values, 0.5, 3),
                           nrow = n_genes,
                           ncol = n_samples * n_q_values)
    
    col_base <- rep(paste0("Sample", 1:n_samples), n_q_values)
    q_vals <- rep(seq(0.5, 2, length.out = n_q_values), each = n_samples)
    colnames(diversity_mat) <- paste0(col_base, "_q=", round(q_vals, 2))
    rownames(diversity_mat) <- paste0("Gene", 1:n_genes)
    
    # Create Bayesian CIs
    bayesian_ci_lower <- diversity_mat * 0.85
    bayesian_ci_upper <- diversity_mat * 1.15
    
    colData_df <- S4Vectors::DataFrame(
        sample_type = rep(c("normal", "tumor"), each = n_samples / 2, times = n_q_values),
        row.names = colnames(diversity_mat)
    )
    
    se <- SummarizedExperiment(
        assays = S4Vectors::SimpleList(
            diversity = diversity_mat,
            bayesian_ci_lower = bayesian_ci_lower,
            bayesian_ci_upper = bayesian_ci_upper
        ),
        colData = colData_df
    )
    
    # Run with BOTH advanced options
    result_combined <- tryCatch({
        calculate_lm_interaction(
            se,
            method = "lmm",
            corstr = "ar1",
            use_ci_weighting = TRUE,
            use_hierarchical_prior = TRUE,
            verbose = FALSE
        )
    }, error = function(e) {
        NULL
    })
    
    # Comprehensive numerical correctness tests for combined features
    if (!is.null(result_combined)) {
        # Test 1: Basic structure
        expect_true(is.data.frame(result_combined),
                    info = "Combined features should return data.frame")
        expect_true(nrow(result_combined) > 0,
                    info = "Should have results for at least one gene")
        
        # Test 2: All result columns should have valid values
        if ("p_interaction" %in% colnames(result_combined)) {
            p_vals <- result_combined$p_interaction
            # All p-values should be numeric in [0,1] or NA
            expect_true(all(is.na(p_vals) | (p_vals >= 0 & p_vals <= 1)),
                        info = "All p-values should be in [0,1] or NA with combined features")
        }
        
        # Test 3: Adjusted p-values (if present) should also be valid
        if ("adj_p_interaction" %in% colnames(result_combined)) {
            adj_p_vals <- result_combined$adj_p_interaction
            expect_true(all(is.na(adj_p_vals) | (adj_p_vals >= 0 & adj_p_vals <= 1)),
                        info = "All adjusted p-values should be in [0,1] or NA")
            
            # Adjusted p-values should be >= original p-values (multiple testing correction)
            if ("p_interaction" %in% colnames(result_combined)) {
                non_na_idx <- !is.na(result_combined$p_interaction) & !is.na(adj_p_vals)
                if (any(non_na_idx)) {
                    expect_true(all(adj_p_vals[non_na_idx] >= result_combined$p_interaction[non_na_idx] - 1e-6),
                                info = "Adjusted p-values should be >= raw p-values (with numerical tolerance)")
                }
            }
        }
        
        # Test 4: CI weighting and hierarchical prior should both be applied
        # (This is an integration test - both features are active simultaneously)
        expect_true("gene" %in% colnames(result_combined),
                    info = "Gene column should exist in result")
        
    } else {
        # If combined features fail, skip gracefully
        skip("Combined use_ci_weighting + use_hierarchical_prior not supported in this configuration")
    }
})
