context("Difference calculation")

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
    expect_length(result, 7)
    expect_equal(mean(result$Pathogenic_mean), 0.65, tolerance = 0.001, scale = 1)
    expect_equal(mean(result$Healthy_mean), 0.25, tolerance = 0.001, scale = 1)
})

context("calculate_difference additional tests")

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
    expect_true("raw_p_values" %in% colnames(res) || "adjusted_p_values" %in% colnames(res))
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
    row_g2 <- res[res$genes == "g2", , drop = FALSE]
    expect_true(nrow(row_g2) == 1)
    expect_true(is.na(row_g2$raw_p_values) || is.na(row_g2$adjusted_p_values))
})


test_that("calculate_fc input validation and pseudocount behavior", {
    mat <- matrix(c(1, 0, -1), nrow = 1)
    samples <- c("A", "A", "B")
    expect_error(calculate_fc(mat, samples = samples[-1], control = "A"), "Length of 'samples' must equal")
    expect_error(calculate_fc(mat, samples = samples, control = "C"), "Control sample type not found")

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
    expect_true(is.matrix(res))
    # result should be 1 row and 2 columns (raw and adjusted)
    expect_equal(dim(res), c(1, 2))
})
