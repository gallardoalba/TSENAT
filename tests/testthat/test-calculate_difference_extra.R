context("calculate_difference extra cases")

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
    # when shuffle used we expect raw_p_values and adjusted_p_values columns
    expect_true(all(c("raw_p_values", "adjusted_p_values") %in% colnames(res)))
})


test_that("Providing multiple samples column names to SummarizedExperiment errors", {
    skip_if_not_installed("SummarizedExperiment")
    mat <- matrix(runif(8), nrow = 2)
    colnames(mat) <- paste0("S", seq_len(ncol(mat)))
    colData_df <- S4Vectors::DataFrame(sample_type = rep(c("A","B","A","B"), length.out = ncol(mat)), row.names = colnames(mat))
    se <- SummarizedExperiment::SummarizedExperiment(assays = S4Vectors::SimpleList(counts = mat), colData = colData_df)

    expect_error(calculate_difference(se, samples = c("sample_type", "foo"), control = "A"), "'samples' must be a single colData column")
})
