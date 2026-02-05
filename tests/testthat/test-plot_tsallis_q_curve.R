library(testthat)

testthat::skip_on_bioc()

context("plot_tsallis_q_curve")

test_that("plot_tsallis_q_curve returns a ggplot object", {
    if (!requireNamespace("SummarizedExperiment",
        quietly = TRUE
    )) {
        skip("SummarizedExperiment required")
    }
    if (!requireNamespace("ggplot2", quietly = TRUE)) skip("ggplot2 required")
    if (!requireNamespace("tidyr", quietly = TRUE)) skip("tidyr required")
    if (!requireNamespace("dplyr", quietly = TRUE)) skip("dplyr required")

    set.seed(1)
    # small synthetic dataset: 30 transcripts, 3 samples
    readcounts <- matrix(rpois(30 * 3, lambda = 10), nrow = 30, ncol = 3)
    colnames(readcounts) <- c("S1_N", "S2_T", "S3_N")
    genes <- rep(paste0("G", 1:10), length.out = nrow(readcounts))

    # Compute a SummarizedExperiment (SE-first API) and map sample metadata
    qvals <- seq(0.01, 0.05, by = 0.01)
    ts_se <- calculate_diversity(readcounts, genes, q = qvals, norm = TRUE)

    # Provide simple coldata matching base sample names (without _q=... suffix)
    coldata_df <- data.frame(
        Sample = c("S1_N", "S2_T", "S3_N"),
        Condition = c("Normal", "Tumor", "Normal"),
        stringsAsFactors = FALSE
    )

    ts_se <- map_coldata_to_se(ts_se, coldata_df)

    p <- plot_tsallis_q_curve(ts_se, genes, q_values = qvals)

    expect_true(inherits(p, "ggplot"))
})

test_that("matrix input is rejected with informative message", {
    set.seed(2)
    readcounts <- matrix(rpois(20 * 2, lambda = 5), nrow = 20, ncol = 2)
    colnames(readcounts) <- c("A_N", "B_T")
    genes <- rep(paste0("G", 1:5), length.out = nrow(readcounts))

    expect_error(
        plot_tsallis_q_curve(readcounts, genes, q_values = seq(0.01, 0.05, by = 0.01)),
        "Matrix or data.frame input is no longer supported"
    )
})
