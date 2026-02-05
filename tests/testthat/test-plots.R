testthat::skip_on_bioc()

# Combined plotting tests: density/violin/MA/top-transcripts, volcano, q-curve

testthat::test_that("plot_diversity_density returns ggplot object", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("tidyr")
    skip_if_not_installed("dplyr")

    library(SummarizedExperiment)
    library(ggplot2)
    library(tidyr)
    library(dplyr)

    # construct minimal SummarizedExperiment
    mat <- matrix(runif(20), nrow = 5, ncol = 4)
    colnames(mat) <- c("S1_N", "S2_T", "S3_N", "S4_T")
    rownames(mat) <- paste0("G", 1:5)
    rowData_df <- DataFrame(genes = rownames(mat))
    colData_df <- DataFrame(samples = colnames(mat))
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rowData_df,
        colData = colData_df
    )

    p <- plot_diversity_density(se)
    expect_s3_class(p, "gg")
})


testthat::test_that("plot_mean_violin returns ggplot object", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("tidyr")
    skip_if_not_installed("dplyr")

    library(SummarizedExperiment)
    library(ggplot2)
    library(tidyr)
    library(dplyr)

    mat <- matrix(runif(20), nrow = 5, ncol = 4)
    colnames(mat) <- c("S1_N", "S2_T", "S3_N", "S4_T")
    rownames(mat) <- paste0("G", 1:5)
    rowData_df <- DataFrame(genes = rownames(mat))
    colData_df <- DataFrame(samples = colnames(mat))
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rowData_df,
        colData = colData_df
    )

    p <- plot_mean_violin(se)
    expect_s3_class(p, "gg")
})


testthat::test_that("plot_ma returns ggplot object", {
    skip_if_not_installed("ggplot2")
    library(ggplot2)

    df <- data.frame(
        Gene = paste0(
            "G",
            1:10
        ),
        A_mean = runif(10),
        B_mean = runif(10),
        log2_fold_change = rnorm(10),
        adjusted_p_values = runif(10)
    )
    p <- plot_ma(df)
    expect_s3_class(p, "gg")
})


testthat::test_that("plot_top_transcripts returns ggplot for synthetic data", {
    skip_on_cran()
    set.seed(42)
    counts <- matrix(rpois(3 * 8, lambda = 20), nrow = 3)
    rownames(counts) <- paste0("tx", 1:3)
    colnames(counts) <- paste0("S", 1:8)
    samples <- c(rep("Normal", 4), rep("Tumor", 4))
    # create simple tx2gene mapping
    tx2 <- data.frame(
        Transcript = rownames(counts),
        Gen = rep(
            "GENE1",
            3
        ),
        stringsAsFactors = FALSE
    )

    p <- plot_top_transcripts(counts,
        gene = "GENE1",
        samples = samples,
        tx2gene = tx2,
        top_n = 2
    )
    testthat::expect_s3_class(p, "ggplot")
})


testthat::test_that("plot_volcano returns a ggplot and annotates top genes", {
    skip_on_cran()
    set.seed(42)
    n <- 20
    df <- data.frame(
        genes = paste0("gene", seq_len(n)),
        mean_difference = rnorm(n),
        adjusted_p_values = p.adjust(runif(n))
    )

    p <- plot_volcano(df,
        x_col = "mean_difference",
        padj_col = "adjusted_p_values",
        top_n = 3
    )
    testthat::expect_s3_class(p, "ggplot")
    # building the plot should not error
    ggplot2::ggplot_build(p)
})


testthat::test_that("plot_tsallis_q_curve returns a ggplot object (SE-first)", {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
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

    p <- plot_tsallis_q_curve(ts_se)
    expect_true(inherits(p, "ggplot"))
})

context("plot_ma selection logic")

test_that("plot_ma prefers mean columns when available", {
  df <- data.frame(A_mean = runif(10), B_mean = runif(10),
                   log2_fold_change = rnorm(10), adjusted_p_values = runif(10))
  p <- plot_ma(df)
  expect_s3_class(p, "ggplot")
})

test_that("plot_ma falls back to median columns when means absent", {
  df <- data.frame(A_median = runif(10), B_median = runif(10),
                   log2_fold_change = rnorm(10), adjusted_p_values = runif(10))
  p <- plot_ma(df)
  expect_s3_class(p, "ggplot")
})

test_that("plot_ma errors on mixed mean/median columns", {
  df <- data.frame(A_mean = runif(10), B_median = runif(10),
                   log2_fold_change = rnorm(10), adjusted_p_values = runif(10))
  expect_error(plot_ma(df), "Could not find two mean or two median columns")
})

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

    p <- plot_tsallis_q_curve(ts_se)

    expect_true(inherits(p, "ggplot"))
})

test_that("matrix input is rejected with informative message", {
    set.seed(2)
    readcounts <- matrix(rpois(20 * 2, lambda = 5), nrow = 20, ncol = 2)
    colnames(readcounts) <- c("A_N", "B_T")
    genes <- rep(paste0("G", 1:5), length.out = nrow(readcounts))

    expect_error(
        plot_tsallis_q_curve(readcounts),
        "plot_tsallis_q_curve requires a SummarizedExperiment from calculate_diversity."
    )
})
