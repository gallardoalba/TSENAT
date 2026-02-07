# Comprehensive testing of all plotting functions
# Tests plot_diversity_density, plot_mean_violin, plot_ma, 
# plot_top_transcripts, plot_volcano, plot_tsallis_q_curve,
# plot_tsallis_gene_profile, plot_tsallis_density_multq, plot_tsallis_violin_multq

library(TSENAT)
skip_on_bioc()

test_that("plot_diversity_density returns ggplot object with valid data", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("tidyr")
    skip_if_not_installed("dplyr")

    library(SummarizedExperiment)
    library(ggplot2)

    # construct minimal SummarizedExperiment
    mat <- matrix(runif(20), nrow = 5, ncol = 4)
    colnames(mat) <- c("S1_N", "S2_T", "S3_N", "S4_T")
    rownames(mat) <- paste0("G", 1:5)
    rowData_df <- S4Vectors::DataFrame(genes = rownames(mat))
    colData_df <- S4Vectors::DataFrame(samples = colnames(mat))
    se <- SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rowData_df,
        colData = colData_df
    )

    p <- plot_diversity_density(se)
    expect_s3_class(p, "gg")
    expect_s3_class(p, "ggplot")
})

test_that("plot_mean_violin returns ggplot object", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("ggplot2")

    library(SummarizedExperiment)
    library(ggplot2)

    mat <- matrix(runif(20), nrow = 5, ncol = 4)
    colnames(mat) <- c("S1_N", "S2_T", "S3_N", "S4_T")
    rownames(mat) <- paste0("G", 1:5)
    rowData_df <- S4Vectors::DataFrame(genes = rownames(mat))
    colData_df <- S4Vectors::DataFrame(samples = colnames(mat))
    se <- SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rowData_df,
        colData = colData_df
    )

    p <- plot_mean_violin(se)
    expect_s3_class(p, "gg")
    expect_s3_class(p, "ggplot")
})

test_that("plot_ma returns ggplot object with mean columns", {
    skip_if_not_installed("ggplot2")
    library(ggplot2)

    df <- data.frame(
        Gene = paste0("G", 1:10),
        A_mean = runif(10),
        B_mean = runif(10),
        log2_fold_change = rnorm(10),
        adjusted_p_values = runif(10)
    )
    p <- plot_ma(df)
    expect_s3_class(p, "gg")
    expect_s3_class(p, "ggplot")
})

test_that("plot_ma returns ggplot object with median columns", {
    skip_if_not_installed("ggplot2")
    library(ggplot2)

    df <- data.frame(
        Gene = paste0("G", 1:10),
        A_median = runif(10),
        B_median = runif(10),
        log2_fold_change = rnorm(10),
        adjusted_p_values = runif(10)
    )
    p <- plot_ma(df)
    expect_s3_class(p, "gg")
})

test_that("plot_ma errors on mixed mean/median columns", {
    skip_if_not_installed("ggplot2")
    library(ggplot2)

    df <- data.frame(
        Gene = paste0("G", 1:10),
        A_mean = runif(10),
        B_median = runif(10),
        log2_fold_change = rnorm(10),
        adjusted_p_values = runif(10)
    )
    expect_error(plot_ma(df), "Could not find two mean or two median columns")
})

test_that("plot_top_transcripts returns ggplot for synthetic data", {
    skip_if_not_installed("ggplot2")
    library(ggplot2)

    set.seed(42)
    counts <- matrix(rpois(3 * 8, lambda = 20), nrow = 3)
    rownames(counts) <- paste0("tx", 1:3)
    colnames(counts) <- paste0("S", 1:8)
    samples <- c(rep("Normal", 4), rep("Tumor", 4))
    # create simple tx2gene mapping
    tx2 <- data.frame(
        Transcript = rownames(counts),
        Gen = rep("GENE1", 3),
        stringsAsFactors = FALSE
    )

    p <- plot_top_transcripts(counts,
        gene = "GENE1",
        samples = samples,
        tx2gene = tx2,
        top_n = 2
    )
    expect_s3_class(p, "ggplot")
})

test_that("plot_volcano returns a ggplot and annotates top genes", {
    skip_if_not_installed("ggplot2")
    library(ggplot2)

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
    expect_s3_class(p, "ggplot")
    # building the plot should not error
    ggplot2::ggplot_build(p)
})

test_that("plot_volcano with custom columns", {
    skip_if_not_installed("ggplot2")
    library(ggplot2)

    set.seed(42)
    n <- 15
    df <- data.frame(
        genes = paste0("gene", seq_len(n)),
        logFC = rnorm(n),
        pval = p.adjust(runif(n))
    )

    p <- plot_volcano(df,
        x_col = "logFC",
        padj_col = "pval",
        top_n = 2
    )
    expect_s3_class(p, "ggplot")
})

test_that("plot_tsallis_q_curve returns ggplot with valid SE", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("tidyr")
    skip_if_not_installed("dplyr")

    library(SummarizedExperiment)
    library(ggplot2)

    set.seed(1)
    readcounts <- matrix(rpois(30 * 3, lambda = 10), nrow = 30, ncol = 3)
    colnames(readcounts) <- c("S1_N", "S2_T", "S3_N")
    genes <- rep(paste0("G", 1:10), length.out = nrow(readcounts))

    qvals <- seq(0.01, 0.05, by = 0.01)
    ts_se <- calculate_diversity(readcounts, genes, q = qvals, norm = TRUE)

    coldata_df <- data.frame(
        Sample = c("S1_N", "S2_T", "S3_N"),
        Condition = c("Normal", "Tumor", "Normal"),
        stringsAsFactors = FALSE
    )

    ts_se <- map_coldata_to_se(ts_se, coldata_df)

    p <- plot_tsallis_q_curve(ts_se)
    expect_true(inherits(p, "ggplot"))
})

test_that("plot_tsallis_gene_profile returns ggplot for single gene", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("tidyr")
    skip_if_not_installed("dplyr")

    library(SummarizedExperiment)
    library(ggplot2)

    set.seed(123)
    readcounts <- matrix(rpois(30 * 4, lambda = 15), nrow = 30, ncol = 4)
    colnames(readcounts) <- c("S1_N", "S2_N", "S3_T", "S4_T")
    genes <- rep(paste0("G", 1:10), length.out = nrow(readcounts))

    qvals <- seq(0.01, 0.05, by = 0.02)
    ts_se <- calculate_diversity(readcounts, genes, q = qvals, norm = TRUE)

    coldata_df <- data.frame(
        Sample = c("S1_N", "S2_N", "S3_T", "S4_T"),
        Condition = c("Normal", "Normal", "Tumor", "Tumor"),
        stringsAsFactors = FALSE
    )

    ts_se <- map_coldata_to_se(ts_se, coldata_df)

    p <- plot_tsallis_gene_profile(ts_se, gene = "G1")
    expect_s3_class(p, "ggplot")
})

test_that("plot_tsallis_density_multq returns ggplot", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("tidyr")
    skip_if_not_installed("dplyr")

    library(SummarizedExperiment)
    library(ggplot2)

    set.seed(456)
    readcounts <- matrix(rpois(25 * 4, lambda = 12), nrow = 25, ncol = 4)
    colnames(readcounts) <- c("S1_N", "S2_N", "S3_T", "S4_T")
    genes <- rep(paste0("G", 1:5), length.out = nrow(readcounts))

    qvals <- seq(0.01, 0.1, by = 0.03)
    ts_se <- calculate_diversity(readcounts, genes, q = qvals, norm = TRUE)

    coldata_df <- data.frame(
        Sample = c("S1_N", "S2_N", "S3_T", "S4_T"),
        Condition = c("Normal", "Normal", "Tumor", "Tumor"),
        stringsAsFactors = FALSE
    )

    ts_se <- map_coldata_to_se(ts_se, coldata_df)

    p <- plot_tsallis_density_multq(ts_se)
    expect_s3_class(p, "ggplot")
})

test_that("plot_tsallis_violin_multq returns ggplot", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("tidyr")
    skip_if_not_installed("dplyr")

    library(SummarizedExperiment)
    library(ggplot2)

    set.seed(789)
    readcounts <- matrix(rpois(20 * 4, lambda = 18), nrow = 20, ncol = 4)
    colnames(readcounts) <- c("S1_N", "S2_N", "S3_T", "S4_T")
    genes <- rep(paste0("G", 1:4), length.out = nrow(readcounts))

    qvals <- c(0.01, 0.05, 0.1)
    ts_se <- calculate_diversity(readcounts, genes, q = qvals, norm = TRUE)

    coldata_df <- data.frame(
        Sample = c("S1_N", "S2_N", "S3_T", "S4_T"),
        Condition = c("Normal", "Normal", "Tumor", "Tumor"),
        stringsAsFactors = FALSE
    )

    ts_se <- map_coldata_to_se(ts_se, coldata_df)

    p <- plot_tsallis_violin_multq(ts_se)
    expect_s3_class(p, "ggplot")
})

test_that("All plot functions produce buildable ggplot objects", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("ggplot2")

    library(SummarizedExperiment)
    library(ggplot2)

    # Create test data
    mat <- matrix(runif(15), nrow = 5, ncol = 3)
    colnames(mat) <- c("S1_N", "S2_T", "S3_N")
    rownames(mat) <- paste0("G", 1:5)
    rowData_df <- S4Vectors::DataFrame(genes = rownames(mat))
    colData_df <- S4Vectors::DataFrame(samples = colnames(mat))
    se <- SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rowData_df,
        colData = colData_df
    )

    # Test that plots can be built without error
    expect_no_error(ggplot_build(plot_diversity_density(se)))
    expect_no_error(ggplot_build(plot_mean_violin(se)))
})
