skip_on_bioc()

context("Additional generate_plots tests")

library(TSENAT)

# plot_ma_expression errors when provided precomputed se missing 'log2_fold_change'

test_that("plot_ma_expression errors when precomputed se lacks log2_fold_change", {
    x <- data.frame(genes = paste0("g", seq_len(3)), mean = runif(3), stringsAsFactors = FALSE)
    bad_se <- data.frame(genes = paste0("g", seq_len(3)), other = rnorm(3), stringsAsFactors = FALSE)
    expect_error(plot_ma_expression(x, se = bad_se), "`se` data.frame must contain 'log2_fold_change'")
})

# plot_ma_expression accepts precomputed fold changes provided as data.frame/matrix with rownames

test_that("plot_ma_expression accepts precomputed fold changes with rownames", {
    x <- data.frame(genes = paste0("g", seq_len(4)), mean = runif(4), stringsAsFactors = FALSE)
    fc <- data.frame(log2_fold_change = rnorm(4))
    rownames(fc) <- paste0("g", seq_len(4))
    p <- plot_ma_expression(x, se = fc)
    expect_s3_class(p, "gg")
    expect_s3_class(p, "ggplot")
})

# plot_top_transcripts writes output_file for single and multiple genes

test_that("plot_top_transcripts writes output files for single and multiple genes", {
    skip_if_not_installed("ggplot2")
    set.seed(10)
    counts <- matrix(rpois(6 * 4, lambda = 10), nrow = 6)
    rownames(counts) <- paste0("tx", 1:6)
    colnames(counts) <- paste0("S", 1:4)
    samples <- c("N", "N", "T", "T")
    tx2 <- data.frame(Transcript = rownames(counts), Gen = rep(c("G1", "G2"), each = 3), stringsAsFactors = FALSE)

    tf1 <- tempfile(fileext = ".png")
    plot_top_transcripts(counts, gene = "G1", samples = samples, tx2gene = tx2, output_file = tf1)
    expect_true(file.exists(tf1) && file.info(tf1)$size > 0)

    tf2 <- tempfile(fileext = ".png")
    plot_top_transcripts(counts, gene = c("G1", "G2"), samples = samples, tx2gene = tx2, output_file = tf2)
    expect_true(file.exists(tf2) && file.info(tf2)$size > 0)
})

# plot_top_transcripts supports metric = 'iqr'

test_that("plot_top_transcripts supports metric 'iqr'", {
    skip_if_not_installed("ggplot2")
    counts <- matrix(rpois(3 * 4, lambda = 5), nrow = 3)
    rownames(counts) <- paste0("tx", 1:3)
    colnames(counts) <- paste0("S", 1:4)
    samples <- c("N", "N", "T", "T")
    tx2 <- data.frame(Transcript = rownames(counts), Gen = rep("G1", 3), stringsAsFactors = FALSE)

    p <- plot_top_transcripts(counts, gene = "G1", samples = samples, tx2gene = tx2, metric = "iqr")
    expect_s3_class(p, "ggplot")
})

# plot_volcano auto-detects x_col when not provided

test_that("plot_volcano auto-detects a numeric x column when x_col is NULL", {
    skip_if_not_installed("ggplot2")
    df <- data.frame(
        genes = paste0("g", seq_len(10)),
        stat = rnorm(10),
        adjusted_p_values = p.adjust(runif(10)),
        stringsAsFactors = FALSE
    )
    p <- plot_volcano(df, x_col = NULL, padj_col = "adjusted_p_values")
    expect_s3_class(p, "ggplot")
})

# .plot_ma_core uses fc_df values when provided; verify y values in plot data correspond to fc_df

test_that(".plot_ma_core uses fc_df values when provided", {
    skip_if_not_installed("ggplot2")
    x <- data.frame(genes = paste0("g", seq_len(5)), mean = runif(5), log2_fold_change = rnorm(5), stringsAsFactors = FALSE)
    fc <- data.frame(genes = x$genes, log2_fold_change = rnorm(5, mean = 5, sd = 0.1), stringsAsFactors = FALSE)

    p <- .plot_ma_core(x, fc_df = fc)
    expect_s3_class(p, "ggplot")
    pb <- ggplot2::ggplot_build(p)
    plotted_y <- pb$data[[1]]$y
    # Reconstruct expected merged df per implementation
    fdf <- as.data.frame(fc, stringsAsFactors = FALSE)
    df <- as.data.frame(x, stringsAsFactors = FALSE)
    df <- merge(df, fdf[, c("genes", "log2_fold_change")], by = "genes", all.x = TRUE, suffixes = c("", ".fc"))
    if ("log2_fold_change.fc" %in% colnames(df)) df$log2_fold_change <- ifelse(!is.na(df$log2_fold_change.fc), df$log2_fold_change.fc, df$log2_fold_change)
    expected_y <- as.numeric(df$log2_fold_change)
    expect_equal(plotted_y, expected_y)
})
