skip_on_bioc()

context("Unit tests for plot_top_transcripts helpers")

library(TSENAT)

# .ptt_select_genes_from_res

test_that(".ptt_select_genes_from_res selects by adjusted_p_values and raw_p_values", {
    res1 <- data.frame(genes = c("A", "B", "C"), adjusted_p_values = c(0.05, 0.01, 0.2), stringsAsFactors = FALSE)
    expect_equal(TSENAT:::.ptt_select_genes_from_res(res1, top_n = 2), c("B", "A"))

    res2 <- data.frame(genes = c("X", "Y", "Z"), raw_p_values = c(0.2, 0.01, 0.05), stringsAsFactors = FALSE)
    expect_equal(TSENAT:::.ptt_select_genes_from_res(res2, top_n = 2), c("Y", "Z"))

    expect_error(TSENAT:::.ptt_select_genes_from_res(NULL, top_n = 2))
    expect_error(TSENAT:::.ptt_select_genes_from_res(data.frame(a = 1), top_n = 2))
})

# .ptt_infer_samples_from_coldata

test_that(".ptt_infer_samples_from_coldata infers samples from data.frame and file path and errors on mismatch", {
    counts <- matrix(1:8, ncol = 4)
    colnames(counts) <- paste0("S", 1:4)

    cdf <- data.frame(sample_type = c("N", "T", "N", "T"), stringsAsFactors = FALSE)
    rownames(cdf) <- colnames(counts)

    samp <- TSENAT:::.ptt_infer_samples_from_coldata(cdf, counts, sample_type_col = "sample_type")
    expect_equal(as.character(samp), as.character(cdf[colnames(counts), "sample_type"]))

    # write as file with sample id column
    tf <- tempfile(fileext = ".tsv")
    dff <- data.frame(sample = colnames(counts), sample_type = c("N", "T", "N", "T"), stringsAsFactors = FALSE)
    utils::write.table(dff, file = tf, sep = "\t", quote = FALSE, row.names = FALSE)

    samp2 <- TSENAT:::.ptt_infer_samples_from_coldata(tf, counts, sample_type_col = "sample_type")
    expect_equal(as.character(samp2), as.character(dff$sample_type))

    # mismatch
    badcdf <- data.frame(other = c("a", "b"))
    expect_error(TSENAT:::.ptt_infer_samples_from_coldata(badcdf, counts, sample_type_col = "sample_type"))
})

# .ptt_read_tx2gene

test_that(".ptt_read_tx2gene reads mapping from data.frame and file and errors on missing columns", {
    mapping <- data.frame(Transcript = c("t1", "t2"), Gen = c("G1", "G1"), stringsAsFactors = FALSE)
    out <- TSENAT:::.ptt_read_tx2gene(mapping)
    expect_equal(out, mapping)

    tf <- tempfile(fileext = ".tsv")
    utils::write.table(mapping, file = tf, sep = "\t", quote = FALSE, row.names = FALSE)
    out2 <- TSENAT:::.ptt_read_tx2gene(tf)
    expect_equal(out2$Transcript, mapping$Transcript)

    expect_error(TSENAT:::.ptt_read_tx2gene(data.frame(a = 1)))
})

# .ptt_make_agg

test_that(".ptt_make_agg returns correct aggregator and label", {
    med <- TSENAT:::.ptt_make_agg("median")
    expect_equal(med$metric_choice, "median")
    expect_equal(med$agg_fun(c(1, 2, NA)), stats::median(c(1, 2, NA), na.rm = TRUE))

    mn <- TSENAT:::.ptt_make_agg("mean")
    expect_equal(mn$agg_fun(c(1, 2, NA)), mean(c(1, 2, NA), na.rm = TRUE))

    varr <- TSENAT:::.ptt_make_agg("variance")
    expect_equal(varr$agg_fun(c(1, 2, 3, NA)), stats::var(c(1, 2, 3, NA), na.rm = TRUE))

    iq <- TSENAT:::.ptt_make_agg("iqr")
    expect_equal(iq$agg_fun(c(1, 2, 3, 4, NA)), stats::IQR(c(1, 2, 3, 4, NA), na.rm = TRUE))

    # check counter side-effect increments
    opt_before <- as.integer(getOption("TSENAT.plot_top_counter", 0))
    TSENAT:::.ptt_make_agg("median")
    expect_true(as.integer(getOption("TSENAT.plot_top_counter", 0)) >= opt_before + 1)
})

# .ptt_build_tx_long & .ptt_aggregate_df_long & .ptt_build_plot_from_summary

test_that("tx long building, aggregation and plot building behave correctly", {
    counts <- matrix(rpois(6 * 2, lambda = 10), nrow = 6)
    rownames(counts) <- paste0("tx", 1:6)
    colnames(counts) <- paste0("S", 1:2)
    mapping <- data.frame(Transcript = rownames(counts), Gen = rep("G1", 6), stringsAsFactors = FALSE)
    samples <- c("N", "T")

    built <- TSENAT:::.ptt_build_tx_long("G1", mapping, counts, samples, top_n = 3)
    expect_true(is.list(built))
    expect_true(all(c("df_long", "txs") %in% names(built)))
    expect_true(length(built$txs) <= 3)
    expect_true(all(c("tx", "sample", "expr", "group") %in% colnames(built$df_long)))

    df_summary <- TSENAT:::.ptt_aggregate_df_long(built$df_long, agg_fun = function(x) mean(x, na.rm = TRUE), pseudocount = 1e-6)
    expect_true(all(c("tx", "group", "expr", "log2expr") %in% colnames(df_summary)))
    expect_true(is.factor(df_summary$tx))

    skip_if_not_installed("ggplot2")
    p <- TSENAT:::.ptt_build_plot_from_summary(df_summary, agg_label_unique = "label")
    expect_s3_class(p, "gg")
})

# .ptt_combine_grid writes to file when output_file provided

test_that(".ptt_combine_grid writes a PNG file when output_file is given", {
    skip_if_not_installed("ggplot2")
    library(ggplot2)

    df <- data.frame(x = 1:3, y = rnorm(3))
    p1 <- ggplot(df, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point()
    p2 <- ggplot(df, ggplot2::aes(x = x, y = -y)) + ggplot2::geom_point()

    tf <- tempfile(fileext = ".png")
    # call grid combiner directly
    TSENAT:::.ptt_combine_grid(list(p1, p2), output_file = tf, agg_label_unique = "agg")
    expect_true(file.exists(tf))
    expect_true(file.info(tf)$size > 0)
})
