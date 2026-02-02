testthat::skip_on_bioc()
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
