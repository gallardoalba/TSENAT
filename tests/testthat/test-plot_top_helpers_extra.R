context("plot_top_helpers extra cases")

library(testthat)

# .ptt_select_genes_from_res
test_that(".ptt_select_genes_from_res errors on NULL or missing genes", {
    expect_error(.ptt_select_genes_from_res(NULL, 3), "Either 'gene' or 'res' must be provided")
    expect_error(.ptt_select_genes_from_res(data.frame(a = 1:3), 2), "must contain a 'genes' column")
})

test_that(".ptt_select_genes_from_res sorts by adjusted or raw p-values and returns unique genes", {
    res <- data.frame(genes = c('g1', 'g2', 'g1', 'g3'), adjusted_p_values = c(0.2, 0.01, 0.05, NA))
    sel <- .ptt_select_genes_from_res(res, top_n = 3)
    expect_true(all(c('g2', 'g1') %in% sel))
    expect_length(unique(sel), length(sel))

    # raw p-values fallback
    res2 <- data.frame(genes = c('a','b','c'), raw_p_values = c(0.5, 0.1, 0.3))
    sel2 <- .ptt_select_genes_from_res(res2, top_n = 2)
    expect_equal(sel2, c('b','c'))
})

# .ptt_infer_samples_from_coldata
test_that(".ptt_infer_samples_from_coldata handles row-named coldata and Sample id column", {
    counts <- matrix(1:6, nrow = 2)
    colnames(counts) <- c('S1','S2','S3')[1:ncol(counts)]

    cdf <- data.frame(sample_type = c('A','B','A'), stringsAsFactors = FALSE)
    rownames(cdf) <- c('S1','S2','S3')
    out <- .ptt_infer_samples_from_coldata(cdf, counts, 'sample_type')
    expect_equal(out, c('A','B','A')[1:ncol(counts)])

    # use Sample id column
    cdf2 <- data.frame(Sample = c('S1','S2','S3'), sample_type = c('A','B','A'), stringsAsFactors = FALSE)
    out2 <- .ptt_infer_samples_from_coldata(cdf2, counts, 'sample_type')
    expect_equal(out2, c('A','B','A')[1:ncol(counts)])

    # mismatched sample ids should error
    cdf_bad <- data.frame(Sample = c('X','Y','Z'), sample_type = c('A','B','A'), stringsAsFactors = FALSE)
    expect_error(.ptt_infer_samples_from_coldata(cdf_bad, counts, 'sample_type'), "coldata sample id column does not match")
    expect_error(.ptt_infer_samples_from_coldata(123, counts, 'sample_type'), "must be a data.frame or path")
})

# .ptt_read_tx2gene
test_that(".ptt_read_tx2gene validates inputs and reads mapping", {
    expect_error(.ptt_read_tx2gene(NULL), "`tx2gene` must be provided")

    bad <- data.frame(X = 1:2)
    expect_error(.ptt_read_tx2gene(bad), "tx2gene must have columns 'Transcript' and 'Gen'")

    good <- data.frame(Transcript = c('t1','t2'), Gen = c('g1','g1'), stringsAsFactors = FALSE)
    out <- .ptt_read_tx2gene(good)
    expect_equal(out, good)

    tf <- tempfile(fileext = '.tsv')
    write.table(good, file = tf, sep = '\t', row.names = FALSE, quote = FALSE)
    outf <- .ptt_read_tx2gene(tf)
    expect_true(is.data.frame(outf))
    unlink(tf)
    expect_error(.ptt_read_tx2gene('no_such_file.tsv'), "tx2gene file not found")
})

# .ptt_make_agg
test_that(".ptt_make_agg returns an aggregation function and label", {
    maj <- .ptt_make_agg('median')
    expect_equal(maj$metric_choice, 'median')
    expect_true(is.function(maj$agg_fun))
    expect_true(grepl('median', maj$agg_label_unique, ignore.case = TRUE))

    m2 <- .ptt_make_agg('iqr')
    expect_equal(m2$metric_choice, 'iqr')
    expect_true(grepl('IQR', m2$agg_label_unique))
})

# .ptt_build_tx_long and .ptt_aggregate_df_long
test_that(".ptt_build_tx_long and aggregation pipeline works and errors appropriately", {
    counts <- matrix(1:12, nrow = 4)
    rownames(counts) <- paste0('tx',1:4)
    colnames(counts) <- paste0('S',1:3)
    mapping <- data.frame(Transcript = paste0('tx',1:4), Gen = c('G1','G1','G2','G2'), stringsAsFactors = FALSE)
    samples <- c('A','B','A')

    expect_error(.ptt_build_tx_long('NOPE', mapping, counts, samples, top_n = NULL), 'No transcripts found')

    res <- .ptt_build_tx_long('G1', mapping, counts, samples, top_n = 1)
    expect_true(is.list(res))
    expect_true(all(c('df_long','txs') %in% names(res)))
    expect_equal(length(unique(res$df_long$tx)), length(res$txs))

    summ <- .ptt_aggregate_df_long(res$df_long, agg_fun = function(x) mean(x, na.rm = TRUE), pseudocount = 0.1)
    expect_true('log2expr' %in% colnames(summ))
    expect_true(is.factor(summ$tx))
})

# .ptt_build_plot_from_summary and combine functions
test_that(".ptt_build_plot_from_summary generates ggplot and combine functions operate", {
    skip_if_not_installed('ggplot2')
    p <- .ptt_build_plot_from_summary(data.frame(tx = factor(c('a','b')), group = c('A','B'), log2expr = c(1,2)), 'Label')
    expect_s3_class(p, 'gg')

    # patchwork combine
    if (rlang::is_installed('patchwork')) {
        skip_if_not_installed('patchwork')
        p2 <- p + p
        combined <- .ptt_combine_patchwork(list(p, p), 'Label')
        expect_true(inherits(combined, 'patchwork'))
    }

    if (rlang::is_installed('cowplot')) {
        skip_if_not_installed('cowplot')
        outp <- .ptt_combine_cowplot(list(p, p), output_file = NULL, agg_label_unique = 'Label')
        expect_true(inherits(outp, 'gtable') || inherits(outp, 'ggplot') || inherits(outp, 'grob'))
    }

    if (rlang::is_installed('grid')) {
        skip_if_not_installed('grid')
        # .ptt_combine_grid returns invisibly NULL when not writing file
        res_grid <- .ptt_combine_grid(list(p, p), output_file = NULL, agg_label_unique = 'Label')
        expect_null(res_grid)
    }
})
