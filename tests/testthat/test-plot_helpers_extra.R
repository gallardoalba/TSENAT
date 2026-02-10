context("plot_helpers extras")

library(testthat)

test_that(".tsenat_format_label handles various inputs", {
  expect_null(.tsenat_format_label(NULL))
  expect_equal(.tsenat_format_label("__FOO_bar  "), "Foo bar")
  expect_equal(.tsenat_format_label(" a "), "A")
  expect_equal(.tsenat_format_label("   "), "")
  expect_equal(.tsenat_format_label("SINGLE"), "Single")
})

test_that(".tsenat_prepare_ma_plot_df handles mean_cols length >=2 and significance detection", {
  df <- data.frame(genes = c('g1', 'g2', 'g3'), meanA = c(1,2,3), meanB = c(1.5, 1.5, 1.5), log2fc = c(0, 1.2, -0.5), padj = c(0.2, 0.01, NA), stringsAsFactors = FALSE)
  res <- .tsenat_prepare_ma_plot_df(df, fold_col = 'log2fc', mean_cols = c('meanA', 'meanB'), x_label = NULL, y_label = 'Log2FC')
  expect_is(res, 'list')
  # when mean_cols length>=2 and x_label is NULL, default to 'meanA vs meanB'
  expect_equal(res$x_label, 'meanA vs meanB')
  expect_true('plot_df' %in% names(res))
  expect_equal(nrow(res$plot_df), 3)
  # gene 2 should be significant (abs(y)>0 and padj<0.05)
  sig <- res$plot_df$significant
  expect_equal(sig, c('non-significant', 'significant', 'non-significant'))
})

test_that(".tsenat_prepare_ma_plot_df handles single mean col and fallback mean/index", {
  df1 <- data.frame(genes = c('g1','g2'), m = c(5,6), fc = c(0,2), stringsAsFactors = FALSE)
  r1 <- .tsenat_prepare_ma_plot_df(df1, fold_col = 'fc', mean_cols = c('m'), x_label = NULL, y_label = NULL)
  expect_equal(r1$x_label, 'm')
  expect_equal(r1$plot_df$x, as.numeric(c(5,6)))

  df2 <- data.frame(genes = c('g1','g2'), mean = c(3,4), fc = c(1, 0), stringsAsFactors = FALSE)
  r2 <- .tsenat_prepare_ma_plot_df(df2, fold_col = 'fc', mean_cols = character(0), x_label = NULL, y_label = NULL)
  expect_equal(r2$x_label, 'Mean')

  df3 <- data.frame(genes = c('g1','g2'), fc = c(1, 2), stringsAsFactors = FALSE)
  r3 <- .tsenat_prepare_ma_plot_df(df3, fold_col = 'fc', mean_cols = character(0), x_label = NULL, y_label = NULL)
  expect_equal(r3$x_label, 'Index')
  expect_equal(r3$plot_df$x, c(1,2))
})


test_that(".tsenat_prepare_volcano_df detects _difference column and formats labels", {
  df <- data.frame(gene = c('a','b','c'), median_difference = c(0.2, -0.5, 0.6), adjusted_p_values = c(0.2, 0.01, 0.001), stringsAsFactors = FALSE)
  res <- .tsenat_prepare_volcano_df(df)
  expect_equal(res$x_col, 'median_difference')
  expect_equal(res$padj_col, 'adjusted_p_values')
  expect_true('df' %in% names(res))
  expect_match(res$x_label_formatted, 'Median')
  expect_match(res$padj_label_formatted, 'Adjusted p values|Adjusted p values')
})

test_that(".tsenat_prepare_volcano_df errors for missing columns and empty data", {
  df <- data.frame(g = 1:3, something = letters[1:3], stringsAsFactors = FALSE)
  # Because 'g' is numeric it will be chosen as x_col but the default padj
  # column 'adjusted_p_values' is missing and an informative error is raised
  expect_error(.tsenat_prepare_volcano_df(df), "Column 'adjusted_p_values' not found")

  df2 <- data.frame(x = c(NA, Inf), adjusted_p_values = c(NA, NA), stringsAsFactors = FALSE)
  expect_error(.tsenat_prepare_volcano_df(df2, x_col = 'x'), "No valid points to plot")

  df3 <- data.frame(x = c(1,2), adj = c(0.01, 0.02), stringsAsFactors = FALSE)
  expect_error(.tsenat_prepare_volcano_df(df3, x_col = 'x', padj_col = 'nope'), "Column 'nope' not found")
})

test_that(".tsenat_prepare_volcano_df handles padj <=0 and signficance logic", {
  df <- data.frame(g = 1:4, value = c(0.2, 0.5, -0.2, 1), adjusted_p_values = c(0, 1e-10, 0.5, 0.001), stringsAsFactors = FALSE)
  res <- .tsenat_prepare_volcano_df(df, x_col = 'value')
  expect_true(all(res$df$padj > 0))
  # label_thresh default 0.1: check significance assignment
  sig <- res$df$significant
  expect_equal(sig, ifelse(abs(res$df$xval) >= 0.1 & res$df$padj < 0.05, 'significant', 'non-significant'))
})