library(testthat)

context("assumptions utilities")

test_that(".format_table_value handles special values and formats", {
  f <- TSENAT:::.format_table_value
  expect_equal(f(0), "0")
  expect_equal(f(NA_real_), as.character(NA_real_))
  expect_match(f(1e-6), "e")
  expect_equal(f(1.23456, decimals = 2), "1.23")
})

test_that(".calculate_assumptions returns placeholder for empty data", {
  calc <- TSENAT:::.calculate_assumptions
  res <- calc(matrix(nrow = 0, ncol = 0))
  expect_s3_class(res, "rank_assumptions")
  expect_true(!is.null(attr(res, "checks")))
  checks <- attr(res, "checks")
  expect_true(!is.null(checks$exchangeability))
  expect_true(!is.null(res$summary_stats))
  expect_equal(res$summary_stats$n_genes, 0)
})

test_that(".compute_concurvity_index skips when mgcv missing or q_values missing", {
  conc <- TSENAT:::.compute_concurvity_index
  # Provide simple matrix and no q_values -> should return error placeholder
  m <- matrix(runif(10), nrow = 5)
  out <- conc(m, q_values = NULL)
  expect_true(grepl("q-values required", out$details) || grepl("mgcv not available", out$status))
})

test_that(".compute_working_correlation_fit returns sensible structure info", {
  wc <- TSENAT:::.compute_working_correlation_fit
  mat <- matrix(rnorm(30), nrow = 5)
  out <- wc(mat, assumed_structure = "exchangeable")
  expect_true(!is.null(out$mean_autocorrelation))
  expect_true(!is.null(out$suitable_structures))
  expect_true(is.character(out$status))
})
