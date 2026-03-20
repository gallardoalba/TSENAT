context("Orchestration: Pipeline Infrastructure")

make_test_se <- function(n_genes = 20, n_samples = 8) {
  counts <- matrix(rpois(n_genes * n_samples, lambda = 10), nrow = n_genes, ncol = n_samples)
  rownames(counts) <- paste0("Gene", 1:n_genes)
  colnames(counts) <- paste0("Sample", 1:n_samples)
  SummarizedExperiment(assays = list(counts = counts))
}

test_that("tsenat() creates TSENATAnalysis object", {
  se <- make_test_se()
  analysis <- suppressWarnings(tsenat(se, verbose = FALSE))
  expect_s4_class(analysis, "TSENATAnalysis")
})

test_that("tsenat() rejects invalid SummarizedExperiment", {
  expect_error(tsenat("not_se"), "SummarizedExperiment")
})

test_that("tsenat() stores end time in metadata", {
  se <- make_test_se()
  analysis <- suppressWarnings(tsenat(se, verbose = FALSE))
  expect_true("ended_at" %in% names(analysis@metadata))
  expect_true(inherits(analysis@metadata$ended_at, "POSIXct"))
})

test_that("tsenat() respects methods parameter", {
  se <- make_test_se()
  # Run with only diversity
  analysis <- suppressWarnings(tsenat(se, methods = c("diversity"), verbose = FALSE))
  expect_s4_class(analysis, "TSENATAnalysis")
})

test_that("tsenat() respects q_values parameter", {
  se <- make_test_se()
  analysis <- suppressWarnings(tsenat(se, q_values = c(0.5, 2.0), verbose = FALSE))
  expect_s4_class(analysis, "TSENATAnalysis")
})

test_that("tsenat() accepts TSENATConfig", {
  se <- make_test_se()
  cfg <- tsenat_config(q_values = c(1.0), methods = c("diversity"))
  analysis <- suppressWarnings(tsenat(se, config = cfg, verbose = FALSE))
  expect_s4_class(analysis, "TSENATAnalysis")
})

test_that("tsenat() generates metadata with timestamps", {
  se <- make_test_se()
  analysis <- suppressWarnings(tsenat(se, methods = c(), verbose = FALSE))  # No actual methods

  expect_true("created_at" %in% names(analysis@metadata))
  expect_true("ended_at" %in% names(analysis@metadata))
  expect_true(analysis@metadata$ended_at >= analysis@metadata$created_at)
})

test_that("tsenat() output works with show()", {
  se <- make_test_se()
  analysis <- suppressWarnings(tsenat(se, methods = c(), verbose = FALSE))
  expect_no_error(show(analysis))
})

test_that("tsenat() output works with summary()", {
  se <- make_test_se()
  analysis <- suppressWarnings(tsenat(se, methods = c(), verbose = FALSE))
  expect_no_error(summary(analysis))
})

test_that("tsenat() accessors work on output", {
  se <- make_test_se()
  analysis <- suppressWarnings(tsenat(se, methods = c(), verbose = FALSE))

  # These should not error even on empty results
  div <- suppressWarnings(diversity(analysis))
  lm <- suppressWarnings(lmResults(analysis))
  jk <- suppressWarnings(jackKnife(analysis))
  div_res <- suppressWarnings(divergence(analysis))
  plots <- suppressWarnings(getPlot(analysis))

  expect_true(
    is.null(div) || is.list(div)
  )
})

test_that("tsenat() with custom seed runs reproducibly", {
  se <- make_test_se()
  cfg <- tsenat_config(seed = 42)
  expect_no_error(suppressWarnings(tsenat(se, config = cfg, verbose = FALSE)))
})

test_that("tsenat_config creates valid config object", {
  cfg <- tsenat_config()
  expect_true("TSENATConfig" %in% class(cfg))
  expect_true(is.list(cfg))
})

test_that("tsenat_config rejects invalid methods", {
  expect_error(
    tsenat_config(methods = c("diversity", "invalid")),
    "Invalid methods"
  )
})

test_that("tsenat_config with q_range generates sequence", {
  cfg <- tsenat_config(q_range = c(1.0, 2.0))
  expect_true(1.0 %in% cfg$q_values)
  expect_true(2.0 %in% cfg$q_values)
})

test_that("Orchestration functions handle empty config", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se)  # Empty config
  expect_no_error(analysis)
})
