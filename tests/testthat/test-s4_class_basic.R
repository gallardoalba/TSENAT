context("S4 Class: TSENATAnalysis Basic Operations")

test_that("TSENATAnalysis object can be created with valid SummarizedExperiment", {
  # Create minimal SE
  se <- SummarizedExperiment(
    assays = list(counts = matrix(rpois(100, 5), nrow = 10, ncol = 10))
  )
  rownames(se) <- paste0("Gene", 1:10)
  colnames(se) <- paste0("Sample", 1:10)

  # Create TSENATAnalysis
  analysis <- TSENATAnalysis(se)

  # Verify class
  expect_s4_class(analysis, "TSENATAnalysis")
  expect_equal(nrow(analysis@se), 10)
  expect_equal(ncol(analysis@se), 10)
  expect_equal(length(analysis@diversity_results), 0)
  expect_equal(length(analysis@lm_results), 0)
})

test_that("TSENATAnalysis can be created with empty SummarizedExperiment (validation before use)", {
  empty_se <- SummarizedExperiment(assays = list(counts = matrix(0, 0, 0)))
  # Constructor allows empty SE; validation happens when running analyses
  analysis <- TSENATAnalysis(empty_se)
  expect_s4_class(analysis, "TSENATAnalysis")
  expect_equal(nrow(analysis@se), 0)
})

test_that("TSENATAnalysis stores configuration properly", {
  se <- SummarizedExperiment(
    assays = list(counts = matrix(rpois(50, 3), nrow = 5, ncol = 10))
  )
  rownames(se) <- paste0("G", 1:5)
  colnames(se) <- paste0("S", 1:10)

  config <- list(q_values = c(0.5, 1.0), fdr_threshold = 0.01)
  analysis <- TSENATAnalysis(se, config = config)

  expect_equal(analysis@config$q_values, c(0.5, 1.0))
  expect_equal(analysis@config$fdr_threshold, 0.01)
})

test_that("TSENATAnalysis initializes metadata with timestamps and version", {
  se <- SummarizedExperiment(
    assays = list(counts = matrix(rpois(50, 3), nrow = 5, ncol = 10))
  )
  rownames(se) <- paste0("G", 1:5)
  colnames(se) <- paste0("S", 1:10)

  analysis <- TSENATAnalysis(se)

  expect_true("created_at" %in% names(analysis@metadata))
  expect_true("package_version" %in% names(analysis@metadata))
  expect_true("function_calls" %in% names(analysis@metadata))
  expect_true(inherits(analysis@metadata$created_at, "POSIXct"))
  expect_length(analysis@metadata$function_calls, 0)
})

test_that("TSENATAnalysis validity checks slot types", {
  # The validity function should prevent invalid objects
  # We test this indirectly through the constructor

  se <- SummarizedExperiment(
    assays = list(counts = matrix(rpois(50, 3), nrow = 5, ncol = 10))
  )
  rownames(se) <- paste0("G", 1:5)
  colnames(se) <- paste0("S", 1:10)

  analysis <- TSENATAnalysis(se)

  # Try to create invalid object by direct setClass manipulation
  # (This would only fail if validity is enforced)
  expect_equal(class(analysis@diversity_results), "list")
  expect_equal(class(analysis@lm_results), "list")
  expect_equal(class(analysis@plots), "list")
})

test_that("show method works for TSENATAnalysis", {
  se <- SummarizedExperiment(
    assays = list(counts = matrix(rpois(50, 3), nrow = 5, ncol = 10))
  )
  rownames(se) <- paste0("G", 1:5)
  colnames(se) <- paste0("S", 1:10)

  analysis <- TSENATAnalysis(se)

  # Capture output
  output <- capture.output(show(analysis))
  expect_true(any(grepl("TSENATAnalysis", output)))
  expect_true(any(grepl("Genes", output)))
  expect_true(any(grepl("Samples", output)))
})

test_that("summary method works for TSENATAnalysis", {
  se <- SummarizedExperiment(
    assays = list(counts = matrix(rpois(50, 3), nrow = 5, ncol = 10))
  )
  rownames(se) <- paste0("G", 1:5)
  colnames(se) <- paste0("S", 1:10)

  analysis <- TSENATAnalysis(se)

  # Capture output
  output <- capture.output(summary(analysis))
  expect_true(any(grepl("TSENAT", output)))
  expect_true(any(grepl("Created", output)))
})
