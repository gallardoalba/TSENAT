context("S4 Methods: Accessor Functions")

# Helper function to create test analysis object
make_test_analysis <- function() {
  se <- SummarizedExperiment(
    assays = list(counts = matrix(rpois(100, 5), nrow = 10, ncol = 10))
  )
  rownames(se) <- paste0("Gene", 1:10)
  colnames(se) <- paste0("Sample", 1:10)
  TSENATAnalysis(se)
}

test_that("diversity accessor returns NULL when no results available", {
  analysis <- make_test_analysis()
  result <- expect_warning(diversity(analysis), "No diversity results found")
  expect_null(result)
})

test_that("diversity accessor returns all results when no q specified", {
  analysis <- make_test_analysis()

  # Manually populate diversity results
  se_q1 <- SummarizedExperiment(
    assays = list(counts = matrix(rnorm(100), nrow = 10))
  )
  analysis@diversity_results$q_1.0 <- se_q1

  result <- diversity(analysis)
  expect_type(result, "list")
  expect_true("q_1.0" %in% names(result))
})

test_that("diversity accessor extracts specific q-value", {
  analysis <- make_test_analysis()

  # Add results for multiple q-values
  se_q1 <- SummarizedExperiment(
    assays = list(counts = matrix(rnorm(100), nrow = 10))
  )
  se_q2 <- SummarizedExperiment(
    assays = list(counts = matrix(rnorm(100), nrow = 10))
  )
  analysis@diversity_results$q_1.0 <- se_q1
  analysis@diversity_results$q_2.0 <- se_q2

  result_q1 <- diversity(analysis, q = 1.0)
  result_q2 <- diversity(analysis, q = 2.0)

  expect_s4_class(result_q1, "SummarizedExperiment")
  expect_s4_class(result_q2, "SummarizedExperiment")
  expect_false(identical(result_q1, result_q2))
})

test_that("diversity accessor errors on missing q-value", {
  analysis <- make_test_analysis()
  analysis@diversity_results$q_1.0 <- SummarizedExperiment(
    assays = list(counts = matrix(rnorm(100), nrow = 10))
  )

  expect_error(diversity(analysis, q = 3.0), "not found")
})

test_that("lmResults accessor returns NULL when no results available", {
  analysis <- make_test_analysis()
  result <- expect_warning(lmResults(analysis), "No LM results found")
  expect_null(result)
})

test_that("lmResults accessor returns all results when no component specified", {
  analysis <- make_test_analysis()

  # Add mock LM results
  analysis@lm_results$main <- list(results = data.frame(gene = 1:5, pval = runif(5)))
  analysis@lm_results$interaction <- list(results = data.frame(gene = 1:5, pval = runif(5)))

  result <- lmResults(analysis)
  expect_type(result, "list")
  expect_equal(length(result), 2)
})

test_that("lmResults accessor extracts specific component", {
  analysis <- make_test_analysis()

  df_main <- data.frame(gene = 1:5, pval = runif(5), estimate = rnorm(5))
  analysis@lm_results$main <- list(results = df_main)

  result <- lmResults(analysis, component = "main")
  expect_type(result, "list")
  expect_true("results" %in% names(result))
})

test_that("jackKnife accessor returns NULL when no results available", {
  analysis <- make_test_analysis()
  result <- expect_warning(jackKnife(analysis), "No jackknife results found")
  expect_null(result)
})

test_that("jackKnife accessor extracts specific q-value", {
  analysis <- make_test_analysis()

  # Add jackknife results
  jk_data <- list(
    confidence_intervals = data.frame(gene = 1:5, ci_lower = runif(5), ci_upper = runif(5) + 1),
    resamples = matrix(rnorm(50), nrow = 5)
  )
  analysis@jackknife_results$q_1.0 <- jk_data

  result <- jackKnife(analysis, q = 1.0)
  expect_type(result, "list")
  expect_true("confidence_intervals" %in% names(result))
})

test_that("divergence accessor returns NULL when no results available", {
  analysis <- make_test_analysis()
  result <- expect_warning(divergence(analysis), "No divergence results found")
  expect_null(result)
})

test_that("divergence accessor returns all results when no component specified", {
  analysis <- make_test_analysis()

  analysis@divergence_results$tsallis <- data.frame(gene = 1:5, div = runif(5))
  analysis@divergence_results$effect_size <- data.frame(gene = 1:5, es = rnorm(5))

  result <- divergence(analysis)
  expect_type(result, "list")
  expect_equal(length(result), 2)
})

test_that("divergence accessor extracts specific component", {
  analysis <- make_test_analysis()

  df <- data.frame(gene = 1:5, tsallis_div = runif(5))
  analysis@divergence_results$tsallis <- df

  result <- divergence(analysis, component = "tsallis")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 5)
})

test_that("getPlot returns NULL when no plots available", {
  analysis <- make_test_analysis()
  result <- expect_warning(getPlot(analysis), "No plots found")
  expect_null(result)
})

test_that("getPlot retrieves specific plot type", {
  analysis <- make_test_analysis()

  library(ggplot2)
  p <- ggplot() + geom_point(aes(x = 1, y = 1))
  analysis@plots$q_curve <- p

  result <- getPlot(analysis, type = "q_curve")
  expect_s3_class(result, "ggplot")
})

test_that("getPlot returns all plots when no type specified", {
  analysis <- make_test_analysis()

  library(ggplot2)
  p1 <- ggplot() + geom_point(aes(x = 1, y = 1))
  p2 <- ggplot() + geom_line(aes(x = c(1, 2), y = c(1, 2)))

  analysis@plots$q_curve <- p1
  analysis@plots$divergence <- p2

  result <- getPlot(analysis)
  expect_type(result, "list")
  expect_equal(length(result), 2)
})

test_that("addPlot stores plot in TSENATAnalysis", {
  analysis <- make_test_analysis()

  library(ggplot2)
  p <- ggplot() + geom_point(aes(x = 1, y = 1))

  analysis <- addPlot(analysis, type = "test_plot", plot = p, replace = FALSE)

  expect_true("test_plot" %in% names(analysis@plots))
  retrieved <- getPlot(analysis, type = "test_plot")
  expect_s3_class(retrieved, "ggplot")
})

test_that("addPlot refuses to overwrite by default", {
  analysis <- make_test_analysis()

  library(ggplot2)
  p1 <- ggplot() + geom_point(aes(x = 1, y = 1))
  p2 <- ggplot() + geom_line(aes(x = c(1, 2), y = c(1, 2)))

  analysis <- addPlot(analysis, type = "test_plot", plot = p1)
  expect_warning(
    addPlot(analysis, type = "test_plot", plot = p2, replace = FALSE),
    "already exists"
  )

  # Should return unchanged object
  expect_equal(length(analysis@plots), 1)
})

test_that("addPlot overwrites when replace=TRUE", {
  analysis <- make_test_analysis()

  library(ggplot2)
  p1 <- ggplot() + geom_point(aes(x = 1, y = 1))
  p2 <- ggplot() + geom_line(aes(x = c(1, 2), y = c(1, 2)))

  analysis <- addPlot(analysis, type = "test_plot", plot = p1)
  analysis <- addPlot(analysis, type = "test_plot", plot = p2, replace = TRUE)

  expect_equal(length(analysis@plots), 1)
  retrieved <- getPlot(analysis, type = "test_plot")
  expect_s3_class(retrieved, "ggplot")
})
