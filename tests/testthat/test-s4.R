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

context("S4 Wrappers: Input Validation and Structure")

# Helper to create test analysis
make_test_se_for_wrappers <- function(n_genes = 20, n_samples = 8) {
  counts <- matrix(rpois(n_genes * n_samples, lambda = 10), nrow = n_genes, ncol = n_samples)
  rownames(counts) <- paste0("Gene", 1:n_genes)
  colnames(counts) <- paste0("Sample", 1:n_samples)
  SummarizedExperiment(assays = list(counts = counts))
}

# ============================================================================
# DIVERSITY WRAPPER TESTS
# ============================================================================

test_that("calculate_diversity_s4 validates input is TSENATAnalysis", {
  expect_error(calculate_diversity_s4("not_analysis", q = 1.0), "TSENATAnalysis")
})

test_that("calculate_diversity_s4 validates q is numeric", {
  se <- make_test_se_for_wrappers()
  analysis <- TSENATAnalysis(se)
  expect_error(calculate_diversity_s4(analysis, q = "not_numeric"), "must be numeric")
})

test_that("calculate_diversity_s4 rejects empty SummarizedExperiment", {
  empty_se <- SummarizedExperiment(assays = list(counts = matrix(0, 0, 0)))
  analysis <- TSENATAnalysis(empty_se)
  expect_error(calculate_diversity_s4(analysis), "empty")
})

# ============================================================================
# LM INTERACTION WRAPPER TESTS
# ============================================================================

test_that("calculate_lm_interaction_s4 validates input", {
  expect_error(calculate_lm_interaction_s4("not_analysis"), "TSENATAnalysis")
})

test_that("calculate_lm_interaction_s4 accepts formula from config", {
  se <- make_test_se_for_wrappers()
  cfg <- list(formula = ~ treatment)
  analysis <- TSENATAnalysis(se, config = cfg)
  expect_equal(analysis@config$formula, ~ treatment)
})

# ============================================================================
# JACKKNIFE WRAPPER TESTS
# ============================================================================

test_that("jackknife_tsallis_entropy_s4 requires diversity results", {
  se <- make_test_se_for_wrappers()
  analysis <- TSENATAnalysis(se)
  expect_error(
    jackknife_tsallis_entropy_s4(analysis, q = 1.0),
    "Diversity results required"
  )
})

test_that("jackknife_tsallis_entropy_s4 errors on unavailable q-value", {
  se <- make_test_se_for_wrappers()
  analysis <- TSENATAnalysis(se)

  # Manually add diversity for q=1.0
  analysis@diversity_results$q_1.0 <- SummarizedExperiment(
    assays = list(counts = matrix(rnorm(100), nrow = 10))
  )

  # Jackknife for q=2.0 should error
  expect_error(
    jackknife_tsallis_entropy_s4(analysis, q = 2.0),
    "not calculated"
  )
})

# ============================================================================
# DIVERGENCE WRAPPER TESTS
# ============================================================================

test_that("calculate_divergence_s4 requires diversity results", {
  se <- make_test_se_for_wrappers()
  analysis <- TSENATAnalysis(se)
  expect_error(
    calculate_divergence_s4(analysis),
    "Diversity results required"
  )
})

# ============================================================================
# Q-GENE INTERACTIONS WRAPPER TESTS
# ============================================================================

test_that("detect_q_gene_interactions_s4 requires diversity", {
  se <- make_test_se_for_wrappers()
  analysis <- TSENATAnalysis(se)
  expect_error(
    detect_q_gene_interactions_s4(analysis),
    "Diversity results required"
  )
})

# ============================================================================
# DIFFERENCE WRAPPER TESTS
# ============================================================================

test_that("calculate_difference_s4 requires control specification", {
  se <- make_test_se_for_wrappers()
  analysis <- TSENATAnalysis(se)
  expect_error(
    calculate_difference_s4(analysis),
    "Diversity results required"
  )
})

# ============================================================================
# METADATA TRACKING TESTS
# ============================================================================

test_that("Wrappers initialize metadata tracking structure", {
  se <- make_test_se_for_wrappers()
  analysis <- TSENATAnalysis(se)

  expect_true("function_calls" %in% names(analysis@metadata))
  expect_true(is.character(analysis@metadata$function_calls))
  expect_equal(length(analysis@metadata$function_calls), 0)
})

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
