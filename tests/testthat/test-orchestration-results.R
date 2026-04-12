# ============================================================================
# TEST SUITE: orchestration_results.R - Results Accessor Functions
# ============================================================================
# Tests for:
# - results() public function
# - .validate_results_params() parameter validation
# - .get_diversity_q_value() q-value extraction
# - .display_diversity_table() display formatting
# - .extract_result_by_type() result extraction
# - .extract_jackknife_result() jackknife extraction
# - .extract_or_compute_switching_tables() lazy computation
# - Helper functions for result processing
#
# Coverage Targets:
# - Lines with hits=0 in cobertura.xml for orchestration_results.R
# - Error conditions and validation
# - Display and formatting functions
# ============================================================================

# Ensure TSENAT is loaded before making TSENAT::: calls
# This prevents lazy-loading during test execution which can cause
# benign stack imbalance warnings in R's namespace system
library(TSENAT)

# Test helpers
.make_test_se_orchr <- function(n_genes = 5, n_samples = 4) {
  counts <- matrix(rpois(n_genes * n_samples, lambda = 50), nrow = n_genes)
  rownames(counts) <- paste0("Gene_", seq_len(n_genes))
  colnames(counts) <- paste0("S", seq_len(n_samples))
  
  coldata <- data.frame(
    sample = colnames(counts),
    condition = rep(c("A", "B"), length.out = n_samples),
    row.names = colnames(counts)
  )
  
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = coldata
  )
}

.make_test_analysis_orchr <- function() {
  se <- .make_test_se_orchr()
  TSENAT::TSENATAnalysis(se, config = list(q = c(1.0, 2.0)))
}

context("orchestration_results: Parameter Validation")

test_that(".validate_results_params rejects invalid TSENATAnalysis object", {
  expect_error(
    TSENAT:::.validate_results_params("not_analysis", type = "diversity", rankBy = "none", format = "auto", filterFDR = NULL),
    "must be a TSENATAnalysis object",
    fixed = TRUE
  )
})

test_that(".validate_results_params rejects invalid rankBy values", {
  analysis <- .make_test_analysis_orchr()
  
  expect_error(
    TSENAT:::.validate_results_params(analysis, type = "diversity", rankBy = "invalid", format = "auto", filterFDR = NULL),
    "must be one of",
    fixed = TRUE
  )
  
  expect_error(
    TSENAT:::.validate_results_params(analysis, type = "diversity", rankBy = "byScore", format = "auto", filterFDR = NULL),
    "must be one of"
  )
})

test_that(".validate_results_params rejects invalid format values", {
  analysis <- .make_test_analysis_orchr()
  
  expect_error(
    TSENAT:::.validate_results_params(analysis, type = "diversity", rankBy = "none", format = "xlsx", filterFDR = NULL),
    "must be one of"
  )
  
  expect_error(
    TSENAT:::.validate_results_params(analysis, type = "diversity", rankBy = "none", format = "invalid_format", filterFDR = NULL),
    "must be one of"
  )
})

test_that(".validate_results_params rejects invalid filterFDR values", {
  analysis <- .make_test_analysis_orchr()
  
  expect_error(
    TSENAT:::.validate_results_params(analysis, type = "diversity", rankBy = "none", format = "auto", filterFDR = -0.1),
    "must be between 0 and 1"
  )
  
  expect_error(
    TSENAT:::.validate_results_params(analysis, type = "diversity", rankBy = "none", format = "auto", filterFDR = 1.5),
    "must be between 0 and 1"
  )
})

test_that(".validate_results_params accepts valid parameters", {
  analysis <- .make_test_analysis_orchr()
  
  # Should not throw errors
  expect_silent(
    TSENAT:::.validate_results_params(analysis, type = "diversity", rankBy = "pvalue", format = "dataframe", filterFDR = 0.05)
  )
  
  expect_silent(
    TSENAT:::.validate_results_params(analysis, type = "lm", rankBy = "qvalue", format = "matrix", filterFDR = NULL)
  )
  
  expect_silent(
    TSENAT:::.validate_results_params(analysis, type = "assumptions", rankBy = "none", format = "list", filterFDR = 0)
  )
})

# ============================================================================
context("orchestration_results: Get Diversity Q-value")

test_that(".get_diversity_q_value returns result for existing q-value (decimal format)", {
  # Create mock diversity results with various q-value formats
  result <- list(
    q_1 = data.frame(gene = "g1", value = 1),
    q_1.5 = data.frame(gene = "g1", value = 2)
  )
  
  # Should find q_1 with various decimal specifications
  expect_equal(
    TSENAT:::.get_diversity_q_value(result, 1.0)$value,
    1
  )
  
  expect_equal(
    TSENAT:::.get_diversity_q_value(result, 1.5)$value,
    2
  )
})

test_that(".get_diversity_q_value handles integer q-values", {
  result <- list(
    q_0 = data.frame(gene = "g1", value = 0),
    q_1 = data.frame(gene = "g1", value = 1),
    q_2 = data.frame(gene = "g1", value = 2)
  )
  
  # Should find integer q-values
  expect_equal(
    TSENAT:::.get_diversity_q_value(result, 0)$value,
    0
  )
  
  expect_equal(
    TSENAT:::.get_diversity_q_value(result, 2)$value,
    2
  )
})

test_that(".get_diversity_q_value handles underscore format for decimals", {
  result <- list(
    q_0_5 = data.frame(gene = "g1", value = 0.5),
    q_1_5 = data.frame(gene = "g1", value = 1.5)
  )
  
  # Underscore conversion: 1.5 -> q_1_5
  expect_equal(
    TSENAT:::.get_diversity_q_value(result, 1.5)$value,
    1.5
  )
})

test_that(".get_diversity_q_value throws error for missing q-value", {
  result <- list(
    q_1_000 = data.frame(gene = "g1", value = 1),
    q_2_000 = data.frame(gene = "g1", value = 2)
  )
  
  expect_error(
    TSENAT:::.get_diversity_q_value(result, 3.0),
    "not found in results"
  )
  
  expect_error(
    TSENAT:::.get_diversity_q_value(result, 0.5),
    "not found in results"
  )
})

test_that(".get_diversity_q_value shows available q-values in error message", {
  result <- list(
    q_0_500 = data.frame(gene = "g1", value = 0.5),
    q_1_500 = data.frame(gene = "g1", value = 1.5)
  )
  
  expect_error(
    TSENAT:::.get_diversity_q_value(result, 2.5),
    "Available q-values"
  )
})

# ============================================================================
context("orchestration_results: Display Diversity Table")

test_that(".display_diversity_table displays formatted output correctly", {
  # Create test data
  entropy_data <- matrix(seq(1, 20), nrow = 5, ncol = 4)
  rownames(entropy_data) <- paste0("Gene_", 1:5)
  colnames(entropy_data) <- c("S1", "S2", "S3", "S4")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = entropy_data),
    colData = data.frame(sample_id = c("S1", "S2", "S3", "S4"))
  )
  
  analysis <- .make_test_analysis_orchr()
  analysis@diversity_results <- list(
    q_1 = se,
    q_1.5 = se,
    q_2 = se
  )
  
  # Capture messages using testthat's capture_messages
  msgs <- testthat::capture_messages({
    TSENAT:::.display_diversity_table(
      analysis = analysis,
      result = list(q_1 = se),
      q = NULL,
      n_genes = 2,
      q_values_table = c(1.0, 1.5, 2.0)
    )
  })
  
  full_msg <- paste(msgs, collapse = "")
  
  # Should contain expected headers and data
  expect_true(grepl("Tsallis entropy across q-spectrum", full_msg))
  expect_true(grepl("Gene", full_msg))
  expect_true(grepl("Gene_1", full_msg))
})

test_that(".display_diversity_table handles empty results gracefully", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = matrix(numeric(0), nrow = 0, ncol = 0))
  )
  
  analysis <- .make_test_analysis_orchr()
  analysis@diversity_results <- list()
  
  # Should not error with empty diversity results
  expect_silent(
    result <- TSENAT:::.display_diversity_table(
      analysis = analysis,
      result = list(),
      q = NULL,
      n_genes = 4,
      q_values_table = c(1.0, 1.5, 2.0)
    )
  )
})

test_that(".display_diversity_table limits displayed genes", {
  entropy_data <- matrix(seq(1, 100), nrow = 20, ncol = 5)
  rownames(entropy_data) <- paste0("Gene_", 1:20)
  colnames(entropy_data) <- paste0("S", 1:5)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = entropy_data)
  )
  
  analysis <- .make_test_analysis_orchr()
  analysis@diversity_results <- list(q_1 = se)
  
  # Capture messages
  msgs <- testthat::capture_messages({
    TSENAT:::.display_diversity_table(
      analysis = analysis,
      result = list(q_1 = se),
      q = NULL,
      n_genes = 3,
      q_values_table = c(1.0)
    )
  })
  
  full_msg <- paste(msgs, collapse = "")
  
  # Should show only 3 genes
  expect_true(grepl("showing 3", full_msg))
})

# ============================================================================
context("orchestration_results: Extract Result by Type")

test_that(".extract_result_by_type returns NULL for empty results", {
  analysis <- .make_test_analysis_orchr()
  
  expect_null(TSENAT:::.extract_result_by_type(analysis, "diversity"))
  expect_null(TSENAT:::.extract_result_by_type(analysis, "divergence"))
  expect_null(TSENAT:::.extract_result_by_type(analysis, "lm"))
  expect_null(TSENAT:::.extract_result_by_type(analysis, "jackknife"))
})

test_that(".extract_result_by_type extracts diversity results", {
  analysis <- .make_test_analysis_orchr()
  
  # Manually create and add diversity results
  entropy_data <- matrix(seq(1, 8), nrow = 2, ncol = 4)
  rownames(entropy_data) <- c("Gene_1", "Gene_2")
  colnames(entropy_data) <- paste0("S", 1:4)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = entropy_data)
  )
  
  analysis@diversity_results <- list(q_1 = se, q_2 = se)
  
  result <- TSENAT:::.extract_result_by_type(analysis, "diversity")
  
  expect_is(result, "list")
  expect_gt(length(result), 0)
})

test_that(".extract_result_by_type extracts lm results with nested structure", {
  analysis <- .make_test_analysis_orchr()
  
  lm_df <- data.frame(
    gene = "g1",
    p_interaction = 0.01,
    adj_p_interaction = 0.05
  )
  
  # Test with nested lm_interaction
  analysis@lm_results <- list(lm_interaction = lm_df)
  
  result <- TSENAT:::.extract_result_by_type(analysis, "lm")
  expect_identical(result, lm_df)
})

test_that(".extract_result_by_type extracts lm results from flat structure", {
  analysis <- .make_test_analysis_orchr()
  
  lm_df <- data.frame(gene = "g1", p_interaction = 0.01)
  analysis@lm_results <- lm_df
  
  result <- TSENAT:::.extract_result_by_type(analysis, "lm")
  expect_identical(result, lm_df)
})

test_that(".extract_result_by_type extracts rank_test results", {
  analysis <- .make_test_analysis_orchr()
  
  rank_df <- data.frame(gene = "g1", p_value = 0.05)
  analysis@rank_test_results <- list(rank_test = rank_df)
  
  result <- TSENAT:::.extract_result_by_type(analysis, "rank_test")
  expect_identical(result, rank_df)
})

test_that(".extract_result_by_type extracts assumptions from metadata", {
  analysis <- .make_test_analysis_orchr()
  
  assumptions_data <- list(
    exchangeability = data.frame(test = "exchangeability"),
    monotonicity = data.frame(test = "monotonicity")
  )
  
  # Directly assign to @metadata slot
  meta <- analysis@metadata
  meta$rankbased_assumptions <- list(result = structure(assumptions_data, checks = TRUE))
  analysis@metadata <- meta
  
  result <- TSENAT:::.extract_result_by_type(analysis, "assumptions")
  
  expect_equal(result, TRUE)
})

test_that(".extract_result_by_type returns NULL for unknown type", {
  analysis <- .make_test_analysis_orchr()
  
  expect_error(
    TSENAT:::.extract_result_by_type(analysis, "unknown_type"),
    "Unknown result type"
  )
})

# ============================================================================
context("orchestration_results: Extract Jackknife Result")

test_that(".extract_jackknife_result returns data.frame directly", {
  analysis <- .make_test_analysis_orchr()
  
  jk_df <- data.frame(
    gene = "g1",
    mean = 0.5,
    ci_lower = 0.3,
    ci_upper = 0.7
  )
  
  analysis@jackknife_results <- jk_df
  
  result <- TSENAT:::.extract_jackknife_result(analysis)
  expect_identical(result, jk_df)
})

test_that(".extract_jackknife_result extracts from nested 'results'", {
  analysis <- .make_test_analysis_orchr()
  
  jk_df <- data.frame(gene = "g1", mean = 0.5)
  analysis@jackknife_results <- list(results = jk_df)
  
  result <- TSENAT:::.extract_jackknife_result(analysis)
  expect_identical(result, jk_df)
})

test_that(".extract_jackknife_result extracts from nested 'ci'", {
  analysis <- .make_test_analysis_orchr()
  
  ci_df <- data.frame(gene = "g1", ci_lower = 0.3, ci_upper = 0.7)
  analysis@jackknife_results <- list(ci = ci_df)
  
  result <- TSENAT:::.extract_jackknife_result(analysis)
  expect_identical(result, ci_df)
})

test_that(".extract_jackknife_result returns NULL for empty slots", {
  analysis <- .make_test_analysis_orchr()
  analysis@jackknife_results <- list()
  
  result <- TSENAT:::.extract_jackknife_result(analysis)
  expect_null(result)
})

# ============================================================================
context("orchestration_results: Extract/Compute Switching Tables")

test_that(".extract_or_compute_switching_tables returns cached tables if present", {
  analysis <- .make_test_analysis_orchr()
  
  cached_tables <- list(
    gene1 = data.frame(comparison = "A_vs_B", switched = TRUE),
    gene2 = data.frame(comparison = "A_vs_B", switched = FALSE)
  )
  
  # Directly assign to @metadata slot
  meta <- analysis@metadata
  meta$switching_tables <- cached_tables
  analysis@metadata <- meta
  
  result <- TSENAT:::.extract_or_compute_switching_tables(analysis)
  expect_identical(result, cached_tables)
})

test_that(".extract_or_compute_switching_tables returns NULL with missing prerequisites", {
  analysis <- .make_test_analysis_orchr()
  analysis@lm_results <- list()
  analysis@jackknife_results <- list()
  
  result <- TSENAT:::.extract_or_compute_switching_tables(analysis)
  expect_null(result)
})

# ============================================================================
context("orchestration_results: Warn Unsupported Parameters")

test_that(".warn_unsupported_params warns for switching_tables with rankBy", {
  expect_warning(
    TSENAT:::.warn_unsupported_params(type = "switching_tables", filterFDR = NULL, rankBy = "pvalue"),
    "rankBy and filterFDR are not supported"
  )
})

test_that(".warn_unsupported_params warns for switching_tables with filterFDR", {
  expect_warning(
    TSENAT:::.warn_unsupported_params(type = "switching_tables", filterFDR = 0.05, rankBy = "none"),
    "rankBy and filterFDR are not supported"
  )
})

test_that(".warn_unsupported_params warns for invalid rankBy usage", {
  expect_warning(
    TSENAT:::.warn_unsupported_params(type = "diversity", filterFDR = NULL, rankBy = "pvalue"),
    "rankBy.*not supported"
  )
})

test_that(".warn_unsupported_params does not warn for valid rankBy on lm", {
  expect_silent(
    TSENAT:::.warn_unsupported_params(type = "lm", filterFDR = NULL, rankBy = "pvalue")
  )
})

test_that(".warn_unsupported_params does not warn for valid rankBy on jackknife", {
  expect_silent(
    TSENAT:::.warn_unsupported_params(type = "jackknife", filterFDR = NULL, rankBy = "qvalue")
  )
})

# ============================================================================
context("orchestration_results: Filter Statistical Results by FDR")

test_that(".filter_statistical_by_fdr filters lm results correctly", {
  result <- data.frame(
    gene = c("g1", "g2", "g3", "g4"),
    adj_p_interaction = c(0.001, 0.05, 0.1, 0.2),
    stringsAsFactors = FALSE
  )
  
  filtered <- TSENAT:::.filter_statistical_by_fdr(result, type = "lm", filterFDR = 0.05)
  
  expect_equal(nrow(filtered), 2)
  expect_equal(filtered$gene, c("g1", "g2"))
})

test_that(".filter_statistical_by_fdr filters rank_test results correctly", {
  result <- data.frame(
    gene = c("g1", "g2", "g3"),
    adj_p_value = c(0.001, 0.06, 0.1),
    stringsAsFactors = FALSE
  )
  
  filtered <- TSENAT:::.filter_statistical_by_fdr(result, type = "rank_test", filterFDR = 0.05)
  
  expect_equal(nrow(filtered), 1)
  expect_equal(filtered$gene, "g1")
})

test_that(".filter_statistical_by_fdr returns full results when filterFDR is NULL", {
  result <- data.frame(
    gene = c("g1", "g2"),
    adj_p_interaction = c(0.001, 0.1)
  )
  
  filtered <- TSENAT:::.filter_statistical_by_fdr(result, type = "lm", filterFDR = NULL)
  
  expect_identical(filtered, result)
})

test_that(".filter_statistical_by_fdr returns full results for non-dataframe", {
  result <- list(something = "else")
  
  filtered <- TSENAT:::.filter_statistical_by_fdr(result, type = "lm", filterFDR = 0.05)
  
  expect_identical(filtered, result)
})

# ============================================================================
context("orchestration_results: Rank Statistical Results")

test_that(".rank_statistical_results returns unchanged when rankBy='none'", {
  result <- data.frame(
    gene = c("g3", "g1", "g2"),
    p_value = c(0.5, 0.001, 0.05)
  )
  
  ranked <- TSENAT:::.rank_statistical_results(result, type = "rank_test", rankBy = "none", n = NA)
  
  expect_identical(ranked, result)
})

test_that(".rank_statistical_results ranks by pvalue", {
  result <- data.frame(
    gene = c("g3", "g1", "g2"),
    p_interaction = c(0.5, 0.001, 0.05),
    stringsAsFactors = FALSE
  )
  
  ranked <- TSENAT:::.rank_statistical_results(result, type = "lm", rankBy = "pvalue", n = NA)
  
  expect_equal(ranked$gene[1], "g1")
  expect_equal(ranked$gene[2], "g2")
  expect_equal(ranked$gene[3], "g3")
})

test_that(".rank_statistical_results limits to top n", {
  result <- data.frame(
    gene = c("g1", "g2", "g3", "g4", "g5"),
    p_value = c(0.001, 0.002, 0.003, 0.004, 0.005),
    stringsAsFactors = FALSE
  )
  
  ranked <- TSENAT:::.rank_statistical_results(result, type = "rank_test", rankBy = "pvalue", n = 2)
  
  expect_equal(nrow(ranked), 2)
  expect_equal(ranked$gene, c("g1", "g2"))
})

# ============================================================================
context("orchestration_results: Main Results Accessor - Error Handling")

test_that("results() throws error for invalid analysis object", {
  expect_error(
    TSENAT::results("not_an_analysis", type = "diversity"),
    "must be a TSENATAnalysis object"
  )
})

test_that("results() throws error for invalid type", {
  analysis <- .make_test_analysis_orchr()
  
  expect_error(
    TSENAT::results(analysis, type = "invalid_type"),
    "Unknown result type"
  )
})

test_that("results() returns NULL for non-computed results", {
  analysis <- .make_test_analysis_orchr()
  
  result <- TSENAT::results(analysis, type = "divergence")
  expect_null(result)
})

# ============================================================================
context("orchestration_results: Main Results Accessor - Diversity Processing")

test_that("results() with display_table=TRUE produces output", {
  entropy_data <- matrix(seq(1, 12), nrow = 3, ncol = 4)
  rownames(entropy_data) <- c("Gene_1", "Gene_2", "Gene_3")
  colnames(entropy_data) <- paste0("S", 1:4)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = entropy_data)
  )
  
  analysis <- .make_test_analysis_orchr()
  analysis@diversity_results <- list(q_1 = se)
  
  # Capture messages while calling results()
  msgs <- testthat::capture_messages({
    result <- TSENAT::results(analysis, type = "diversity", display_table = TRUE, n_genes = 2)
  })
  
  full_msg <- paste(msgs, collapse = "")
  
  # Check that message was produced
  expect_true(grepl("Tsallis entropy", full_msg))
  # Result should be a list when q is not specified (returns all q-values as list)
  expect_is(result, "list")
})

test_that("results() extracts specific q-value from diversity", {
  entropy_data <- matrix(seq(1, 8), nrow = 2, ncol = 4)
  rownames(entropy_data) <- c("Gene_1", "Gene_2")
  colnames(entropy_data) <- paste0("S", 1:4)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = entropy_data)
  )
  
  analysis <- .make_test_analysis_orchr()
  analysis@diversity_results <- list(
    q_1 = se,
    q_2 = se
  )
  
  # Request specific q-value
  result <- TSENAT::results(analysis, type = "diversity", q = 1.0)
  
  expect_is(result, "SummarizedExperiment")
  expect_equal(nrow(result), 2)
})

# ============================================================================
context("orchestration_results: Main Results Accessor - Statistical Results")

test_that("results() applies rankBy to lm results", {
  analysis <- .make_test_analysis_orchr()
  
  lm_df <- data.frame(
    gene = c("g3", "g1", "g2"),
    p_interaction = c(0.5, 0.001, 0.05),
    adj_p_interaction = c(0.6, 0.01, 0.1),
    stringsAsFactors = FALSE
  )
  
  analysis@lm_results <- list(lm_interaction = lm_df)
  
  # Get ranked results
  result <- TSENAT::results(analysis, type = "lm", rankBy = "pvalue")
  
  expect_equal(result$gene[1], "g1")
  expect_equal(nrow(result), 3)
})

test_that("results() applies n parameter to limit results", {
  analysis <- .make_test_analysis_orchr()
  
  lm_df <- data.frame(
    gene = c("g1", "g2", "g3", "g4", "g5"),
    p_interaction = c(0.005, 0.01, 0.02, 0.1, 0.5),
    stringsAsFactors = FALSE
  )
  
  analysis@lm_results <- list(lm_interaction = lm_df)
  
  result <- TSENAT::results(analysis, type = "lm", rankBy = "pvalue", n = 2)
  
  expect_equal(nrow(result), 2)
  expect_equal(result$gene, c("g1", "g2"))
})

test_that("results() applies FDR filter", {
  analysis <- .make_test_analysis_orchr()
  
  rank_df <- data.frame(
    gene = c("g1", "g2", "g3", "g4"),
    adj_p_value = c(0.001, 0.05, 0.1, 0.5),
    stringsAsFactors = FALSE
  )
  
  analysis@rank_test_results <- list(rank_test = rank_df)
  
  result <- TSENAT::results(analysis, type = "rank_test", filterFDR = 0.05)
  
  expect_equal(nrow(result), 2)
  expect_equal(result$gene, c("g1", "g2"))
})

# ============================================================================
context("orchestration_results: Lazy Computation of Switching Tables")

test_that("results() with type='switching_tables' returns cached if available", {
  analysis <- .make_test_analysis_orchr()
  
  cached <- list(gene1 = data.frame(test = 1))
  
  # Directly assign to @metadata slot
  meta <- analysis@metadata
  meta$switching_tables <- cached
  analysis@metadata <- meta
  
  result <- TSENAT::results(analysis, type = "switching_tables")
  
  expect_identical(result, cached)
})

test_that("results() warns about unsupported rankBy for switching_tables", {
  analysis <- .make_test_analysis_orchr()
  
  # Directly assign to @metadata slot
  meta <- analysis@metadata
  meta$switching_tables <- list(test = 1)
  analysis@metadata <- meta
  
  expect_warning(
    TSENAT::results(analysis, type = "switching_tables", rankBy = "pvalue"),
    "not supported"
  )
})
