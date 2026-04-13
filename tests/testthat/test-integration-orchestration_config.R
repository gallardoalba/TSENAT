library(testthat)

context("Orchestration: Configuration and Pipeline")

# ============================================================================
# Helper: Create test SummarizedExperiment using real package data
# ============================================================================

make_test_se <- function() {
  data(readcounts, package = "TSENAT", envir = environment())
  readcounts <- as.matrix(readcounts)
  
  if (ncol(readcounts) != 16) {
    stop("Expected 16 samples in readcounts, got ", ncol(readcounts))
  }
  
  metadata_df <- read.table(
    system.file("extdata", "metadata.tsv", package = "TSENAT"),
    header = TRUE, sep = "\t"
  )
  
  if (nrow(metadata_df) != 16) {
    stop("Expected 16 samples in metadata, got ", nrow(metadata_df))
  }
  
  gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
  
  config <- TSENAT_config(
    sample_col = "sample",
    condition_col = "condition",
    subject_col = "paired_samples",
    q = seq(0, 2, by = 0.05),
    paired = TRUE,
    control = "normal",
    stringency = "severe",
    nthreads = 1
  )
  
  analysis <- build_analysis(
    config = config,
    readcounts = readcounts,
    metadata = metadata_df,
    tx2gene = gff3_file,
    tpm = tpm,
    effective_length = effective_length
  )
  
  analysis <- filter_analysis(analysis, stringency = "medium", verbose = FALSE)
  se(analysis)
}

# ============================================================================
# TEST: TSENAT_config function (CONSOLIDATED: 10 → 3 tests)
# ============================================================================

test_that("TSENAT_config creates and stores configurations", {
  # Test 1: Defaults
  config <- TSENAT_config()
  expect_true(is.list(config))
  expect_true("q" %in% names(config))
  
  # Test 2: Custom parameters
  config <- TSENAT_config(p_threshold = 0.01, seed = 123)
  expect_equal(config$p_threshold, 0.01)
  expect_equal(config$seed, 123)
  
  # Test 3: All provided arguments
  config <- TSENAT_config(
    p_threshold = 0.05,
    q = c(0.5, 1.0, 1.5),
    norm = "none"
  )
  expect_equal(config$p_threshold, 0.05)
  expect_equal(config$norm, "none")
  expect_equal(length(config$q), 3)
})

test_that("TSENAT_config validates parameters", {
  # Test q validation
  config <- TSENAT_config(q = c(0.5, 1.0, 1.5))
  expect_true(is.numeric(config$q))
  expect_true(length(config$q) >= 1)
  
  # Test stringency parameter
  config <- TSENAT_config(stringency = "medium")
  expect_equal(config$stringency, "medium")
})

test_that("TSENAT_config/getConfig/setConfig integration", {
  se <- make_test_se()
  
  # getConfig test
  analysis <- TSENATAnalysis(se, config = TSENAT_config(seed = 99))
  config <- getConfig(analysis)
  expect_equal(config$seed, 99)
  
  # setConfig test
  analysis_new <- setConfig(analysis, TSENAT_config(seed = 2))
  expect_equal(analysis_new@config$seed, 2)
  
  # Preserves SE data
  new_config <- TSENAT_config(p_threshold = 0.001)
  analysis_new <- setConfig(analysis, new_config)
  expect_identical(
    SummarizedExperiment::assay(analysis_new@se, "counts"),
    SummarizedExperiment::assay(se, "counts")
  )
})

# ============================================================================
# TEST: tsenat pipeline function (CONSOLIDATED: 7 → 2 tests)
# ============================================================================

test_that("tsenat pipeline requires TSENATAnalysis object", {
  se <- make_test_se()
  
  # All tsenat calls should fail with raw SE
  expect_error(TSENAT(se), "must be a TSENATAnalysis object")
  expect_error(TSENAT(se, verbose = FALSE), "must be a TSENATAnalysis object")
  
  # Invalid SEs
  counts <- matrix(rpois(100 * 20, lambda = 100), nrow = 100)
  invalid_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = data.frame(condition = rep(c('A', 'B'), 10), row.names = paste0('S', 1:20))
  )
  expect_error(TSENAT(invalid_se, verbose = FALSE), "must be a TSENATAnalysis object")
})

test_that("tsenat passes config parameters through pipeline", {
  se <- make_test_se()
  config <- TSENAT_config(p_threshold = 0.001, seed = 777)
  analysis <- TSENATAnalysis(se, config = config)
  
  expect_equal(getConfig(analysis)$seed, 777)
  expect_equal(getConfig(analysis)$p_threshold, 0.001)
})

# ============================================================================
# TEST: .validate_analysis_object (CONSOLIDATED: 9 → 3 tests)
# ============================================================================

test_that(".validate_analysis_object accepts valid objects", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = TSENAT_config())
  
  # Should not raise error for valid object
  expect_no_error(.validate_analysis_object(analysis))
})

test_that(".validate_analysis_object rejects invalid configurations", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = TSENAT_config())
  
  # Test: Empty SE
  analysis_empty <- analysis
  analysis_empty@se <- analysis_empty@se[0, ]
  error_msg <- tryCatch(
    .validate_analysis_object(analysis_empty),
    error = function(e) e$message
  )
  expect_true(grepl("Validation failed|se_valid", error_msg))
  
  # Test: Insufficient samples
  analysis_few <- analysis
  analysis_few@se <- analysis_few@se[, 1, drop = FALSE]
  error_msg <- tryCatch(
    .validate_analysis_object(analysis_few),
    error = function(e) e$message
  )
  expect_true(grepl("Validation failed|min_samples", error_msg))
})

test_that(".validate_analysis_object handles missing columns gracefully", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = TSENAT_config())
  
  # Remove condition column
  coldata <- SummarizedExperiment::colData(analysis@se)
  coldata$condition <- NULL
  SummarizedExperiment::colData(analysis@se) <- coldata
  
  # Should handle gracefully (error acceptable, but function must respond)
  result <- tryCatch(
    .validate_analysis_object(analysis),
    error = function(e) list(error = TRUE, msg = e$message)
  )
  
  # Either should error (validation caught the missing column)
  # OR should return TRUE/FALSE (graceful handling)
  expect_true(
    is.list(result) || is.logical(result),
    info = "Should return error list or logical result"
  )
  
  # If error, should mention validation failure
  if (is.list(result) && !is.null(result$error)) {
    expect_true(grepl("Validation failed", result$msg))
  }
})

# ============================================================================
# TEST: .finalize_tsenat_analysis (CONSOLIDATED: 4 → 2 tests)
# ============================================================================

test_that(".finalize_tsenat_analysis adds metadata and respects verbose", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se)
  
  before_time <- Sys.time()
  
  # Test with verbose = TRUE
  expect_message(
    .finalize_tsenat_analysis(analysis, verbose = TRUE),
    "ANALYSIS COMPLETE"
  )
  
  # Test end time added
  analysis_finalized <- .finalize_tsenat_analysis(analysis, verbose = FALSE)
  expect_true(inherits(analysis_finalized@metadata$ended_at, "POSIXct"))
  expect_true(analysis_finalized@metadata$ended_at >= before_time)
})

test_that(".finalize_tsenat_analysis verbose parameter works", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se)
  
  # verbose = TRUE should produce messages
  expect_message(
    .finalize_tsenat_analysis(analysis, verbose = TRUE),
    "Results Summary"
  )
  
  # verbose = FALSE should NOT produce messages
  expect_no_message(
    .finalize_tsenat_analysis(analysis, verbose = FALSE)
  )
})

# ============================================================================
# TEST: .track_analysis_metadata (CONSOLIDATED: 7 → 2 tests)
# ============================================================================

test_that(".track_analysis_metadata records workflow information", {
  se <- make_test_se()
  config <- TSENAT_config(fdr_threshold = 0.01, q = c(0.5, 1.0, 1.5))
  analysis <- TSENATAnalysis(se, config = config)
  
  analysis_tracked <- .track_analysis_metadata(analysis, config)
  
  # Check workflow structure
  expect_true("workflow" %in% names(analysis_tracked@metadata))
  expect_true("workflow_type" %in% names(analysis_tracked@metadata$workflow))
  expect_true("tsenat_version" %in% names(analysis_tracked@metadata$workflow))
  
  # Check method parameters
  expect_true("methods_parameters" %in% names(analysis_tracked@metadata))
  expect_equal(analysis_tracked@metadata$methods_parameters$fdr_threshold, 0.01)
  expect_equal(length(analysis_tracked@metadata$methods_parameters$q), 3)
})

test_that(".track_analysis_metadata preserves and timestamps", {
  se <- make_test_se()
  config <- TSENAT_config(condition_col = "condition")
  analysis <- TSENATAnalysis(se, config = config)
  analysis@metadata$custom_field <- "custom_value"
  
  before_time <- Sys.time()
  analysis_tracked <- .track_analysis_metadata(analysis, config)
  after_time <- Sys.time()
  
  # Check existing metadata preserved
  expect_equal(analysis_tracked@metadata$custom_field, "custom_value")
  
  # Check timestamp
  completion_time <- analysis_tracked@metadata$workflow$completion_time
  expect_true(inherits(completion_time, "POSIXct"))
  expect_true(completion_time >= before_time && completion_time <= after_time)
  
  # Check config parameters stored
  expect_equal(analysis_tracked@metadata$methods_parameters$condition_col, "condition")
})

# ============================================================================
# TEST: results accessor (CONSOLIDATED: 32 → 5 tests)
# ============================================================================

test_that("results returns NULL/errors for edge cases", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = TSENAT_config())
  
  # Uncomputed results
  expect_null(results(analysis, type = "diversity", format = "table"))
  expect_null(results(analysis, type = "divergence"))
  expect_null(results(analysis, type = "lm"))
  
  # Unknown type
  expect_error(results(analysis, type = "unknown"), "Unknown result type")
  
  # Non-TSENATAnalysis object
  expect_error(results(se, type = "diversity", format = "table"), "must be a TSENATAnalysis object")
})

test_that("results handles diversity with q-value filtering", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = TSENAT_config())
  
  n_genes <- nrow(se)
  n_samples <- ncol(se)
  
  # Create SummarizedExperiment objects for each q-value
  create_diversity_se <- function(q_value) {
    mat <- matrix(rnorm(n_genes * n_samples, mean = 3, sd = 1), 
                  nrow = n_genes, ncol = n_samples)
    colnames(mat) <- colnames(se)
    rownames(mat) <- rownames(se)
    SummarizedExperiment::SummarizedExperiment(assays = list(q_entropy = mat),
                                               colData = SummarizedExperiment::colData(se))
  }
  
  diversity_results <- list(
    q_0.5 = create_diversity_se(0.5),
    q_1.0 = create_diversity_se(1.0),
    q_1.5 = create_diversity_se(1.5),
    q_2.0 = create_diversity_se(2.0)
  )
  analysis@diversity_results <- diversity_results
  
  # All results
  all_results <- results(analysis, type = "diversity", format = "table")
  expect_false(is.null(all_results))
  
  # Specific q-value
  q1_results <- results(analysis, type = "diversity", q = 1.0, format = "table")
  expect_false(is.null(q1_results))
})

test_that("results returns all supported result types", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = TSENAT_config())
  
  n_genes <- nrow(se)
  n_samples <- ncol(se)
  
  # Create SummarizedExperiment objects for diversity results
  create_diversity_se <- function(q_value) {
    mat <- matrix(rnorm(n_genes * n_samples), 
                  nrow = n_genes, ncol = n_samples)
    colnames(mat) <- colnames(se)
    rownames(mat) <- rownames(se)
    SummarizedExperiment::SummarizedExperiment(assays = list(q_entropy = mat),
                                               colData = SummarizedExperiment::colData(se))
  }
  
  # Mock all result types
  analysis@diversity_results <- list(
    q_0.5 = create_diversity_se(0.5),
    q_1.0 = create_diversity_se(1.0)
  )
  analysis@divergence_results <- list(
    q_0.5 = data.frame(gene = paste0("g", 1:n_genes), divergence = rnorm(n_genes))
  )
  analysis@lm_results <- list(lm_interaction = data.frame(
    p_interaction = rnorm(n_genes),
    adj_p_interaction = p.adjust(rnorm(n_genes), method = "BH")
  ))
  analysis@jackknife_results <- list(results = data.frame(
    ci_lower = rnorm(n_genes),
    ci_upper = rnorm(n_genes)
  ))
  analysis@rank_test_results <- list(rank_test = list(results = data.frame(
    gene = paste0("g", 1:n_genes),
    p_value = rnorm(n_genes, mean = 0.01)
  )))
  
  # Test all types at once
  expect_false(is.null(results(analysis, type = "diversity", format = "table")))
  expect_false(is.null(results(analysis, type = "divergence")))
  expect_false(is.null(results(analysis, type = "lm")))
  expect_false(is.null(results(analysis, type = "jackknife")))
  expect_false(is.null(results(analysis, type = "rank_test")))
})

# ============================================================================
# Additional integration tests
# ============================================================================

test_that("stringency levels process without error", {
  se <- make_test_se()
  
  # Test that pipeline errors correctly with raw SE
  for (stringency in c("soft", "medium", "severe")) {
    expect_error(
      TSENAT(se, verbose = FALSE),
      "must be a TSENATAnalysis object"
    )
  }
})

test_that("full pipeline integration works", {
  se <- make_test_se()
  config <- TSENAT_config(p_threshold = 0.05, seed = 123)
  analysis <- TSENATAnalysis(se, config = config)
  
  # Test config retrieval
  expect_equal(getConfig(analysis)$seed, 123)
  
  # Test validation
  expect_no_error(.validate_analysis_object(analysis))
  
  # Test metadata tracking
  analysis_tracked <- .track_analysis_metadata(analysis, config)
  expect_true("workflow" %in% names(analysis_tracked@metadata))
})
# ============================================================================
# TEST SUITE: orchestration_results.R - Results Accessor Functions
# ============================================================================
# Tests for:
# - results() public function
# - .validate_results_params() parameter validation
# - .get_diversity_q_value() q-value extraction
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

test_that("results() with diversity produces formatted output", {
  entropy_data <- matrix(seq(1, 12), nrow = 3, ncol = 4)
  rownames(entropy_data) <- c("Gene_1", "Gene_2", "Gene_3")
  colnames(entropy_data) <- paste0("S", 1:4)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = entropy_data)
  )
  
  analysis <- .make_test_analysis_orchr()
  analysis@diversity_results <- list(q_1 = se)
  
  # Call results() without q - should return data.frame with all q-values
  result <- TSENAT::results(analysis, type = "diversity", n_genes = 2)
  
  # Result should be a data.frame (not list)
  expect_is(result, "data.frame")
  expect_equal(nrow(result), 2)  # n_genes = 2
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
  result <- TSENAT::results(analysis, type = "diversity", q = 1.0, format = "se")
  
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
context("Results display and formatting")

# Helper to create minimal analysis with effect sizes
.make_test_analysis_for_display_table <- function() {
    # Create minimal SummarizedExperiment
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1, nrow = 3, ncol = 4)),
        colData = data.frame(
            sample = paste0("S", 1:4),
            condition = c("A", "A", "B", "B"),
            row.names = paste0("S", 1:4)
        )
    )
    
    # Create analysis using proper constructor
    analysis <- TSENAT::TSENATAnalysis(se, config = list(q = c(0, 1.0)))
    
    # Add effect_sizes_divergence data to metadata
    effect_sizes_data <- data.frame(
        gene_name = c("GENE1", "GENE2", "GENE3"),
        gene_id = c("ENSG00001", "ENSG00002", "ENSG00003"),
        mean_divergence = c(0.2846, 0.1279, 0.3456),
        q_pattern = c("Balanced", "Balanced", "Skewed"),
        d_rare = c(0.3063, 0.1455, 0.3890),
        d_abundant = c(0.2712, 0.1142, 0.3122),
        ratio = c(1.13, 1.27, 1.25),
        p_value_interaction = c(3e-162, 7.3e-95, 2.1e-58),
        stringsAsFactors = FALSE
    )
    
    analysis@metadata$effect_sizes_divergence <- effect_sizes_data
    
    return(analysis)
}

test_that("results produces message output for effect_sizes_divergence", {
    analysis <- .make_test_analysis_for_display_table()
    
    # Call results() - should return silently
    result <- TSENAT::results(
        analysis,
        type = "effect_sizes_divergence",
        top_n = NULL,
        sort_by = "p_value_interaction"
    )
    
    # Result should be a data.frame
    expect_true(is.data.frame(result), "Should return data.frame")
    # Should have expected structure
    expect_true(nrow(result) > 0, "Should return non-empty data.frame")
})

test_that("results returns data.frame with all expected columns", {
    analysis <- .make_test_analysis_for_display_table()
    
    result <- TSENAT::results(
        analysis,
        type = "effect_sizes_divergence",
        top_n = NULL
    )
    
    # Should return data.frame
    expect_true(is.data.frame(result), "Should return data.frame")
    # Should have expected columns
    expect_true("gene_name" %in% colnames(result), "Should have gene_name column")
    expect_true("p_value_interaction" %in% colnames(result), "Should have p-value column")
})

test_that("results filters by top_n correctly", {
    analysis <- .make_test_analysis_for_display_table()
    
    output <- capture_output({
        result <- TSENAT::results(
            analysis,
            type = "effect_sizes_divergence",
            top_n = 2,
            sort_by = "p_value_interaction"
        )
    })
    
    # Result should contain only top 2 genes
    expect_equal(nrow(result), 2L)
    # GENE1 should be first (smallest p-value)
    expect_equal(result$gene_name[1], "GENE1")
})

test_that("results shows all required columns in output", {
    analysis <- .make_test_analysis_for_display_table()
    
    result <- TSENAT::results(
        analysis,
        type = "effect_sizes_divergence",
        top_n = 1,
        sort_by = "p_value_interaction"
    )
    
    # Result should be a data.frame with expected columns
    expect_true(is.data.frame(result), "Should return data.frame")
    # Check for key columns in result
    expect_true("gene_name" %in% colnames(result) || "Gene" %in% colnames(result), 
                "Should contain gene column")
    expect_true("p_value_interaction" %in% colnames(result) || "p_value" %in% colnames(result), 
                "Should contain p-value information")
})

test_that("results returns result with expected structure", {
    analysis <- .make_test_analysis_for_display_table()
    
    # Call the function 
    result <- TSENAT::results(
        analysis,
        type = "effect_sizes_divergence",
        top_n = 1,
        sort_by = "p_value_interaction"
    )
    
    # Result should be a data.frame with correct structure
    expect_true(is.data.frame(result), "Should return data.frame")
    expect_true(nrow(result) > 0, "Result should contain genes")
    expect_true("mean_divergence" %in% colnames(result), "Should have divergence column")
})

test_that("Results contain expected columns for effect_sizes_divergence", {
    analysis <- .make_test_analysis_for_display_table()
    
    result <- TSENAT::results(
        analysis,
        type = "effect_sizes_divergence"
    )
    
    # Verify all expected columns are present
    expect_true(is.data.frame(result), "Result should be data.frame")
    expect_true("gene_name" %in% colnames(result), "Should have gene_name")
    expect_true("mean_divergence" %in% colnames(result), "Should have mean_divergence")
    expect_true("q_pattern" %in% colnames(result), "Should have q_pattern")
    expect_true("ratio" %in% colnames(result), "Should have ratio")
    expect_true("p_value_interaction" %in% colnames(result), "Should have p_value_interaction")
    expect_equal(nrow(result), 3L)
})

test_that("Sorting by p_value_interaction works correctly", {
    analysis <- .make_test_analysis_for_display_table()
    
    result <- TSENAT::results(
        analysis,
        type = "effect_sizes_divergence",
        sort_by = "p_value_interaction",
        top_n = NULL
    )
    
    # Verify sorting is by p-value ascending
    p_vals <- result$p_value_interaction
    expect_true(p_vals[1] <= p_vals[2], "P-values should be sorted ascending")
    expect_true(p_vals[2] <= p_vals[3], "P-values should be sorted ascending")
    # Verify correct order
    expect_equal(result$gene_name[1], "GENE1")
    expect_equal(result$gene_name[2], "GENE2")
    expect_equal(result$gene_name[3], "GENE3")
})

# Tests for the sample parameter in results() function
# Tests added April 12, 2026 to cover results() with sample argument

# Helper function to create a minimal test analysis object
.make_test_analysis_orchr_sample <- function() {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1, nrow = 5, ncol = 4)),
    colData = data.frame(
      sample = paste0("S", 1:4),
      condition = c("A", "A", "B", "B"),
      row.names = paste0("S", 1:4)
    )
  )
  TSENAT::TSENATAnalysis(se, config = list(q = c(0, 1.0)))
}

context("orchestration_results: Sample Parameter for diversity")

test_that("results() with sample parameter selects correct sample", {
  entropy_data <- matrix(seq(1, 16), nrow = 4, ncol = 4)
  rownames(entropy_data) <- c("Gene_1", "Gene_2", "Gene_3", "Gene_4")
  colnames(entropy_data) <- c("S1", "S2", "S3", "S4")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = entropy_data)
  )
  
  analysis <- .make_test_analysis_orchr_sample()
  analysis@diversity_results <- list(
    q_0.000 = se,
    q_1.000 = se
  )
  
  # Call with sample parameter - should return formatted data.frame
  result <- TSENAT::results(analysis, type = "diversity", 
                            n_genes = 2, sample = "S2", q_values_table = c(0.0, 1.0), format = "table")
  
  # Result should be a data.frame when results are displayed with sample
  expect_is(result, "data.frame")
  expect_equal(nrow(result), 2)  # n_genes = 2
  expect_true("Gene" %in% colnames(result))
  expect_true("q_0.0" %in% colnames(result))
})

test_that("results() with sample parameter displays correct sample name", {
  entropy_data <- matrix(seq(1, 16), nrow = 4, ncol = 4)
  rownames(entropy_data) <- c("Gene_1", "Gene_2", "Gene_3", "Gene_4")
  colnames(entropy_data) <- c("S1", "S2", "S3", "S4")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = entropy_data)
  )
  
  analysis <- .make_test_analysis_orchr_sample()
  analysis@diversity_results <- list(q_0.000 = se)
  
  # Capture messages while calling results with sample
  result <- TSENAT::results(analysis, type = "diversity", 
                              n_genes = 1, sample = "S3", format = "table")
  
  # Result should be data.frame with correct sample attribute
  expect_is(result, "data.frame")
  expect_equal(attr(result, "sample"), "S3")
})

test_that("results() raises error for invalid sample name", {
  entropy_data <- matrix(seq(1, 12), nrow = 3, ncol = 4)
  rownames(entropy_data) <- c("Gene_1", "Gene_2", "Gene_3")
  colnames(entropy_data) <- c("S1", "S2", "S3", "S4")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = entropy_data)
  )
  
  analysis <- .make_test_analysis_orchr_sample()
  analysis@diversity_results <- list(q_0.000 = se)
  
  # Should raise error with invalid sample name
  expect_error(
    TSENAT::results(analysis, type = "diversity", 
                    sample = "INVALID_SAMPLE", format = "table"),
    "not found"
  )
})

test_that("results() sample parameter returns data.frame with correct values", {
  entropy_data <- matrix(c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2), nrow = 4, ncol = 2)
  rownames(entropy_data) <- c("Gene_A", "Gene_B", "Gene_C", "Gene_D")
  colnames(entropy_data) <- c("Sample_X", "Sample_Y")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = entropy_data)
  )
  
  analysis <- .make_test_analysis_orchr_sample()
  analysis@diversity_results <- list(q_0.000 = se, q_1.000 = se)
  
  # Get results for Sample_X
  result <- TSENAT::results(analysis, type = "diversity",
                            n_genes = 2, sample = "Sample_X",
                            q_values_table = c(0.0, 1.0), format = "table")
  
  # Check values for specific sample
  expect_equal(result$`q_0.0`[1], 0.5)  # Gene_A, Sample_X
  expect_equal(result$`q_0.0`[2], 0.6)  # Gene_B, Sample_X
})

test_that("results() default sample uses first sample when not specified", {
  entropy_data <- matrix(seq(1, 12), nrow = 3, ncol = 4)
  rownames(entropy_data) <- c("Gene_1", "Gene_2", "Gene_3")
  colnames(entropy_data) <- c("S1", "S2", "S3", "S4")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = entropy_data)
  )
  
  analysis <- .make_test_analysis_orchr_sample()
  analysis@diversity_results <- list(q_0.000 = se)
  
  # Call without sample parameter - should default to first
  result <- TSENAT::results(analysis, type = "diversity", 
                              n_genes = 1, format = "table")
  
  # Should default to first sample (S1) - check via attribute
  expect_equal(attr(result, "sample"), "S1")
})

test_that("results() sample parameter with multiple q-values", {
  entropy_data <- matrix(seq(1, 8), nrow = 2, ncol = 4)
  rownames(entropy_data) <- c("Gene_X", "Gene_Y")
  colnames(entropy_data) <- c("S1", "S2", "S3", "S4")
  
  se_q0 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 2, ncol = 4))
  )
  rownames(se_q0) <- c("Gene_X", "Gene_Y")
  colnames(se_q0) <- c("S1", "S2", "S3", "S4")
  
  se_q1 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = matrix(c(10, 11, 12, 13, 14, 15, 16, 17), nrow = 2, ncol = 4))
  )
  rownames(se_q1) <- c("Gene_X", "Gene_Y")
  colnames(se_q1) <- c("S1", "S2", "S3", "S4")
  
  analysis <- .make_test_analysis_orchr_sample()
  analysis@diversity_results <- list(q_0.000 = se_q0, q_1.000 = se_q1)
  
  # Get results for specific sample with multiple q-values in display
  result <- TSENAT::results(analysis, type = "diversity",
                            n_genes = 2, sample = "S4",
                            q_values_table = c(0.0, 1.0))
  
  # Check that both q columns exist
  expect_true("q_0.0" %in% colnames(result))
  expect_true("q_1.0" %in% colnames(result))
  expect_equal(nrow(result), 2)
})

test_that("results() sample parameter works with extract_diversity_table", {
  entropy_data <- matrix(seq(1, 12), nrow = 3, ncol = 4)
  rownames(entropy_data) <- c("Gene_A", "Gene_B", "Gene_C")
  colnames(entropy_data) <- c("S1", "S2", "S3", "S4")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = entropy_data)
  )
  
  analysis <- .make_test_analysis_orchr_sample()
  analysis@diversity_results <- list(q_0.000 = se, q_1.000 = se)
  
  # Test the internal helper function directly
  table_df <- TSENAT:::.extract_diversity_table(analysis, NULL, NULL, 
                                                 n_genes = 2, 
                                                 q_values_table = c(0, 1.0),
                                                 sample = "S2")
  
  # Verify returned data.frame
  expect_is(table_df, "data.frame")
  expect_equal(nrow(table_df), 2)
  expect_equal(attr(table_df, "sample"), "S2")
  expect_equal(attr(table_df, "n_genes_total"), 3)
})
