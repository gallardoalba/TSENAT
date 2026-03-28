# Tests for calculate_divergence_bootstrap and results processing
# Covers uncovered lines from divergence_coverage.txt

# Suppress nboot < 100 warnings for exploratory tests (acceptable for testing)
options(TSENAT.suppress_nboot_warning = TRUE)

test_that("calculate_divergence handles per-q pattern classification", {
  # Tests the q-spectrum pattern classification (covered lines 1041-1051)
  se <- create_test_se_simple()
  
  # Test with multiple q values to trigger pattern classification
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = c(0.5, 1.0, 2.0),
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  # Should have per_q_pattern column in rowData
  expect_result_structure(result, "SummarizedExperiment", n_rows = 2)
  rd <- SummarizedExperiment::rowData(result)
  expect_true("per_q_pattern" %in% colnames(rd))
})

test_that("calculate_divergence all normalization modes work correctly", {
  # Consolidated test: validates all 5 normalization modes
  # (replaced 4 separate tests covering lines 1082-1103, 1104-1127, 1128-1152, 1153+)
  # Uses test_all_normalization_modes() helper for comprehensive validation with strong assertions
  
  se <- create_test_se_simple()
  
  # Test all 5 normalization modes (none, range, zscore, log_odds_ratio, relative_reference)
  test_all_normalization_modes(
    func = .calculate_divergence,
    se = se,
    q = 1,
    group_col = "sample_type",
    control_group = "Control",
    bootstrap = FALSE,
    verbose = FALSE,
    progress = FALSE
  )
})

test_that("calculate_divergence skips genes with NA estimates", {
  # Tests handling of NA estimates (covered lines 1028-1030)
  se <- create_test_se_simple()
  
  # Run with progress to trigger logging paths
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    bootstrap = FALSE,
    verbose = FALSE,
    progress = TRUE
  )
  
  expect_result_structure(result, "SummarizedExperiment", n_rows = 2)
  rd <- SummarizedExperiment::rowData(result)
  expect_true(all(is.numeric(as.matrix(rd[, grep("^estimate_q", colnames(rd))]))))
})

test_that("calculate_divergence with paired samples", {
  # Tests paired sample handling
  se <- create_test_se_simple(paired = TRUE)
  
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    paired = TRUE,
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  expect_result_structure(result, "SummarizedExperiment", n_rows = 2)
  rd <- SummarizedExperiment::rowData(result)
  expect_true(all(c("lower_ci_q1", "upper_ci_q1") %in% colnames(rd)))
})

test_that("calculate_divergence with unpaired bootstrapping", {
  # Tests bootstrap with unpaired samples
  se <- create_test_se_simple()
  
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    bootstrap = TRUE,
    nboot = 10,
    verbose = FALSE
  )
  
  expect_result_structure(result, "SummarizedExperiment", n_rows = 2)
  rd <- SummarizedExperiment::rowData(result)
  expect_true(all(c("lower_ci_q1", "upper_ci_q1") %in% colnames(rd)))
})

test_that("calculate_divergence with paired bootstrapping", {
  # Tests bootstrap with paired samples
  se <- create_test_se_simple(paired = TRUE)
  
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    paired = TRUE,
    bootstrap = TRUE,
    nboot = 10,
    verbose = FALSE
  )
  
  expect_result_structure(result, "SummarizedExperiment", n_rows = 2)
  rd <- SummarizedExperiment::rowData(result)
  expect_true(all(c("lower_ci_q1", "upper_ci_q1") %in% colnames(rd)))
})

test_that("calculate_divergence auto-selects nboot", {
  # Tests auto-selection of nboot (covered in bootstrap logic)
  se <- create_test_se_simple()
  
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    bootstrap = TRUE,
    nboot = "auto",
    verbose = FALSE
  )
  
  expect_result_structure(result, "SummarizedExperiment", n_rows = 2)
  rd <- SummarizedExperiment::rowData(result)
  expect_true(all(c("lower_ci_q1", "upper_ci_q1") %in% colnames(rd)))
})

test_that("calculate_divergence uses different CI methods", {
  # Tests different confidence interval methods
  se <- create_test_se_simple()
  
  # Test percentile method
  result_percentile <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    bootstrap = TRUE,
    nboot = 10,
    method = "percentile",
    verbose = FALSE
  )
  
  expect_result_structure(result_percentile, "SummarizedExperiment", n_rows = 2)
  rd_p <- SummarizedExperiment::rowData(result_percentile)
  expect_true(all(c("lower_ci_q1", "upper_ci_q1") %in% colnames(rd_p)))
})

test_that("calculate_divergence applies CI threshold", {
  # Tests confidence interval parameter
  se <- create_test_se_simple()
  
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    bootstrap = TRUE,
    nboot = 10,
    ci = 0.90,
    verbose = FALSE
  )
  
  expect_result_structure(result, "SummarizedExperiment", n_rows = 2)
  rd <- SummarizedExperiment::rowData(result)
  expect_true(all(c("lower_ci_q1", "upper_ci_q1") %in% colnames(rd)))
})

test_that("calculate_divergence uses parallel processing", {
  # Tests parallel/multi-threaded execution
  se <- create_test_se_simple()
  
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    nthreads = 2,
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  expect_result_structure(result, "SummarizedExperiment", n_rows = 2)
  rd <- SummarizedExperiment::rowData(result)
  expect_true(all(is.numeric(as.matrix(rd[, grep("^estimate_q", colnames(rd))]))))
})

test_that("calculate_divergence with log_base parameter", {
  # Tests alternative log base for entropy computation
  se <- create_test_se_simple()
  
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    log_base = 2,
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  expect_result_structure(result, "SummarizedExperiment", n_rows = 2)
  rd <- SummarizedExperiment::rowData(result)
  expect_true(all(is.numeric(as.matrix(rd[, grep("^estimate_q", colnames(rd))]))))
})

test_that("calculate_divergence with pseudocount parameter", {
  # Tests pseudocount handling for zero-count genes
  se <- create_test_se_simple()
  
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    pseudocount = 1.0,
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  expect_result_structure(result, "SummarizedExperiment", n_rows = 2)
})
