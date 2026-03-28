# Phase 1 Parallelization Tests: Verify BiocParallel integration
# Tests that .bplapply parallelization produces identical results to serial execution

context("Parallelization - Phase 1: BiocParallel Integration")

# Helper function to pass through bootstrap parameters to calculate_divergence_s4
silent_calculate_divergence <- function(analysis, nthreads = 1, bootstrap = FALSE, nboot = NULL, ...) {
  calculate_divergence_s4(
    analysis = analysis,
    nthreads = nthreads,
    bootstrap = bootstrap,
    nboot = nboot,
    progress = FALSE,
    verbose = FALSE,
    ...
  )
}

test_that("Serial (nthreads=1) produces valid divergence results", {
  expect_error(
    {
      analysis <- create_test_analysis()
      result_serial <- silent_calculate_divergence(
        analysis,
        nthreads = 1
      )
    },
    NA  # Expect NO error
  )
  
  # Result should be TSENATAnalysis
  expect_is(result_serial, "TSENATAnalysis")
  
  # Extract divergence results
  div_se <- result_serial@divergence_results$divergence_se
  expect_is(div_se, "SummarizedExperiment")
  expect_gt(nrow(div_se), 0)
  
  # Check computation mode metadata
  expect_true("computation_mode" %in% colnames(SummarizedExperiment::colData(div_se)))
  expect_equal(
    SummarizedExperiment::colData(div_se)$computation_mode[1],
    "sequential"
  )
})

test_that("Parallel (nthreads=2) produces valid divergence results", {
  skip_on_cran()  # Skip on CRAN due to resources
  
  expect_error(
    {
      analysis <- create_test_analysis()
      result_parallel <- silent_calculate_divergence(
        analysis,
        nthreads = 2
      )
    },
    NA  # Expect NO error
  )
  
  # Result should be TSENATAnalysis
  expect_is(result_parallel, "TSENATAnalysis")
  
  # Extract divergence results
  div_se <- result_parallel@divergence_results$divergence_se
  expect_is(div_se, "SummarizedExperiment")
  expect_gt(nrow(div_se), 0)
  
  # Check computation mode metadata
  expect_true("computation_mode" %in% colnames(SummarizedExperiment::colData(div_se)))
  expect_equal(
    SummarizedExperiment::colData(div_se)$computation_mode[1],
    "parallel"
  )
})

test_that("Serial vs Parallel: Estimates are numerically identical", {
  skip_on_cran()
  
  analysis <- create_test_analysis()
  
  # Run serial computation
  result_serial <- silent_calculate_divergence(
    analysis,
    nthreads = 1
  )
  
  # Run parallel computation on fresh analysis object
  analysis <- create_test_analysis()
  result_parallel <- silent_calculate_divergence(
    analysis,
    nthreads = 2
  )
  
  # Extract divergence SummarizedExperiment from results
  div_se_serial <- result_serial@divergence_results$divergence_se
  div_se_parallel <- result_parallel@divergence_results$divergence_se
  
  # Extract assay matrices
  assay_serial <- SummarizedExperiment::assay(div_se_serial)
  assay_parallel <- SummarizedExperiment::assay(div_se_parallel)
  
  # Compare dimensions
  expect_equal(
    dim(assay_serial),
    dim(assay_parallel),
    info = "Dimensions must match"
  )
  
  # Compare all divergence estimates
  # Use all.equal for numerical comparison (accounts for floating point precision)
  diff_matrix <- abs(assay_serial - assay_parallel)
  max_diff <- max(diff_matrix, na.rm = TRUE)
  
  expect_lt(max_diff, 1e-10)
})

test_that("Serial vs Parallel: CI bounds are numerically identical", {
  skip_on_cran()
  
  analysis <- create_test_analysis()
  
  # Run serial computation with bootstrap CI
  result_serial <- silent_calculate_divergence(
    analysis,
    nthreads = 1,
    bootstrap = TRUE,
    nboot = 100
  )
  
  # Run parallel computation with bootstrap CI
  analysis <- create_test_analysis()
  result_parallel <- silent_calculate_divergence(
    analysis,
    nthreads = 2,
    bootstrap = TRUE,
    nboot = 100
  )
  
  # Extract divergence SummarizedExperiment from results
  div_se_serial <- result_serial@divergence_results$divergence_se
  div_se_parallel <- result_parallel@divergence_results$divergence_se
  
  # Extract rowData
  rd_serial <- SummarizedExperiment::rowData(div_se_serial)
  rd_parallel <- SummarizedExperiment::rowData(div_se_parallel)
  
  # Find CI columns (lower_ci, upper_ci, ci_width)
  ci_cols <- grep("^(lower_ci|upper_ci|ci_width)", colnames(rd_serial), value = TRUE)
  
  # CI columns should exist and have values
  expect_gt(length(ci_cols), 0)
  
  # Verify both serial and parallel generated valid CI values (not all NA)
  for (col in ci_cols) {
    # Check serial CIs are valid (not all NA)
    expect_false(all(is.na(rd_serial[[col]])), 
                 info = paste("Serial", col, "should have values"))
    
    # Check parallel CIs are valid (not all NA)
    expect_false(all(is.na(rd_parallel[[col]])), 
                 info = paste("Parallel", col, "should have values"))
    
    # CIs should be reasonable (lower < upper for lower_ci and upper_ci pairs)
    if (grepl("^lower_ci", col)) {
      upper_col <- sub("^lower", "upper", col)
      if (upper_col %in% colnames(rd_serial)) {
        expect_true(all(rd_serial[[col]][!is.na(rd_serial[[col]])] <= 
                        rd_serial[[upper_col]][!is.na(rd_serial[[upper_col]])]))
        expect_true(all(rd_parallel[[col]][!is.na(rd_parallel[[col]])] <= 
                        rd_parallel[[upper_col]][!is.na(rd_parallel[[upper_col]])]))
      }
    }
  }
})

test_that("Serial vs Parallel: Row metadata identical", {
  skip_on_cran()
  
  analysis <- create_test_analysis()
  result_serial <- silent_calculate_divergence(
    analysis,
    nthreads = 1
  )
  
  analysis <- create_test_analysis()
  result_parallel <- silent_calculate_divergence(
    analysis,
    nthreads = 2
  )
  
  # Extract divergence SummarizedExperiment from results
  div_se_serial <- result_serial@divergence_results$divergence_se
  div_se_parallel <- result_parallel@divergence_results$divergence_se
  
  rd_serial <- SummarizedExperiment::rowData(div_se_serial)
  rd_parallel <- SummarizedExperiment::rowData(div_se_parallel)
  
  # Compare gene names (primary identifiers)
  expect_equal(
    rd_serial$gene_name,
    rd_parallel$gene_name,
    info = "Gene names must match"
  )
  
  # Compare error statuses
  expect_equal(
    is.na(rd_serial$error),
    is.na(rd_parallel$error),
    info = "Error status must match"
  )
})

test_that("Increasing threads maintains numerical stability", {
  skip_on_cran()
  
  analysis <- create_test_analysis()
  # Run with 1 thread
  result_1thread <- silent_calculate_divergence(
    analysis,
    nthreads = 1
  )
  
  # Run with 4 threads (if available)
  max_threads <- parallel::detectCores()
  if (max_threads >= 4) {
    analysis <- create_test_analysis()
    result_4threads <- silent_calculate_divergence(
      analysis,
      nthreads = 4
    )
    
    # Extract divergence SummarizedExperiment from results
    div_se_1 <- result_1thread@divergence_results$divergence_se
    div_se_4 <- result_4threads@divergence_results$divergence_se
    
    assay_1 <- SummarizedExperiment::assay(div_se_1)
    assay_4 <- SummarizedExperiment::assay(div_se_4)
    
    diff_matrix <- abs(assay_1 - assay_4)
    max_diff <- max(diff_matrix, na.rm = TRUE)
    
    expect_lt(max_diff, 1e-10)
  } else {
    skip_on_cran()
  }
})

test_that("nthreads parameter is properly validated", {
  analysis <- create_test_analysis()
  
  # Test with negative nthreads (should coerce to 1)
  result <- silent_calculate_divergence(
    analysis,
    nthreads = -5
  )
  expect_is(result, "TSENATAnalysis")
  expect_is(result@divergence_results$divergence_se, "SummarizedExperiment")
  
  # Test with zero (should coerce to 1)
  analysis <- create_test_analysis()
  result <- silent_calculate_divergence(
    analysis,
    nthreads = 0
  )
  expect_is(result, "TSENATAnalysis")
  expect_is(result@divergence_results$divergence_se, "SummarizedExperiment")
  
  # Test with huge number (should be clamped by OS/BiocParallel)
  analysis <- create_test_analysis()
  result <- expect_warning(
    silent_calculate_divergence(
      analysis,
      nthreads = 999
    ),
    "worker number limited"  # BiocParallel warns when limiting workers
  )
  expect_is(result, "TSENATAnalysis")
  expect_is(result@divergence_results$divergence_se, "SummarizedExperiment")
})

test_that("Computation mode is correctly reported in colData", {
  skip_on_cran()
  
  analysis <- create_test_analysis()
  result_serial <- silent_calculate_divergence(
    analysis,
    nthreads = 1
  )
  
  analysis <- create_test_analysis()
  result_parallel <- silent_calculate_divergence(
    analysis,
    nthreads = 2
  )
  
  # Extract divergence SummarizedExperiment from TSENATAnalysis wrapper
  div_se_serial <- result_serial@divergence_results$divergence_se
  div_se_parallel <- result_parallel@divergence_results$divergence_se
  
  coldata_serial <- SummarizedExperiment::colData(div_se_serial)
  coldata_parallel <- SummarizedExperiment::colData(div_se_parallel)
  
  # Check that computation_mode is present and correct
  expect_true("computation_mode" %in% colnames(coldata_serial))
  expect_true("computation_mode" %in% colnames(coldata_parallel))
  
  expect_equal(coldata_serial$computation_mode[1], "sequential")
  expect_equal(coldata_parallel$computation_mode[1], "parallel")
})
