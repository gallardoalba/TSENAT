# Tests for calculate_lm_interaction_s4 S4 wrapper function

# Setup: Local helper to create test analysis with optional diversity results
.create_test_analysis <- function(
    precompute_diversity = TRUE, 
    q_values = c(0.5, 0.75, 1.0, 1.5, 2.0),
    include_divergence = TRUE) {
  set.seed(42)
  
  # Always create base analysis with diversity first
  analysis <- create_test_analysis(
    q_values = q_values,
    include_divergence = include_divergence,
    include_lm_results = FALSE
  )
  
  # If precompute_diversity = FALSE, clear the diversity results
  if (!precompute_diversity) {
    analysis@diversity_results <- list()
  }
  
  return(analysis)
}

test_that("S4 Wrappers: calculate_lm_interaction_s4 accepts new arguments", {
  # Create analysis with diversity results (default behavior)
  analysis <- .create_test_analysis()
  
  # Ensure diversity results exist
  if (length(analysis@diversity_results) == 0) {
    skip_on_cran()
  }
  
  # Test each argument individually
  args_to_test <- list(
    list(condition_col = "condition"),
    list(method = "lmm"),
    list(paired = FALSE),
    list(pcorr = "BH"),
    list(return_model_data = TRUE),
    list(nthreads = 1)
  )
  
  for (args in args_to_test) {
    arg_string <- paste(names(args), collapse = ", ")
    # Test that arguments are accepted
    result <- tryCatch({
      do.call(calculate_lm_interaction_s4, 
              c(list(analysis = analysis, verbose = FALSE), args))
    }, error = function(e) {
      list(error = paste("Error:", e$message))
    })
    
    # Should accept arguments without syntax errors
    expect_true(!is.list(result) || !("error" %in% names(result)),
                info = paste("Failed for argument:", arg_string,
                            "Error:", if(is.list(result) && "error" %in% names(result)) result$error else "None"))
  }
})

# Tests for .sync_coldata_from_diversity helper
test_that("Helper: .sync_coldata_from_diversity syncs colData from diversity results", {
  analysis <- .create_test_analysis()
  
  if (length(analysis@diversity_results) == 0) {
    skip("No diversity results available")
  }
  
  # Get original colData
  orig_coldata <- SummarizedExperiment::colData(analysis@se)
  
  # Sync colData
  result <- .sync_coldata_from_diversity(analysis, verbose = FALSE)
  
  # Should return a TSENATAnalysis object
  expect_s4_class(result, "TSENATAnalysis")
  
  # colData should be updated
  new_coldata <- SummarizedExperiment::colData(result@se)
  expect_equal(nrow(new_coldata), nrow(orig_coldata))
})

test_that("Helper: .sync_coldata_from_diversity handles empty diversity results", {
  # Create minimal analysis without diversity (set q_values to empty vector)
  analysis <- .create_test_analysis(precompute_diversity = FALSE)
  
  # Empty diversity results should return unchanged analysis
  result <- .sync_coldata_from_diversity(analysis, verbose = FALSE)
  expect_s4_class(result, "TSENATAnalysis")
  expect_equal(length(result@diversity_results), 0)
})

# Tests for .extract_lm_params helper
test_that("Helper: .extract_lm_params extracts and resolves parameters", {
  analysis <- .create_test_analysis()
  
  # Extract parameters with defaults
  params <- .extract_lm_params(
    analysis,
    condition_col = "condition",
    method = "gam",
    pcorr = "BY",
    verbose = FALSE
  )
  
  expect_is(params, "list")
  expect_equal(params$condition_col, "condition")
  expect_equal(params$method, "gam")
  expect_equal(params$pcorr, "BY")
  expect_equal(params$paired, FALSE)
})

test_that("Helper: .extract_lm_params resolves from config", {
  analysis <- .create_test_analysis()
  analysis@config$method <- "lmm"
  analysis@config$pcorr <- "Hochberg"
  
  # Extract with NULL parameters (should use config)
  params <- .extract_lm_params(
    analysis,
    condition_col = NULL,
    method = NULL,
    pcorr = NULL,
    verbose = FALSE
  )
  
  expect_is(params, "list")
  expect_equal(params$method, "lmm")
  expect_equal(params$pcorr, "Hochberg")
})

# Tests for .combine_diversity_results_for_lm helper
test_that("Helper: .combine_diversity_results_for_lm combines multi-q results", {
  analysis <- .create_test_analysis()
  
  if (length(analysis@diversity_results) < 2) {
    skip("Need multiple q-values for this test")
  }
  
  # Combine diversity results
  combined_se <- .combine_diversity_results_for_lm(analysis@diversity_results)
  
  expect_s4_class(combined_se, "SummarizedExperiment")
  expect_gt(ncol(combined_se), 0)
  
  # Should have combined assays and colData
  expect_equal(ncol(combined_se), nrow(SummarizedExperiment::colData(combined_se)))
  
  # Assay should be named 'diversity'
  expect_true("diversity" %in% SummarizedExperiment::assayNames(combined_se))
})

test_that("Helper: .combine_diversity_results_for_lm adds q-value suffixes", {
  analysis <- .create_test_analysis()
  
  if (length(analysis@diversity_results) < 2) {
    skip("Need multiple q-values for this test")
  }
  
  combined_se <- .combine_diversity_results_for_lm(analysis@diversity_results)
  
  # Check that column names have q-value suffixes
  col_names <- colnames(combined_se)
  has_q_suffix <- any(grepl("_q=", col_names))
  expect_true(has_q_suffix, info = "Column names should have _q= suffix")
})

# Tests for .build_lm_args helper
test_that("Helper: .build_lm_args builds argument list correctly", {
  analysis <- .create_test_analysis()
  
  if (length(analysis@diversity_results) == 0) {
    skip("No diversity results available")
  }
  
  combined_se <- .combine_diversity_results_for_lm(analysis@diversity_results)
  
  params <- list(
    condition_col = "condition",
    method = "gam",
    paired = FALSE,
    pcorr = "BH"
  )
  
  args <- .build_lm_args(combined_se, params, return_model_data = TRUE, verbose = FALSE)
  
  expect_is(args, "list")
  expect_true("se" %in% names(args))
  expect_true("condition_col" %in% names(args))
  expect_true("method" %in% names(args))
  expect_true("return_model_data" %in% names(args))
  expect_equal(args$return_model_data, TRUE)
})

test_that("Helper: .build_lm_args respects NULL parameters", {
  analysis <- .create_test_analysis()
  
  if (length(analysis@diversity_results) == 0) {
    skip("No diversity results available")
  }
  
  combined_se <- .combine_diversity_results_for_lm(analysis@diversity_results)
  
  params <- list(
    condition_col = NULL,
    method = NULL,
    paired = FALSE,
    pcorr = NULL
  )
  
  args <- .build_lm_args(combined_se, params, return_model_data = FALSE, verbose = FALSE)
  
  # NULL parameters should not be in args (or should be handled gracefully)
  expect_is(args, "list")
  expect_true("se" %in% names(args))
})

# Tests for .validate_and_extract_lm_result helper
test_that("Helper: .validate_and_extract_lm_result validates data frame", {
  # Create a valid result
  result_df <- data.frame(
    gene = c("gene1", "gene2"),
    adj_p_interaction = c(0.01, 0.05),
    row.names = NULL
  )
  
  extracted <- .validate_and_extract_lm_result(result_df)
  
  expect_is(extracted, "list")
  expect_true("results" %in% names(extracted))
  expect_true("model_data" %in% names(extracted))
  expect_equal(nrow(extracted$results), 2)
})

test_that("Helper: .validate_and_extract_lm_result extracts from list result", {
  result_list <- list(
    results = data.frame(
      gene = c("gene1", "gene2"),
      adj_p_interaction = c(0.01, 0.05)
    ),
    model_data = list(some_model = "data")
  )
  
  extracted <- .validate_and_extract_lm_result(result_list)
  
  expect_is(extracted, "list")
  expect_equal(nrow(extracted$results), 2)
  expect_is(extracted$model_data, "list")
})

test_that("Helper: .validate_and_extract_lm_result handles empty results", {
  empty_df <- data.frame(
    gene = character(0),
    adj_p_interaction = numeric(0)
  )
  
  # Expect a warning about empty results (expected behavior)
  extracted <- expect_warning(
    .validate_and_extract_lm_result(empty_df),
    "Result is empty"
  )
  
  expect_is(extracted, "list")
  expect_equal(nrow(extracted$results), 0)
})

test_that("Helper: .validate_and_extract_lm_result detects missing columns", {
  incomplete_df <- data.frame(
    gene = c("gene1", "gene2"),
    p_value = c(0.01, 0.05)
  )
  
  expect_error(
    .validate_and_extract_lm_result(incomplete_df),
    "Missing columns"
  )
})

# Tests for .store_lm_results_in_analysis helper
test_that("Helper: .store_lm_results_in_analysis stores results in analysis", {
  analysis <- .create_test_analysis()
  
  lm_results_df <- data.frame(
    gene = c("gene1", "gene2"),
    adj_p_interaction = c(0.01, 0.05)
  )
  
  result <- .store_lm_results_in_analysis(analysis, lm_results_df)
  
  expect_s4_class(result, "TSENATAnalysis")
  expect_true("lm_interaction" %in% names(result@lm_results))
  expect_equal(nrow(result@lm_results$lm_interaction), 2)
})

test_that("Helper: .store_lm_results_in_analysis stores model_data if provided", {
  analysis <- .create_test_analysis()
  
  lm_results_df <- data.frame(
    gene = c("gene1", "gene2"),
    adj_p_interaction = c(0.01, 0.05)
  )
  
  model_data <- list(models = "some_models", metadata = "test")
  
  result <- .store_lm_results_in_analysis(analysis, lm_results_df, model_data)
  
  expect_s4_class(result, "TSENATAnalysis")
  expect_true("lm_interaction_model_data" %in% names(result@lm_results))
  expect_is(result@lm_results$lm_interaction_model_data, "list")
})

test_that("Helper: .store_lm_results_in_analysis tracks function call", {
  analysis <- .create_test_analysis()
  
  lm_results_df <- data.frame(
    gene = c("gene1"),
    adj_p_interaction = c(0.01)
  )
  
  result <- .store_lm_results_in_analysis(analysis, lm_results_df)
  
  # Check that function call was tracked
  expect_true("calculate_lm_interaction" %in% result@metadata$function_calls)
})
