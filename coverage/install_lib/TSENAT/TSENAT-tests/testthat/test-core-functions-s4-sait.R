# Tests for calculate_sait S4 wrapper function

# Skip entire test file on Bioconductor due to long runtime (10.37s)
skip_on_bioc()

# Setup: Local helper to create test analysis with optional diversity and SAIT results
.create_test_analysis <- function(
    precompute_diversity = TRUE, 
    q_values = c(0.5, 0.75, 1.0, 1.5, 2.0),
    include_divergence = TRUE,
    include_sait_results = FALSE,
    seed = 42,
    verbose = FALSE) {
  set.seed(seed)
  
  # Always create base analysis with diversity first
  analysis <- create_test_analysis(
    q_values = q_values,
    include_divergence = include_divergence,
    include_sait_results = FALSE,
    seed = seed,
    verbose = verbose
  )
  
  # If precompute_diversity = FALSE, clear the diversity results
  if (!precompute_diversity) {
    analysis@diversity_results <- list()
  }

  if (include_sait_results) {
    analysis@sait_results <- list(
      sait_interaction = data.frame(
        gene = character(0),
        adj_p_interaction = numeric(0),
        stringsAsFactors = FALSE
      )
    )
  }
  
  return(analysis)
}

test_that("S4 Wrappers: calculate_sait accepts new arguments", {
  # Create analysis with diversity results (default behavior)
  analysis <- .create_test_analysis()
  
  # Ensure diversity results exist
  if (length(analysis@diversity_results) == 0) {
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
      do.call(calculate_sait, 
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

# Tests for .extract_sait_params helper
test_that("Helper: .extract_sait_params extracts and resolves parameters", {
  analysis <- .create_test_analysis()
  
  # Extract parameters with defaults
  params <- .extract_sait_params(
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

test_that("S4 wrapper: calculate_sait fails fast on invalid method arguments", {
  analysis <- .create_test_analysis()

  expect_error(
    calculate_sait(analysis, method = "definitely_not_a_real_method", verbose = FALSE),
    "should be one of|not found|Available columns",
    fixed = FALSE
  )
})

test_that("Helper: .extract_sait_params resolves from config", {
  analysis <- .create_test_analysis()
  analysis@config$method <- "lmm"
  analysis@config$pcorr <- "Hochberg"
  
  # Extract with NULL parameters (should use config)
  params <- .extract_sait_params(
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

test_that("calculate_sait fails when condition_col cannot be auto-detected", {
  analysis <- .create_test_analysis(include_sait_results = TRUE)
  analysis@config$condition_col <- NULL

  cd <- SummarizedExperiment::colData(analysis@se)
  empty_cd <- S4Vectors::DataFrame(
    sample_id = cd$sample_id,
    row.names = rownames(cd)
  )
  SummarizedExperiment::colData(analysis@se) <- empty_cd

  # Ensure the diversity results also do not supply a condition column after sync
  for (i in seq_along(analysis@diversity_results)) {
    se_q <- analysis@diversity_results[[i]]
    if (is(se_q, "SummarizedExperiment")) {
      SummarizedExperiment::colData(se_q) <- empty_cd
      analysis@diversity_results[[i]] <- se_q
    }
  }

  expect_error(
    calculate_sait(analysis, verbose = FALSE),
    "condition_col is required for SAIT interaction analysis"
  )
})

# Tests for .combine_diversity_results_for_sait helper
test_that("Helper: .combine_diversity_results_for_sait combines multi-q results", {
  analysis <- .create_test_analysis()
  
  if (length(analysis@diversity_results) < 2) {
    skip("Need multiple q-values for this test")
  }
  
  # Combine diversity results
  combined_se <- .combine_diversity_results_for_sait(analysis@diversity_results)
  
  expect_s4_class(combined_se, "SummarizedExperiment")
  expect_gt(ncol(combined_se), 0)
  
  # Should have combined assays and colData
  expect_equal(ncol(combined_se), nrow(SummarizedExperiment::colData(combined_se)))
  
  # Assay should be named 'diversity'
  expect_true("diversity" %in% SummarizedExperiment::assayNames(combined_se))
})

test_that("Helper: .combine_diversity_results_for_sait adds q-value suffixes", {
  analysis <- .create_test_analysis()
  
  if (length(analysis@diversity_results) < 2) {
    skip("Need multiple q-values for this test")
  }
  
  combined_se <- .combine_diversity_results_for_sait(analysis@diversity_results)
  
  # Check that column names have q-value suffixes
  col_names <- colnames(combined_se)
  has_q_suffix <- any(grepl("_q=", col_names))
  expect_true(has_q_suffix, info = "Column names should have _q= suffix")
})

# Tests for .build_sait_args helper
test_that("Helper: .build_sait_args builds argument list correctly", {
  analysis <- .create_test_analysis()
  
  if (length(analysis@diversity_results) == 0) {
    skip("No diversity results available")
  }
  
  combined_se <- .combine_diversity_results_for_sait(analysis@diversity_results)
  
  params <- list(
    condition_col = "condition",
    method = "gam",
    paired = FALSE,
    pcorr = "BH"
  )
  
  args <- .build_sait_args(combined_se, params, return_model_data = TRUE, verbose = FALSE)
  
  expect_is(args, "list")
  expect_true("se" %in% names(args))
  expect_true("condition_col" %in% names(args))
  expect_true("method" %in% names(args))
  expect_true("return_model_data" %in% names(args))
  expect_equal(args$return_model_data, TRUE)
})

test_that("Helper: .build_sait_args respects NULL parameters", {
  analysis <- .create_test_analysis()
  
  if (length(analysis@diversity_results) == 0) {
    skip("No diversity results available")
  }
  
  combined_se <- .combine_diversity_results_for_sait(analysis@diversity_results)
  
  params <- list(
    condition_col = NULL,
    method = NULL,
    paired = FALSE,
    pcorr = NULL
  )
  
  args <- .build_sait_args(combined_se, params, return_model_data = FALSE, verbose = FALSE)
  
  # NULL parameters should not be in args (or should be handled gracefully)
  expect_is(args, "list")
  expect_true("se" %in% names(args))
})

# Tests for .validate_and_extract_sait_result helper
test_that("Helper: .validate_and_extract_sait_result validates data frame", {
  # Create a valid result
  result_df <- data.frame(
    gene = c("gene1", "gene2"),
    adj_p_interaction = c(0.01, 0.05),
    row.names = NULL
  )
  
  extracted <- .validate_and_extract_sait_result(result_df)
  
  expect_is(extracted, "list")
  expect_true("results" %in% names(extracted))
  expect_true("model_data" %in% names(extracted))
  expect_equal(nrow(extracted$results), 2)
})

test_that("Helper: .validate_and_extract_sait_result extracts from list result", {
  result_list <- list(
    results = data.frame(
      gene = c("gene1", "gene2"),
      adj_p_interaction = c(0.01, 0.05)
    ),
    model_data = list(some_model = "data")
  )
  
  extracted <- .validate_and_extract_sait_result(result_list)
  
  expect_is(extracted, "list")
  expect_equal(nrow(extracted$results), 2)
  expect_is(extracted$model_data, "list")
})

test_that("Helper: .validate_and_extract_sait_result handles empty results", {
  empty_df <- data.frame(
    gene = character(0),
    adj_p_interaction = numeric(0)
  )
  
  # Expect a warning about empty results (expected behavior)
  extracted <- expect_warning(
    .validate_and_extract_sait_result(empty_df),
    "Result is empty"
  )
  
  expect_is(extracted, "list")
  expect_equal(nrow(extracted$results), 0)
})

test_that("Helper: .validate_and_extract_sait_result detects missing columns", {
  incomplete_df <- data.frame(
    gene = c("gene1", "gene2"),
    p_value = c(0.01, 0.05)
  )
  
  expect_error(
    .validate_and_extract_sait_result(incomplete_df),
    "Missing columns"
  )
})

# Tests for .store_sait_results_in_analysis helper
test_that("Helper: .store_sait_results_in_analysis stores results in analysis", {
  analysis <- .create_test_analysis()
  
  sait_results_df <- data.frame(
    gene = c("gene1", "gene2"),
    adj_p_interaction = c(0.01, 0.05)
  )
  
  result <- .store_sait_results_in_analysis(analysis, sait_results_df)
  
  expect_s4_class(result, "TSENATAnalysis")
  expect_true("sait_interaction" %in% names(result@sait_results))
  expect_equal(nrow(result@sait_results$sait_interaction), 2)
})

test_that("Helper: .store_sait_results_in_analysis stores model_data if provided", {
  analysis <- .create_test_analysis()
  
  sait_results_df <- data.frame(
    gene = c("gene1", "gene2"),
    adj_p_interaction = c(0.01, 0.05)
  )
  
  model_data <- list(models = "some_models", metadata = "test")
  
  result <- .store_sait_results_in_analysis(analysis, sait_results_df, model_data)
  
  expect_s4_class(result, "TSENATAnalysis")
  expect_true("sait_interaction_model_data" %in% names(result@sait_results))
  expect_is(result@sait_results$sait_interaction_model_data, "list")
})

test_that("Helper: .store_sait_results_in_analysis tracks function call", {
  analysis <- .create_test_analysis()
  
  sait_results_df <- data.frame(
    gene = c("gene1"),
    adj_p_interaction = c(0.01)
  )
  
  result <- .store_sait_results_in_analysis(analysis, sait_results_df)
  
  # Check that function call was tracked
  expect_true("calculate_sait" %in% result@metadata$function_calls)
})

# ============================================================================
# TESTS FOR REFACTORED HELPERS: ._sait_resolve_params, ._sait_execute_and_store
# ============================================================================

test_that("._sait_resolve_params resolves all parameters from config", {
    analysis <- .create_test_analysis()
    analysis@config$fdr_threshold <- 0.1
    analysis@config$verbose <- TRUE
    analysis@config$paired <- FALSE
    analysis@config$return_model_data <- FALSE
    
    SummarizedExperiment::colData(analysis@se)$condition <- c("A", "A", "B", "B")
    analysis@config$condition_col <- "condition"
    
    resolved <- TSENAT:::._sait_resolve_params(
        analysis, fdr_threshold = NULL, formula = NULL,
        condition_col = NULL, method = "lmm", paired = NULL,
        subject_col = NULL, nthreads = NULL, multicorr = NULL,
        corstr = NULL, pcorr = NULL, verbose = NULL,
        return_model_data = NULL, output_file = NULL
    )
    
    expect_equal(resolved$fdr_threshold, 0.1)
    expect_true(resolved$verbose)
    expect_false(resolved$paired)
    expect_false(resolved$return_model_data)
    expect_equal(resolved$params$method, "lmm")
})

test_that("._sait_execute_and_store returns analysis with results", {
    analysis <- .create_test_analysis()
    SummarizedExperiment::colData(analysis@se)$condition <- c("A", "A", "B", "B")
    
    # Add diversity results: one SE per q-value, each with 4 genes x 4 samples
    n_genes <- 4
    n_samples <- 4
    q_vals <- c(0.1, 0.5, 1.0, 1.5, 2.0)
    div_results <- list()
    for (q in q_vals) {
        div_mat <- matrix(rnorm(n_genes * n_samples, mean = 0.8, sd = 0.1), nrow = n_genes)
        colnames(div_mat) <- paste0("S", 1:n_samples)
        rownames(div_mat) <- paste0("Gene", 1:n_genes)
        q_key <- paste0("q_", sprintf("%.2f", q))
        div_results[[q_key]] <- SummarizedExperiment::SummarizedExperiment(
            assays = list(diversity = div_mat),
            colData = data.frame(condition = rep(c("A", "B"), each = 2),
                                 row.names = paste0("S", 1:n_samples))
        )
    }
    analysis@diversity_results <- div_results
    
    resolved <- list(
        fdr_threshold = NULL, formula = NULL, output_file = NULL,
        verbose = FALSE, paired = FALSE, return_model_data = TRUE,
        params = list(
            condition_col = "condition", method = "lmm", subject_col = NULL,
            nthreads = 1, multicorr = "hochberg", corstr = NULL, pcorr = "BH",
            paired = FALSE
        )
    )
    
    result <- TSENAT:::._sait_execute_and_store(analysis, resolved)
    
    expect_s4_class(result, "TSENATAnalysis")
    expect_true("sait_interaction" %in% names(result@sait_results))
})

test_that("._sait_execute_and_store handles empty results gracefully", {
    analysis <- .create_test_analysis()
    SummarizedExperiment::colData(analysis@se)$condition <- c("A", "A", "B", "B")
    
    # Add diversity results: one SE per q-value
    n_genes <- 4
    n_samples <- 4
    q_vals <- c(0.1, 0.5, 1.0, 1.5, 2.0)
    div_results <- list()
    for (q in q_vals) {
        div_mat <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)
        colnames(div_mat) <- paste0("S", 1:n_samples)
        rownames(div_mat) <- paste0("Gene", 1:n_genes)
        q_key <- paste0("q_", sprintf("%.2f", q))
        div_results[[q_key]] <- SummarizedExperiment::SummarizedExperiment(
            assays = list(diversity = div_mat),
            colData = data.frame(condition = rep(c("A", "B"), each = 2),
                                 row.names = paste0("S", 1:n_samples))
        )
    }
    analysis@diversity_results <- div_results
    
    resolved <- list(
        fdr_threshold = NULL, formula = NULL, output_file = NULL,
        verbose = FALSE, paired = FALSE, return_model_data = FALSE,
        params = list(
            condition_col = "condition", method = "lmm", subject_col = NULL,
            nthreads = 1, multicorr = "hochberg", corstr = NULL, pcorr = "BH",
            paired = FALSE
        )
    )
    
    result <- TSENAT:::._sait_execute_and_store(analysis, resolved)
    
    expect_s4_class(result, "TSENATAnalysis")
})
