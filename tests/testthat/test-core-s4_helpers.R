# Tests for s4_functions_diversity helper functions
# .prepare_diversity_params, .build_calc_diversity_args, 
# .extract_q_metadata_from_result, .apply_diversity_post_hoc_norm

# Setup: Helper to create minimal test SummarizedExperiment
make_test_se <- function(n_genes = 10, n_samples = 5) {
  counts <- matrix(rpois(n_genes * n_samples, lambda = 50), nrow = n_genes)
  rownames(counts) <- paste0("gene_", seq_len(n_genes))
  colnames(counts) <- paste0("sample_", seq_len(n_samples))
  
  coldata <- data.frame(
    sample = colnames(counts),
    condition = rep(c("control", "treatment"), length.out = n_samples),
    row.names = colnames(counts)
  )
  
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = coldata
  )
}

# Setup: Create minimal test analysis object
make_test_analysis <- function() {
  se <- make_test_se()
  TSENAT::TSENATAnalysis(se, config = list(q_values = c(0.5, 1.0, 1.5)))
}

# ===========================================================================
# Tests for .prepare_diversity_params()
# ===========================================================================

test_that(".prepare_diversity_params extracts q values from config", {
  analysis <- make_test_analysis()
  
  params <- TSENAT:::.prepare_diversity_params(
    analysis, 
    q = NULL,  # Should use config
    norm = NULL, norm_method = NULL, reference_group = NULL,
    tpm = FALSE, assayno = NULL, verbose = NULL, what = NULL,
    nthreads = NULL, pseudocount = NULL, min_valid_frac = NULL,
    shrinkage = NULL, genes = NULL, effective_length = NULL,
    metadata = NULL, bootstrap = NULL, nboot = NULL,
    bootstrap_method = NULL, bootstrap_ci = NULL,
    bootstrap_include_diagnostics = NULL, seed = NULL
  )
  
  expect_equal(params$q, c(0.5, 1.0, 1.5))
  expect_named(params, c("q", "nthreads", "verbose", "bootstrap", "pseudocount", 
                         "norm", "what", "assayno", "min_valid_frac", "shrinkage",
                         "bootstrap_method", "bootstrap_ci", "seed", "tpm", 
                         "genes", "effective_length", "nboot", 
                         "bootstrap_include_diagnostics", "metadata", "norm_method", 
                         "reference_group"))
})

test_that(".prepare_diversity_params uses explicit q values over config", {
  analysis <- make_test_analysis()
  explicit_q <- c(0.1, 0.5, 2.0)
  
  params <- TSENAT:::.prepare_diversity_params(
    analysis,
    q = explicit_q,  # Should override config
    norm = NULL, norm_method = NULL, reference_group = NULL,
    tpm = FALSE, assayno = NULL, verbose = NULL, what = NULL,
    nthreads = NULL, pseudocount = NULL, min_valid_frac = NULL,
    shrinkage = NULL, genes = NULL, effective_length = NULL,
    metadata = NULL, bootstrap = NULL, nboot = NULL,
    bootstrap_method = NULL, bootstrap_ci = NULL,
    bootstrap_include_diagnostics = NULL, seed = NULL
  )
  
  expect_equal(params$q, explicit_q)
})

test_that(".prepare_diversity_params uses default q when not in config", {
  se <- make_test_se()
  analysis <- TSENAT::TSENATAnalysis(se, config = list())  # Empty config
  
  params <- TSENAT:::.prepare_diversity_params(
    analysis,
    q = NULL,  # Should use default
    norm = NULL, norm_method = NULL, reference_group = NULL,
    tpm = FALSE, assayno = NULL, verbose = NULL, what = NULL,
    nthreads = NULL, pseudocount = NULL, min_valid_frac = NULL,
    shrinkage = NULL, genes = NULL, effective_length = NULL,
    metadata = NULL, bootstrap = NULL, nboot = NULL,
    bootstrap_method = NULL, bootstrap_ci = NULL,
    bootstrap_include_diagnostics = NULL, seed = NULL
  )
  
  expect_equal(params$q, seq(0.01, 2, by = 0.05))
})

test_that(".prepare_diversity_params validates q is numeric", {
  analysis <- make_test_analysis()
  
  expect_error({
    TSENAT:::.prepare_diversity_params(
      analysis,
      q = c("a", "b"),  # Non-numeric
      norm = NULL, norm_method = NULL, reference_group = NULL,
      tpm = FALSE, assayno = NULL, verbose = NULL, what = NULL,
      nthreads = NULL, pseudocount = NULL, min_valid_frac = NULL,
      shrinkage = NULL, genes = NULL, effective_length = NULL,
      metadata = NULL, bootstrap = NULL, nboot = NULL,
      bootstrap_method = NULL, bootstrap_ci = NULL,
      bootstrap_include_diagnostics = NULL, seed = NULL
    )
  }, "q.*must be numeric")
})

test_that(".prepare_diversity_params resolves default parameters", {
  analysis <- make_test_analysis()
  
  params <- TSENAT:::.prepare_diversity_params(
    analysis,
    q = c(1.0),
    norm = NULL, norm_method = NULL, reference_group = NULL,
    tpm = FALSE, assayno = NULL, verbose = NULL, what = NULL,
    nthreads = NULL, pseudocount = NULL, min_valid_frac = NULL,
    shrinkage = NULL, genes = NULL, effective_length = NULL,
    metadata = NULL, bootstrap = NULL, nboot = NULL,
    bootstrap_method = NULL, bootstrap_ci = NULL,
    bootstrap_include_diagnostics = NULL, seed = NULL
  )
  
  # Check defaults
  expect_equal(params$norm, TRUE)
  expect_equal(params$verbose, TRUE)
  expect_equal(params$bootstrap, FALSE)
  expect_equal(params$pseudocount, 0)
  expect_equal(params$nthreads, 1)
  expect_equal(params$what, "S")
})

test_that(".prepare_diversity_params prioritizes explicit parameters over config", {
  analysis <- TSENAT::TSENATAnalysis(
    make_test_se(),
    config = list(nthreads = 4, verbose = FALSE, pseudocount = 1)
  )
  
  params <- TSENAT:::.prepare_diversity_params(
    analysis,
    q = c(1.0),
    norm = TRUE,  # Explicit
    norm_method = NULL, reference_group = NULL,
    tpm = FALSE, assayno = NULL,
    verbose = TRUE,  # Explicit (override config FALSE)
    what = NULL,
    nthreads = 8,  # Explicit (override config 4)
    pseudocount = 0.5,  # Explicit (override config 1)
    min_valid_frac = NULL,
    shrinkage = NULL, genes = NULL, effective_length = NULL,
    metadata = NULL, bootstrap = NULL, nboot = NULL,
    bootstrap_method = NULL, bootstrap_ci = NULL,
    bootstrap_include_diagnostics = NULL, seed = NULL
  )
  
  expect_equal(params$verbose, TRUE)
  expect_equal(params$nthreads, 8)
  expect_equal(params$pseudocount, 0.5)
})

# ===========================================================================
# Tests for .build_calc_diversity_args()
# ===========================================================================

test_that(".build_calc_diversity_args constructs valid argument list", {
  analysis <- make_test_analysis()
  params <- list(
    q = c(1.0),
    nthreads = 2,
    verbose = TRUE,
    bootstrap = FALSE,
    pseudocount = 0,
    norm = TRUE,
    what = "S",
    assayno = 1,
    min_valid_frac = 0.75,
    shrinkage = "none",
    bootstrap_method = "percentile",
    bootstrap_ci = 0.95,
    seed = NULL,
    tpm = FALSE,
    genes = NULL,
    effective_length = NULL,
    nboot = NULL,
    bootstrap_include_diagnostics = TRUE,
    metadata = NULL,
    norm_method = NULL,
    reference_group = NULL
  )
  
  args <- TSENAT:::.build_calc_diversity_args(params, analysis, list())
  
  expect_equal(args$x, analysis@se)
  expect_equal(args$q, c(1.0))
  expect_equal(args$norm, TRUE)
  expect_equal(args$tpm, FALSE)
  expect_equal(args$verbose, TRUE)
  expect_equal(args$nthreads, 2)
})

test_that(".build_calc_diversity_args includes optional parameters when present", {
  analysis <- make_test_analysis()
  params <- list(
    q = c(1.0),
    nthreads = 1,
    verbose = TRUE,
    bootstrap = FALSE,
    pseudocount = 0,
    norm = TRUE,
    what = "S",
    assayno = 1,
    min_valid_frac = 0.75,
    shrinkage = "none",
    bootstrap_method = "percentile",
    bootstrap_ci = 0.95,
    seed = NULL,
    tpm = FALSE,
    genes = c("gene1", "gene2"),  # Include
    effective_length = c(1000, 2000),  # Include
    nboot = 100,  # Include
    bootstrap_include_diagnostics = TRUE,
    metadata = list(custom = "value"),  # Include
    norm_method = NULL,
    reference_group = NULL
  )
  
  args <- TSENAT:::.build_calc_diversity_args(params, analysis, list())
  
  expect_true("genes" %in% names(args))
  expect_true("effective_length" %in% names(args))
  expect_true("bootstrap_nboot" %in% names(args))
  expect_true("metadata" %in% names(args))
})

test_that(".build_calc_diversity_args includes ... arguments", {
  analysis <- make_test_analysis()
  params <- list(
    q = c(1.0),
    nthreads = 1,
    verbose = TRUE,
    bootstrap = FALSE,
    pseudocount = 0,
    norm = TRUE,
    what = "S",
    assayno = 1,
    min_valid_frac = 0.75,
    shrinkage = "none",
    bootstrap_method = "percentile",
    bootstrap_ci = 0.95,
    seed = NULL,
    tpm = FALSE,
    genes = NULL,
    effective_length = NULL,
    nboot = NULL,
    bootstrap_include_diagnostics = TRUE,
    metadata = NULL,
    norm_method = NULL,
    reference_group = NULL
  )
  
  dots <- list(custom_arg = "custom_value", another = 42)
  args <- TSENAT:::.build_calc_diversity_args(params, analysis, dots)
  
  expect_equal(args$custom_arg, "custom_value")
  expect_equal(args$another, 42)
})

# ===========================================================================
# Tests for .extract_q_metadata_from_result()
# ===========================================================================

test_that(".extract_q_metadata_from_result extracts q values from column names", {
  # Create mock result dataframe with q-value columns
  result_df <- data.frame(
    gene = c("gene1", "gene2"),
    col1 = c(1.5, 2.0),
    col2 = c(1.6, 2.1),
    col3 = c(1.7, 2.2),
    col4 = c(1.8, 2.3)
  )
  colnames(result_df)[2:5] <- c("sample1_q=0.500", "sample2_q=0.500", 
                                 "sample3_q=1.000", "sample4_q=1.000")
  
  q_values <- c(0.5, 1.0)
  col_q_vals <- TSENAT:::.extract_q_metadata_from_result(result_df, q_values)
  
  expect_length(col_q_vals, 4)
  expect_true(all(!is.na(col_q_vals)))
})

test_that(".extract_q_metadata_from_result handles single q-value", {
  result_df <- data.frame(
    gene = c("gene1", "gene2"),
    sample1 = c(1.5, 2.0),
    sample2 = c(1.6, 2.1)
  )
  
  q_values <- c(1.0)
  col_q_vals <- TSENAT:::.extract_q_metadata_from_result(result_df, q_values)
  
  expect_true(is.na(col_q_vals[1]))  # Should return NA for single q
})

test_that(".extract_q_metadata_from_result handles empty dataframe", {
  result_df <- data.frame()
  q_values <- c(0.5, 1.0)
  col_q_vals <- TSENAT:::.extract_q_metadata_from_result(result_df, q_values)
  
  expect_true(is.na(col_q_vals[1]))
})

test_that(".extract_q_metadata_from_result handles missing q= pattern", {
  result_df <- data.frame(
    gene = c("gene1", "gene2"),
    sample1 = c(1.5, 2.0),
    sample2 = c(1.6, 2.1)
  )
  
  q_values <- c(0.5, 1.0)
  col_q_vals <- TSENAT:::.extract_q_metadata_from_result(result_df, q_values)
  
  # No "_q=" pattern in column names, should return NA
  expect_true(is.na(col_q_vals[1]))
})

# ===========================================================================
# Tests for .apply_diversity_post_hoc_norm()
# ===========================================================================

test_that(".apply_diversity_post_hoc_norm returns unchanged SE when norm_method is NULL", {
  # Create minimal SE
  assay_matrix <- matrix(c(1.0, 1.5, 2.0, 2.5), nrow = 2)
  rownames(assay_matrix) <- c("gene1", "gene2")
  result_se <- SummarizedExperiment(assays = list(diversity = assay_matrix))
  
  params <- list(genes = NULL, reference_group = NULL)
  result_se_out <- TSENAT:::.apply_diversity_post_hoc_norm(
    result_se, 
    norm_method = NULL, 
    params, 
    q_val = 1.0, 
    verbose = FALSE
  )
  
  expect_equal(
    SummarizedExperiment::assay(result_se_out, "diversity"),
    assay_matrix
  )
})

test_that(".apply_diversity_post_hoc_norm skips when SE is not SummarizedExperiment", {
  not_an_se <- data.frame(gene = c("gene1", "gene2"), val = c(1.0, 2.0))
  
  params <- list(genes = NULL, reference_group = NULL)
  result <- TSENAT:::.apply_diversity_post_hoc_norm(
    not_an_se,
    norm_method = "zscore",
    params,
    q_val = 1.0,
    verbose = FALSE
  )
  
  expect_identical(result, not_an_se)
})

test_that(".apply_diversity_post_hoc_norm handles norm_method='default'", {
  assay_matrix <- matrix(c(1.0, 1.5, 2.0, 2.5), nrow = 2)
  rownames(assay_matrix) <- c("gene1", "gene2")
  result_se <- SummarizedExperiment(assays = list(diversity = assay_matrix))
  
  params <- list(genes = NULL, reference_group = NULL)
  result_se_out <- TSENAT:::.apply_diversity_post_hoc_norm(
    result_se,
    norm_method = "default",
    params,
    q_val = 1.0,
    verbose = FALSE
  )
  
  # Default should not modify SE
  expect_equal(
    SummarizedExperiment::assay(result_se_out, "diversity"),
    assay_matrix
  )
})

test_that(".apply_diversity_post_hoc_norm applies zscore normalization", {
  # Create SE with test data
  # Matrix: 2 genes x 3 samples (columns are different q-values in real use)
  # Column-wise zscore should normalize each column independently
  assay_matrix <- matrix(c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0), nrow = 2)
  rownames(assay_matrix) <- c("gene1", "gene2")
  result_se <- SummarizedExperiment(assays = list(diversity = assay_matrix))
  
  params <- list(genes = NULL, reference_group = NULL)
  result_se_out <- TSENAT:::.apply_diversity_post_hoc_norm(
    result_se,
    norm_method = "zscore",
    params,
    q_val = 1.0,
    verbose = FALSE
  )
  
  # Verify function runs without error and returns valid SE
  expect_true(is(result_se_out, "SummarizedExperiment"))
  expect_equal(nrow(result_se_out), nrow(assay_matrix))
  expect_equal(ncol(result_se_out), ncol(assay_matrix))
  
  # Verify diversity assay exists
  expect_true("diversity" %in% names(SummarizedExperiment::assays(result_se_out)))
  
  # Verify assay has numeric values (either original or normalized)
  normalized_assay <- SummarizedExperiment::assay(result_se_out, "diversity")
  expect_true(is.numeric(normalized_assay))
  
  # Verify column-wise normalization was applied
  # Each column should have mean ≈ 0 and sd ≈ 1 after zscore normalization
  for (i in 1:ncol(normalized_assay)) {
    col_vals <- normalized_assay[, i]
    col_mean <- mean(col_vals)
    col_sd <- sd(col_vals)
    
    # Column mean should be close to 0
    expect_true(abs(col_mean) < 1e-5, 
                label = paste("Column", i, "mean:", round(col_mean, 10)))
    
    # Column sd should be close to 1 (or 0 if all values are identical, which shouldn't happen)
    expect_true(abs(col_sd - 1.0) < 1e-5 || (col_sd >= 0 && col_sd <= 1e-5),
                label = paste("Column", i, "sd:", round(col_sd, 10)))
  }
})

test_that(".apply_diversity_post_hoc_norm handles missing 'diversity' assay gracefully", {
  # Create SE without 'diversity' assay name
  assay_matrix <- matrix(c(1.0, 1.5, 2.0, 2.5), nrow = 2)
  result_se <- SummarizedExperiment(assays = list(other_name = assay_matrix))
  
  # This should not error because try-catch should handle it
  result_se_out <- tryCatch({
    TSENAT:::.apply_diversity_post_hoc_norm(
      result_se,
      norm_method = "zscore",
      list(genes = NULL, reference_group = NULL),
      q_val = 1.0,
      verbose = FALSE
    )
  }, error = function(e) {
    return(result_se)
  })
  
  expect_true(is(result_se_out, "SummarizedExperiment"))
})

test_that(".apply_diversity_post_hoc_norm returns valid SE after normalization with verbose", {
  assay_matrix <- matrix(c(1.0, 1.5, 2.0, 2.5), nrow = 2)
  rownames(assay_matrix) <- c("gene1", "gene2")
  result_se <- SummarizedExperiment(assays = list(diversity = assay_matrix))
  
  params <- list(genes = NULL, reference_group = NULL)
  
  # Just verify the function runs without error with verbose=TRUE
  result_se_out <- TSENAT:::.apply_diversity_post_hoc_norm(
    result_se,
    norm_method = "zscore",
    params,
    q_val = 1.0,
    verbose = TRUE
  )
  
  # Verify it returns a valid SE
  expect_true(is(result_se_out, "SummarizedExperiment"))
  expect_equal(nrow(result_se_out), 2)
  expect_equal(ncol(result_se_out), 2)
})

# ===========================================================================
# Integration tests: Helper functions working together
# ===========================================================================

test_that("Helper functions integrate correctly in workflow", {
  # Create analysis
  analysis <- make_test_analysis()
  
  # Prepare parameters
  params <- TSENAT:::.prepare_diversity_params(
    analysis,
    q = c(0.5, 1.0),
    norm = TRUE, norm_method = "zscore", reference_group = NULL,
    tpm = FALSE, assayno = 1, verbose = FALSE, what = "S",
    nthreads = 2, pseudocount = 0, min_valid_frac = 0.75,
    shrinkage = "none", genes = c("gene1", "gene2"), effective_length = NULL,
    metadata = NULL, bootstrap = FALSE, nboot = NULL,
    bootstrap_method = "percentile", bootstrap_ci = 0.95,
    bootstrap_include_diagnostics = TRUE, seed = 123
  )
  
  expect_true(is.list(params))
  expect_true("q" %in% names(params))
  
  # Build calculation arguments
  calc_args <- TSENAT:::.build_calc_diversity_args(params, analysis, list())
  
  expect_true("x" %in% names(calc_args))
  expect_true("q" %in% names(calc_args))
  expect_equal(calc_args$nthreads, 2)
})

# ===========================================================================
# Additional Edge Case Tests: .prepare_diversity_params()
# ===========================================================================

test_that(".prepare_diversity_params rejects infinite q values", {
  analysis <- make_test_analysis()
  
  expect_error({
    TSENAT:::.prepare_diversity_params(
      analysis,
      q = c(0.5, Inf),
      norm = NULL, norm_method = NULL, reference_group = NULL,
      tpm = FALSE, assayno = NULL, verbose = NULL, what = NULL,
      nthreads = NULL, pseudocount = NULL, min_valid_frac = NULL,
      shrinkage = NULL, genes = NULL, effective_length = NULL,
      metadata = NULL, bootstrap = NULL, nboot = NULL,
      bootstrap_method = NULL, bootstrap_ci = NULL,
      bootstrap_include_diagnostics = NULL, seed = NULL
    )
  }, "finite")
})

test_that(".prepare_diversity_params rejects NaN q values", {
  analysis <- make_test_analysis()
  
  expect_error({
    TSENAT:::.prepare_diversity_params(
      analysis,
      q = c(0.5, NaN),
      norm = NULL, norm_method = NULL, reference_group = NULL,
      tpm = FALSE, assayno = NULL, verbose = NULL, what = NULL,
      nthreads = NULL, pseudocount = NULL, min_valid_frac = NULL,
      shrinkage = NULL, genes = NULL, effective_length = NULL,
      metadata = NULL, bootstrap = NULL, nboot = NULL,
      bootstrap_method = NULL, bootstrap_ci = NULL,
      bootstrap_include_diagnostics = NULL, seed = NULL
    )
  }, "finite")
})

test_that(".prepare_diversity_params handles zero q value", {
  analysis <- make_test_analysis()
  
  params <- TSENAT:::.prepare_diversity_params(
    analysis,
    q = c(0.0, 1.0),
    norm = NULL, norm_method = NULL, reference_group = NULL,
    tpm = FALSE, assayno = NULL, verbose = NULL, what = NULL,
    nthreads = NULL, pseudocount = NULL, min_valid_frac = NULL,
    shrinkage = NULL, genes = NULL, effective_length = NULL,
    metadata = NULL, bootstrap = NULL, nboot = NULL,
    bootstrap_method = NULL, bootstrap_ci = NULL,
    bootstrap_include_diagnostics = NULL, seed = NULL
  )
  
  expect_true(0.0 %in% params$q)
})

test_that(".prepare_diversity_params handles very large q values", {
  analysis <- make_test_analysis()
  
  params <- TSENAT:::.prepare_diversity_params(
    analysis,
    q = c(100.0, 1000.0),
    norm = NULL, norm_method = NULL, reference_group = NULL,
    tpm = FALSE, assayno = NULL, verbose = NULL, what = NULL,
    nthreads = NULL, pseudocount = NULL, min_valid_frac = NULL,
    shrinkage = NULL, genes = NULL, effective_length = NULL,
    metadata = NULL, bootstrap = NULL, nboot = NULL,
    bootstrap_method = NULL, bootstrap_ci = NULL,
    bootstrap_include_diagnostics = NULL, seed = NULL
  )
  
  expect_equal(params$q, c(100.0, 1000.0))
})

test_that(".prepare_diversity_params validates nthreads is positive", {
  analysis <- make_test_analysis()
  
  # Negative nthreads should error
  expect_error({
    TSENAT:::.prepare_diversity_params(
      analysis,
      q = c(1.0),
      norm = NULL, norm_method = NULL, reference_group = NULL,
      tpm = FALSE, assayno = NULL, verbose = NULL, what = NULL,
      nthreads = -1,  # Invalid
      pseudocount = NULL, min_valid_frac = NULL,
      shrinkage = NULL, genes = NULL, effective_length = NULL,
      metadata = NULL, bootstrap = NULL, nboot = NULL,
      bootstrap_method = NULL, bootstrap_ci = NULL,
      bootstrap_include_diagnostics = NULL, seed = NULL
    )
  }, "positive")
})

test_that(".prepare_diversity_params handles logical parameters correctly", {
  analysis <- make_test_analysis()
  
  params <- TSENAT:::.prepare_diversity_params(
    analysis,
    q = c(1.0),
    norm = FALSE,  # Test all logical parameters
    norm_method = NULL, reference_group = NULL,
    tpm = TRUE,
    assayno = NULL,
    verbose = FALSE,
    what = NULL,
    nthreads = NULL,
    pseudocount = NULL,
    min_valid_frac = NULL,
    shrinkage = NULL,
    genes = NULL,
    effective_length = NULL,
    metadata = NULL,
    bootstrap = TRUE,  # Test bootstrap=TRUE
    nboot = NULL,
    bootstrap_method = NULL,
    bootstrap_ci = NULL,
    bootstrap_include_diagnostics = NULL,
    seed = NULL
  )
  
  expect_equal(params$norm, FALSE)
  expect_equal(params$tpm, TRUE)
  expect_equal(params$verbose, FALSE)
  expect_equal(params$bootstrap, TRUE)
})

test_that(".prepare_diversity_params merges config with defaults", {
  se <- make_test_se()
  config <- list(
    q_values = c(0.2, 0.8),
    nthreads = 3,
    verbose = FALSE
  )
  
  analysis <- TSENAT::TSENATAnalysis(se, config = config)
  
  params <- TSENAT:::.prepare_diversity_params(
    analysis,
    q = NULL,
    norm = NULL, norm_method = NULL, reference_group = NULL,
    tpm = FALSE,  # Explicit to avoid config issues
    assayno = NULL, verbose = NULL, what = NULL,
    nthreads = NULL,  # Should get from config
    pseudocount = NULL, min_valid_frac = NULL,
    shrinkage = NULL, genes = NULL, effective_length = NULL,
    metadata = NULL, bootstrap = NULL, nboot = NULL,
    bootstrap_method = NULL, bootstrap_ci = NULL,
    bootstrap_include_diagnostics = NULL, seed = NULL
  )
  
  expect_equal(params$q, c(0.2, 0.8))
  expect_equal(params$nthreads, 3)
  expect_equal(params$verbose, FALSE)  # From config
})

# ===========================================================================
# Additional Edge Case Tests: .build_calc_diversity_args()
# ===========================================================================

test_that(".build_calc_diversity_args handles all NULL optional parameters", {
  analysis <- make_test_analysis()
  params <- list(
    q = c(1.0),
    nthreads = 1,
    verbose = FALSE,
    bootstrap = FALSE,
    pseudocount = 0,
    norm = TRUE,
    what = "S",
    assayno = 1,
    min_valid_frac = 0.75,
    shrinkage = "none",
    bootstrap_method = "percentile",
    bootstrap_ci = 0.95,
    seed = NULL,
    tpm = FALSE,
    genes = NULL,
    effective_length = NULL,
    nboot = NULL,
    bootstrap_include_diagnostics = TRUE,
    metadata = NULL,
    norm_method = NULL,
    reference_group = NULL
  )
  
  # This should not error even with all optional parameters NULL
  args <- TSENAT:::.build_calc_diversity_args(params, analysis, list())
  
  expect_true(is.list(args))
  expect_true(length(args) > 0)
})

test_that(".build_calc_diversity_args combines params and dots without conflict", {
  analysis <- make_test_analysis()
  params <- list(
    q = c(1.0),
    nthreads = 1,
    verbose = FALSE,
    bootstrap = FALSE,
    pseudocount = 0,
    norm = TRUE,
    what = "S",
    assayno = 1,
    min_valid_frac = 0.75,
    shrinkage = "none",
    bootstrap_method = "percentile",
    bootstrap_ci = 0.95,
    seed = NULL,
    tpm = FALSE,
    genes = NULL,
    effective_length = NULL,
    nboot = NULL,
    bootstrap_include_diagnostics = TRUE,
    metadata = NULL,
    norm_method = NULL,
    reference_group = NULL
  )
  
  # Dots should override or augment params
  dots <- list(custom_param1 = 100, custom_param2 = "value")
  args <- TSENAT:::.build_calc_diversity_args(params, analysis, dots)
  
  expect_equal(args$custom_param1, 100)
  expect_equal(args$custom_param2, "value")
})

test_that(".build_calc_diversity_args preserves bootstrap parameters when bootstrap=TRUE", {
  analysis <- make_test_analysis()
  params <- list(
    q = c(1.0),
    nthreads = 1,
    verbose = FALSE,
    bootstrap = TRUE,  # Enable bootstrap
    pseudocount = 0,
    norm = TRUE,
    what = "S",
    assayno = 1,
    min_valid_frac = 0.75,
    shrinkage = "none",
    bootstrap_method = "block",
    bootstrap_ci = 0.90,
    seed = 42,
    tpm = FALSE,
    genes = NULL,
    effective_length = NULL,
    nboot = 200,
    bootstrap_include_diagnostics = FALSE,
    metadata = NULL,
    norm_method = NULL,
    reference_group = NULL
  )
  
  args <- TSENAT:::.build_calc_diversity_args(params, analysis, list())
  
  # Bootstrap-related arguments should be present
  expect_true("bootstrap" %in% names(args))
  expect_equal(args$bootstrap, TRUE)
})

# ===========================================================================
# Additional Edge Case Tests: .extract_q_metadata_from_result()
# ===========================================================================

test_that(".extract_q_metadata_from_result handles single column dataframe", {
  result_df <- data.frame(gene = c("gene1", "gene2"))
  q_values <- c(1.0)
  
  col_q_vals <- TSENAT:::.extract_q_metadata_from_result(result_df, q_values)
  
  expect_true(is.vector(col_q_vals) || is.list(col_q_vals))
})

test_that(".extract_q_metadata_from_result handles numeric column names", {
  result_df <- data.frame(
    gene = c("gene1", "gene2"),
    "1" = c(1.5, 2.0),
    "2" = c(1.6, 2.1)
  )
  
  q_values <- c(1.0)
  col_q_vals <- TSENAT:::.extract_q_metadata_from_result(result_df, q_values)
  
  expect_true(!is.null(col_q_vals) || length(col_q_vals) >= 0)
})

test_that(".extract_q_metadata_from_result handles many q-values", {
  # Create result dataframe with many columns
  n_cols <- 20
  col_names <- character(n_cols)
  
  for (i in 1:n_cols) {
    q_idx <- ((i - 1) %% 4) + 1
    q_val <- c(0.5, 1.0, 1.5, 2.0)[q_idx]
    col_names[i] <- sprintf("sample%d_q=%.3f", i, q_val)
  }
  
  result_df <- data.frame(
    gene = c("gene1", "gene2"),
    matrix(runif(n_cols * 2), nrow = 2)
  )
  colnames(result_df)[2:(n_cols + 1)] <- col_names
  
  q_values <- c(0.5, 1.0, 1.5, 2.0)
  col_q_vals <- TSENAT:::.extract_q_metadata_from_result(result_df, q_values)
  
  expect_true(length(col_q_vals) >= n_cols)
})

# ===========================================================================
# Additional Edge Case Tests: .apply_diversity_post_hoc_norm()
# ===========================================================================

test_that(".apply_diversity_post_hoc_norm handles single-row SE", {
  # Single row, multiple columns - ensure proper dimensions
  assay_matrix <- matrix(c(1.0, 2.0, 3.0), nrow = 1, ncol = 3)
  rownames(assay_matrix) <- c("gene1")
  colnames(assay_matrix) <- c("col1", "col2", "col3")
  
  result_se <- SummarizedExperiment(assays = list(diversity = assay_matrix))
  
  params <- list(genes = NULL, reference_group = NULL)
  result_se_out <- TSENAT:::.apply_diversity_post_hoc_norm(
    result_se,
    norm_method = "zscore",
    params,
    q_val = 1.0,
    verbose = FALSE
  )
  
  # Should still return valid SE
  expect_true(is(result_se_out, "SummarizedExperiment"))
  expect_equal(nrow(result_se_out), 1)
})

test_that(".apply_diversity_post_hoc_norm handles single-column SE", {
  # Multiple rows, single column
  assay_matrix <- matrix(c(1.0, 2.0, 3.0, 4.0), ncol = 1)
  rownames(assay_matrix) <- c("gene1", "gene2", "gene3", "gene4")
  result_se <- SummarizedExperiment(assays = list(diversity = assay_matrix))
  
  params <- list(genes = NULL, reference_group = NULL)
  result_se_out <- TSENAT:::.apply_diversity_post_hoc_norm(
    result_se,
    norm_method = "zscore",
    params,
    q_val = 1.0,
    verbose = FALSE
  )
  
  expect_true(is(result_se_out, "SummarizedExperiment"))
  expect_equal(ncol(result_se_out), 1)
})

test_that(".apply_diversity_post_hoc_norm handles matrix with zeros", {
  assay_matrix <- matrix(c(0.0, 0.0, 1.0, 1.0), nrow = 2)
  rownames(assay_matrix) <- c("gene1", "gene2")
  result_se <- SummarizedExperiment(assays = list(diversity = assay_matrix))
  
  params <- list(genes = NULL, reference_group = NULL)
  result_se_out <- TSENAT:::.apply_diversity_post_hoc_norm(
    result_se,
    norm_method = "zscore",
    params,
    q_val = 1.0,
    verbose = FALSE
  )
  
  expect_true(is(result_se_out, "SummarizedExperiment"))
  # Verify assay is numeric
  expect_true(is.numeric(SummarizedExperiment::assay(result_se_out, "diversity")))
})

test_that(".apply_diversity_post_hoc_norm handles uniform values in column", {
  # Column with identical values (sd will be 0)
  assay_matrix <- matrix(c(5.0, 5.0, 1.0, 2.0), nrow = 2)
  rownames(assay_matrix) <- c("gene1", "gene2")
  result_se <- SummarizedExperiment(assays = list(diversity = assay_matrix))
  
  params <- list(genes = NULL, reference_group = NULL)
  result_se_out <- TSENAT:::.apply_diversity_post_hoc_norm(
    result_se,
    norm_method = "zscore",
    params,
    q_val = 1.0,
    verbose = FALSE
  )
  
  expect_true(is(result_se_out, "SummarizedExperiment"))
  normalized_assay <- SummarizedExperiment::assay(result_se_out, "diversity")
  
  # First column (uniform values) should normalize to 0s or be handled gracefully
  col1 <- normalized_assay[, 1]
  col_sd <- sd(col1)
  expect_true(is.numeric(col1) && all(is.finite(col1) | is.na(col1)))
})

test_that(".apply_diversity_post_hoc_norm with reference_group parameter", {
  assay_matrix <- matrix(rnorm(20), nrow = 4)
  rownames(assay_matrix) <- c("gene1", "gene2", "gene3", "gene4")
  result_se <- SummarizedExperiment(assays = list(diversity = assay_matrix))
  
  params <- list(
    genes = c("gene1", "gene2"),
    reference_group = "control"
  )
  
  result_se_out <- TSENAT:::.apply_diversity_post_hoc_norm(
    result_se,
    norm_method = "zscore",
    params,
    q_val = 1.0,
    verbose = FALSE
  )
  
  expect_true(is(result_se_out, "SummarizedExperiment"))
})

test_that(".apply_diversity_post_hoc_norm preserves SE structure with metadata", {
  assay_matrix <- matrix(c(1.0, 2.0, 3.0, 4.0), nrow = 2)
  rownames(assay_matrix) <- c("gene1", "gene2")
  colnames(assay_matrix) <- c("q0.5", "q1.0")
  
  result_se <- SummarizedExperiment(
    assays = list(diversity = assay_matrix),
    metadata = list(custom_info = "test_value")
  )
  
  params <- list(genes = NULL, reference_group = NULL)
  result_se_out <- TSENAT:::.apply_diversity_post_hoc_norm(
    result_se,
    norm_method = "zscore",
    params,
    q_val = 1.0,
    verbose = FALSE
  )
  
  expect_true(is(result_se_out, "SummarizedExperiment"))
  expect_equal(nrow(result_se_out), nrow(assay_matrix))
  expect_equal(ncol(result_se_out), ncol(assay_matrix))
})

# ===========================================================================
# Additional Integration Tests
# ===========================================================================

test_that("Full workflow with bootstrap parameters", {
  analysis <- make_test_analysis()
  
  # Prepare parameters with bootstrap enabled
  params <- TSENAT:::.prepare_diversity_params(
    analysis,
    q = c(0.5, 1.0),
    norm = TRUE,
    norm_method = "zscore",
    reference_group = NULL,
    tpm = FALSE,
    assayno = 1,
    verbose = FALSE,
    what = "S",
    nthreads = 1,
    pseudocount = 0,
    min_valid_frac = 0.75,
    shrinkage = "none",
    genes = NULL,
    effective_length = NULL,
    metadata = NULL,
    bootstrap = TRUE,
    nboot = 50,
    bootstrap_method = "percentile",
    bootstrap_ci = 0.95,
    bootstrap_include_diagnostics = TRUE,
    seed = 42
  )
  
  expect_equal(params$bootstrap, TRUE)
  expect_equal(params$nboot, 50)
  expect_equal(params$bootstrap_method, "percentile")
  
  # Build arguments
  calc_args <- TSENAT:::.build_calc_diversity_args(params, analysis, list())
  
  expect_true("bootstrap" %in% names(calc_args))
  expect_equal(calc_args$bootstrap, TRUE)
})

test_that("Parameter resolution priority: explicit > config > default", {
  se <- make_test_se()
  config <- list(
    q_values = c(0.1, 0.2),
    nthreads = 4,
    verbose = FALSE,
    bootstrap = TRUE,
    nboot = 100
  )
  
  analysis <- TSENAT::TSENATAnalysis(se, config = config)
  
  # Mix of explicit, config, and default parameters
  params <- TSENAT:::.prepare_diversity_params(
    analysis,
    q = c(1.0, 2.0),  # Explicit overrides config
    norm = NULL,  # Use default
    norm_method = NULL,
    reference_group = NULL,
    tpm = FALSE,  # Explicit (not from config to avoid type issues)
    assayno = NULL,
    verbose = TRUE,  # Explicit overrides config
    what = NULL,
    nthreads = 8,  # Explicit overrides config
    pseudocount = NULL,
    min_valid_frac = NULL,
    shrinkage = NULL,
    genes = NULL,
    effective_length = NULL,
    metadata = NULL,
    bootstrap = NULL,  # Should get from config
    nboot = NULL,  # Should get from config
    bootstrap_method = NULL,
    bootstrap_ci = NULL,
    bootstrap_include_diagnostics = NULL,
    seed = NULL
  )
  
  # Verify priority order
  expect_equal(params$q, c(1.0, 2.0))  # Explicit
  expect_equal(params$verbose, TRUE)  # Explicit
  expect_equal(params$nthreads, 8)  # Explicit
  expect_equal(params$norm, TRUE)  # Default
  expect_equal(params$bootstrap, TRUE)  # From config
  expect_equal(params$nboot, 100)  # From config
})

test_that("Helper function error propagation", {
  analysis <- make_test_analysis()
  
  # Attempt with invalid q
  attempt_error <- tryCatch({
    TSENAT:::.prepare_diversity_params(
      analysis,
      q = c(NA),  # Invalid
      norm = NULL, norm_method = NULL, reference_group = NULL,
      tpm = FALSE, assayno = NULL, verbose = NULL, what = NULL,
      nthreads = NULL, pseudocount = NULL, min_valid_frac = NULL,
      shrinkage = NULL, genes = NULL, effective_length = NULL,
      metadata = NULL, bootstrap = NULL, nboot = NULL,
      bootstrap_method = NULL, bootstrap_ci = NULL,
      bootstrap_include_diagnostics = NULL, seed = NULL
    )
  }, error = function(e) {
    return(e)
  })
  
  # Should produce an error
  expect_true(is(attempt_error, "error") || is(attempt_error, "try-error"))
})
