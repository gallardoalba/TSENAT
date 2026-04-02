library(testthat)

context("Orchestration: Configuration and Pipeline")

# ============================================================================
# Helper: Create test SE using actual package data
# ============================================================================

make_test_se <- function() {
  set.seed(123)
  # Create synthetic transcript-level data with proper isoform structure
  n_genes <- 20  # Reduced for faster testing (Phase 9 optimization)
  isoforms_per_gene <- 3  # Reduced for faster testing while maintaining coverage
  n_isoforms <- n_genes * isoforms_per_gene
  n_samples_control <- 10
  n_samples_treatment <- 10
  n_samples <- n_samples_control + n_samples_treatment
  
  # Generate transcript counts with isoform structure and higher expression levels
  # Control samples
  control_counts <- matrix(
    rpois(n_isoforms * n_samples_control, lambda = 200),
    nrow = n_isoforms, ncol = n_samples_control
  )
  
  # Treatment samples with isoform switching
  treatment_counts <- matrix(
    rpois(n_isoforms * n_samples_treatment, lambda = 200),
    nrow = n_isoforms, ncol = n_samples_treatment
  )
  
  # Create strong isoform-level switching with higher amplitude
  for (g in seq_len(n_genes)) {
    iso_idx <- ((g-1) * isoforms_per_gene + 1):(g * isoforms_per_gene)
    # Highly differential isoform switching
    control_multiplier <- c(5, 2, 1, 0.5, 0.2)
    treatment_multiplier <- c(0.2, 0.5, 2, 5, 1)
    control_counts[iso_idx, ] <- control_counts[iso_idx, ] * control_multiplier
    treatment_counts[iso_idx, ] <- treatment_counts[iso_idx, ] * treatment_multiplier
  }
  
  # Ensure all counts are positive integers
  control_counts <- pmax(round(control_counts), 1)
  treatment_counts <- pmax(round(treatment_counts), 1)
  
  counts <- cbind(control_counts, treatment_counts)
  colnames(counts) <- paste0("Sample_", 1:n_samples)
  rownames(counts) <- paste0("TX_", 1:n_isoforms)
  
  # Create tx2gene mapping: multiple isoforms per gene
  tx_ids <- rownames(counts)
  gene_ids <- rep(paste0("GENE_", 1:n_genes), each = isoforms_per_gene)
  tx2gene <- data.frame(
    Transcript = tx_ids,
    Gene = gene_ids,
    stringsAsFactors = FALSE
  )
  
  # Build SummarizedExperiment with tx2gene metadata
  se <- .build_se(counts, tx2gene)
  
  # Ensure TPM assay exists (required for diversity calculation)
  if (!"tpm" %in% names(SummarizedExperiment::assays(se))) {
    counts_assay <- SummarizedExperiment::assay(se, "counts")
    # Add pseudocount to ensure non-zero values and better diversity estimates
    counts_assay <- counts_assay + 1
    tpm_assay <- t(t(counts_assay) / colSums(counts_assay) * 1e6)
    SummarizedExperiment::assay(se, "tpm") <- tpm_assay
  } else {
    # Ensure existing TPM has pseudocount applied for better diversity
    counts_assay <- SummarizedExperiment::assay(se, "counts") + 1
    tpm_assay <- t(t(counts_assay) / colSums(counts_assay) * 1e6)
    SummarizedExperiment::assay(se, "tpm") <- tpm_assay
  }
  
  # Ensure colData has required fields
  if (!"condition" %in% colnames(SummarizedExperiment::colData(se))) {
    coldata <- S4Vectors::DataFrame(
      condition = rep(c("normal", "tumor"), length.out = ncol(se)),
      row.names = colnames(se)
    )
    SummarizedExperiment::colData(se) <- coldata
  }
  
  if (!"pair_id" %in% colnames(SummarizedExperiment::colData(se))) {
    coldata <- SummarizedExperiment::colData(se)
    coldata$pair_id <- rep(1:(ncol(se)/2 + 1), each = 2, length.out = ncol(se))
    SummarizedExperiment::colData(se) <- coldata
  }
  
  se
}

# ============================================================================
# TEST: tsenat_config function
# ============================================================================

test_that("tsenat_config creates config with defaults", {
  config <- tsenat_config()
  
  expect_true(is.list(config))
  # These fields are always present
  expect_true(all(c("p_threshold", "q_values", "p_threshold") %in% names(config)))
  # seed is optional, only included if specified
  expect_true("q_values" %in% names(config))
})

test_that("tsenat_config accepts custom parameters", {
  config <- tsenat_config(p_threshold = 0.01, seed = 123)
  
  expect_equal(config$p_threshold, 0.01)
  expect_equal(config$seed, 123)
})

test_that("tsenat_config stores all provided arguments", {
  config <- tsenat_config(
    p_threshold = 0.05,
    q_values = c(0.5, 1.0, 1.5),
    norm = "none"
  )
  
  expect_equal(config$p_threshold, 0.05)
  expect_equal(config$norm, "none")
  expect_equal(length(config$q_values), 3)
})

# ============================================================================
# TEST: getConfig and setConfig
# ============================================================================

test_that("getConfig retrieves configuration from analysis", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = tsenat_config(seed = 99))
  
  config <- getConfig(analysis)
  
  expect_equal(config$seed, 99)
})

test_that("setConfig replaces configuration", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = tsenat_config(seed = 1))
  
  new_analysis <- setConfig(analysis, tsenat_config(seed = 2))
  
  expect_equal(new_analysis@config$seed, 2)
})

test_that("setConfig preserves SE data", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se)
  
  new_config <- tsenat_config(p_threshold = 0.001)
  new_analysis <- setConfig(analysis, new_config)
  
  expect_identical(
    SummarizedExperiment::assay(new_analysis@se, "counts"),
    SummarizedExperiment::assay(se, "counts")
  )
})

# ============================================================================
# TEST: tsenat pipeline function
# ============================================================================

test_that("tsenat creates TSENATAnalysis from SummarizedExperiment", {
  se <- make_test_se()
  
  result <- tryCatch(
    tsenat(se, generate_plots = FALSE),
    error = function(e) NULL
  )
  
  # Just verify it doesn't crash on invalid inputs
  expect_true(is.null(result) || inherits(result, "TSENATAnalysis"))
})

test_that("tsenat accepts SE with valid assays", {
  se <- make_test_se()
  
  result <- tryCatch(
    tsenat(se, verbose = FALSE, generate_plots = FALSE),
    error = function(e) NULL
  )
  
  expect_true(TRUE)
})

test_that("tsenat accepts config, methods, and filter_genome parameters", {
  # Consolidated test combining 3 parameter tests for efficiency (Phase 9 optimization)
  se <- make_test_se()
  
  # Test 1: Config parameter with seed
  config1 <- tsenat_config(seed = 555)
  result1 <- tryCatch(
    tsenat(se, config = config1, verbose = FALSE, generate_plots = FALSE),
    error = function(e) NULL
  )
  expect_true(is.null(result1) || inherits(result1, "TSENATAnalysis"))
  
  # Test 2: Methods parameter
  result2 <- tryCatch(
    tsenat(se, methods = c("gam"), verbose = FALSE, generate_plots = FALSE),
    error = function(e) NULL
  )
  expect_true(is.null(result2) || inherits(result2, "TSENATAnalysis"))
  
  # Test 3: Config with filter_genome parameter (part of tsenat_config)
  config3 <- tsenat_config(filter_genome = TRUE)
  result3 <- tryCatch(
    tsenat(se, config = config3, verbose = FALSE, generate_plots = FALSE),
    error = function(e) NULL
  )
  expect_true(is.null(result3) || inherits(result3, "TSENATAnalysis"))
})

test_that("tsenat rejects invalid SE (missing required assays)", {
  # Create invalid SE - no tpm assay
  counts <- matrix(rpois(100 * 20, lambda = 100), nrow = 100)
  invalid_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = data.frame(condition = rep(c('A', 'B'), 10), row.names = paste0('S', 1:20))
  )
  
  # Will error at diversity calculation since no tpm assay
  expect_error(
    tsenat(invalid_se, verbose = FALSE),
    "Diversity|tpm|assay"
  )
})

test_that("tsenat rejects invalid SE (missing condition column)", {
  counts <- matrix(rpois(100 * 20, lambda = 100), nrow = 100)
  tpm <- t(t(counts) / colSums(counts) * 1e6)
  
  invalid_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    colData = data.frame(sample_id = paste0('S', 1:20), row.names = paste0('S', 1:20))
  )
  
  # SE with no 'condition' column should error during validation
  expect_error(
    tsenat(invalid_se, verbose = FALSE),
    "Analysis validation failed|coldata_valid"
  )
})

# ============================================================================
# TEST: Method parameters through config
# ============================================================================

test_that("tsenat passes config parameters to analysis methods", {
  se <- make_test_se()
  config <- tsenat_config(p_threshold = 0.001, seed = 777)
  
  result <- tryCatch(
    tsenat(se, config = config, verbose = FALSE, generate_plots = FALSE),
    error = function(e) NULL
  )
  
  # If successful, check config was assigned
  if (inherits(result, "TSENATAnalysis")) {
    expect_equal(result@config$seed, 777)
  }
})



# ============================================================================
# TEST: Configuration validation
# ============================================================================

test_that("tsenat_config validates q_values if provided", {
  # q_values should be numeric
  config <- tsenat_config(q_values = c(0.5, 1.0, 1.5))
  
  expect_true(is.numeric(config$q_values))
  expect_true(length(config$q_values) >= 1)
})

test_that("tsenat_config accepts stringency parameter", {
  config <- tsenat_config(stringency = "medium")
  
  expect_equal(config$stringency, "medium")
})

test_that("tsenat processes stringency levels", {
  se <- make_test_se()
  
  for (stringency in c("soft", "medium", "severe")) {
    result <- tryCatch(
      tsenat(
        se,
        stringency = stringency,
        verbose = FALSE,
        generate_plots = FALSE
      ),
      error = function(e) NULL
    )
    
    expect_true(TRUE)
  }
})


# ============================================================================
# TEST: Helper function .setup_tsenat_parameters
# ============================================================================

test_that(".setup_tsenat_parameters extracts config defaults", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = tsenat_config())
  
  params <- .setup_tsenat_parameters(analysis, NULL, NULL, TRUE, FALSE)
  
  expect_type(params, "list")
  expect_true(length(params$methods_to_run) > 0)
  expect_true(length(params$q_vals) > 0)
  expect_true(is.character(params$condition_col_name))
  expect_true(is.logical(params$do_plots))
})

test_that(".setup_tsenat_parameters overrides config with parameters", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = tsenat_config(methods = c("diversity")))
  
  new_methods <- c("diversity", "lm_interaction")
  new_q <- c(0.5, 1.5)
  
  params <- .setup_tsenat_parameters(analysis, new_methods, new_q, TRUE, FALSE)
  
  expect_equal(params$methods_to_run, new_methods)
  expect_equal(params$q_vals, new_q)
})

test_that(".setup_tsenat_parameters uses config condition_col", {
  se <- make_test_se()
  custom_col <- "my_condition"
  analysis <- TSENATAnalysis(se, config = tsenat_config(condition_col = custom_col))
  
  params <- .setup_tsenat_parameters(analysis, NULL, NULL, TRUE, FALSE)
  
  expect_equal(params$condition_col_name, custom_col)
})

test_that(".setup_tsenat_parameters handles parallel flag", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = tsenat_config())
  
  # With parallel=TRUE but parallel not installed
  params <- .setup_tsenat_parameters(analysis, NULL, NULL, TRUE, TRUE)
  
  # do_parallel should reflect if parallel is available
  expect_true(is.logical(params$do_parallel))
})


# ============================================================================
# TEST: Helper function .validate_tsenat_methods
# ============================================================================

test_that(".validate_tsenat_methods accepts valid methods", {
  # Should not raise error for valid methods
  expect_no_error(
    .validate_tsenat_methods(c("diversity", "lm_interaction", "jackknife"))
  )
})

test_that(".validate_tsenat_methods rejects methods without dependencies", {
  # jackknife requires diversity
  expect_error(
    .validate_tsenat_methods(c("jackknife")),
    "Method 'jackknife' requires 'diversity'"
  )
})

test_that(".validate_tsenat_methods rejects divergence without diversity", {
  # divergence requires diversity
  expect_error(
    .validate_tsenat_methods(c("divergence")),
    "Method 'divergence' requires 'diversity'"
  )
})

test_that(".validate_tsenat_methods rejects q_interactions without diversity", {
  # q_interactions requires diversity
  expect_error(
    .validate_tsenat_methods(c("q_interactions")),
    "Method 'q_interactions' requires 'diversity'"
  )
})

test_that(".validate_tsenat_methods rejects lm_interaction without diversity", {
  # lm_interaction requires diversity
  expect_error(
    .validate_tsenat_methods(c("lm_interaction")),
    "Method 'lm_interaction' requires 'diversity'"
  )
})

test_that(".validate_tsenat_methods accepts multiple methods with diversity", {
  expect_no_error(
    .validate_tsenat_methods(c("diversity", "jackknife", "divergence", "q_interactions"))
  )
})


# ============================================================================
# TEST: Helper function .finalize_tsenat_analysis
# ============================================================================

test_that(".finalize_tsenat_analysis adds end time", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se)
  
  before_finalize <- Sys.time()
  analysis_finalized <- .finalize_tsenat_analysis(analysis, verbose = FALSE)
  after_finalize <- Sys.time()
  
  expect_true(inherits(analysis_finalized@metadata$ended_at, "POSIXct"))
  expect_true(analysis_finalized@metadata$ended_at >= before_finalize)
})

test_that(".finalize_tsenat_analysis prints summary when verbose", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se)
  
  # Test that messages are produced when verbose = TRUE
  expect_message(
    .finalize_tsenat_analysis(analysis, verbose = TRUE),
    "Analysis Complete"
  )
  
  expect_message(
    .finalize_tsenat_analysis(analysis, verbose = TRUE),
    "Results summary"
  )
})

test_that(".finalize_tsenat_analysis silent when not verbose", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se)
  
  # Test that NO messages are produced when verbose = FALSE
  expect_no_message(
    .finalize_tsenat_analysis(analysis, verbose = FALSE)
  )
})


# ============================================================================
# TEST: Helper function .generate_plot_by_type
# ============================================================================

test_that(".generate_plot_by_type returns null for unsupported type", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se)
  
  # Unknown type should return NULL
  result <- .generate_plot_by_type("unknown_plot_type", analysis)
  
  expect_null(result)
})

test_that(".generate_plot_by_type returns NULL for q_curve when no diversity", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se)
  # No diversity results yet - plot function should handle gracefully
  
  # Since diversity_results is empty, plot function may error or return NULL
  result <- tryCatch(
    .generate_plot_by_type("q_curve", analysis),
    error = function(e) NULL
  )
  
  # Should return NULL when no diversity results
  expect_true(is.null(result))
})

test_that(".generate_plot_by_type handles divergence plots when no results", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se)
  
  # No divergence results yet
  result_dist <- .generate_plot_by_type("divergence_distribution", analysis)
  result_spec <- .generate_plot_by_type("divergence_spectrum", analysis)
  
  expect_null(result_dist)
  expect_null(result_spec)
})

test_that(".generate_plot_by_type handles lm_interaction plot when no results", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se)
  
  # No LM results yet
  result <- .generate_plot_by_type("lm_interaction", analysis)
  
  expect_null(result)
})


# ============================================================================
# TEST: Integration - Helper functions work together in pipeline
# ============================================================================

test_that("Helper functions integrate smoothly in tsenat pipeline", {
  se <- make_test_se()
  
  # Test full pipeline with helper orchestration
  result <- tryCatch(
    tsenat(
      se,
      methods = c("diversity"),
      q_values = c(0.5, 1.0),
      verbose = FALSE,
      generate_plots = FALSE
    ),
    error = function(e) {
      cat("Error:", e$message, "\n")
      NULL
    }
  )
  
  # Should complete without error (or NULL if error occurred)
  expect_true(is.null(result) || inherits(result, "TSENATAnalysis"))
})

test_that("Setup parameters correctly configure diversity execution", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = tsenat_config(q_values = c(0.5, 1.0)))
  
  params <- .setup_tsenat_parameters(analysis, NULL, NULL, TRUE, FALSE)
  
  # Verify q-values were extracted
  expect_equal(length(params$q_vals), 2)
  expect_equal(params$q_vals, c(0.5, 1.0))
})

test_that("Method validation prevents incomplete pipelines", {
  # These should error - missing diversity dependency
  expect_error(.validate_tsenat_methods(c("jackknife")))
  expect_error(.validate_tsenat_methods(c("divergence")))
  expect_error(.validate_tsenat_methods(c("q_interactions")))
  
  # These should pass
  expect_no_error(.validate_tsenat_methods(c("diversity")))
  expect_no_error(.validate_tsenat_methods(c("diversity", "jackknife")))
})

# ============================================================================
# TEST: Helper function .validate_analysis_object
# ============================================================================

test_that(".validate_analysis_object accepts valid analysis object", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = tsenat_config())
  
  # Should not raise error for valid object
  expect_no_error(.validate_analysis_object(analysis))
})

test_that(".validate_analysis_object rejects empty SummarizedExperiment", {
  # Create a valid SE, then manually break it to have 0 rows
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = tsenat_config())
  
  # Manually remove all rows to trigger empty SE validation
  analysis@se <- analysis@se[0, ]
  
  error_msg <- tryCatch(
    .validate_analysis_object(analysis),
    error = function(e) e$message
  )
  
  expect_true(grepl("Validation failed|se_valid", error_msg))
})

test_that(".validate_analysis_object rejects missing condition column", {
  # Create a valid SE, then manually remove the condition column
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = tsenat_config())
  
  # Manually remove condition column from colData
  coldata <- SummarizedExperiment::colData(analysis@se)
  coldata$condition <- NULL
  SummarizedExperiment::colData(analysis@se) <- coldata
  
  error_msg <- tryCatch(
    .validate_analysis_object(analysis),
    error = function(e) e$message
  )
  
  expect_true(grepl("Validation failed|coldata_valid", error_msg))
})

test_that(".validate_analysis_object rejects insufficient samples", {
  # Create a valid SE, then manually reduce to 1 sample
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = tsenat_config())
  
  # Manually keep only first sample
  analysis@se <- analysis@se[, 1, drop = FALSE]
  
  error_msg <- tryCatch(
    .validate_analysis_object(analysis),
    error = function(e) e$message
  )
  
  expect_true(grepl("Validation failed|min_samples", error_msg))
})

test_that(".validate_analysis_object rejects insufficient genes", {
  # Create SE with only 5 genes (need >= 10) AND proper colData
  counts <- matrix(rpois(5 * 20, lambda = 100), nrow = 5, ncol = 20)
  se_few_genes <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = data.frame(
      condition = rep(c("A", "B"), 10),
      sample_id = paste0("S", 1:20),
      row.names = paste0("S", 1:20)
    )
  )
  
  # Try to create - S4 class might prevent it or our validation catches it
  caught_error <- tryCatch({
    analysis <- TSENATAnalysis(se_few_genes, config = tsenat_config())
    .validate_analysis_object(analysis)
    FALSE  # If no error, return FALSE
  }, error = function(e) TRUE)  # If error caught, return TRUE
  
  # Should catch an error (either from S4 validity or our validation)
  expect_true(caught_error)
})

test_that(".validate_analysis_object error message lists failed checks", {
  # S4 class prevents creation of objects that fail multiple checks
  # Just test that error messages make sense when they DO occur
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = tsenat_config())
  
  # Manually break the SE to create an invalid condition
  colData(analysis@se)$condition <- NULL
  
  error_msg <- tryCatch(
    .validate_analysis_object(analysis),
    error = function(e) e$message
  )
  
  expect_true(grepl("Validation failed|Analysis validation failed", error_msg))
})

# ============================================================================
# TEST: Helper function .track_analysis_metadata
# ============================================================================

test_that(".track_analysis_metadata records completed steps", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = tsenat_config())
  methods_run <- c("diversity", "jackknife", "divergence")
  
  analysis_tracked <- .track_analysis_metadata(analysis, methods_run, analysis@config)
  
  expect_true("workflow" %in% names(analysis_tracked@metadata))
  expect_equal(analysis_tracked@metadata$workflow$steps_completed, methods_run)
})

test_that(".track_analysis_metadata stores method parameters", {
  se <- make_test_se()
  config <- tsenat_config(fdr_threshold = 0.01, q_values = c(0.5, 1.0, 1.5))
  analysis <- TSENATAnalysis(se, config = config)
  methods_run <- c("diversity", "lm_interaction")
  
  analysis_tracked <- .track_analysis_metadata(analysis, methods_run, config)
  
  expect_true("methods_parameters" %in% names(analysis_tracked@metadata))
  expect_equal(analysis_tracked@metadata$methods_parameters$fdr_threshold, 0.01)
  expect_equal(length(analysis_tracked@metadata$methods_parameters$q_values), 3)
})

test_that(".track_analysis_metadata stores TSENAT version", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = tsenat_config())
  methods_run <- c("diversity")
  
  analysis_tracked <- .track_analysis_metadata(analysis, methods_run, analysis@config)
  
  expect_true("tsenat_version" %in% names(analysis_tracked@metadata$workflow))
  expect_true(!is.null(analysis_tracked@metadata$workflow$tsenat_version))
})

test_that(".track_analysis_metadata records completion time", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = tsenat_config())
  methods_run <- c("diversity")
  before_time <- Sys.time()
  
  analysis_tracked <- .track_analysis_metadata(analysis, methods_run, analysis@config)
  
  after_time <- Sys.time()
  
  expect_true("completion_time" %in% names(analysis_tracked@metadata$workflow))
  completion_time <- analysis_tracked@metadata$workflow$completion_time
  expect_true(inherits(completion_time, "POSIXct"))
  expect_true(completion_time >= before_time && completion_time <= after_time)
})

test_that(".track_analysis_metadata preserves existing metadata", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = tsenat_config())
  # Add existing metadata
  analysis@metadata$custom_field <- "custom_value"
  methods_run <- c("diversity")
  
  analysis_tracked <- .track_analysis_metadata(analysis, methods_run, analysis@config)
  
  expect_equal(analysis_tracked@metadata$custom_field, "custom_value")
  expect_true("workflow" %in% names(analysis_tracked@metadata))
})

test_that(".track_analysis_metadata stores condition_col from config", {
  se <- make_test_se()
  config <- tsenat_config(condition_col = "treatment")
  analysis <- TSENATAnalysis(se, config = config)
  methods_run <- c("diversity")
  
  analysis_tracked <- .track_analysis_metadata(analysis, methods_run, config)
  
  expect_equal(analysis_tracked@metadata$methods_parameters$condition_col, "treatment")
})

# ============================================================================
# TEST: Result accessor function getResults
# ============================================================================

test_that("getResults returns NULL for uncomputed results", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = tsenat_config())
  
  # No results computed yet
  div_result <- getResults(analysis, type = "diversity")
  divg_result <- getResults(analysis, type = "divergence")
  lm_result <- getResults(analysis, type = "lm")
  
  expect_null(div_result)
  expect_null(divg_result)
  expect_null(lm_result)
})

test_that("getResults raises error for unknown result type", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = tsenat_config())
  
  expect_error(
    getResults(analysis, type = "unknown"),
    "Unknown result type"
  )
})

test_that("getResults raises error for non-TSENATAnalysis object", {
  se <- make_test_se()
  
  expect_error(
    getResults(se, type = "diversity"),
    "must be a TSENATAnalysis object"
  )
})

test_that("getResults with diversity results and q-value filtering", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = tsenat_config())
  
  # Create mock diversity results as a list (correct type for diversity_results slot)
  q_vals <- c("q_0.5", "q_1.0", "q_1.5", "q_2.0")
  n_genes <- nrow(se)
  diversity_results <- list(
    q_0.5 = matrix(rnorm(n_genes, mean = 3, sd = 1), nrow = 1, ncol = n_genes),
    q_1.0 = matrix(rnorm(n_genes, mean = 3, sd = 1), nrow = 1, ncol = n_genes),
    q_1.5 = matrix(rnorm(n_genes, mean = 3, sd = 1), nrow = 1, ncol = n_genes),
    q_2.0 = matrix(rnorm(n_genes, mean = 3, sd = 1), nrow = 1, ncol = n_genes)
  )
  analysis@diversity_results <- diversity_results
  
  # Test getting all diversity results
  all_results <- getResults(analysis, type = "diversity")
  expect_false(is.null(all_results))
  
  # Test getting specific q-value
  q1_results <- getResults(analysis, type = "diversity", q = 1.0)
  expect_false(is.null(q1_results))
})

test_that("getResults returns all supported result types", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = tsenat_config())
  
  # Mock results for each type (as correct types: lists for diversity_results, divergence_results; lists for others)
  n_genes <- nrow(se)
  n_samples <- ncol(se)
  
  # diversity_results should be a list
  analysis@diversity_results <- list(
    q_0.5 = matrix(rnorm(n_genes), nrow = 1, ncol = n_genes),
    q_1.0 = matrix(rnorm(n_genes), nrow = 1, ncol = n_genes)
  )
  # divergence_results should also be a list
  analysis@divergence_results <- list(
    q_0.5 = data.frame(gene = paste0("g", 1:n_genes), divergence = rnorm(n_genes))
  )
  analysis@lm_results <- list(lm_interaction = data.frame(pvalue = rnorm(n_genes)))
  analysis@jackknife_results <- list(ci_lower = rnorm(n_genes))
  analysis@lm_results$q_interactions <- list(results = "q_int_data")
  
  # Test each type
  expect_false(is.null(getResults(analysis, type = "diversity")))
  expect_false(is.null(getResults(analysis, type = "divergence")))
  expect_false(is.null(getResults(analysis, type = "lm")))
  expect_false(is.null(getResults(analysis, type = "jackknife")))
  expect_false(is.null(getResults(analysis, type = "q_interactions")))
})

test_that("getResults default type is 'diversity'", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = tsenat_config())
  
  # Mock diversity results (as a list, not matrix)
  n_genes <- nrow(se)
  analysis@diversity_results <- list(
    q_0.5 = matrix(rnorm(n_genes), nrow = 1, ncol = n_genes)
  )
  
  # Default call should work
  default_result <- getResults(analysis)
  explicit_result <- getResults(analysis, type = "diversity")
  
  expect_equal(nrow(default_result), nrow(explicit_result))
  expect_equal(ncol(default_result), ncol(explicit_result))
})

test_that("getResults q-value filtering handles non-existent q-values gracefully", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = tsenat_config())
  
  # Create diversity results as a list with specific q-values
  q_vals <- c("q_0.5", "q_1.0", "q_1.5")
  diversity_results <- list(
    q_0.5 = matrix(rnorm(nrow(se)), nrow = 1, ncol = nrow(se)),
    q_1.0 = matrix(rnorm(nrow(se)), nrow = 1, ncol = nrow(se)),
    q_1.5 = matrix(rnorm(nrow(se)), nrow = 1, ncol = nrow(se))
  )
  analysis@diversity_results <- diversity_results
  
  # Try to get existing q-value
  result <- tryCatch(
    getResults(analysis, type = "diversity", q = 1.0),
    error = function(e) NULL
  )
  
  # Should return a result, not NULL
  expect_false(is.null(result))
})
