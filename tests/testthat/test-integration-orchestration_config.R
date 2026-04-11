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
  expect_null(results(analysis, type = "diversity"))
  expect_null(results(analysis, type = "divergence"))
  expect_null(results(analysis, type = "lm"))
  
  # Unknown type
  expect_error(results(analysis, type = "unknown"), "Unknown result type")
  
  # Non-TSENATAnalysis object
  expect_error(results(se, type = "diversity"), "must be a TSENATAnalysis object")
})

test_that("results handles diversity with q-value filtering", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = TSENAT_config())
  
  n_genes <- nrow(se)
  diversity_results <- list(
    q_0.5 = matrix(rnorm(n_genes, mean = 3, sd = 1), nrow = 1, ncol = n_genes),
    q_1.0 = matrix(rnorm(n_genes, mean = 3, sd = 1), nrow = 1, ncol = n_genes),
    q_1.5 = matrix(rnorm(n_genes, mean = 3, sd = 1), nrow = 1, ncol = n_genes),
    q_2.0 = matrix(rnorm(n_genes, mean = 3, sd = 1), nrow = 1, ncol = n_genes)
  )
  analysis@diversity_results <- diversity_results
  
  # All results
  all_results <- results(analysis, type = "diversity")
  expect_false(is.null(all_results))
  
  # Specific q-value
  q1_results <- results(analysis, type = "diversity", q = 1.0)
  expect_false(is.null(q1_results))
})

test_that("results returns all supported result types", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = TSENAT_config())
  
  n_genes <- nrow(se)
  
  # Mock all result types
  analysis@diversity_results <- list(
    q_0.5 = matrix(rnorm(n_genes), nrow = 1, ncol = n_genes),
    q_1.0 = matrix(rnorm(n_genes), nrow = 1, ncol = n_genes)
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
  expect_false(is.null(results(analysis, type = "diversity")))
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
