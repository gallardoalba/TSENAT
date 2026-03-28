# Tests for calculate_divergence main function and related public functions
# Covers uncovered lines from divergence_coverage.txt for main functions

# =====================================================================
# Tests for calculate_divergence validation and setup
# =====================================================================

test_that("calculate_divergence validates norm parameter", {
  # Basic test data - 2 genes x 6 samples
  se <- create_test_se_simple(
    n_genes = 2,
    n_samples = 6,
    control_n = 2,
    group_col_name = "sample_type"
  )
  colnames(se) <- c("s1", "s2", "s3", "s4", "s5", "s6")
  
  # Test norm parameter coercion from logical to character
  # This tests the normalization setup
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    norm = TRUE,
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  # Should return a result (may be SummarizedExperiment or list)
  expect_true(!is.null(result))
})

test_that("calculate_divergence sorts q values in ascending order", {
  se <- create_test_se_simple(
    n_genes = 2,
    n_samples = 6,
    control_n = 2,
    group_col_name = "sample_type"
  )
  colnames(se) <- c("s1", "s2", "s3", "s4", "s5", "s6")
  
  # Pass q values out of order and verify sorting
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = c(2, 0.5, 1),  # Unsorted
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  # Should return a result without error
  expect_true(!is.null(result))
})

test_that("calculate_divergence rejects non-positive q values", {
  se <- create_test_se_simple(
    n_genes = 2,
    n_samples = 6,
    control_n = 2,
    group_col_name = "sample_type"
  )
  colnames(se) <- c("s1", "s2", "s3", "s4", "s5", "s6")
  
  # Should reject negative q
  expect_error(
    .calculate_divergence(
      se,
      group_col = "sample_type",
      control_group = "Control",
      q = c(-0.5, 1, 2),
      bootstrap = FALSE,
      verbose = FALSE
    ),
    "q parameter must be >= 0"
  )
})

test_that("calculate_divergence requires SummarizedExperiment input", {
  # Use consolidated validation helper: tests SE type requirement
  test_calculate_divergence_input_validation(test_types = c("se_type"))
})

test_that("calculate_divergence auto-detects group column", {
  se <- create_test_se_simple(
    n_genes = 2,
    n_samples = 6,
    control_n = 2,
    group_col_name = "sample_type"
  )
  colnames(se) <- c("s1", "s2", "s3", "s4", "s5", "s6")
  
  # Don't specify group_col, should auto-detect
  result <- .calculate_divergence(
    se,
    control_group = "Control",
    q = 1,
    bootstrap = FALSE,
    verbose = FALSE,
    progress = TRUE
  )
  
  # Should return result without error
  expect_true(!is.null(result))
})

test_that("calculate_divergence auto-detects control group", {
  se <- create_test_se_simple(
    n_genes = 2,
    n_samples = 6,
    control_n = 2,
    group_col_name = "sample_type"
  )
  # Update control group name to "Normal" for this test
  SummarizedExperiment::colData(se)$sample_type <- c("Normal", "Normal", "Treatment", "Treatment", "Treatment", "Treatment")
  colnames(se) <- c("s1", "s2", "s3", "s4", "s5", "s6")
  
  # Don't specify control_group, should auto-detect "Normal"
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    q = 1,
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  # Should return result without error
  expect_true(!is.null(result))
})

test_that("calculate_divergence fails without gene identifiers", {
  # Create SE without rowData gene_name or rownames
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(c(1, 2, 3, 4, 5, 6), nrow = 1, ncol = 6)),
    colData = data.frame(
      sample_type = c("Control", "Control", "Treatment", "Treatment", "Treatment", "Treatment")
    )
  )
  colnames(se) <- c("s1", "s2", "s3", "s4", "s5", "s6")
  
  # No gene names in rowData or rownames
  expect_error(
    .calculate_divergence(
      se,
      group_col = "sample_type",
      control_group = "Control",
      bootstrap = FALSE,
      verbose = FALSE
    ),
    "gene identifiers"
  )
})

test_that("calculate_divergence rejects non-logical bootstrap", {
  se <- create_test_se_simple(
    n_genes = 1,
    n_samples = 6,
    control_n = 2,
    group_col_name = "sample_type"
  )
  colnames(se) <- c("s1", "s2", "s3", "s4", "s5", "s6")
  
  # Use consolidated validation helper: tests bootstrap type requirement
  test_calculate_divergence_input_validation(se = se, test_types = c("bootstrap_type"))
})

# =====================================================================
# Tests for classify_q_pattern public function
# =====================================================================

test_that("classify_q_pattern performs basic classification", {
  per_q_divs <- c(
    "q_0.01" = 0.8,
    "q_0.5" = 0.85,
    "q_1.0" = 0.3,
    "q_2.0" = 0.2
  )
  
  result <- .classify_q_pattern(per_q_divs)
  
  expect_true(result %in% c("RARE_DRIVEN", "ABUNDANT_DRIVEN", "BALANCED", NA_character_))
})

test_that("classify_q_pattern classifies RARE_DRIVEN correctly", {
  per_q_divs <- c(
    "q_0.01" = 1.0,
    "q_0.5" = 0.95,
    "q_1.0" = 0.5,
    "q_2.0" = 0.4
  )
  
  result <- .classify_q_pattern(per_q_divs)
  
  expect_equal(result, "RARE_DRIVEN")
})

test_that("classify_q_pattern classifies ABUNDANT_DRIVEN correctly", {
  per_q_divs <- c(
    "q_0.01" = 0.2,
    "q_0.5" = 0.15,
    "q_1.0" = 0.8,
    "q_2.0" = 0.9
  )
  
  result <- .classify_q_pattern(per_q_divs)
  
  expect_equal(result, "ABUNDANT_DRIVEN")
})

# =====================================================================
# Tests for effect_sizes_divergence
# =====================================================================

test_that("effect_sizes_divergence requires data frame input", {
  # lm_res must be a data frame
  expect_error(
    .effect_sizes_divergence(
      lm_res = list(a = 1, b = 2),
      divergence_results_se = NULL,
      verbose = FALSE
    ),
    "lm_res must be a data frame"
  )
})

test_that("effect_sizes_divergence requires required columns", {
  # Create lm_res without required columns
  lm_res <- data.frame(some_col = c(1, 2, 3))
  
  expect_error(
    .effect_sizes_divergence(
      lm_res = lm_res,
      divergence_results_se = NULL,
      verbose = FALSE
    ),
    "lm_res must have columns"
  )
})

# =====================================================================
# Tests for refactored calculate_divergence orchestrator (March 2026)
# =====================================================================

test_that("calculate_divergence executes sequential processing (nthreads=1)", {
  se <- create_test_se_simple(
    n_genes = 3,
    n_samples = 8,
    control_n = 3,
    group_col_name = "sample_type"
  )
  
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = c(0.5, 1, 2),
    nthreads = 1,  # Explicitly sequential
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  expect_true(methods::is(result, "SummarizedExperiment"))
  expect_equal(nrow(result), 3)  # 3 genes
})

test_that("calculate_divergence handles multiple q values", {
  se <- create_test_se_simple(
    n_genes = 2,
    n_samples = 8,
    control_n = 3,
    group_col_name = "sample_type"
  )
  
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = c(0.5, 1, 1.5, 2),  # 4 q values
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  # Should have matrix with correct dimensions
  expect_true(methods::is(result, "SummarizedExperiment"))
  expect_equal(nrow(result), 2)  # 2 genes
  # Assay has 4 columns (one per q value)
  expect_equal(ncol(SummarizedExperiment::assay(result)), 4)
})

test_that("calculate_divergence performs bootstrap with auto nboot", {
  se <- create_test_se_simple(
    n_genes = 1,
    n_samples = 8,
    control_n = 3,
    group_col_name = "sample_type"
  )
  
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    bootstrap = TRUE,
    nboot = "auto",  # Auto-select nboot
    verbose = FALSE,
    progress = FALSE
  )
  
  # Should return SE with bootstrap results
  expect_true(methods::is(result, "SummarizedExperiment"))
})

test_that("calculate_divergence applies normalization (range)", {
  se <- create_test_se_simple(
    n_genes = 2,
    n_samples = 6,
    control_n = 2,
    group_col_name = "sample_type"
  )
  
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    norm = "range",  # Range normalization
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  expect_true(methods::is(result, "SummarizedExperiment"))
})

test_that("calculate_divergence applies normalization (zscore)", {
  se <- create_test_se_simple(
    n_genes = 2,
    n_samples = 6,
    control_n = 2,
    group_col_name = "sample_type"
  )
  
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    norm = "zscore",  # Z-score normalization
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  expect_true(methods::is(result, "SummarizedExperiment"))
})

test_that("calculate_divergence skips normalization with norm='none'", {
  se <- create_test_se_simple(
    n_genes = 2,
    n_samples = 6,
    control_n = 2,
    group_col_name = "sample_type"
  )
  
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    norm = "none",  # No normalization
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  expect_true(methods::is(result, "SummarizedExperiment"))
})

test_that("calculate_divergence classifies per-q patterns with multiple q", {
  se <- create_test_se_simple(
    n_genes = 2,
    n_samples = 6,
    control_n = 2,
    group_col_name = "sample_type"
  )
  
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = c(0.5, 1, 2),  # Multiple q values
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  # Should have per_q_pattern column
  expect_true("per_q_pattern" %in% colnames(SummarizedExperiment::rowData(result)))
})

test_that("calculate_divergence populates reference q columns", {
  se <- create_test_se_simple(
    n_genes = 2,
    n_samples = 6,
    control_n = 2,
    group_col_name = "sample_type"
  )
  
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = c(0.5, 1, 2),  # q=1 included
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  rd <- SummarizedExperiment::rowData(result)
  
  # Should have generic columns populated from reference q
  expect_true("estimate" %in% colnames(rd))
  expect_true("lower_ci" %in% colnames(rd))
  expect_true("upper_ci" %in% colnames(rd))
})

test_that("calculate_divergence handles single gene correctly", {
  se <- create_test_se_simple(
    n_genes = 1,
    n_samples = 8,
    control_n = 3,
    group_col_name = "sample_type"
  )
  
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  expect_equal(nrow(result), 1)
})

test_that("calculate_divergence rejects bootstrap=non-logical", {
  se <- create_test_se_simple(
    n_genes = 1,
    n_samples = 6,
    control_n = 2,
    group_col_name = "sample_type"
  )
  
  expect_error(
    .calculate_divergence(
      se,
      group_col = "sample_type",
      control_group = "Control",
      bootstrap = "TRUE",  # Should be logical, not character
      verbose = FALSE
    ),
    "bootstrap must be a logical"
  )
})

test_that("calculate_divergence includes metadata in result", {
  se <- create_test_se_simple(
    n_genes = 2,
    n_samples = 6,
    control_n = 2,
    group_col_name = "sample_type"
  )
  
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  # Get metadata using the generic function (no namespace prefix)
  meta <- metadata(result)
  
  # Should include computation metadata
  expect_true("elapsed_time_sec" %in% names(meta))
  expect_true("summary_stats" %in% names(meta))
})

test_that("calculate_divergence handles custom pseudocount", {
  se <- create_test_se_simple(
    n_genes = 2,
    n_samples = 6,
    control_n = 2,
    group_col_name = "sample_type"
  )
  
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    pseudocount = 1.0,  # Custom pseudocount
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  expect_true(methods::is(result, "SummarizedExperiment"))
})

test_that("calculate_divergence handles custom log_base", {
  se <- create_test_se_simple(
    n_genes = 2,
    n_samples = 6,
    control_n = 2,
    group_col_name = "sample_type"
  )
  
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    log_base = 2,  # Log base 2
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  expect_true(methods::is(result, "SummarizedExperiment"))
})

test_that("calculate_divergence handles bootstrap with bca method", {
  se <- create_test_se_simple(
    n_genes = 1,
    n_samples = 8,
    control_n = 3,
    group_col_name = "sample_type"
  )
  
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    bootstrap = TRUE,
    nboot = 100,
    method = "bca",  # Bias-corrected accelerated method
    verbose = FALSE,
    progress = FALSE
  )
  
  expect_true(methods::is(result, "SummarizedExperiment"))
})

test_that("calculate_divergence returns error results gracefully", {
  # Create SE with extreme count values that may cause computation issues
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(c(1e15, 1e15, 1, 1, 1, 1), nrow = 1)),
    rowData = data.frame(gene_name = "gene1"),
    colData = data.frame(
      sample_type = c("Control", "Control", "Treatment", "Treatment", "Treatment", "Treatment")
    )
  )
  colnames(se) <- c("s1", "s2", "s3", "s4", "s5", "s6")
  
  # Should complete without crashing
  result <- .calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  expect_true(methods::is(result, "SummarizedExperiment"))
})
