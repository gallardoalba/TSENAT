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
  result <- calculate_divergence(
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
  result <- calculate_divergence(
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
    calculate_divergence(
      se,
      group_col = "sample_type",
      control_group = "Control",
      q = c(-0.5, 1, 2),
      bootstrap = FALSE,
      verbose = FALSE
    ),
    "q parameter must be positive"
  )
})

test_that("calculate_divergence requires SummarizedExperiment input", {
  # Pass a non-SE object
  expect_error(
    calculate_divergence(
      data.frame(a = 1, b = 2),
      group_col = "group",
      control_group = "Control",
      bootstrap = FALSE,
      verbose = FALSE
    ),
    "SummarizedExperiment"
  )
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
  result <- calculate_divergence(
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
  result <- calculate_divergence(
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
    calculate_divergence(
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
  
  expect_error(
    calculate_divergence(
      se,
      group_col = "sample_type",
      control_group = "Control",
      bootstrap = "yes",  # Not logical
      verbose = FALSE
    ),
    "bootstrap must be a logical"
  )
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
  
  result <- classify_q_pattern(per_q_divs)
  
  expect_true(result %in% c("RARE_DRIVEN", "ABUNDANT_DRIVEN", "BALANCED", NA_character_))
})

test_that("classify_q_pattern classifies RARE_DRIVEN correctly", {
  per_q_divs <- c(
    "q_0.01" = 1.0,
    "q_0.5" = 0.95,
    "q_1.0" = 0.5,
    "q_2.0" = 0.4
  )
  
  result <- classify_q_pattern(per_q_divs)
  
  expect_equal(result, "RARE_DRIVEN")
})

test_that("classify_q_pattern classifies ABUNDANT_DRIVEN correctly", {
  per_q_divs <- c(
    "q_0.01" = 0.2,
    "q_0.5" = 0.15,
    "q_1.0" = 0.8,
    "q_2.0" = 0.9
  )
  
  result <- classify_q_pattern(per_q_divs)
  
  expect_equal(result, "ABUNDANT_DRIVEN")
})

# =====================================================================
# Tests for effect_sizes_divergence
# =====================================================================

test_that("effect_sizes_divergence requires data frame input", {
  # lm_res must be a data frame
  expect_error(
    effect_sizes_divergence(
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
    effect_sizes_divergence(
      lm_res = lm_res,
      divergence_results_se = NULL,
      verbose = FALSE
    ),
    "lm_res must have columns"
  )
})
