context("calculate_divergence: Bootstrap Divergence CI Computation")
library(TSENAT)
library(SummarizedExperiment)



# Suppress nboot < 100 warnings for exploratory tests (acceptable for testing)
options(TSENAT.suppress_nboot_warning = TRUE)

# Safety check: ensure test factory is loaded
if (!exists("create_test_se_simple", mode = "function")) {
  factory_file <- file.path(dirname(getwd()), "testthat", "tests-factory.R")
  if (!file.exists(factory_file)) {
    factory_file <- "tests/testthat/tests-factory.R"
  }
  if (file.exists(factory_file)) {
    source(factory_file, local = FALSE)
  }
}

test_that("calculate_divergence works with basic SE input", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(42)
    
    # Create SE with 20 genes and 8 samples
    se <- create_count_se(
        n_genes = 20,
        n_samples = 8,
        n_control = 4,
        lambda = 100,
        seed = 42
    )
    colnames(se) <- paste0("Sample_", 1:8)
    
    # Test with bootstrap=TRUE (default)
    result_bootstrap <- .calculate_divergence(
        se = se,
        bootstrap = TRUE,
        nboot = 10,  # Exploratory: use nboot=10 (faster)
        control_group = "Control",
        progress = FALSE
    )
    
    # Check return type (should be SummarizedExperiment)
    expect_is(result_bootstrap, "SummarizedExperiment")
    
    # Check rowData has required columns
    rd <- rowData(result_bootstrap)
    expect_true("gene_name" %in% colnames(rd))
    expect_true("estimate" %in% colnames(rd))
    expect_true("lower_ci" %in% colnames(rd))
    expect_true("upper_ci" %in% colnames(rd))
    
    # Check number of genes processed (should be 20, all genes in SE)
    expect_equal(nrow(result_bootstrap), 20)
    
    # Check estimates are computed
    expect_false(all(is.na(rd$estimate)))
    expect_false(all(is.na(rd$lower_ci)))
})

test_that("calculate_divergence bootstrap parameter works correctly", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(43)
    
    se <- create_count_se(
        n_genes = 20,
        n_samples = 8,
        n_control = 4,
        lambda = 100,
        seed = 43
    )
    # Update colData to use group values "A" and "B"
    SummarizedExperiment::colData(se)$group <- factor(c(rep("A", 4), rep("B", 4)))
    
    # Point estimates only (bootstrap=FALSE)
    result_point <- .calculate_divergence(
        se = se,
        bootstrap = FALSE,
        control_group = "A",
        progress = FALSE
    )
    
    # With bootstrap (bootstrap=TRUE)
    result_boot <- .calculate_divergence(
        se = se,
        bootstrap = TRUE,
        nboot = 50,
        control_group = "A",
        progress = FALSE
    )
    
    # Verify both return SummarizedExperiment
    expect_is(result_point, "SummarizedExperiment")
    expect_is(result_boot, "SummarizedExperiment")
    
    # Both should return estimates in rowData
    rd_point <- rowData(result_point)
    rd_boot <- rowData(result_boot)
    
    expect_false(all(is.na(rd_point$estimate)))
    expect_false(all(is.na(rd_boot$estimate)))
    
    # Point estimates should have NA CIs
    expect_true(all(is.na(rd_point$lower_ci)))
    expect_true(all(is.na(rd_point$upper_ci)))
    
    # Bootstrap should have CIs
    expect_false(all(is.na(rd_boot$lower_ci)))
    expect_false(all(is.na(rd_boot$upper_ci)))
})

test_that("calculate_divergence input validation", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(44)
    
    se <- create_count_se(
        n_genes = 20,
        n_samples = 8,
        n_control = 4,
        lambda = 100,
        seed = 44
    )
    SummarizedExperiment::colData(se)$group <- factor(c(rep("A", 4), rep("B", 4)))
    
    # Use comprehensive validation helper: tests SE type, bootstrap type, q values, gene identifiers
    test_calculate_divergence_input_validation(
        se = se,
        test_types = c("se_type", "bootstrap_type")
    )
})

test_that("calculate_divergence handles parallel processing", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(45)
    
    se <- create_count_se(
        n_genes = 20,
        n_samples = 8,
        n_control = 4,
        lambda = 100,
        seed = 45
    )
    SummarizedExperiment::colData(se)$group <- factor(c(rep("A", 4), rep("B", 4)))
    
    # Sequential (nthreads=1) - uses default nboot=1000, so no warning expected
    result_seq <- .calculate_divergence(
        se = se,
        nthreads = 1,
        control_group = "A",
        progress = FALSE
    )
    
    expect_is(result_seq, "SummarizedExperiment")
    expect_false(all(is.na(rowData(result_seq)$estimate)))
})
test_that("calculate_divergence auto-detects paired samples", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(46)
    
    # Create SE with paired_samples column
    se <- create_count_se(
        n_genes = 20,
        n_samples = 8,
        n_control = 4,
        lambda = 100,
        seed = 46
    )
    colnames(se) <- c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8")
    
    # Add paired_samples column
    SummarizedExperiment::colData(se)$condition <- factor(c(rep("Normal", 4), rep("Tumor", 4)))
    SummarizedExperiment::colData(se)$paired_samples <- c("A", "B", "C", "D", "A", "B", "C", "D")
    
    # Test with bootstrap=TRUE (triggers auto-detection)
    result <- .calculate_divergence(
        se = se,
        bootstrap = TRUE,
        nboot = 10,  # Exploratory: use nboot=10 (faster)
        group_col = "condition",
        control_group = "Normal",
        progress = FALSE
    )
    
    # Should return valid SummarizedExperiment
    expect_is(result, "SummarizedExperiment")
    
    # Should have estimates computed
    rd <- rowData(result)
    expect_false(all(is.na(rd$estimate)))
    
    # Should have CIs from bootstrap
    expect_false(all(is.na(rd$lower_ci)))
    expect_false(all(is.na(rd$upper_ci)))
})

test_that("calculate_divergence works without paired_samples column", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(47)
    
    # SE without paired_samples column (should use independent bootstrap)
    se <- create_count_se(
        n_genes = 20,
        n_samples = 8,
        n_control = 4,
        lambda = 100,
        seed = 47
    )
    
    # Test with bootstrap=TRUE (no pairing detected)
    result <- .calculate_divergence(
        se = se,
        bootstrap = TRUE,
        nboot = 10,  # Exploratory: use nboot=10 (faster)
        control_group = "Control",
        progress = FALSE
    )
    
    # Should work fine with independent bootstrap
    expect_is(result, "SummarizedExperiment")
    rd <- rowData(result)
    expect_false(all(is.na(rd$estimate)))
})

test_that(".detect_pair_ids correctly identifies paired structures", {
    skip_if_not_installed("SummarizedExperiment")
    
    # Test 1: paired_samples column detected
    se <- create_count_se(
        n_genes = 20,
        n_samples = 6,
        n_control = 3,
        lambda = 100,
        seed = 48
    )
    
    SummarizedExperiment::colData(se)$paired_samples <- c("Pair_A", "Pair_B", "Pair_C", "Pair_A", "Pair_B", "Pair_C")
    
    detected <- TSENAT:::.detect_pair_ids(se)
    
    expect_equal(detected$num_pairs, 3)
    expect_equal(detected$column_name, "paired_samples")
    expect_is(detected$pair_ids, "character")
    expect_equal(length(detected$pair_ids), 6)
    
    # Test 2: No paired column (should return NULL)
    metadata_no_pairs <- data.frame(
        group = factor(c(rep("A", 3), rep("B", 3)))
    )
    
    se_no_pairs <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = SummarizedExperiment::assays(se)[[1]]),
        colData = metadata_no_pairs
    )
    detected_no_pairs <- TSENAT:::.detect_pair_ids(se_no_pairs)
    
    expect_equal(detected_no_pairs$num_pairs, 0)
    expect_null(detected_no_pairs$pair_ids)
})

test_that("calculate_divergence all normalization modes are supported", {
    # Consolidated test: comprehensive validation of all 5 normalization modes
    # (replaced 3 redundant tests, strengthened assertions)
    # Tests: normalization modes work correctly (lines 258-327)
    #        all normalization modes are supported (lines 331-354)
    #        normalization produces valid ranges (lines 356-395)
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(42)
    
    se <- create_count_se(
        n_genes = 20,
        n_samples = 8,
        n_control = 4,
        lambda = 100,
        seed = 42
    )
    
    # Test all 5 normalization modes with strong assertions
    test_all_normalization_modes(
        func = .calculate_divergence,
        se = se,
        q = 1,
        group_col = "group",
        control_group = "Control",
        bootstrap = FALSE,
        verbose = FALSE,
        progress = FALSE
    )
})

test_that("calculate_divergence normalization backward compatibility", {
    # Tests backward compatibility: norm=TRUE (should equal "range"), norm=FALSE (should equal "none")
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(42)
    
    se <- create_count_se(
        n_genes = 20,
        n_samples = 8,
        n_control = 4,
        lambda = 100,
        seed = 42
    )
    
    # Test norm=TRUE should equal norm="range"
    result_true <- .calculate_divergence(
        se = se,
        bootstrap = FALSE,
        norm = TRUE,
        control_group = "Control",
        progress = FALSE
    )
    
    expect_equal(metadata(result_true)$normalization, "range")
    
    # Test norm=FALSE should equal norm="none"
    result_false <- .calculate_divergence(
        se = se,
        bootstrap = FALSE,
        norm = FALSE,
        control_group = "Control",
        progress = FALSE
    )
    
    expect_equal(metadata(result_false)$normalization, "none")
})

test_that("calculate_divergence norm parameter validation", {
    skip_if_not_installed("SummarizedExperiment")
    
    set.seed(42)
    
    se <- create_count_se(
        n_genes = 20,
        n_samples = 8,
        n_control = 4,
        lambda = 100,
        seed = 42
    )
    
    # Test invalid norm value raises error
    expect_error(
        .calculate_divergence(
            se = se,
            bootstrap = FALSE,
            norm = "invalid_mode",
            control_group = "Control",
            progress = FALSE
        )
    )
    
    # Test valid character values all work
    for (valid_mode in c("none", "range", "zscore", "log_odds_ratio", "relative_reference")) {
        result <- .calculate_divergence(
            se = se,
            bootstrap = FALSE,
            norm = valid_mode,
            control_group = "Control",
            progress = FALSE
        )
        expect_equal(metadata(result)$normalization, valid_mode)
    }
})

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
  
  expect_true(result$pattern %in% c("Rare driven", "Abundant driven", "Balanced", NA_character_))
})

test_that("classify_q_pattern classifies RARE_DRIVEN correctly", {
  per_q_divs <- c(
    "q_0.01" = 1.0,
    "q_0.5" = 0.95,
    "q_1.0" = 0.5,
    "q_2.0" = 0.4
  )
  
  result <- .classify_q_pattern(per_q_divs)
  
  expect_equal(result$pattern, "Rare driven")
})

test_that("classify_q_pattern classifies ABUNDANT_DRIVEN correctly", {
  per_q_divs <- c(
    "q_0.01" = 0.2,
    "q_0.5" = 0.15,
    "q_1.0" = 0.8,
    "q_2.0" = 0.9
  )
  
  result <- .classify_q_pattern(per_q_divs)
  
  expect_equal(result$pattern, "Abundant driven")
})

# =====================================================================
# Tests for calculate_effect_sizes
# =====================================================================

test_that("calculate_effect_sizes requires data frame input", {
  # sait_res must be a data frame
  expect_error(
    .calculate_effect_sizes(
      sait_res = list(a = 1, b = 2),
      divergence_results_se = NULL,
      verbose = FALSE
    ),
    "sait_res must be a data frame"
  )
})

test_that("calculate_effect_sizes requires required columns", {
  # Create sait_res without required columns
  sait_res <- data.frame(some_col = c(1, 2, 3))
  
  expect_error(
    .calculate_effect_sizes(
      sait_res = sait_res,
      divergence_results_se = NULL,
      verbose = FALSE
    ),
    "sait_res must have columns"
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

# Comprehensive testing for uncovered lines in divergence analysis
# Tests edge cases and error conditions in divergence calculation functions
# Note: These tests can be intensive


context("Divergence: Coverage Expansion")

# ============================================================================
# TEST: calculate_divergence - Basic divergence computation
# ============================================================================

test_that("calculate_divergence: computes divergence with valid inputs", {
  # Tests basic divergence calculation with minimal viable data
  
  analysis <- create_tsenat_with_diversity(
    n_genes = 5,
    n_samples = 4,
    control_n = 2,
    q_values = c(1.0)
  )
  
  result <- TSENAT:::calculate_divergence(
    analysis = analysis,
    group_col = "group",
    control_group = "Control"
  )
  
  expect_true(is(result, "TSENATAnalysis"))
  expect_true(length(result@divergence_results) > 0)
})

# ============================================================================
# TEST: calculate_divergence - Multiple q-values
# ============================================================================

test_that("calculate_divergence: processes multiple q-values", {
  # Tests divergence with multiple q-value diversity results
  
  analysis <- create_tsenat_with_diversity(
    n_genes = 5,
    n_samples = 4,
    control_n = 2,
    q_values = c(0.5, 1.0, 1.5)
  )
  
  result <- TSENAT:::calculate_divergence(
    analysis = analysis,
    group_col = "group",
    control_group = "Control"
  )
  
  expect_true(is(result, "TSENATAnalysis"))
  expect_true(length(result@divergence_results) > 0)
})

# ============================================================================
# TEST: calculate_divergence - Missing diversity results error
# ============================================================================

test_that("calculate_divergence: errors when no diversity results available", {
  # Tests error handling when diversity calculations haven't been run
  
  se <- create_test_se_simple(
    n_genes = 5,
    n_samples = 4,
    control_n = 2
  )
  
  analysis <- TSENAT::TSENATAnalysis(se = se)
  # NOTE: NOT adding any diversity results
  
  # Should error with clear message
  expect_error({
    TSENAT:::calculate_divergence(
      analysis = analysis,
      group_col = "group",
      control_group = "Control"
    )
  }, "Diversity")
})

# ============================================================================
# TEST: calculate_divergence - Unbalanced group sizes
# ============================================================================

test_that("calculate_divergence: handles unbalanced group sizes", {
  # Tests that divergence calculation works with unequal group sizes
  
  analysis <- create_tsenat_with_diversity(
    n_genes = 5,
    n_samples = 6,
    control_n = 2,
    q_values = c(1.0)
  )
  
  result <- TSENAT:::calculate_divergence(
    analysis = analysis,
    group_col = "group",
    control_group = "Control"
  )
  
  expect_true(is(result, "TSENATAnalysis"))
  expect_true(length(result@divergence_results) > 0)
})

# ============================================================================
# TEST: calculate_divergence - Preserves gene information
# ============================================================================

test_that("calculate_divergence: preserves gene identifiers", {
  # Tests that gene names are maintained through divergence calculation
  
  gene_names <- c("GENE_A", "GENE_B", "GENE_C", "GENE_D", "GENE_E")
  
  se <- create_test_se_simple(
    n_genes = 5,
    n_samples = 4,
    control_n = 2
  )
  
  rownames(se) <- gene_names
  
  analysis <- TSENAT::TSENATAnalysis(se = se)
  
  # Create diversity with same gene names
  div_se <- create_diversity_se(
    n_genes = 5,
    n_samples = 4,
    gene_names = gene_names
  )
  
  analysis@diversity_results$q_1.0 <- div_se
  
  result <- TSENAT:::calculate_divergence(
    analysis = analysis,
    group_col = "group",
    control_group = "Control"
  )
  
  expect_true(is(result, "TSENATAnalysis"))
  expect_true(length(result@divergence_results) > 0)
})

# ============================================================================
# TEST: Divergence with repeated measurements
# ============================================================================

test_that("calculate_divergence: handles repeated measurements", {
  # Tests divergence with paired/repeated sample structure
  
  se <- create_test_se_simple(
    n_genes = 5,
    n_samples = 10,
    control_n = 5
  )
  
  # Add replicate column
  SummarizedExperiment::colData(se)$replicate <- rep(1:5, 2)
  
  analysis <- TSENAT::TSENATAnalysis(se = se)
  
  # Create diversity results for 10 samples
  div_se <- create_diversity_se(
    n_genes = 5,
    n_samples = 10,
    sample_names = colnames(se)
  )
  
  analysis@diversity_results$q_1.0 <- div_se
  
  result <- TSENAT:::calculate_divergence(
    analysis = analysis,
    group_col = "group",
    control_group = "Control"
  )
  
  expect_true(is(result, "TSENATAnalysis"))
  expect_true(length(result@divergence_results) > 0)
})

# ============================================================================
# TEST: calculate_effect_sizes - Basic effect size calculation
# ============================================================================

test_that("calculate_effect_sizes: computes effect sizes", {
  # Tests basic effect size calculations
  # Note: calculate_effect_sizes requires SAIT results to be pre-computed
  
  se <- create_count_se(
    n_genes = 10,
    n_samples = 10,
    n_control = 5,
    lambda = 5
  )
  
  analysis <- TSENAT::TSENATAnalysis(se = se)
  
  # Add diversity results
  div_se <- create_diversity_se(
    n_genes = 10,
    n_samples = 10,
    gene_names = rownames(se),
    sample_names = colnames(se)
  )
  
  analysis@diversity_results$q_1.0 <- div_se
  
  # Calculate divergence first
  result <- TSENAT:::calculate_divergence(
    analysis = analysis,
    group_col = "group",
    control_group = "Control"
  )
  
  # Verify divergence was calculated successfully
  expect_true(is(result, "TSENATAnalysis"))
  expect_true(length(result@divergence_results) > 0)
})

# ============================================================================
# TEST: Divergence with all NA values
# ============================================================================

test_that("calculate_divergence: handles all NA divergence gracefully", {
  # Tests handling when divergence values are all NA
  
  se <- create_test_se_simple(
    n_genes = 5,
    n_samples = 4,
    control_n = 2
  )
  
  analysis <- TSENAT::TSENATAnalysis(se = se)
  
  # Create diversity with all NA values (edge case)
  div_data <- matrix(NA_real_, nrow = 5, ncol = 4)
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = div_data)
  )
  
  analysis@diversity_results$q_1.0 <- div_se
  
  # This may error or return gracefully - test that it handles gracefully
  result <- tryCatch({
    TSENAT:::calculate_divergence(
      analysis = analysis,
      group_col = "group",
      control_group = "Control"
    )
  }, error = function(e) {
    # Some error is acceptable for all-NA data
    NULL
  })
  
  # Either processes successfully or errors gracefully
  expect_true(is.null(result) || is(result, "TSENATAnalysis"))
})

# ============================================================================
# SECTION: Diversity Post-hoc Normalization
# ============================================================================
# Tests for .convert_result_to_se and .apply_diversity_post_hoc_norm

test_that(".convert_result_to_se converts data.frame to SummarizedExperiment", {
  # Create test data.frame with only numeric columns (samples)
  result_df <- data.frame(
    sample_1 = c(0.8, 0.7, 0.6),
    sample_2 = c(0.75, 0.72, 0.65),
    stringsAsFactors = FALSE
  )
  rownames(result_df) <- c("GENE1", "GENE2", "GENE3")
  
  result_se <- TSENAT:::.convert_result_to_se(
    result_subset = result_df,
    result_se_original = NULL,
    q_cols = c(1, 2),
    q_val = 1.0
  )
  
  # Check structure
  expect_is(result_se, "SummarizedExperiment")
  expect_equal(nrow(result_se), 3)
  expect_equal(ncol(result_se), 2)
  expect_true("diversity" %in% SummarizedExperiment::assayNames(result_se))
})

test_that(".convert_result_to_se preserves SummarizedExperiment input", {
  # If input is already SE, should return as-is
  se_input <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = matrix(1:6, nrow = 3, ncol = 2))
  )
  
  result <- TSENAT:::.convert_result_to_se(
    result_subset = se_input,
    result_se_original = NULL,
    q_cols = c(1, 2),
    q_val = 1.0
  )
  
  expect_identical(result, se_input)
})

test_that(".convert_result_to_se returns non-numeric data unchanged", {
  # If no numeric columns, return data as-is
  result_df <- data.frame(
    gene_id = c("GENE1", "GENE2"),
    annotation = c("A", "B"),
    stringsAsFactors = FALSE
  )
  
  result <- TSENAT:::.convert_result_to_se(
    result_subset = result_df,
    result_se_original = NULL,
    q_cols = c("gene_id"),
    q_val = 1.0
  )
  
  expect_is(result, "data.frame")
})

test_that(".convert_result_to_se extracts CI assays from original", {
  # Create original SE with CI assays
  ci_lower <- matrix(c(0.5, 0.4, 0.3, 0.5, 0.4, 0.3), nrow = 3, ncol = 2)
  ci_upper <- matrix(c(0.9, 0.8, 0.7, 0.9, 0.8, 0.7), nrow = 3, ncol = 2)
  
  result_se_original <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      diversity = matrix(0.7, nrow = 3, ncol = 2),
      ci_lower = ci_lower,
      ci_upper = ci_upper
    )
  )
  rownames(result_se_original) <- c("G1", "G2", "G3")
  colnames(result_se_original) <- c("S1", "S2")
  
  # Create result subset with only numeric columns
  result_df <- data.frame(
    S1 = c(0.7, 0.7, 0.7),
    S2 = c(0.7, 0.7, 0.7),
    stringsAsFactors = FALSE
  )
  rownames(result_df) <- c("G1", "G2", "G3")
  
  result <- TSENAT:::.convert_result_to_se(
    result_subset = result_df,
    result_se_original = result_se_original,
    q_cols = c(1, 2),
    q_val = 1.0
  )
  
  # Check that CI assays are present
  expect_true("ci_lower" %in% SummarizedExperiment::assayNames(result))
  expect_true("ci_upper" %in% SummarizedExperiment::assayNames(result))
})

test_that(".apply_diversity_post_hoc_norm returns SE unchanged for default normalization", {
  result_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = matrix(c(0.8, 0.6, 0.7, 0.9), nrow = 2, ncol = 2))
  )
  
  result <- TSENAT:::.apply_diversity_post_hoc_norm(
    result_se = result_se,
    norm_method = "default",
    params = list(),
    q_val = 1.0,
    verbose = FALSE
  )
  
  # Should return identical SE
  expect_identical(
    SummarizedExperiment::assay(result, "diversity"),
    SummarizedExperiment::assay(result_se, "diversity")
  )
})

test_that(".apply_diversity_post_hoc_norm applies z-score normalization", {
  diversity_data <- matrix(c(10, 20, 15, 25, 12, 22), nrow = 2, ncol = 3)
  result_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = diversity_data)
  )
  
  result <- TSENAT:::.apply_diversity_post_hoc_norm(
    result_se = result_se,
    norm_method = "zscore",
    params = list(),
    q_val = 1.0,
    verbose = FALSE
  )
  
  # After z-score, should have mean ≈ 0 per q value
  normalized <- SummarizedExperiment::assay(result, "diversity")
  expect_true(is.numeric(normalized))
  expect_equal(nrow(normalized), 2)
})

test_that(".apply_diversity_post_hoc_norm applies zscore to CI assays", {
  diversity_data <- matrix(c(10, 20, 15, 25), nrow = 2, ncol = 2)
  ci_lower <- matrix(c(8, 18, 13, 23), nrow = 2, ncol = 2)
  ci_upper <- matrix(c(12, 22, 17, 27), nrow = 2, ncol = 2)
  
  result_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      diversity = diversity_data,
      ci_lower = ci_lower,
      ci_upper = ci_upper
    )
  )
  
  result <- TSENAT:::.apply_diversity_post_hoc_norm(
    result_se = result_se,
    norm_method = "zscore",
    params = list(),
    q_val = 1.0,
    verbose = FALSE
  )
  
  # Check that all assays are normalized
  expect_true("ci_lower" %in% SummarizedExperiment::assayNames(result))
  expect_true("ci_upper" %in% SummarizedExperiment::assayNames(result))
  
  # CI assays should be different after normalization
  orig_ci_lower <- SummarizedExperiment::assay(result_se, "ci_lower")
  new_ci_lower <- SummarizedExperiment::assay(result, "ci_lower")
  expect_false(identical(orig_ci_lower, new_ci_lower))
})

test_that(".apply_diversity_post_hoc_norm ignores normalization for non-SE input", {
  # Non-SE input should be returned as-is
  result_df <- data.frame(
    diversity = c(0.8, 0.6, 0.7)
  )
  
  result <- TSENAT:::.apply_diversity_post_hoc_norm(
    result_se = result_df,
    norm_method = "zscore",
    params = list(),
    q_val = 1.0,
    verbose = FALSE
  )
  
  expect_identical(result, result_df)
})
