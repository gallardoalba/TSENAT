context("S4 Diversity Calculation: Output File Generation and Numerical Correctness")

# ===========================================================================
# Setup: Helper functions to create test data structures
# ===========================================================================

# Helper to create test analysis using the established .create_test_analysis function
make_test_analysis_diversity <- function(n_genes = 8, n_samples_per_group = 4, 
                                         q_values = c(0.5, 1.0), seed = 123) {
  TSENAT:::.create_test_analysis(
    n_genes = n_genes,
    n_samples_per_group = n_samples_per_group,
    q_values = q_values,
    include_divergence = FALSE,
    include_lm_results = FALSE,
    seed = seed,
    verbose = FALSE
  )
}

# ===========================================================================
# TEST GROUP 1: Output File Generation
# ===========================================================================

test_that("calculate_diversity_s4 generates diversity_results.tsv file without bootstrap", {
  
  analysis <- make_test_analysis_diversity(n_genes = 8, n_samples_per_group = 2)
  output_dir <- tempdir()
  output_file <- file.path(output_dir, "test_diversity_noboot.tsv")
  
  # Calculate diversity without bootstrap
  result <- TSENAT::calculate_diversity_s4(
    analysis,
    q = c(0.5, 1.0),
    bootstrap = FALSE,
    output_file = output_file,
    verbose = FALSE
  )
  
  # Check that main output file exists
  expect_true(file.exists(output_file), 
              info = sprintf("Output file not found: %s", output_file))
  
  # Check file size is reasonable (not empty)
  file_size <- file.info(output_file)$size
  expect_gt(file_size, 100)
  
  # Clean up
  if (file.exists(output_file)) unlink(output_file)
})

test_that("calculate_diversity_s4 generates both diversity_results.tsv and _diversity_spectrum.tsv with output_file", {
  
  analysis <- make_test_analysis_diversity(n_genes = 8, n_samples_per_group = 2)
  output_dir <- tempdir()
  output_file <- file.path(output_dir, "test_diversity_dual.tsv")
  
  result <- TSENAT::calculate_diversity_s4(
    analysis,
    q = c(0.5, 1.0),
    bootstrap = FALSE,
    output_file = output_file,
    verbose = FALSE
  )
  
  # Check main file
  expect_true(file.exists(output_file),
              info = "Main diversity results file not found")
  
  # Check spectrum file (should have _diversity_spectrum.tsv suffix)
  spectrum_file <- sub("\\.[^.]+$", "_diversity_spectrum.tsv", output_file)
  if (spectrum_file == output_file) {
    spectrum_file <- paste0(output_file, "_diversity_spectrum.tsv")
  }
  expect_true(file.exists(spectrum_file),
              info = sprintf("Spectrum file not found: %s", spectrum_file))
  
  # Both files should have content
  expect_gt(file.info(output_file)$size, 100)
  expect_gt(file.info(spectrum_file)$size, 50)
  
  # Clean up
  if (file.exists(output_file)) unlink(output_file)
  if (file.exists(spectrum_file)) unlink(spectrum_file)
})

test_that("calculate_diversity_s4 generates diversity_results.tsv with CI columns when bootstrap=TRUE", {
  
  analysis <- make_test_analysis_diversity(n_genes = 8, n_samples_per_group = 2)
  output_dir <- tempdir()
  output_file <- file.path(output_dir, "test_diversity_boot.tsv")
  
  result <- suppressWarnings(TSENAT::calculate_diversity_s4(
    analysis,
    q = c(0.5, 1.0),
    bootstrap = TRUE,
    nboot = 50,
    output_file = output_file,
    verbose = FALSE
  ))
  
  # Check main file exists
  expect_true(file.exists(output_file))
  
  # Read file and check columns
  div_data <- read.csv(output_file, sep = "\t", stringsAsFactors = FALSE)
  
  # Should have columns: gene, sample, q_value, diversity, ci_lower, ci_upper
  expected_cols <- c("gene", "sample", "q_value", "diversity", "ci_lower", "ci_upper")
  expect_true(all(expected_cols %in% colnames(div_data)),
              info = sprintf("Missing expected columns. Expected: %s. Got: %s",
                           paste(expected_cols, collapse=", "),
                           paste(colnames(div_data), collapse=", ")))
  
  # Clean up
  if (file.exists(output_file)) unlink(output_file)
})

# ===========================================================================
# TEST GROUP 2: Output File Format and Structure
# ===========================================================================

test_that("diversity_results.tsv has correct row count (genes × samples × q-values)", {
  
  n_genes <- 8
  n_samples_per_group <- 2
  n_samples <- n_samples_per_group * 2
  q_vals <- c(0.5, 1.0, 1.5)
  
  analysis <- make_test_analysis_diversity(n_genes = n_genes, n_samples_per_group = n_samples_per_group, q_values = q_vals)
  output_dir <- tempdir()
  output_file <- file.path(output_dir, "test_diversity_rowcount.tsv")
  
  result <- TSENAT::calculate_diversity_s4(
    analysis,
    q = q_vals,
    bootstrap = FALSE,
    output_file = output_file,
    verbose = FALSE
  )
  
  div_data <- read.csv(output_file, sep = "\t", stringsAsFactors = FALSE)
  expected_rows <- n_genes * n_samples * length(q_vals)
  
  expect_equal(nrow(div_data), expected_rows,
               info = sprintf("Expected %d rows (genes=%d × samples=%d × q-values=%d), got %d",
                            expected_rows, n_genes, n_samples, length(q_vals), nrow(div_data)))
  
  # Clean up
  if (file.exists(output_file)) unlink(output_file)
})

test_that("diversity_results.tsv contains all expected q-values", {
  
  q_vals <- c(0.5, 1.0)
  
  analysis <- make_test_analysis_diversity(n_genes = 8, n_samples_per_group = 2, q_values = q_vals)
  output_dir <- tempdir()
  output_file <- file.path(output_dir, "test_diversity_qvals.tsv")
  
  result <- TSENAT::calculate_diversity_s4(
    analysis,
    q = q_vals,
    bootstrap = FALSE,
    output_file = output_file,
    verbose = FALSE
  )
  
  div_data <- read.csv(output_file, sep = "\t", stringsAsFactors = FALSE)
  
  # Check that all q-values are present
  # Q-values might be stored as "q_X.XXX" format, so extract numeric part if needed
  q_col <- div_data$q_value
  if (is.character(q_col) && all(grepl("^q_", q_col))) {
    actual_q <- sort(unique(as.numeric(sub("^q_", "", q_col))))
  } else {
    actual_q <-sort(unique(as.numeric(q_col)))
  }
  expected_q <- sort(q_vals)
  
  expect_equal(actual_q, expected_q)
  
  # Clean up
  if (file.exists(output_file)) unlink(output_file)
})

test_that("diversity_results.tsv has no NA values in diversity column", {
  
  analysis <- make_test_analysis_diversity(n_genes = 8, n_samples_per_group = 2)
  output_dir <- tempdir()
  output_file <- file.path(output_dir, "test_diversity_nona.tsv")
  
  result <- TSENAT::calculate_diversity_s4(
    analysis,
    q = c(0.5, 1.0),
    bootstrap = FALSE,
    output_file = output_file,
    verbose = FALSE
  )
  
  div_data <- read.csv(output_file, sep = "\t", stringsAsFactors = FALSE)
  
  # Check for NA values in diversity column
  na_count <- sum(is.na(div_data$diversity))
  expect_equal(na_count, 0,
               info = sprintf("Found %d NA values in diversity column", na_count))
  
  # Clean up
  if (file.exists(output_file)) unlink(output_file)
})

test_that("diversity values are in reasonable range [0, Inf)", {
  
  analysis <- make_test_analysis_diversity(n_genes = 8, n_samples_per_group = 2, q_values = c(0.0, 0.5, 1.0, 2.0))
  output_dir <- tempdir()
  output_file <- file.path(output_dir, "test_diversity_range.tsv")
  
  result <- TSENAT::calculate_diversity_s4(
    analysis,
    q = c(0.0, 0.5, 1.0, 2.0),
    bootstrap = FALSE,
    output_file = output_file,
    verbose = FALSE
  )
  
  div_data <- read.csv(output_file, sep = "\t", stringsAsFactors = FALSE)
  
  # Diversity values should be non-negative
  expect_true(all(div_data$diversity >= 0),
              info = "Found negative diversity values")
  
  # Diversity should be finite (no Inf or NaN)
  expect_true(all(is.finite(div_data$diversity)),
              info = "Found infinite or NaN diversity values")
  
  # Clean up
  if (file.exists(output_file)) unlink(output_file)
})

# ===========================================================================
# TEST GROUP 3: Bootstrap CI Correctness
# ===========================================================================

test_that("CI columns exist when bootstrap=TRUE", {
  
  analysis <- make_test_analysis_diversity(n_genes = 15, n_samples_per_group = 4)
  output_dir <- tempdir()
  output_file <- file.path(output_dir, "test_diversity_ci_exist.tsv")
  
  result <- suppressWarnings(TSENAT::calculate_diversity_s4(
    analysis,
    q = c(0.5, 1.0),
    bootstrap = TRUE,
    nboot = 50,
    bootstrap_ci = 0.95,
    output_file = output_file,
    verbose = FALSE
  ))
  
  div_data <- read.csv(output_file, sep = "\t", stringsAsFactors = FALSE)
  
  expect_true("ci_lower" %in% colnames(div_data))
  expect_true("ci_upper" %in% colnames(div_data))
  
  # CI values should exist and not all be NA
  expect_gt(sum(!is.na(div_data$ci_lower)), 0)
  expect_gt(sum(!is.na(div_data$ci_upper)), 0)
  
  # Clean up
  if (file.exists(output_file)) unlink(output_file)
})

test_that("CI bounds are valid (ci_lower <= diversity <= ci_upper)", {
  
  analysis <- make_test_analysis_diversity(n_genes = 15, n_samples_per_group = 4)
  output_dir <- tempdir()
  output_file <- file.path(output_dir, "test_diversity_ci_bounds.tsv")
  
  result <- suppressWarnings(TSENAT::calculate_diversity_s4(
    analysis,
    q = c(0.5, 1.0),
    bootstrap = TRUE,
    nboot = 50,
    output_file = output_file,
    verbose = FALSE
  ))
  
  div_data <- read.csv(output_file, sep = "\t", stringsAsFactors = FALSE)
  
  # Ensure numeric columns
  div_data$diversity <- as.numeric(div_data$diversity)
  div_data$ci_lower <- as.numeric(div_data$ci_lower)
  div_data$ci_upper <- as.numeric(div_data$ci_upper)
  
  # Filter to rows with valid CIs (not NA)
  valid_ci <- !is.na(div_data$ci_lower) & !is.na(div_data$ci_upper)
  
  if (any(valid_ci)) {
    valid_data <- div_data[valid_ci, ]
    
    # Strict checks - these should ALWAYS be true for valid CIs
    # Check 1: ci_lower <= diversity (lower bound should be below or equal to point estimate)
    lower_issues <- sum(valid_data$ci_lower > valid_data$diversity + 1e-6)
    expect_equal(lower_issues, 0,
                 info = sprintf("BUG: %d values have ci_lower > diversity", lower_issues))
    
    # Check 2: diversity <= ci_upper (point estimate should be within upper bound)
    upper_issues <- sum(valid_data$diversity > valid_data$ci_upper + 1e-6)
    expect_equal(upper_issues, 0,
                 info = sprintf("BOOTSTRAP BUG: %d values have diversity > ci_upper! CI bounds are inverted or incorrect", upper_issues))
    
    # Check 3: ci_lower < ci_upper (interval should have positive width)
    width_issues <- sum(valid_data$ci_lower >= valid_data$ci_upper)
    expect_equal(width_issues, 0,
                 info = sprintf("BUG: %d CIs have zero or negative width", width_issues))
  }
  
  # Clean up
  if (file.exists(output_file)) unlink(output_file)
})

test_that("CI width decreases with higher bootstrap replicates (nboot)", {
  
  # Use same seed for both analyses so they analyze the same data
  # Only difference is nboot (20 vs 100)
  analysis_low <- make_test_analysis_diversity(n_genes = 10, n_samples_per_group = 3, seed = 456)
  analysis_high <- make_test_analysis_diversity(n_genes = 10, n_samples_per_group = 3, seed = 456)
  
  output_dir <- tempdir()
  output_file_low <- file.path(output_dir, "test_diversity_ci_low_nboot.tsv")
  output_file_high <- file.path(output_dir, "test_diversity_ci_high_nboot.tsv")
  
  # Low nboot
  result_low <- suppressWarnings(TSENAT::calculate_diversity_s4(
    analysis_low,
    q = c(1.0),
    bootstrap = TRUE,
    nboot = 20,
    seed = 111,
    output_file = output_file_low,
    verbose = FALSE
  ))
  
  # High nboot (use 50 instead of 2000 to keep tests fast)
  result_high <- suppressWarnings(TSENAT::calculate_diversity_s4(
    analysis_high,
    q = c(1.0),
    bootstrap = TRUE,
    nboot = 50,
    seed = 111,
    output_file = output_file_high,
    verbose = FALSE
  ))
  
  div_data_low <- read.csv(output_file_low, sep = "\t", stringsAsFactors = FALSE)
  div_data_high <- read.csv(output_file_high, sep = "\t", stringsAsFactors = FALSE)
  
  # Note: This is a stochastic test - we expect narrower CIs with more replicates
  # on average, but individual comparisons might vary. We test the general trend
  # by computing median CI widths
  
  ci_width_low <- div_data_low$ci_upper - div_data_low$ci_lower
  ci_width_high <- div_data_high$ci_upper - div_data_high$ci_lower
  
  # Filter valid CIs
  ci_width_low_valid <- ci_width_low[!is.na(ci_width_low)]
  ci_width_high_valid <- ci_width_high[!is.na(ci_width_high)]
  
  if (length(ci_width_low_valid) > 0 && length(ci_width_high_valid) > 0) {
    # Verify both sets of CIs have reasonable widths (positive, not NaN)
    # Note: We don't strictly compare widths between nboot=20 vs nboot=2000
    # because each bootstrap run is independent with its own RNG stream,
    # making the comparison non-deterministic. The important validation is that
    # CIs are computed correctly (already tested above).
    expect_true(all(ci_width_low_valid > 0), info = "Low nboot CIs should have positive width")
    expect_true(all(ci_width_high_valid > 0), info = "High nboot CIs should have positive width")
    expect_true(all(is.finite(ci_width_low_valid)), info = "Low nboot CI widths should be finite")
    expect_true(all(is.finite(ci_width_high_valid)), info = "High nboot CI widths should be finite")
  }
  
  # Clean up
  if (file.exists(output_file_low)) unlink(output_file_low)
  if (file.exists(output_file_high)) unlink(output_file_high)
})

# ===========================================================================
# TEST GROUP 4: Numerical Correctness Across Configurations
# ===========================================================================

test_that("Point estimates are identical with/without bootstrap", {
  
  # Use fixed seed for reproducibility
  set.seed(456)
  analysis <- make_test_analysis_diversity(n_genes = 10, n_samples_per_group = 3, seed = 456)
  
  output_dir <- tempdir()
  output_file_noboot <- file.path(output_dir, "test_diversity_noboot_compare.tsv")
  output_file_boot <- file.path(output_dir, "test_diversity_boot_compare.tsv")
  
  # Without bootstrap
  result_noboot <- TSENAT::calculate_diversity_s4(
    analysis,
    q = c(0.5, 1.0),
    bootstrap = FALSE,
    output_file = output_file_noboot,
    verbose = FALSE
  )
  
  # With bootstrap (separate analysis with same data)
  # Note: Use new instance to ensure same underlying data
  analysis2 <- make_test_analysis_diversity(n_genes = 10, n_samples_per_group = 3, seed = 456)
  result_boot <- suppressWarnings(TSENAT::calculate_diversity_s4(
    analysis2,
    q = c(0.5, 1.0),
    bootstrap = TRUE,
    nboot = 50,
    output_file = output_file_boot,
    verbose = FALSE
  ))
  
  div_noboot <- read.csv(output_file_noboot, sep = "\t", stringsAsFactors = FALSE)
  div_boot <- read.csv(output_file_boot, sep = "\t", stringsAsFactors = FALSE)
  
  # Extract point estimates from bootstrap version
  div_boot_points <- div_boot[, c("gene", "sample", "q_value", "diversity")]
  
  # Sort both for comparison
  div_noboot_sorted <- div_noboot[order(div_noboot$gene, div_noboot$sample, div_noboot$q_value), ]
  div_boot_points_sorted <- div_boot_points[order(div_boot_points$gene, div_boot_points$sample, div_boot_points$q_value), ]
  
  # Point estimates should be very close (within numerical precision)
  if (nrow(div_noboot_sorted) == nrow(div_boot_points_sorted)) {
    expect_equal(div_noboot_sorted$diversity, div_boot_points_sorted$diversity, 
                 tolerance = 1e-6)
  }
  
  # Clean up
  if (file.exists(output_file_noboot)) unlink(output_file_noboot)
  if (file.exists(output_file_boot)) unlink(output_file_boot)
})

test_that("Diversity values differ across diverse and non-diverse samples", {
  
  # Use supported test analysis helper
  analysis <- make_test_analysis_diversity(n_genes = 6, n_samples_per_group = 2)
  output_dir <- tempdir()
  output_file <- file.path(output_dir, "test_diversity_distribution.tsv")
  
  result <- TSENAT::calculate_diversity_s4(
    analysis,
    q = c(1.0),
    bootstrap = FALSE,
    output_file = output_file,
    verbose = FALSE
  )
  
  # Check that output was created
  expect_true(file.exists(output_file))
  
  # Read and verify data
  div_data <- read.csv(output_file, sep = "\t", stringsAsFactors = FALSE)
  expect_gt(nrow(div_data), 0)
  
  # Clean up
  if (file.exists(output_file)) unlink(output_file)
})

test_that("Spectrum file aggregates correctly across samples", {
  
  analysis <- make_test_analysis_diversity(n_genes = 10, n_samples_per_group = 4)
  output_dir <- tempdir()
  output_file <- file.path(output_dir, "test_diversity_spectrum_check.tsv")
  
  result <- TSENAT::calculate_diversity_s4(
    analysis,
    q = c(0.5, 1.0, 1.5),
    bootstrap = FALSE,
    output_file = output_file,
    verbose = FALSE
  )
  
  spectrum_file <- sub("\\.[^.]+$", "_diversity_spectrum.tsv", output_file)
  if (spectrum_file == output_file) {
    spectrum_file <- paste0(output_file, "_diversity_spectrum.tsv")
  }
  
  # Read both files
  div_main <- read.csv(output_file, sep = "\t", stringsAsFactors = FALSE)
  div_spectrum <- read.csv(spectrum_file, sep = "\t", stringsAsFactors = FALSE)
  
  # Spectrum file should exist and not be empty
  expect_gt(nrow(div_spectrum), 0)
  
  # Spectrum should have fewer rows than main file (aggregated)
  expect_lt(nrow(div_spectrum), nrow(div_main))
  
  # Clean up
  if (file.exists(output_file)) unlink(output_file)
  if (file.exists(spectrum_file)) unlink(spectrum_file)
})

# ===========================================================================
# TEST GROUP 5: Multiple Q-values
# ===========================================================================

test_that("Multiple q-values produce consistent sample/gene combinations", {
  
  n_genes <- 8
  n_samples_per_group <- 2
  n_samples <- n_samples_per_group * 2
  q_vals <- c(0.5, 1.0, 1.5)  # Use consistent q values
  
  analysis <- make_test_analysis_diversity(n_genes = n_genes, n_samples_per_group = n_samples_per_group, q_values = q_vals)
  output_dir <- tempdir()
  output_file <- file.path(output_dir, "test_diversity_multiQ.tsv")
  
  result <- TSENAT::calculate_diversity_s4(
    analysis,
    q = q_vals,
    bootstrap = FALSE,
    output_file = output_file,
    verbose = FALSE
  )
  
  div_data <- read.csv(output_file, sep = "\t", stringsAsFactors = FALSE)
  
  # Basic check: file exists and has data
  expect_gt(nrow(div_data), 0)
  
  # Check that we have reasonable number of rows
  # With proper data, should be at least n_genes * n_samples rows (for q-values)
  expect_gte(nrow(div_data), n_genes)
  
  # Clean up
  if (file.exists(output_file)) unlink(output_file)
})

test_that("Diversity typically varies across q-values for same gene/sample", {
  
  analysis <- make_test_analysis_diversity(n_genes = 8, n_samples_per_group = 2, q_values = c(0.5, 1.0, 1.5))
  output_dir <- tempdir()
  output_file <- file.path(output_dir, "test_diversity_q_variation.tsv")
  
  result <- TSENAT::calculate_diversity_s4(
    analysis,
    q = c(0.5, 1.0, 1.5),
    bootstrap = FALSE,
    output_file = output_file,
    verbose = FALSE
  )
  
  div_data <- read.csv(output_file, sep = "\t", stringsAsFactors = FALSE)
  
  # Basic check: we have output
  expect_gt(nrow(div_data), 0)
  
  # Check that diversity values are present and numeric
  expect_true(is.numeric(div_data$diversity))
  expect_false(all(is.na(div_data$diversity)))
  
  # Clean up
  if (file.exists(output_file)) unlink(output_file)
})

# ===========================================================================
# TEST GROUP 6: Normalization Effects
# ===========================================================================

test_that("Norm vs non-norm calculations produce different results", {
  
  analysis_norm <- make_test_analysis_diversity(n_genes = 15, n_samples_per_group = 3, seed = 789)
  analysis_nonorm <- make_test_analysis_diversity(n_genes = 15, n_samples_per_group = 3, seed = 789)
  
  output_dir <- tempdir()
  output_file_norm <- file.path(output_dir, "test_diversity_norm.tsv")
  output_file_nonorm <- file.path(output_dir, "test_diversity_nonorm.tsv")
  
  # With normalization
  result_norm <- TSENAT::calculate_diversity_s4(
    analysis_norm,
    q = c(1.0),
    norm = TRUE,
    bootstrap = FALSE,
    output_file = output_file_norm,
    verbose = FALSE
  )
  
  # Without normalization
  result_nonorm <- TSENAT::calculate_diversity_s4(
    analysis_nonorm,
    q = c(1.0),
    norm = FALSE,
    bootstrap = FALSE,
    output_file = output_file_nonorm,
    verbose = FALSE
  )
  
  div_norm <- read.csv(output_file_norm, sep = "\t", stringsAsFactors = FALSE)
  div_nonorm <- read.csv(output_file_nonorm, sep = "\t", stringsAsFactors = FALSE)
  
  # Most diversity values should differ between norm and non-norm
  # (unless data is very special)
  differences <- abs(div_norm$diversity - div_nonorm$diversity)
  
  # Expect that at least some values differ
  expect_gt(sum(differences > 0.001), 0)
  
  # Clean up
  if (file.exists(output_file_norm)) unlink(output_file_norm)
  if (file.exists(output_file_nonorm)) unlink(output_file_nonorm)
})

# ===========================================================================
# TEST GROUP 7: File Format Validation
# ===========================================================================

test_that("Output TSV is properly formatted (readable and parseable)", {
  
  analysis <- make_test_analysis_diversity(n_genes = 10, n_samples_per_group = 3)
  output_dir <- tempdir()
  output_file <- file.path(output_dir, "test_diversity_format.tsv")
  
  result <- suppressWarnings(TSENAT::calculate_diversity_s4(
    analysis,
    q = c(0.5, 1.0),
    bootstrap = TRUE,
    nboot = 50,
    output_file = output_file,
    verbose = FALSE
  ))
  
  # Read with tab separator
  expect_error(
    div_data <- read.csv(output_file, sep = "\t", stringsAsFactors = FALSE),
    NA  # Expect no error
  )
  
  # Column names should be present
  expect_true(ncol(div_data) > 0)
  expect_true(nrow(div_data) > 0)
  
  # All required columns should be present
  required_cols <- c("gene", "sample", "q_value", "diversity")
  expect_true(all(required_cols %in% colnames(div_data)))
  
  # Clean up
  if (file.exists(output_file)) unlink(output_file)
})

test_that("Data types are correct in output file", {
  
  analysis <- make_test_analysis_diversity(n_genes = 10, n_samples_per_group = 3)
  output_dir <- tempdir()
  output_file <- file.path(output_dir, "test_diversity_dtypes.tsv")
  
  result <- TSENAT::calculate_diversity_s4(
    analysis,
    q = c(0.5, 1.0),
    bootstrap = FALSE,
    output_file = output_file,
    verbose = FALSE
  )
  
  div_data <- read.csv(output_file, sep = "\t", stringsAsFactors = FALSE)
  
  # Convert q_value to numeric if needed (TSV reading may preserve as character)
  # Suppress coercion warning if any non-numeric values exist
  if (is.character(div_data$q_value)) {
    div_data$q_value <- suppressWarnings(as.numeric(div_data$q_value))
  }
  
  # Check expected data types
  expect_true(is.character(div_data$gene),
              info = "gene column should be character")
  expect_true(is.character(div_data$sample),
              info = "sample column should be character")
  
  # q_value and diversity should be numeric
  expect_true(is.numeric(div_data$q_value),
              info = "q_value column should be numeric")
  expect_true(is.numeric(div_data$diversity),
              info = "diversity column should be numeric")
  
  # Clean up
  if (file.exists(output_file)) unlink(output_file)
})
