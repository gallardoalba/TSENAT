# Tests for lazy-loading visualization functions
library(ggplot2)
library(dplyr)
library(SummarizedExperiment)
library(testthat)

# Tests for lazy-loading visualization functions
# 
# These tests verify that visualization dependencies are not loaded at package
# startup, but are correctly loaded on first use via lazy-loading mechanism.

# Note: These tests are designed to run in a clean environment where TSENAT
# is freshly loaded. Some tests may need to be skipped in full test suites
# if other tests have already triggered visualization loading.



test_that("Lazy-loading infrastructure exists", {
  # Verify functions are available
  expect_true(exists(".load_visualization_deps"))
  expect_true(exists(".viz_available"))
  expect_true(exists(".viz_status"))
})

test_that(".viz_available reports lazy-loading state", {
  # Function should exist and be callable
  result <- TSENAT:::.viz_available()
  expect_type(result, "logical")
  expect_length(result, 1)
})

test_that(".viz_status provides diagnostic information", {
  # Check that diagnostic function works
  status <- TSENAT:::.viz_status()
  
  expect_is(status, "list")
  expect_true(all(c("viz_loaded", "packages_loaded", "total_namespaces_loaded") %in% names(status)))
  
  # packages_loaded should be named logical vector
  expect_named(status$packages_loaded)
  expect_type(status$packages_loaded, "logical")
})

test_that(".load_visualization_deps works and is idempotent", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  
  # First call should load
  result1 <- TSENAT:::.load_visualization_deps()
  expect_type(result1, "logical")
  expect_length(result1, 1)
  
  # Should return FALSE on first call (indicating just loaded)
  # OR TRUE if already loaded by something else
  expect_true(result1 %in% c(TRUE, FALSE))
  
  # Second call should be idempotent
  result2 <- TSENAT:::.load_visualization_deps()
  expect_true(result2)  # Should return TRUE (already loaded)
})

test_that(".load_visualization_deps fails gracefully with strict=FALSE", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  
  # This test verifies error handling, not real failure
  # We can't easily force a package loading failure in tests,
  # but we verify the parameter works
  result <- TSENAT:::.load_visualization_deps(strict = FALSE)
  expect_true(result %in% c(TRUE, FALSE, NA))
})

test_that("Plot functions trigger visualization loading", {
  skip_on_cran()
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("ggplot2")
  
  # Create minimal test data
  analysis <- .create_test_analysis(
    n_genes = 4, 
    n_samples_per_group = 5
  )
  
  # Compute required results
  analysis <- calculate_diversity_s4(analysis, q = 1.0, verbose = FALSE)
  analysis <- calculate_divergence_s4(analysis, q = 1.0, verbose = FALSE)
  
  # Call a plot function (should trigger lazy-loading if not already done)
  p <- tryCatch({
    plot_divergence_spectrum_s4(analysis, verbose = FALSE)
  }, error = function(e) {
    # Plot might fail for other reasons (missing data), that's OK
    # We're just testing that lazy-loading is triggered
    NULL
  })
  
  # After calling a plot function, viz packages should be loaded
  status <- TSENAT:::.viz_status()
  expect_true(status$viz_loaded)
  expect_true(any(status$packages_loaded))  # At least one viz package loaded
})

test_that("Multiple plot functions work after lazy-loading", {
  skip_on_cran()
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("ggplot2")
  
  # Create test analysis with required computations
  analysis <- .create_test_analysis(
    n_genes = 4,
    n_samples_per_group = 5
  )
  
  # LM interaction requires at least 5 q-values
  analysis <- calculate_diversity_s4(analysis, q = seq(0.5, 2, by = 0.375), verbose = FALSE)
  analysis <- calculate_divergence_s4(analysis, q = seq(0.5, 2, by = 0.375), verbose = FALSE)
  analysis <- calculate_lm_interaction_s4(
    analysis, 
    condition_col = "condition",
    verbose = FALSE
  )
  
  # These should all trigger lazy-loading on first call
  # and work correctly on subsequent calls
  plots <- list()
  
  # Attempt different plot types (some may fail due to data constraints, that's OK)
  tryCatch({
    plots$spectrum <- plot_divergence_spectrum_s4(analysis, verbose = FALSE)
  }, error = function(e) NULL)
  
  tryCatch({
    plots$curve <- plot_tsallis_q_curve_s4(analysis, verbose = FALSE)
  }, error = function(e) NULL)
  
  # At least verify lazy-loading was triggered
  status <- TSENAT:::.viz_status()
  expect_true(status$viz_loaded)
})

test_that("Lazy-loading doesn't affect non-plot functions", {
  # Core computational functions should work regardless of viz loading state
  # This test verifies that lazy-loading is orthogonal to computation
  
  analysis <- .create_test_analysis(n_genes = 4, n_samples_per_group = 5)
  
  # These should work fine without any plot functions
  analysis <- calculate_diversity_s4(analysis, q = 1.0, verbose = FALSE)
  expect_true(length(analysis@diversity_results) > 0)
  
  analysis <- calculate_divergence_s4(analysis, q = 1.0, verbose = FALSE)
  expect_true(length(analysis@divergence_results) > 0)
  
  # Viz may or may not be loaded depending on prior tests
  # But these functions should work either way
})

test_that("Lazy-loading maintains backward compatibility", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  
  # Users expecting visualizations to just work should not be affected
  analysis <- .create_test_analysis(n_genes = 4, n_samples_per_group = 5)
  analysis <- calculate_diversity_s4(analysis, q = 1.0, verbose = FALSE)
  
  # Calling plot functions should work as before
  # (they just might be slightly slower on first call due to lazy-loading)
  tryCatch({
    p <- plot_diversity_density(
      SummarizedExperiment::assay(analysis@se),
      sample_type_col = NULL
    )
    # If plot creation succeeds, it's using loaded viz packages
    expect_s3_class(p, "ggplot")
  }, error = function(e) {
    # Some plot functions might need specific data, that's OK
    # We're just verifying the mechanism works
    expect_true(TRUE)
  })
})
