# =============================================================================
# Jackknife Test Consolidation Helpers
# =============================================================================
# Reusable test helpers for jackknife_tsallis_entropy to eliminate redundancy
# across test-jackknife_diagnostics.R and test-jackknife-coverage.R
# =============================================================================

#' Test Input Type Handling for Jackknife Functions
#'
#' @param func_name Function name as string (e.g., "jackknife_tsallis_entropy")
#' @param test_vec Numeric vector for testing (optional)
#' @param test_mat Numeric matrix for testing (optional)
#' @param extra_args List of additional arguments to pass to function
#'
#' @keywords internal
test_jackknife_input_types <- function(
    func_name = "jackknife_tsallis_entropy",
    test_vec = c(100, 50, 30, 20),
    test_mat = matrix(c(100, 50, 30, 20, 80, 40, 20, 10), nrow = 2, byrow = TRUE),
    extra_args = list(q = 1, print_results = FALSE)
) {
  skip_if_not_installed("TSENAT")
  
  # Test vector input
  vec_result <- do.call(func_name, c(list(x = test_vec), extra_args))
  expect_true(!is.null(vec_result), info = "Vector input should return result")
  
  # Test matrix input
  mat_result <- do.call(func_name, c(list(x = test_mat), extra_args))
  expect_true(!is.null(mat_result), info = "Matrix input should return result")
  expect_is(mat_result, "tsenat_jackknife_list", info = "Matrix should return list")
  expect_equal(length(mat_result), nrow(test_mat), info = "Matrix result should have nrow(x) elements")
  
  # Test data.frame input
  df_result <- do.call(func_name, c(list(x = as.data.frame(test_mat)), extra_args))
  expect_true(!is.null(df_result), info = "Data frame input should return result")
  expect_is(df_result, "tsenat_jackknife_list", info = "Data frame should return list")
}

#' Test Parameter Validation for Jackknife Functions
#'
#' @param func_name Function name as string
#' @param valid_counts Valid count vector to use as baseline
#' @param valid_args List of other valid arguments
#'
#' @keywords internal
test_jackknife_parameter_validation <- function(
    func_name = "jackknife_tsallis_entropy",
    valid_counts = c(100, 50, 30, 20),
    valid_args = list(q = 1, print_results = FALSE)
) {
  skip_if_not_installed("TSENAT")
  
  # Test rejection of negative counts
  expect_error(
    do.call(func_name, c(list(x = c(100, -50, 30, 20)), valid_args)),
    info = "Should reject negative counts"
  )
  
  # Test rejection of NA values
  expect_error(
    do.call(func_name, c(list(x = c(100, NA, 30, 20)), valid_args)),
    info = "Should reject NA values"
  )
  
  # Test rejection of invalid q parameter
  expect_error(
    do.call(func_name, c(list(x = valid_counts, q = -1), valid_args)),
    info = "Should reject negative q"
  )
  
  # Test minimum size requirement (n >= 2)
  expect_error(
    do.call(func_name, c(list(x = c(100)), valid_args)),
    info = "Should require at least 2 transcripts"
  )
}

#' Test Multi-q Support for Resampling Functions
#'
#' @param func_name Function name (e.g., "jackknife_tsallis_entropy")
#' @param test_data Data to use (vector, matrix, or SE)
#' @param q_vector Vector of q values to test (e.g., c(1, 2, 3))
#' @param extra_args List of additional arguments
#'
#' @keywords internal
test_multiq_support <- function(
    func_name = "jackknife_tsallis_entropy",
    test_data = c(100, 50, 30, 20),
    q_vector = c(0.5, 1, 2),
    extra_args = list(print_results = FALSE)
) {
  skip_if_not_installed("TSENAT")
  
  # Test single q (backward compatibility)
  single_q_result <- do.call(func_name, c(list(x = test_data, q = q_vector[1]), extra_args))
  expect_is(single_q_result, "tsenat_jackknife", info = "Single q should return tsenat_jackknife")
  
  # Test multiple q values
  multi_q_result <- do.call(func_name, c(list(x = test_data, q = q_vector), extra_args))
  expect_is(multi_q_result, "tsenat_jackknife_list", info = "Multiple q should return list")
  expect_equal(length(multi_q_result), length(q_vector), 
               info = "Result should have length equal to q vector")
  
  # Verify different q values produce different estimates
  estimates <- sapply(multi_q_result, function(r) r$estimate)
  estimate_range <- max(estimates) - min(estimates)
  expect_true(estimate_range > 0, 
              info = "Different q values should produce different estimates")
}

#' Assert Jackknife Result Structure
#'
#' @param result Result from jackknife function
#' @param n_transcripts Expected number of transcripts
#' @param check_influence If TRUE, verify influence field exists
#'
#' @keywords internal
assert_jackknife_result_valid <- function(result, n_transcripts = NULL, check_influence = TRUE) {
  expect_is(result, "tsenat_jackknife", info = "Result should be tsenat_jackknife object")
  
  # Check required fields
  expect_true(!is.null(result$estimate), info = "estimate field required")
  expect_true(!is.null(result$jackknife_estimates), info = "jackknife_estimates field required")
  expect_true(!is.null(result$jackknife_se), info = "jackknife_se field required")
  
  # Check types
  expect_is(result$estimate, "numeric", info = "estimate must be numeric")
  expect_is(result$jackknife_se, "numeric", info = "jackknife_se must be numeric")
  expect_is(result$jackknife_estimates, "numeric", info = "jackknife_estimates must be numeric")
  
  # Check sizes match if specified
  if (!is.null(n_transcripts)) {
    expect_equal(length(result$jackknife_estimates), n_transcripts,
                 info = "jackknife_estimates length should match n_transcripts")
  }
  
  # Check influence optionally
  if (check_influence) {
    expect_true(!is.null(result$influence), info = "influence field expected")
    expect_is(result$influence, "numeric", info = "influence should be numeric")
  }
  
  invisible(result)
}

#' Assert Jackknife List Structure
#'
#' @param result Result from jackknife with matrix/multi-gene input
#' @param expected_length Expected number of genes/rows
#'
#' @keywords internal
assert_jackknife_list_valid <- function(result, expected_length = NULL) {
  expect_is(result, "tsenat_jackknife_list", info = "Result should be tsenat_jackknife_list")
  
  # Check it's a list
  expect_true(is.list(result), info = "Result should be a list")
  
  # Check each element is valid jackknife result
  for (i in seq_along(result)) {
    assert_jackknife_result_valid(result[[i]], check_influence = FALSE)
  }
  
  # Check length if specified
  if (!is.null(expected_length)) {
    expect_equal(length(result), expected_length,
                 info = "List length should match expected_length")
  }
  
  invisible(result)
}
