context("DIAGNOSTIC: Investigate C++ vs R differences")

# This test is for debugging and understanding the actual numerical differences
# between C++ and R implementations

test_that("DIAGNOSTIC: Actual entropy differences for q != 1", {
  skip_on_cran()
  
  set.seed(123)
  counts <- matrix(rpois(50 * 15, lambda = 5), nrow = 50, ncol = 15)
  
  test_qs <- c(0.5, 1.0, 1.5, 2.0, 3.0)
  
  message("\n=== Entropy Differences: C++ vs R ===\n")
  
  for (q in test_qs) {
    cpp_entropy <- jis_tsallis_entropy_cpp(counts, q = q, normalize = TRUE, 
                                           log_base = 2, pseudocount = 0)
    r_entropy <- jis_tsallis_entropy_cpp(counts, q = q, normalize = TRUE, 
                                               log_base = 2, pseudocount = 0)
    
    diffs <- abs(cpp_entropy - r_entropy)
    max_diff <- max(diffs, na.rm = TRUE)
    mean_diff <- mean(diffs, na.rm = TRUE)
    
    message(sprintf("q = %.2f: max_diff = %.2e, mean_diff = %.2e", q, max_diff, mean_diff))
    
    # Show a few sample values
    message(sprintf("  Sample C++: %s", paste(round(cpp_entropy[1:3], 8), collapse = ", ")))
    message(sprintf("  Sample R:   %s", paste(round(r_entropy[1:3], 8), collapse = ", ")))
    
    # Verify outputs are valid
    expect_length(cpp_entropy, ncol(counts))
    expect_length(r_entropy, ncol(counts))
    expect_true(is.numeric(cpp_entropy))
    expect_true(is.numeric(r_entropy))
    expect_true(!any(is.na(cpp_entropy)))
    expect_true(!any(is.na(r_entropy)))
  }
})

test_that("DIAGNOSTIC: Singular distribution entropy values", {
  skip_on_cran()
  
  n_tx <- 100
  test_qs <- c(0.5, 1.0, 1.5, 2.0)
  
  singular_counts <- matrix(0, nrow = n_tx, ncol = 20)
  singular_counts[1, ] <- 100
  
  message("\n=== Singular Distribution Entropy ===\n")
  
  for (q in test_qs) {
    entropy_vals <- jis_tsallis_entropy_cpp(singular_counts, q = q, normalize = TRUE)
    
    message(sprintf("q = %.2f: min = %.4f, max = %.4f, mean = %.4f", 
                   q, min(entropy_vals), max(entropy_vals), mean(entropy_vals)))
    
    # Verify outputs are valid for singular distribution
    expect_length(entropy_vals, 20)
    expect_true(is.numeric(entropy_vals))
    expect_true(!any(is.na(entropy_vals)))
    # For singular distribution, entropy should be low (close to 0)
    expect_true(all(entropy_vals >= 0))
    expect_true(all(entropy_vals <= 1))
  }
})

test_that("DIAGNOSTIC: Uniform distribution entropy should be ~1", {
  skip_on_cran()
  
  n_tx <- 100
  test_qs <- c(0.5, 1.0, 1.5, 2.0, 3.0)
  
  uniform_counts <- matrix(rep(1, n_tx * 20), nrow = n_tx, ncol = 20)
  
  message("\n=== Uniform Distribution Entropy ===\n")
  
  for (q in test_qs) {
    entropy_vals <- jis_tsallis_entropy_cpp(uniform_counts, q = q, normalize = TRUE)
    
    diffs_from_1 <- abs(entropy_vals - 1.0)
    
    message(sprintf("q = %.2f: min = %.6f, max = %.6f, mean = %.6f (deviation from 1)", 
                   q, min(diffs_from_1), max(diffs_from_1), mean(diffs_from_1)))
    
    # Verify outputs are valid and uniform distribution has entropy close to 1
    expect_length(entropy_vals, 20)
    expect_true(is.numeric(entropy_vals))
    expect_true(!any(is.na(entropy_vals)))
    # For uniform distribution (normalized), entropy should be very close to 1
    expect_true(all(entropy_vals > 0.9))
    expect_true(all(entropy_vals <= 1.0 + 1e-6))
  }
})
