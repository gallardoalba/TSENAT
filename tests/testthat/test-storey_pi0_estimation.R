library(TSENAT)

context("storey_pi0_estimation: Storey Prior Probability Estimation")

test_that("estimate_storey_pi0: lambda method works correctly", {
  # Create test data: mixture of 400 nulls and 100 signals
  set.seed(42)
  n_null <- 400
  n_signal <- 100
  pvalues <- c(
    runif(n_null),           # Null: uniform [0,1]
    rbeta(n_signal, 0.5, 1)  # Signal: skewed toward small p
  )
  
  # Estimate π₀
  result <- estimate_storey_pi0(pvalues, lambda = 0.5, pi0_method = "lambda")
  
  expect_is(result, "list")
  expect_true("pi0" %in% names(result))
  expect_true("lambda" %in% names(result))
  expect_true("pi0_method" %in% names(result))
  
  # π₀ should be close to 0.80 (400/500)
  expect_gt(result$pi0, 0.7)
  expect_lt(result$pi0, 0.9)
  expect_equal(result$pi0_method, "lambda")
  expect_equal(result$lambda, 0.5)
})

test_that("estimate_storey_pi0: handles edge cases", {
  # All nulls
  pvalues_all_null <- runif(100)
  result_all <- estimate_storey_pi0(pvalues_all_null, lambda = 0.5)
  expect_equal(result_all$pi0, 1.0)
  
  # All signals (very small p-values)
  set.seed(123)
  pvalues_all_signal <- rbeta(100, 0.1, 1)  # Heavily skewed toward 0
  result_all_sig <- estimate_storey_pi0(pvalues_all_signal, lambda = 0.5)
  expect_lt(result_all_sig$pi0, 0.3)  # Most p-values in [0, 0.5]
})

test_that("estimate_storey_pi0: rejects invalid inputs", {
  expect_error(estimate_storey_pi0(c(1.5, 0.5, 0.3)))  # p > 1
  expect_error(estimate_storey_pi0(c(-0.1, 0.5, 0.3)))  # p < 0
  expect_error(estimate_storey_pi0(numeric(0)))  # Empty
  expect_error(estimate_storey_pi0(c(0.5), lambda = 1.5))  # lambda >= 1
  expect_error(estimate_storey_pi0(c(0.5), lambda = -0.1))  # lambda < 0
})

test_that("compute_storey_qvalues: produces valid q-values", {
  set.seed(42)
  pvalues <- c(runif(400), rbeta(100, 0.5, 1))
  
  qvalues <- compute_storey_qvalues(pvalues, pi0 = 0.8)
  
  # Check properties
  expect_length(qvalues, length(pvalues))
  expect_true(all(qvalues >= 0, na.rm = TRUE))
  expect_true(all(qvalues <= 1, na.rm = TRUE))
  expect_equal(sum(is.na(qvalues)), 0)
})

test_that("compute_storey_qvalues: enforces monotonicity", {
  set.seed(42)
  pvalues <- sort(runif(100))  # Sort to make monotonicity obvious
  
  qvalues <- compute_storey_qvalues(pvalues, pi0 = NULL, robust = TRUE)
  
  # Check: if p_i < p_j then q_i <= q_j (allowing for numerical tolerance)
  for (i in 1:(length(pvalues) - 1)) {
    expect_lte(qvalues[i], qvalues[i + 1] + 1e-10)
  }
})

test_that("compute_storey_qvalues: handles NAs correctly", {
  pvalues_with_na <- c(0.01, NA, 0.05, 0.1, NA, 0.3)
  
  qvalues <- compute_storey_qvalues(pvalues_with_na, pi0 = 0.8, na.rm = TRUE)
  
  # Check NAs preserved in same positions
  expect_true(is.na(qvalues[2]))
  expect_true(is.na(qvalues[5]))
  expect_false(is.na(qvalues[1]))
  expect_false(is.na(qvalues[3]))
})
