context("calculate_fc defensive checks")

library(testthat)

test_that("calculate_fc errors when control missing or not found", {
  x <- matrix(runif(8), nrow = 2)
  samples <- c("A","A","B","B")
  expect_error(calculate_fc(x, samples, NULL), "`control` must be provided to calculate_fc")
  expect_error(calculate_fc(x, samples, "C"), "Control sample type not found in samples\\.")
  # mismatched samples length
  expect_error(calculate_fc(x, samples[-1], "A"), "Length of 'samples' must equal number of columns in 'x'")
})
