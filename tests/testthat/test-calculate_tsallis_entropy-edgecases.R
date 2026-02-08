test_that("calculate_tsallis_entropy handles zero-sum and q=1 correctly", {
  x_uniform <- c(1, 1, 1)
  # Uniform distribution normalized entropy should be 1 for any q when norm=TRUE
  s_unif <- calculate_tsallis_entropy(x_uniform, q = c(0.5, 1, 2), norm = TRUE, what = "S")
  expect_equal(as.numeric(s_unif), rep(1, 3))

  # Zero-sum input returns NA
  x_zero <- c(0, 0, 0)
  s_zero <- calculate_tsallis_entropy(x_zero, q = c(0.5, 1, 2), norm = TRUE, what = "S")
  expect_true(all(is.na(s_zero)))

  # D at q = 1 equals exp(Shannon) when using natural log base
  x <- c(10, 5, 0)
  p <- x / sum(x)
  sh <- -sum(ifelse(p > 0, p * log(p), 0))
  expected_D1 <- exp(sh)
  D1 <- calculate_tsallis_entropy(x, q = 1, what = "D")
  expect_equal(as.numeric(D1), expected_D1)
})

test_that("calculate_diversity rejects non-positive q", {
  mat <- matrix(1, nrow = 3, ncol = 2)
  genes <- letters[1:3]
  expect_error(calculate_diversity(mat, genes = genes, q = 0), "q")
})
