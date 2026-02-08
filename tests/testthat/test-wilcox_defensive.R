context("Wilcoxon defensive behavior")

test_that("wilcoxon handles all-NA and constant rows without error and returns matrix", {
  # two groups of 3 samples each
  samples <- c(rep("A", 3), rep("B", 3))
  # build matrix with 3 rows: all-NA, constant, variable
  m <- matrix(nrow = 3, ncol = 6)
  m[1, ] <- NA_real_
  m[2, ] <- rep(5, 6)
  m[3, ] <- c(1,2,3,4,5,6)

  res <- wilcoxon(m, samples, pcorr = "none")
  expect_true(is.matrix(res) || is.data.frame(res))
  # two columns: raw and adjusted
  expect_equal(ncol(res), 2)
  # raw p-values produced and NA-handling yields numeric outputs
  raw <- as.numeric(res[,1])
  expect_length(raw, 3)
  # all-NA row should produce p-value of 1 (by convention used elsewhere)
  expect_true(is.finite(raw[1]))
})
