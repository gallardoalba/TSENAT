test_that("wilcoxon() returns named p-value columns", {
  mat <- matrix(runif(20), nrow = 5)
  samples <- rep(c("A", "B"), each = 5)
  res <- wilcoxon(mat, samples, pcorr = "none", paired = FALSE, exact = FALSE)
  expect_true(is.matrix(res) || is.data.frame(res))
  expect_true(all(c("raw_p_values", "adjusted_p_values") %in% colnames(res)))
})
