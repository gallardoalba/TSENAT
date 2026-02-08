test_that("label_shuffling() returns named p-value columns", {
  mat <- matrix(c(
    0.1, 0.2, 0.3, 0.4,
    1.0, 1.2, 0.9, 1.1
  ), nrow = 2, byrow = TRUE)
  colnames(mat) <- paste0("S", 1:4)
  samples <- c("A", "A", "B", "B")
  set.seed(1)
  res <- label_shuffling(mat, samples, control = "A", method = "mean", randomizations = 10, pcorr = "none")
  expect_true(is.matrix(res) || is.data.frame(res))
  expect_true(all(c("raw_p_values", "adjusted_p_values") %in% colnames(res)))
})
