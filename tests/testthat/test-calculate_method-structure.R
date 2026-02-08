test_that("calculate_method returns expected column names for multiple q", {
  # small synthetic transcript matrix: 4 transcripts, 2 genes (2 transcripts each), 2 samples
  mat <- matrix(c(
    10, 5,
    0, 0,
    2, 8,
    3, 7
  ), nrow = 4, byrow = TRUE)
  colnames(mat) <- c("S1", "S2")
  genes <- c("g1", "g1", "g2", "g2")
  res <- calculate_method(mat, genes, norm = TRUE, q = c(0.5, 1), what = "S")
  # expected column names: S1_q=0.5, S1_q=1, S2_q=0.5, S2_q=1
  expect_true(all(c("S1_q=0.5", "S1_q=1", "S2_q=0.5", "S2_q=1") %in% colnames(res)))
  expect_equal(res$Gene, unique(genes))
})

test_that("calculate_method handles missing sample names by producing columns nevertheless", {
  mat <- matrix(rep(1, 6), nrow = 3)
  colnames(mat) <- NULL
  genes <- letters[1:3]
  res <- calculate_method(mat, genes, q = 1, what = "S")
  # should return a data.frame with columns (Gene + one column)
  expect_true(is.data.frame(res))
  expect_true(ncol(res) >= 2)
})
