test_that("label_shuffling handles constant null distribution correctly", {
  # Build a simple matrix where permutation yields identical nulls
  # Two groups of 2 samples each: group A (control), group B (case)
  mat <- matrix(c(
    0.1, 0.1, 0.1, 0.1,  # gene1: identical across samples -> permuted log2fc all zero
    1.0, 1.0, 2.0, 2.0     # gene2: different between groups -> non-zero observed
  ), nrow = 2, byrow = TRUE)
  colnames(mat) <- paste0("S", 1:4)
  samples <- c("A", "A", "B", "B")

  # gene1 observed log2fc should be 0; nulls will be all 0 -> pval should be 1/(n+1)
  # gene2 observed log2fc should be non-zero -> pval should reflect permutations

  # Use a small number of permutations for reproducibility
  set.seed(42)
  res <- label_shuffling(mat, samples, control = "A", method = "mean", randomizations = 10, pcorr = "none")
  expect_true(is.matrix(res) || is.data.frame(res))
  raw <- as.numeric(res[,1])

  # For gene1: obs = 0, nulls all zero -> count = (all nulls |null| >= |obs|) = n_perm
  # empirical p = (n_perm + 1) / (n_perm + 1) = 1
  expect_equal(raw[1], 1)

  # For gene2: expect p between 0 and 1
  expect_true(raw[2] > 0 && raw[2] <= 1)
})
