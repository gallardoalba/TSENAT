context("label_shuffling paired permutations")

test_that("paired permutation requires even number of samples", {
    x <- matrix(runif(6), nrow = 2)
    samples <- c("A", "B", "A")
    expect_error(
        label_shuffling(x, samples, control = "A", method = "mean", paired = TRUE),
        "Paired permutation requires an even number of samples"
    )
})

test_that("test_differential forwards paired to label_shuffling and is reproducible", {
    set.seed(123)
    x <- matrix(rnorm(8), nrow = 2)
    samples <- rep(c("Normal", "Tumor"), times = 2)

    # direct call with same internal seed (default 123 in test_differential)
    set.seed(123)
    direct <- label_shuffling(x, samples, control = "Normal", method = "mean", randomizations = 50, pcorr = "none", paired = TRUE)

    wrapped <- test_differential(x, samples, control = "Normal", method = "shuffle", fc_method = "mean", randomizations = 50, pcorr = "none", paired = TRUE)

    expect_equal(direct, wrapped)
})

context("label_shuffling signflip exact and sampled")

test_that("signflip exact enumeration runs and returns valid p-values", {
    set.seed(42)
    x <- matrix(rnorm(8), nrow = 2)
    samples <- rep(c("Normal", "Tumor"), times = 2) # two pairs
    # total combinations = 2^2 = 4
    res <- label_shuffling(x, samples, control = "Normal", method = "mean", randomizations = 4, paired = TRUE, paired_method = "signflip")
    expect_true(is.matrix(res))
    expect_equal(ncol(res), 2)
    expect_equal(nrow(res), nrow(x))
    expect_true(all(res >= 0 & res <= 1, na.rm = TRUE))
})

test_that("signflip sampled returns same shape and in-range p-values", {
    set.seed(42)
    x <- matrix(rnorm(8), nrow = 2)
    samples <- rep(c("Normal", "Tumor"), times = 2)
    res <- label_shuffling(x, samples, control = "Normal", method = "mean", randomizations = 10, paired = TRUE, paired_method = "signflip")
    expect_true(is.matrix(res))
    expect_equal(ncol(res), 2)
    expect_equal(nrow(res), nrow(x))
    expect_true(all(res >= 0 & res <= 1, na.rm = TRUE))
})

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
