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
