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
