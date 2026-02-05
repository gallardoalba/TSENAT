context("test_differential wrapper")

test_that("test_differential wilcoxon matches direct wilcoxon", {
    x <- matrix(c(
        1, 2, 3, 6,
        5, 6, 7, 9
    ), nrow = 2, byrow = TRUE)
    samples <- c("Normal", "Tumor", "Normal", "Tumor")

    res_wil <- wilcoxon(x, samples, paired = TRUE, exact = TRUE)
    res_wrap <- test_differential(
        x,
        samples,
        method = "wilcoxon",
        paired = TRUE,
        exact = TRUE
    )

    expect_equal(res_wil, res_wrap)
})

test_that("test_differential shuffle is reproducible with default seed", {
    x <- matrix(c(
        1, 2, 3, 6,
        5, 6, 7, 9
    ), nrow = 2, byrow = TRUE)
    samples <- c("Normal", "Tumor", "Normal", "Tumor")

    # default seed is set inside the wrapper; two successive calls should match
    res_a <- test_differential(
        x,
        samples,
        control = "Normal",
        method = "shuffle",
        randomizations = 50
    )
    res_b <- test_differential(
        x,
        samples,
        control = "Normal",
        method = "shuffle",
        randomizations = 50
    )
    expect_equal(res_a, res_b)

    # different seed should usually change the permutation results
    res_c <- test_differential(
        x,
        samples,
        control = "Normal",
        method = "shuffle",
        randomizations = 50,
        seed = 124
    )
    expect_false(all(res_a == res_c))
})

test_that("test_differential shuffle requires control argument", {
    x <- matrix(runif(8), nrow = 2)
    samples <- rep(c("A", "B"), length.out = ncol(x))
    expect_error(
        test_differential(x, samples, method = "shuffle"),
        "`control` must be provided when method = 'shuffle'"
    )
})
