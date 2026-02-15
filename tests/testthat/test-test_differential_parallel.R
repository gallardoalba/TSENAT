context("test_differential parallelization: wilcoxon")

test_that("test_differential wilcoxon with nthreads=1", {
    x <- matrix(rnorm(200), nrow = 20)
    samples <- c(rep("A", 5), rep("B", 5))
    
    result <- test_differential(
        x, samples,
        method = "wilcoxon",
        nthreads = 1
    )
    
    expect_is(result, c("matrix", "array"))
    expect_equal(ncol(result), 2)
    expect_equal(nrow(result), 20)
    expect_true(all(colnames(result) == c("raw_p_values", "adjusted_p_values")))
})

test_that("test_differential wilcoxon with nthreads=2", {
    x <- matrix(rnorm(200), nrow = 20)
    samples <- c(rep("A", 5), rep("B", 5))
    
    result_serial <- test_differential(
        x, samples,
        method = "wilcoxon",
        nthreads = 1
    )
    
    result_parallel <- test_differential(
        x, samples,
        method = "wilcoxon",
        nthreads = 2
    )
    
    # Serial and parallel should produce identical results
    expect_equal(result_serial, result_parallel, tolerance = 1e-10)
})

context("test_differential parallelization: label_shuffling")

test_that("test_differential shuffle with nthreads=1", {
    x <- matrix(rnorm(200), nrow = 20)
    samples <- c(rep("Ctrl", 5), rep("Case", 5))
    
    set.seed(42)
    result <- test_differential(
        x, samples,
        control = "Ctrl",
        method = "shuffle",
        randomizations = 100,
        nthreads = 1
    )
    
    expect_is(result, c("matrix", "array"))
    expect_equal(ncol(result), 2)
    expect_equal(nrow(result), 20)
    expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
})

test_that("test_differential shuffle with nthreads=2", {
    x <- matrix(rnorm(200), nrow = 20)
    samples <- c(rep("Ctrl", 5), rep("Case", 5))
    
    set.seed(42)
    result_serial <- test_differential(
        x, samples,
        control = "Ctrl",
        method = "shuffle",
        randomizations = 100,
        nthreads = 1
    )
    
    set.seed(42)
    result_parallel <- test_differential(
        x, samples,
        control = "Ctrl",
        method = "shuffle",
        randomizations = 100,
        nthreads = 2
    )
    
    # Serial and parallel should produce identical results with same seed
    expect_equal(result_serial, result_parallel, tolerance = 1e-10)
})

context("test_differential parallelization: paired tests")

test_that("test_differential paired wilcoxon with nthreads", {
    x <- matrix(rnorm(60), nrow = 10)
    samples <- rep(c("A", "B"), 3)
    
    result_serial <- test_differential(
        x, samples,
        method = "wilcoxon",
        paired = TRUE,
        nthreads = 1
    )
    
    result_parallel <- test_differential(
        x, samples,
        method = "wilcoxon",
        paired = TRUE,
        nthreads = 2
    )
    
    expect_is(result_serial, c("matrix", "array"))
    expect_equal(result_serial, result_parallel, tolerance = 1e-10)
    expect_equal(nrow(result_serial), 10)
})

test_that("test_differential paired shuffle with nthreads and swap", {
    x <- matrix(rnorm(60), nrow = 10)
    samples <- rep(c("Ctrl", "Case"), 3)
    
    set.seed(42)
    result_serial <- test_differential(
        x, samples,
        control = "Ctrl",
        method = "shuffle",
        paired = TRUE,
        paired_method = "swap",
        randomizations = 100,
        nthreads = 1
    )
    
    set.seed(42)
    result_parallel <- test_differential(
        x, samples,
        control = "Ctrl",
        method = "shuffle",
        paired = TRUE,
        paired_method = "swap",
        randomizations = 100,
        nthreads = 2
    )
    
    expect_is(result_serial, c("matrix", "array"))
    expect_equal(result_serial, result_parallel, tolerance = 1e-10)
    expect_equal(nrow(result_serial), 10)
    expect_true(all(result_serial >= 0 & result_serial <= 1, na.rm = TRUE))
})

test_that("test_differential paired shuffle with nthreads and signflip", {
    x <- matrix(rnorm(60), nrow = 10)
    samples <- rep(c("Ctrl", "Case"), 3)
    
    set.seed(42)
    result_serial <- test_differential(
        x, samples,
        control = "Ctrl",
        method = "shuffle",
        paired = TRUE,
        paired_method = "signflip",
        randomizations = 50,
        nthreads = 1
    )
    
    set.seed(42)
    result_parallel <- test_differential(
        x, samples,
        control = "Ctrl",
        method = "shuffle",
        paired = TRUE,
        paired_method = "signflip",
        randomizations = 50,
        nthreads = 2
    )
    
    expect_is(result_serial, c("matrix", "array"))
    expect_equal(result_serial, result_parallel, tolerance = 1e-10)
    expect_equal(nrow(result_serial), 10)
    expect_true(all(result_serial >= 0 & result_serial <= 1, na.rm = TRUE))
})

test_that("test_differential shuffle reproducibility with nthreads", {
    x <- matrix(rnorm(200), nrow = 20)
    samples <- c(rep("A", 5), rep("B", 5))
    
    set.seed(123)
    result1 <- test_differential(
        x, samples,
        control = "A",
        method = "shuffle",
        randomizations = 100,
        seed = 999,
        nthreads = 1
    )
    
    set.seed(123)
    result2 <- test_differential(
        x, samples,
        control = "A",
        method = "shuffle",
        randomizations = 100,
        seed = 999,
        nthreads = 1
    )
    
    # Results should be identical with the same seed
    expect_equal(result1, result2, tolerance = 1e-10)
})
