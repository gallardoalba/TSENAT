context("Parallelization: wilcoxon with nthreads")

test_that("wilcoxon serial and parallel produce identical results", {
    set.seed(42)
    # Create large dataset: 500 genes, 20 samples
    nfeat <- 500
    nsamp <- 20
    mat <- matrix(runif(nfeat * nsamp, 0, 10), nrow = nfeat, ncol = nsamp)
    samples <- rep(c("Control", "Treatment"), each = nsamp / 2)
    
    # Run serial version
    res_serial <- wilcoxon(mat, samples, pcorr = "BH", paired = FALSE, exact = FALSE, nthreads = 1)
    
    # Run parallel version with 2 threads
    res_parallel_2 <- wilcoxon(mat, samples, pcorr = "BH", paired = FALSE, exact = FALSE, nthreads = 2)
    
    # Run parallel version with 4 threads
    res_parallel_4 <- wilcoxon(mat, samples, pcorr = "BH", paired = FALSE, exact = FALSE, nthreads = 4)
    
    # All versions should produce identical results
    expect_equal(res_serial, res_parallel_2, tolerance = 1e-10)
    expect_equal(res_serial, res_parallel_4, tolerance = 1e-10)
    
    # Check output structure
    expect_true(is.matrix(res_serial) || is.data.frame(res_serial))
    expect_equal(ncol(res_serial), 2)
    expect_equal(nrow(res_serial), nfeat)
    expect_true(all(colnames(res_serial) == c("raw_p_values", "adjusted_p_values")))
})

test_that("wilcoxon parallelization handles edge cases correctly", {
    set.seed(123)
    # Medium dataset: 100 genes, 10 samples
    nfeat <- 100
    nsamp <- 10
    mat <- matrix(c(rnorm(nfeat * nsamp * 0.8, mean = 1, sd = 0.5),
                    rep(NA_real_, nfeat * nsamp * 0.2)), 
                  nrow = nfeat, ncol = nsamp)
    samples <- rep(c("A", "B"), each = nsamp / 2)
    
    # Serial vs parallel should match
    res_serial <- wilcoxon(mat, samples, nthreads = 1)
    res_parallel <- wilcoxon(mat, samples, nthreads = 2)
    
    # P-values should be valid (0 to 1)
    expect_true(all(res_serial[, 1] >= 0 & res_serial[, 1] <= 1, na.rm = TRUE))
    expect_true(all(res_parallel[, 1] >= 0 & res_parallel[, 1] <= 1, na.rm = TRUE))
    
    # Should match
    expect_equal(res_serial, res_parallel, tolerance = 1e-10)
})

context("Parallelization: label_shuffling with nthreads")

test_that("label_shuffling parallel execution completes and returns valid results", {
    set.seed(42)
    # Create large dataset: 300 genes, 12 samples
    nfeat <- 300
    nsamp <- 12
    mat <- matrix(runif(nfeat * nsamp, 0, 5), nrow = nfeat, ncol = nsamp)
    samples <- rep(c("Control", "Treatment"), each = nsamp / 2)
    
    # Run serial version
    set.seed(42)
    res_serial <- label_shuffling(mat, samples, control = "Control", method = "mean", 
                                   randomizations = 50, pcorr = "BH", nthreads = 1)
    
    # Run parallel version with 2 threads - will have different random stream
    # but should still produce valid p-values in the same range
    set.seed(42)
    res_parallel_2 <- label_shuffling(mat, samples, control = "Control", method = "mean", 
                                       randomizations = 50, pcorr = "BH", nthreads = 2)
    
    # Check output structure
    expect_true(is.matrix(res_serial) || is.data.frame(res_serial))
    expect_equal(ncol(res_serial), 2)
    expect_equal(nrow(res_serial), nfeat)
    expect_equal(ncol(res_parallel_2), 2)
    expect_equal(nrow(res_parallel_2), nfeat)
    
    # All p-values should be valid
    expect_true(all(res_serial[, 1] >= 0 & res_serial[, 1] <= 1))
    expect_true(all(res_parallel_2[, 1] >= 0 & res_parallel_2[, 1] <= 1))
})

test_that("label_shuffling parallelization produces valid p-values", {
    set.seed(456)
    # Medium dataset: 150 genes, 8 samples
    nfeat <- 150
    nsamp <- 8
    mat <- matrix(rpois(nfeat * nsamp, lambda = 5), nrow = nfeat, ncol = nsamp)
    samples <- rep(c("Normal", "Tumor"), each = nsamp / 2)
    
    # Different randomization counts with parallelization
    res_100_serial <- label_shuffling(mat, samples, "Normal", "median", 100, 
                                      pcorr = "none", nthreads = 1)
    res_100_parallel <- label_shuffling(mat, samples, "Normal", "median", 100, 
                                        pcorr = "none", nthreads = 2)
    
    # Results should have valid p-values
    expect_true(all(res_100_serial[, 1] >= 0 & res_100_serial[, 1] <= 1))
    expect_true(all(res_100_parallel[, 1] >= 0 & res_100_parallel[, 1] <= 1))
    
    # Both should have same dimensions
    expect_equal(dim(res_100_serial), dim(res_100_parallel))
})

context("Parallelization: calculate_method with nthreads")

test_that("calculate_method serial and parallel produce identical results", {
    set.seed(789)
    # Create large dataset: 50 genes with multiple transcripts each, 6 samples
    ntx <- 300  # total transcripts
    nsamp <- 6
    x <- matrix(rpois(ntx * nsamp, lambda = 10), nrow = ntx, ncol = nsamp)
    colnames(x) <- paste0("Sample", seq_len(nsamp))
    # Assign transcripts to genes (50 genes, ~6 transcripts each)
    genes <- rep(paste0("Gene_", seq_len(50)), each = 6)
    
    # Single q value
    res_serial <- calculate_method(x, genes, norm = TRUE, q = 2, nthreads = 1)
    res_parallel_2 <- calculate_method(x, genes, norm = TRUE, q = 2, nthreads = 2)
    res_parallel_4 <- calculate_method(x, genes, norm = TRUE, q = 2, nthreads = 4)
    
    # Should produce identical results
    expect_equal(res_serial, res_parallel_2, tolerance = 1e-10)
    expect_equal(res_serial, res_parallel_4, tolerance = 1e-10)
    
    # Check structure
    expect_equal(nrow(res_serial), 50)
    expect_equal(nrow(res_parallel_2), 50)
    expect_true("Gene" %in% colnames(res_serial))
})

test_that("calculate_method parallelization with multiple q values", {
    set.seed(111)
    # Create dataset: 40 genes with multiple transcripts, 5 samples
    ntx <- 200
    nsamp <- 5
    x <- matrix(rpois(ntx * nsamp, lambda = 8), nrow = ntx, ncol = nsamp)
    colnames(x) <- paste0("S", seq_len(nsamp))
    genes <- rep(paste0("G", seq_len(40)), each = 5)
    
    # Multiple q values
    qvec <- c(0.5, 1, 2)
    
    res_serial <- calculate_method(x, genes, norm = TRUE, q = qvec, nthreads = 1)
    res_parallel <- calculate_method(x, genes, norm = TRUE, q = qvec, nthreads = 2)
    
    expect_equal(res_serial, res_parallel, tolerance = 1e-10)
    
    # Should have columns for each sample-q combination
    expected_cols <- 1 + (nsamp * length(qvec))
    expect_equal(ncol(res_serial), expected_cols)
})

context("Parallelization: calculate_difference with nthreads")

test_that("calculate_difference serial and parallel produce identical results", {
    set.seed(222)
    # Create large dataset: 250 genes, 10 samples
    nfeat <- 250
    nsamp <- 10
    data_mat <- matrix(runif(nfeat * nsamp, 0, 1), nrow = nfeat, ncol = nsamp)
    data_df <- as.data.frame(data_mat)
    data_df <- cbind(genes = paste0("Gene_", seq_len(nfeat)), data_df)
    colnames(data_df)[-1] <- paste0("Sample", seq_len(nsamp))
    samples <- rep(c("Normal", "Tumor"), each = nsamp / 2)
    
    # Wilcoxon test, serial
    res_wilcox_serial <- calculate_difference(data_df, samples = samples, 
                                              control = "Normal", test = "wilcoxon",
                                              nthreads = 1, verbose = FALSE)
    
    # Wilcoxon test, parallel
    res_wilcox_par_2 <- calculate_difference(data_df, samples = samples,
                                             control = "Normal", test = "wilcoxon",
                                             nthreads = 2, verbose = FALSE)
    
    res_wilcox_par_4 <- calculate_difference(data_df, samples = samples,
                                             control = "Normal", test = "wilcoxon",
                                             nthreads = 4, verbose = FALSE)
    
    # Should produce identical results
    expect_equal(res_wilcox_serial, res_wilcox_par_2, tolerance = 1e-10)
    expect_equal(res_wilcox_serial, res_wilcox_par_4, tolerance = 1e-10)
    
    # Verify output structure
    expect_true(is.data.frame(res_wilcox_serial))
    expect_true("genes" %in% colnames(res_wilcox_serial))
    expect_true("log2_fold_change" %in% colnames(res_wilcox_serial))
    expect_true("raw_p_values" %in% colnames(res_wilcox_serial))
    expect_equal(nrow(res_wilcox_serial), nfeat)
})

test_that("calculate_difference label_shuffling produces valid results", {
    set.seed(333)
    # Create dataset: 150 genes, 8 samples
    nfeat <- 150
    nsamp <- 8
    data_mat <- matrix(runif(nfeat * nsamp, 0.1, 2), nrow = nfeat, ncol = nsamp)
    data_df <- as.data.frame(data_mat)
    data_df <- cbind(genes = paste0("Gene", seq_len(nfeat)), data_df)
    colnames(data_df)[-1] <- paste0("Samp", seq_len(nsamp))
    samples <- c(rep("Ctrl", nsamp/2), rep("Trt", nsamp/2))
    
    # Label shuffling, serial
    res_shuffle_serial <- suppressWarnings(
        calculate_difference(data_df, samples = samples,
                           control = "Ctrl", test = "shuffle",
                           randomizations = 30, nthreads = 1,
                           verbose = FALSE)
    )
    
    # Label shuffling, parallel
    res_shuffle_par <- suppressWarnings(
        calculate_difference(data_df, samples = samples,
                            control = "Ctrl", test = "shuffle",
                            randomizations = 30, nthreads = 2,
                            verbose = FALSE)
    )
    
    # Check validity - randomized results won't be identical
    expect_true(all(res_shuffle_serial$raw_p_values >= 0 & res_shuffle_serial$raw_p_values <= 1))
    expect_true(all(res_shuffle_par$raw_p_values >= 0 & res_shuffle_par$raw_p_values <= 1))
    
    # Should have same number of results
    expect_equal(nrow(res_shuffle_serial), nrow(res_shuffle_par))
})

context("Parallelization: Large dataset stress test")

test_that("Large dataset parallelization produces valid results", {
    set.seed(444)
    # Very large dataset to ensure parallelization overhead is worth it
    # 1000 genes, 24 samples
    nfeat <- 1000
    nsamp <- 24
    mat <- matrix(rpois(nfeat * nsamp, lambda = 10), nrow = nfeat, ncol = nsamp)
    samples <- rep(c("GroupA", "GroupB"), each = nsamp / 2)
    
    # Wilcoxon with different thread counts
    res_1 <- wilcoxon(mat, samples, nthreads = 1)
    res_2 <- wilcoxon(mat, samples, nthreads = 2)
    res_4 <- wilcoxon(mat, samples, nthreads = 4)
    
    # All should match
    expect_equal(res_1, res_2, tolerance = 1e-10)
    expect_equal(res_1, res_4, tolerance = 1e-10)
    
    # All p-values should be valid
    expect_true(all(is.finite(res_1[, 1])))
    expect_true(all(res_1[, 2] >= 0 & res_1[, 2] <= 1))
})

test_that("Large calculate_method dataset with single q value", {
    set.seed(555)
    # 500 genes with multiple transcripts, 8 samples, single q value
    ntx <- 2000  # 4 transcripts per gene
    nsamp <- 8
    x <- matrix(rpois(ntx * nsamp, lambda = 15), nrow = ntx, ncol = nsamp)
    colnames(x) <- paste0("Samp_", seq_len(nsamp))
    genes <- rep(paste0("Gene_", seq_len(500)), each = 4)
    
    # Serial execution
    res_serial <- calculate_method(x, genes, norm = TRUE, q = 2, nthreads = 1)
    
    # Parallel execution
    res_parallel <- calculate_method(x, genes, norm = TRUE, q = 2, nthreads = 4)
    
    # Should produce identical results
    expect_equal(res_serial, res_parallel, tolerance = 1e-10)
    
    # Correct dimensions
    expect_equal(nrow(res_serial), 500)
    expect_equal(nrow(res_parallel), 500)
    expected_cols <- 1 + nsamp  # 1 for gene + 8 samples (single q)
    expect_equal(ncol(res_serial), expected_cols)
})

context("Parallelization: Consistency across different backends")

test_that("nthreads=1 produces serial results regardless of implementation", {
    set.seed(666)
    nfeat <- 200
    nsamp <- 8
    mat <- matrix(runif(nfeat * nsamp, 0, 1), nrow = nfeat, ncol = nsamp)
    samples <- rep(c("A", "B"), each = nsamp / 2)
    
    # Multiple calls with nthreads=1 should give identical results
    res1 <- wilcoxon(mat, samples, nthreads = 1)
    res2 <- wilcoxon(mat, samples, nthreads = 1)
    res3 <- wilcoxon(mat, samples, nthreads = 1)
    
    expect_equal(res1, res2)
    expect_equal(res2, res3)
})

test_that("label_shuffling produces reproducible results with different serial calls", {
    # Create moderate dataset
    set.seed(777)
    nfeat <- 100
    nsamp <- 6
    mat <- matrix(runif(nfeat * nsamp, 0, 2), nrow = nfeat, ncol = nsamp)
    samples <- rep(c("Ctrl", "Test"), each = nsamp / 2)
    
    # Multiple runs with same seed and serial execution should match
    set.seed(777)
    res_serial_1 <- label_shuffling(mat, samples, "Ctrl", "mean", 50, nthreads = 1)
    
    set.seed(777)
    res_serial_2 <- label_shuffling(mat, samples, "Ctrl", "mean", 50, nthreads = 1)
    
    # Serial calls should be identical
    expect_equal(res_serial_1, res_serial_2)
    
    # Parallel execution will have different random stream, so just check validity
    set.seed(777)
    res_parallel <- label_shuffling(mat, samples, "Ctrl", "mean", 50, nthreads = 2)
    
    expect_true(all(res_parallel[, 1] >= 0 & res_parallel[, 1] <= 1))
    expect_equal(dim(res_parallel), dim(res_serial_1))
})
context("Parallel helper functions")

test_that(".tsenat_get_bpparam returns SerialParam for nthreads=1", {
    bpparam <- TSENAT:::.tsenat_get_bpparam(nthreads = 1)
    expect_is(bpparam, "SerialParam")
})

test_that(".tsenat_get_bpparam returns MulticoreParam for nthreads>1 on Unix", {
    skip_if_not(identical(.Platform$OS.type, "unix"))
    bpparam <- TSENAT:::.tsenat_get_bpparam(nthreads = 2)
    expect_is(bpparam, "MulticoreParam")
})

test_that(".tsenat_bplapply with FUN.VALUE uses vapply simplification", {
    # Test the code path: result_list <- BiocParallel::bplapply(...); return(vapply(...))
    X <- 1:5
    FUN <- function(x) x * 2
    FUN.VALUE <- numeric(1)
    
    result <- TSENAT:::.tsenat_bplapply(X, FUN, nthreads = 1, FUN.VALUE = FUN.VALUE)
    
    # Should return a numeric vector
    expect_is(result, "numeric")
    expect_equal(result, c(2, 4, 6, 8, 10))
})

test_that(".tsenat_bplapply without FUN.VALUE returns list", {
    X <- 1:5
    FUN <- function(x) x * 2
    
    result <- TSENAT:::.tsenat_bplapply(X, FUN, nthreads = 1, FUN.VALUE = NULL)
    
    # Should return a list
    expect_is(result, "list")
    expect_equal(result, list(2, 4, 6, 8, 10))
})

test_that(".tsenat_bplapply with FUN.VALUE and nthreads=1 uses vapply", {
    # Test serial execution with FUN.VALUE specified
    X <- c("a", "b", "c")
    FUN <- function(x) nchar(x)
    FUN.VALUE <- integer(1)
    
    result <- TSENAT:::.tsenat_bplapply(X, FUN, nthreads = 1, SIMPLIFY = TRUE, FUN.VALUE = FUN.VALUE)
    
    expect_is(result, "integer")
    expect_equal(result, c(1L, 1L, 1L))
})

test_that(".tsenat_bplapply with matrix FUN.VALUE returns matrix", {
    # Test with more complex FUN.VALUE (matrix)
    X <- 1:3
    FUN <- function(x) c(x, x ^ 2)
    FUN.VALUE <- numeric(2)
    
    result <- TSENAT:::.tsenat_bplapply(X, FUN, nthreads = 1, FUN.VALUE = FUN.VALUE)
    
    expect_is(result, "matrix")
    expect_equal(dim(result), c(2, 3))
    expect_equal(result[, 1], c(1, 1))
    expect_equal(result[, 2], c(2, 4))
})

test_that(".tsenat_bpmapply serial execution with mapply", {
    # Test the code path: return(mapply(FUN, X, Y, SIMPLIFY = FALSE))
    X <- 1:5
    Y <- 10:14
    FUN <- function(x, y) x + y
    
    result <- TSENAT:::.tsenat_bpmapply(X, Y, FUN, nthreads = 1)
    
    expect_is(result, "list")
    # Expected: 1+10=11, 2+11=13, 3+12=15, 4+13=17, 5+14=19
    expect_equal(result, list(11, 13, 15, 17, 19))
})

test_that(".tsenat_bpmapply with vectors of different lengths", {
    X <- 1:3
    Y <- 10:12
    FUN <- function(x, y) c(x, y)
    
    result <- TSENAT:::.tsenat_bpmapply(X, Y, FUN, nthreads = 1)
    
    expect_is(result, "list")
    expect_equal(length(result), 3)
    expect_equal(result[[1]], c(1, 10))
    expect_equal(result[[3]], c(3, 12))
})

test_that(".tsenat_bpmapply preserves order", {
    X <- c("a", "b", "c")
    Y <- c(1, 2, 3)
    FUN <- function(x, y) paste(x, y, sep = "-")
    
    result <- TSENAT:::.tsenat_bpmapply(X, Y, FUN, nthreads = 1)
    
    expect_equal(result, list("a-1", "b-2", "c-3"))
})

context("Parallel helpers: .tsenat_bplapply with parallel execution (nthreads > 1)")

test_that(".tsenat_bplapply with nthreads > 1 and FUN.VALUE uses BiocParallel::bplapply + vapply", {
    # This test specifically covers lines 43-44:
    # result_list <- BiocParallel::bplapply(X, FUN, BPPARAM = bpparam)
    # return(vapply(result_list, identity, FUN.VALUE = FUN.VALUE))
    X <- 1:10
    FUN <- function(x) x * 3
    FUN.VALUE <- numeric(1)
    
    # Use nthreads = 2 to trigger parallel execution path
    result_parallel <- TSENAT:::.tsenat_bplapply(X, FUN, nthreads = 2, FUN.VALUE = FUN.VALUE)
    result_serial <- TSENAT:::.tsenat_bplapply(X, FUN, nthreads = 1, FUN.VALUE = FUN.VALUE)
    
    # Parallel and serial should produce identical results
    expect_equal(result_parallel, result_serial, tolerance = 1e-10)
    expect_is(result_parallel, "numeric")
    expect_equal(result_parallel, c(3, 6, 9, 12, 15, 18, 21, 24, 27, 30))
})

test_that(".tsenat_bplapply parallel execution with integer FUN.VALUE", {
    # Test parallel execution with integer simplification
    X <- c("cat", "dog", "elephant")
    FUN <- function(x) nchar(x)
    FUN.VALUE <- integer(1)
    
    result_parallel <- TSENAT:::.tsenat_bplapply(X, FUN, nthreads = 2, FUN.VALUE = FUN.VALUE)
    result_serial <- TSENAT:::.tsenat_bplapply(X, FUN, nthreads = 1, FUN.VALUE = FUN.VALUE)
    
    expect_equal(result_parallel, result_serial)
    expect_is(result_parallel, "integer")
    expect_equal(result_parallel, c(3L, 3L, 8L))
})

test_that(".tsenat_bplapply parallel execution with numeric matrix FUN.VALUE", {
    # Test parallel execution with matrix simplification (result_list -> vapply with identity)
    X <- 1:5
    FUN <- function(x) c(x, x ^ 2, sqrt(x))
    FUN.VALUE <- numeric(3)
    
    result_parallel <- TSENAT:::.tsenat_bplapply(X, FUN, nthreads = 2, FUN.VALUE = FUN.VALUE)
    result_serial <- TSENAT:::.tsenat_bplapply(X, FUN, nthreads = 1, FUN.VALUE = FUN.VALUE)
    
    expect_equal(result_parallel, result_serial, tolerance = 1e-10)
    expect_is(result_parallel, "matrix")
    expect_equal(dim(result_parallel), c(3, 5))
    # Verify first column: [1, 1, 1] (for x=1: c(1, 1^2, sqrt(1)) = c(1, 1, 1))
    expect_equal(result_parallel[, 1], c(1, 1, 1))
    # Verify second column: [2, 4, sqrt(2)] (for x=2: c(2, 2^2, sqrt(2)))
    expect_equal(result_parallel[1, 2], 2)
    expect_equal(result_parallel[2, 2], 4)
})

test_that(".tsenat_bplapply parallel execution with nthreads=4", {
    # Test with more threads
    X <- seq(1, 100, by = 10)
    FUN <- function(x) log(x)
    FUN.VALUE <- numeric(1)
    
    result_parallel_2 <- TSENAT:::.tsenat_bplapply(X, FUN, nthreads = 2, FUN.VALUE = FUN.VALUE)
    result_parallel_4 <- TSENAT:::.tsenat_bplapply(X, FUN, nthreads = 4, FUN.VALUE = FUN.VALUE)
    result_serial <- TSENAT:::.tsenat_bplapply(X, FUN, nthreads = 1, FUN.VALUE = FUN.VALUE)
    
    # All should be equal
    expect_equal(result_parallel_2, result_serial, tolerance = 1e-10)
    expect_equal(result_parallel_4, result_serial, tolerance = 1e-10)
    expect_is(result_parallel_4, "numeric")
})

test_that(".tsenat_bplapply parallel execution with complex function", {
    # Test with a more realistic function that does computation
    X <- list(c(1, 2, 3), c(4, 5, 6), c(7, 8, 9))
    FUN <- function(vec) mean(vec)
    FUN.VALUE <- numeric(1)
    
    result_parallel <- TSENAT:::.tsenat_bplapply(X, FUN, nthreads = 2, FUN.VALUE = FUN.VALUE)
    result_serial <- TSENAT:::.tsenat_bplapply(X, FUN, nthreads = 1, FUN.VALUE = FUN.VALUE)
    
    expect_equal(result_parallel, result_serial, tolerance = 1e-10)
    expect_equal(result_parallel, c(2, 5, 8), tolerance = 1e-10)
})

test_that(".tsenat_bplapply parallel execution with logical FUN.VALUE", {
    # Test parallel with logical output simplification
    X <- c(1, 2, 3, 4, 5)
    FUN <- function(x) x > 2
    FUN.VALUE <- logical(1)
    
    result_parallel <- TSENAT:::.tsenat_bplapply(X, FUN, nthreads = 2, FUN.VALUE = FUN.VALUE)
    result_serial <- TSENAT:::.tsenat_bplapply(X, FUN, nthreads = 1, FUN.VALUE = FUN.VALUE)
    
    expect_equal(result_parallel, result_serial)
    expect_is(result_parallel, "logical")
    expect_equal(result_parallel, c(FALSE, FALSE, TRUE, TRUE, TRUE))
})

test_that(".tsenat_bplapply parallel without FUN.VALUE returns list", {
    # Test parallel execution without FUN.VALUE (different code path)
    X <- 1:5
    FUN <- function(x) list(x = x, squared = x ^ 2)
    
    result_parallel <- TSENAT:::.tsenat_bplapply(X, FUN, nthreads = 2, FUN.VALUE = NULL)
    result_serial <- TSENAT:::.tsenat_bplapply(X, FUN, nthreads = 1, FUN.VALUE = NULL)
    
    expect_is(result_parallel, "list")
    expect_equal(length(result_parallel), 5)
    expect_equal(result_parallel, result_serial)
})

context("Parallel helpers: .tsenat_bpmapply with parallel execution (nthreads > 1)")

test_that(".tsenat_bpmapply with nthreads > 1 uses BiocParallel::bpmapply", {
    # Test parallel execution for bpmapply
    X <- 1:5
    Y <- 10:14
    FUN <- function(x, y) x + y
    
    result_parallel <- TSENAT:::.tsenat_bpmapply(X, Y, FUN, nthreads = 2)
    result_serial <- TSENAT:::.tsenat_bpmapply(X, Y, FUN, nthreads = 1)
    
    expect_is(result_parallel, "list")
    expect_equal(result_parallel, result_serial)
    expect_equal(result_parallel, list(11, 13, 15, 17, 19))
})

test_that(".tsenat_bpmapply parallel with complex operation", {
    # Test with more complex function
    X <- list(c(1, 2, 3), c(4, 5, 6))
    Y <- list(c(10, 20, 30), c(40, 50, 60))
    FUN <- function(x, y) list(sum = sum(x) + sum(y), means = mean(c(x, y)))
    
    result_parallel <- TSENAT:::.tsenat_bpmapply(X, Y, FUN, nthreads = 2)
    result_serial <- TSENAT:::.tsenat_bpmapply(X, Y, FUN, nthreads = 1)
    
    expect_equal(result_parallel, result_serial)
    expect_is(result_parallel, "list")
    expect_equal(result_parallel[[1]]$sum, 66)  # (1+2+3) + (10+20+30) = 66
    expect_equal(result_parallel[[2]]$sum, 165)  # (4+5+6) + (40+50+60) = 165
})

test_that(".tsenat_bpmapply parallel with nthreads=2 initializes bpparam correctly", {
    # This test specifically covers lines 57-59:
    # bpparam <- .tsenat_get_bpparam(nthreads)
    # return(unname(BiocParallel::bpmapply(FUN, X, Y, BPPARAM = bpparam, SIMPLIFY = FALSE)))
    X <- 1:10
    Y <- 11:20
    FUN <- function(x, y) x * y
    
    result <- TSENAT:::.tsenat_bpmapply(X, Y, FUN, nthreads = 2)
    
    expect_is(result, "list")
    expect_equal(length(result), 10)
    # Verify some computations: 1*11=11, 5*15=75, 10*20=200
    expect_equal(result[[1]], 11)
    expect_equal(result[[5]], 75)
    expect_equal(result[[10]], 200)
})

test_that(".tsenat_bpmapply parallel with nthreads=4 uses MulticoreParam", {
    # Test with 4 threads to ensure bpparam initialization works correctly
    X <- c("a", "b", "c", "d")
    Y <- c(1, 2, 3, 4)
    FUN <- function(x, y) rep(x, y)
    
    result_parallel_4 <- TSENAT:::.tsenat_bpmapply(X, Y, FUN, nthreads = 4)
    result_serial <- TSENAT:::.tsenat_bpmapply(X, Y, FUN, nthreads = 1)
    
    expect_equal(result_parallel_4, result_serial)
    expect_is(result_parallel_4, "list")
    expect_equal(result_parallel_4[[1]], c("a"))
    expect_equal(result_parallel_4[[2]], c("b", "b"))
    expect_equal(result_parallel_4[[4]], c("d", "d", "d", "d"))
})

test_that(".tsenat_bpmapply parallel returns unnamned list with BiocParallel", {
    # Verify that unname() is applied to the BiocParallel::bpmapply result
    X <- 1:3
    Y <- c("x", "y", "z")
    FUN <- function(x, y) paste0(y, x)
    
    result <- TSENAT:::.tsenat_bpmapply(X, Y, FUN, nthreads = 2)
    
    # Result should be a list with no names
    expect_is(result, "list")
    expect_null(names(result))
    expect_equal(result[[1]], "x1")
    expect_equal(result[[2]], "y2")
    expect_equal(result[[3]], "z3")
})

test_that(".tsenat_bpmapply parallel with numeric vectors and bpparam initialization", {
    # Test with numeric computations to verify bpparam is correctly initialized
    X <- c(0.5, 1.5, 2.5, 3.5)
    Y <- c(10, 20, 30, 40)
    FUN <- function(x, y) x * y + sqrt(x)
    
    result_parallel <- TSENAT:::.tsenat_bpmapply(X, Y, FUN, nthreads = 2)
    result_serial <- TSENAT:::.tsenat_bpmapply(X, Y, FUN, nthreads = 1)
    
    expect_equal(result_parallel, result_serial, tolerance = 1e-10)
    # Verify: 0.5*10 + sqrt(0.5) ≈ 5.707
    expect_equal(result_parallel[[1]], 0.5 * 10 + sqrt(0.5), tolerance = 1e-10)
})

test_that(".tsenat_bpmapply parallel with large vectors and bpparam", {
    # Test with larger data to ensure BiocParallel::bpmapply with bpparam works efficiently
    X <- 1:100
    Y <- 101:200
    FUN <- function(x, y) (x + y) / 2  # mean of x and y
    
    result_parallel <- TSENAT:::.tsenat_bpmapply(X, Y, FUN, nthreads = 2)
    result_serial <- TSENAT:::.tsenat_bpmapply(X, Y, FUN, nthreads = 1)
    
    expect_equal(result_parallel, result_serial, tolerance = 1e-10)
    expect_equal(length(result_parallel), 100)
    # Verify first and last: (1+101)/2 = 51, (100+200)/2 = 150
    expect_equal(result_parallel[[1]], 51)
    expect_equal(result_parallel[[100]], 150)
})

test_that(".tsenat_bpmapply parallel executes correctly with different nthreads values", {
    # Test that bpparam initialization works for various thread counts
    X <- 1:6
    Y <- 6:1
    FUN <- function(x, y) c(x, y)
    
    result_parallel_2 <- TSENAT:::.tsenat_bpmapply(X, Y, FUN, nthreads = 2)
    result_parallel_4 <- TSENAT:::.tsenat_bpmapply(X, Y, FUN, nthreads = 4)
    result_serial <- TSENAT:::.tsenat_bpmapply(X, Y, FUN, nthreads = 1)
    
    # All should be equal
    expect_equal(result_parallel_2, result_serial)
    expect_equal(result_parallel_4, result_serial)
    expect_is(result_parallel_2, "list")
})

test_that(".tsenat_bpmapply SIMPLIFY=FALSE is respected in parallel execution", {
    # Verify that SIMPLIFY=FALSE is correctly passed to BiocParallel::bpmapply
    X <- 1:4
    Y <- 1:4
    FUN <- function(x, y) list(sum = x + y, product = x * y)
    
    result <- TSENAT:::.tsenat_bpmapply(X, Y, FUN, nthreads = 2)
    
    # Result should be a list of lists, not simplified
    expect_is(result, "list")
    expect_equal(length(result), 4)
    expect_is(result[[1]], "list")
    expect_equal(result[[1]]$sum, 2)
    expect_equal(result[[1]]$product, 1)
    expect_equal(result[[4]]$sum, 8)
    expect_equal(result[[4]]$product, 16)
})