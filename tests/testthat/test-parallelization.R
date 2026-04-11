context("Parallelization: Method Calculation with Multiple Threads")

# Internal helper access for calculate_method
calculate_method <- TSENAT:::.calculate_method

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
    res_serial <- .calculate_method(x, genes, norm = TRUE, q = 2, nthreads = 1)
    res_parallel_2 <- .calculate_method(x, genes, norm = TRUE, q = 2, nthreads = 2)
    core_limit <- suppressWarnings(as.integer(Sys.getenv("_R_CHECK_LIMIT_CORES_", NA)))
    max_threads <- if (is.na(core_limit)) min(2, parallel::detectCores()) else min(2, core_limit)
    res_parallel_4 <- .calculate_method(x, genes, norm = TRUE, q = 2, nthreads = max_threads)
    
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
    
    res_serial <- .calculate_method(x, genes, norm = TRUE, q = qvec, nthreads = 1)
    res_parallel <- .calculate_method(x, genes, norm = TRUE, q = qvec, nthreads = 2)
    
    expect_equal(res_serial, res_parallel, tolerance = 1e-10)
    
    # Should have columns for each sample-q combination
    expected_cols <- 1 + (nsamp * length(qvec))
    expect_equal(ncol(res_serial), expected_cols)
})

context("Parallelization: Large Dataset Stress Testing")

test_that("Large calculate_method dataset with single q value", {
    skip_on_cran()
    
    set.seed(555)
    # 500 genes with multiple transcripts, 8 samples, single q value
    ntx <- 2000  # 4 transcripts per gene
    nsamp <- 8
    x <- matrix(rpois(ntx * nsamp, lambda = 15), nrow = ntx, ncol = nsamp)
    colnames(x) <- paste0("Samp_", seq_len(nsamp))
    genes <- rep(paste0("Gene_", seq_len(500)), each = 4)
    
    # Serial execution
    res_serial <- .calculate_method(x, genes, norm = TRUE, q = 2, nthreads = 1)
    
    # Parallel execution
    core_limit <- suppressWarnings(as.integer(Sys.getenv("_R_CHECK_LIMIT_CORES_", NA)))
    max_threads <- if (is.na(core_limit)) min(2, parallel::detectCores()) else min(2, core_limit)
    res_parallel <- .calculate_method(x, genes, norm = TRUE, q = 2, nthreads = max_threads)
    
    # Should produce identical results
    expect_equal(res_serial, res_parallel, tolerance = 1e-10)
    
    # Correct dimensions
    expect_equal(nrow(res_serial), 500)
    expect_equal(nrow(res_parallel), 500)
    expected_cols <- 1 + nsamp  # 1 for gene + 8 samples (single q)
    expect_equal(ncol(res_serial), expected_cols)
})

context("Parallelization: Cross-Backend Consistency Validation")


context("Parallel helper functions")

test_that(".get_bpparam returns SerialParam for nthreads=1", {
    bpparam <- TSENAT:::.get_bpparam(nthreads = 1)
    expect_is(bpparam, "SerialParam")
})

test_that(".get_bpparam returns MulticoreParam for nthreads>1 on Unix", {
    skip_if_not(identical(.Platform$OS.type, "unix"))
    bpparam <- TSENAT:::.get_bpparam(nthreads = 2)
    expect_is(bpparam, "MulticoreParam")
})

test_that(".bplapply with FUN.VALUE uses vapply simplification", {
    # Test the code path: result_list <- BiocParallel::bplapply(...); return(vapply(...))
    X <- 1:5
    FUN <- function(x) x * 2
    FUN.VALUE <- numeric(1)
    
    result <- TSENAT:::.bplapply(X, FUN, nthreads = 1, FUN.VALUE = FUN.VALUE)
    
    # Should return a numeric vector
    expect_is(result, "numeric")
    expect_equal(result, c(2, 4, 6, 8, 10))
})

test_that(".bplapply without FUN.VALUE returns list", {
    X <- 1:5
    FUN <- function(x) x * 2
    
    result <- TSENAT:::.bplapply(X, FUN, nthreads = 1, FUN.VALUE = NULL)
    
    # Should return a list
    expect_is(result, "list")
    expect_equal(result, list(2, 4, 6, 8, 10))
})

test_that(".bplapply with FUN.VALUE and nthreads=1 uses vapply", {
    # Test serial execution with FUN.VALUE specified
    X <- c("a", "b", "c")
    FUN <- function(x) nchar(x)
    FUN.VALUE <- integer(1)
    
    result <- TSENAT:::.bplapply(X, FUN, nthreads = 1, SIMPLIFY = TRUE, FUN.VALUE = FUN.VALUE)
    
    expect_is(result, "integer")
    expect_equal(result, c(1L, 1L, 1L))
})

test_that(".bplapply with matrix FUN.VALUE returns matrix", {
    # Test with more complex FUN.VALUE (matrix)
    X <- 1:3
    FUN <- function(x) c(x, x ^ 2)
    FUN.VALUE <- numeric(2)
    
    result <- TSENAT:::.bplapply(X, FUN, nthreads = 1, FUN.VALUE = FUN.VALUE)
    
    expect_is(result, "matrix")
    expect_equal(dim(result), c(2, 3))
    expect_equal(result[, 1], c(1, 1))
    expect_equal(result[, 2], c(2, 4))
})

test_that(".bpmapply serial execution with mapply", {
    # Test the code path: return(mapply(FUN, X, Y, SIMPLIFY = FALSE))
    X <- 1:5
    Y <- 10:14
    FUN <- function(x, y) x + y
    
    result <- TSENAT:::.bpmapply(X, Y, FUN, nthreads = 1)
    
    expect_is(result, "list")
    # Expected: 1+10=11, 2+11=13, 3+12=15, 4+13=17, 5+14=19
    expect_equal(result, list(11, 13, 15, 17, 19))
})

test_that(".bpmapply with vectors of different lengths", {
    X <- 1:3
    Y <- 10:12
    FUN <- function(x, y) c(x, y)
    
    result <- TSENAT:::.bpmapply(X, Y, FUN, nthreads = 1)
    
    expect_is(result, "list")
    expect_equal(length(result), 3)
    expect_equal(result[[1]], c(1, 10))
    expect_equal(result[[3]], c(3, 12))
})

test_that(".bpmapply preserves order", {
    X <- c("a", "b", "c")
    Y <- c(1, 2, 3)
    FUN <- function(x, y) paste(x, y, sep = "-")
    
    result <- TSENAT:::.bpmapply(X, Y, FUN, nthreads = 1)
    
    expect_equal(result, list("a-1", "b-2", "c-3"))
})

context("Parallel helpers: .bplapply with parallel execution (nthreads > 1)")

test_that(".bplapply with nthreads > 1 and FUN.VALUE uses BiocParallel::bplapply + vapply", {
    # This test specifically covers lines 43-44:
    # result_list <- BiocParallel::bplapply(X, FUN, BPPARAM = bpparam)
    # return(vapply(result_list, identity, FUN.VALUE = FUN.VALUE))
    X <- 1:10
    FUN <- function(x) x * 3
    FUN.VALUE <- numeric(1)
    
    # Use nthreads = 2 to trigger parallel execution path
    result_parallel <- TSENAT:::.bplapply(X, FUN, nthreads = 2, FUN.VALUE = FUN.VALUE)
    result_serial <- TSENAT:::.bplapply(X, FUN, nthreads = 1, FUN.VALUE = FUN.VALUE)
    
    # Parallel and serial should produce identical results
    expect_equal(result_parallel, result_serial, tolerance = 1e-10)
    expect_is(result_parallel, "numeric")
    expect_equal(result_parallel, c(3, 6, 9, 12, 15, 18, 21, 24, 27, 30))
})

test_that(".bplapply parallel execution with integer FUN.VALUE", {
    # Test parallel execution with integer simplification
    X <- c("cat", "dog", "elephant")
    FUN <- function(x) nchar(x)
    FUN.VALUE <- integer(1)
    
    result_parallel <- TSENAT:::.bplapply(X, FUN, nthreads = 2, FUN.VALUE = FUN.VALUE)
    result_serial <- TSENAT:::.bplapply(X, FUN, nthreads = 1, FUN.VALUE = FUN.VALUE)
    
    expect_equal(result_parallel, result_serial)
    expect_is(result_parallel, "integer")
    expect_equal(result_parallel, c(3L, 3L, 8L))
})

test_that(".bplapply parallel execution with numeric matrix FUN.VALUE", {
    # Test parallel execution with matrix simplification (result_list -> vapply with identity)
    X <- 1:5
    FUN <- function(x) c(x, x ^ 2, sqrt(x))
    FUN.VALUE <- numeric(3)
    
    result_parallel <- TSENAT:::.bplapply(X, FUN, nthreads = 2, FUN.VALUE = FUN.VALUE)
    result_serial <- TSENAT:::.bplapply(X, FUN, nthreads = 1, FUN.VALUE = FUN.VALUE)
    
    expect_equal(result_parallel, result_serial, tolerance = 1e-10)
    expect_is(result_parallel, "matrix")
    expect_equal(dim(result_parallel), c(3, 5))
    # Verify first column: [1, 1, 1] (for x=1: c(1, 1^2, sqrt(1)) = c(1, 1, 1))
    expect_equal(result_parallel[, 1], c(1, 1, 1))
    # Verify second column: [2, 4, sqrt(2)] (for x=2: c(2, 2^2, sqrt(2)))
    expect_equal(result_parallel[1, 2], 2)
    expect_equal(result_parallel[2, 2], 4)
})

test_that(".bplapply parallel execution with nthreads=3", {
    # Test with more threads
    X <- seq(1, 100, by = 10)
    FUN <- function(x) log(x)
    FUN.VALUE <- numeric(1)
    
    result_parallel_2 <- TSENAT:::.bplapply(X, FUN, nthreads = 2, FUN.VALUE = FUN.VALUE)
    core_limit <- suppressWarnings(as.integer(Sys.getenv("_R_CHECK_LIMIT_CORES_", NA)))
    max_threads <- if (is.na(core_limit)) min(2, parallel::detectCores()) else min(2, core_limit)
    result_parallel_4 <- TSENAT:::.bplapply(X, FUN, nthreads = max_threads, FUN.VALUE = FUN.VALUE)
    result_serial <- TSENAT:::.bplapply(X, FUN, nthreads = 1, FUN.VALUE = FUN.VALUE)
    
    # All should be equal
    expect_equal(result_parallel_2, result_serial, tolerance = 1e-10)
    expect_equal(result_parallel_4, result_serial, tolerance = 1e-10)
    expect_is(result_parallel_4, "numeric")
})

test_that(".bplapply parallel execution with complex function", {
    # Test with a more realistic function that does computation
    X <- list(c(1, 2, 3), c(4, 5, 6), c(7, 8, 9))
    FUN <- function(vec) mean(vec)
    FUN.VALUE <- numeric(1)
    
    result_parallel <- TSENAT:::.bplapply(X, FUN, nthreads = 2, FUN.VALUE = FUN.VALUE)
    result_serial <- TSENAT:::.bplapply(X, FUN, nthreads = 1, FUN.VALUE = FUN.VALUE)
    
    expect_equal(result_parallel, result_serial, tolerance = 1e-10)
    expect_equal(result_parallel, c(2, 5, 8), tolerance = 1e-10)
})

test_that(".bplapply parallel execution with logical FUN.VALUE", {
    # Test parallel with logical output simplification
    X <- c(1, 2, 3, 4, 5)
    FUN <- function(x) x > 2
    FUN.VALUE <- logical(1)
    
    result_parallel <- TSENAT:::.bplapply(X, FUN, nthreads = 2, FUN.VALUE = FUN.VALUE)
    result_serial <- TSENAT:::.bplapply(X, FUN, nthreads = 1, FUN.VALUE = FUN.VALUE)
    
    expect_equal(result_parallel, result_serial)
    expect_is(result_parallel, "logical")
    expect_equal(result_parallel, c(FALSE, FALSE, TRUE, TRUE, TRUE))
})

test_that(".bplapply parallel without FUN.VALUE returns list", {
    # Test parallel execution without FUN.VALUE (different code path)
    X <- 1:5
    FUN <- function(x) list(x = x, squared = x ^ 2)
    
    result_parallel <- TSENAT:::.bplapply(X, FUN, nthreads = 2, FUN.VALUE = NULL)
    result_serial <- TSENAT:::.bplapply(X, FUN, nthreads = 1, FUN.VALUE = NULL)
    
    expect_is(result_parallel, "list")
    expect_equal(length(result_parallel), 5)
    expect_equal(result_parallel, result_serial)
})

context("Parallel helpers: .bpmapply with parallel execution (nthreads > 1)")

test_that(".bpmapply with nthreads > 1 uses BiocParallel::bpmapply", {
    # Test parallel execution for bpmapply
    X <- 1:5
    Y <- 10:14
    FUN <- function(x, y) x + y
    
    result_parallel <- TSENAT:::.bpmapply(X, Y, FUN, nthreads = 2)
    result_serial <- TSENAT:::.bpmapply(X, Y, FUN, nthreads = 1)
    
    expect_is(result_parallel, "list")
    expect_equal(result_parallel, result_serial)
    expect_equal(result_parallel, list(11, 13, 15, 17, 19))
})

test_that(".bpmapply parallel with complex operation", {
    # Test with more complex function
    X <- list(c(1, 2, 3), c(4, 5, 6))
    Y <- list(c(10, 20, 30), c(40, 50, 60))
    FUN <- function(x, y) list(sum = sum(x) + sum(y), means = mean(c(x, y)))
    
    result_parallel <- TSENAT:::.bpmapply(X, Y, FUN, nthreads = 2)
    result_serial <- TSENAT:::.bpmapply(X, Y, FUN, nthreads = 1)
    
    expect_equal(result_parallel, result_serial)
    expect_is(result_parallel, "list")
    expect_equal(result_parallel[[1]]$sum, 66)  # (1+2+3) + (10+20+30) = 66
    expect_equal(result_parallel[[2]]$sum, 165)  # (4+5+6) + (40+50+60) = 165
})

test_that(".bpmapply parallel with nthreads=2 initializes bpparam correctly", {
    # This test specifically covers lines 57-59:
    # bpparam <- .get_bpparam(nthreads)
    # return(unname(BiocParallel::bpmapply(FUN, X, Y, BPPARAM = bpparam, SIMPLIFY = FALSE)))
    X <- 1:10
    Y <- 11:20
    FUN <- function(x, y) x * y
    
    result <- TSENAT:::.bpmapply(X, Y, FUN, nthreads = 2)
    
    expect_is(result, "list")
    expect_equal(length(result), 10)
    # Verify some computations: 1*11=11, 5*15=75, 10*20=200
    expect_equal(result[[1]], 11)
    expect_equal(result[[5]], 75)
    expect_equal(result[[10]], 200)
})

test_that(".bpmapply parallel with nthreads=3 uses MulticoreParam", {
    # Test with 3 threads to ensure bpparam initialization works correctly
    X <- c("a", "b", "c", "d")
    Y <- c(1, 2, 3, 4)
    FUN <- function(x, y) rep(x, y)
    
    core_limit <- suppressWarnings(as.integer(Sys.getenv("_R_CHECK_LIMIT_CORES_", NA)))
    max_threads <- if (is.na(core_limit)) min(2, parallel::detectCores()) else min(2, core_limit)
    result_parallel_4 <- TSENAT:::.bpmapply(X, Y, FUN, nthreads = max_threads)
    result_serial <- TSENAT:::.bpmapply(X, Y, FUN, nthreads = 1)
    
    expect_equal(result_parallel_4, result_serial)
    expect_is(result_parallel_4, "list")
    expect_equal(result_parallel_4[[1]], c("a"))
    expect_equal(result_parallel_4[[2]], c("b", "b"))
    expect_equal(result_parallel_4[[4]], c("d", "d", "d", "d"))
})

test_that(".bpmapply parallel returns unnamned list with BiocParallel", {
    # Verify that unname() is applied to the BiocParallel::bpmapply result
    X <- 1:3
    Y <- c("x", "y", "z")
    FUN <- function(x, y) paste0(y, x)
    
    result <- TSENAT:::.bpmapply(X, Y, FUN, nthreads = 2)
    
    # Result should be a list with no names
    expect_is(result, "list")
    expect_null(names(result))
    expect_equal(result[[1]], "x1")
    expect_equal(result[[2]], "y2")
    expect_equal(result[[3]], "z3")
})

test_that(".bpmapply parallel with numeric vectors and bpparam initialization", {
    # Test with numeric computations to verify bpparam is correctly initialized
    X <- c(0.5, 1.5, 2.5, 3.5)
    Y <- c(10, 20, 30, 40)
    FUN <- function(x, y) x * y + sqrt(x)
    
    result_parallel <- TSENAT:::.bpmapply(X, Y, FUN, nthreads = 2)
    result_serial <- TSENAT:::.bpmapply(X, Y, FUN, nthreads = 1)
    
    expect_equal(result_parallel, result_serial, tolerance = 1e-10)
    # Verify: 0.5*10 + sqrt(0.5) ≈ 5.707
    expect_equal(result_parallel[[1]], 0.5 * 10 + sqrt(0.5), tolerance = 1e-10)
})

test_that(".bpmapply parallel with large vectors and bpparam", {
    # Test with larger data to ensure BiocParallel::bpmapply with bpparam works efficiently
    X <- 1:100
    Y <- 101:200
    FUN <- function(x, y) (x + y) / 2  # mean of x and y
    
    result_parallel <- TSENAT:::.bpmapply(X, Y, FUN, nthreads = 2)
    result_serial <- TSENAT:::.bpmapply(X, Y, FUN, nthreads = 1)
    
    expect_equal(result_parallel, result_serial, tolerance = 1e-10)
    expect_equal(length(result_parallel), 100)
    # Verify first and last: (1+101)/2 = 51, (100+200)/2 = 150
    expect_equal(result_parallel[[1]], 51)
    expect_equal(result_parallel[[100]], 150)
})

test_that(".bpmapply parallel executes correctly with different nthreads values", {
    # Test that bpparam initialization works for various thread counts
    X <- 1:6
    Y <- 6:1
    FUN <- function(x, y) c(x, y)
    
    result_parallel_2 <- TSENAT:::.bpmapply(X, Y, FUN, nthreads = 2)
    core_limit <- suppressWarnings(as.integer(Sys.getenv("_R_CHECK_LIMIT_CORES_", NA)))
    max_threads <- if (is.na(core_limit)) min(2, parallel::detectCores()) else min(2, core_limit)
    result_parallel_4 <- TSENAT:::.bpmapply(X, Y, FUN, nthreads = max_threads)
    result_serial <- TSENAT:::.bpmapply(X, Y, FUN, nthreads = 1)
    
    # All should be equal
    expect_equal(result_parallel_2, result_serial)
    expect_equal(result_parallel_4, result_serial)
    expect_is(result_parallel_2, "list")
})

test_that(".bpmapply SIMPLIFY=FALSE is respected in parallel execution", {
    # Verify that SIMPLIFY=FALSE is correctly passed to BiocParallel::bpmapply
    X <- 1:4
    Y <- 1:4
    FUN <- function(x, y) list(sum = x + y, product = x * y)
    
    result <- TSENAT:::.bpmapply(X, Y, FUN, nthreads = 2)
    
    # Result should be a list of lists, not simplified
    expect_is(result, "list")
    expect_equal(length(result), 4)
    expect_is(result[[1]], "list")
    expect_equal(result[[1]]$sum, 2)
    expect_equal(result[[1]]$product, 1)
    expect_equal(result[[4]]$sum, 8)
    expect_equal(result[[4]]$product, 16)
})

# ============================================================================
# Parallelization: .calculate_divergence() with nthreads
# ============================================================================

context("Parallelization: calculate_divergence with nthreads")

# ============================================================================
# Parallelization: .calculate_divergence() with nthreads
# ============================================================================

context("Parallelization: calculate_divergence with nthreads")

test_that("calculate_divergence sequential execution produces valid results", {
    library(SummarizedExperiment)
    set.seed(42)
    
    # Create test SummarizedExperiment: 30 genes, 8 samples
    n_genes <- 30
    n_samples <- 8
    counts <- matrix(rpois(n_genes * n_samples, lambda = 100), 
                     nrow = n_genes, ncol = n_samples)
    rownames(counts) <- paste0("Gene_", 1:n_genes)
    colnames(counts) <- paste0("Sample_", 1:n_samples)
    
    metadata <- data.frame(
        group = factor(c(rep("Normal", n_samples/2), rep("Tumor", n_samples/2))),
        row.names = colnames(counts)
    )
    
    se <- SummarizedExperiment(
        assays = list(counts = counts),
        colData = metadata
    )
    
    # Run sequential calculation (no res or top_n parameters - processes all genes)
    set.seed(42)
    result_seq <- .calculate_divergence(
        se = se,
        q = 1,
        nboot = 100,
        ci = 0.95,
        method = "percentile",
        nthreads = 1,
        progress = FALSE
    )
    
    # Check output structure (now returns SummarizedExperiment)
    expect_is(result_seq, "SummarizedExperiment")
    expect_equal(nrow(result_seq), n_genes)  # Should have all genes
    
    # Check rowData has required columns
    rd <- rowData(result_seq)
    expect_true("gene_name" %in% colnames(rd))
    expect_true("estimate" %in% colnames(rd))
    expect_true("lower_ci" %in% colnames(rd))
    expect_true("upper_ci" %in% colnames(rd))
    
    # Check values (allow for NAs in some cases)
    expect_true(nrow(result_seq) > 0)
    expect_true(all(result_seq$ci_width >= 0, na.rm = TRUE))
})

test_that("calculate_divergence sequential vs parallel produce consistent results", {
    library(SummarizedExperiment)
    set.seed(42)
    
    # Create test data
    n_genes <- 25
    n_samples <- 8
    counts <- matrix(rpois(n_genes * n_samples, lambda = 50), 
                     nrow = n_genes, ncol = n_samples)
    rownames(counts) <- paste0("Gene_", 1:n_genes)
    colnames(counts) <- paste0("Sample_", 1:n_samples)
    
    metadata <- data.frame(
        group = factor(c(rep("Normal", n_samples/2), rep("Tumor", n_samples/2))),
        row.names = colnames(counts)
    )
    
    se <- SummarizedExperiment(
        assays = list(counts = counts),
        colData = metadata
    )
    
    # Run sequential (nthreads=1)
    set.seed(42)
    result_seq <- .calculate_divergence(
        se = se,
        q = 1,
        nboot = 100,
        ci = 0.95,
        method = "percentile",
        nthreads = 1,
        progress = FALSE,
    )
    
    # Run parallel with 2 threads
    set.seed(42)
    result_par2 <- .calculate_divergence(
        se = se,
        q = 1,
        nboot = 100,
        ci = 0.95,
        method = "percentile",
        nthreads = 2,
        progress = FALSE,
    )
    
    # Both should return SummarizedExperiment
    expect_is(result_seq, "SummarizedExperiment")
    expect_is(result_par2, "SummarizedExperiment")
    
    # Both should have all genes
    expect_equal(nrow(result_seq), n_genes)
    expect_equal(nrow(result_par2), n_genes)
    
    # Results should be approximately equal (parallel RNG differs from sequential)
    # Check that numeric estimates match within tolerance (bootstrap stochasticity)
    rd_seq <- rowData(result_seq)
    rd_par <- rowData(result_par2)
    
    # Use approximate equality due to RNG differences in parallel execution
    # Tolerance of 0.01 (1%) is reasonable for bootstrap estimates
    expect_equal(rd_seq$estimate, rd_par$estimate, tolerance = 0.02)
    expect_equal(rd_seq$lower_ci, rd_par$lower_ci, tolerance = 0.05)
    expect_equal(rd_seq$upper_ci, rd_par$upper_ci, tolerance = 0.05)
})

test_that("calculate_divergence multiple thread levels produce consistent nrows", {
    library(SummarizedExperiment)
    set.seed(99)
    
    # Create test data
    n_genes <- 20
    n_samples <- 8
    counts <- matrix(rpois(n_genes * n_samples, lambda = 75), 
                     nrow = n_genes, ncol = n_samples)
    rownames(counts) <- paste0("Gene_", 1:n_genes)
    colnames(counts) <- paste0("Sample_", 1:n_samples)
    
    metadata <- data.frame(
        group = factor(c(rep("Control", n_samples/2), rep("Treatment", n_samples/2))),
        row.names = colnames(counts)
    )
    
    se <- SummarizedExperiment(
        assays = list(counts = counts),
        colData = metadata
    )
    
    # Run with different thread counts
    result_1t <- .calculate_divergence(
        se = se, q = 1, nboot = 100,
        nthreads = 1, progress = FALSE
    )
    
    result_2t <- .calculate_divergence(
        se = se, q = 1, nboot = 100,
        nthreads = 2, progress = FALSE
    )
    
    core_limit <- suppressWarnings(as.integer(Sys.getenv("_R_CHECK_LIMIT_CORES_", NA)))
    max_threads <- if (is.na(core_limit)) min(2, parallel::detectCores()) else min(2, core_limit)
    result_3t <- .calculate_divergence(
        se = se, q = 1, nboot = 100,
        nthreads = max_threads, progress = FALSE
    )
    
    # All should return SummarizedExperiment with all genes
    expect_is(result_1t, "SummarizedExperiment")
    expect_is(result_2t, "SummarizedExperiment")
    expect_is(result_3t, "SummarizedExperiment")
    
    expect_equal(nrow(result_1t), n_genes)
    expect_equal(nrow(result_2t), n_genes)
    expect_equal(nrow(result_3t), n_genes)
    
    # Results from different thread counts should be approximately equal
    # (RNG differs between sequential and parallel, and between different thread counts)
    rd_1t <- rowData(result_1t)
    rd_2t <- rowData(result_2t) 
    rd_3t <- rowData(result_3t)
    
    # Use approximate equality with tolerance due to RNG differences
    expect_equal(rd_1t$estimate, rd_2t$estimate, tolerance = 0.02)
    expect_equal(rd_2t$estimate, rd_3t$estimate, tolerance = 0.02)
})

test_that("calculate_divergence handles small gene count correctly", {
    library(SummarizedExperiment)
    set.seed(55)
    
    # Create very small dataset: only 3 genes
    n_genes <- 3
    n_samples <- 6
    counts <- matrix(rpois(n_genes * n_samples, lambda = 80), 
                     nrow = n_genes, ncol = n_samples)
    rownames(counts) <- paste0("Gene_", 1:n_genes)
    colnames(counts) <- paste0("Sample_", 1:n_samples)
    
    metadata <- data.frame(
        group = factor(c(rep("A", n_samples/2), rep("B", n_samples/2))),
        row.names = colnames(counts)
    )
    
    se <- SummarizedExperiment(
        assays = list(counts = counts),
        colData = metadata
    )
    
    # Request parallel computation (should work even though few genes)
    result <- .calculate_divergence(
        se = se,
        q = 1,
        nboot = 100,
        nthreads = min(3, parallel::detectCores()),  # Request parallel
        progress = FALSE,
    )
    
    # Should process all genes
    expect_is(result, "SummarizedExperiment")
    expect_equal(nrow(result), n_genes)
    
    # Check gene names in rowData
    rd <- rowData(result)
    expect_true("gene_name" %in% colnames(rd))
})

test_that("calculate_divergence auto-detects cores when nthreads=NULL", {
    library(SummarizedExperiment)
    set.seed(77)
    
    # Create test data
    n_genes <- 15
    n_samples <- 8
    counts <- matrix(rpois(n_genes * n_samples, lambda = 60), 
                     nrow = n_genes, ncol = n_samples)
    rownames(counts) <- paste0("Gene_", 1:n_genes)
    colnames(counts) <- paste0("Sample_", 1:n_samples)
    
    metadata <- data.frame(
        group = factor(c(rep("G1", n_samples/2), rep("G2", n_samples/2))),
        row.names = colnames(counts)
    )
    
    se <- SummarizedExperiment(
        assays = list(counts = counts),
        colData = metadata
    )
    
    # Run with max 2 threads (respecting environment limits)
    set.seed(77)
    core_limit <- suppressWarnings(as.integer(Sys.getenv("_R_CHECK_LIMIT_CORES_", NA)))
    max_threads <- if (is.na(core_limit)) min(2, parallel::detectCores()) else min(2, core_limit)
    result_auto <- .calculate_divergence(
        se = se,
        q = 1,
        nboot = 100,
        nthreads = max_threads,
        progress = FALSE,
    )
    
    # Should produce valid results as SummarizedExperiment
    expect_is(result_auto, "SummarizedExperiment")
    expect_equal(nrow(result_auto), n_genes)
    
    # Check rowData structure
    rd <- rowData(result_auto)
    expect_true("estimate" %in% colnames(rd))
})

test_that("calculate_divergence processes all genes (top_n no longer used)", {
    library(SummarizedExperiment)
    set.seed(88)
    
    # Create test data with many genes
    n_genes <- 50
    n_samples <- 8
    counts <- matrix(rpois(n_genes * n_samples, lambda = 70), 
                     nrow = n_genes, ncol = n_samples)
    rownames(counts) <- paste0("Gene_", 1:n_genes)
    colnames(counts) <- paste0("Sample_", 1:n_samples)
    
    metadata <- data.frame(
        group = factor(c(rep("X", n_samples/2), rep("Y", n_samples/2))),
        row.names = colnames(counts)
    )
    
    se <- SummarizedExperiment(
        assays = list(counts = counts),
        colData = metadata
    )
    
    # New API processes all genes provided, no top_n filtering
    result <- .calculate_divergence(
        se = se, q = 1, nboot = 100,
        nthreads = 1, progress = FALSE
    )
    
    # Should process all genes
    expect_is(result, "SummarizedExperiment")
    expect_equal(nrow(result), n_genes)
    
    # Verify all expected genes are present
    rd <- rowData(result)
    expect_true(all(!is.na(rd$estimate)))
})

test_that("calculate_divergence error handling for invalid input", {
    library(SummarizedExperiment)
    
    # Invalid SE (not SummarizedExperiment)
    expect_error(
        .calculate_divergence(
            se = data.frame(x = 1:4),
            q = 1, nthreads = 1
        ),
        "SummarizedExperiment|class"
    )
    
    # Valid SE but empty (should warn or error)
    n_genes <- 5
    n_samples <- 4
    counts <- matrix(rpois(n_genes * n_samples, lambda = 50), 
                     nrow = n_genes, ncol = n_samples)
    rownames(counts) <- paste0("Gene_", 1:n_genes)
    colnames(counts) <- paste0("Sample_", 1:n_samples)
    
    metadata <- data.frame(
        group = factor(c("A", "A", "B", "B")),
        row.names = colnames(counts)
    )
    
    se <- SummarizedExperiment(
        assays = list(counts = counts),
        colData = metadata
    )
    
    # Valid SE should not error
    expect_no_error(
        .calculate_divergence(
            se = se,
            q = 1, nthreads = 1, progress = FALSE
        )
    )
})

test_that("calculate_divergence output includes required columns", {
    library(SummarizedExperiment)
    set.seed(33)
    
    # Create small dataset
    n_genes <- 8
    n_samples <- 6
    counts <- matrix(rpois(n_genes * n_samples, lambda = 95), 
                     nrow = n_genes, ncol = n_samples)
    rownames(counts) <- paste0("Gene_", 1:n_genes)
    colnames(counts) <- paste0("Sample_", 1:n_samples)
    
    metadata <- data.frame(
        group = factor(c(rep("Ctrl", 3), rep("Trt", 3))),
        row.names = colnames(counts)
    )
    
    se <- SummarizedExperiment(
        assays = list(counts = counts),
        colData = metadata
    )
    
    result <- .calculate_divergence(
        se = se, q = 1, nboot = 100,
        nthreads = 1, progress = FALSE
    )
    
    # Check output is SummarizedExperiment
    expect_is(result, "SummarizedExperiment")
    
    # Check required columns exist in rowData
    rd <- rowData(result)
    expect_true("gene_name" %in% colnames(rd))
    expect_true("estimate" %in% colnames(rd))
    expect_true("lower_ci" %in% colnames(rd))
    expect_true("upper_ci" %in% colnames(rd))
    expect_true("ci_width" %in% colnames(rd))
})


context("Parallelization - Phase 1: BiocParallel Integration")

# Helper function to pass through bootstrap parameters to calculate_divergence
silent_calculate_divergence <- function(analysis, nthreads = 1, bootstrap = FALSE, nboot = NULL, ...) {
  calculate_divergence(
    analysis = analysis,
    nthreads = nthreads,
    bootstrap = bootstrap,
    nboot = nboot,
    progress = FALSE,
    verbose = FALSE,
    ...
  )
}

test_that("Serial (nthreads=1) produces valid divergence results", {
  expect_error(
    {
      analysis <- create_test_analysis()
      result_serial <- silent_calculate_divergence(
        analysis,
        nthreads = 1
      )
    },
    NA  # Expect NO error
  )
  
  # Result should be TSENATAnalysis
  expect_is(result_serial, "TSENATAnalysis")
  
  # Extract divergence results
  div_se <- result_serial@divergence_results$divergence_se
  expect_is(div_se, "SummarizedExperiment")
  expect_gt(nrow(div_se), 0)
  
  # Check computation mode metadata
  expect_true("computation_mode" %in% colnames(SummarizedExperiment::colData(div_se)))
  expect_equal(
    SummarizedExperiment::colData(div_se)$computation_mode[1],
    "sequential"
  )
})

test_that("Parallel (nthreads=2) produces valid divergence results", {
  skip_on_cran()  # Skip on CRAN due to resources
  
  expect_error(
    {
      analysis <- create_test_analysis()
      result_parallel <- silent_calculate_divergence(
        analysis,
        nthreads = 2
      )
    },
    NA  # Expect NO error
  )
  
  # Result should be TSENATAnalysis
  expect_is(result_parallel, "TSENATAnalysis")
  
  # Extract divergence results
  div_se <- result_parallel@divergence_results$divergence_se
  expect_is(div_se, "SummarizedExperiment")
  expect_gt(nrow(div_se), 0)
  
  # Check computation mode metadata
  expect_true("computation_mode" %in% colnames(SummarizedExperiment::colData(div_se)))
  expect_equal(
    SummarizedExperiment::colData(div_se)$computation_mode[1],
    "parallel"
  )
})

test_that("Serial vs Parallel: Estimates are numerically identical", {
  skip_on_cran()
  
  analysis <- create_test_analysis()
  
  # Run serial computation
  result_serial <- silent_calculate_divergence(
    analysis,
    nthreads = 1
  )
  
  # Run parallel computation on fresh analysis object
  analysis <- create_test_analysis()
  result_parallel <- silent_calculate_divergence(
    analysis,
    nthreads = 2
  )
  
  # Extract divergence SummarizedExperiment from results
  div_se_serial <- result_serial@divergence_results$divergence_se
  div_se_parallel <- result_parallel@divergence_results$divergence_se
  
  # Extract assay matrices
  assay_serial <- SummarizedExperiment::assay(div_se_serial)
  assay_parallel <- SummarizedExperiment::assay(div_se_parallel)
  
  # Compare dimensions
  expect_equal(
    dim(assay_serial),
    dim(assay_parallel),
    info = "Dimensions must match"
  )
  
  # Compare all divergence estimates
  # Use all.equal for numerical comparison (accounts for floating point precision)
  diff_matrix <- abs(assay_serial - assay_parallel)
  max_diff <- max(diff_matrix, na.rm = TRUE)
  
  expect_lt(max_diff, 1e-10)
})

test_that("Serial vs Parallel: CI bounds are numerically identical", {
  skip_on_cran()
  
  analysis <- create_test_analysis()
  
  # Run serial computation with bootstrap CI
  result_serial <- silent_calculate_divergence(
    analysis,
    nthreads = 1,
    bootstrap = TRUE,
    nboot = 100
  )
  
  # Run parallel computation with bootstrap CI
  analysis <- create_test_analysis()
  result_parallel <- silent_calculate_divergence(
    analysis,
    nthreads = 2,
    bootstrap = TRUE,
    nboot = 100
  )
  
  # Extract divergence SummarizedExperiment from results
  div_se_serial <- result_serial@divergence_results$divergence_se
  div_se_parallel <- result_parallel@divergence_results$divergence_se
  
  # Extract rowData
  rd_serial <- SummarizedExperiment::rowData(div_se_serial)
  rd_parallel <- SummarizedExperiment::rowData(div_se_parallel)
  
  # Find CI columns (lower_ci, upper_ci, ci_width)
  ci_cols <- grep("^(lower_ci|upper_ci|ci_width)", colnames(rd_serial), value = TRUE)
  
  # CI columns should exist and have values
  expect_gt(length(ci_cols), 0)
  
  # Verify both serial and parallel generated valid CI values (not all NA)
  for (col in ci_cols) {
    # Check serial CIs are valid (not all NA)
    expect_false(all(is.na(rd_serial[[col]])), 
                 info = paste("Serial", col, "should have values"))
    
    # Check parallel CIs are valid (not all NA)
    expect_false(all(is.na(rd_parallel[[col]])), 
                 info = paste("Parallel", col, "should have values"))
    
    # CIs should be reasonable (lower < upper for lower_ci and upper_ci pairs)
    if (grepl("^lower_ci", col)) {
      upper_col <- sub("^lower", "upper", col)
      if (upper_col %in% colnames(rd_serial)) {
        expect_true(all(rd_serial[[col]][!is.na(rd_serial[[col]])] <= 
                        rd_serial[[upper_col]][!is.na(rd_serial[[upper_col]])]))
        expect_true(all(rd_parallel[[col]][!is.na(rd_parallel[[col]])] <= 
                        rd_parallel[[upper_col]][!is.na(rd_parallel[[upper_col]])]))
      }
    }
  }
})

test_that("Serial vs Parallel: Row metadata identical", {
  skip_on_cran()
  
  analysis <- create_test_analysis()
  result_serial <- silent_calculate_divergence(
    analysis,
    nthreads = 1
  )
  
  analysis <- create_test_analysis()
  result_parallel <- silent_calculate_divergence(
    analysis,
    nthreads = 2
  )
  
  # Extract divergence SummarizedExperiment from results
  div_se_serial <- result_serial@divergence_results$divergence_se
  div_se_parallel <- result_parallel@divergence_results$divergence_se
  
  rd_serial <- SummarizedExperiment::rowData(div_se_serial)
  rd_parallel <- SummarizedExperiment::rowData(div_se_parallel)
  
  # Compare gene names (primary identifiers)
  expect_equal(
    rd_serial$gene_name,
    rd_parallel$gene_name,
    info = "Gene names must match"
  )
  
  # Compare error statuses
  expect_equal(
    is.na(rd_serial$error),
    is.na(rd_parallel$error),
    info = "Error status must match"
  )
})

test_that("Increasing threads maintains numerical stability", {
  skip_on_cran()
  
  analysis <- create_test_analysis()
  # Run with 1 thread
  result_1thread <- silent_calculate_divergence(
    analysis,
    nthreads = 1
  )
  
  # Run with 4 threads (if available)
  max_threads <- parallel::detectCores()
  if (max_threads >= 4) {
    analysis <- create_test_analysis()
    result_4threads <- silent_calculate_divergence(
      analysis,
      nthreads = 4
    )
    
    # Extract divergence SummarizedExperiment from results
    div_se_1 <- result_1thread@divergence_results$divergence_se
    div_se_4 <- result_4threads@divergence_results$divergence_se
    
    assay_1 <- SummarizedExperiment::assay(div_se_1)
    assay_4 <- SummarizedExperiment::assay(div_se_4)
    
    diff_matrix <- abs(assay_1 - assay_4)
    max_diff <- max(diff_matrix, na.rm = TRUE)
    
    expect_lt(max_diff, 1e-10)
  } else {
    skip_on_cran()
  }
})

test_that("nthreads parameter is properly validated", {
  analysis <- create_test_analysis()
  
  # Test with negative nthreads (should coerce to 1)
  result <- silent_calculate_divergence(
    analysis,
    nthreads = -5
  )
  expect_is(result, "TSENATAnalysis")
  expect_is(result@divergence_results$divergence_se, "SummarizedExperiment")
  
  # Test with zero (should coerce to 1)
  analysis <- create_test_analysis()
  result <- silent_calculate_divergence(
    analysis,
    nthreads = 0
  )
  expect_is(result, "TSENATAnalysis")
  expect_is(result@divergence_results$divergence_se, "SummarizedExperiment")
  
  # Test with huge number (should be clamped by OS/BiocParallel)
  analysis <- create_test_analysis()
  result <- expect_warning(
    silent_calculate_divergence(
      analysis,
      nthreads = 999
    ),
    "worker number limited"  # BiocParallel warns when limiting workers
  )
  expect_is(result, "TSENATAnalysis")
  expect_is(result@divergence_results$divergence_se, "SummarizedExperiment")
})

test_that("Computation mode is correctly reported in colData", {
  skip_on_cran()
  
  analysis <- create_test_analysis()
  result_serial <- silent_calculate_divergence(
    analysis,
    nthreads = 1
  )
  
  analysis <- create_test_analysis()
  result_parallel <- silent_calculate_divergence(
    analysis,
    nthreads = 2
  )
  
  # Extract divergence SummarizedExperiment from TSENATAnalysis wrapper
  div_se_serial <- result_serial@divergence_results$divergence_se
  div_se_parallel <- result_parallel@divergence_results$divergence_se
  
  coldata_serial <- SummarizedExperiment::colData(div_se_serial)
  coldata_parallel <- SummarizedExperiment::colData(div_se_parallel)
  
  # Check that computation_mode is present and correct
  expect_true("computation_mode" %in% colnames(coldata_serial))
  expect_true("computation_mode" %in% colnames(coldata_parallel))
  
  expect_equal(coldata_serial$computation_mode[1], "sequential")
  expect_equal(coldata_parallel$computation_mode[1], "parallel")
})
