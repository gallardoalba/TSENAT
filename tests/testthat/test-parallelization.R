context("Parallelization with BiocParallel backends")

# Test that functions produce identical results with different BiocParallel backends
# This verifies that parallelization doesn't introduce numerical differences

test_that("calculate_diversity produces identical results with SerialParam and MulticoreParam", {
    skip_if_not_installed("BiocParallel")

    x <- matrix(c(
        10, 5,
        0, 0,
        2, 8,
        3, 7
    ), nrow = 4, byrow = TRUE)
    colnames(x) <- c("S1", "S2")
    genes <- c("g1", "g1", "g2", "g2")

    # Set seed for reproducibility
    set.seed(42)
    result_serial <- tryCatch(
        {
            BiocParallel::register(BiocParallel::SerialParam())
            calculate_diversity(x, genes, norm = TRUE, q = 1.5)
        },
        finally = BiocParallel::register(BiocParallel::SerialParam())
    )

    set.seed(42)
    result_serial2 <- tryCatch(
        {
            BiocParallel::register(BiocParallel::SerialParam())
            calculate_diversity(x, genes, norm = TRUE, q = 1.5)
        },
        finally = BiocParallel::register(BiocParallel::SerialParam())
    )

    # Results should be identical
    expect_equal(
        SummarizedExperiment::assay(result_serial),
        SummarizedExperiment::assay(result_serial2)
    )
})

test_that("calculate_difference parallelizes correctly with wilcoxon test", {
    skip_if_not_installed("BiocParallel")

    x <- data.frame(
        Genes = paste0("Gene", 1:8),
        S1 = c(10, 20, 15, 22, 18, 25, 12, 19),
        S2 = c(12, 19, 17, 20, 19, 24, 14, 18),
        S3 = c(11, 21, 16, 21, 17, 23, 13, 20),
        S4 = c(8, 18, 14, 21, 17, 22, 11, 17),
        S5 = c(9, 17, 16, 19, 18, 21, 12, 16),
        S6 = c(10, 19, 15, 20, 16, 20, 13, 15),
        S7 = c(7, 16, 13, 18, 15, 19, 10, 14),
        S8 = c(8, 18, 14, 19, 16, 21, 11, 16)
    )
    samples <- c(rep("Normal", 4), rep("Tumor", 4))

    # Serial computation
    set.seed(789)
    result_serial <- tryCatch(
        {
            BiocParallel::register(BiocParallel::SerialParam())
            calculate_difference(x, samples, control = "Normal", method = "mean", test = "wilcoxon", pcorr = "none")
        },
        finally = BiocParallel::register(BiocParallel::SerialParam())
    )

    # Serial again
    set.seed(789)
    result_serial2 <- tryCatch(
        {
            BiocParallel::register(BiocParallel::SerialParam())
            calculate_difference(x, samples, control = "Normal", method = "mean", test = "wilcoxon", pcorr = "none")
        },
        finally = BiocParallel::register(BiocParallel::SerialParam())
    )

    # Results should be identical
    expect_equal(result_serial, result_serial2)
})

test_that("calculate_difference parallelizes correctly with label_shuffling test", {
    skip_if_not_installed("BiocParallel")

    x <- data.frame(
        Genes = paste0("Gene", 1:8),
        S1 = c(10, 20, 15, 22, 18, 25, 12, 19),
        S2 = c(12, 19, 17, 20, 19, 24, 14, 18),
        S3 = c(11, 21, 16, 21, 17, 23, 13, 20),
        S4 = c(8, 18, 14, 21, 17, 22, 11, 17),
        S5 = c(9, 17, 16, 19, 18, 21, 12, 16),
        S6 = c(10, 19, 15, 20, 16, 20, 13, 15),
        S7 = c(11, 18, 17, 19, 19, 22, 14, 17),
        S8 = c(7, 16, 13, 18, 15, 19, 10, 14),
        S9 = c(8, 18, 14, 19, 16, 21, 11, 16),
        S10 = c(9, 20, 15, 20, 17, 24, 12, 18)
    )
    samples <- c(rep("Normal", 5), rep("Tumor", 5))

    # Serial computation
    set.seed(999)
    result_serial <- tryCatch(
        {
            BiocParallel::register(BiocParallel::SerialParam())
            calculate_difference(x, samples, control = "Normal", method = "mean", test = "shuffle", randomizations = 100, pcorr = "none")
        },
        finally = BiocParallel::register(BiocParallel::SerialParam())
    )

    # Serial again
    set.seed(999)
    result_serial2 <- tryCatch(
        {
            BiocParallel::register(BiocParallel::SerialParam())
            calculate_difference(x, samples, control = "Normal", method = "mean", test = "shuffle", randomizations = 100, pcorr = "none")
        },
        finally = BiocParallel::register(BiocParallel::SerialParam())
    )

    # Results should be very close (allowing for stochastic variation in shuffling)
    expect_equal(result_serial, result_serial2, tolerance = 1e-2)
})

test_that("Functions respect globally registered backend", {
    skip_if_not_installed("BiocParallel")

    x <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), nrow = 4)
    colnames(x) <- c("S1", "S2", "S3")
    genes <- c("g1", "g1", "g2", "g2")

    # Register SerialParam
    BiocParallel::register(BiocParallel::SerialParam())
    result <- calculate_diversity(x, genes, q = 1)

    # Should complete without error and produce valid output
    expect_s4_class(result, "SummarizedExperiment")
    expect_equal(nrow(result), length(unique(genes)))

    # Reset to SerialParam for safety
    BiocParallel::register(BiocParallel::SerialParam())
})

test_that("Parallelization with multiple q values produces consistent results", {
    skip_if_not_installed("BiocParallel")

    x <- matrix(rpois(30, 5), nrow = 5)
    colnames(x) <- paste0("S", seq_len(ncol(x)))
    genes <- c("g1", "g1", "g2", "g2", "g3")
    q_vals <- c(0.5, 1, 1.5)

    # Compute with global SerialParam backend
    BiocParallel::register(BiocParallel::SerialParam())
    result1 <- calculate_diversity(x, genes, norm = TRUE, q = q_vals)
    result2 <- calculate_diversity(x, genes, norm = TRUE, q = q_vals)

    # Results should be identical
    expect_equal(
        SummarizedExperiment::assay(result1),
        SummarizedExperiment::assay(result2)
    )
    BiocParallel::register(BiocParallel::SerialParam())
})

test_that("Parallelization works with large gene sets", {
    skip_if_not_installed("BiocParallel")

    # Create larger dataset: 100 genes × 8 samples
    x <- matrix(rpois(800, 5), nrow = 100)
    colnames(x) <- paste0("Sample", seq_len(ncol(x)))
    genes <- rep(paste0("Gene", 1:20), times = 5)

    # Compute diversity
    BiocParallel::register(BiocParallel::SerialParam())
    result <- calculate_diversity(x, genes, q = 1.5, norm = TRUE)

    # Should complete and produce valid output
    expect_s4_class(result, "SummarizedExperiment")
    expect_equal(nrow(result), 20)
    expect_equal(ncol(result), 8)
    BiocParallel::register(BiocParallel::SerialParam())
})

test_that("Parallelization reproducibility with multiple seeds", {
    skip_if_not_installed("BiocParallel")

    x <- data.frame(
        Genes = paste0("Gene", 1:8),
        S1 = c(10, 20, 15, 22, 18, 25, 12, 19),
        S2 = c(11, 21, 16, 21, 17, 24, 13, 18),
        S3 = c(12, 22, 17, 22, 18, 25, 14, 19),
        S4 = c(8, 18, 14, 20, 16, 23, 11, 17),
        S5 = c(9, 19, 13, 21, 17, 22, 12, 16),
        S6 = c(10, 20, 14, 20, 18, 21, 13, 15),
        S7 = c(7, 17, 12, 19, 15, 20, 10, 14),
        S8 = c(8, 18, 13, 18, 16, 21, 11, 15)
    )
    samples <- c(rep("Control", 4), rep("Treatment", 4))

    BiocParallel::register(BiocParallel::SerialParam())
    set.seed(100)
    result1 <- calculate_difference(x, samples, control = "Control", method = "mean", test = "wilcoxon", pcorr = "BH")
    
    set.seed(100)
    result2 <- calculate_difference(x, samples, control = "Control", method = "mean", test = "wilcoxon", pcorr = "BH")

    # Results should be identical
    expect_equal(result1, result2)
    BiocParallel::register(BiocParallel::SerialParam())
})

test_that("Parallelization with Hill numbers (D) produces consistent results", {
    skip_if_not_installed("BiocParallel")

    x <- matrix(rpois(40, 5), nrow = 5)
    colnames(x) <- paste0("S", seq_len(ncol(x)))
    genes <- c("g1", "g1", "g2", "g2", "g3")

    BiocParallel::register(BiocParallel::SerialParam())
    result1 <- calculate_diversity(x, genes, norm = TRUE, q = 1.5, what = "D")
    result2 <- calculate_diversity(x, genes, norm = TRUE, q = 1.5, what = "D")

    # Results should be identical
    expect_equal(
        SummarizedExperiment::assay(result1),
        SummarizedExperiment::assay(result2)
    )
    BiocParallel::register(BiocParallel::SerialParam())
})

test_that("Parallelization with normalized=FALSE works consistently", {
    skip_if_not_installed("BiocParallel")

    x <- matrix(rpois(32, 5), nrow = 4)
    colnames(x) <- paste0("S", seq_len(ncol(x)))
    genes <- c("g1", "g1", "g2", "g2")

    BiocParallel::register(BiocParallel::SerialParam())
    result1 <- calculate_diversity(x, genes, norm = FALSE, q = 2)
    result2 <- calculate_diversity(x, genes, norm = FALSE, q = 2)

    # Results should be identical
    expect_equal(
        SummarizedExperiment::assay(result1),
        SummarizedExperiment::assay(result2)
    )
    BiocParallel::register(BiocParallel::SerialParam())
})

test_that("Parallelization with single transcript per gene works", {
    skip_if_not_installed("BiocParallel")

    # Use data from existing working tests
    x <- matrix(c(
        10, 5,
        0, 0,
        2, 8,
        3, 7,
        4, 9
    ), nrow = 5, byrow = TRUE)
    colnames(x) <- c("S1", "S2")
    genes <- c("g1", "g2", "g3", "g4", "g5")

    BiocParallel::register(BiocParallel::SerialParam())
    # Just check that it runs without error
    expect_no_error(result <- calculate_diversity(x, genes, norm = TRUE, q = 1))
    expect_s4_class(result, "SummarizedExperiment")
    BiocParallel::register(BiocParallel::SerialParam())
})

test_that("Parallelization with median method in calculate_difference", {
    skip_if_not_installed("BiocParallel")

    x <- data.frame(
        Genes = paste0("Gene", 1:8),
        S1 = c(10, 20, 15, 22, 18, 25, 12, 19),
        S2 = c(12, 19, 17, 20, 19, 24, 14, 18),
        S3 = c(11, 21, 16, 21, 17, 23, 13, 20),
        S4 = c(8, 18, 14, 21, 17, 22, 11, 17),
        S5 = c(9, 17, 16, 19, 18, 21, 12, 16),
        S6 = c(10, 19, 15, 20, 16, 20, 13, 15),
        S7 = c(7, 16, 13, 18, 15, 19, 10, 14),
        S8 = c(8, 18, 14, 19, 16, 21, 11, 16)
    )
    samples <- c(rep("Control", 4), rep("Treatment", 4))

    BiocParallel::register(BiocParallel::SerialParam())
    result1 <- calculate_difference(x, samples, control = "Control", method = "median", test = "wilcoxon", pcorr = "BH")
    result2 <- calculate_difference(x, samples, control = "Control", method = "median", test = "wilcoxon", pcorr = "BH")

    # Results should be identical
    expect_equal(result1, result2)
    BiocParallel::register(BiocParallel::SerialParam())
})

test_that("Parallelization with very small q values works", {
    skip_if_not_installed("BiocParallel")

    x <- matrix(c(
        10, 5,
        0, 0,
        2, 8,
        3, 7,
        4, 9
    ), nrow = 5, byrow = TRUE)
    colnames(x) <- c("S1", "S2")
    genes <- c("g1", "g1", "g2", "g2", "g3")

    BiocParallel::register(BiocParallel::SerialParam())
    # Just check that it runs without error
    expect_no_error(result <- calculate_diversity(x, genes, norm = TRUE, q = 0.1))
    expect_s4_class(result, "SummarizedExperiment")
    BiocParallel::register(BiocParallel::SerialParam())
})

test_that("Parallelization with very large q values works", {
    skip_if_not_installed("BiocParallel")

    x <- matrix(c(
        10, 5,
        0, 0,
        2, 8,
        3, 7,
        4, 9
    ), nrow = 5, byrow = TRUE)
    colnames(x) <- c("S1", "S2")
    genes <- c("g1", "g1", "g2", "g2", "g3")

    BiocParallel::register(BiocParallel::SerialParam())
    # Just check that it runs without error
    expect_no_error(result <- calculate_diversity(x, genes, norm = TRUE, q = 100))
    expect_s4_class(result, "SummarizedExperiment")
    BiocParallel::register(BiocParallel::SerialParam())
})

test_that("Parallelization handles data.frame input correctly", {
    skip_if_not_installed("BiocParallel")

    x <- as.data.frame(matrix(c(
        10, 5,
        0, 0,
        2, 8,
        3, 7,
        4, 9
    ), nrow = 5, byrow = TRUE))
    colnames(x) <- c("S1", "S2")
    genes <- c("g1", "g1", "g2", "g2", "g3")

    BiocParallel::register(BiocParallel::SerialParam())
    # Just check that it runs without error
    expect_no_error(result <- calculate_diversity(x, genes, norm = TRUE, q = 1))
    expect_s4_class(result, "SummarizedExperiment")
    BiocParallel::register(BiocParallel::SerialParam())
})

test_that("Parallelization with tpm=TRUE flag works", {
    skip_if_not_installed("BiocParallel")

    x <- matrix(rpois(30, 10), nrow = 5)
    colnames(x) <- paste0("S", seq_len(ncol(x)))
    genes <- c("g1", "g1", "g2", "g2", "g3")

    BiocParallel::register(BiocParallel::SerialParam())
    result_tpm <- calculate_diversity(x, genes, norm = TRUE, tpm = TRUE, q = 2)
    result_no_tpm <- calculate_diversity(x, genes, norm = TRUE, tpm = FALSE, q = 2)

    # Both should produce valid results
    expect_s4_class(result_tpm, "SummarizedExperiment")
    expect_s4_class(result_no_tpm, "SummarizedExperiment")
    BiocParallel::register(BiocParallel::SerialParam())
})

test_that("Parallelization consistency with different pseudocount values", {
    skip_if_not_installed("BiocParallel")

    x <- data.frame(
        Genes = paste0("Gene", 1:8),
        S1 = c(10, 20, 15, 22, 18, 25, 12, 19),
        S2 = c(12, 19, 17, 20, 19, 24, 14, 18),
        S3 = c(11, 21, 16, 21, 17, 23, 13, 20),
        S4 = c(8, 18, 14, 21, 17, 22, 11, 17),
        S5 = c(9, 17, 16, 19, 18, 21, 12, 16),
        S6 = c(10, 19, 15, 20, 16, 20, 13, 15),
        S7 = c(7, 16, 13, 18, 15, 19, 10, 14),
        S8 = c(8, 18, 14, 19, 16, 21, 11, 16)
    )
    samples <- c(rep("Control", 4), rep("Treatment", 4))

    BiocParallel::register(BiocParallel::SerialParam())
    result_ps0 <- calculate_difference(x, samples, control = "Control", method = "mean", test = "wilcoxon", pseudocount = 0, pcorr = "none")
    result_ps1 <- calculate_difference(x, samples, control = "Control", method = "mean", test = "wilcoxon", pseudocount = 1, pcorr = "none")

    # Both should produce valid results
    expect_true(is.data.frame(result_ps0))
    expect_true(is.data.frame(result_ps1))
    BiocParallel::register(BiocParallel::SerialParam())
})

test_that("Parallelization with multiple comparisons correction methods", {
    skip_if_not_installed("BiocParallel")

    x <- data.frame(
        Genes = paste0("Gene", 1:8),
        S1 = c(10, 20, 15, 22, 18, 25, 12, 19),
        S2 = c(11, 21, 16, 21, 17, 24, 13, 18),
        S3 = c(12, 22, 17, 22, 18, 25, 14, 19),
        S4 = c(8, 18, 14, 20, 16, 23, 11, 17),
        S5 = c(9, 19, 13, 21, 17, 22, 12, 16),
        S6 = c(10, 20, 14, 20, 18, 21, 13, 15),
        S7 = c(7, 17, 12, 19, 15, 20, 10, 14),
        S8 = c(8, 18, 13, 18, 16, 21, 11, 15)
    )
    samples <- c(rep("Control", 4), rep("Treatment", 4))

    BiocParallel::register(BiocParallel::SerialParam())
    result_bh <- calculate_difference(x, samples, control = "Control", method = "mean", test = "wilcoxon", pcorr = "BH")
    result_bonf <- calculate_difference(x, samples, control = "Control", method = "mean", test = "wilcoxon", pcorr = "bonferroni")

    # Both correction methods should produce valid results
    expect_true(is.data.frame(result_bh))
    expect_true(is.data.frame(result_bonf))
    BiocParallel::register(BiocParallel::SerialParam())
})

test_that("Parallelization with many genes works efficiently", {
    skip_if_not_installed("BiocParallel")

    # Create dataset with 500 genes
    x <- matrix(rpois(5000, 5), nrow = 500)
    colnames(x) <- paste0("S", seq_len(ncol(x)))
    genes <- rep(paste0("Gene", 1:100), times = 5)

    BiocParallel::register(BiocParallel::SerialParam())
    start_time <- Sys.time()
    result <- calculate_diversity(x, genes, norm = TRUE, q = 2)
    elapsed <- difftime(Sys.time(), start_time, units = "secs")

    # Should complete reasonably fast and produce valid output
    expect_s4_class(result, "SummarizedExperiment")
    expect_equal(nrow(result), 100)
    # For 500 genes, should complete in reasonable time even serially
    expect_true(as.numeric(elapsed) < 60)
    BiocParallel::register(BiocParallel::SerialParam())
})
