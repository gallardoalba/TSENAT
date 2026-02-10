skip_on_bioc()

context("Extra diversity function tests")

library(TSENAT)

# calculate_tsallis_entropy argument validation and edge cases

test_that("calculate_tsallis_entropy validates inputs", {
    expect_error(calculate_tsallis_entropy("notnum", q = 2), "x must be numeric")
    expect_error(calculate_tsallis_entropy(c(1, 2, 3), q = "a"), "q must be numeric")
    expect_error(calculate_tsallis_entropy(c(1, 2, 3), q = c(-1, 2)), "q must be greater than 0")
})

test_that("calculate_tsallis_entropy handles zero-sum vectors and returns NA", {
    x <- c(0, 0, 0)
    expect_true(all(is.na(calculate_tsallis_entropy(x, q = 1, what = "S"))))
    expect_true(all(is.na(calculate_tsallis_entropy(x, q = 1, what = "D"))))
    both <- calculate_tsallis_entropy(x, q = c(0.5, 1, 2), what = "both")
    expect_true(all(is.na(both$S)))
    expect_true(all(is.na(both$D)))
})

test_that("calculate_tsallis_entropy computes expected values for simple distributions", {
    # single-dominant distribution -> entropy 0, diversity 1 for all q
    x <- c(10, 0, 0)
    S <- calculate_tsallis_entropy(x, q = c(0.5, 1, 2), norm = FALSE, what = "S")
    D <- calculate_tsallis_entropy(x, q = c(0.5, 1, 2), what = "D")
    expect_equal(as.numeric(S), rep(0, 3))
    expect_equal(as.numeric(D), rep(1, 3))

    # uniform distribution p = (1/3,1/3,1/3) with norm = TRUE should yield S in [0,1]
    x2 <- c(1, 1, 1)
    S_unif <- calculate_tsallis_entropy(x2, q = c(0.5, 1, 2), norm = TRUE, what = "S")
    expect_true(all(S_unif >= 0 & S_unif <= 1))

    # q=1 should match Shannon entropy normalization when norm=TRUE
    S_q1 <- calculate_tsallis_entropy(x2, q = 1, norm = TRUE, what = "S")
    # For uniform distribution, Shannon entropy = log(n)/log(n) = 1 when normalized
    expect_equal(as.numeric(S_q1), 1)
})

# .tsenat_prepare_diversity_input behaviours

test_that(".tsenat_prepare_diversity_input accepts data.frame and emits matrices", {
    df <- data.frame(S1 = c(1, 2), S2 = c(3, 4))
    res <- TSENAT:::.tsenat_prepare_diversity_input(df)
    expect_true(is.matrix(res$x))
    expect_null(res$se_assay_mat)
})

test_that(".tsenat_prepare_diversity_input handles tximport-like lists and tpm flag", {
    lst <- list(counts = matrix(1:4, nrow = 2), abundance = matrix(5:8, nrow = 2), other = 1, other2 = 2)
    res_counts <- TSENAT:::.tsenat_prepare_diversity_input(lst, tpm = FALSE)
    expect_true(is.matrix(res_counts$x))
    expect_equal(res_counts$x[1, 1], 1)

    res_abun <- TSENAT:::.tsenat_prepare_diversity_input(lst, tpm = TRUE)
    expect_equal(res_abun$x[1, 1], 5)
})

test_that(".tsenat_prepare_diversity_input warns/messages for tpm non-list inputs", {
    mat <- matrix(1:6, nrow = 3)
    expect_message(TSENAT:::.tsenat_prepare_diversity_input(mat, tpm = TRUE, verbose = TRUE), "tpm as a logical argument is only interpreted")
})

test_that(".tsenat_prepare_diversity_input handles SummarizedExperiment metadata readcounts and tx2gene mapping", {
    # Construct SE with metadata readcounts and tx2gene
    rc <- matrix(1:6, nrow = 3)
    rownames(rc) <- paste0("tx", 1:3)
    tx2 <- data.frame(Transcript = rownames(rc), Gen = c("g1", "g1", "g2"), stringsAsFactors = FALSE)

    se <- SummarizedExperiment::SummarizedExperiment(assays = S4Vectors::SimpleList(diversity = rc))
    S4Vectors::metadata(se)$readcounts <- rc
    S4Vectors::metadata(se)$tx2gene <- tx2

    res <- TSENAT:::.tsenat_prepare_diversity_input(se)
    expect_true(is.matrix(res$x))
    expect_equal(res$genes, c("g1", "g1", "g2"))
    expect_true(!is.null(res$se_assay_mat))
})

# invalid input types

test_that(".tsenat_prepare_diversity_input errors on unsupported input types", {
    expect_error(TSENAT:::.tsenat_prepare_diversity_input(12345), "Input data type is not supported")
})

# invalid assayno should error

test_that(".tsenat_prepare_diversity_input errors on invalid assayno for SummarizedExperiment", {
    rc <- matrix(1:4, nrow = 2)
    se <- SummarizedExperiment::SummarizedExperiment(assays = S4Vectors::SimpleList(a = rc))
    expect_error(TSENAT:::.tsenat_prepare_diversity_input(se, assayno = 2), "Please provide a valid assay number")
})
