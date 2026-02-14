context("Tsallis entropy calculations")

test_that("Tsallis entropy calculation is mathematically correct", {
    # Mathematical reference:
    # H_q(p) = (1 - sum(p_i^q)) / (q-1)
    # For q=1, H_1(p) = -sum(p_i * log2(p_i))
    read_counts <- c(0, 0, 5, 4, 1)
    p <- read_counts / sum(read_counts)

    # q = 2 (unnormalized)
    q2 <- 2
    manual_q2 <- (1 - sum(p^q2)) / (q2 - 1)
    tsallis_q2 <- calculate_tsallis_entropy(read_counts, q = q2, norm = FALSE)
    expect_equal(tsallis_q2, manual_q2, tolerance = 1e-8)

    # q = 2 (normalized)
    max_tsallis_q2 <- (1 - length(read_counts)^(1 - q2)) / (q2 - 1)
    manual_q2_norm <- manual_q2 / max_tsallis_q2
    tsallis_q2_norm <- calculate_tsallis_entropy(read_counts, q = q2, norm = TRUE)
    expect_equal(tsallis_q2_norm, manual_q2_norm, tolerance = 1e-8)
    expect_true(tsallis_q2_norm <= 1 && tsallis_q2_norm >= 0)

    # q = 1 (Shannon, unnormalized) -- use natural log by default
    manual_shannon <- -sum(ifelse(p > 0, p * log(p), 0))
    tsallis_q1 <- calculate_tsallis_entropy(read_counts, q = 1, norm = FALSE)
    expect_equal(tsallis_q1, manual_shannon, tolerance = 1e-8)

    # q = 1 (Shannon, normalized)
    manual_shannon_norm <- manual_shannon / log(length(read_counts))
    tsallis_q1_norm <- calculate_tsallis_entropy(read_counts, q = 1, norm = TRUE)
    expect_equal(tsallis_q1_norm, manual_shannon_norm, tolerance = 1e-8)
    expect_true(tsallis_q1_norm <= 1 && tsallis_q1_norm >= 0)

    # q = 1.5 (unnormalized)
    q15 <- 1.5
    manual_q15 <- (1 - sum(p^q15)) / (q15 - 1)
    tsallis_q15 <- calculate_tsallis_entropy(read_counts, q = q15, norm = FALSE)
    expect_equal(tsallis_q15, manual_q15, tolerance = 1e-8)

    # Vector q
    qvec <- c(1, 1.5, 2)
    tsallis_vec <- calculate_tsallis_entropy(read_counts, q = qvec, norm = FALSE)
    manual_vec <- vapply(qvec, function(qi) {
        if (abs(qi - 1) < .Machine$double.eps^0.5) {
            -sum(ifelse(p > 0, p * log(p), 0))
        } else {
            (1 - sum(p^qi)) / (qi - 1)
        }
    }, numeric(1))
    expect_equal(as.numeric(tsallis_vec),
        as.numeric(manual_vec),
        tolerance = 1e-8
    )
    expect_named(tsallis_vec, paste0("q=", qvec))

    # Edge cases
    expect_true(is.nan(calculate_tsallis_entropy(c(1), q = 2)))
    expect_true(is.na(calculate_tsallis_entropy(c(0, 0), q = 2)))
    expect_error(calculate_tsallis_entropy(read_counts, q = 0))
    expect_error(calculate_tsallis_entropy(read_counts, q = -1))
})

context("diversity_helpers extras")

library(testthat)

# .tsenat_calc_S: q ~= 1 and q != 1, normalized and not
test_that(".tsenat_calc_S computes Shannon and Tsallis correctly", {
    p <- c(0.5, 0.5)
    # Shannon with base 2: entropy = 1; normalized dividing by log2(2)=1 -> still 1
    s1 <- .tsenat_calc_S(p = p, q = 1, tol = 1e-8, n = 2, log_base = 2, norm = TRUE)
    expect_equal(s1, 1)
    # Tsallis q=2: S_2 = (1 - sum(p^2)) / (2-1) = 1 - (0.25 + 0.25) = 0.5
    s2 <- .tsenat_calc_S(p = p, q = 2, tol = 1e-8, n = 2, log_base = 2, norm = FALSE)
    expect_equal(s2, 0.5)
})

# .tsenat_calc_D: q close to 1 and other q
test_that(".tsenat_calc_D computes Hill numbers for q=1 and q!=1", {
    p <- c(0.5, 0.5)
    d1 <- .tsenat_calc_D(p = p, q = 1, tol = 1e-8, log_base = 2)
    # For q=1, sh = 1 (base 2), D1 = (log_base)^sh = 2^1 = 2
    expect_equal(d1, 2)
    d2 <- .tsenat_calc_D(p = p, q = 2, tol = 1e-8, log_base = 2)
    # For q=2, spq = sum(p^2)=0.5, Dq = spq^(1/(1-2)) = 0.5^( -1) = 2
    expect_equal(d2, 2)
})

# Input preparation errors and conversion
test_that(".tsenat_prepare_diversity_input rejects unsupported input types", {
    expect_error(.tsenat_prepare_diversity_input(1:5), "Input data type is not supported")
})

test_that(".tsenat_prepare_diversity_input handles data.frame conversion and provided genes", {
    df <- data.frame(a = 1:3, b = 2:4)
    res <- .tsenat_prepare_diversity_input(df, genes = c("g1", "g2", "g3"))
    expect_true(is.matrix(res$x))
    expect_equal(res$genes, c("g1", "g2", "g3"))
})

test_that(".tsenat_prepare_diversity_input handles tximport-like lists and tpm flag", {
    counts <- matrix(1:6, nrow = 3)
    abundance <- matrix(7:12, nrow = 3)
    # tximport-like lists are typically length 4 and contain named elements
    xlist <- list(counts = counts, abundance = abundance, txOut = TRUE, other = NULL)
    # default tpm = FALSE uses counts
    r1 <- .tsenat_prepare_diversity_input(xlist, genes = c("g1", "g2", "g3"))
    expect_true(is.matrix(r1$x))
    expect_equal(r1$x[1, 1], counts[1, 1])

    # tpm = TRUE uses abundance
    r2 <- .tsenat_prepare_diversity_input(xlist, genes = c("g1", "g2", "g3"), tpm = TRUE)
    expect_equal(r2$x[1, 1], abundance[1, 1])

    # improper list should error
    expect_error(.tsenat_prepare_diversity_input(list(foo = 1)), "cannot find any expression data")
})

test_that(".tsenat_prepare_diversity_input handles DGEList-like objects and messages when verbose", {
    counts <- matrix(rpois(6, lambda = 10), nrow = 3)
    dge <- list(counts = counts)
    class(dge) <- "DGEList"
    expect_message(.tsenat_prepare_diversity_input(dge, genes = c("g1", "g2", "g3"), verbose = TRUE), "DGEList contains transcript-level")
    expect_message(.tsenat_prepare_diversity_input(dge, genes = c("g1", "g2", "g3"), verbose = TRUE, tpm = TRUE), "tpm as a logical argument")
})

test_that(".tsenat_prepare_diversity_input handles SummarizedExperiment variants and tx2gene mapping", {
    mat <- matrix(1:6, nrow = 3)
    rownames(mat) <- paste0("tx", 1:3)
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat))
    # when genes not provided, should use rownames
    res <- .tsenat_prepare_diversity_input(se, genes = NULL)
    expect_true(is.matrix(res$x))
    expect_equal(res$genes, rownames(mat))

    # when metadata contains readcounts and tx2gene, prefer metadata mapping
    md <- list(readcounts = mat, tx2gene = data.frame(Transcript = paste0("tx", 1:3), Gen = c("gA", "gA", "gB"), stringsAsFactors = FALSE))
    se2 <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat), metadata = md)
    res2 <- .tsenat_prepare_diversity_input(se2, genes = NULL)
    expect_true(is.matrix(res2$x))
    expect_equal(res2$genes, c("gA", "gA", "gB"))

    # invalid assay number should error
    expect_error(.tsenat_prepare_diversity_input(se, genes = NULL, assayno = 10), "provide a valid assay number")
})

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
