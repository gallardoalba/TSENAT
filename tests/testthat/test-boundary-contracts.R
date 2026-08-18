context("AUDIT3 boundary contracts")

# Regression tests for the audit3 boundary conditions:
# A. q=0 entropy is independent of pseudocount (raw-support policy)
# B. non-natural log_base is rejected for multi-q spectra
# C. duplicated q within subject x condition is rejected deterministically
# D. TPM + effective-length double normalization is rejected

test_that("A: q=0 entropy uses raw support, independent of pseudocount", {
    x <- c(100, 0, 0)

    s0_no_pc <- .calculate_tsallis_entropy(x, q = 0, what = "S", pseudocount = 0)
    s0_pc <- .calculate_tsallis_entropy(x, q = 0, what = "S", pseudocount = 0.5)

    # Raw support is 1 observed isoform -> S_0 = 0, regardless of pseudocount
    expect_equal(unname(s0_pc), 0)
    expect_equal(unname(s0_pc), unname(s0_no_pc))

    # Hill D_0 = effective richness = raw support count = 1
    d0_pc <- .calculate_tsallis_entropy(x, q = 0, what = "D", pseudocount = 0.5)
    expect_equal(unname(d0_pc), 1)

    # .entropy_single follows the same policy
    e_no_pc <- .entropy_single(x, q = 0, norm = FALSE, pseudocount = 0)
    e_pc <- .entropy_single(x, q = 0, norm = FALSE, pseudocount = 0.5)
    expect_equal(e_pc, 0)
    expect_equal(e_pc, e_no_pc)

    # audit final2 exact regression: pseudocount = 0 vs 1 must agree at q = 0
    expect_equal(
        .calculate_tsallis_entropy(c(100, 0, 0), q = 0, pseudocount = 0),
        .calculate_tsallis_entropy(c(100, 0, 0), q = 0, pseudocount = 1),
        tolerance = 1e-12
    )
})

test_that("B: log_base != exp(1) is rejected for multi-q spectra", {
    x <- c(10, 5, 0)
    expect_error(
        .calculate_tsallis_entropy(x, q = c(1, 2), log_base = 2),
        "log_base"
    )
    expect_error(
        .calculate_tsallis_entropy(x, q = c(0.5, 1.5), log_base = 10),
        "log_base"
    )
    # single-q with non-natural base remains allowed (q=1 Shannon convention)
    expect_true(is.finite(.calculate_tsallis_entropy(x, q = 1, log_base = 2, norm = FALSE)))
})

test_that("C: duplicated q within subject x condition is rejected deterministically", {
    df <- data.frame(
        entropy = rnorm(8),
        q = c(0, 1, 2, 3, 0, 1, 2, 3),
        subject = factor(rep(1:2, each = 4)),
        condition = factor(rep(c("A", "B"), each = 4))
    )
    # duplicate q=3 within subject 1 x condition A block
    df_dup <- rbind(
        df,
        data.frame(
            entropy = 0,
            q = 3,
            subject = factor(1, levels = levels(df$subject)),
            condition = factor("A", levels = levels(df$condition))
        )
    )

    expect_error(
        .build_ar1_cor(df_dup, grid_col = "time_idx"),
        "Duplicated q values"
    )

    # clean data builds fine
    expect_true(is.list(.build_ar1_cor(df, grid_col = "time_idx")))
})

test_that("D: tpm = TRUE with effective_length is rejected (no double normalization)", {
    skip_if_not_installed("SummarizedExperiment")

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            counts = matrix(seq_len(12), nrow = 4),
            tpm = matrix(stats::runif(12), nrow = 4)
        ),
        colData = S4Vectors::DataFrame(sample = paste0("S", 1:3),
            condition = c("A", "A", "B"))
    )
    S4Vectors::metadata(se)$effective_length <- c(100, 200, 300, 400)
    rownames(se) <- paste0("TX", 1:4)
    genes <- c("G1", "G1", "G2", "G2")

    expect_error(
        .calculate_diversity(se, genes = genes, q = 1, tpm = TRUE, verbose = FALSE),
        "effective_length"
    )
})

test_that("E: pseudocount is added AFTER effective-length normalization", {
    x <- c(100, 0, 0)
    qs <- c(0.5, 1, 2)

    # With zero-count transcripts, swapping their effective lengths must not
    # change the entropy: the pseudocount is constant on the
    # effective-abundance scale. (It would vary under the old
    # pseudocount-before-length order, where the prior shrinks by 1/L_eff.)
    e1 <- .calculate_tsallis_entropy(x, q = qs, what = "S", norm = FALSE,
        pseudocount = 1, effective_length = c(1, 10, 100))
    e2 <- .calculate_tsallis_entropy(x, q = qs, what = "S", norm = FALSE,
        pseudocount = 1, effective_length = c(1, 100, 10))
    expect_equal(e1, e2, tolerance = 1e-12)

    # Manual q=1 (Shannon) check: effective abundances then pseudocount
    a <- c(100/1, 0/10, 0/100) + 1
    p <- a/sum(a)
    expect_equal(unname(e1["q=1"]), -sum(p * log(p)), tolerance = 1e-10)
})

test_that("F: pseudocount = 'auto' with TPM input is rejected", {
    skip_if_not_installed("SummarizedExperiment")

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            counts = matrix(stats::rpois(12, lambda = 5), nrow = 4),
            tpm = matrix(stats::runif(12, 0, 50), nrow = 4)
        ),
        colData = S4Vectors::DataFrame(sample = paste0("S", 1:3))
    )
    rownames(se) <- paste0("TX", 1:4)
    genes <- c("G1", "G1", "G2", "G2")

    expect_error(
        .calculate_diversity(se, genes = genes, q = 1, tpm = TRUE,
            pseudocount = "auto", verbose = FALSE),
        "'auto'"
    )
})

test_that("G: pseudocount = 'auto' is estimated from the resolved counts of a tximport-style list", {
    # Each gene gets >= 2 transcripts: with a single transcript the normalized
    # Shannon entropy is NaN by design (0/0), independent of pseudocount
    # resolution.
    counts <- matrix(stats::rpois(24, lambda = 10), nrow = 6)
    rownames(counts) <- paste0("TX", 1:6)
    colnames(counts) <- paste0("S", 1:4)
    lst <- list(
        counts = counts,
        abundance = counts/rowSums(counts) * 1e6,
        length = rep(1000, 6),
        countsFromAbundance = "no"
    )
    genes <- c("G1", "G1", "G2", "G2", "G3", "G3")

    res <- .calculate_diversity(lst, genes = genes, q = 1,
        pseudocount = "auto", min_valid_frac = 0, verbose = FALSE)

    expect_s4_class(res, "SummarizedExperiment")
    expect_equal(nrow(res), 3L)  # 3 genes (min_valid_frac = 0 keeps sparse genes)
    expect_true(all(is.finite(SummarizedExperiment::assay(res))))
})
