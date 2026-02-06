library(testthat)

context("calculate_fc pseudocount behavior")

test_that("scale-aware pseudocount replaces non-positive group summaries", {
    # construct a small matrix with two genes and four samples (2 groups)
    # samples: Normal,Normal,Tumor,Tumor
    samples <- c("Normal", "Normal", "Tumor", "Tumor")
    # gene1: Normal means = 0, Tumor means = 4
    # gene2: Normal means = 2, Tumor means = 4
    mat <- matrix(c(
        0, 0, 4, 4,  # gene1
        2, 2, 4, 4   # gene2
    ), byrow = TRUE, nrow = 2)
    rownames(mat) <- c("gene1", "gene2")

    res_auto <- calculate_fc(mat, samples, control = "Normal", method = "mean", pseudocount = 0)
    # Note: `calculate_fc` returns log2(non-control / control). With
    # control = "Normal" the non-control (Tumor) values are 4 and 4.
    # After automatic pseudocount selection (half of smallest positive = 1),
    # gene1: log2(4 / 1) = 2; gene2: log2(4 / 2) = 1.
    expect_equal(as.numeric(res_auto$log2_fold_change[1]), 2, tolerance = 1e-8)
    expect_equal(as.numeric(res_auto$log2_fold_change[2]), 1, tolerance = 1e-8)
})

test_that("explicit pseudocount is honored and differs from automatic choice", {
    samples <- c("Normal", "Normal", "Tumor", "Tumor")
    mat <- matrix(c(
        0, 0, 4, 4,
        2, 2, 4, 4
    ), byrow = TRUE, nrow = 2)
    rownames(mat) <- c("gene1", "gene2")

    res_explicit <- calculate_fc(mat, samples, control = "Normal", method = "mean", pseudocount = 1e-6)
    # explicit tiny pseudocount should lead to a very large positive log2fc
    # (Tumor / near-zero Normal). Expect log2fc >> 10
    expect_true(as.numeric(res_explicit$log2_fold_change[1]) > 10)
    # gene2 unaffected by small pseudocount (still ~ 1)
    expect_equal(as.numeric(res_explicit$log2_fold_change[2]), 1, tolerance = 1e-8)
})

test_that("calculate_difference forwards pseudocount to calculate_fc", {
    samples <- c("Normal", "Normal", "Tumor", "Tumor")
    mat <- matrix(c(
        0, 0, 4, 4,
        2, 2, 4, 4
    ), byrow = TRUE, nrow = 2)
    rownames(mat) <- c("gene1", "gene2")

    # build data.frame as calculate_difference expects (genes col + sample cols)
    df <- as.data.frame(mat)
    df <- cbind(genes = rownames(mat), df)
    # call calculate_difference which should forward pseudocount
    # small sample sizes trigger a warning; assert it and capture the result
    res <- testthat::expect_warning(
        calculate_difference(df, samples = samples, control = "Normal", method = "mean", test = "wilcoxon", paired = FALSE, pseudocount = 0),
        "Low sample size for wilcoxon"
    )
    # find matching rows and compare log2_fold_change
    expect_equal(as.numeric(res$log2_fold_change[1]), 2, tolerance = 1e-8)
    expect_equal(as.numeric(res$log2_fold_change[2]), 1, tolerance = 1e-8)
})
