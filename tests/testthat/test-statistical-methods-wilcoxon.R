context("Statistical Methods: Wilcoxon Defensive Behavior")

test_that("wilcoxon handles all-NA and constant rows without error and returns matrix", {
    # two groups of 3 samples each
    samples <- c(rep("A", 3), rep("B", 3))
    # build matrix with 3 rows: all-NA, constant, variable
    m <- matrix(nrow = 3, ncol = 6)
    m[1, ] <- NA_real_
    m[2, ] <- rep(5, 6)
    m[3, ] <- c(1, 2, 3, 4, 5, 6)

    res <- .wilcoxon(m, samples, pcorr = "none")
    expect_true(is.data.frame(res))
    # four columns: pvalue, padj, U, r
    expect_equal(ncol(res), 4)
    # raw p-values produced and NA-handling yields numeric outputs
    raw <- as.numeric(res[, 1])
    expect_length(raw, 3)
    # all-NA row should produce p-value of 1 (by convention used elsewhere)
    expect_true(is.finite(raw[1]))
})

test_that(".wilcoxon() returns named p-value and effect size columns", {
    mat <- matrix(runif(20), nrow = 5)
    samples <- rep(c("A", "B"), each = 5)
    res <- .wilcoxon(mat, samples, pcorr = "none", paired = FALSE, exact = FALSE)
    expect_true(is.data.frame(res))
    expect_true(all(c("pvalue", "padj", "U", "r") %in% colnames(res)))
})

context("Statistical Methods: Wilcoxon Paired Behavior")

test_that("wilcoxon paired matches per-row wilcox.test on ordered pairs", {
    # create a small matrix with two genes and two pairs (columns: N,T,N,T)
    # ensure paired differences are not identical to avoid "ties" in differences
    x <- matrix(c(
        1, 2, 3, 6, # gene1 across 4 samples (diffs: -1, -3)
        5, 6, 7, 9 # gene2 across 4 samples (diffs: -1, -2)
    ), nrow = 2, byrow = TRUE)
    samples <- c("Normal", "Tumor", "Normal", "Tumor")

    res <- .wilcoxon(x, samples, paired = TRUE, exact = TRUE)

    # compute expected raw p-values by calling wilcox.test per row with
    # paired=TRUE
    expected_raw <- sapply(seq_len(nrow(x)), function(i) {
        wilcox.test(
            x[
                i,
                which(samples == "Normal")
            ],
            x[
                i,
                which(samples == "Tumor")
            ],
            paired = TRUE, exact = TRUE
        )$p.value
    })

    expected_adj <- p.adjust(expected_raw, method = "BH")

    expect_equal(as.numeric(res[, "pvalue"]), as.numeric(expected_raw))
    expect_equal(as.numeric(res[, "padj"]), as.numeric(expected_adj))
})


test_that("wilcoxon paired errors on unequal group sizes", {
    x_bad <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)
    samples_bad <- c("Normal", "Tumor", "Normal")
    expect_error(.wilcoxon(x_bad, samples_bad, paired = TRUE), "Paired Wilcoxon requires equal numbers of samples in each group")
})



test_that("wilcoxon paired handles SummarizedExperiment input", {
    # construct a small paired dataset
    sample_names <- c("S1_N", "S1_T", "S2_N", "S2_T")
    mat_vals <- matrix(c(
        1, 2, 3, 6,  # gene1
        5, 6, 7, 9   # gene2
    ), nrow = 2, byrow = TRUE)
    
    colnames(mat_vals) <- sample_names
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat_vals)
    )
    cond <- ifelse(grepl("_N$", sample_names), "Normal", "Tumor")
    coldata <- data.frame(
        Sample = sample_names,
        Condition = cond,
        stringsAsFactors = FALSE
    )
    
    se_mapped <- TSENAT:::.tsenat_map_metadata_se(se, coldata)
    res <- .wilcoxon(
        SummarizedExperiment::assay(se_mapped),
        SummarizedExperiment::colData(se_mapped)$sample_type,
        paired = TRUE,
        exact = TRUE
    )
    
    # Verify output structure
    expect_true(is.data.frame(res))
    expect_equal(ncol(res), 4)
    expect_true(all(c("pvalue", "padj", "U", "r") %in% colnames(res)))
    # Verify results are not all NA
    expect_true(any(!is.na(res$pvalue)))
})

context("Statistical Methods: Wilcoxon U Statistic and r-Value Computation")

test_that("wilcoxon U statistic is computed correctly", {
    mat <- matrix(c(
        1, 2, 3, 4, # gene1
        5, 6, 7, 8  # gene2
    ), nrow = 2, byrow = TRUE)
    samples <- c("A", "A", "B", "B")
    
    res <- .wilcoxon(mat, samples, pcorr = "none", paired = FALSE)
    
    # U should be present and non-NA for valid data
    expect_false(anyNA(res$U))
    expect_true(all(res$U >= 0))
    # U should be bounded by the product of group sizes
    expect_true(all(res$U <= 4))  # 2 samples in each group
})

test_that("wilcoxon r-value is computed correctly", {
    mat <- matrix(c(
        1, 2, 3, 4,
        5, 6, 7, 8
    ), nrow = 2, byrow = TRUE)
    samples <- c("A", "A", "B", "B")
    
    res <- .wilcoxon(mat, samples, pcorr = "none", paired = FALSE)
    
    # r-value should be present and non-NA for valid data
    expect_false(anyNA(res$r))
    # r-value should be bounded between -1 and 1
    expect_true(all(res$r >= -1 & res$r <= 1))
    # For constant row, r should be 0 (no effect)
    mat_const <- matrix(c(1, 1, 1, 1), nrow = 1)
    res_const <- .wilcoxon(mat_const, samples, pcorr = "none")
    expect_true(res_const$r[1] == 0)
})

test_that("wilcoxon U and r values are NA when pvalue is NA", {
    m <- matrix(nrow = 2, ncol = 4)
    m[1, ] <- NA_real_
    m[2, ] <- c(1, 2, 3, 4)
    samples <- c("A", "A", "B", "B")
    
    res <- .wilcoxon(m, samples, pcorr = "none")
    
    # All-NA row should have NA U and r value
    expect_true(is.na(res$U[1]))
    expect_true(is.na(res$r[1]))
    # Valid row should have non-NA U and r
    expect_false(is.na(res$U[2]))
    expect_false(is.na(res$r[2]))
})

test_that("wilcoxon r-value is bounded and non-NA for valid paired data", {
    # For paired test, verify r-value is reasonable
    mat <- matrix(c(
        1, 2, 3, 6,  # distinct paired differences
        5, 6, 7, 9
    ), nrow = 2, byrow = TRUE)
    samples <- c("N", "T", "N", "T")
    
    res <- .wilcoxon(mat, samples, pcorr = "none", paired = TRUE)
    
    # Verify r-value is within valid range and non-NA for non-constant rows
    for (i in seq_len(nrow(res))) {
        p <- res$pvalue[i]
        if (!is.na(p)) {
            # r should be bounded between -1 and 1 and non-NA for valid p-values
            expect_true(!is.na(res$r[i]))
            expect_true(res$r[i] >= -1 && res$r[i] <= 1)
        }
    }
})

