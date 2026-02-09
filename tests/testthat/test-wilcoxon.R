context("Wilcoxon defensive behavior")

test_that("wilcoxon handles all-NA and constant rows without error and returns matrix", {
    # two groups of 3 samples each
    samples <- c(rep("A", 3), rep("B", 3))
    # build matrix with 3 rows: all-NA, constant, variable
    m <- matrix(nrow = 3, ncol = 6)
    m[1, ] <- NA_real_
    m[2, ] <- rep(5, 6)
    m[3, ] <- c(1, 2, 3, 4, 5, 6)

    res <- wilcoxon(m, samples, pcorr = "none")
    expect_true(is.matrix(res) || is.data.frame(res))
    # two columns: raw and adjusted
    expect_equal(ncol(res), 2)
    # raw p-values produced and NA-handling yields numeric outputs
    raw <- as.numeric(res[, 1])
    expect_length(raw, 3)
    # all-NA row should produce p-value of 1 (by convention used elsewhere)
    expect_true(is.finite(raw[1]))
})

test_that("wilcoxon() returns named p-value columns", {
    mat <- matrix(runif(20), nrow = 5)
    samples <- rep(c("A", "B"), each = 5)
    res <- wilcoxon(mat, samples, pcorr = "none", paired = FALSE, exact = FALSE)
    expect_true(is.matrix(res) || is.data.frame(res))
    expect_true(all(c("raw_p_values", "adjusted_p_values") %in% colnames(res)))
})

context("Wilcoxon paired behavior")

test_that("wilcoxon paired matches per-row wilcox.test on ordered pairs", {
    # create a small matrix with two genes and two pairs (columns: N,T,N,T)
    # ensure paired differences are not identical to avoid "ties" in differences
    x <- matrix(c(
        1, 2, 3, 6, # gene1 across 4 samples (diffs: -1, -3)
        5, 6, 7, 9 # gene2 across 4 samples (diffs: -1, -2)
    ), nrow = 2, byrow = TRUE)
    samples <- c("Normal", "Tumor", "Normal", "Tumor")

    res <- wilcoxon(x, samples, paired = TRUE, exact = TRUE)

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

    expect_equal(as.numeric(res[, 1]), as.numeric(expected_raw))
    expect_equal(as.numeric(res[, 2]), as.numeric(expected_adj))
})


test_that("wilcoxon paired errors on unequal group sizes", {
    x_bad <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)
    samples_bad <- c("Normal", "Tumor", "Normal")
    expect_error(wilcoxon(x_bad, samples_bad, paired = TRUE), "Paired Wilcoxon requires equal numbers of samples in each group")
})


test_that("wilcoxon paired is invariant to shuffled SE when mapped with paired=TRUE", {
    # construct a small paired dataset with two bases (S1, S2)
    sample_names <- c("S1_N", "S1_T", "S2_N", "S2_T")
    # create values so paired differences are non-tied
    mat_vals <- matrix(c(
        1, 2, 3, 6,  # gene1
        5, 6, 7, 9   # gene2
    ), nrow = 2, byrow = TRUE)

    cols <- paste0(sample_names, "_q=0.1")
    colnames(mat_vals) <- cols

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat_vals)
    )
    # base (ordered) coldata
    cond <- ifelse(grepl("_N$", sample_names), "Normal", "Tumor")
    coldata_base <- data.frame(
        Sample = sample_names,
        Condition = cond,
        stringsAsFactors = FALSE
    )

    se_mapped_base <- map_metadata(se, coldata_base, paired = TRUE)
    res_base <- wilcoxon(
        SummarizedExperiment::assay(se_mapped_base),
        SummarizedExperiment::colData(se_mapped_base)$sample_type,
        paired = TRUE,
        exact = TRUE
    )

    # shuffled SE and shuffled coldata rows
    perm <- c(2, 3, 1, 4) # deterministic shuffle
    se_shuffled <- se[, perm]
    coldata_shuffled <- coldata_base[perm, , drop = FALSE]

    se_mapped_shuffled <- map_metadata(se_shuffled,
        coldata_shuffled, paired = TRUE)
    res_shuffled <- wilcoxon(
        SummarizedExperiment::assay(se_mapped_shuffled),
        SummarizedExperiment::colData(se_mapped_shuffled)$sample_type,
        paired = TRUE,
        exact = TRUE
    )

    expect_equal(as.numeric(res_base[, 1]), as.numeric(res_shuffled[, 1]))
    expect_equal(as.numeric(res_base[, 2]), as.numeric(res_shuffled[, 2]))
})
