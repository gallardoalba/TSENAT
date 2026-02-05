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

    se_mapped_base <- map_coldata_to_se(se, coldata_base, paired = TRUE)
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

    se_mapped_shuffled <- map_coldata_to_se(se_shuffled,
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
