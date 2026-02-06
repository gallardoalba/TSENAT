library(testthat)

context("calculate_lm_interaction lmm pvalue options")

test_that("lmm returns LRT and Satterthwaite p-values and respects pvalue arg", {
    qvec <- seq(0.01, 0.1, by = 0.01)
    sample_names <- rep(c("S1_N", "S2_T"), each = length(qvec))
    coln <- paste0(sample_names, "_q=", qvec)

    set.seed(5)
    gene1_vals <- c(qvec * 1, qvec * 2) + rnorm(length(coln), sd = 1e-3)
    gene2_vals <- c(qvec * 1, qvec * 1) + rnorm(length(coln), sd = 1e-3)
    mat <- rbind(g1 = gene1_vals, g2 = gene2_vals)
    colnames(mat) <- coln
    rownames(mat) <- c("g1", "g2")

    rd <- data.frame(
        genes = rownames(mat),
        row.names = rownames(mat),
        stringsAsFactors = FALSE
    )
    cd <- data.frame(
        samples = sample_names,
        row.names = coln,
        stringsAsFactors = FALSE
    )

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )

    res_both <- calculate_lm_interaction(se, sample_type_col = "samples",
        method = "lmm", pvalue = "both", min_obs = 8)

    expect_s3_class(res_both, "data.frame")
    expect_true(all(c("gene", "p_interaction", "p_lrt", "p_satterthwaite", "adj_p_interaction") %in% colnames(res_both)))

    res_satt <- calculate_lm_interaction(se, sample_type_col = "samples",
        method = "lmm", pvalue = "satterthwaite", min_obs = 8)

    expect_true(all(c("gene", "p_interaction", "p_lrt", "p_satterthwaite") %in% colnames(res_satt)))
    # when pvalue = 'satterthwaite' the primary p_interaction should equal p_satterthwaite when available
    expect_true(all(is.na(res_satt$p_satterthwaite) | abs(res_satt$p_interaction - res_satt$p_satterthwaite) < 1e-12))

    res_lrt <- calculate_lm_interaction(se, sample_type_col = "samples",
        method = "lmm", pvalue = "lrt", min_obs = 8)

    expect_true(all(c("gene", "p_interaction", "p_lrt") %in% colnames(res_lrt)))
    # when pvalue = 'lrt' primary p_interaction should equal p_lrt when available
    expect_true(all(is.na(res_lrt$p_lrt) | abs(res_lrt$p_interaction - res_lrt$p_lrt) < 1e-12))
})
