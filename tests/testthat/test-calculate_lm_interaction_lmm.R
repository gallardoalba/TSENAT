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

    res_se_both <- calculate_lm_interaction(se, sample_type_col = "samples",
        method = "lmm", pvalue = "both", min_obs = 8)

    expect_s4_class(res_se_both, "SummarizedExperiment")
    rd_both <- SummarizedExperiment::rowData(res_se_both)
    expect_true(all(c("p_interaction", "p_lrt", "p_satterthwaite", "adj_p_interaction") %in% colnames(rd_both)))

    res_se_satt <- calculate_lm_interaction(se, sample_type_col = "samples",
        method = "lmm", pvalue = "satterthwaite", min_obs = 8)

    rd_satt <- SummarizedExperiment::rowData(res_se_satt)
    expect_true(all(c("p_interaction", "p_lrt", "p_satterthwaite") %in% colnames(rd_satt)))
    # when pvalue = 'satterthwaite' the primary p_interaction should equal p_satterthwaite when available
    mask <- !is.na(rd_satt$p_satterthwaite)
    expect_true(all(is.na(rd_satt$p_satterthwaite) | abs(rd_satt$p_interaction[mask] - rd_satt$p_satterthwaite[mask]) < 1e-8))

    res_se_lrt <- calculate_lm_interaction(se, sample_type_col = "samples",
        method = "lmm", pvalue = "lrt", min_obs = 8)

    rd_lrt <- SummarizedExperiment::rowData(res_se_lrt)
    expect_true(all(c("p_interaction", "p_lrt") %in% colnames(rd_lrt)))
    mask2 <- !is.na(rd_lrt$p_lrt)
    expect_true(all(is.na(rd_lrt$p_lrt) | abs(rd_lrt$p_interaction[mask2] - rd_lrt$p_lrt[mask2]) < 1e-8))
})
