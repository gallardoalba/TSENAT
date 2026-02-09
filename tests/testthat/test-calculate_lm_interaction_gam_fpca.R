library(testthat)

context("calculate_lm_interaction gam and fpca methods")

test_that("gam method attaches p_interaction to rowData when mgcv available", {
    skip_if_not_installed("mgcv")

    qvec <- seq(0.01, 0.1, by = 0.01)
    # define three sample IDs, each measured across all q values
    sample_ids <- rep(c("S1", "S2"), each = length(qvec))
    coln <- paste0(sample_ids, "_q=", rep(qvec, times = 2))

    set.seed(7)
    gene1_vals <- c(qvec * 1, qvec * 2) + rnorm(length(coln), sd = 1e-3)
    gene2_vals <- c(qvec * 1, qvec * 1) + rnorm(length(coln), sd = 1e-3)
    mat <- rbind(g1 = gene1_vals, g2 = gene2_vals)
    colnames(mat) <- coln
    rownames(mat) <- c("g1", "g2")

    rd <- data.frame(genes = rownames(mat), row.names = rownames(mat), stringsAsFactors = FALSE)
    # sample-type mapping: first sample is Normal, second is Tumor
    cd <- data.frame(samples = sample_ids, sample_type = rep(c("Normal", "Tumor"), each = length(qvec)), row.names = coln, stringsAsFactors = FALSE)

    se <- SummarizedExperiment::SummarizedExperiment(assays = list(diversity = mat), rowData = rd, colData = cd)

    res_se <- calculate_lm_interaction(se, sample_type_col = "sample_type", method = "gam", min_obs = 8)
    expect_s4_class(res_se, "SummarizedExperiment")
    rd_out <- SummarizedExperiment::rowData(res_se)
    expect_true("p_interaction" %in% colnames(rd_out))
})

test_that("fpca method attaches p_interaction to rowData", {
    # Use 4 samples and 5 q values to ensure adequate FPCA matrix coverage
    qvec <- seq(0.01, 0.05, by = 0.01)
    sample_ids <- rep(c("S1", "S2", "S3", "S4"), each = length(qvec))
    coln <- paste0(sample_ids, "_q=", rep(qvec, times = 4))

    set.seed(9)
    # create clear differences so PCA can pick up group separation
    gene1_vals <- c(qvec * 1, qvec * 2, qvec * 1.5, qvec * 2.2) + rnorm(length(coln), sd = 1e-4)
    gene2_vals <- c(qvec * 1, qvec * 1, qvec * 1, qvec * 1) + rnorm(length(coln), sd = 1e-4)
    mat <- rbind(g1 = gene1_vals, g2 = gene2_vals)
    colnames(mat) <- coln
    rownames(mat) <- c("g1", "g2")

    rd <- data.frame(genes = rownames(mat), row.names = rownames(mat), stringsAsFactors = FALSE)
    cd <- data.frame(samples = sample_ids, sample_type = rep(c("Normal", "Tumor", "Normal", "Tumor"), each = length(qvec)), row.names = coln, stringsAsFactors = FALSE)

    se <- SummarizedExperiment::SummarizedExperiment(assays = list(diversity = mat), rowData = rd, colData = cd)

    # lower min_obs so test is robust; 4 samples with full q coverage should pass
    res_se <- calculate_lm_interaction(se, sample_type_col = "sample_type", method = "fpca", min_obs = 2)
    expect_s4_class(res_se, "SummarizedExperiment")
    rd_out <- SummarizedExperiment::rowData(res_se)
    expect_true("p_interaction" %in% colnames(rd_out))
})

test_that("paired lmm with subject_col attaches results when lme4 available", {
    skip_if_not_installed("lme4")

    qvec <- seq(0.01, 0.05, by = 0.01)
    sample_ids <- rep(c("S1", "S2", "S3"), each = length(qvec))
    coln <- paste0(sample_ids, "_q=", rep(qvec, times = 3))

    set.seed(11)
    gene1_vals <- c(qvec * 1, qvec * 2, qvec * 1.2) + rnorm(length(coln), sd = 1e-3)
    mat <- rbind(g1 = gene1_vals)
    colnames(mat) <- coln
    rownames(mat) <- c("g1")

    rd <- data.frame(genes = rownames(mat), row.names = rownames(mat), stringsAsFactors = FALSE)
    # build colData with sample identifiers per column and subject ids (sample_base)
    sample_ids <- rep(c("S1", "S2", "S3"), each = length(qvec))
    cd <- data.frame(samples = sample_ids, sample_type = rep(c("Normal", "Tumor", "Normal"), each = length(qvec)), sample_base = sample_ids, row.names = coln, stringsAsFactors = FALSE)

    se <- SummarizedExperiment::SummarizedExperiment(assays = list(diversity = mat), rowData = rd, colData = cd)

    res_se <- calculate_lm_interaction(se, sample_type_col = "sample_type", method = "lmm", subject_col = "sample_base", min_obs = 3)
    expect_s4_class(res_se, "SummarizedExperiment")
    rd_out <- SummarizedExperiment::rowData(res_se)
    expect_true("p_interaction" %in% colnames(rd_out))
})
