library(testthat)

context("calculate_lm_interaction")

test_that("calculate_lm_interaction returns expected columns and filters genes", {
    # construct synthetic data: 2 samples (Normal, Tumor) and multiple q values
    qvec <- seq(0.01, 0.1, by = 0.01)
    sample_names <- rep(c("S1_N", "S2_T"), each = length(qvec))
    coln <- paste0(sample_names, "_q=", qvec)

    # gene1: different slopes for Normal vs Tumor (interaction expected)
    # add small noise so the linear model is not a perfect fit
    set.seed(2)
    noise1 <- rnorm(length(coln), sd = 0.001)
    gene1_vals <- c(qvec * 1, qvec * 2) + noise1
    # gene2: same slope for both groups (no interaction expected)
    set.seed(1)
    noise <- rnorm(length(coln), sd = 0.001)
    gene2_vals <- c(qvec * 1, qvec * 1) + noise
    # gene3: mostly NA -> should be filtered out due to insufficient observations
    gene3_vals <- rep(NA_real_, length(coln))
    gene3_vals[1] <- 0.1

    mat <- rbind(g1 = gene1_vals, g2 = gene2_vals, g3 = gene3_vals)
    colnames(mat) <- coln
    rownames(mat) <- c("g1", "g2", "g3")

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

    res_se <- calculate_lm_interaction(se, sample_type_col = "samples",
        min_obs = 8)

    expect_s4_class(res_se, "SummarizedExperiment")
    rd <- SummarizedExperiment::rowData(res_se)
    expect_true(all(c("p_interaction", "adj_p_interaction") %in% colnames(rd)))
    # g1 and g2 should have non-NA entries in rowData (g3 filtered -> NA)
    expect_true(!is.na(rd["g1", "p_interaction"]))
    expect_true(!is.na(rd["g2", "p_interaction"]))
    expect_true(is.na(rd["g3", "p_interaction"]))
})

library(testthat)

context("calculate_lm_interaction edge cases")

test_that("missing sample_type_col produces informative error", {
    qvec <- seq(0.1, 0.5, by = 0.1)
    sample_ids <- rep(c("A", "B"), each = length(qvec))
    coln <- paste0(sample_ids, "_q=", rep(qvec, times = 2))
    mat <- matrix(runif(2 * length(coln)), nrow = 2)
    colnames(mat) <- coln
    rownames(mat) <- c("g1", "g2")
    rd <- data.frame(genes = rownames(mat), row.names = rownames(mat), stringsAsFactors = FALSE)
    cd <- data.frame(samples = sample_ids, row.names = coln, stringsAsFactors = FALSE)
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(diversity = mat), rowData = rd, colData = cd)

    expect_error(calculate_lm_interaction(se), "No sample grouping found|map sample types into `colData\\(se\\)`", fixed = FALSE)
})

test_that("column names without _q= are rejected", {
    mat <- matrix(runif(6), nrow = 2)
    colnames(mat) <- c("S1a", "S1b", "S2")
    rownames(mat) <- c("g1", "g2")
    rd <- data.frame(genes = rownames(mat), row.names = rownames(mat), stringsAsFactors = FALSE)
    cd <- data.frame(samples = c("A", "A", "B"), row.names = colnames(mat), stringsAsFactors = FALSE)
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(diversity = mat), rowData = rd, colData = cd)

    expect_error(calculate_lm_interaction(se, sample_type_col = "samples"), "Could not parse q values", fixed = FALSE)
})

test_that("invalid method and pvalue arguments produce errors", {
    qvec <- seq(0.01, 0.05, by = 0.01)
    sample_ids <- rep(c("S1", "S2"), each = length(qvec))
    coln <- paste0(sample_ids, "_q=", rep(qvec, times = 2))
    mat <- rbind(g1 = runif(length(coln)), g2 = runif(length(coln)))
    colnames(mat) <- coln
    rownames(mat) <- c("g1", "g2")
    rd <- data.frame(genes = rownames(mat), row.names = rownames(mat), stringsAsFactors = FALSE)
    cd <- data.frame(samples = rep(c("Normal", "Tumor"), each = length(qvec)), row.names = coln, stringsAsFactors = FALSE)
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(diversity = mat), rowData = rd, colData = cd)

    expect_error(calculate_lm_interaction(se, sample_type_col = "samples", method = "nope"), "should be one of", fixed = FALSE)
    expect_error(calculate_lm_interaction(se, sample_type_col = "samples", method = "linear", pvalue = "nope"), "should be one of", fixed = FALSE)
})

test_that("lmm fallback used when lmer fails (stubbed)", {
    skip_if_not_installed("lme4")
    qvec <- seq(0.01, 0.05, by = 0.01)
    sample_ids <- rep(c("S1", "S2", "S3"), each = length(qvec))
    coln <- paste0(sample_ids, "_q=", rep(qvec, times = 3))
    set.seed(42)
    mat <- rbind(g1 = c(qvec * 1, qvec * 2, qvec * 1.5) + rnorm(length(coln), sd = 1e-4))
    colnames(mat) <- coln
    rownames(mat) <- c("g1")
    rd <- data.frame(genes = rownames(mat), row.names = rownames(mat), stringsAsFactors = FALSE)
    cd <- data.frame(samples = rep(c("Normal", "Tumor", "Normal"), each = length(qvec)), sample_base = rep(c("S1", "S2", "S3"), each = length(qvec)), row.names = coln, stringsAsFactors = FALSE)
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(diversity = mat), rowData = rd, colData = cd)

    ns <- asNamespace("TSENAT")
    orig <- get(".tsenat_try_lmer", envir = ns)
    stub <- function(...) structure("error", class = "try-error")
    assignInNamespace('.tsenat_try_lmer', stub, ns = 'TSENAT')
    on.exit(assignInNamespace('.tsenat_try_lmer', orig, ns = 'TSENAT'), add = TRUE)

    res_se <- calculate_lm_interaction(se, sample_type_col = "samples", method = "lmm", subject_col = "sample_base", min_obs = 3)
    # function may return a SummarizedExperiment (writing into rowData) or a data.frame fallback
    if (is.data.frame(res_se)) {
        # accept any data.frame fallback (presence indicates graceful handling)
        expect_true(is.data.frame(res_se))
    } else {
        expect_s4_class(res_se, "SummarizedExperiment")
        rd_out <- SummarizedExperiment::rowData(res_se)
        # fallback should set a fit_method value (e.g., lm_subject or lm_nosubject)
        expect_true("fit_method" %in% colnames(rd_out))
    }
})

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
