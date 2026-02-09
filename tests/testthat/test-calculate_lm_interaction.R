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
