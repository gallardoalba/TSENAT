context("plot_tsallis_gene_profile extra cases")

library(testthat)

# vector-of-genes input returns a named list of ggplots
test_that("plot_tsallis_gene_profile accepts vector of genes and returns list of ggplots", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("ggplot2")

    set.seed(101)
    readcounts <- matrix(rpois(30 * 4, lambda = 15), nrow = 30, ncol = 4)
    colnames(readcounts) <- c("S1_N", "S2_N", "S3_T", "S4_T")
    genes <- rep(paste0("G", 1:10), length.out = nrow(readcounts))

    qvals <- seq(0.01, 0.05, by = 0.02)
    ts_se <- calculate_diversity(readcounts, genes, q = qvals, norm = TRUE)

    coldata_df <- data.frame(
        Sample = c("S1_N", "S2_N", "S3_T", "S4_T"),
        Condition = c("Normal", "Normal", "Tumor", "Tumor"),
        stringsAsFactors = FALSE
    )

    ts_se <- map_metadata(ts_se, coldata_df)

    plots <- plot_tsallis_gene_profile(ts_se, gene = c("G1", "G2"))
    expect_type(plots, "list")
    expect_equal(length(plots), 2)
    expect_true(all(c("G1", "G2") %in% names(plots)))
    lapply(plots, function(p) expect_s3_class(p, "ggplot"))
})

# lm_res supplied with gene = NULL picks top n_top genes
test_that("plot_tsallis_gene_profile uses lm_res when gene is NULL and returns up to n_top plots", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("ggplot2")

    set.seed(202)
    # Build a SummarizedExperiment with a controlled interaction for first few genes
    base_samples <- paste0("S", 1:6)
    qvals <- c(0.1, 0.5)
    # column names like S1_q=0.1, S1_q=0.5, S2_q=0.1, S2_q=0.5, ...
    sample_cols <- as.character(t(outer(base_samples, qvals, function(s, q) paste0(s, "_q=", q))))

    n_genes <- 50
    mat <- matrix(rnorm(n_genes * length(sample_cols), sd = 0.1), nrow = n_genes, ncol = length(sample_cols))
    rownames(mat) <- paste0("G", seq_len(n_genes))
    colnames(mat) <- sample_cols

    # create a strong q:group interaction for first 5 genes: tumor samples have slope with q
    conds <- rep(c("Normal", "Tumor"), length.out = length(base_samples))
    for (g in 1:5) {
        for (j in seq_along(sample_cols)) {
            base_idx <- ((j - 1) %/% length(qvals)) + 1
            qv <- qvals[((j - 1) %% length(qvals)) + 1]
            cond <- conds[base_idx]
            mat[g, j] <- 0.2 + (cond == "Tumor") * (1.0 * qv) + rnorm(1, 0, 0.05)
        }
    }

    se <- SummarizedExperiment::SummarizedExperiment(assays = list(diversity = mat))
    coldata_df <- data.frame(
        Sample = base_samples,
        Condition = conds,
        stringsAsFactors = FALSE
    )
    ts_se <- map_metadata(se, coldata_df)

    # Run linear lm interaction (fast) to get lm_res
    lm_res <- calculate_lm_interaction(ts_se, sample_type_col = "sample_type", method = "linear", pvalue = "lrt", min_obs = 2)
    if (nrow(lm_res) == 0) {
        # fallback synthetic lm_res: ensure top genes exist for plotting
        pvals <- runif(50, min = 0.01, max = 1)
        pvals[1:5] <- runif(5, min = 0, max = 1e-6)
        lm_res <- data.frame(gene = paste0("G", seq_len(50)), p_interaction = pvals, adj_p_interaction = stats::p.adjust(pvals, method = "BH"), stringsAsFactors = FALSE)
    }
    # request top 5 instead of default 10 for speed/stability
    plots <- plot_tsallis_gene_profile(ts_se, gene = NULL, lm_res = lm_res, n_top = 5)
    expect_type(plots, "list")
    expect_lte(length(plots), 5)
    # each element must be ggplot
    lapply(plots, function(p) expect_s3_class(p, "ggplot"))
})
