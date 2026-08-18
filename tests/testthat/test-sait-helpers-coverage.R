context("SAIT helpers: additional branch coverage")

# ============================================================================
# .test_residual_normality — extraction/edge branches
# ============================================================================

test_that(".test_residual_normality handles insufficient residuals", {
    m <- stats::lm(y ~ x, data = data.frame(x = 1:2, y = c(1, 2)))
    res <- TSENAT:::.test_residual_normality(m, "gam")
    expect_equal(res$test_status, "error")
    expect_equal(res$n_residuals, 2)
    expect_match(res$report, "Insufficient residuals")
})

test_that(".test_residual_normality handles non-gaussian gam (pearson residuals)", {
    set.seed(1)
    d <- data.frame(x = seq_len(10), y = exp(rnorm(10)))
    m <- suppressWarnings(mgcv::gam(y ~ s(x), family = Gamma(link = "log"), data = d))
    res <- TSENAT:::.test_residual_normality(m, "gam")
    expect_true(res$test_status %in% c("pass", "fail"))
    expect_true(is.finite(res$shapiro_p_value) || is.na(res$shapiro_p_value))
})

test_that(".test_residual_normality emits verbose message on successful extraction", {
    set.seed(1)
    d <- data.frame(x = seq_len(10), y = rnorm(10))
    m <- mgcv::gam(y ~ s(x), data = d)
    expect_message(
        TSENAT:::.test_residual_normality(m, "gam", verbose = TRUE),
        "residuals appear"
    )
})

test_that(".test_residual_normality emits verbose message when extraction fails", {
    # Accessing $family on an atomic vector throws, which hits the
    # error-extraction branch of the tryCatch.
    expect_message(
        TSENAT:::.test_residual_normality(1:5, "gam", verbose = TRUE),
        "Could not extract residuals"
    )
})

test_that(".test_residual_normality handles shapiro.test execution failure", {
    # Constant residuals make shapiro.test() error -> execution-failed branch.
    m <- stats::lm(y ~ 1, data = data.frame(y = rep(3, 4)))
    res <- TSENAT:::.test_residual_normality(m, "gam")
    expect_equal(res$test_status, "error")
    expect_match(res$report, "Shapiro-Wilk test execution failed")
})

# ============================================================================
# .is_bounded_0_1 — empty and non-finite branches
# ============================================================================

test_that(".is_bounded_0_1 returns FALSE for empty input", {
    expect_false(TSENAT:::.is_bounded_0_1(numeric(0)))
})

test_that(".is_bounded_0_1 returns FALSE for all non-finite input", {
    expect_false(TSENAT:::.is_bounded_0_1(c(NA, NaN, Inf, -Inf)))
})

# ============================================================================
# .detect_heteroscedasticity — non-finite group variance + verbose
# ============================================================================

test_that(".detect_heteroscedasticity sets group var ratio NA for singleton groups", {
    df <- data.frame(entropy = c(0.5, 0.6), q = c(1, 2), group = c("A", "B"))
    res <- TSENAT:::.detect_heteroscedasticity(df, c(1, 2), c("A", "B"), verbose = TRUE)
    expect_true(is.na(res$var_ratio_group))
})

# ============================================================================
# .estimate_variance_weights — NULL / verbose / uniform / residual branches
# ============================================================================

test_that(".estimate_variance_weights returns NULL when OLS fails (power)", {
    df <- data.frame(q = 1:5, group = rep("a", 5))  # no 'entropy' column
    expect_null(
        TSENAT:::.estimate_variance_weights(df, q_vals = 1:5, method = "power")
    )
})

test_that(".estimate_variance_weights emits verbose power message", {
    set.seed(1)
    df <- data.frame(
        entropy = c(0.8, 0.7, 0.5, 0.4, 0.3, 0.2),
        q = rep(c(0.5, 1, 2), 2),
        group = rep(c("A", "B"), each = 3)
    )
    expect_message(
        TSENAT:::.estimate_variance_weights(df, q_vals = c(0.5, 1, 2), method = "power", verbose = TRUE),
        "power parameter"
    )
})

test_that(".estimate_variance_weights falls back to uniform weights", {
    df <- data.frame(q = 1:5)  # no 'entropy' -> residual OLS fails
    res <- TSENAT:::.estimate_variance_weights(df, q_vals = 1:5, method = "residual")
    expect_equal(res$method, "uniform")
    expect_equal(res$weights, rep(1, nrow(df)))
})

test_that(".estimate_variance_weights computes residual weights with group", {
    set.seed(1)
    df <- data.frame(
        entropy = c(0.8, 0.7, 0.5, 0.4, 0.3, 0.2),
        q = rep(c(0.5, 1, 2), 2),
        group = rep(c("A", "B"), each = 3)
    )
    res <- TSENAT:::.estimate_variance_weights(df, q_vals = c(0.5, 1, 2), method = "residual")
    expect_equal(res$method, "residual")
    expect_length(res$weights, nrow(df))
})

# ============================================================================
# .validate_sait_interaction_input — wy_randomizations + paired detection
# ============================================================================

test_that(".validate_sait_interaction_input defaults wy_randomizations to 1000", {
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = matrix(1:6, nrow = 2))
    )
    res <- TSENAT:::.validate_sait_interaction_input(
        method = "gam", pvalue = "Wald", corstr = "ar1", regularization = "pca",
        multicorr = "hochberg", pcorr = "BH", storey = FALSE,
        wy_randomizations = NULL, paired = FALSE, subject_col = NULL,
        se = se, verbose = FALSE
    )
    expect_equal(res$wy_randomizations, 1000)
})

test_that(".validate_sait_interaction_input rejects non-numeric wy_randomizations", {
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = matrix(1:6, nrow = 2))
    )
    expect_error(
        TSENAT:::.validate_sait_interaction_input(
            method = "gam", pvalue = "Wald", corstr = "ar1", regularization = "pca",
            multicorr = "hochberg", pcorr = "BH", storey = FALSE,
            wy_randomizations = TRUE, paired = FALSE, subject_col = NULL,
            se = se, verbose = FALSE
        ),
        "wy_randomizations must be numeric"
    )
})

test_that(".validate_sait_interaction_input auto-detects paired_samples", {
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = matrix(1:4, nrow = 2)),
        colData = S4Vectors::DataFrame(paired_samples = c("A", "B"), row.names = c("S1", "S2"))
    )
    expect_message(
        res <- TSENAT:::.validate_sait_interaction_input(
            method = "gam", pvalue = "Wald", corstr = "ar1", regularization = "pca",
            multicorr = "hochberg", pcorr = "BH", storey = FALSE,
            wy_randomizations = NULL, paired = TRUE, subject_col = NULL,
            se = se, verbose = TRUE
        ),
        "paired_samples"
    )
    expect_equal(res$subject_col, "paired_samples")
})

test_that(".validate_sait_interaction_input auto-detects sample_base", {
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = matrix(1:4, nrow = 2)),
        colData = S4Vectors::DataFrame(sample_base = c("A", "B"), row.names = c("S1", "S2"))
    )
    expect_message(
        res <- TSENAT:::.validate_sait_interaction_input(
            method = "gam", pvalue = "Wald", corstr = "ar1", regularization = "pca",
            multicorr = "hochberg", pcorr = "BH", storey = FALSE,
            wy_randomizations = NULL, paired = TRUE, subject_col = NULL,
            se = se, verbose = TRUE
        ),
        "sample_base"
    )
    expect_equal(res$subject_col, "sample_base")
})

test_that(".validate_sait_interaction_input errors for paired without pairing column", {
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = matrix(1:4, nrow = 2)),
        colData = S4Vectors::DataFrame(condition = c("A", "B"), row.names = c("S1", "S2"))
    )
    expect_error(
        TSENAT:::.validate_sait_interaction_input(
            method = "gam", pvalue = "Wald", corstr = "ar1", regularization = "pca",
            multicorr = "hochberg", pcorr = "BH", storey = FALSE,
            wy_randomizations = NULL, paired = TRUE, subject_col = NULL,
            se = se, verbose = FALSE
        ),
        "paired=TRUE requires"
    )
})

# ============================================================================
# .parse_sample_metadata — validation branches
# ============================================================================

test_that(".parse_sample_metadata errors when assay is missing", {
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:6, nrow = 2))
    )
    # assay() itself throws when the assay is absent, before the internal
    # is.null() guard is reached.
    expect_error(
        TSENAT:::.parse_sample_metadata(se, "condition", "diversity", FALSE),
        "diversity"
    )
})

test_that(".parse_sample_metadata errors when assay has no column names", {
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = matrix(1:6, nrow = 2))
    )
    expect_error(
        TSENAT:::.parse_sample_metadata(se, "condition", "diversity", FALSE),
        "No column names"
    )
})

test_that(".parse_sample_metadata errors when some columns lack '_q='", {
    mat <- matrix(1:4, nrow = 2)
    colnames(mat) <- c("S1_q=1", "S2")  # second column missing '_q='
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(diversity = mat))
    expect_error(
        TSENAT:::.parse_sample_metadata(se, "condition", "diversity", FALSE),
        "Some column names are missing '_q='"
    )
})

test_that(".parse_sample_metadata emits verbose parsing message", {
    qvec <- c(0.5, 1.0)
    coln <- c("S1_q=0.5", "S1_q=1.0", "S2_q=0.5", "S2_q=1.0")
    mat <- matrix(1:8, nrow = 2)
    colnames(mat) <- coln
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        colData = S4Vectors::DataFrame(condition = rep(c("A", "B"), each = 2), row.names = coln)
    )
    expect_message(
        TSENAT:::.parse_sample_metadata(se, "condition", "diversity", verbose = TRUE),
        "parsed samples and groups"
    )
})

# ============================================================================
# .adjust_pvalues_multicorr — verbose + validation branches
# ============================================================================

test_that(".adjust_pvalues_multicorr verbose hochberg", {
    expect_message(
        TSENAT:::.adjust_pvalues_multicorr(c(0.01, 0.02), multicorr = "hochberg", verbose = TRUE),
        "Hochberg"
    )
})

test_that(".adjust_pvalues_multicorr verbose BH", {
    expect_message(
        TSENAT:::.adjust_pvalues_multicorr(c(0.01, 0.02), multicorr = "bh", verbose = TRUE),
        "Benjamini-Hochberg"
    )
})

test_that(".adjust_pvalues_multicorr verbose Benjamini-Yekutieli", {
    expect_message(
        TSENAT:::.adjust_pvalues_multicorr(c(0.01, 0.02), multicorr = "benjamini-yekutieli", verbose = TRUE),
        "Benjamini-Yekutieli"
    )
})

test_that(".adjust_pvalues_multicorr errors on missing WY block_col", {
    p_values <- c(0.01, 0.02)
    metadata <- list(group_vec = c("A", "B", "A", "B"), q_vals = c(1, 2))
    mat <- matrix(1:8, nrow = 2, dimnames = list(
        c("g1", "g2"), c("S1_q=1", "S1_q=2", "S2_q=1", "S2_q=2")
    ))
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        colData = S4Vectors::DataFrame(condition = rep(c("A", "B"), each = 2), row.names = colnames(mat))
    )
    expect_error(
        TSENAT:::.adjust_pvalues_multicorr(
            p_values, multicorr = "westfall-young", wy_randomizations = 10,
            metadata = metadata, mat = mat, rownames_mat = rownames(mat), se = se,
            block_col = "bogus"
        ),
        "block_col"
    )
})

test_that(".adjust_pvalues_multicorr errors on missing WY strata_col", {
    p_values <- c(0.01, 0.02)
    metadata <- list(group_vec = c("A", "B", "A", "B"), q_vals = c(1, 2))
    mat <- matrix(1:8, nrow = 2, dimnames = list(
        c("g1", "g2"), c("S1_q=1", "S1_q=2", "S2_q=1", "S2_q=2")
    ))
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        colData = S4Vectors::DataFrame(condition = rep(c("A", "B"), each = 2), row.names = colnames(mat))
    )
    expect_error(
        TSENAT:::.adjust_pvalues_multicorr(
            p_values, multicorr = "westfall-young", wy_randomizations = 10,
            metadata = metadata, mat = mat, rownames_mat = rownames(mat), se = se,
            strata_col = "bogus"
        ),
        "strata_col"
    )
})

test_that(".adjust_pvalues_multicorr applies Storey verbose message", {
    skip_if_not_installed("fdrtool")
    expect_message(
        TSENAT:::.adjust_pvalues_multicorr(
            c(0.01, 0.02, 0.5), multicorr = "bh", verbose = TRUE, storey = TRUE
        ),
        "Storey"
    )
})

# ============================================================================
# .map_gene_annotations — verbose and fallback branches
# ============================================================================

test_that(".map_gene_annotations emits verbose annotation message", {
    rd <- S4Vectors::DataFrame(gene_name = c("A", "B"), row.names = c("g1", "g2"))
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = matrix(1:4, nrow = 2)),
        rowData = rd
    )
    res <- data.frame(gene = c("g1", "g2"))
    expect_message(
        TSENAT:::.map_gene_annotations(res, se, verbose = TRUE),
        "Gene annotations"
    )
})

test_that(".map_gene_annotations uses rownames as IDs when no id column", {
    rd <- S4Vectors::DataFrame(gene_name = c("A", "B"), row.names = c("g1", "g2"))
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = matrix(1:4, nrow = 2)),
        rowData = rd
    )
    res <- data.frame(gene = c("g1", "g2"))
    out <- TSENAT:::.map_gene_annotations(res, se, verbose = FALSE)
    expect_equal(out$gene_id, c("g1", "g2"))
    expect_equal(out$gene_name, c("A", "B"))
})

test_that(".map_gene_annotations emits mapping fallback message for unmapped genes", {
    rd <- S4Vectors::DataFrame(gene_name = c("A"), row.names = c("g1"))
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = matrix(1:2, nrow = 1)),
        rowData = rd
    )
    res <- data.frame(gene = c("g1", "g_unknown"))
    expect_message(
        TSENAT:::.map_gene_annotations(res, se, verbose = TRUE),
        "Gene mapping"
    )
})

test_that(".map_gene_annotations emits message when rowData lacks gene_name", {
    rd <- S4Vectors::DataFrame(some_col = c("x"), row.names = c("g1"))
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = matrix(1:2, nrow = 1)),
        rowData = rd
    )
    res <- data.frame(gene = c("g1"))
    expect_message(
        TSENAT:::.map_gene_annotations(res, se, verbose = TRUE),
        "gene_name column not found"
    )
})
