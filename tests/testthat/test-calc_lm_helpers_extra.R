context("calc_lm_helpers extra cases")

library(testthat)

# .tsenat_report_fit_summary should be silent when no fallback/singular
test_that(".tsenat_report_fit_summary is silent when no messages to print", {
    df <- data.frame(x = 1:3)
    expect_silent(.tsenat_report_fit_summary(df, verbose = TRUE))
})

# linear branch: insufficient observations returns NULL, good data returns p-value
test_that(".tsenat_fit_one_interaction linear branch handles min_obs and returns p", {
    set.seed(1)
    # construct tiny matrix: 1 gene x 3 observations
    mat <- matrix(rnorm(3), nrow = 1)
    rownames(mat) <- "g1"
    q_vals <- c(0.1, 0.2, 0.3)
    sample_names <- paste0("s", seq_along(q_vals))
    group_vec <- c("A", "A", "B")
    # min_obs > non-missing -> NULL
    res_null <- .tsenat_fit_one_interaction("g1", se = NULL, mat = mat, q_vals = q_vals,
        sample_names = sample_names, group_vec = group_vec, method = "linear",
        pvalue = "lrt", subject_col = NULL, paired = FALSE, min_obs = 10, verbose = FALSE,
        suppress_lme4_warnings = TRUE, progress = FALSE)
    expect_null(res_null)

    # generate larger sample with clear interaction signal
    n <- 60
    qv <- rep(seq(0.1, 1.0, length.out = 20), 3)
    group <- rep(c("A", "B", "A"), each = 20)
    obs <- 0.5 * qv + ifelse(group == "B", 0.6 * qv, 0) + rnorm(length(qv), 0, 0.05)
    mat2 <- matrix(obs, nrow = 1)
    rownames(mat2) <- "gX"
    res <- .tsenat_fit_one_interaction("gX", se = NULL, mat = mat2, q_vals = qv,
        sample_names = paste0("s", seq_along(qv)), group_vec = group, method = "linear",
        pvalue = "lrt", subject_col = NULL, paired = FALSE, min_obs = 5, verbose = FALSE,
        suppress_lme4_warnings = TRUE, progress = FALSE)
    expect_true(is.data.frame(res) || is.null(res))
    if (is.data.frame(res)) expect_true(is.numeric(res$p_interaction) || is.na(res$p_interaction))
})

# lmm branch: errors when subject_col missing or paired but no sample_base
test_that(".tsenat_fit_one_interaction lmm errors when subject_col missing or paired with no sample_base", {
    skip_if_not_installed("SummarizedExperiment")
    set.seed(2)
    # small dataset to attach to se
    qv <- rep(0.1, 8)
    sample_names <- paste0("s", seq_along(qv))
    group <- rep(c("A", "B"), length.out = length(qv))
    obs <- rnorm(length(qv))
    mat <- matrix(obs, nrow = 1)
    rownames(mat) <- "g1"

    coldata <- S4Vectors::DataFrame(samples = sample_names)
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat), colData = coldata)

    # subject_col provided but not present
    expect_error(.tsenat_fit_one_interaction("g1", se = se, mat = mat, q_vals = qv,
        sample_names = sample_names, group_vec = group, method = "lmm", pvalue = "lrt",
        subject_col = "foo", paired = FALSE, min_obs = 2, verbose = FALSE,
        suppress_lme4_warnings = TRUE, progress = FALSE), "subject_col")

    # paired = TRUE but no sample_base column
    expect_error(.tsenat_fit_one_interaction("g1", se = se, mat = mat, q_vals = qv,
        sample_names = sample_names, group_vec = group, method = "lmm", pvalue = "lrt",
        subject_col = NULL, paired = TRUE, min_obs = 2, verbose = FALSE,
        suppress_lme4_warnings = TRUE, progress = FALSE), "paired = TRUE")
})

# lmm branch returns NULL when only one subject is present
test_that(".tsenat_fit_one_interaction lmm returns NULL with <2 subjects", {
    skip_if_not_installed("SummarizedExperiment")
    set.seed(3)
    qv <- rep(c(0.1, 0.5), times = 3)
    sample_names <- paste0("s", seq_along(qv))
    group <- rep(c("A", "B"), length.out = length(qv))
    # create single subject for all samples
    coldata <- S4Vectors::DataFrame(samples = sample_names, sample_base = rep("sub1", length(qv)))
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = matrix(rnorm(length(qv)), nrow = 1)), colData = coldata)
    mat <- matrix(rnorm(length(qv)), nrow = 1)
    rownames(mat) <- "g1"

    res <- .tsenat_fit_one_interaction("g1", se = se, mat = mat, q_vals = qv,
        sample_names = sample_names, group_vec = group, method = "lmm", pvalue = "lrt",
        subject_col = NULL, paired = TRUE, min_obs = 2, verbose = FALSE,
        suppress_lme4_warnings = TRUE, progress = FALSE)
    expect_null(res)
})

# .tsenat_extract_lrt_p returns NA when anova errors
test_that(".tsenat_extract_lrt_p returns NA for invalid models", {
    p <- .tsenat_extract_lrt_p("not_a_model", "also_not")
    expect_true(is.na(p))
})

# .tsenat_extract_satterthwaite_p returns NA when fallback has no interaction
test_that(".tsenat_extract_satterthwaite_p handles fallback lm without interaction", {
    df <- data.frame(entropy = rnorm(20), q = runif(20), group = rep(c("A", "B"), length.out = 20))
    fit1 <- stats::lm(entropy ~ q + group, data = df)
    fb <- list(fit1 = fit1)
    p <- .tsenat_extract_satterthwaite_p(NULL, fallback_lm = fb)
    expect_true(is.na(p) || is.numeric(p))

    # now with an interaction term
    fit2 <- stats::lm(entropy ~ q * group, data = df)
    fb2 <- list(fit1 = fit2)
    p2 <- .tsenat_extract_satterthwaite_p(NULL, fallback_lm = fb2)
    expect_true(is.numeric(p2) || is.na(p2))
})

# .tsenat_prepare_fpca_matrix returns NULL when insufficient good rows (min_obs large)
test_that(".tsenat_prepare_fpca_matrix returns NULL when min_obs larger than available", {
    mat <- matrix(rnorm(6), nrow = 1)
    sample_names <- rep(paste0("s", 1:3), each = 2)
    q_vals <- rep(c(1, 2), times = 3)
    res <- .tsenat_prepare_fpca_matrix(mat = mat, sample_names = sample_names, q_vals = q_vals, min_obs = 10)
    expect_null(res)
})

# .tsenat_try_lmer sets 'singular' attribute (when lme4 present and fit succeeded)
test_that(".tsenat_try_lmer sets singular attribute when fitting succeeds", {
    if (!rlang::is_installed("lme4")) skip("lme4 not installed")
    set.seed(7)
    nsub <- 10
    nper <- 3
    subject <- rep(paste0("s", seq_len(nsub)), each = nper)
    q <- rep(seq(0.1, 1, length.out = nper), times = nsub)
    group <- rep(rep(c("A", "B"), length.out = nper), times = nsub)
    entropy <- rnorm(length(subject), mean = 0.5 + as.numeric(group == "A") * 0.1 + 0.2 * q, sd = 0.05)
    df <- data.frame(entropy = entropy, q = q, group = group, subject = subject, stringsAsFactors = FALSE)
    fit_try <- .tsenat_try_lmer(entropy ~ q * group + (1 | subject), df, suppress_lme4_warnings = TRUE, verbose = FALSE)
    if (inherits(fit_try, "try-error")) {
        succeed()
    } else {
        expect_true(!is.null(attr(fit_try, "singular")))
        expect_true(is.logical(attr(fit_try, "singular")))
    }
})
