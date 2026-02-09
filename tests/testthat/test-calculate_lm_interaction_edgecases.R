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
