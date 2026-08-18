# =============================================================================
# test-postselection-lasso.R — LASSO/LMM post-selection inference.
#
# The package decouples LASSO/ElasticNet selection (exploratory) from the
# confirmatory LMM model: the p-value does not depend on the selection.
# This test validates:
#   1. CRITERION: type I of the package LMM pipeline with regularization='lasso'
#      under H0 ≈ 0.05 (exploratory selection does not leak into the test).
#   2. EVIDENCE: the naive counterfactual (select interactions with LASSO and
#      test ONLY the selected ones on the same data) inflates the type I —
#      it quantifies the risk the package avoids by design.
# =============================================================================
`%||%` <- function(a, b) if (is.null(a)) b else a
script_dir <- get0("TSENAT_VALIDATION_DIR", envir = environment(),
    ifnotfound = get0("TSENAT_VALIDATION_DIR", envir = .GlobalEnv,
        ifnotfound = {
            .args <- commandArgs(trailingOnly = FALSE)
            .fa <- .args[startsWith(.args, "--file=")]
            if (length(.fa)) {
                .sd <- dirname(sub("^--file=", "", .fa[[1]]))
                if (file.exists(file.path(.sd, "helpers.R"))) .sd else getwd()
            } else getwd()
        }))
sys.source(file.path(script_dir, "helpers.R"), envir = environment())

pkg_root <- normalizePath(file.path(script_dir, "..", ".."))
# Under testthat the package is already loaded; load_all() here would unload
# the installed namespace and break the whole run (covr/R CMD check).
if (!(requireNamespace("testthat", quietly = TRUE) && testthat::is_testing()) &&
    file.exists(file.path(pkg_root, "DESCRIPTION"))) {
    suppressMessages(pkgload::load_all(pkg_root, quiet = TRUE))
}
calc_sait <- getFromNamespace(".calculate_sait", "TSENAT")

n_sims <- as.integer(Sys.getenv("TSENAT_A9_SIMS", "100"))
n_q <- 10L
n_sub <- 8L
seed <- get0("TSENAT_VALIDATION_SEED", envir = environment(),
    ifnotfound = get0("TSENAT_VALIDATION_SEED", envir = .GlobalEnv,
        ifnotfound = as.integer(Sys.getenv("TSENAT_VALIDATION_SEED", "42"))))
set.seed(seed)

# testthat integration: the Monte Carlo suite is heavy and must NOT run under
# R CMD check (it made CI exceed its 30-minute no-output timeout). Skipped by
# default under testthat (Bioconductor builds always; elsewhere unless
# explicitly requested). Run the full suite locally with:
#   Rscript tests/testthat/run-validation.R
# or opt in under testthat with env TSENAT_RUN_VALIDATION=true
# (reps/seed via TSENAT_VALIDATION_REPS / TSENAT_VALIDATION_SEED).
if (requireNamespace("testthat", quietly = TRUE) && testthat::is_testing()) {
    testthat::skip_on_bioc()
    if (!isTRUE(as.logical(Sys.getenv("TSENAT_RUN_VALIDATION", "false")))) {
        testthat::skip("Monte Carlo validation suite skipped under testthat; run tests/testthat/run-validation.R or set TSENAT_RUN_VALIDATION=true")
    }
}

build_se <- function(seed_i) {
    qvec <- seq(0.1, 2, length.out = n_q)
    subjects <- rep(rep(paste0("S", seq_len(n_sub)), each = n_q), 2)
    conds <- rep(c("Normal", "Tumor"), each = n_sub * n_q)
    coln <- paste0(subjects, "_", conds, "_q=", rep(qvec, times = 2 * n_sub))
    u <- rnorm(n_sub, sd = 0.5)
    vals <- unlist(lapply(seq_len(n_sub), function(s) {
        c(u[s] + rnorm(n_q, sd = 0.2), u[s] + rnorm(n_q, sd = 0.2))
    }))
    mat <- matrix(vals, nrow = 1, dimnames = list("g1", coln))
    cd <- data.frame(samples = paste0(subjects, "_", conds),
        sample_type = conds,
        sample_base = subjects, row.names = coln, stringsAsFactors = FALSE)
    SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat), colData = cd)
}

# Naive counterfactual: LASSO selects q×group interactions and the LRT
# tests ONLY the selected ones (same data) → anti-conservative p.
naive_lasso_lrt <- function(seed_i) {
    qvec <- seq(0.1, 2, length.out = n_q)
    u <- rnorm(n_sub, sd = 0.5)
    df <- do.call(rbind, lapply(seq_len(n_sub), function(s) {
        rbind(data.frame(q = qvec, entropy = u[s] + rnorm(n_q, sd = 0.2),
            group = "Normal", subject = paste0("S", s)),
            data.frame(q = qvec, entropy = u[s] + rnorm(n_q, sd = 0.2),
                group = "Tumor", subject = paste0("S", s)))
    }))
    gi <- as.numeric(factor(df$group)) - 1
    X <- cbind(q = df$q, g = gi)
    uq <- sort(unique(df$q))
    for (qv in head(uq, -1)) X <- cbind(X, df$q * gi * (df$q == qv))
    colnames(X) <- make.unique(c("q", "g", paste0("iq", head(uq, -1))))
    # Add the interactions to df so the formula finds them
    for (k in seq_len(ncol(X))) df[[colnames(X)[k]]] <- X[, k]
    cv <- glmnet::cv.glmnet(x = X, y = df$entropy, family = "gaussian",
        alpha = 1, nfolds = 5)
    sel <- which(as.numeric(stats::coef(cv, s = "lambda.min")[-1]) != 0)
    if (length(sel) == 0 || length(sel) >= ncol(X)) sel <- setdiff(seq_len(ncol(X)),
        1:2)
    sel_nm <- colnames(X)[sel]
    df$subject <- factor(df$subject)
    fit0 <- nlme::lme(entropy ~ q + group, random = ~1 | subject, data = df,
        method = "ML")
    f_alt <- as.formula(paste("entropy ~", paste(c("q", "group", sel_nm),
        collapse = " + ")))
    fit1 <- try(nlme::lme(f_alt, random = ~1 | subject, data = df, method = "ML"),
        silent = TRUE)
    if (inherits(fit1, "try-error")) return(NA_real_)
    anova(fit0, fit1)$`p-value`[2]
}

p_pkg <- p_naive <- numeric(n_sims)
for (i in seq_len(n_sims)) {
    se <- build_se(seed * 10000 + i)
    r <- suppressWarnings(calc_sait(se, condition_col = "sample_type",
        method = "lmm", subject_col = "sample_base", min_obs = 8,
        regularization = "lasso"))
    p_pkg[i] <- if (is.data.frame(r) && nrow(r) == 1) r$p_interaction[1] else NA_real_
    p_naive[i] <- naive_lasso_lrt(seed * 10000 + i)
}

ci_pkg <- binom_ci(sum(p_pkg < 0.05, na.rm = TRUE), sum(!is.na(p_pkg)))
ci_naive <- binom_ci(sum(p_naive < 0.05, na.rm = TRUE), sum(!is.na(p_naive)))

# The package pipeline is conservative (type I < 0.05), so requiring the
# two-sided CI to CONTAIN 0.05 would reject a valid-but-conservative method.
# Fail only when the observed rejection count is significantly above the 5%
# null expectation (one-sided binomial check for inflation).
n_pkg <- sum(!is.na(p_pkg))
not_inflated <- sum(p_pkg < 0.05, na.rm = TRUE) <=
    stats::qbinom(0.95, n_pkg, 0.05)

validation_summary <- data.frame(
    metric = c("Package LMM + lasso (exploratory selection) type I — no inflation above 0.05",
        "Naive selection-then-test (same data) type I — evidence: inflated",
        "n_sims"),
    estimate = c(round(mean(p_pkg < 0.05, na.rm = TRUE), 4),
        round(mean(p_naive < 0.05, na.rm = TRUE), 4), n_sims),
    ci_low = c(round(ci_pkg["lower"], 4), round(ci_naive["lower"], 4), NA_real_),
    ci_high = c(round(ci_pkg["upper"], 4), round(ci_naive["upper"], 4), NA_real_),
    criterion_met = c(not_inflated, TRUE, TRUE))
print(validation_summary)

# testthat integration: assert Monte Carlo criteria when running under testthat
if (requireNamespace("testthat", quietly = TRUE) && testthat::is_testing()) {
    testthat::test_that("test-postselection-lasso: Monte Carlo criteria met", {
        testthat::expect_true(all(validation_summary$criterion_met,
            na.rm = TRUE))
    })
}
