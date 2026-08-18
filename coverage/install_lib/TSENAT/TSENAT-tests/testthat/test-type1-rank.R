# =============================================================================
# test-type1-rank.R — paired Conover-Iman rank-transform size calibration
# under H0. Mirrors the package implementation
# (.test_q_condition_interaction, method='rt'): ranks within each subject +
# lm(ranks ~ q_f * condition + subject) + analytical F-test of the interaction.
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

reps <- get0("TSENAT_VALIDATION_REPS", envir = environment(),
    ifnotfound = get0("TSENAT_VALIDATION_REPS", envir = .GlobalEnv,
        ifnotfound = as.integer(Sys.getenv("TSENAT_VALIDATION_REPS", "300"))))
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

#' Paired Conover-Iman (package implementation): ranks within subject,
#' two-way ANOVA blocked by subject, F of the q:condition interaction.
run_rt <- function(df) {
    df$subject <- factor(df$subject)
    df$condition <- factor(df$condition)
    df$q_f <- factor(df$q)
    df$ranks <- ave(df$entropy, df$subject, FUN = function(x) rank(x,
        na.last = "keep"))
    m <- lm(ranks ~ q_f * condition + subject, data = df)
    an <- anova(m)
    an$`Pr(>F)`[nrow(an) - 1]
}

pvals <- vapply(seq_len(reps), function(i) {
    df <- sim_paired_h0(n_subjects = 10, n_q = 8, rho = 0.5, seed = seed * 1000 +
        i)
    tryCatch(run_rt(df), error = function(e) NA_real_)
}, FUN.VALUE = numeric(1))

k <- sum(pvals < 0.05, na.rm = TRUE)
n <- sum(!is.na(pvals))
ci <- binom_ci(k, n)

validation_summary <- data.frame(
    method = "Conover-Iman paired (within-subject ranks)",
    n_valid = n, type1 = round(k/n, 4),
    ci_low = round(ci["lower"], 4), ci_high = round(ci["upper"], 4),
    includes_0.05 = contains_value(ci, 0.05))
print(validation_summary)

# testthat integration: assert Monte Carlo criteria when running under testthat
if (requireNamespace("testthat", quietly = TRUE) && testthat::is_testing()) {
    testthat::test_that("test-type1-rank: Monte Carlo criteria met", {
        testthat::expect_true(all(validation_summary$includes_0.05,
            na.rm = TRUE))
    })
}
