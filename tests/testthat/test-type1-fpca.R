# =============================================================================
# test-type1-fpca.R — paired FPCA size calibration under H0.
# Package approach: Hotelling T² on the scores of the first principal
# components of the per-subject differences d_s(q) = A_s(q) - B_s(q);
# H0: the mean difference function is 0.
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

#' Paired FPCA (mirror of the package, sait_fpca.R): PCA on the POOLED
#' curves (subjects x conditions), scores per subject×condition, per-subject
#' score differences, one-sample Hotelling T² on the differences.
run_fpca <- function(df, max_pc = 3) {
    df$subject <- factor(df$subject)
    df$condition <- factor(df$condition)
    qs <- sort(unique(df$q))
    # Pooled matrix: one row per (subject, condition), q columns
    wide <- do.call(rbind, lapply(levels(df$subject), function(s) {
        do.call(rbind, lapply(levels(df$condition), function(cc) {
            v <- df$entropy[df$subject == s & df$condition == cc]
            v[order(qs)]
        }))
    }))  # (n_sub*2) x n_q
    n_sub <- length(levels(df$subject))
    n_q <- ncol(wide)
    K <- min(max_pc, n_q, n_sub - 1)
    pca <- prcomp(wide, center = TRUE, scale. = FALSE, rank. = K)
    scores <- pca$x  # (n_sub*2) x K
    # Per-subject difference: row A - row B
    D <- scores[seq(1, nrow(scores), by = 2), , drop = FALSE] -
        scores[seq(2, nrow(scores), by = 2), , drop = FALSE]
    n <- nrow(D)
    dbar <- colMeans(D)
    S <- cov(D)
    t2 <- if (K == 1) n * dbar^2/S[1, 1] else as.numeric(n * t(dbar) %*%
        solve(S) %*% dbar)
    if (!is.finite(t2) || t2 <= 0) return(1)
    f_stat <- t2 * (n - K)/(K * (n - 1))
    stats::pf(f_stat, df1 = K, df2 = n - K, lower.tail = FALSE)
}

pvals <- vapply(seq_len(reps), function(i) {
    df <- sim_paired_h0(n_subjects = 10, n_q = 8, rho = 0.5, seed = seed * 2000 +
        i)
    tryCatch(run_fpca(df), error = function(e) NA_real_)
}, FUN.VALUE = numeric(1))

k <- sum(pvals < 0.05, na.rm = TRUE)
n <- sum(!is.na(pvals))
ci <- binom_ci(k, n)

validation_summary <- data.frame(
    method = "FPCA paired (Hotelling T2 on PC scores)",
    n_valid = n, type1 = round(k/n, 4),
    ci_low = round(ci["lower"], 4), ci_high = round(ci["upper"], 4),
    includes_0.05 = contains_value(ci, 0.05))
print(validation_summary)

# testthat integration: assert Monte Carlo criteria when running under testthat
if (requireNamespace("testthat", quietly = TRUE) && testthat::is_testing()) {
    testthat::test_that("test-type1-fpca: Monte Carlo criteria met", {
        testthat::expect_true(all(validation_summary$includes_0.05,
            na.rm = TRUE))
    })
}
