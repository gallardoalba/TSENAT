# =============================================================================
# test-qic-auto-selection.R — corstr selection via QIC.
# corstr='auto' selects the structure on the SAME inference data; the
# recommendation is to pre-specify for confirmatory analysis (the package
# default is corstr='ar1') and use 'auto' as exploratory/sensitivity. This
# test measures the type I of the GEE pipeline with corstr='auto' under H0
# (evidence) and verifies that the pre-specified default controls (criterion).
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
select_gee <- getFromNamespace(".select_gee_correlation", "TSENAT")

n_sims <- as.integer(Sys.getenv("TSENAT_QIC_SIMS", "200"))
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

run_gee_corstr <- function(df, corstr) {
    df$subject <- factor(df$subject)
    df$group <- factor(df$condition)
    df <- df[order(df$subject, df$group, df$q), ]
    if (!requireNamespace("geepack", quietly = TRUE)) return(NA_real_)
    fit <- geepack::geeglm(entropy ~ q * group, id = subject, data = df,
        family = gaussian(), corstr = corstr)
    coefs <- stats::coef(fit)
    ia <- which(grepl("^q:group", names(coefs)))
    if (length(ia) == 0) return(NA_real_)
    V <- stats::vcov(fit)[ia, ia, drop = FALSE]
    w <- as.numeric(t(coefs[ia]) %*% solve(V) %*% coefs[ia])
    n_clusters <- length(unique(df$subject))
    f_stat <- w/length(ia)
    stats::pf(f_stat, df1 = length(ia), df2 = n_clusters - length(ia),
        lower.tail = FALSE)
}

run_gee_auto <- function(df) {
    # QIC selection on the same data (mirrors the package's corstr='auto')
    sel <- select_gee(df, entropy ~ q + group, entropy ~ q * group,
        subject = df$subject, criteria = "qic")
    best <- sel$best_corstr
    if (is.null(best) || is.na(best)) return(NA_real_)
    run_gee_corstr(df, best)
}

p_ar1 <- p_auto <- numeric(n_sims)
for (i in seq_len(n_sims)) {
    df <- sim_paired_h0(n_subjects = 10, n_q = 8, rho = 0.5,
        seed = seed * 6000 + i)
    p_ar1[i] <- suppressWarnings(run_gee_corstr(df, "ar1"))
    p_auto[i] <- suppressWarnings(run_gee_auto(df))
}

ci_ar1 <- binom_ci(sum(p_ar1 < 0.05, na.rm = TRUE), sum(!is.na(p_ar1)))
ci_auto <- binom_ci(sum(p_auto < 0.05, na.rm = TRUE), sum(!is.na(p_auto)))

validation_summary <- data.frame(
    metric = c("GEE corstr='ar1' (pre-specified default) type I — CI must include 0.05",
        "GEE corstr='auto' (QIC on same data) type I — evidence",
        "n_sims"),
    estimate = c(round(mean(p_ar1 < 0.05, na.rm = TRUE), 4),
        round(mean(p_auto < 0.05, na.rm = TRUE), 4), n_sims),
    ci_low = c(round(ci_ar1["lower"], 4), round(ci_auto["lower"], 4), NA_real_),
    ci_high = c(round(ci_ar1["upper"], 4), round(ci_auto["upper"], 4), NA_real_),
    criterion_met = c(contains_value(ci_ar1, 0.05), TRUE, TRUE))
print(validation_summary)

# testthat integration: assert Monte Carlo criteria when running under testthat
if (requireNamespace("testthat", quietly = TRUE) && testthat::is_testing()) {
    testthat::test_that("test-qic-auto-selection: Monte Carlo criteria met", {
        testthat::expect_true(all(validation_summary$criterion_met,
            na.rm = TRUE))
    })
}
