# =============================================================================
# test-bootstrap-coverage.R — Tsallis entropy bootstrap CI coverage.
#
# Design: the package bootstrap resamples reads (multinomial) from the
# observed count vector. Two criteria:
#   1. Variance calibration: sd(bootstrap dist.) / sd(real sampling dist.
#      of the estimator) must include 1 (interval [0.85, 1.15]).
#   2. Population coverage at sufficient depth (N=100000, B=2000): the 95%
#      percentile CI must cover H(p_true) at a frequency that includes 0.95.
# Evidence: coverage at moderate depth (N=5000), which is lower — read-level
# resampling treats p_hat as the population and undercovers H(p_true) at
# typical depths; it improves with depth.
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

calc_H <- getFromNamespace(".calculate_tsallis_entropy", "TSENAT")
calc_boot <- getFromNamespace(".calculate_tsallis_entropy_bootstrap", "TSENAT")

n_sims <- as.integer(Sys.getenv("TSENAT_BOOT_SIMS", "200"))
n_boot <- as.integer(Sys.getenv("TSENAT_BOOT_NBOOT", "2000"))
n_mc <- as.integer(Sys.getenv("TSENAT_BOOT_MC", "100"))
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

cov_pop <- sd_ratio <- cov_bca <- numeric(n_sims)
for (i in seq_len(n_sims)) {
    set.seed(seed * 1000 + i)  # per-iteration seed: stable against RNG consumption changes
    p_true <- rgamma(200, shape = 0.5)
    p_true <- p_true/sum(p_true)
    x <- as.numeric(rmultinom(1, size = 100000, prob = p_true))

    r <- calc_boot(x, q = 2, norm = TRUE, nboot = n_boot, ci = 0.95,
        method = "percentile", what = "S", include_diagnostics = FALSE,
        verbose = FALSE, show_messages = FALSE, pseudocount = 0)

    H_pop <- as.numeric(calc_H(p_true, q = 2, norm = TRUE, what = "S"))
    cov_pop[i] <- H_pop >= r$lower_ci && H_pop <= r$upper_ci

    # Evidence: BCa CI (bias/skewness correction) on the same data
    r_bca <- calc_boot(x, q = 2, norm = TRUE, nboot = n_boot, ci = 0.95,
        method = "bca", what = "S", include_diagnostics = FALSE,
        verbose = FALSE, show_messages = FALSE, pseudocount = 0)
    cov_bca[i] <- H_pop >= r_bca$lower_ci && H_pop <= r_bca$upper_ci

    # Real sampling SD of the estimator (multinomial MC) vs bootstrap SD
    hs <- replicate(n_mc, as.numeric(calc_H(
        as.numeric(rmultinom(1, 100000, p_true))/100000,
        q = 2, norm = TRUE, what = "S")))
    sd_ratio[i] <- stats::sd(r$bootstrap_dist, na.rm = TRUE)/stats::sd(hs)
}

# Evidencia: cobertura a profundidad moderada (N=5000)
cov_shallow <- numeric(n_sims)
for (i in seq_len(n_sims)) {
    set.seed(seed * 2000 + i)
    p_true <- rgamma(200, shape = 0.5)
    p_true <- p_true/sum(p_true)
    x <- as.numeric(rmultinom(1, size = 5000, prob = p_true))
    r <- calc_boot(x, q = 2, norm = TRUE, nboot = n_boot, ci = 0.95,
        method = "percentile", what = "S", include_diagnostics = FALSE,
        verbose = FALSE, show_messages = FALSE, pseudocount = 0)
    H_pop <- as.numeric(calc_H(p_true, q = 2, norm = TRUE, what = "S"))
    cov_shallow[i] <- H_pop >= r$lower_ci && H_pop <= r$upper_ci
}

ci_pop <- binom_ci(sum(cov_pop), n_sims)
ci_shallow <- binom_ci(sum(cov_shallow), n_sims)
ci_bca <- binom_ci(sum(cov_bca), n_sims)
sd_ci <- binom_ci(sum(sd_ratio >= 0.85 & sd_ratio <= 1.15), n_sims)

validation_summary <- data.frame(
    metric = c("Population coverage H(p_true), N=100000 (CI must include 0.95)",
        "BCa coverage H(p_true), N=100000 (evidence)",
        "Fraction with sd(boot)/sd(sampling) in [0.85,1.15] (>=0.90)",
        "Mean sd(boot)/sd(sampling) (evidence: ~1)",
        "Population coverage H(p_true), N=5000 (evidence: depth-limited)",
        "n_sims", "n_boot"),
    estimate = c(round(mean(cov_pop), 4), round(mean(cov_bca), 4),
        round(mean(sd_ratio >= 0.85 & sd_ratio <= 1.15), 4),
        round(mean(sd_ratio), 4),
        round(mean(cov_shallow), 4), n_sims, n_boot),
    ci_low = c(round(ci_pop["lower"], 4), round(ci_bca["lower"], 4),
        round(sd_ci["lower"], 4),
        NA_real_, round(ci_shallow["lower"], 4), NA_real_, NA_real_),
    ci_high = c(round(ci_pop["upper"], 4), round(ci_bca["upper"], 4),
        round(sd_ci["upper"], 4),
        NA_real_, round(ci_shallow["upper"], 4), NA_real_, NA_real_),
    criterion_met = c(contains_value(ci_pop, 0.95), TRUE,
        sd_ci["lower"] >= 0.90,
        TRUE, TRUE, TRUE, TRUE))
print(validation_summary)

# testthat integration: assert Monte Carlo criteria when running under testthat
if (requireNamespace("testthat", quietly = TRUE) && testthat::is_testing()) {
    testthat::test_that("test-bootstrap-coverage: Monte Carlo criteria met", {
        testthat::expect_true(all(validation_summary$criterion_met,
            na.rm = TRUE))
    })
}
