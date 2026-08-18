# =============================================================================
# test-effect-sizes.R — divergence estimator bias and coverage.
# calculate_effect_sizes() reports the per-q divergence estimates of the
# LMM/divergence pipeline.
#
# Design: p_A ~ Dirichlet, p_B = mixture of p_A and another Dirichlet (known
# divergence). K replicates per group ~ Multinomial(N=50000, p). Package
# estimator: .tsallis_divergence_scalar on pooled counts with pseudocount 0.5.
# Truth (same estimand, count scale): .tsallis_divergence_scalar(p*N0,
# p*N0, pc=0.5) with N0 = K*50000 (the pseudocount acts on counts, not on
# probabilities). Criteria: (1) mean bias ≈ 0 (CI includes 0); (2) percentile
# bootstrap CI coverage includes 0.95 at K=32; evidence: coverage at K=8
# (lower — read resampling over pooled counts ignores part of the
# between-replicate variability).
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
ts_div <- getFromNamespace(".tsallis_divergence_scalar", "TSENAT")
ts_div_boot <- getFromNamespace(".calculate_divergence_bootstrap", "TSENAT")

n_sims <- as.integer(Sys.getenv("TSENAT_EFF_SIMS", "150"))
n_boot <- as.integer(Sys.getenv("TSENAT_EFF_NBOOT", "500"))
n_iso <- 100L
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

one_sim <- function(K, depth = 50000) {
    p_a <- rgamma(n_iso, shape = 1)
    p_a <- p_a/sum(p_a)
    p_shift <- rgamma(n_iso, shape = 1)
    p_shift <- p_shift/sum(p_shift)
    p_b <- (0.75 * p_a + 0.25 * p_shift)
    p_b <- p_b/sum(p_b)
    x <- rowSums(replicate(K, rmultinom(1, depth, p_a)))
    y <- rowSums(replicate(K, rmultinom(1, depth, p_b)))
    N0 <- K * depth
    truth <- ts_div(p_a * N0, p_b * N0, q_val = 1, pseudocount = 0.5)
    est <- ts_div(x, y, q_val = 1, pseudocount = 0.5)
    r <- ts_div_boot(x, y, q = 1, nboot = n_boot, ci = 0.95,
        method = "percentile", pseudocount = 0.5)
    c(truth = truth, est = est, lower = r$lower_ci, upper = r$upper_ci)
}

# Main criterion: K = 32 replicates per group
bias <- cov <- truth_vec <- numeric(n_sims)
for (i in seq_len(n_sims)) {
    r <- one_sim(K = 32)
    bias[i] <- r["est"] - r["truth"]
    cov[i] <- r["truth"] >= r["lower"] && r["truth"] <= r["upper"]
    truth_vec[i] <- r["truth"]
}

# Evidence: K = 8 replicates (lower pooled depth)
bias8 <- cov8 <- numeric(n_sims)
for (i in seq_len(n_sims)) {
    r <- one_sim(K = 8)
    bias8[i] <- r["est"] - r["truth"]
    cov8[i] <- r["truth"] >= r["lower"] && r["truth"] <= r["upper"]
}

ci_cov <- binom_ci(sum(cov), n_sims)
ci_cov8 <- binom_ci(sum(cov8), n_sims)

rel_bias_crit <- abs(mean(bias)) < 0.01 * mean(truth_vec)
validation_summary <- data.frame(
    metric = c("Relative bias |mean bias|/mean(truth) < 0.01 (K=32)",
        "Bootstrap CI coverage of truth (K=32) — must include 0.95",
        "Bootstrap CI coverage (K=8) — evidence: shallower",
        "n_sims", "n_boot"),
    estimate = c(round(mean(bias)/mean(truth_vec), 5), round(mean(cov), 4),
        round(mean(cov8), 4), n_sims, n_boot),
    ci_low = c(round((mean(bias) - 1.96 * stats::sd(bias)/sqrt(n_sims))/mean(truth_vec), 5),
        round(ci_cov["lower"], 4), round(ci_cov8["lower"], 4), NA_real_, NA_real_),
    ci_high = c(round((mean(bias) + 1.96 * stats::sd(bias)/sqrt(n_sims))/mean(truth_vec), 5),
        round(ci_cov["upper"], 4), round(ci_cov8["upper"], 4), NA_real_, NA_real_),
    criterion_met = c(rel_bias_crit, contains_value(ci_cov, 0.95), TRUE, TRUE,
        TRUE))
print(validation_summary)

# testthat integration: assert Monte Carlo criteria when running under testthat
if (requireNamespace("testthat", quietly = TRUE) && testthat::is_testing()) {
    testthat::test_that("test-effect-sizes: Monte Carlo criteria met", {
        testthat::expect_true(all(validation_summary$criterion_met,
            na.rm = TRUE))
    })
}
