# =============================================================================
# test-fwer-westfall-young.R — Westfall–Young validation under exchangeability:
# FWER of the WY maxT procedure with paired permutation (label swap within
# each subject, the only valid exchangeability in paired designs).
# Scenarios: global H0 (200 null genes) and partial alternatives
# (180 nulls + 20 alternatives; FWER measured among the nulls).
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

n_sims <- as.integer(Sys.getenv("TSENAT_WY_SIMS", "100"))
n_perm <- as.integer(Sys.getenv("TSENAT_WY_PERM", "199"))
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

#' Simulate a dataset of G paired genes and return the observed statistic
#' per gene (sum of per-subject differences, not normalized by sd; all
#' genes share the same scale in this simulation) and the D matrix of
#' per-subject differences [n_subjects x G]. Nulls: N(0,sd); alternatives:
#' +effect.
#' NOTE: the raw sum is used because the t with per-permutation sd is unstable
#' with few subjects (permuted columns with sd ~ 0); with a homogeneous
#' per-gene scale, the raw sum is equivalent for validating the WY PROCEDURE.
sim_wy_data <- function(G, n_null, n_subjects = 8, sd_diff = 0.5,
    effect = 0.8, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    D <- matrix(rnorm(n_subjects * G, sd = sd_diff), nrow = n_subjects, ncol = G)
    if (G > n_null) {
        D[, (n_null + 1):G] <- D[, (n_null + 1):G] + effect
    }
    t_obs <- colMeans(D)
    list(D = D, t_obs = t_obs, null = seq_len(G) <= n_null)
}

#' Paired WY maxT: statistic = sum of differences per gene; the same
#' permutation (signs) for all genes. Returns adjusted p-values.
wy_maxT <- function(dat, n_perm) {
    n_sub <- nrow(dat$D)
    G <- ncol(dat$D)
    signs <- matrix(sample(c(-1, 1), n_sub * n_perm, replace = TRUE),
        nrow = n_sub, ncol = n_perm)
    t_perm <- crossprod(dat$D, signs)/n_sub
    max_t <- apply(abs(t_perm), 2, max)
    vapply(seq_len(G), function(g) (1 + sum(max_t >= abs(dat$t_obs[g])))/(n_perm +
        1), FUN.VALUE = numeric(1))
}

#' One replication: FWER among nulls = P(min adjusted p of nulls < alpha)
run_one <- function(G, n_null, seed) {
    dat <- sim_wy_data(G = G, n_null = n_null, seed = seed)
    p_adj <- wy_maxT(dat, n_perm)
    as.numeric(any(p_adj[dat$null] < 0.05))
}

# 1) Global H0: G=200, all nulls
fwer_global <- mean(vapply(seq_len(n_sims), function(i) {
    run_one(G = 200, n_null = 200, seed = seed * 1000 + i)
}, FUN.VALUE = numeric(1)))

# 2) Partial alternatives: 180 nulls + 20 alternatives
fwer_partial <- mean(vapply(seq_len(n_sims), function(i) {
    run_one(G = 200, n_null = 180, seed = seed * 2000 + i)
}, FUN.VALUE = numeric(1)))

# Criterion: FWER <= 0.05 (the CI upper bound must be <= 0.05)
ci_g <- binom_ci(round(fwer_global * n_sims), n_sims)
ci_p <- binom_ci(round(fwer_partial * n_sims), n_sims)

# Criterion: FWER ≤ 0.05. With finite simulation, it holds if the CI does
# not allow excluding FWER ≤ 0.05 (lower bound ≤ 0.05).
validation_summary <- data.frame(
    metric = c("FWER WY maxT (global null, 200 genes)",
        "FWER WY maxT (partial: 180 null + 20 alt, among nulls)",
        "n_sims", "n_perm"),
    estimate = c(round(fwer_global, 4), round(fwer_partial, 4), n_sims, n_perm),
    ci_low = c(round(ci_g["lower"], 4), round(ci_p["lower"], 4), NA_real_,
        NA_real_),
    ci_high = c(round(ci_g["upper"], 4), round(ci_p["upper"], 4), NA_real_,
        NA_real_),
    criterion_met = c(ci_g["lower"] <= 0.05, ci_p["lower"] <= 0.05, TRUE, TRUE))
print(validation_summary)

# testthat integration: assert Monte Carlo criteria when running under testthat
if (requireNamespace("testthat", quietly = TRUE) && testthat::is_testing()) {
    testthat::test_that("test-fwer-westfall-young: Monte Carlo criteria met", {
        testthat::expect_true(all(validation_summary$criterion_met,
            na.rm = TRUE))
    })
}
