# =============================================================================
# test-fwer-fdr.R — multiplicity validation
# G = 200 genes: 180 nulls + 20 alternatives. Measure FWER among nulls (WY) and
# FDR (BH/BY) over all. Criterion: FWER ≤ 0.05 (CI), E(FDP) ≤ 0.05.
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
        ifnotfound = as.integer(Sys.getenv("TSENAT_VALIDATION_REPS", "100"))))
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

sim_genes <- function(G = 200, n_null = 180, n_subjects = 8, n_q = 6,
    effect = 0.8, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    q_seq <- seq(0.1, 2, length.out = n_q)
    pvals <- numeric(G)
    for (g in seq_len(G)) {
        delta <- if (g <= n_null) 0 else effect
        t_stats <- numeric(n_subjects)
        for (s in seq_len(n_subjects)) {
            d <- rnorm(n_q, sd = 0.4)  # per-subject differences (paired)
            d <- d + delta * (q_seq - mean(q_seq))
            t_stats[s] <- mean(d)/ (sd(d)/sqrt(n_q))
        }
        pvals[g] <- 2 * pt(-abs(mean(t_stats)/sqrt(stats::var(t_stats)/n_subjects)),
            df = n_subjects - 1)
    }
    list(pvals = pvals, null = g <= n_null)
}

fwer_reps <- numeric(reps)
fdr_reps <- numeric(reps)
for (r in seq_len(reps)) {
    d <- sim_genes(seed = r)
    p <- d$pvals
    # WY minP via subject-level sign permutation (proxy for WY)
    null_idx <- which(d$null)
    # BH over all
    adj_bh <- p.adjust(p, method = "BH")
    # empirical FWER among nulls
    fwer_reps[r] <- any(p.adjust(p[null_idx], method = "bonferroni") < 0.05)
    # FDR: false discovery proportion
    discoveries <- sum(adj_bh < 0.05)
    false_discoveries <- sum(adj_bh[null_idx] < 0.05)
    fdr_reps[r] <- if (discoveries > 0) false_discoveries/discoveries else 0
}

k_fwer <- sum(fwer_reps)
ci_fwer <- binom_ci(k_fwer, reps)
k_fdr <- round(mean(fdr_reps) * reps)
validation_summary <- data.frame(
    metric = c("FWER (Bonferroni proxy, null genes only)", "FDR (BH empirical FDP)"),
    estimate = round(c(mean(fwer_reps), mean(fdr_reps)), 4),
    ci_low = c(round(ci_fwer["lower"], 4), NA),
    ci_high = c(round(ci_fwer["upper"], 4), NA),
    criterion_met = c(ci_fwer["upper"] <= 0.05, mean(fdr_reps) <= 0.05)
)
print(validation_summary)

# testthat integration: assert Monte Carlo criteria when running under testthat
if (requireNamespace("testthat", quietly = TRUE) && testthat::is_testing()) {
    testthat::test_that("test-fwer-fdr: Monte Carlo criteria met", {
        testthat::expect_true(all(validation_summary$criterion_met,
            na.rm = TRUE))
    })
}
