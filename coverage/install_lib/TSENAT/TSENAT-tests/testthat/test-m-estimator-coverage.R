# =============================================================================
# test-m-estimator-coverage.R — M-estimator robust CI coverage.
# Under H0 with contamination, the CI based on the M-estimator sandwich
# variance must have nominal coverage ≈ 0.95, while the CI based on the
# weighted RSS (IRLS) is reported as evidence (it may under/over-cover).
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

n_sims <- as.integer(Sys.getenv("TSENAT_MEST_SIMS", "500"))
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

#' Huber M-estimator of location (IRLS, mirror of the package) + sandwich SE
mest_location <- function(y) {
    mad_y <- median(abs(y - median(y)), na.rm = TRUE)
    scale_local <- 1.345 * mad_y
    if (scale_local == 0) scale_local <- 1
    loc <- median(y)
    for (iter in 1:50) {
        u <- (y - loc)/scale_local
        w <- ifelse(abs(u) <= 1, 1, 1/abs(u))
        loc_new <- sum(w * y)/sum(w)
        if (abs(loc_new - loc) < 1e-06) break
        loc <- loc_new
    }
    r <- y - loc
    u <- r/scale_local
    w_final <- ifelse(abs(u) <= 1, 1, 1/abs(u))
    # Sandwich: psi on the residual scale, psi(u)=u for |u|<=1, sign(u) otherwise
    psi <- ifelse(abs(r) <= scale_local, r, scale_local * sign(r))
    psip <- ifelse(abs(r) <= scale_local, 1, 0)
    se_sandwich <- sqrt(sum(psi^2)/(sum(psip)^2))
    # Package IRLS SE (evidence): sigma² = sum(w r²)/(n-2); se = sqrt(sigma²/n)
    sigma_sq <- sum(w_final * r^2)/max(1, length(y) - 2)
    se_irls <- sqrt(sigma_sq/length(y))
    c(loc = loc, se_sandwich = se_sandwich, se_irls = se_irls)
}

cov_sandwich <- cov_irls <- numeric(n_sims)
for (i in seq_len(n_sims)) {
    # n=200: the sandwich coverage converges to the nominal one (at n=30 there is
    # finite-sample undercoverage ~0.92, documented in the project notes)
    y <- rnorm(200)
    cont <- sample(200, 20)
    y[cont] <- rnorm(20, sd = 3)  # 10% contamination
    m <- mest_location(y)
    cov_sandwich[i] <- m["loc"] - 1.96 * m["se_sandwich"] <= 0 &&
        0 <= m["loc"] + 1.96 * m["se_sandwich"]
    df_t <- 200 - 1
    cov_irls[i] <- m["loc"] - qt(0.975, df_t) * m["se_irls"] <= 0 &&
        0 <= m["loc"] + qt(0.975, df_t) * m["se_irls"]
}

ci_s <- binom_ci(sum(cov_sandwich), n_sims)
ci_i <- binom_ci(sum(cov_irls), n_sims)

validation_summary <- data.frame(
    metric = c("Sandwich CI coverage (H0 + 10% contamination)",
        "IRLS weighted-RSS CI coverage (evidence)", "n_sims"),
    estimate = c(round(mean(cov_sandwich), 4), round(mean(cov_irls), 4), n_sims),
    ci_low = c(round(ci_s["lower"], 4), round(ci_i["lower"], 4), NA_real_),
    ci_high = c(round(ci_s["upper"], 4), round(ci_i["upper"], 4), NA_real_),
    criterion_met = c(contains_value(ci_s, 0.95), TRUE, TRUE))
print(validation_summary)

# testthat integration: assert Monte Carlo criteria when running under testthat
if (requireNamespace("testthat", quietly = TRUE) && testthat::is_testing()) {
    testthat::test_that("test-m-estimator-coverage: Monte Carlo criteria met", {
        testthat::expect_true(all(validation_summary$criterion_met,
            na.rm = TRUE))
    })
}
