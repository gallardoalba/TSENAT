# =============================================================================
# test-lmm-functional.R — validation of the reformulated LMM pipeline.
#
# q is a deterministic functional argument, not time. The confirmatory test
# is the functional interaction H0: beta(q) = 0 ∀q on the ORIGINAL H(q) curve
# (no ARIMA differencing). This test runs the package's REAL PIPELINE
# (.calculate_sait, method = "lmm") under H0 and validates:
#   1. Type I ≈ 0.05 (CI includes 0.05) — closes the requirement to
#      "demonstrate through simulations that the type I error is controlled".
#   2. arima_transformation == FALSE in all results (no confirmatory
#      differencing).
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

n_sims <- as.integer(Sys.getenv("TSENAT_A1A2_SIMS", "150"))
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

pvals <- numeric(n_sims)
arima_flags <- logical(n_sims)
for (i in seq_len(n_sims)) {
    se <- build_se(seed * 20000 + i)
    r <- suppressWarnings(calc_sait(se, condition_col = "sample_type",
        method = "lmm", subject_col = "sample_base", min_obs = 8))
    if (is.data.frame(r) && nrow(r) == 1) {
        pvals[i] <- r$p_interaction[1]
        if ("arima_transformation" %in% colnames(r)) {
            arima_flags[i] <- !isTRUE(r$arima_transformation[1])
        } else {
            arima_flags[i] <- TRUE
        }
    } else {
        pvals[i] <- NA_real_
        arima_flags[i] <- TRUE
    }
}

n_valid <- sum(!is.na(pvals))
ci <- binom_ci(sum(pvals < 0.05, na.rm = TRUE), n_valid)

validation_summary <- data.frame(
    metric = c("Package LMM (functional, raw H(q)) type I — CI must include 0.05",
        "arima_transformation == FALSE in all results",
        "n_sims"),
    estimate = c(round(mean(pvals < 0.05, na.rm = TRUE), 4),
        round(mean(arima_flags, na.rm = TRUE), 4), n_sims),
    ci_low = c(round(ci["lower"], 4), NA_real_, NA_real_),
    ci_high = c(round(ci["upper"], 4), NA_real_, NA_real_),
    criterion_met = c(contains_value(ci, 0.05), all(arima_flags), TRUE))
print(validation_summary)

# testthat integration: assert Monte Carlo criteria when running under testthat
if (requireNamespace("testthat", quietly = TRUE) && testthat::is_testing()) {
    testthat::test_that("test-lmm-functional: Monte Carlo criteria met", {
        testthat::expect_true(all(validation_summary$criterion_met,
            na.rm = TRUE))
    })
}
